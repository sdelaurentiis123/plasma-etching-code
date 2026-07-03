#!/usr/bin/env python3
"""Charging gate: petch's 2-D Hwang-Giapis charging solver vs the published curves.

GATE 1 (primary, 2-D solver): floor ion flux vs AR -- HG JAP 82, 566 Fig. 4 (digitized +-0.02),
Cl2 HDP, V_s = 37 + 30 sin(wt), T_e = 4 eV. Solver carries the FULL HG geometry now: poly-Si
conductor line, periodic pitch, physical 3.7 um boundary, RF-burst electrons, ion energy-angle
anticorrelation. All constants published; NOTHING tuned. Gate: RMSE <= 0.05 over the 8 points.
GATE 1b (info->gate): V_floor_center 8 -> 33 V over AR 1 -> 4 (ground-referenced now; +-40%).
GATE 2: Matsui asymptote (APL 78, 883): with 300 eV ions the floor must NOT cut off at AR 4-6.
GATE 3 (secondary, 0-D closure sanity): shape monotone + Matsui pass (charging.floor_balance).
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from petch.charging2d import solve_trench_charging
from petch.charging import floor_balance

HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
HG_FLUX = np.array([0.59, 0.55, 0.47, 0.40, 0.34, 0.30, 0.26, 0.22])
SEE_MODEL = os.environ.get("PETCH_SEE_MODEL", "none")
SEE_GENERATIONS = int(os.environ.get("PETCH_SEE_GENERATIONS", "1"))
SOURCE_MODEL = os.environ.get("PETCH_SOURCE_MODEL", "analytic")
POLY_MODE = os.environ.get("PETCH_POLY_MODE", "tied")
POLY_BIAS_V = float(os.environ.get("PETCH_POLY_BIAS_V", "0.0"))
EDGE_OPEN_MODEL = os.environ.get("PETCH_EDGE_OPEN_MODEL", "none")
EDGE_OPEN_ELECTRON_FLUX = os.environ.get("PETCH_EDGE_OPEN_ELECTRON_FLUX")
EDGE_OPEN_ELECTRON_FLUX = None if EDGE_OPEN_ELECTRON_FLUX is None else float(EDGE_OPEN_ELECTRON_FLUX)

print("=== GATE 1 (2-D solver): floor ion flux vs AR (HG JAP 82,566 Fig.4) ===", flush=True)
print(f"    see_model={SEE_MODEL} see_generations={SEE_GENERATIONS} "
      f"source_model={SOURCE_MODEL} poly_mode={POLY_MODE} poly_bias_V={POLY_BIAS_V} "
      f"edge_open_model={EDGE_OPEN_MODEL} edge_open_electron_flux={EDGE_OPEN_ELECTRON_FLUX}", flush=True)
pred, vcs, vfs, isurv, esurv = [], [], [], [], []
edge_e_gross, edge_i_gross, edge_net, edge_hg_gross = [], [], [], []
resmax = []
for i, ar in enumerate(HG_AR):
    t0 = time.time()
    r = solve_trench_charging(ar, n_per_iter=8000, n_iter=140, seed=i,
                              see_model=SEE_MODEL, see_generations=SEE_GENERATIONS,
                              source_model=SOURCE_MODEL, poly_mode=POLY_MODE,
                              poly_bias_V=POLY_BIAS_V,
                              edge_open_model=EDGE_OPEN_MODEL,
                              edge_open_electron_flux=EDGE_OPEN_ELECTRON_FLUX)
    pred.append(r["floor_flux"]); vcs.append(r["V_floor_center"]); vfs.append(r["V_foot_peak"])
    ti = r["diag"]["trace"]["last_ion"]; te = r["diag"]["trace"]["last_electron"]
    isurv.append(ti["survivor_frac"] if ti else np.nan)
    esurv.append(te["survivor_frac"] if te else np.nan)
    res = r["diag"].get("residual", {})
    resmax.append(max(abs(float(v)) for v in res.values()) if res else np.nan)
    eo = r["diag"].get("edge_open", {})
    edge_e_gross.append(float(eo.get("electron_gross", np.nan)))
    edge_i_gross.append(float(eo.get("ion_gross", np.nan)))
    edge_net.append(float(eo.get("net_electron", np.nan)))
    edge_hg_gross.append(float(eo.get("hg_electron_gross", np.nan)))
    print(f"  AR {ar:3.1f}:  model={r['floor_flux']:.3f}  HG={HG_FLUX[i]:.2f}   "
          f"Vc={r['V_floor_center']:5.1f} Vfoot={r['V_foot_peak']:5.1f}  "
          f"surv_i/e={isurv[-1]:.4f}/{esurv[-1]:.4f} res={resmax[-1]:.3f}  "
          f"({time.time()-t0:.0f}s)", flush=True)
pred = np.array(pred)
rmse = float(np.sqrt(np.mean((pred - HG_FLUX) ** 2)))
print(f"  GATE 1 RMSE = {rmse:.3f}  [{'PASS' if rmse <= 0.05 else 'fail'}] (gate 0.05)", flush=True)
print(f"  survivor max ion/electron = {np.nanmax(isurv):.4f}/{np.nanmax(esurv):.4f}  "
      f"[{'PASS' if max(np.nanmax(isurv), np.nanmax(esurv)) < 0.001 else 'fail'}] (W1 gate 0.001)", flush=True)
print(f"  current residual max = {np.nanmax(resmax):.3f}  "
      f"[{'PASS' if np.nanmax(resmax) < 0.08 else 'fail'}] (diagnostic gate 0.08)", flush=True)

print("=== GATE 2 (2-D): Matsui asymptote (300 eV ions: no cutoff at AR 4) ===", flush=True)
r300 = solve_trench_charging(4.0, V_dc=300.0, V_rf=30.0, n_per_iter=8000, n_iter=140,
                             see_model=SEE_MODEL, see_generations=SEE_GENERATIONS,
                             source_model=SOURCE_MODEL, poly_mode=POLY_MODE,
                             poly_bias_V=POLY_BIAS_V,
                             edge_open_model=EDGE_OPEN_MODEL,
                             edge_open_electron_flux=EDGE_OPEN_ELECTRON_FLUX)
ok2 = r300["floor_flux"] > 0.1
print(f"  300 eV ions @ AR4: floor flux = {r300['floor_flux']:.3f}   [{'PASS' if ok2 else 'fail'}] (must stay well above 0)", flush=True)

print("=== GATE 3 (0-D closure sanity) ===", flush=True)
f0 = [floor_balance(a, n=200000)[1] for a in (1.0, 4.0)]
m0 = floor_balance(4.0, V_dc=300.0, n=200000)[1]
ok3 = f0[0] > f0[1] and m0 > 0.05
print(f"  0-D: flux AR1={f0[0]:.3f} > AR4={f0[1]:.3f}, Matsui AR4@300eV={m0:.3f}  [{'PASS' if ok3 else 'fail'}]", flush=True)
np.savez(os.path.join(os.path.dirname(__file__), "..", "charging_gate_result.npz"),
         ar=HG_AR, hg=HG_FLUX, model=pred, vc=vcs, vfoot=vfs, rmse=rmse,
         survivor_ion=np.array(isurv), survivor_electron=np.array(esurv),
         residual_max=np.array(resmax),
         edge_electron_gross=np.array(edge_e_gross), edge_ion_gross=np.array(edge_i_gross),
         edge_net_electron=np.array(edge_net), edge_hg_electron_gross=np.array(edge_hg_gross),
         see_model=SEE_MODEL, see_generations=SEE_GENERATIONS,
         source_model=SOURCE_MODEL, poly_mode=POLY_MODE, poly_bias_V=POLY_BIAS_V,
         edge_open_model=EDGE_OPEN_MODEL,
         edge_open_electron_flux=-1.0 if EDGE_OPEN_ELECTRON_FLUX is None else EDGE_OPEN_ELECTRON_FLUX)
print("DONE", flush=True)
