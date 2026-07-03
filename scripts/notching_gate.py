#!/usr/bin/env python3
"""Notching mechanism gates: the deflected-ion population + poly-line potential vs Hwang-Giapis.

The notch is dug by ions deflected off the charged floor into the poly-Si sidewall foot (HG JVST B
15,70; JAP 82,566). The solver now carries HG's geometry: PR mask (insulator) over a 0.3 µm poly-Si
CONDUCTOR line (one floating equipotential, explicit charge redistribution) on oxide, periodic
line/space pitch, sheath boundary at the physical 3.7 µm, RF-burst electrons (residual Child-law
barrier), ion energy-angle anticorrelation. All constants published; NOTHING tuned.

  GATE A: mean impact energy of ions on the poly sidewall RISES 15 -> 27.5 eV over AR 1 -> 4
          (HG JAP 82,566; +-30%).
  GATE B: deflected-ion flux to the poly is ~AR-INDEPENDENT (max/min <= 2 over AR >= 1.6).
  GATE C: the poly-line floating potential rises ~6 -> 39 V over AR 1 -> 4 (+-30%) — the
          V_m - V_p gap that accelerates the notch ions (NEW: conductor physics).

Writes notching_gate_result.npz for viz/notching.png.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from petch.charging2d import solve_trench_charging

HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
HG_FOOT_E = np.array([15.0, 16.5, 17.5, 20.0, 23.0, 25.0, 26.5, 27.5])   # eV, poly-sidewall ions
HG_VPOLY = np.array([6.0, 9.0, 15.0, 20.0, 27.0, 31.0, 36.0, 39.0])       # V (6->39, Fig. 6 trend)
SEE_MODEL = os.environ.get("PETCH_SEE_MODEL", "none")
SEE_GENERATIONS = int(os.environ.get("PETCH_SEE_GENERATIONS", "1"))
SOURCE_MODEL = os.environ.get("PETCH_SOURCE_MODEL", "analytic")
POLY_MODE = os.environ.get("PETCH_POLY_MODE", "tied")
POLY_BIAS_V = float(os.environ.get("PETCH_POLY_BIAS_V", "0.0"))

Em, Fl, Vp, Vc, Fx, isurv, esurv = [], [], [], [], [], [], []
print(f"see_model={SEE_MODEL} see_generations={SEE_GENERATIONS} "
      f"source_model={SOURCE_MODEL} poly_mode={POLY_MODE} poly_bias_V={POLY_BIAS_V}", flush=True)
for i, ar in enumerate(HG_AR):
    r = solve_trench_charging(ar, n_per_iter=8000, n_iter=110, seed=3 + i,
                              see_model=SEE_MODEL, see_generations=SEE_GENERATIONS,
                              source_model=SOURCE_MODEL, poly_mode=POLY_MODE,
                              poly_bias_V=POLY_BIAS_V)
    Em.append(r["foot_ion_Emean"]); Fl.append(r["foot_ion_flux"])
    Vp.append(r["V_poly"]); Vc.append(r["V_floor_center"]); Fx.append(r["floor_flux"])
    ti = r["diag"]["trace"]["last_ion"]; te = r["diag"]["trace"]["last_electron"]
    isurv.append(ti["survivor_frac"] if ti else np.nan)
    esurv.append(te["survivor_frac"] if te else np.nan)
    print(f"  AR {ar:3.1f}:  footE={Em[-1]:5.1f} eV (HG {HG_FOOT_E[i]:.1f})   footFlux={Fl[-1]:.3f}   "
          f"Vpoly={Vp[-1]:5.1f} (HG {HG_VPOLY[i]:.0f})   Vc={Vc[-1]:5.1f}   floorFlux={Fx[-1]:.3f}   "
          f"surv_i/e={isurv[-1]:.4f}/{esurv[-1]:.4f}", flush=True)
Em = np.array(Em); Fl = np.array(Fl); Vp = np.array(Vp); Vc = np.array(Vc); Fx = np.array(Fx)

relE = np.abs(Em - HG_FOOT_E) / HG_FOOT_E
okA = bool((relE <= 0.30).all()) and Em[-1] > Em[0]
print(f"GATE A (foot ion energy 15->27.5 eV, +-30%): max rel err = {relE.max()*100:.0f}%  "
      f"rising = {Em[-1] > Em[0]}   [{'PASS' if okA else 'fail'}]")
sel = HG_AR >= 1.6
ratio = Fl[sel].max() / max(Fl[sel].min(), 1e-9)
okB = bool(ratio <= 2.0)
print(f"GATE B (foot flux ~AR-independent, AR>=1.6): max/min = {ratio:.2f}  [{'PASS' if okB else 'fail'}]")
relP = np.abs(Vp - HG_VPOLY) / HG_VPOLY
okC = bool((relP <= 0.30).all()) and Vp[-1] > Vp[0]
print(f"GATE C (poly potential 6->39 V, +-30%): max rel err = {relP.max()*100:.0f}%  "
      f"rising = {Vp[-1] > Vp[0]}   [{'PASS' if okC else 'fail'}]")
print(f"W1 survivor gate: max ion/electron = {np.nanmax(isurv):.4f}/{np.nanmax(esurv):.4f}  "
      f"[{'PASS' if max(np.nanmax(isurv), np.nanmax(esurv)) < 0.001 else 'fail'}] (gate 0.001)")
np.savez(os.path.join(os.path.dirname(__file__), "..", "notching_gate_result.npz"),
         ar=HG_AR, foot_E=Em, hg_foot_E=HG_FOOT_E, foot_flux=Fl, vpoly=Vp, hg_vpoly=HG_VPOLY,
         vc=Vc, floor_flux=Fx, okA=okA, okB=okB, okC=okC,
         survivor_ion=np.array(isurv), survivor_electron=np.array(esurv),
         see_model=SEE_MODEL, see_generations=SEE_GENERATIONS,
         source_model=SOURCE_MODEL, poly_mode=POLY_MODE, poly_bias_V=POLY_BIAS_V)
print("DONE")
