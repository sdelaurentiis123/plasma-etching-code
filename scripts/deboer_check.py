#!/usr/bin/env python3
"""petch vs the de Boer/Blauw cryo SF6/O2 DRIE RIE-lag EXPERIMENT (real wafer ARDE), with CURRENT
params (cal_F / Ysp_scale removed 2026-06-18 -- deboer_final.py is stale). Two configs:
  A) FAITHFUL = the exact config that matches ViennaPS (belen + viennaps yields + ion reflection).
  B) CRYO-like = collimated ions (sub-deg IADF) + high passivation, best-effort with current knobs.
Both seed-averaged, compared to experiment normARDE [1.0,0.43,0.29,0.20] @ AR [0,10,20,40]. cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W, DX, NS = 2.0, 0.25, 160
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)


def curve(par_over, flags_over, seeds=(0, 1)):
    arN = np.linspace(0.5, 30, 60); acc = np.zeros_like(arN); cnt = np.zeros_like(arN)
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1; p.update(par_over)
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         neutral_transport="mc", sampling="sobol", warm_start_coverage=True, **flags_over)
        g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                           reinit_method="fsm", verbose=False, record_depth_every=4, seed_offset=sd * 100, **GEO)
        h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
        tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
        r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
        nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
        samp = np.interp(arN, ar, nr, left=1.0, right=np.nan)
        ok = ~np.isnan(samp); acc[ok] += samp[ok]; cnt[ok] += 1
    return arN, np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan), float(dd.max() / W)


print(f"device={t3.DEVICE}  de Boer empirical check (seed-avg)...\n", flush=True)
arF, nrF, armaxF = curve(dict(), dict(ion_reflection=True))
arC, nrC, armaxC = curve(dict(ion_ang_sigma=np.deg2rad(0.8), betaE=0.9), dict(ion_reflection=True))
rF = np.interp(EXP_AR, arF, nrF); rC = np.interp(EXP_AR, arC, nrC)
rmseF = float(np.sqrt(np.nanmean((rF - EXP_R) ** 2))); rmseC = float(np.sqrt(np.nanmean((rC - EXP_R) ** 2)))
np.savez("/root/deboer_check.npz", arF=arF, nrF=nrF, arC=arC, nrC=nrC, exp_ar=EXP_AR, exp_r=EXP_R)
print("================  petch vs de Boer cryo SF6/O2 wafer experiment  ================", flush=True)
print(f"  de Boer EXPERIMENT      @AR 0/10/20/40 = {EXP_R}", flush=True)
print(f"  petch FAITHFUL (=ViennaPS cfg)         = {np.round(rF,3)}  RMSE {rmseF:.3f}  (reached AR~{armaxF:.0f})", flush=True)
print(f"  petch CRYO-like (0.8deg IADF, betaE.9) = {np.round(rC,3)}  RMSE {rmseC:.3f}  (reached AR~{armaxC:.0f})", flush=True)
print("\n  NOTE: petch+ViennaPS are ballistic; de Boer is a real cryo wafer. Gap at high AR is the", flush=True)
print("  structural reaction-model limit (missing coverage-independent etch channel), not a knob.", flush=True)
print("done", flush=True)
