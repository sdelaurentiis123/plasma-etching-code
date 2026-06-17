#!/usr/bin/env python3
"""Match de Boer with de-Boer PROCESS PARAMETERS (not ViennaPS defaults). The two gaps each have a
physical cause and a physical knob:
  - KNEE too deep (AR10=0.88 vs 0.43): petch defaults are etchant-RICH (cal_F=12 -> 1800x excess ->
    coverage saturated until ~98% neutral depletion -> knee at AR15-20, like ViennaPS). de Boer cryo
    DRIE is etchant-STARVED/ion-driven -> early knee. Knob: lower cal_F (etchant excess).
  - FLOOR too low (AR40=0.07 vs 0.20): ion IADF=2.5deg -> ions hit a sidewall by AR~11 -> no floor.
    de Boer needs a sub-degree IADF core. Knob: ion_ang_sigma.
Sweep cal_F (knee) x ion_ang_sigma (floor) with a modest sputter floor; find the best de Boer match.
RMSE<~0.06 -> petch REPRODUCES de Boer with physical process params. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 160


def curve(calF, sig_deg, beta=0.9, ysp=4.0):
    p = dict(petch.PAR)
    p['rate_scale'] = 0.30; p['betaE'] = beta; p['Ysp_scale'] = ysp; p['periodic_y'] = 1
    p['cal_F'] = calF; p['ion_ang_sigma'] = np.deg2rad(sig_deg)
    fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol",
                     warm_start_coverage=True)
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=30000, n_neu=30000,
                       reinit_method="fsm", verbose=False, record_depth_every=5, **GEO)
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
    ours = np.interp(EXP_AR, ar, nr)
    return ours, float(np.sqrt(np.mean((ours - EXP_R) ** 2)))


print(f"device={t3.DEVICE}  de Boer {EXP_R} @ AR {EXP_AR}  (process params: cal_F x ion_ang_sigma)\n", flush=True)
print(f"  {'cal_F':>6} {'sig':>5}   {'@AR 0/10/20/40':>26}   RMSE", flush=True)
best = None
for calF in [1.0, 2.0, 4.0]:
    for sig in [0.5, 1.0]:
        ours, rmse = curve(calF, sig)
        tag = "  <--" if (best is None or rmse < best[3]) else ""
        print(f"  {calF:6.1f} {sig:5.1f}   {np.round(ours,3)!s:>26}   {rmse:.3f}{tag}", flush=True)
        if best is None or rmse < best[3]:
            best = (calF, sig, ours, rmse)
print(f"\n  BEST cal_F={best[0]} sig={best[1]}deg -> {np.round(best[2],3)}  RMSE={best[3]:.3f}", flush=True)
print(f"  de Boer target                              {EXP_R}", flush=True)
