#!/usr/bin/env python3
"""Refine the de Boer match around the cal_F~3 / sig~0.8 optimum (coarse sweep hit RMSE 0.086 at
cal_F=4/sig=1.0). Tune the mid-AR (AR20 came in low: knee drops to floor too fast) by softening the
neutral sticking and lifting the sputter floor a touch. Goal RMSE < 0.06. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 160


def curve(calF, sig_deg, beta, ysp):
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
    return np.interp(EXP_AR, ar, nr), float(np.sqrt(np.mean((np.interp(EXP_AR, ar, nr) - EXP_R) ** 2)))


print(f"device={t3.DEVICE}  de Boer {EXP_R}  (refine cal_F/sig/betaE/Ysp)\n", flush=True)
print(f"  {'cal_F':>5} {'sig':>4} {'beta':>4} {'ysp':>4}   {'@AR 0/10/20/40':>26}   RMSE", flush=True)
best = None
grid = [(3.0, 0.7, 0.80, 8), (3.0, 0.8, 0.78, 12), (3.5, 0.8, 0.80, 10), (3.0, 0.6, 0.82, 14),
        (2.5, 0.7, 0.80, 10), (3.5, 0.7, 0.78, 8), (4.0, 0.8, 0.80, 12), (3.0, 0.9, 0.80, 14)]
for calF, sig, beta, ysp in grid:
    ours, rmse = curve(calF, sig, beta, ysp)
    tag = "  <--" if (best is None or rmse < best[-1]) else ""
    print(f"  {calF:5.1f} {sig:4.1f} {beta:4.2f} {ysp:4.0f}   {np.round(ours,3)!s:>26}   {rmse:.3f}{tag}", flush=True)
    if best is None or rmse < best[-1]:
        best = (calF, sig, beta, ysp, ours, rmse)
print(f"\n  BEST cal_F={best[0]} sig={best[1]} betaE={best[2]} Ysp={best[3]} -> {np.round(best[4],3)}  RMSE={best[5]:.3f}", flush=True)
print(f"  de Boer target                                          {EXP_R}", flush=True)
