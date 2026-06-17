#!/usr/bin/env python3
"""Tighten the KNEE match to de Boer (ignore the high-AR floor frontier for now). The knee is still
right-shifted (petch hits 0.43 at ~AR13 vs experiment AR9; petch too high at AR5). Push more etchant-
starved (lower cal_F) + faster near-top depletion (higher betaE). SEED-AVERAGED inside the sweep (2
seeds) so we optimize signal not noise. Metric = RMSE over the KNEE band AR {5,8,10,14,18} only.
PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

# experiment, knee band only (exclude the AR>20 floor frontier)
KNEE_AR = np.array([5.0, 8.0, 10.0, 14.0, 18.0])
KNEE_R = np.interp(KNEE_AR, [0, 10, 20, 40], [1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 160


def curve(calF, beta, sig=0.8, ysp=10.0, seeds=(0, 1)):
    acc = np.zeros(len(KNEE_AR)); cnt = 0
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1
        p['cal_F'] = calF; p['betaE'] = beta; p['ion_ang_sigma'] = np.deg2rad(sig); p['Ysp_scale'] = ysp
        fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol",
                         warm_start_coverage=True)
        g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                           reinit_method="fsm", verbose=False, record_depth_every=4, seed_offset=sd * 100, **GEO)
        h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
        tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
        r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
        nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
        acc += np.interp(KNEE_AR, ar, nr); cnt += 1
    ours = acc / cnt
    return ours, float(np.sqrt(np.mean((ours - KNEE_R) ** 2)))


print(f"device={t3.DEVICE}  de Boer KNEE band {np.round(KNEE_R,3)} @ AR {KNEE_AR}  (2-seed avg)\n", flush=True)
print(f"  {'cal_F':>5} {'beta':>4}   {'@AR 5/8/10/14/18':>30}   kneeRMSE", flush=True)
best = None
for calF in [1.5, 2.0, 2.5, 3.0]:
    for beta in [0.85, 0.93]:
        ours, rmse = curve(calF, beta)
        tag = "  <--" if (best is None or rmse < best[-1]) else ""
        print(f"  {calF:5.1f} {beta:4.2f}   {np.round(ours,3)!s:>30}   {rmse:.3f}{tag}", flush=True)
        if best is None or rmse < best[-1]:
            best = (calF, beta, ours, rmse)
print(f"\n  BEST cal_F={best[0]} betaE={best[1]} -> {np.round(best[2],3)}  kneeRMSE={best[3]:.3f}", flush=True)
print(f"  de Boer knee band                       {np.round(KNEE_R,3)}", flush=True)
