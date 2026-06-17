#!/usr/bin/env python3
"""Two-component fit to de Boer: the experiment's ARDE = a steep NEUTRAL (chemical) knee + a flat
ION-SPUTTER floor. petch has both terms (belen.py) but the ViennaPS-default physical-sputter floor is
~7% of the mouth rate while de Boer's high-AR tail is ~20%. Sweep (Ysp_scale = sputter floor lift)
x (betaE = neutral sticking = knee position) on MC transport and find the best match to
1.0/0.43/0.29/0.20. If RMSE drops below ~0.06 -> the two-component model REPRODUCES de Boer.
PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 160


def curve(beta, ysp):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta; p['Ysp_scale'] = ysp; p['periodic_y'] = 1
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


print(f"device={t3.DEVICE}  de Boer {EXP_R} @ AR {EXP_AR}  (2-knob: betaE x Ysp_scale, MC)\n", flush=True)
print(f"  {'betaE':>6} {'Ysp':>5}   {'@AR 0/10/20/40':>26}   RMSE", flush=True)
best = None
for beta in [0.80, 0.88, 0.94]:
    for ysp in [1.0, 4.0, 8.0, 14.0]:
        ours, rmse = curve(beta, ysp)
        tag = "  <--" if (best is None or rmse < best[3]) else ""
        print(f"  {beta:6.2f} {ysp:5.1f}   {np.round(ours,3)!s:>26}   {rmse:.3f}{tag}", flush=True)
        if best is None or rmse < best[3]:
            best = (beta, ysp, ours, rmse)
print(f"\n  BEST betaE={best[0]} Ysp_scale={best[1]} -> {np.round(best[2],3)}  RMSE={best[3]:.3f}", flush=True)
print(f"  de Boer target                                {EXP_R}", flush=True)
