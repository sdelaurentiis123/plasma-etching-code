#!/usr/bin/env python3
"""Calibration vs physics: sweep wall sticking (betaE) with periodic-y and see how close the trench ARDE
gets to de Boer. If a betaE snaps onto 1.0/0.43/0.29/0.20 it's mostly calibration; if none does, the
residual is the ballistic-vs-Knudsen transport shape. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0; GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 140


def curve(beta, nt="radiosity"):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta; p['periodic_y'] = 1
    fl = petch.Flags(coverage_sticking=True, neutral_transport=nt, sampling="sobol", warm_start_coverage=True)
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=30000, n_neu=30000,
                       reinit_method="fsm", verbose=False, record_depth_every=6, **GEO)
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    t = st / NS * GEO['t_end']; r = np.gradient(dd, t); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
    ours = np.interp(EXP_AR, ar, nr)
    return ours, float(np.sqrt(np.mean((ours - EXP_R) ** 2)))


print(f"device={t3.DEVICE}  de Boer {EXP_R} @ AR {EXP_AR}  (radiosity + periodic-y)\n", flush=True)
print(f"  {'betaE':>6}   {'@AR 0/10/20/40':>26}   RMSE", flush=True)
best = None
for b in [0.47, 0.65, 0.8, 0.9, 0.97]:
    ours, rmse = curve(b)
    print(f"  {b:6.2f}   {np.round(ours,3)!s:>26}   {rmse:.3f}", flush=True)
    if best is None or rmse < best[2]:
        best = (b, ours, rmse)
print(f"\n  BEST betaE={best[0]} RMSE={best[2]:.3f}.  If RMSE stays >~0.15 at every betaE -> not calibration,", flush=True)
print(f"  it's the ballistic-vs-Knudsen transport shape (needs the conductance/molecular-flow model).", flush=True)
