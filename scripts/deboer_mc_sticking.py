#!/usr/bin/env python3
"""DECISIVE TEST for the de Boer gap. petch's MC neutral transport (_trace3d_cov_rr) is already
free-molecular Knudsen MC (sticking S=bare*beta at every wall hit, diffuse cosine re-emission, RR,
no bounce cap). So the question is NOT "ballistic vs Knudsen" -- it's whether the STICKING coefficient
is the missing calibration. The old betaE sweep ran on RADIOSITY (which over-couples and barely
responds to beta at depth). This sweeps betaE on MC transport.

If some betaE lands MC on de Boer 1.0/0.43/0.29/0.20 -> we match with calibrated sticking, no new
module. If MC saturates too gentle at every betaE -> genuinely missing physics (report which).
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


def curve(beta):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta; p['periodic_y'] = 1
    fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol",
                     warm_start_coverage=True)
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=30000, n_neu=30000,
                       reinit_method="fsm", verbose=False, record_depth_every=5, **GEO)
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
    ours = np.interp(EXP_AR, ar, nr)
    return ours, float(np.sqrt(np.mean((ours - EXP_R) ** 2))), ar.max()


print(f"device={t3.DEVICE}  de Boer {EXP_R} @ AR {EXP_AR}  (MC transport + periodic-y)\n", flush=True)
print(f"  {'betaE':>6}   {'@AR 0/10/20/40':>26}   {'RMSE':>6}  AR_max", flush=True)
best = None
for b in [0.70, 0.85, 0.93, 0.97, 0.99]:
    ours, rmse, armax = curve(b)
    print(f"  {b:6.2f}   {np.round(ours,3)!s:>26}   {rmse:6.3f}  {armax:.1f}", flush=True)
    if best is None or rmse < best[2]:
        best = (b, ours, rmse)
print(f"\n  BEST betaE={best[0]} RMSE={best[2]:.3f}", flush=True)
print(f"  if RMSE<~0.08 at some betaE -> MC+sticking REPRODUCES de Boer (match achieved, calibrate sticking).", flush=True)
print(f"  if RMSE stays >~0.15 everywhere -> MC saturates too gentle; missing physics, not transport.", flush=True)
