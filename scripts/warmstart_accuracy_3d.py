#!/usr/bin/env python3
"""WAVE 2c GATE: is WARM-START coverage accuracy-neutral? The coverage fixed point converges to the
SAME point whether seeded from bare=ones (cold, n_fp=4) or from the previous step's coverage (warm,
n_fp=1-2) -- the front moves <1 cell/step so coverage barely moves. If the depth TRAJECTORY matches
cold n_fp=4, warm-start is a free 2-4x flux win (fewer neutral launches) with ZERO physics change.

Needs the CLEAN regime (30k rays, dx=0.25) -- at low ray count MC noise (~1 cell) swamps the warm-vs-
cold signal. Run on a box: PETCH_DEVICE=cuda python this. Truth = cold n_fp=6 (well converged).
Geometry is tuned to NOT bottom out (final depth < sub_top) so the trajectory is informative.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)   # deep substrate: do NOT bottom out
NS, NRAY = 20, 30000


def traj(label, warm, nfp):
    p = dict(petch.PAR); p['n_fp'] = nfp
    fl = petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=warm)
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=NRAY, n_neu=NRAY, reinit_method="fsm", verbose=False,
                       record_depth_every=2, **GEO)
    h = g['depth_history']
    st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    launches = 2 * nfp * NS    # neutral MC launches over the run
    print(f"  {label:22s} depth_final={dd[-1]:5.2f}um  neutral_launches={launches}", flush=True)
    return st, dd


print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS} rays={NRAY}\n", flush=True)
# warm up JIT (fsm + coupled kernels) so timings/first-step aren't penalized
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                   n_ion=3000, n_neu=3000, reinit_method="fsm", verbose=False, **GEO)
st0, base = traj("COLD n_fp=8 (truth)", False, 8)
runs = {
    "COLD n_fp=6": traj("COLD n_fp=6", False, 6),
    "COLD n_fp=4": traj("COLD n_fp=4", False, 4),
    "WARM n_fp=2": traj("WARM n_fp=2", True, 2),
    "WARM n_fp=1": traj("WARM n_fp=1", True, 1),
    "COLD n_fp=2 (ctrl)": traj("COLD n_fp=2 (ctrl)", False, 2),
}
print("\n  max |depth - cold_nfp4| over the trajectory:", flush=True)
for k, (st, dd) in runs.items():
    d = np.interp(st0, st, dd)
    err = np.abs(d - base)
    tag = "PASS (accuracy-neutral)" if err.mean() < 0.30 else "drifts"
    nfp_k = int(k.split("n_fp=")[1][0])
    print(f"    {k:22s} max {err.max():.3f}um  mean {err.mean():.3f}um  ({8.0/nfp_k:.1f}x fewer launches than truth) -> {tag}", flush=True)
print("\n  If WARM n_fp=1/2 PASS but COLD n_fp=2 drifts -> warm-start gives full accuracy at 2-4x less flux work.")
