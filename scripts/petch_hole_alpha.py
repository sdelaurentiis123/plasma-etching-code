#!/usr/bin/env python3
"""Does lowering flux_smooth_alpha (which fixed the TRENCH gap) REGRESS holes? Run a petch HOLE ARDE
(depth vs AR over durations) at alpha=1.0 (full smooth) vs alpha=0.0 (none). If the hole ARDE gets
STEEPER at alpha=0 (lower normalized rate at depth), holes prefer full smoothing -> the trench and hole
smoothing optima CONFLICT -> a single global alpha can't serve both (need radius-based smoothing).
Same geometry/resolution as the trench test (dx=0.04). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W = 0.04, 0.5
SUB = 7.0
DURS = [0.4, 0.8, 1.3, 1.9]
GEO = dict(Lx=1.5, Ly=1.5, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=True)


def arde(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(DURS)
    return armid, rate / rate[0]


def hole_curve(alpha, nsm, seeds=(0, 1, 2)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['flux_smooth'] = nsm; p['flux_smooth_alpha'] = alpha
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.max_depth_3d(
            t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl, n_ion=40000,
                           n_neu=40000, reinit_method="fsm", verbose=False, seed_offset=sd * 100, **GEO))
            for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


hole_curve(1.0, 1, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  petch HOLE ARDE: alpha sensitivity (does low alpha regress holes?)\n", flush=True)
for alpha, nsm, lab in [(1.0, 1, "full smooth"), (0.0, 0, "no smooth")]:
    ar, nr = hole_curve(alpha, nsm)
    print(f"  alpha={alpha} ({lab:11s}): AR {np.round(ar,2)}  nr {np.round(nr,3)}", flush=True)
print("\n  if no-smooth nr is LOWER at depth -> holes get steeper without smoothing -> conflict with the", flush=True)
print("  trench (which wanted low alpha). One global alpha can't serve both -> radius-based smoothing.", flush=True)
