#!/usr/bin/env python3
"""Decisive test for the radius-based smoothing fix: does ONE fixed disk radius match ViennaPS on BOTH
the TRENCH (which wanted weak smoothing) and the HOLE (which needs smoothing or it stalls)? Compare
flux_smooth_mode='radius' (sweep radius in dx units) against the edge (1-ring) baseline, for both
geometries. ViennaPS trench ref nr=1.0/0.862/0.732 @ AR 3.73/6.1/8.58. Holes: radius should keep the
hole reaching high AR like edge-full-smooth (no stall). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.4, 0.8, 1.3, 1.9]
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def arde(dep):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, np.diff(dep) / np.diff(DURS) / (np.diff(dep) / np.diff(DURS))[0]


def curve(hole, mode, radius, alpha=1.0, nsm=1, seeds=(0, 1, 2)):
    Lxy = 1.5
    GEO = dict(Lx=Lxy, Ly=Lxy if hole else 0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=hole)
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['flux_smooth'] = nsm; p['flux_smooth_alpha'] = alpha
        p['flux_smooth_mode'] = mode; p['flux_smooth_radius'] = radius
        if not hole:
            p['periodic_y'] = 1
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        dfun = t3.max_depth_3d if hole else t3.center_depth_3d
        deps = np.array([dfun(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl,
                         n_ion=40000, n_neu=40000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


curve(False, 'edge', 0, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  radius-smoothing unified-fix test\n", flush=True)
print(f"  ViennaPS TRENCH ref: nr {np.round(VPS_NR,3)} @ AR {VPS_AR}\n", flush=True)
print("  TRENCH (want nr -> ViennaPS 0.862/0.732):", flush=True)
ar, nr = curve(False, 'edge', 0)
print(f"    edge full-smooth : nr@vpsAR {np.round(np.interp(VPS_AR,ar,nr),3)}  gap {np.sqrt(np.mean((np.interp(VPS_AR,ar,nr)-VPS_NR)**2)):.3f}", flush=True)
for R in [1.2, 1.7, 2.5]:
    ar, nr = curve(False, 'radius', R)
    pp = np.interp(VPS_AR, ar, nr)
    print(f"    radius={R:>3} dx     : nr@vpsAR {np.round(pp,3)}  gap {np.sqrt(np.mean((pp-VPS_NR)**2)):.3f}", flush=True)
print("\n  HOLE (must NOT stall; edge full-smooth reaches AR~10):", flush=True)
ar, nr = curve(True, 'edge', 0)
print(f"    edge full-smooth : AR {np.round(ar,1)}  nr {np.round(nr,3)}  maxAR {ar.max():.1f}", flush=True)
for R in [1.2, 1.7, 2.5]:
    ar, nr = curve(True, 'radius', R)
    print(f"    radius={R:>3} dx     : AR {np.round(ar,1)}  nr {np.round(nr,3)}  maxAR {ar.max():.1f}", flush=True)
print("\n  WIN if some radius gives trench gap<~0.02 AND hole maxAR ~ edge-full (no stall).", flush=True)
