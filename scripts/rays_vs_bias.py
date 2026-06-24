#!/usr/bin/env python3
"""Is the fudge-free trench steepness UNDER-SAMPLING or a real transport BIAS? Run fudge-free petch
(Fflux=1800, cal_F gone) at rising ray counts. If the deep-floor nr climbs toward ViennaPS (0.862/0.732)
as rays rise -> it was under-sampling (fix = deep-floor importance sampling / more rays). If it plateaus
too steep -> real transport bias (keep hunting). Extended durations (fudge-free etches slow). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.6, 1.4, 2.6, 4.2]
GEO = dict(Lx=1.5, Ly=0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def arde(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    r = np.diff(dep) / np.diff(DURS); return armid, r / r[0]


def curve(nrays, seeds=(0, 1)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1   # Fflux=1800 default (fudge-free)
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=nrays, n_neu=nrays, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


curve(40000, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  fudge-free (Fflux=1800) trench: under-sampling or bias?\n", flush=True)
print(f"  ViennaPS ref: nr 0.862/0.732 @ AR 6.1/8.58\n", flush=True)
print(f"  {'n_rays':>8}   {'nr @vpsAR':>22}   gapRMSE", flush=True)
for nr_n in [40000, 120000, 300000]:
    ar, nr = curve(nr_n)
    pp = np.interp(VPS_AR, ar, nr)
    gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
    print(f"  {nr_n:8d}   {np.round(pp,3)!s:>22}   {gap:.3f}", flush=True)
print(f"\n  ViennaPS target           {np.round(VPS_NR,3)}", flush=True)
print("  climbs toward ViennaPS w/ rays -> under-sampling. plateaus steep -> real transport bias.", flush=True)
