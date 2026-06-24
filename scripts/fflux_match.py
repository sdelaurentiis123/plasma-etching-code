#!/usr/bin/env python3
"""cal_F removal overshot: fudge-free petch (Fflux=1800) is TOO STEEP (nr 0.44 vs ViennaPS 0.73). cal_F
just scaled the etchant flux, so sweep the physical Fflux (= 1800 * factor) to find what lands petch on
the converged ViennaPS trench ref (nr 1.0/0.862/0.732 @ AR 3.73/6.10/8.58). Tells us the real
normalization gap: if petch needs Fflux >> 1800 to match, petch under-delivers neutral flux to the deep
floor vs ViennaPS by that factor. Extended durations so even low Fflux spans the AR range. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.5, 1.1, 1.9, 3.0]          # extended so low-Fflux still reaches AR~8-9
GEO = dict(Lx=1.5, Ly=0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def arde(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    r = np.diff(dep) / np.diff(DURS); return armid, r / r[0]


def curve(fflux, seeds=(0, 1, 2)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['Fflux'] = fflux
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=40000, n_neu=40000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


curve(1800.0, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  Fflux sweep to match ViennaPS trench (ref nr 1.0/0.862/0.732)\n", flush=True)
print(f"  {'Fflux':>7} {'(xViennaPS)':>11}   {'petch nr @vpsAR':>22}   gapRMSE", flush=True)
best = None
for k in [1, 2, 4, 6, 9, 12]:
    ar, nr = curve(1800.0 * k)
    pp = np.interp(VPS_AR, ar, nr)
    gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
    tag = "  <--" if (best is None or gap < best[-1]) else ""
    print(f"  {1800*k:7.0f} {('x'+str(k)):>11}   {np.round(pp,3)!s:>22}   {gap:.3f}{tag}", flush=True)
    if best is None or gap < best[-1]:
        best = (k, pp, gap)
print(f"\n  BEST Fflux=x{best[0]} ({1800*best[0]:.0f}) -> {np.round(best[1],3)}  gapRMSE={best[2]:.3f}", flush=True)
print(f"  ViennaPS uses Fflux=1800 (x1). If petch needs x{best[0]} to match -> petch under-delivers deep", flush=True)
print(f"  neutral flux ~{best[0]}x vs ViennaPS = the real normalization/transport gap cal_F was masking.", flush=True)
