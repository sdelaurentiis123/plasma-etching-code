#!/usr/bin/env python3
"""Verify the FUDGE-FREE petch (cal_F removed) matches ViennaPS. Compares the normalized trench ARDE to
the CONVERGED ViennaPS reference measured earlier (dx=0.04, raysPerPoint converged): depths
1.297/2.437/3.666/4.918 -> nr 1.0/0.862/0.732 @ AR 3.73/6.10/8.58. Also runs the hole to see where the
un-fudged hole ARDE lands. Normalized ARDE (rate_scale cancels), seed-averaged. PETCH_DEVICE=cuda."""
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
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    r = np.diff(dep) / np.diff(DURS)
    return armid, r / r[0]


def curve(hole, seeds=(0, 1, 2)):
    Lxy = 1.5
    GEO = dict(Lx=Lxy, Ly=Lxy if hole else 0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=hole)
    accd = None; topdepth = []
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1
        if not hole:
            p['periodic_y'] = 1
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        dfun = t3.max_depth_3d if hole else t3.center_depth_3d
        deps = np.array([dfun(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl,
                         n_ion=40000, n_neu=40000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    dep = accd / len(seeds)
    return dep, arde(dep)


curve(False, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  FUDGE-FREE petch (cal_F removed) vs ViennaPS\n", flush=True)
print(f"  ViennaPS trench ref (converged): nr {np.round(VPS_NR,3)} @ AR {VPS_AR}\n", flush=True)

dep, (ar, nr) = curve(False)
pp = np.interp(VPS_AR, ar, nr)
gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
print(f"  TRENCH fudge-free: depths {np.round(dep,3)}", flush=True)
print(f"    AR {np.round(ar,2)}  nr {np.round(nr,3)}", flush=True)
print(f"    nr @ViennaPS AR = {np.round(pp,3)}  vs {np.round(VPS_NR,3)}  -> gapRMSE {gap:.3f}", flush=True)
print(f"    {'MATCH (<0.03)' if gap < 0.03 else ('CLOSE (<0.05)' if gap < 0.05 else 'still off')}\n", flush=True)

deph, (arh, nrh) = curve(True)
print(f"  HOLE fudge-free: depths {np.round(deph,3)}", flush=True)
print(f"    AR {np.round(arh,2)}  nr {np.round(nrh,3)}  (honest un-fudged hole ARDE)", flush=True)
print(f"\n  trench gap to ViennaPS = {gap:.3f}. rate_scale only sets absolute depth (cancels in ARDE).", flush=True)
