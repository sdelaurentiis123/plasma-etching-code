#!/usr/bin/env python3
"""FIX petch's over-gentle ARDE: sweep flux_smooth_alpha (lateral floor-feeding strength) to steepen
petch onto the ViennaPS reference. ViennaPS points measured this run (dx=0.04, same trench), seed-
averaged petch. Hypothesis: petch's flux smoothing over-feeds the HARC floor -> gentler ARDE; lower
alpha -> steeper -> matches ViennaPS. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
SUB = 7.0
DURS = [0.4, 0.8, 1.3, 1.9]                     # match the ViennaPS reference points
GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)

# ViennaPS reference (this box, CPU_TRIANGLE, dx=0.04)
VPS_DEP = np.array([1.297, 2.437, 3.666, 4.918])


def arde(dep, dur):
    dep = np.asarray(dep, float); dur = np.asarray(dur, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dur)
    return armid, rate / rate[0]


vps_ar, vps_nr = arde(VPS_DEP, DURS)


def petch_curve(alpha, n_smooth, seeds=(0, 1, 2)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
        p['flux_smooth'] = n_smooth; p['flux_smooth_alpha'] = alpha
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(
            t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl, n_ion=40000,
                           n_neu=40000, reinit_method="fsm", verbose=False, seed_offset=sd * 100, **GEO))
            for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds), DURS)


petch_curve(1.0, 1, seeds=(0,))  # warm
print(f"device={t3.DEVICE}", flush=True)
print(f"  ViennaPS ref: AR {np.round(vps_ar,2)}  nr {np.round(vps_nr,3)}\n", flush=True)
print(f"  {'alpha':>5} {'nsm':>3}   {'petch nr @ViennaPS AR':>24}   gapRMSE", flush=True)
best = None
for alpha, nsm in [(1.0, 1), (0.6, 1), (0.35, 1), (0.15, 1), (0.0, 0)]:
    ar, nr = petch_curve(alpha, nsm)
    pp = np.interp(vps_ar, ar, nr)
    gap = float(np.sqrt(np.mean((pp - vps_nr) ** 2)))
    tag = "  <--" if (best is None or gap < best[-1]) else ""
    print(f"  {alpha:5.2f} {nsm:3d}   {np.round(pp,3)!s:>24}   {gap:.3f}{tag}", flush=True)
    if best is None or gap < best[-1]:
        best = (alpha, nsm, pp, gap)
print(f"\n  BEST alpha={best[0]} nsm={best[1]} -> {np.round(best[2],3)}  gapRMSE={best[3]:.3f}", flush=True)
print(f"  ViennaPS target                 {np.round(vps_nr,3)}", flush=True)
print(f"  (default alpha=1.0 row = the over-gentle bias; lower alpha should shrink gapRMSE)", flush=True)
