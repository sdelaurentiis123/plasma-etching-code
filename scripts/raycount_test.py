#!/usr/bin/env python3
"""GROUNDED test: ViennaRay's smoothFlux formula == petch's exactly (self w=1 + normal-dot-weighted
neighbors, normalized; radius ~1 diskRadius). So the ~0.05 trench gap is NOT the smoothing -- it's that
petch's raw MC flux to the deep floor is more UNDER-SAMPLED, so the (identical) smoothing boosts petch's
starved deep floor more -> gentler. TEST: raise n_neu; if petch steepens toward ViennaPS (0.73) the gap
is sampling. no-smooth@high-rays = the true physical ARDE (smoothing-bias-free). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.4, 0.8, 1.3, 1.9]
GEO = dict(Lx=1.5, Ly=0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def curve(nrays, nsm, seeds=(0, 1)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['flux_smooth'] = nsm
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(
            t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl, n_ion=nrays,
                           n_neu=nrays, reinit_method="fsm", verbose=False, seed_offset=sd * 100, **GEO))
            for dr in DURS])
        accd = deps if accd is None else accd + deps
    dep = accd / len(seeds)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    nr = (np.diff(dep) / np.diff(DURS)); nr = nr / nr[0]
    return armid, nr


curve(40000, 1, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  ViennaPS trench ref nr {np.round(VPS_NR,3)} @ AR {VPS_AR}\n", flush=True)
print(f"  {'n_rays':>7} {'smooth':>7}   {'petch nr @vpsAR':>22}   gapRMSE", flush=True)
for nr_n, nsm, lab in [(40000, 1, 'full'), (160000, 1, 'full'), (400000, 1, 'full'),
                       (40000, 0, 'none'), (400000, 0, 'none')]:
    ar, nr = curve(nr_n, nsm)
    pp = np.interp(VPS_AR, ar, nr)
    gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
    print(f"  {nr_n:7d} {lab:>7}   {np.round(pp,3)!s:>22}   {gap:.3f}", flush=True)
print(f"\n  ViennaPS target          {np.round(VPS_NR,3)}", flush=True)
print("  if full-smooth gap SHRINKS as rays rise -> gap is under-sampling (smoothing over-boosts starved", flush=True)
print("  deep floor). if no-smooth@400k ~ ViennaPS -> ViennaPS = the well-sampled physical ARDE.", flush=True)
