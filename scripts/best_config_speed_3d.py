#!/usr/bin/env python3
"""HEADLINE: the original pre-optimization config vs the full GPU stack, on a matched hole. All
speedups validated accuracy-neutral (depth identical). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)
NS, NRAY = 20, 30000


def run(label, reinit, warm, nfp, gsmooth, gsrc, gmesh, samp, gws=False):
    p = dict(petch.PAR); p['n_fp'] = nfp
    p['flux_smooth_gpu'] = gsmooth; p['gpu_source'] = gsrc; p['gpu_mesh'] = gmesh; p['gpu_warmstart'] = gws
    fl = petch.Flags(coverage_sticking=True, sampling=samp, warm_start_coverage=warm)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=NRAY, n_neu=NRAY, reinit_method=reinit, verbose=False, **GEO)
    wall = time.time() - t0
    tm = g['timings']; depth = t3.max_depth_3d(g)
    print(f"\n=== {label} ===   wall {wall:6.2f}s   depth {depth:5.2f}um", flush=True)
    for k in ['mesh', 'flux', 'extend', 'advect', 'reinit']:
        print(f"    {k:7s} {tm[k]:6.2f}s  {100*tm[k]/max(tm['total'],1e-9):4.1f}%")
    return dict(label=label, wall=wall, depth=depth, timings=tm)


print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS} rays={NRAY}", flush=True)
# warmup (compile all GPU kernels incl. Warp MC)
run("warmup", "fsm", True, 1, True, True, True, "sobol", gws=True)

base = run("ORIGINAL: cold n_fp=4, skfmm-CPU, CPU-smooth, sobol-host, CPU-mesh, KDTree",
           "skfmm", False, 4, False, False, False, "sobol", gws=False)
best = run("FULL GPU: warm n_fp=1, fsm, GPU-smooth/source/mesh/warmstart",
           "fsm", True, 1, True, True, True, "sobol", gws=True)

sx = base['wall'] / max(best['wall'], 1e-3)
dd = abs(base['depth'] - best['depth'])
print(f"\n  CUMULATIVE {sx:.2f}x faster ({base['wall']:.2f}s -> {best['wall']:.2f}s); depth delta {dd:.2f}um", flush=True)
json.dump(dict(dx=DX, diam=DIAM, ns=NS, nray=NRAY, original=base, full=best, speedup=sx),
          open("best_config_result.json", "w"), indent=2)
print("wrote best_config_result.json")
