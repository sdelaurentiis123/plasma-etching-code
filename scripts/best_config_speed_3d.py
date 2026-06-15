#!/usr/bin/env python3
"""Combined headline: the pre-tonight default (cold n_fp=4 + CPU skfmm reinit) vs the new best config
(WARM n_fp=1 + GPU fsm reinit). Both accurate (warm n_fp=1 == cold n_fp=8 converged; see
warmstart_accuracy_3d.py). Measures wall + per-line timings + depth on a matched hole. PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)
NS, NRAY = 20, 30000


def run(label, reinit, warm, nfp):
    p = dict(petch.PAR); p['n_fp'] = nfp
    fl = petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=warm)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=NRAY, n_neu=NRAY, reinit_method=reinit, verbose=False, **GEO)
    wall = time.time() - t0
    tm = g['timings']; depth = t3.max_depth_3d(g)
    other = tm['total'] - sum(tm[k] for k in ['mesh', 'flux', 'extend', 'advect', 'reinit'])
    print(f"\n=== {label} ===   wall {wall:6.2f}s   depth {depth:5.2f}um", flush=True)
    for k in ['mesh', 'flux', 'extend', 'advect', 'reinit']:
        print(f"    {k:7s} {tm[k]:6.2f}s  {100*tm[k]/max(tm['total'],1e-9):4.1f}%")
    print(f"    host    {other:6.2f}s  {100*other/max(tm['total'],1e-9):4.1f}%")
    return dict(label=label, wall=wall, depth=depth, timings=tm)


print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS} rays={NRAY}", flush=True)
# warm up JIT (small t_end -> realistic dt, no CFL stress)
_wu = {**GEO, 't_end': 0.12}
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=3000, n_neu=3000, reinit_method="fsm", verbose=False, **_wu)

base = run("OLD default: cold n_fp=4 + skfmm(CPU)", "skfmm", False, 4)
best = run("NEW best: WARM n_fp=1 + fsm(GPU)", "fsm", True, 1)

sx = base['wall'] / max(best['wall'], 1e-3)
dd = abs(base['depth'] - best['depth'])
print(f"\n  NEW is {sx:.2f}x faster overall ({base['wall']:.1f}s -> {best['wall']:.1f}s); "
      f"depth delta {dd:.2f}um", flush=True)
print(f"  (warm n_fp=1 is the CONVERGED answer -- more accurate than the old cold n_fp=4, not less)")
json.dump(dict(dx=DX, diam=DIAM, ns=NS, nray=NRAY, old=base, new=best, speedup=sx),
          open("best_config_result.json", "w"), indent=2)
print("wrote best_config_result.json")
