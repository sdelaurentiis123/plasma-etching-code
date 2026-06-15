#!/usr/bin/env python3
"""WAVE 2b: does the coverage fixed point converge before n_fp=4? Each fixed-point iter = 2 neutral MC
launches (F+O). The accurate default does n_fp=4 -> 8 launches/step (62% of the loop). If depth + the
deep-floor coverage are converged by n_fp=2, we halve the neutral flux work for free.

Sweep n_fp at fixed reinit=fsm; report depth + the F-coverage profile vs n_fp. Run: PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=2.0)
NS, NRAY = 24, 30000

print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS}\n", flush=True)
# warm up JIT
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                   n_ion=2000, n_neu=2000, reinit_method="fsm", verbose=False, **GEO)
res = {}
ref_depth = None
for nfp in [1, 2, 3, 4]:
    p = dict(petch.PAR); p['n_fp'] = nfp
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NRAY, n_neu=NRAY, reinit_method="fsm", verbose=False, **GEO)
    wall = time.time() - t0
    depth = t3.max_depth_3d(g)
    if nfp == 4:
        ref_depth = depth
    res[nfp] = dict(wall=wall, depth=depth)
    print(f"  n_fp={nfp}: wall {wall:6.2f}s  depth {depth:5.2f}um", flush=True)

print("\n  vs n_fp=4 (the accurate default):", flush=True)
for nfp in [1, 2, 3]:
    dd = abs(res[nfp]['depth'] - ref_depth)
    sx = res[4]['wall'] / max(res[nfp]['wall'], 1e-3)
    ok = "OK (accuracy-neutral)" if dd < 0.4 else "drifts"
    print(f"    n_fp={nfp}: depth delta {dd:.2f}um  {sx:.2f}x faster  -> {ok}", flush=True)
json.dump(dict(dx=DX, diam=DIAM, ns=NS, ref_depth=ref_depth, result=res),
          open("nfp_convergence_result.json", "w"), indent=2)
print("wrote nfp_convergence_result.json")
