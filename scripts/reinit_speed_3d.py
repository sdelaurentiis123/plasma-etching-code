#!/usr/bin/env python3
"""WAVE 2 SPEED: FSM (GPU Jacobi Godunov-Eikonal) vs skfmm-narrow (CPU) reinit, on the full 3D loop.

Correctness is already proven (reinit_correctness_3d.py, run on CPU: |grad|~1, depth delta 0.000).
This measures the wall-clock win on a GPU: FSM keeps reinit on-device (no per-step CPU round-trip =
the ~40% bottleneck). Reports total wall, the reinit line of the timings dict, and depth parity.
Run on a GPU box: PETCH_DEVICE=cuda python scripts/reinit_speed_3d.py
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
NS, NRAY = 30, 30000


def run(meth):
    p = dict(petch.PAR)
    # warm-up compile not timed: one tiny step
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NRAY, n_neu=NRAY, reinit_method=meth, verbose=False, **GEO)
    wall = time.time() - t0
    tm = g['timings']
    return wall, tm, t3.max_depth_3d(g)


print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS} rays={NRAY}\n", flush=True)
# warm up the JIT (compile both reinit kernels) with a 2-step throwaway
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                   n_ion=2000, n_neu=2000, reinit_method="fsm", verbose=False, **GEO)
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                   n_ion=2000, n_neu=2000, reinit_method="skfmm", verbose=False, **GEO)

res = {}
for meth in ["skfmm", "fsm"]:
    wall, tm, depth = run(meth)
    res[meth] = dict(wall=wall, reinit=tm['reinit'], total=tm['total'], depth=depth,
                     pct=100 * tm['reinit'] / max(tm['total'], 1e-9))
    print(f"  {meth:6s}  wall {wall:6.2f}s   reinit {tm['reinit']:6.2f}s ({res[meth]['pct']:4.1f}%)   "
          f"depth {depth:5.2f}um", flush=True)

sx = res['skfmm']['wall'] / max(res['fsm']['wall'], 1e-3)
rx = res['skfmm']['reinit'] / max(res['fsm']['reinit'], 1e-3)
dd = abs(res['skfmm']['depth'] - res['fsm']['depth'])
print(f"\n  FSM reinit: {rx:.1f}x faster on the reinit line, {sx:.2f}x faster overall; "
      f"depth delta {dd:.2f}um ({'accuracy-neutral' if dd < 0.5 else 'CHECK'})", flush=True)
json.dump(dict(dx=DX, diam=DIAM, ns=NS, nray=NRAY, result=res, reinit_speedup=rx,
               total_speedup=sx, depth_delta=dd), open("reinit_speed_result.json", "w"), indent=2)
print("wrote reinit_speed_result.json")
