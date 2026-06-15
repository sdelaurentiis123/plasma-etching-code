#!/usr/bin/env python3
"""Minimal petch example: etch a 3D hole and print the profile.

On an NVIDIA GPU this runs the full GPU pipeline (~14x faster than ViennaPS) automatically.
On CPU / Apple Silicon the SAME code runs the portable numpy/skimage path (slower).

    PETCH_DEVICE=cuda python examples/etch_hole.py     # GPU
    PETCH_DEVICE=cpu  python examples/etch_hole.py     # CPU (or just omit; defaults to cpu)
"""
import time
import petch
from petch import threed as t3

# The fast + accurate config is ONE flag set. warm_start_coverage makes the coverage fixed point
# converge in 1 iteration; on CUDA the GPU speedups (marching cubes, on-device flux, warm-start,
# source generation) all auto-enable. On CPU they auto-fall-back. No other tuning needed.
flags = petch.Flags(coverage_sticking=True, warm_start_coverage=True, sampling="sobol")

t0 = time.time()
geo = t3.run_etch_3d(
    trench_width=6.0, dx=0.25, n_steps=40,
    Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0,
    par=dict(petch.PAR, rate_scale=0.07),
    flags=flags, n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False,
)
wall = time.time() - t0

print(f"device      : {t3.DEVICE}")
print(f"wall clock  : {wall:.2f} s   (ViennaPS-GPU on the same hole: ~22 s)")
print(f"center depth: {t3.center_depth_3d(geo):.2f} µm")
print(f"max depth   : {t3.max_depth_3d(geo):.2f} µm")
print(f"mesh faces  : {len(t3.extract_mesh_3d(geo['phi'], 0.25)[1])}")
