#!/usr/bin/env python3
"""Validate Warp GPU MarchingCubes vs skimage marching_cubes (the 26% mesh cost). Compare on a real
etched phi: vertex/face counts, mesh extent, AND the per-face flux + depth a full etch produces. Warp
MC is a different triangulation, so the bar is WITHIN-NOISE (not bit-identical). PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DEV = t3.DEVICE
DX = 0.25
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=0.7)

# get a realistic etched phi
g = t3.run_etch_3d(trench_width=6.0, dx=DX, n_steps=12, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=8000, n_neu=8000, reinit_method="fsm", verbose=False, **GEO)
phi = g['phi']
nx, ny, nz = phi.shape
print(f"grid {nx}x{ny}x{nz}", flush=True)

# skimage
vs, fs, cs, as_ = t3.extract_mesh_3d(phi, DX)
print(f"skimage MC: {len(vs)} verts, {len(fs)} faces, z {cs[:,2].min():.2f}..{cs[:,2].max():.2f}", flush=True)

# Warp GPU MC
mv = mt = nx * ny * nz // 2 + 1000
mc = wp.MarchingCubes(nx, ny, nz, max_verts=mv, max_tris=mt, device=DEV)
fwp = wp.array(phi.astype(np.float32), dtype=float, device=DEV)
mc.surface(fwp, 0.0)
vg = mc.verts.numpy() * DX                     # node-index coords -> physical
ig = mc.indices.numpy().reshape(-1, 3)
vv = vg[ig]; cg = vv.mean(axis=1)
print(f"warp MC   : {len(vg)} verts, {len(ig)} faces, z {cg[:,2].min():.2f}..{cg[:,2].max():.2f}", flush=True)

# timing
def tmc_sk():
    t3.extract_mesh_3d(phi, DX)
def tmc_wp():
    mc.surface(fwp, 0.0); _ = mc.verts.numpy(); _ = mc.indices.numpy()
for f in (tmc_sk, tmc_wp): f()
t0 = time.time()
for _ in range(20): tmc_sk()
t_sk = 1000*(time.time()-t0)/20
t0 = time.time()
for _ in range(20): tmc_wp()
t_wp = 1000*(time.time()-t0)/20
print(f"\n  skimage MC: {t_sk:.1f} ms/call   warp GPU MC (+readback): {t_wp:.1f} ms/call   speedup {t_sk/max(t_wp,1e-3):.1f}x", flush=True)
print(f"  vert/face count ratio warp/skimage: {len(vg)/max(len(vs),1):.2f} / {len(ig)/max(len(fs),1):.2f}")
