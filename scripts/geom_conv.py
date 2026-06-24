#!/usr/bin/env python3
"""Does the smooth-SDF geometry make the trench width grid-INDEPENDENT? Build make_trench_3d at
dx=0.25/0.15/0.10, extract the mesh, measure the mask-wall x-position (true = W/2). Binary carve
quantizes it per-dx (the dx-drift source); smooth SDF should land on W/2 at every dx."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
from petch import threed as t3

W, XE, YE = 0.5, 1.5, 1.0
MASKTOP = 0.3 + 0.16            # sub_top + mask_th region to sample the wall


def wall_x(dx):
    sub_top = 3.3; mask_th = 0.16
    geo = t3.make_trench_3d(XE, YE, sub_top + mask_th + 0.3, dx, W, mask_th, sub_top, hole=False)
    verts, faces, cen, areas = t3.extract_mesh_3d(geo['phi'], dx)
    v = verts[faces]; nrm = np.cross(v[:, 1]-v[:, 0], v[:, 2]-v[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
    cx = cen[:, 0] - XE/2; cz = cen[:, 2]
    # mask-wall facets: near-vertical (|nrm_x|>0.5), in the mask z-band, |x| near W/2
    wall = (np.abs(nrm[:, 0]) > 0.5) & (cz > sub_top) & (cz < sub_top + mask_th) & (np.abs(np.abs(cx) - W/2) < 0.2)
    if wall.sum() == 0:
        return np.nan, 0
    return float(np.abs(cx[wall]).mean()), int(wall.sum())


print(f"device={t3.DEVICE}   true wall x = W/2 = {W/2:.4f}\n", flush=True)
print("dx     wall_x   err(cells)   nfacets", flush=True)
for dx in [0.25, 0.15, 0.10]:
    wx, n = wall_x(dx)
    print(f"{dx:.2f}   {wx:.4f}   {abs(wx-W/2)/dx:.3f}        {n}", flush=True)
print("\n  grid-INDEPENDENT => wall_x ~ 0.2500 at every dx (was: binary carve quantizes it per-dx).", flush=True)
