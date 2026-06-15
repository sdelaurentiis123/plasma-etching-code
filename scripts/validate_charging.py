#!/usr/bin/env python3
"""Validate the surface-charging model vs Hwang-Giapis 1997: the floor ion CURRENT drops with aspect
ratio (~60% by AR~4) because electrons (diffuse) are more shadowed than ions (directional) -> the floor
floats positive -> ion throttling. We test the RAW physics m_e/m_i and the implied charging factor on a
single DEEP hole: every face has a local AR = (sub_top - z)/width, so one hole gives the whole curve.
Run on a box: PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DEV = t3.DEVICE
DX, W = 0.25, 4.0
GEO = dict(Lx=12, Ly=12, Lz=40, mask_th=2, sub_top=34, hole=True, t_end=2.4)

# etch a DEEP hole (charging off) to get a high-AR cavity
g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=40, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, **GEO)
phi = g['phi']; sub_top = GEO['sub_top']
verts, faces, centroids, areas = t3.extract_mesh_3d(phi, DX)
F = len(faces); A = np.maximum(areas, 0.3 * np.median(areas))
A_src = GEO['Lx'] * GEO['Ly']; z_src = GEO['Lz'] - DX
mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV),
               indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
rng = np.random.default_rng(1)
N = 60000


def trace_ion():
    o, d = t3._source3d('ion', N, GEO['Lx'], GEO['Ly'], z_src, petch.PAR['ion_ang_sigma'], 'sobol', rng, 11)
    fi = wp.zeros(F, dtype=float, device=DEV); ai = wp.zeros(F, dtype=float, device=DEV)
    wp.launch(t3._trace3d, dim=N, device=DEV, inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEV),
              wp.array(d, dtype=wp.vec3, device=DEV), 1.0, 0, 11, 0, 0.34, 0.8, fi, ai])
    return np.clip((fi.numpy() / A) / (N / A_src), 0.0, 1.5)


def trace_elec():
    o, d = t3._source3d('neutral', N, GEO['Lx'], GEO['Ly'], z_src, 0.04, 'sobol', rng, 13)   # cosine (diffuse)
    fe = wp.zeros(F, dtype=float, device=DEV)
    wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, wp.array(o, dtype=wp.vec3, device=DEV),
              wp.array(d, dtype=wp.vec3, device=DEV), wp.array(np.ones(F, np.float32), dtype=float, device=DEV),
              1.0, 13, fe])
    return np.clip((fe.numpy() / A) / (N / A_src), 0.0, 8.0)


m_i = trace_ion(); m_e = trace_elec()
# faces inside the hole (near the axis), local AR = depth-below-mask / width
x, y, z = centroids[:, 0], centroids[:, 1], centroids[:, 2]
r = np.sqrt((x - GEO['Lx']/2)**2 + (y - GEO['Ly']/2)**2)
inside = (r < W * 0.7) & (z < sub_top - 0.5)
ar = (sub_top - z) / W
ratio = np.clip(m_e / np.maximum(m_i, 1e-6), 0.0, 1.0)

print(f"device={DEV}  deep hole AR_max~{ar[inside].max():.1f}  faces_inside={inside.sum()}\n", flush=True)
print("  Hwang-Giapis: floor ion current ~1.0/0.72/0.50/0.42 at AR 0/2/4/6 (electrons shadowed -> floor throttled)\n")
print("   AR   <m_e/m_i>   implied f_charge(alpha=0.85)   n_faces")
HG = {0: 1.0, 2: 0.72, 4: 0.50, 6: 0.42, 8: 0.38}
for arc in [0, 2, 4, 6, 8]:
    sel = inside & (np.abs(ar - arc) < 1.0)
    if sel.sum() < 3:
        continue
    rr = ratio[sel].mean()
    fch = 1.0 - 0.85 * (1.0 - rr)
    print(f"  {arc:3d}   {rr:6.3f}      {fch:6.3f}   (HG {HG[arc]:.2f})     {sel.sum()}", flush=True)
print("\n  if f_charge(AR) tracks the HG column, the charging model reproduces the measured floor-current rolloff.")
