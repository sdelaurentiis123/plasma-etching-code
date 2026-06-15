#!/usr/bin/env python3
"""Validate surface charging vs Hwang-Giapis 1997: the FLOOR ion current drops with aspect ratio
(~60% by AR~4). Correct test = etch holes to INCREASING depth, and for each measure the charging factor
at the floor vs the floor's AR, anchored to the open field (most-exposed faces). Run on a box: cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DEV = t3.DEVICE
DX, W = 0.25, 4.0
GEO = dict(Lx=12, Ly=12, Lz=40, mask_th=2, sub_top=34, hole=True, t_end=None)
N = 60000


def floor_fcharge(phi):
    verts, faces, c, areas = t3.extract_mesh_3d(phi, DX)
    F = len(faces); A = np.maximum(areas, 0.3 * np.median(areas))
    A_src = GEO['Lx'] * GEO['Ly']; z_src = GEO['Lz'] - DX
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
    rng = np.random.default_rng(1)
    oi, di = t3._source3d('ion', N, GEO['Lx'], GEO['Ly'], z_src, petch.PAR['ion_ang_sigma'], 'sobol', rng, 11)
    fi = wp.zeros(F, dtype=float, device=DEV); ai = wp.zeros(F, dtype=float, device=DEV)
    wp.launch(t3._trace3d, dim=N, device=DEV, inputs=[mesh.id, wp.array(oi, dtype=wp.vec3, device=DEV),
              wp.array(di, dtype=wp.vec3, device=DEV), 1.0, 0, 11, 0, 0.34, 0.8, fi, ai])
    m_i = np.clip((fi.numpy() / A) / (N / A_src), 0.0, 1.5)
    oe, de = t3._source3d('neutral', N, GEO['Lx'], GEO['Ly'], z_src, 0.04, 'sobol', rng, 13)
    fe = wp.zeros(F, dtype=float, device=DEV)
    wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, wp.array(oe, dtype=wp.vec3, device=DEV),
              wp.array(de, dtype=wp.vec3, device=DEV), wp.array(np.ones(F, np.float32), dtype=float, device=DEV),
              1.0, 13, fe])
    m_e = np.clip((fe.numpy() / A) / (N / A_src), 0.0, 8.0)
    ref_i = np.percentile(m_i[m_i > 1e-6], 90); ref_e = np.percentile(m_e[m_e > 1e-6], 90)
    sh_i = np.clip(m_i / ref_i, 0, 1); sh_e = np.clip(m_e / ref_e, 0, 1)
    ratio = np.clip(sh_e / np.maximum(sh_i, 1e-3), 0, 1)
    x, y, z = c[:, 0], c[:, 1], c[:, 2]
    r = np.sqrt((x - GEO['Lx']/2)**2 + (y - GEO['Ly']/2)**2)
    zf = z[(r < W * 0.55)].min() if (r < W * 0.55).any() else GEO['sub_top']
    floor = (r < W * 0.55) & (z < zf + 1.2)                 # deepest central faces = the floor
    ar = (GEO['sub_top'] - z[floor].mean()) / W
    return ar, float(ratio[floor].mean())


print(f"device={DEV}  charging vs Hwang-Giapis (floor current vs AR)\n", flush=True)
print("  HG: floor ion current ~1.0/0.72/0.50/0.42 at AR 0/2/4/6\n")
print("   steps   AR_floor   <shadow_e/shadow_i>   f_charge(a=0.85)", flush=True)
pts = []
for ns in [4, 8, 14, 20, 28, 36]:
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=ns, par=dict(petch.PAR, n_fp=1),
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, t_end=ns * 0.06, **{k: v for k, v in GEO.items() if k != 't_end'})
    ar, rr = floor_fcharge(g['phi'])
    fch = 1.0 - 0.85 * (1.0 - rr)
    pts.append((ar, rr, fch))
    print(f"   {ns:3d}     {ar:5.1f}      {rr:6.3f}              {fch:6.3f}", flush=True)
print("\n  if f_charge falls ~1.0->0.4 as AR_floor goes 0->4, the charging model matches Hwang-Giapis.")
