"""ISOLATION TEST (petch side): instantaneous bottom etch rate on a CLEAN vertical trench carved to a
fixed depth D -- one flux+rate solve, NO etch dynamics (advection/reinit). Compare to vps_inst.py on the
identical geometry to separate flux-physics differences from etch-dynamics differences. cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import skfmm
import petch
from petch import threed as t3

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
MASK = 2 * DX
SUBTOP = 5.2                      # substrate top z (deep domain so the floor never clamps)
LZ = SUBTOP + MASK + 0.3
DS = [1.0, 1.85, 2.7, 3.55, 4.3]  # depths -> AR ~2.0, 3.7, 5.4, 7.1, 8.6
FLG = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                  sampling="sobol", neutral_transport="mc", ion_reflection=True)


def clean_trench(D):
    nx, ny, nz = round(XE / DX), round(YE / DX), round(LZ / DX)
    xs, ys, zs = np.arange(nx) * DX, np.arange(ny) * DX, np.arange(nz) * DX
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    opening = np.abs(X - XE / 2) < W / 2
    carve = opening & (Z >= SUBTOP - D) & (Z < SUBTOP)        # clean vertical trench of depth D
    solid = (Z < SUBTOP) & ~carve
    mask = ((Z >= SUBTOP) & (Z < SUBTOP + MASK)) & ~opening
    solid = solid | mask
    phi = skfmm.distance(np.where(solid, 1.0, -1.0), dx=DX)
    return dict(xs=xs, ys=ys, zs=zs, dx=DX, phi=phi, mask=mask, Lx=XE, Ly=YE, Lz=LZ,
                sub_top=SUBTOP, trench_width=W, hole=False)


def inst_rate(D):
    g = clean_trench(D)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['n_fp'] = 8
    verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
    mesh = t3.wp.Mesh(points=t3.wp.array(verts, dtype=t3.wp.vec3, device=t3.DEVICE),
                      indices=t3.wp.array(faces.flatten(), dtype=t3.wp.int32, device=t3.DEVICE))
    m_i, m_F, m_O, cos_i, bare = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, g, p,
                                                       n_ion=150000, n_neu=150000, seed=1,
                                                       sampling="sobol", flags=FLG)
    is_mask = t3.faces_in_mask(centroids, g, MASK, W, hole=False)
    V = t3.surface_rate(m_i, m_F, m_O, cos_i, is_mask, p, flags=FLG)
    V = np.nan_to_num(V, nan=0.0, posinf=0.0, neginf=0.0)
    cz = centroids[:, 2]; cx = centroids[:, 0]
    flo = (cz < cz.min() + 1.5 * DX) & (np.abs(cx - XE / 2) < 0.10) & (V > 0)
    return float(V[flo].mean()), float((SUBTOP - cz.min()) / W), int(flo.sum())


print("petch INSTANTANEOUS bottom rate on clean vertical trenches (no etch dynamics)\n", flush=True)
rows = [inst_rate(D) for D in DS]
v0 = rows[0][0]
print("  AR    nr(=v/v0)   floor_faces", flush=True)
for v, ar, n in rows:
    print(f"  {ar:4.1f}   {v / v0:.3f}        {n}", flush=True)
print("\n  compare to ViennaPS instantaneous (vps_inst.py) on the SAME clean trenches.", flush=True)
print("  petch ETCHED nr @AR8.6 was 0.602; ViennaPS ETCHED 0.731.", flush=True)
import json
json.dump([(ar, v / v0) for v, ar, n in rows], open("/root/petch_inst.json", "w"))
