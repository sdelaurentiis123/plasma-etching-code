"""petch floor flux by channel (etchant m_F vs ion-enhanced ionFlux*ie_n) vs AR on clean trenches --
to compare per-channel rolloff against ViennaPS (vps_flux.py) and pinpoint which channel is deficient."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import skfmm
import petch
from petch import threed as t3

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
MASK = 2 * DX; SUBTOP = 5.2; LZ = SUBTOP + MASK + 0.3
DS = [1.0, 1.85, 2.7, 3.55, 4.3]
FLG = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                  sampling="sobol", neutral_transport="mc", ion_reflection=True)


def clean_trench(D):
    nx, ny, nz = round(XE / DX), round(YE / DX), round(LZ / DX)
    xs, ys, zs = np.arange(nx) * DX, np.arange(ny) * DX, np.arange(nz) * DX
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    opening = np.abs(X - XE / 2) < W / 2
    carve = opening & (Z >= SUBTOP - D) & (Z < SUBTOP)
    solid = (Z < SUBTOP) & ~carve
    mask = ((Z >= SUBTOP) & (Z < SUBTOP + MASK)) & ~opening
    phi = skfmm.distance(np.where(solid | mask, 1.0, -1.0), dx=DX)
    return dict(xs=xs, ys=ys, zs=zs, dx=DX, phi=phi, mask=(mask), Lx=XE, Ly=YE, Lz=LZ,
                sub_top=SUBTOP, trench_width=W, hole=False)


def floor_flux(D):
    g = clean_trench(D)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['n_fp'] = 16
    verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
    mesh = t3.wp.Mesh(points=t3.wp.array(verts, dtype=t3.wp.vec3, device=t3.DEVICE),
                      indices=t3.wp.array(faces.flatten(), dtype=t3.wp.int32, device=t3.DEVICE))
    m_i, m_F, m_O, cos_i, bare = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, g, p, n_ion=150000,
                                                       n_neu=150000, seed=1, sampling="sobol", flags=FLG)
    ie_n = p['_ion_yield'][1]; ienh = p['ionFlux'] * ie_n
    cz = centroids[:, 2]; cx = centroids[:, 0]
    flo = (cz < cz.min() + 1.5 * DX) & (np.abs(cx - XE / 2) < 0.10)
    fld = (np.abs(cx - XE / 2) > 0.5) & (cz > SUBTOP - 0.1)
    return (m_F[flo].mean(), ienh[flo].mean(), m_F[fld].mean() if fld.any() else 1, ienh[fld].mean() if fld.any() else 1)


print("petch floor flux by channel vs AR (clean trench)\n", flush=True)
rows = [floor_flux(D) for D in DS]
e0, i0 = rows[0][0], rows[0][1]
print("  AR   etchant(m_F)   ionEnhanced   | field_etch field_ienh", flush=True)
for (e, ie, fe, fie), D in zip(rows, DS):
    print(f"  {D/W:4.1f}   {e/e0:.3f}         {ie/i0:.3f}        |  {fe:.2f}  {fie:.2f}", flush=True)
print("\n  compare per-channel rolloff to ViennaPS (vps_flux.py).", flush=True)
