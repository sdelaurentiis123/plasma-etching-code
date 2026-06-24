"""Is petch's clean-trench floor-flux deficit a COVERAGE-CONVERGENCE artifact? Sweep n_fp (cold fixed-
point iterations) on the identical clean vertical trenches. If deep-AR nr RISES with n_fp toward ViennaPS
(1, 0.924, 0.841, 0.745, ~0.73), the cold fixed point was under-converged (sidewall coverage too low ->
neutrals over-stick -> floor starved). cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import skfmm
import petch
from petch import threed as t3

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
MASK = 2 * DX
SUBTOP = 5.2
LZ = SUBTOP + MASK + 0.3
DS = [1.0, 1.85, 2.7, 3.55, 4.3]            # AR 2.0, 3.7, 5.4, 7.1, 8.6
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
    solid = solid | mask
    phi = skfmm.distance(np.where(solid, 1.0, -1.0), dx=DX)
    return dict(xs=xs, ys=ys, zs=zs, dx=DX, phi=phi, mask=mask, Lx=XE, Ly=YE, Lz=LZ,
                sub_top=SUBTOP, trench_width=W, hole=False)


def inst_rate(D, n_fp):
    g = clean_trench(D)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['n_fp'] = n_fp
    verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
    mesh = t3.wp.Mesh(points=t3.wp.array(verts, dtype=t3.wp.vec3, device=t3.DEVICE),
                      indices=t3.wp.array(faces.flatten(), dtype=t3.wp.int32, device=t3.DEVICE))
    m_i, m_F, m_O, cos_i, bare = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, g, p, n_ion=150000,
                                                       n_neu=150000, seed=1, sampling="sobol", flags=FLG)
    is_mask = t3.faces_in_mask(centroids, g, MASK, W, hole=False)
    V = np.nan_to_num(t3.surface_rate(m_i, m_F, m_O, cos_i, is_mask, p, flags=FLG), nan=0.0, posinf=0.0)
    cz = centroids[:, 2]; cx = centroids[:, 0]
    flo = (cz < cz.min() + 1.5 * DX) & (np.abs(cx - XE / 2) < 0.10) & (V > 0)
    # also report mean neutral F flux at the floor to see if the lever is neutral coverage
    return float(V[flo].mean()), float(m_F[flo].mean())


print("petch clean-trench nr vs n_fp (ViennaPS: 1, 0.924, 0.841, 0.745, ~0.73)\n", flush=True)
for n_fp in [4, 8, 16, 30]:
    vs = [inst_rate(D, n_fp) for D in DS]
    v0 = vs[0][0]
    nr = [round(v / v0, 3) for v, mf in vs]
    mF = [round(mf, 3) for v, mf in vs]
    print(f"  n_fp={n_fp:2d}: nr {nr}   floor_mF {mF}", flush=True)
print("\ndone", flush=True)
