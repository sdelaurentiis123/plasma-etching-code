#!/usr/bin/env python3
"""TIME-INDEPENDENT ARDE: build STATIC clean trenches at a range of aspect ratios (smooth SDF), compute
the INSTANTANEOUS normalized bottom etch rate from the flux+chemistry directly -- no etching, no time, no
early-stop, no depth metric. nr(AR) = V_floor / V_openfield is a pure geometric property. Report it across
dx to show grid-convergence. This is how a research tool should characterize ARDE. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np, warp as wp
import skfmm
import petch
from petch import threed as t3

W, XE, YE = 0.5, 2.0, 1.0


def clean_trench_phi(dx, D, sub_top):
    """Smooth-SDF clean trench of depth D in a flat substrate (phi>0 solid). Floor at z=sub_top-D."""
    Lz = sub_top + 0.4
    nx, ny, nz = round(XE/dx), round(YE/dx), round(Lz/dx)
    xs, ys, zs = np.arange(nx)*dx, np.arange(ny)*dx, np.arange(nz)*dx
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    r = np.abs(X - XE/2)
    # graphics convention f<0 inside solid: substrate minus trench-gas-box
    f_sub = Z - sub_top                                  # <0 below sub_top (substrate)
    f_box = np.maximum(r - W/2, (sub_top - D) - Z)       # <0 inside trench box (r<W/2 & z>sub_top-D)
    f = np.maximum(f_sub, -f_box)                        # substrate AND NOT box
    phi = skfmm.distance(-f, dx=dx)                      # phi>0 solid, sub-cell zero contour
    geo = dict(xs=xs, ys=ys, zs=zs, dx=dx, phi=phi, Lx=XE, Ly=YE, Lz=Lz, sub_top=sub_top,
               trench_width=W, hole=False, mask=np.zeros_like(phi, bool))
    return geo


def bottom_rate(dx, D, sub_top, n_rays=200000):
    geo = clean_trench_phi(dx, D, sub_top)
    verts, faces, cen, areas = t3.extract_mesh_3d(geo['phi'], dx)
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=t3.DEVICE),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=t3.DEVICE))
    par = dict(petch.PAR); fl = petch.Flags(chemistry="belen", yield_angular="viennaps",
                                            coverage_sticking=True, sampling="sobol", ion_reflection=True)
    bi = None
    for _it in range(5):                                # converge the coupled-coverage fixed point
        m_i, m_F, m_O, cos_i, cov = t3.mc_flux_3d_coupled(mesh, verts, faces, areas, geo, par,
                                                          n_ion=n_rays, n_neu=n_rays, seed=_it,
                                                          sampling="sobol", flags=fl, bare_init=bi)
        bi = cov
    is_mask = np.zeros(len(faces), bool)
    V = t3.surface_rate(m_i, m_F, m_O, cos_i, is_mask, par, flags=fl)
    cx = cen[:, 0]-XE/2; cz = cen[:, 2]
    floor = (np.abs(cx) < 0.15) & (cz < (sub_top - D) + 3*dx)         # trench floor facets
    field = cz > sub_top - 0.5*dx                                     # open field (substrate top, r>W/2)
    vf = np.nanmean(V[floor]) if floor.sum() else np.nan
    v0 = np.nanmean(V[field]) if field.sum() else np.nan
    return vf, v0, vf/v0 if v0 else np.nan, int(floor.sum())


SUB = 8.0
ARs = [2, 4, 6, 8, 10]
print(f"device={t3.DEVICE}   STATIC instantaneous ARDE nr(AR)=V_floor/V_field\n", flush=True)
for dx in [0.04, 0.025]:
    print(f"dx={dx}:", flush=True)
    for AR in ARs:
        vf, v0, nr, nfl = bottom_rate(dx, AR*W, SUB)
        print(f"   AR {AR:2d}  nr={nr:.3f}   (floor faces {nfl})", flush=True)
print("\n  Time-independent, deterministic. ViennaPS deep-AR (AR8.6) ~ 0.73.", flush=True)
