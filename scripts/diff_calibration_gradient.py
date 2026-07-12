#!/usr/bin/env python3
"""Phase 1a: exact CALIBRATION gradients via implicit differentiation of the radiosity fixed point.

The near-term differentiable-engine moat is gradient-based CALIBRATION of the declared uncertain
surface parameters (sticking, yields) -- calibrate on structure N, predict N+1, faster and more
data-efficient than a derivative-free search. Those parameters enter the diffuse-radiosity operator
M(s) = I - (1-s) B SMOOTHLY (unlike GEOMETRY parameters, which move ray-hit/visibility boundaries and
are the genuinely discontinuous, harder case). So the calibration gradient is exact:

    floor = c^T H,   M(s) H = D,   M = I - (1-s) B  ->  dM/ds = +B
    d(floor)/ds = -c^T M^{-1} B H = -(M^{-T} c)^T (B H)          [one adjoint linear solve]

This matches central finite difference to ~1e-7, i.e. it is exact, not an approximation. It is the
implicit-function-theorem adjoint through the (linear) radiosity fixed point -- the same structure that
will wrap the (nonlinear) charging fixed point later. GEOMETRY/shape gradients remain the open hard
problem (discontinuity treatment); this closes the calibration half. Run:
python scripts/diff_calibration_gradient.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from deboer_arde_static import make_rectangular_trench_geometry_3d, thermal_neutral_qmc_boundary_state
from petch.boundary_transport_3d import estimate_diffuse_form_factors_3d, gather_boundary_state_ballistic_3d
from petch.feature_step_3d import _face_material_ids, _surface_gas_normals
from petch.threed import extract_mesh_3d


def radiosity_system(ar, *, opening=0.10, dx=0.02, mask=0.05):
    """Build the diffuse-radiosity exchange operator B, direct flux D, floor functional c."""
    etched = ar * opening
    substrate_top = etched + max(4.0 * dx, 0.05)
    g = make_rectangular_trench_geometry_3d(
        cell_width=2.0 * opening, cell_length=max(6.0 * dx, 0.06),
        domain_height=substrate_top + mask + max(6.0 * dx, 0.06), dx=dx,
        opening_width=opening, mask_thickness=mask, substrate_top=substrate_top, etched_depth=etched)
    verts, faces, centroids, areas = extract_mesh_3d(g.phi, g.dx)
    gn = _surface_gas_normals(verts, faces, centroids, g)
    domain = (np.asarray(g.phi.shape) - 1) * g.dx
    src_z = float(domain[2])
    flux0 = 1.0e20
    boundary = thermal_neutral_qmc_boundary_state("F", 19.0, 0.05, flux0, log2_samples=16,
                                                  reference_plane_m=src_z * g.mesh_length_unit_m)
    tr = gather_boundary_state_ballistic_3d(
        boundary, {"F": "neutral_reactant"}, verts, faces, areas, centroids, gn,
        source_bounds=(0.0, float(domain[0]), 0.0, float(domain[1])), source_z=src_z,
        mesh_length_unit_m=g.mesh_length_unit_m, mesh_origin_m=g.mesh_origin_m,
        face_quadrature_points=3, periodic_lateral=True, domain_size=domain, ray_offset=1e-3 * g.dx)
    direct = np.asarray(tr.surface_fluxes.neutral_flux_m2_s["F"], float)
    ff = estimate_diffuse_form_factors_3d(verts, faces, centroids, gn, rays_per_face=64,
                                          domain_size=domain, periodic_lateral=True, ray_offset=1e-3 * g.dx)
    pa = areas * g.mesh_length_unit_m ** 2
    mat = _face_material_ids(centroids, g)
    up = gn[:, 2] > 0.5
    floor_z = centroids[(mat == 1) & up, 2].min()
    floor = (mat == 1) & up & (centroids[:, 2] <= floor_z + g.dx)
    B = sparse.coo_matrix((ff.transfer_fraction * pa[ff.source_face] / pa[ff.target_face],
                           (ff.target_face, ff.source_face)), shape=(len(faces), len(faces))).tocsr()
    c = np.zeros(len(faces)); c[floor] = pa[floor] / pa[floor].sum()
    return B, direct, c, flux0


def floor_transmission(B, direct, c, s):
    M = (sparse.eye(B.shape[0], format="csr") - (1.0 - s) * B).tocsc()
    H = spsolve(M, direct)
    return float(c @ H), M, H


def sticking_gradient(B, c, M, H):
    """Exact d(floor flux)/ds via one adjoint solve: -(M^{-T} c)^T (B H)."""
    lam = spsolve(M.T, c)
    return float(-lam @ (B @ H))


def main():
    print("Phase 1a: exact calibration gradient d(floor flux)/d(sticking) -- adjoint vs central FD")
    print(f"{'AR':>4} {'s':>5} {'floorT':>9} {'analytic':>13} {'central_FD':>13} {'rel_err':>10}")
    for ar in (1.0, 4.0):
        B, direct, c, flux0 = radiosity_system(ar)
        for s in (0.05, 0.2):
            f, M, H = floor_transmission(B, direct, c, s)
            g_an = sticking_gradient(B, c, M, H) / flux0
            h = 1e-4
            g_fd = (floor_transmission(B, direct, c, s + h)[0]
                    - floor_transmission(B, direct, c, s - h)[0]) / (2 * h) / flux0
            rel = abs(g_an - g_fd) / max(abs(g_fd), 1e-30)
            print(f"{ar:>4.0f} {s:>5.2f} {f/flux0:>9.4f} {g_an:>13.5e} {g_fd:>13.5e} {rel:>10.2e}",
                  flush=True)


if __name__ == "__main__":
    main()
