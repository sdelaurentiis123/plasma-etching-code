#!/usr/bin/env python3
"""Phase 1 moat DEMO: data-efficient gradient calibration via the exact radiosity adjoint.

The differentiable-engine wedge is calibrating MANY declared surface parameters from limited data,
cheaply. The adjoint (reverse mode) computes the FULL gradient w.r.t. all K parameters in ONE extra
linear solve, independent of K; a derivative-free / finite-difference search needs O(K) solves per
gradient. So gradient calibration cost is ~flat in K while derivative-free grows linearly -- the
data-efficiency advantage that beats a derivative-free baseline (e.g. Krueger 2024) as the parameter
count rises.

Concretely: recover a K-band, depth-varying effective reaction-probability (sticking) map on the
trench walls from the floor incident-flux observations. Scalar loss L(theta)=1/2||H[floor]-target||^2;
gradient dL/dtheta via one adjoint solve `lam = M^{-T} residual`, `dL/dr = -H (B^T lam)`, summed per
band. Compares scipy BFGS with the analytic adjoint gradient vs BFGS with a 2-point finite-difference
gradient, counting radiosity solves. Run: python scripts/diff_calibration_demo.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy import sparse
from scipy.optimize import minimize
from scipy.sparse.linalg import spsolve

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from deboer_arde_static import make_rectangular_trench_geometry_3d, thermal_neutral_qmc_boundary_state
from petch.boundary_transport_3d import estimate_diffuse_form_factors_3d, gather_boundary_state_ballistic_3d
from petch.feature_step_3d import _face_material_ids, _surface_gas_normals
from petch.threed import extract_mesh_3d

SOLVES = [0]


def setup(ar=4.0, opening=0.10, dx=0.02, mask=0.05):
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
    return B, direct, floor, centroids[:, 2], flux0


def band_assign(z, K):
    """Assign each face to one of K depth bands (parameter -> faces map)."""
    lo, hi = z.min(), z.max()
    idx = np.clip(((z - lo) / (hi - lo + 1e-12) * K).astype(int), 0, K - 1)
    return idx


def solve_H(B, direct, r):
    # M = I - B diag(1-r): incident on target from reflected fraction (1-r) on source via B
    SOLVES[0] += 1
    M = (sparse.eye(B.shape[0], format="csr") - B @ sparse.diags(1.0 - r)).tocsc()
    return spsolve(M, direct), M


def make_problem(B, direct, floor, z, K):
    band = band_assign(z, K)

    def r_of(theta):
        return np.clip(theta[band], 0.0, 1.0)

    def obs(theta):
        H, _ = solve_H(B, direct, r_of(theta))
        return H[floor]

    def loss(theta):
        return 0.5 * float(np.sum((obs(theta) - target) ** 2))

    def grad(theta):
        r = r_of(theta)
        H, M = solve_H(B, direct, r)
        resid = np.zeros(B.shape[0]); resid[floor] = H[floor] - target
        SOLVES[0] += 1
        lam = spsolve(M.T.tocsc(), resid)
        dL_dr = -H * (B.T @ lam)                 # per-face gradient wrt reaction r
        return np.array([dL_dr[band == k].sum() for k in range(K)])

    theta_true = 0.06 + 0.14 * np.linspace(0, 1, K) ** 2   # depth-varying ground truth
    H_true, _ = solve_H(B, direct, np.clip(theta_true[band], 0, 1))
    target = H_true[floor]
    return loss, grad, theta_true, band


def main():
    B, direct, floor, z, flux0 = setup(ar=4.0)
    print("Phase 1 moat demo: recover a K-band sticking map from floor fluxes.")
    print("Adjoint gives all K gradients in ONE solve; finite-difference needs O(K).")
    print(f"{'K':>4} {'adjoint_solves':>15} {'FD_solves':>12} {'ratio':>7} {'adj_err':>9} {'fd_err':>9}")
    for K in (2, 6, 12):
        loss, grad, theta_true, band = make_problem(B, direct, floor, z, K)
        x0 = np.full(K, 0.10)
        SOLVES[0] = 0
        ra = minimize(loss, x0, jac=grad, method="L-BFGS-B",
                      bounds=[(0.0, 1.0)] * K, options={"maxiter": 60, "ftol": 1e-14, "gtol": 1e-10})
        adj_solves = SOLVES[0]
        adj_err = float(np.linalg.norm(ra.x - theta_true) / np.linalg.norm(theta_true))
        SOLVES[0] = 0
        rf = minimize(loss, x0, jac="2-point", method="L-BFGS-B",
                      bounds=[(0.0, 1.0)] * K, options={"maxiter": 60, "ftol": 1e-14, "gtol": 1e-10})
        fd_solves = SOLVES[0]
        fd_err = float(np.linalg.norm(rf.x - theta_true) / np.linalg.norm(theta_true))
        print(f"{K:>4} {adj_solves:>15} {fd_solves:>12} {fd_solves / max(adj_solves, 1):>7.1f} "
              f"{adj_err:>9.2e} {fd_err:>9.2e}", flush=True)
    print("\nadjoint solve-count is ~flat in K; finite-difference grows ~O(K) -> the calibration wedge.")


if __name__ == "__main__":
    main()
