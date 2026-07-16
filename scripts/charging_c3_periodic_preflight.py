#!/usr/bin/env python3
"""Standalone local/GPU preflight for the corrected C3 periodic charge operator."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import _geometry_and_poisson  # noqa: E402
from petch.charging_coupled_3d import _validate_periodic_topology_3d  # noqa: E402
from petch.charging_poisson import EPS0  # noqa: E402
from petch.charging_poisson_3d import (  # noqa: E402
    CompatibleQ1SurfaceChargeProjector3D,
    NodalPoissonSystem3D,
    lump_triangle_sheet_charge_3d,
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--warm-proposal", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-cuda", action="store_true")
    args = parser.parse_args()

    # Exact manufactured solution: a uniform sheet over a periodic parallel-plate capacitor.
    cell_shape = (3, 4, 7)
    spacing = np.asarray((13e-9, 19e-9, 23e-9))
    epsilon = np.full(cell_shape, 3.9)
    fixed = np.zeros(tuple(np.asarray(cell_shape) + 1), dtype=bool)
    fixed[:, :, -1] = True
    analytic = NodalPoissonSystem3D(
        epsilon, spacing, fixed, periodic_axes=(0, 1))
    sigma = 6.2e-4
    x_weight = np.full(cell_shape[0] + 1, spacing[0]); x_weight[[0, -1]] *= 0.5
    y_weight = np.full(cell_shape[1] + 1, spacing[1]); y_weight[[0, -1]] *= 0.5
    charge = np.zeros(analytic.shape)
    charge[:, :, 0] = sigma * x_weight[:, None] * y_weight[None, :]
    voltage, analytic_diagnostics = analytic.solve(charge)
    expected_surface = sigma * cell_shape[2] * spacing[2] / (EPS0 * 3.9)
    expected = np.linspace(expected_surface, 0.0, cell_shape[2] + 1)
    analytic_error_v = float(np.max(np.abs(voltage - expected[None, None, :])))
    _require(analytic_error_v < 2e-10, "periodic analytic capacitor gate failed")

    geometry, poisson = _geometry_and_poisson(0.25)
    _require(poisson.periodic_axes == (0, 1), "real-trench Poisson topology is not periodic")
    _require(np.array_equal(geometry.phi[0], geometry.phi[-1]),
             "real-trench geometry is not x-periodic")
    _require(np.array_equal(geometry.phi[:, 0], geometry.phi[:, -1]),
             "real-trench geometry is not y-periodic")
    random_charge = np.random.default_rng(20260714).normal(
        scale=2e-18, size=poisson.shape)
    random_charge[:, :, -1] = 0.0
    random_voltage, random_diagnostics = poisson.solve(random_charge)
    random_x_seam_v = float(np.max(np.abs(random_voltage[0] - random_voltage[-1])))
    random_y_seam_v = float(np.max(np.abs(random_voltage[:, 0] - random_voltage[:, -1])))
    _require(random_x_seam_v == 0.0 and random_y_seam_v == 0.0,
             "periodic prolongation is not bitwise continuous")
    _validate_periodic_topology_3d(poisson, True)
    mismatch_refused = False
    try:
        _validate_periodic_topology_3d(poisson, False)
    except ValueError:
        mismatch_refused = True
    _require(mismatch_refused, "particle/field topology mismatch was not refused")

    warm_path = args.warm_proposal.resolve()
    with np.load(warm_path) as warm:
        warm_sigma = np.asarray(warm["sigma_c_per_m2"], dtype=float)
        warm_face_charge = np.asarray(warm["face_charge_c"], dtype=float)
        warm_node_charge = np.asarray(warm["charge_node_c"], dtype=float)
        warm_potential = np.asarray(warm["potential_v"], dtype=float)
        vertices = np.asarray(warm["vertices"])
        faces = np.asarray(warm["faces"], dtype=int)
        areas = np.asarray(warm["areas"], dtype=float)
        archived_axes = tuple(int(value) for value in np.asarray(
            warm["poisson_periodic_axes"]).ravel())
    _require(archived_axes == poisson.periodic_axes,
             "warm proposal periodic topology does not match engine")
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    face_inventory_error_c = float(np.max(np.abs(
        warm_face_charge - warm_sigma * physical_area)))
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson, vertices, faces, grid_spacing=geometry.dx,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    reconstructed_node = poisson.canonicalize_charge(lump_triangle_sheet_charge_3d(
        poisson.shape, vertices, faces, warm_sigma, grid_spacing=geometry.dx,
        coordinate_length_unit_m=geometry.mesh_length_unit_m))
    node_scale = max(float(np.sum(np.abs(warm_node_charge))), np.finfo(float).tiny)
    warm_node_error = float(
        np.sum(np.abs(reconstructed_node - warm_node_charge)) / node_scale)
    solved_warm_potential, warm_diagnostics = poisson.solve(warm_node_charge)
    warm_potential_scale = max(float(np.linalg.norm(warm_potential)), np.finfo(float).tiny)
    warm_potential_error = float(
        np.linalg.norm(solved_warm_potential - warm_potential) / warm_potential_scale)
    warm_x_seam_v = float(np.max(np.abs(solved_warm_potential[0] - solved_warm_potential[-1])))
    warm_y_seam_v = float(np.max(np.abs(
        solved_warm_potential[:, 0] - solved_warm_potential[:, -1])))
    warm_null_fraction = projector.unresolved_fraction(warm_face_charge)
    _require(face_inventory_error_c < 2e-29, "warm face inventory is inconsistent")
    _require(warm_node_error < 5e-13, "warm face/node projection is inconsistent")
    _require(warm_potential_error < 5e-13, "warm saved/recomputed potential is inconsistent")
    _require(warm_x_seam_v == 0.0 and warm_y_seam_v == 0.0,
             "warm proposal contains a periodic field seam")
    _require(warm_null_fraction < 5e-13,
             "warm proposal retains periodic-Q1-null face inventory")

    import warp as wp
    cuda_available = bool(wp.is_cuda_available())
    if args.require_cuda:
        _require(cuda_available, "CUDA was required but Warp reports no CUDA device")
    device = str(wp.get_device("cuda:0")) if cuda_available else str(wp.get_device("cpu"))

    artifact = dict(
        schema="petch.charging.c3.periodic-preflight.v1",
        passed=True,
        analytic_parallel_plate_max_error_v=analytic_error_v,
        analytic_poisson_charge_balance_c=analytic_diagnostics.charge_balance_c,
        random_charge_x_seam_v=random_x_seam_v,
        random_charge_y_seam_v=random_y_seam_v,
        random_poisson_charge_balance_c=random_diagnostics.charge_balance_c,
        mismatch_refused=mismatch_refused,
        geometry_periodic_x=True, geometry_periodic_y=True,
        poisson_periodic_axes=list(poisson.periodic_axes),
        poisson_independent_node_shape=list(poisson.reduced_shape),
        periodic_face_coupling_rank=projector.rank,
        periodic_face_coupling_nullity=projector.nullity,
        warm_face_inventory_max_error_c=face_inventory_error_c,
        warm_node_relative_l1_error=warm_node_error,
        warm_potential_relative_l2_error=warm_potential_error,
        warm_x_seam_v=warm_x_seam_v, warm_y_seam_v=warm_y_seam_v,
        warm_periodic_null_fraction=warm_null_fraction,
        warm_poisson_charge_balance_c=warm_diagnostics.charge_balance_c,
        cuda_available=cuda_available, selected_device=device,
        provenance=dict(
            warm_proposal_sha256=_hash(warm_path),
            charging_poisson_sha256=_hash(ROOT / "src/petch/charging_poisson_3d.py"),
            charging_coupled_sha256=_hash(ROOT / "src/petch/charging_coupled_3d.py"),
            charging_coevolution_sha256=_hash(
                ROOT / "src/petch/charging_coevolution_3d.py"),
            script_sha256=_hash(Path(__file__).resolve())))
    encoded = json.dumps(artifact, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
