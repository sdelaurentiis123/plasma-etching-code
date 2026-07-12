import numpy as np

import petch.charging_nodal_fixed_point as nodal_fixed_point_module
from petch.adaptive_quadrature import AdaptiveQuadratureResult
from petch.boundary_state import (
    PlasmaBoundaryState,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
)
from petch.charging_nodal_fixed_point import (
    _confidence_separated_log_ratio,
    solve_boundary_state_charging_nodal,
)
from petch.charging_backward import _gas_faces
from petch.charging_nodal import material_face_nodes


def _balanced_boundary():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    velocity = [[0.0, 0.0, 1.0]]
    ion = SpeciesBoundaryState(
        "ion", 1, 40.0, 1.0, velocity, [1.0], density_model=density)
    electron = SpeciesBoundaryState(
        "electron", -1, 5.485799e-4, 1.0, velocity, [1.0], density_model=density)
    return PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)


def test_confidence_separated_residual_has_exact_limit_and_unresolved_band():
    residual, ilo, ihi, elo, ehi = _confidence_separated_log_ratio(
        np.array([1.0, 2.0, 0.1, 2.0]),
        np.array([1.0, 1.0, 0.2, 1.0]),
        np.array([0.1, 0.1, 0.01, 0.0]),
        np.array([0.1, 0.1, 0.01, 0.0]), 2.0)

    assert residual[0] == 0.0
    assert np.isclose(residual[1], np.log(1.8 / 1.2))
    assert np.isclose(residual[2], np.log(0.12 / 0.18))
    assert np.isclose(residual[3], np.log(2.0))
    assert ilo[0] <= ehi[0] and elo[0] <= ihi[0]


def test_nodal_surface_fixed_point_balances_on_actual_boundary_vertices():
    nx, nz = 8, 6
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    result = solve_boundary_state_charging_nodal(
        solid, np.zeros_like(solid, dtype=int), _balanced_boundary(),
        n_iter=2, min_iter=1, balance_tol=1e-12,
        n_face_position=2, field_sweeps=100)
    assert result["surface_discretization"] == "boundary_nodal"
    assert result["iterations"] == 1
    assert result["interval_balance_final"]["max_abs_log_ratio"] == 0.0
    assert np.allclose(result["ion_current"], result["electron_current"], atol=1e-14)
    assert result["dielectric_nodes"].shape == (nx + 1, 2)
    assert np.allclose(result["boundary_nodal_voltage"], 0.0, atol=1e-14)


def test_nodal_charging_returns_last_evaluated_state_not_unassessed_step():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    velocity = [[0.0, 0.0, 1.0]]
    boundary = PlasmaBoundaryState((
        SpeciesBoundaryState(
            "ion", 1, 40.0, 2.0, velocity, [1.0], density_model=density),
        SpeciesBoundaryState(
            "electron", -1, 5.485799e-4, 1.0, velocity, [1.0], density_model=density),
    ), reference_plane_m=0.0)
    solid = np.zeros((6, 5), dtype=bool); solid[:, -1] = True
    result = solve_boundary_state_charging_nodal(
        solid, np.zeros_like(solid, dtype=int), boundary,
        n_iter=1, min_iter=1, balance_tol=None, beta=0.5,
        response_energy_eV=4.0, field_sweeps=50, trust_region=False)

    assert np.all(result["boundary_nodal_voltage"] == 0.0)
    assert np.all(result["surface_voltage"][solid] == 0.0)
    assert np.isclose(result["balance_final"]["max_abs_log_ratio"], np.log(2.0))
    assert np.allclose(result["ion_current"], 2.0 * result["electron_current"])
    assert np.all(result["ion_current"] > 0.0)


def test_endpoint_current_basis_removes_face_constant_alternating_null_mode():
    solid = np.zeros((20, 18), dtype=bool)
    solid[5, 1:] = True; solid[14, 1:] = True; solid[5:15, 15:] = True
    cells, normals = _gas_faces(solid, solid)
    endpoints = [material_face_nodes(cell, normal) for cell, normal in zip(cells, normals)]
    nodes = sorted({node for pair in endpoints for node in pair})
    index = {node: row for row, node in enumerate(nodes)}
    face_constant = np.zeros((len(nodes), len(endpoints)))
    endpoint_basis = np.zeros((len(nodes), 2 * len(endpoints)))
    for face, pair in enumerate(endpoints):
        for local, node in enumerate(pair):
            face_constant[index[node], face] = 0.5
            endpoint_basis[index[node], 2 * face + local] = 1.0
    assert np.linalg.matrix_rank(face_constant) < len(nodes)
    assert np.linalg.matrix_rank(endpoint_basis) == len(nodes)


def test_nodal_fixed_point_warm_starts_both_transport_estimators(monkeypatch):
    calls = []

    def hybrid(boundary, species_name, potential, solid, cells, normals, **kwargs):
        calls.append((
            species_name,
            kwargs["adjoint_options"].get("initial_log2_samples"),
            kwargs["forward_options"].get("initial_log2_samples")))
        face_count = len(cells); unique_cells = list(dict.fromkeys(cells))
        replicates = np.ones((4, face_count))
        endpoints = np.full((4, face_count, 2), 0.5)
        adjoint = AdaptiveQuadratureResult(
            np.ones(face_count), np.zeros(face_count), replicates,
            np.full(face_count, 9), 1.0, 0.0, True, 1, 1,
            auxiliary_mean=endpoints.mean(axis=0), auxiliary_replicates=endpoints)
        forward = AdaptiveQuadratureResult(
            np.ones(face_count), np.zeros(face_count), replicates,
            np.full(face_count, 11), 1.0, 0.0, True, 1, 1,
            auxiliary_mean=endpoints.mean(axis=0), auxiliary_replicates=endpoints)
        return dict(
            selected_face_mean=np.ones(face_count),
            selected_face_stderr=np.zeros(face_count),
            selected_face_replicates=replicates,
            selected_endpoint_stderr=np.zeros((face_count, 2)),
            selected_endpoint_replicates=endpoints,
            unique_cells=np.asarray(unique_cells),
            method=np.full(len(unique_cells), "adjoint"), converged=True,
            adjoint=adjoint, forward=forward)

    monkeypatch.setattr(
        nodal_fixed_point_module, "bidirectional_boundary_state_cell_flux", hybrid)
    solid = np.zeros((4, 4), dtype=bool); solid[:, -1] = True
    result = solve_boundary_state_charging_nodal(
        solid, np.zeros_like(solid, dtype=int), _balanced_boundary(),
        n_iter=2, min_iter=1, balance_tol=None, field_sweeps=20,
        trust_region=False, adaptive_quadrature=dict(
            bidirectional=True, base_log2=6, max_log2=12, n_replicates=4,
            forward_options=dict(base_log2=8, max_log2=14, n_replicates=4)))

    assert calls[0][1] is None and calls[0][2] is None
    assert calls[1][1] is None and calls[1][2] is None
    assert np.all(calls[2][1] == 9) and np.all(calls[2][2] == 11)
    assert np.all(calls[3][1] == 9) and np.all(calls[3][2] == 11)
    assert np.all(result["adaptive_levels"]["ion"] == 9)
    assert np.all(result["forward_adaptive_levels"]["electron"] == 11)
