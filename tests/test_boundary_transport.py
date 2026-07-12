import numpy as np
import pytest

import petch.boundary_transport as boundary_transport_module
import petch.charging_backward as charging_backward_module
from petch.adaptive_quadrature import AdaptiveQuadratureResult

from petch.boundary_state import (
    PlasmaBoundaryState,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
    instantaneous_sinusoidal_ion_boundary_state,
    maxwellian_electron_boundary_state,
)
from petch.boundary_transport import (
    adaptive_adjoint_boundary_state_face_flux,
    adaptive_forward_boundary_state_cell_flux,
    bidirectional_boundary_state_cell_flux,
    adjoint_boundary_state_face_flux,
    adjoint_boundary_state_floor_flux,
    boundary_launches_2d,
    forward_boundary_state_cell_flux_qmc,
    trace_boundary_state_floor_flux,
)
from petch.charging2d import _build_edge_array_geometry
from petch.charging_backward import self_consistent_backward, solve_boundary_state_charging


def test_boundary_launcher_preserves_probability_flux_and_joint_energy():
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 2e19,
        velocity_sqrt_eV=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        weight=np.array([0.25, 0.75]),
    )
    launches = boundary_launches_2d(species, 2.0, 6.0, 8)
    assert np.isclose(launches.normalized_weight.sum(), 1.0)
    assert np.isclose(launches.flux_weight_m2_s.sum(), 2e19)
    launched_energy = launches.vx ** 2 + launches.vy ** 2 + launches.vz ** 2
    assert np.isclose(np.dot(launches.normalized_weight, launched_energy), species.mean_energy_eV)
    assert launches.x.min() > 2.0 and launches.x.max() < 6.0


def test_boundary_state_transport_matches_open_vertical_ion_flux():
    nx, nz = 24, 18
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    target = np.zeros_like(solid); target[:, -1] = True
    V = np.zeros((nx + 1, nz + 1))
    boundary = instantaneous_sinusoidal_ion_boundary_state(
        37.0, 30.0, 4.0, 40.0, 3e19, n_phase=64, ion_name="Ar+")
    result = trace_boundary_state_floor_flux(
        boundary, "Ar+", V, solid, target, n_position=64)
    assert np.isclose(result["normalized_flux"], 1.0, atol=1e-12)
    assert np.isclose(result["absolute_flux_m2_s"], 3e19, rtol=1e-12)


def test_forward_qmc_cell_flux_matches_open_surface_cell_measure():
    nx, nz = 16, 10
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-1e-9, 1e-9]), np.array([-1e-9, 1e-9]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "vertical", 0, 40.0, 2e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(x, nz - 1) for x in range(nx)]
    result = forward_boundary_state_cell_flux_qmc(
        boundary, "vertical", np.zeros((nx + 1, nz + 1)), solid, cells,
        log2_samples=10, seed=31)
    assert np.allclose(result["per_cell"], 1.0, atol=0.03)
    assert np.isclose(result["normalized_total"], nx, atol=1e-12)


def test_same_transport_adapter_accepts_neutral_reactive_species():
    nx, nz = 16, 12
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    target = np.zeros_like(solid); target[:, -1] = True
    neutral = SpeciesBoundaryState(
        "CF2", 0, 50.0, 7e20,
        velocity_sqrt_eV=np.array([[0.0, 0.0, 0.2], [0.05, 0.0, 0.3]]),
        weight=np.array([0.6, 0.4]), provenance={"source": "reactor"})
    boundary = PlasmaBoundaryState((neutral,), reference_plane_m=0.0)
    result = trace_boundary_state_floor_flux(
        boundary, "CF2", np.zeros((nx + 1, nz + 1)), solid, target, n_position=32,
        max_steps=20000)
    assert np.isclose(result["normalized_flux"], 1.0, atol=1e-12)
    assert np.isclose(result["absolute_flux_m2_s"], 7e20, rtol=1e-12)


@pytest.mark.parametrize("aspect_ratio", [1, 4, 16])
def test_same_boundary_transport_engine_spans_aspect_ratio_ladder(aspect_ratio):
    width = 8; depth = aspect_ratio * width
    nx = 3 * width; nz = depth + 2
    left, right, floor = width, 2 * width, depth
    solid = np.zeros((nx, nz), dtype=bool)
    solid[left - 1, :floor + 1] = True
    solid[right, :floor + 1] = True
    solid[left - 1:right + 1, floor] = True
    target = np.zeros_like(solid); target[left:right, floor] = True
    vertical = SpeciesBoundaryState("test", 0, 40.0, 1e19, [[0.0, 0.0, 1.0]], [1.0])
    boundary = PlasmaBoundaryState((vertical,), reference_plane_m=0.0)
    result = trace_boundary_state_floor_flux(
        boundary, "test", np.zeros((nx + 1, nz + 1)), solid, target,
        n_position=3 * width, max_steps=1000 * nz)
    assert np.isclose(result["normalized_flux"], 1.0, atol=1e-12)


@pytest.mark.parametrize("charge_number,name", [(1, "Ar+"), (0, "CF2")])
def test_same_boundary_density_drives_ion_and_neutral_adjoint(charge_number, name):
    nx, nz = 12, 10
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    V = np.zeros((nx + 1, nz + 1))
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.5, 0.5]), np.array([-0.5, 0.5]), np.array([0.5, 1.5])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        name, charge_number, 40.0, 2e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(x, nz - 1) for x in range(nx)]
    result = adjoint_boundary_state_floor_flux(boundary, name, V, solid, cells, n_face_position=4)
    assert np.isclose(result["normalized_flux"], 1.0, atol=1e-12)
    assert np.isclose(result["absolute_flux_m2_s"], 2e19, rtol=1e-12)


def test_adjoint_preserves_phase_label_when_scoring_plasma_exit():
    class PhaseRequiredDensity:
        def log_flux_density(self, velocity_sqrt_eV, phase_rad=None, position_m=None):
            if phase_rad is None:
                raise ValueError("phase label was lost")
            velocity_sqrt_eV = np.asarray(velocity_sqrt_eV)
            phase_rad = np.asarray(phase_rad)
            return np.log((1.0 + 0.25 * np.cos(phase_rad)) / (2.0 * np.pi)) + np.zeros(
                velocity_sqrt_eV.shape[:-1])

    nx, nz = 8, 6
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    phase = np.array([0.25, 2.25])
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 1e19, [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], [0.5, 0.5],
        phase_rad=phase, density_model=PhaseRequiredDensity())
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    result = adjoint_boundary_state_floor_flux(
        boundary, "ion", np.zeros((nx + 1, nz + 1)), solid,
        [(x, nz - 1) for x in range(nx)], n_face_position=2)
    assert np.isclose(result["normalized_flux"], 1.0, atol=1e-12)


def test_adaptive_adjoint_refines_generic_face_phase_space_with_error_evidence():
    nx, nz = 12, 10
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.5, 0.5]), np.array([-0.5, 0.5]), np.array([0.5, 1.5])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 1e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(x, nz - 1) for x in range(nx)]
    result = adaptive_adjoint_boundary_state_face_flux(
        boundary, "ion", np.zeros((nx + 1, nz + 1)), solid, cells,
        [(0.0, -1.0)] * len(cells), base_log2=4, max_log2=8, n_replicates=3,
        absolute_tolerance=2e-3, relative_tolerance=0.0,
        element_absolute_tolerance=2e-3, n_face_position=2)
    assert result.converged
    assert np.isclose(result.total_mean, 1.0, atol=2e-3)


def test_bidirectional_estimator_selects_by_uncertainty_not_region_name():
    nx, nz = 12, 10
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.2, 0.2]), np.array([-0.2, 0.2]), np.array([0.8, 1.2])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "test", 1, 40.0, 1e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(x, nz - 1) for x in range(nx)]
    result = bidirectional_boundary_state_cell_flux(
        boundary, "test", np.zeros((nx + 1, nz + 1)), solid, cells,
        [(0.0, -1.0)] * len(cells),
        adjoint_options=dict(base_log2=4, max_log2=7, n_replicates=4),
        forward_options=dict(base_log2=7, max_log2=11, n_replicates=4),
        element_absolute_tolerance=0.04, element_relative_tolerance=0.0)
    assert result["converged"]
    assert np.allclose(result["per_face"], 1.0, atol=0.04)
    assert set(result["method"]) <= {"forward", "adjoint"}
    assert result["method_within_tolerance"].all()
    assert np.all(np.isfinite(result["estimator_discrepancy_sigma"]))
    assert result["forward_cell_mean"].shape == (len(cells),)
    assert result["adjoint_cell_mean"].shape == (len(cells),)


def test_bidirectional_hysteresis_never_retains_an_uncertified_estimator(monkeypatch):
    mean_adjoint = 0.1402
    stderr_adjoint = 0.0338
    delta = stderr_adjoint * np.sqrt(3.0)
    adjoint_replicates = np.array([
        [mean_adjoint - delta], [mean_adjoint - delta],
        [mean_adjoint + delta], [mean_adjoint + delta]])
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([mean_adjoint]),
        element_stderr=np.array([stderr_adjoint]),
        element_replicates=adjoint_replicates,
        log2_samples=np.array([12]), total_mean=mean_adjoint,
        total_stderr=stderr_adjoint, converged=True, rounds=1, evaluations=1)
    mean_forward = 0.1562
    stderr_forward = 0.0269
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([mean_forward]),
        element_stderr=np.array([stderr_forward]),
        element_replicates=np.array([[mean_forward]] * 4),
        log2_samples=np.array([12]), total_mean=mean_forward,
        total_stderr=stderr_forward, converged=True, rounds=1, evaluations=1)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)

    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], method_hint=np.array(["adjoint"]),
        element_absolute_tolerance=0.01, element_relative_tolerance=0.15,
        switch_factor=2.0)

    assert result["converged"]
    assert result["method"][0] == "forward"
    assert np.isclose(result["per_face"][0], mean_forward)


def test_bidirectional_frozen_estimator_refuses_instead_of_switching_map(monkeypatch):
    mean_adjoint = 0.1402
    stderr_adjoint = 0.0338
    delta = stderr_adjoint * np.sqrt(3.0)
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([mean_adjoint]),
        element_stderr=np.array([stderr_adjoint]),
        element_replicates=np.array([
            [mean_adjoint - delta], [mean_adjoint - delta],
            [mean_adjoint + delta], [mean_adjoint + delta]]),
        log2_samples=np.array([12]), total_mean=mean_adjoint,
        total_stderr=stderr_adjoint, converged=True, rounds=1, evaluations=1)
    mean_forward = 0.1562
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([mean_forward]), element_stderr=np.array([0.0269]),
        element_replicates=np.array([[mean_forward]] * 4),
        log2_samples=np.array([12]), total_mean=mean_forward,
        total_stderr=0.0269, converged=True, rounds=1, evaluations=1)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)

    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], method_hint=np.array(["adjoint"]),
        freeze_method_hint=True, element_absolute_tolerance=0.01,
        element_relative_tolerance=0.15, switch_factor=2.0)

    assert result["method_hint_frozen"]
    assert result["method"][0] == "adjoint"
    assert not result["method_within_tolerance"][0]
    assert not result["cell_converged"][0]
    assert not result["converged"]
    assert np.isclose(result["per_face"][0], mean_adjoint)


def test_bidirectional_adjoint_zero_is_not_exact_when_forward_sees_flux(monkeypatch):
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([0.0]), element_stderr=np.array([0.0]),
        element_replicates=np.zeros((4, 1)), log2_samples=np.array([12]),
        total_mean=0.0, total_stderr=0.0, converged=True, rounds=1, evaluations=1)
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([0.1]), element_stderr=np.array([0.01]),
        element_replicates=np.array([[0.09], [0.09], [0.11], [0.11]]),
        log2_samples=np.array([12]), total_mean=0.1, total_stderr=0.01,
        converged=True, rounds=1, evaluations=1)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)

    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], method_hint=np.array(["adjoint"]),
        element_absolute_tolerance=0.02, element_relative_tolerance=0.15)

    assert result["converged"]
    assert result["adjoint_zero_unresolved"][0]
    assert result["adjoint_support_unresolved"][0]
    assert result["method"][0] == "forward"
    assert np.isclose(result["per_face"][0], 0.1)


def test_bidirectional_adjoint_zero_keeps_forward_zero_hit_upper_bound(monkeypatch):
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([0.0]), element_stderr=np.array([0.0]),
        element_replicates=np.zeros((4, 1)), log2_samples=np.array([12]),
        total_mean=0.0, total_stderr=0.0, converged=True, rounds=1, evaluations=1)
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([0.0]), element_stderr=np.array([0.005]),
        element_replicates=np.zeros((4, 1)), log2_samples=np.array([12]),
        total_mean=0.0, total_stderr=0.005, converged=True, rounds=1, evaluations=1)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)

    result = bidirectional_boundary_state_cell_flux(
        object(), "electron", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], method_hint=np.array(["adjoint"]),
        element_absolute_tolerance=0.01, element_relative_tolerance=0.0)

    assert result["converged"]
    assert result["adjoint_zero_unresolved"][0]
    assert result["method"][0] == "forward"
    assert result["per_face"][0] == 0.0
    assert result["per_face_stderr"][0] == 0.005


def test_bidirectional_refuses_inconsistent_estimators_that_both_claim_precision(monkeypatch):
    def estimate(mean):
        delta = 0.001 * np.sqrt(3.0)
        return AdaptiveQuadratureResult(
            element_mean=np.array([mean]), element_stderr=np.array([0.001]),
            element_replicates=np.array([
                [mean - delta], [mean - delta], [mean + delta], [mean + delta]]),
            log2_samples=np.array([12]), total_mean=mean, total_stderr=0.001,
            converged=True, rounds=1, evaluations=1)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: estimate(0.10))
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: estimate(0.20))
    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], element_absolute_tolerance=0.01,
        element_relative_tolerance=0.0, consistency_sigma=5.0)

    assert result["method_within_tolerance"][0]
    assert not result["estimator_consistent"][0]
    assert not result["cell_converged"][0]
    assert not result["converged"]


def test_bidirectional_cross_refines_a_missed_adjoint_mode(monkeypatch):
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([0.05]), element_stderr=np.array([0.001]),
        element_replicates=np.array([[0.048], [0.050], [0.050], [0.052]]),
        log2_samples=np.array([12]), total_mean=0.05, total_stderr=0.001,
        converged=True, rounds=1, evaluations=1)
    levels_seen = []

    def adjoint_estimate(*args, **kwargs):
        levels = np.asarray(kwargs.get("initial_log2_samples", [6]), dtype=int)
        level = int(levels[0]); levels_seen.append(level)
        mean = 0.05 if level >= 10 else 0.0
        stderr = 0.001 if level >= 10 else 0.0
        delta = stderr * np.sqrt(3.0)
        replicates = np.array([
            [mean - delta], [mean - delta], [mean + delta], [mean + delta]])
        return AdaptiveQuadratureResult(
            element_mean=np.array([mean]), element_stderr=np.array([stderr]),
            element_replicates=replicates, log2_samples=np.array([level]),
            total_mean=mean, total_stderr=stderr, converged=True,
            rounds=1, evaluations=1)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        adjoint_estimate)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)
    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)],
        adjoint_options={"base_log2": 6, "max_log2": 12},
        element_absolute_tolerance=0.01, element_relative_tolerance=0.0)

    assert levels_seen == [6, 8, 10]
    assert result["cross_refinement_rounds"] == 2
    assert result["converged"]
    assert result["estimator_consistent"][0]
    assert np.allclose(result["selected_endpoint_mean"][0], [0.025, 0.025])


def test_bidirectional_cross_refines_forward_and_uses_final_endpoint_ensemble(monkeypatch):
    adjoint_replicates = np.array([[0.098], [0.100], [0.100], [0.102]])
    adjoint_endpoints = np.repeat(
        np.array([[[0.04, 0.06]]]), adjoint_replicates.shape[0], axis=0)
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([0.1]), element_stderr=np.array([0.001]),
        element_replicates=adjoint_replicates, log2_samples=np.array([12]),
        total_mean=0.1, total_stderr=0.001, converged=True, rounds=1, evaluations=1,
        auxiliary_mean=adjoint_endpoints.mean(axis=0),
        auxiliary_replicates=adjoint_endpoints)
    levels_seen = []

    def forward_estimate(*args, **kwargs):
        level = int(np.asarray(kwargs.get("initial_log2_samples", [12])).max())
        levels_seen.append(level)
        mean = 0.2 if level == 12 else 0.1
        delta = 0.001 * np.sqrt(3.0)
        replicates = np.array([
            [mean - delta], [mean - delta], [mean + delta], [mean + delta]])
        endpoint = np.array([0.10, 0.10]) if level == 12 else np.array([0.03, 0.07])
        endpoints = np.repeat(endpoint[None, None, :], 4, axis=0)
        return AdaptiveQuadratureResult(
            element_mean=np.array([mean]), element_stderr=np.array([0.001]),
            element_replicates=replicates, log2_samples=np.array([level]),
            total_mean=mean, total_stderr=0.001, converged=True, rounds=1, evaluations=1,
            auxiliary_mean=endpoints.mean(axis=0), auxiliary_replicates=endpoints)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        forward_estimate)
    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], method_hint=np.array(["forward"]),
        adjoint_options={"max_log2": 12},
        forward_options={"base_log2": 12, "max_log2": 14},
        element_absolute_tolerance=0.01, element_relative_tolerance=0.0)

    assert levels_seen == [12, 13]
    assert result["forward_cross_refinement_rounds"] == 1
    assert result["estimator_consistent"][0]
    assert result["method"][0] == "forward"
    assert np.allclose(result["selected_endpoint_mean"][0], [0.03, 0.07])


def test_bidirectional_refines_nonoverlapping_forward_and_adjoint_support(monkeypatch):
    forward = AdaptiveQuadratureResult(
        element_mean=np.array([0.034]), element_stderr=np.array([0.009]),
        element_replicates=np.array([[0.016], [0.034], [0.034], [0.052]]),
        log2_samples=np.array([12]), total_mean=0.034, total_stderr=0.009,
        converged=True, rounds=1, evaluations=1)
    levels_seen = []

    def adjoint_estimate(*args, **kwargs):
        level = int(np.asarray(kwargs.get("initial_log2_samples", [6])).max())
        levels_seen.append(level)
        mean = 0.034 if level >= 10 else 0.0
        stderr = 0.004 if level >= 10 else 1e-31
        delta = stderr * np.sqrt(3.0)
        return AdaptiveQuadratureResult(
            element_mean=np.array([mean]), element_stderr=np.array([stderr]),
            element_replicates=np.array([
                [mean - delta], [mean - delta], [mean + delta], [mean + delta]]),
            log2_samples=np.array([level]), total_mean=mean, total_stderr=stderr,
            converged=True, rounds=1, evaluations=1)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        adjoint_estimate)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: forward)
    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(1.0, 0.0)],
        adjoint_options={"base_log2": 6, "max_log2": 12},
        element_absolute_tolerance=0.01, element_relative_tolerance=0.15,
        consistency_sigma=5.0, support_sigma=2.0, support_ratio=0.5)

    assert levels_seen == [6, 8, 10]
    assert result["cross_refinement_rounds"] == 2
    assert not result["adjoint_support_unresolved"][0]
    assert result["converged"]


def test_bidirectional_support_guard_does_not_reject_modest_multicell_disagreement(monkeypatch):
    def estimate(mean, stderr):
        delta = stderr * np.sqrt(3.0)
        return AdaptiveQuadratureResult(
            element_mean=np.array([mean]), element_stderr=np.array([stderr]),
            element_replicates=np.array([
                [mean - delta], [mean - delta], [mean + delta], [mean + delta]]),
            log2_samples=np.array([14]), total_mean=mean, total_stderr=stderr,
            converged=True, rounds=1, evaluations=1)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: estimate(0.297, 0.0005))
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        lambda *args, **kwargs: estimate(0.364, 0.021))
    result = bidirectional_boundary_state_cell_flux(
        object(), "electron", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0)], [(0.0, -1.0)], element_absolute_tolerance=0.05,
        element_relative_tolerance=0.0, consistency_sigma=5.0,
        support_sigma=2.0, support_ratio=0.5)

    assert result["estimator_discrepancy_sigma"][0] < 5.0
    assert not result["adjoint_support_unresolved"][0]
    assert result["converged"]


def test_bidirectional_refines_forward_uncertainty_after_pooling_cell_faces(monkeypatch):
    adjoint = AdaptiveQuadratureResult(
        element_mean=np.array([0.05, 0.05]),
        element_stderr=np.array([0.001, 0.001]),
        element_replicates=np.array([
            [0.049, 0.049], [0.050, 0.050], [0.050, 0.050], [0.051, 0.051]]),
        log2_samples=np.array([12, 12]), total_mean=0.05, total_stderr=0.001,
        converged=True, rounds=1, evaluations=1)
    levels_seen = []

    def forward_estimate(*args, **kwargs):
        level = int(np.asarray(kwargs.get("initial_log2_samples", [12])).max())
        levels_seen.append(level)
        face_stderr = 0.009 if level == 12 else 0.004
        delta = face_stderr * np.sqrt(3.0)
        replicates = np.array([
            [0.05 - delta, 0.05 - delta], [0.05 - delta, 0.05 - delta],
            [0.05 + delta, 0.05 + delta], [0.05 + delta, 0.05 + delta]])
        return AdaptiveQuadratureResult(
            element_mean=np.array([0.05, 0.05]),
            element_stderr=np.array([face_stderr, face_stderr]),
            element_replicates=replicates, log2_samples=np.array([level, level]),
            total_mean=0.05, total_stderr=face_stderr,
            converged=True, rounds=1, evaluations=1)

    monkeypatch.setattr(
        boundary_transport_module, "adaptive_adjoint_boundary_state_face_flux",
        lambda *args, **kwargs: adjoint)
    monkeypatch.setattr(
        boundary_transport_module, "adaptive_forward_boundary_state_cell_flux",
        forward_estimate)
    result = bidirectional_boundary_state_cell_flux(
        object(), "ion", np.zeros((2, 2)), np.zeros((1, 1), dtype=bool),
        [(0, 0), (0, 0)], [(1.0, 0.0), (-1.0, 0.0)],
        forward_options={"base_log2": 12, "max_log2": 14},
        element_absolute_tolerance=0.005, element_relative_tolerance=0.1)

    assert levels_seen == [12, 13]
    assert result["forward_pool_refinement_rounds"] == 1
    assert result["converged"]


def test_charging_warm_starts_bidirectional_levels_from_accepted_iteration(monkeypatch):
    calls = []

    def hybrid(boundary, species_name, potential, solid, cells, normals, **kwargs):
        calls.append((
            species_name,
            kwargs["adjoint_options"].get("initial_log2_samples"),
            kwargs["forward_options"].get("initial_log2_samples")))
        face_count = len(cells)
        unique_cells = list(dict.fromkeys(cells))
        replicates = np.ones((4, face_count))
        adjoint = AdaptiveQuadratureResult(
            element_mean=np.ones(face_count), element_stderr=np.zeros(face_count),
            element_replicates=replicates, log2_samples=np.full(face_count, 9),
            total_mean=1.0, total_stderr=0.0, converged=True, rounds=1, evaluations=1)
        forward = AdaptiveQuadratureResult(
            element_mean=np.ones(face_count), element_stderr=np.zeros(face_count),
            element_replicates=replicates, log2_samples=np.full(face_count, 11),
            total_mean=1.0, total_stderr=0.0, converged=True, rounds=1, evaluations=1)
        return dict(
            per_face=np.ones(face_count), per_face_stderr=np.zeros(face_count),
            unique_cells=np.asarray(unique_cells),
            method=np.full(len(unique_cells), "adjoint"), converged=True,
            adjoint=adjoint, forward=forward)

    monkeypatch.setattr(
        charging_backward_module, "bidirectional_boundary_state_cell_flux", hybrid)
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    ion = SpeciesBoundaryState(
        "ion", 1, 40.0, 1.0, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    electron = SpeciesBoundaryState(
        "electron", -1, 5.4858e-4, 1.0, [[0.0, 0.0, 1.0]], [1.0],
        density_model=density)
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    solid = np.zeros((4, 4), dtype=bool); solid[:, -1] = True
    result = solve_boundary_state_charging(
        solid, np.zeros_like(solid, dtype=int), boundary,
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


def test_charging_returns_last_evaluated_state_not_an_unassessed_final_step():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    velocity = [[0.0, 0.0, 1.0]]
    boundary = PlasmaBoundaryState((
        SpeciesBoundaryState(
            "ion", 1, 40.0, 2.0, velocity, [1.0], density_model=density),
        SpeciesBoundaryState(
            "electron", -1, 5.4858e-4, 1.0, velocity, [1.0], density_model=density),
    ), reference_plane_m=0.0)
    solid = np.zeros((6, 5), dtype=bool); solid[:, -1] = True
    result = solve_boundary_state_charging(
        solid, np.zeros_like(solid, dtype=int), boundary,
        n_iter=1, min_iter=1, balance_tol=None, beta=0.5,
        response_energy_eV=4.0, field_sweeps=50, trust_region=False)

    assert np.all(result["surface_voltage"][solid] == 0.0)
    assert np.isclose(result["balance_final"]["max_abs_log_ratio"], np.log(2.0))
    assert np.allclose(result["ion_current"], 2.0)
    assert np.allclose(result["electron_current"], 1.0)


def test_adaptive_forward_zero_hits_have_nonzero_confidence_bound():
    nx, nz = 64, 8
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-1e-9, 1e-9]), np.array([-1e-9, 1e-9]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "vertical", 0, 40.0, 1e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    # At 2^4 source samples this particular one-cell target receives no hits for the fixed scrambles.
    result = adaptive_forward_boundary_state_cell_flux(
        boundary, "vertical", np.zeros((nx + 1, nz + 1)), solid, [(0, nz - 1)],
        base_log2=4, max_log2=4, n_replicates=3,
        absolute_tolerance=1.0, relative_tolerance=0.0,
        element_absolute_tolerance=1e-6)
    assert result.element_mean[0] == 0.0
    assert result.element_stderr[0] > 0.0
    assert not result.converged


def test_adaptive_forward_rare_hits_keep_counting_uncertainty_floor():
    nx, nz = 128, 8
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-1e-9, 1e-9]), np.array([-1e-9, 1e-9]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "vertical", 0, 40.0, 1e19, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    result = adaptive_forward_boundary_state_cell_flux(
        boundary, "vertical", np.zeros((nx + 1, nz + 1)), solid, [(0, nz - 1)],
        base_log2=10, max_log2=10, n_replicates=4,
        absolute_tolerance=1.0, relative_tolerance=0.0,
        element_absolute_tolerance=1.0)
    probability = result.element_mean[0] / nx
    binomial_floor = nx * np.sqrt(probability * (1.0 - probability) / (4 * 2 ** 10))
    assert result.element_mean[0] > 0.0
    assert result.element_stderr[0] >= binomial_floor - 1e-15


def test_forward_transport_resolves_exact_oriented_hit_face_without_losing_cell_current():
    nx, nz = 10, 8
    solid = np.zeros((nx, nz), dtype=bool); solid[4, 5] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-1e-9, 1e-9]), np.array([-1e-9, 1e-9]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "vertical", 0, 40.0, 1.0, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    potential = np.zeros((nx + 1, nz + 1))
    cell = forward_boundary_state_cell_flux_qmc(
        boundary, "vertical", potential, solid, [(4, 5)], log2_samples=10, seed=31)
    normals = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    faces = forward_boundary_state_cell_flux_qmc(
        boundary, "vertical", potential, solid, [(4, 5)] * 4,
        normals=normals, log2_samples=10, seed=31)
    assert np.isclose(faces["per_cell"].sum(), cell["per_cell"][0])
    assert np.isclose(faces["per_cell"][2], cell["per_cell"][0])
    assert np.allclose(faces["per_face_endpoint"].sum(axis=1), faces["per_cell"])
    assert np.count_nonzero(faces["per_cell"]) == 1


def test_arbitrary_face_adjoint_has_correct_wall_flux_jacobian():
    nx, nz = 128, 12
    solid = np.zeros((nx, nz), dtype=bool); solid[100, :] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([0.5, 1.5]), np.array([-0.5, 0.5]), np.array([1.0, 2.0])),
        np.ones((1, 1, 1)))
    vx = np.linspace(0.5, 1.5, 64, endpoint=False) + 0.5 / 64
    vz = np.linspace(1.0, 2.0, 128, endpoint=False) + 0.5 / 128
    xx, zz = np.meshgrid(vx, vz, indexing="ij")
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 2e19, np.column_stack((xx.ravel(), np.zeros(xx.size), zz.ravel())),
        np.ones(xx.size), density_model=density)
    # The surface proposal is expressed in local (tangent, out-of-plane, inward-normal) coordinates.
    # For the left-facing wall n=(-1,0), local tangent=-global vz and local normal=global vx.
    local_density = RectilinearVelocityHistogramDensity(
        (np.array([-2.0, -1.0]), np.array([-0.5, 0.5]), np.array([0.5, 1.5])),
        np.ones((1, 1, 1)))
    proposal = SpeciesBoundaryState(
        "wall-proposal", 1, 40.0, 1.0,
        np.column_stack((-zz.ravel(), np.zeros(xx.size), xx.ravel())),
        np.ones(xx.size), density_model=local_density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(100, z) for z in range(2, 9)]
    result = adjoint_boundary_state_face_flux(
        boundary, "ion", np.zeros((nx + 1, nz + 1)), solid, cells,
        [(-1.0, 0.0)] * len(cells), proposal_species=proposal,
        n_face_position=2, max_steps=1000, want_energy=True)
    # The finite tensor midpoint rule must reproduce its own Liouville Jacobian exactly; its continuum
    # limit is E[vx/vz] = 1*ln(2).
    expected = float(np.mean(xx / zz))
    assert np.isclose(result["normalized_flux"], expected, atol=1e-12)
    assert np.allclose(result["per_face_endpoint"].sum(axis=1), result["per_face"], atol=1e-14)
    assert np.isclose(expected, np.log(2.0), rtol=1.1e-3)
    expected_energy = float(np.sum((xx / zz) * (xx * xx + zz * zz)) / np.sum(xx / zz))
    assert np.isclose(result["mean_impact_energy_eV"], expected_energy, atol=1e-12)


def test_self_consistent_charging_consumes_unified_boundary_state_without_source_branches():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.5, 0.5]), np.array([-0.5, 0.5]), np.array([1.0, 2.0])),
        np.ones((1, 1, 1)))
    velocity = np.array([
        [-0.25, 0.0, 1.25], [0.25, 0.0, 1.25],
        [-0.25, 0.0, 1.75], [0.25, 0.0, 1.75],
    ])
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, velocity, np.ones(4), density_model=density)
    electron = SpeciesBoundaryState(
        "e-", -1, 5.4858e-4, 1e19, velocity, np.ones(4), density_model=density)
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    geometry = _build_edge_array_geometry(1.0, W=16, mouth=20)
    result = self_consistent_backward(
        geometry, n_iter=1, n_wall=4, n_floor=3, sweeps=100,
        boundary_state=boundary, ion_species="Ar+", electron_species="e-",
        n_face_position=2)
    # Identical countercharged phase-space measures have identical collisionless trajectories at V=0,
    # so every material capacitor is exactly at its current-balance fixed point.
    assert result["balance_preupdate"]["max_abs_log_ratio"] == 0.0
    assert np.all(result["Vs"] == 0.0)
    assert result["field_final"]["max_abs"] == 0.0


def test_general_charging_solver_uses_only_material_grid_components_and_boundary_state():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.5, 0.5]), np.array([-0.5, 0.5]), np.array([1.0, 2.0])),
        np.ones((1, 1, 1)))
    velocity = np.array([
        [-0.25, 0.0, 1.25], [0.25, 0.0, 1.25],
        [-0.25, 0.0, 1.75], [0.25, 0.0, 1.75],
    ])
    boundary = PlasmaBoundaryState((
        SpeciesBoundaryState("positive_a", 1, 40.0, 0.4e19, velocity, np.ones(4), density_model=density),
        SpeciesBoundaryState("positive_b", 1, 20.0, 0.6e19, velocity, np.ones(4), density_model=density),
        SpeciesBoundaryState("negative", -1, 5.4858e-4, 1e19, velocity, np.ones(4), density_model=density),
    ), reference_plane_m=0.0)
    solid = np.zeros((24, 18), dtype=bool)
    solid[5, 7:] = True; solid[18, 7:] = True; solid[5:19, 15:] = True
    conductor_ids = np.zeros_like(solid, dtype=int)
    conductor_ids[5, 10:14] = 1; conductor_ids[18, 10:14] = 1
    result = solve_boundary_state_charging(
        solid, conductor_ids, boundary, ion_species=("positive_a", "positive_b"),
        electron_species="negative",
        n_iter=2, min_iter=1, n_face_position=2, field_sweeps=100)
    assert result["balance_final"]["max_abs_log_ratio"] <= 3 * np.finfo(float).eps
    assert np.allclose(result["surface_voltage"], 0.0, atol=2e-15)
    assert abs(result["conductor_voltage"][1]) <= 2e-15
    assert set(result["species_current"]) == {"positive_a", "positive_b", "negative"}
    assert result["current_scale_m2_s"] == 1e19


def test_general_charging_rejects_material_on_grounded_plasma_reference_plane():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    velocity = [[0.0, 0.0, 1.0]]
    boundary = PlasmaBoundaryState((
        SpeciesBoundaryState("ion", 1, 40.0, 1.0, velocity, [1.0], density_model=density),
        SpeciesBoundaryState("electron", -1, 5.4858e-4, 1.0, velocity, [1.0], density_model=density),
    ), reference_plane_m=0.0)
    solid = np.zeros((8, 6), dtype=bool); solid[3, 0] = True
    with pytest.raises(ValueError, match="gas-only top row"):
        solve_boundary_state_charging(solid, np.zeros_like(solid, dtype=int), boundary)


def test_charging_trust_region_rolls_back_merit_increase():
    nx, nz = 8, 6
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    conductor_ids = np.zeros_like(solid, dtype=int)
    ion_density = RectilinearVelocityHistogramDensity(
        (np.array([-0.2, 0.2]), np.array([-0.2, 0.2]), np.array([1.0, 2.0])),
        np.ones((1, 1, 1)))
    ion = SpeciesBoundaryState(
        "ion", 1, 40.0, 1.0, [[0.0, 0.0, 1.25], [0.0, 0.0, 1.75]], [1.0, 1.0],
        density_model=ion_density)
    electron = maxwellian_electron_boundary_state(
        1.0, 1.0, n_transverse=3, n_normal=6).get("electron")
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    initial = np.zeros_like(solid, dtype=float); initial[solid] = 1.0
    result = solve_boundary_state_charging(
        solid, conductor_ids, boundary, initial_surface_voltage=initial,
        n_iter=20, min_iter=2, balance_tol=0.1, beta=0.8, dVmax=2.0,
        field_sweeps=100, trust_region=True)
    accepted_merit = np.array([
        item["max_abs_log_ratio"] for item in result["interval_balance_history"]])
    assert result["rejected_steps"] >= 1
    assert np.all(np.diff(accepted_merit) <= 1e-14)
    assert result["interval_balance_final"]["max_abs_log_ratio"] <= 0.1
    accelerated = solve_boundary_state_charging(
        solid, conductor_ids, boundary, initial_surface_voltage=initial,
        n_iter=20, min_iter=2, balance_tol=0.1, beta=0.8, dVmax=2.0,
        field_sweeps=100, trust_region=True, nonlinear_update="anderson")
    assert accelerated["interval_balance_final"]["max_abs_log_ratio"] <= 0.1
    assert accelerated["iterations"] <= result["iterations"]


def _uniform_box_species(name, charge, flux, vx_edges, vz_edges, nx=8, nz=16):
    density = RectilinearVelocityHistogramDensity(
        (np.asarray(vx_edges), np.array([-0.5, 0.5]), np.asarray(vz_edges)),
        np.ones((1, 1, 1)))
    vx = np.linspace(vx_edges[0], vx_edges[1], nx, endpoint=False) + (vx_edges[1] - vx_edges[0]) / (2 * nx)
    vz = np.linspace(vz_edges[0], vz_edges[1], nz, endpoint=False) + (vz_edges[1] - vz_edges[0]) / (2 * nz)
    xx, zz = np.meshgrid(vx, vz, indexing="ij")
    velocity = np.column_stack((xx.ravel(), np.zeros(xx.size), zz.ravel()))
    return SpeciesBoundaryState(name, charge, 40.0, flux, velocity, np.ones(xx.size), density_model=density)


def test_unified_forward_adjoint_reciprocity_in_nonuniform_field_with_separate_proposal():
    nx, nz = 40, 30
    left, right, floor = 10, 30, 25
    solid = np.zeros((nx, nz), dtype=bool)
    solid[left - 1, :floor + 1] = True; solid[right, :floor + 1] = True
    solid[left - 1:right + 1, floor] = True
    target = np.zeros_like(solid); target[left:right, floor] = True
    ii, jj = np.meshgrid(np.arange(nx + 1), np.arange(nz + 1), indexing="ij")
    # Harmonic bilinear potential: both Ex and Ez vary spatially.
    V = 0.015 * jj + 0.0008 * (ii - nx / 2) * jj
    physical = _uniform_box_species("ion", 1, 2e19, (-0.4, 0.4), (5.0, 7.0), nx=12, nz=24)
    # Proposal nodes align both physical histogram edges; misaligned midpoint quadrature converges only
    # first order at the discontinuous support boundary and is tested separately by density gates.
    proposal = _uniform_box_species("proposal", 1, 1.0, (-1.0, 1.0), (4.0, 8.0), nx=20, nz=32)
    boundary = PlasmaBoundaryState((physical,), reference_plane_m=0.0)
    forward = trace_boundary_state_floor_flux(
        boundary, "ion", V, solid, target, n_position=160, max_steps=200 * nz)
    cells = [(x, floor) for x in range(left, right)]
    backward = adjoint_boundary_state_floor_flux(
        boundary, "ion", V, solid, cells, proposal_species=proposal,
        n_face_position=16, max_steps=200 * nz)
    assert np.isclose(backward["normalized_flux"], forward["normalized_flux"], rtol=0.01, atol=0.002), (
        backward["normalized_flux"], forward["normalized_flux"])


def test_bidirectional_nonuniform_field_refuses_state_dependent_timestep():
    nx, nz = 6, 5
    solid = np.zeros((nx, nz), dtype=bool); solid[:, -1] = True
    density = RectilinearVelocityHistogramDensity(
        (np.array([-0.1, 0.1]), np.array([-0.1, 0.1]), np.array([0.9, 1.1])),
        np.ones((1, 1, 1)))
    species = SpeciesBoundaryState(
        "ion", 1, 40.0, 1.0, [[0.0, 0.0, 1.0]], [1.0], density_model=density)
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    potential = np.tile(np.arange(nz + 1), (nx + 1, 1)).astype(float)
    cells = [(index, nz - 1) for index in range(nx)]
    normals = [(0.0, -1.0)] * nx

    with pytest.raises(ValueError, match="reversible adjoint map"):
        bidirectional_boundary_state_cell_flux(
            boundary, "ion", potential, solid, cells, normals,
            adjoint_options=dict(base_log2=4, max_log2=4, n_replicates=4),
            forward_options=dict(base_log2=4, max_log2=4, n_replicates=4))

    result = bidirectional_boundary_state_cell_flux(
        boundary, "ion", potential, solid, cells, normals,
        element_absolute_tolerance=10.0, consistency_sigma=1e6,
        adjoint_options=dict(
            base_log2=4, max_log2=4, n_replicates=4, fixed_dt=0.02),
        forward_options=dict(
            base_log2=4, max_log2=4, n_replicates=4, fixed_dt=0.02))
    assert result["converged"]
