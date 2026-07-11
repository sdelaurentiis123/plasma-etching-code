import numpy as np
import pytest

from petch.boundary_state import (
    PlasmaBoundaryState,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
    instantaneous_sinusoidal_ion_boundary_state,
)
from petch.boundary_transport import (
    adjoint_boundary_state_face_flux,
    adjoint_boundary_state_floor_flux,
    boundary_launches_2d,
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
    boundary = PlasmaBoundaryState((species,), reference_plane_m=0.0)
    cells = [(100, z) for z in range(2, 9)]
    result = adjoint_boundary_state_face_flux(
        boundary, "ion", np.zeros((nx + 1, nz + 1)), solid, cells,
        [(-1.0, 0.0)] * len(cells), n_face_position=2, max_steps=1000, want_energy=True)
    # The finite tensor midpoint rule must reproduce its own Liouville Jacobian exactly; its continuum
    # limit is E[vx/vz] = 1*ln(2).
    expected = float(np.mean(xx / zz))
    assert np.isclose(result["normalized_flux"], expected, atol=1e-12)
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
    RectilinearVelocityHistogramDensity,
