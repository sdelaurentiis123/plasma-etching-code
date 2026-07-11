import numpy as np
import pytest

from petch.boundary_state import (
    PlasmaBoundaryState,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
    instantaneous_sinusoidal_ion_boundary_state,
)
from petch.boundary_transport import (
    adjoint_boundary_state_floor_flux,
    boundary_launches_2d,
    trace_boundary_state_floor_flux,
)


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
