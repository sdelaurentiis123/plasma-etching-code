import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.boundary_transport_3d import trace_boundary_state_field_3d
from petch.charging_poisson import EPS0
from petch.charging_poisson_3d import (
    NodalPoissonSystem3D,
    assemble_q1_stiffness_3d,
    lump_triangle_sheet_charge_3d,
)


def _uniform_surface_charge(cell_shape, sigma, spacing):
    nx, ny, nz = cell_shape
    hx, hy, _ = spacing
    x_weight = np.full(nx + 1, hx); x_weight[[0, -1]] *= 0.5
    y_weight = np.full(ny + 1, hy); y_weight[[0, -1]] *= 0.5
    charge = np.zeros((nx + 1, ny + 1, nz + 1))
    charge[:, :, 0] = sigma * x_weight[:, None] * y_weight[None, :]
    return charge


def test_q1_3d_stiffness_is_symmetric_and_annihilates_constant_voltage():
    stiffness = assemble_q1_stiffness_3d(
        np.ones((2, 3, 4)), spacing_m=(11e-9, 17e-9, 23e-9))

    asymmetry = stiffness - stiffness.T
    assert asymmetry.nnz == 0 or np.max(np.abs(asymmetry.data)) < 1e-20
    assert np.max(np.abs(stiffness @ np.ones(stiffness.shape[0]))) < 1e-20


def test_q1_3d_poisson_refuses_fully_constrained_grid_without_unknowns():
    with pytest.raises(ValueError, match="leave at least one free node"):
        NodalPoissonSystem3D(np.ones((1, 1, 1)), 1e-9, np.ones((2, 2, 2), dtype=bool))


def test_q1_3d_poisson_reproduces_uniform_parallel_plate_capacitance():
    cell_shape = (3, 4, 7); spacing = np.array([13e-9, 19e-9, 23e-9])
    epsilon = np.full(cell_shape, 3.9)
    fixed = np.zeros(tuple(np.asarray(cell_shape) + 1), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, spacing, fixed)
    sigma = 6.2e-4
    voltage, diagnostics = system.solve(_uniform_surface_charge(cell_shape, sigma, spacing))

    expected_surface = sigma * cell_shape[2] * spacing[2] / (EPS0 * 3.9)
    expected_line = np.linspace(expected_surface, 0.0, cell_shape[2] + 1)
    assert np.allclose(voltage, expected_line[None, None, :], rtol=2e-12, atol=2e-12)
    assert diagnostics.max_abs_free_charge_residual_c < 1e-25
    assert abs(diagnostics.charge_balance_c) < 1e-24


def test_q1_3d_poisson_reproduces_series_dielectric_capacitance_and_displacement():
    cell_shape = (2, 3, 8); spacing = np.array([10e-9, 14e-9, 21e-9])
    epsilon = np.ones(cell_shape); epsilon[:, :, 4:] = 4.0
    fixed = np.zeros(tuple(np.asarray(cell_shape) + 1), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, spacing, fixed)
    sigma = 2.4e-4
    voltage, _ = system.solve(_uniform_surface_charge(cell_shape, sigma, spacing))

    expected_surface = sigma * spacing[2] * (4.0 / 1.0 + 4.0 / 4.0) / EPS0
    mean_voltage = voltage.mean(axis=(0, 1))
    displacement = epsilon.mean(axis=(0, 1)) * (mean_voltage[:-1] - mean_voltage[1:]) / spacing[2]
    assert np.allclose(voltage[:, :, 0], expected_surface, rtol=2e-12, atol=2e-12)
    assert np.allclose(displacement, displacement[0], rtol=2e-12, atol=2e-12)


def test_q1_3d_poisson_closes_global_gauss_law():
    epsilon = np.ones((3, 2, 4)); epsilon[:, :, 2:] = 3.9
    fixed = np.zeros((4, 3, 5), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, 18e-9, fixed)
    charge = np.zeros((4, 3, 5)); charge[1, 1, 1] = 2.5e-18; charge[2, 1, 3] = -0.4e-18
    _, diagnostics = system.solve(charge)

    assert np.isclose(diagnostics.specified_charge_c, charge.sum())
    assert np.isclose(
        diagnostics.dirichlet_reaction_charge_c, -charge.sum(), rtol=2e-12, atol=1e-29)
    assert abs(diagnostics.charge_balance_c) < 1e-28


def test_triangle_sheet_projection_conserves_charge_and_first_moment():
    vertices = np.array([
        [0.10, 0.20, 0.30],
        [0.80, 0.25, 0.40],
        [0.20, 0.90, 0.60],
    ])
    faces = np.array([[0, 1, 2]]); sigma = np.array([3.7e-4]); unit_m = 2.5e-9
    charge = lump_triangle_sheet_charge_3d(
        (3, 3, 3), vertices, faces, sigma, grid_spacing=1.0,
        coordinate_length_unit_m=unit_m)

    edge_a = (vertices[1] - vertices[0]) * unit_m
    edge_b = (vertices[2] - vertices[0]) * unit_m
    expected_charge = sigma[0] * 0.5 * np.linalg.norm(np.cross(edge_a, edge_b))
    grid = np.stack(np.meshgrid(np.arange(3), np.arange(3), np.arange(3), indexing="ij"))
    deposited_centroid = np.array([
        np.sum(charge * grid[axis]) / charge.sum() for axis in range(3)])

    assert np.isclose(charge.sum(), expected_charge, rtol=1e-14)
    assert np.allclose(deposited_centroid, vertices.mean(axis=0), rtol=1e-13, atol=1e-14)


def test_triangle_sheet_projection_accepts_float32_vertices_on_grid_endpoint():
    # Marching cubes emits float32 coordinates.  The endpoint rounds upward before division by the
    # float64 grid spacing: float32(0.3) / 0.01 > 30, although the intended vertex is on the grid.
    vertices = np.array([
        [0.30, 0.10, 0.10],
        [0.30, 0.11, 0.10],
        [0.30, 0.10, 0.11],
    ], dtype=np.float32)
    faces = np.array([[0, 1, 2]])

    charge = lump_triangle_sheet_charge_3d(
        (31, 31, 31), vertices, faces, np.array([1.0]), grid_spacing=0.01)

    assert np.isclose(charge.sum(), 0.00005, rtol=2e-6)
    assert np.count_nonzero(charge[-1]) > 0


def test_q1_3d_diagonal_capacitance_is_positive_and_symmetric():
    epsilon = np.full((3, 3, 4), 3.9)
    fixed = np.zeros((4, 4, 5), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, 20e-9, fixed)
    capacitance = system.diagonal_capacitance(np.array([[1, 1, 0], [2, 2, 0]]))

    assert np.all(capacitance > 0.0)
    assert np.allclose(capacitance[0], capacitance[1], rtol=1e-12)


def test_support_response_inverse_reproduces_requested_voltage_step_exactly():
    epsilon = np.ones((3, 3, 5))
    fixed = np.zeros((4, 4, 6), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, 20e-9, fixed)
    nodes = np.array([[1, 1, 1], [2, 1, 1], [1, 2, 2], [2, 2, 2]])
    requested = np.array([-3.0, 1.5, 2.0, -0.75])
    response = system.voltage_response(nodes)
    charge = np.zeros(system.shape)
    charge[tuple(nodes.T)] = np.linalg.solve(response, requested)

    voltage, _ = system.solve(charge)

    assert np.allclose(voltage[tuple(nodes.T)], requested, rtol=3e-13, atol=3e-13)


def test_solved_dielectric_sheet_potential_decelerates_ion_by_electrostatic_work():
    cell_shape = (1, 1, 10); spacing_m = np.full(3, 0.1e-6)
    fixed = np.zeros((2, 2, 11), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(np.ones(cell_shape), spacing_m, fixed)
    target_surface_voltage = 10.0
    sigma = target_surface_voltage * EPS0 / (cell_shape[2] * spacing_m[2])
    potential, _ = system.solve(_uniform_surface_charge(cell_shape, sigma, spacing_m))

    vertices = np.array([
        [0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
        [0.1, 0.1, 0.0], [0.0, 0.1, 0.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    areas = np.full(2, 0.005)
    ion = SpeciesBoundaryState(
        "Ar+", 1, 40.0, 1e19, [[0.0, 0.0, np.sqrt(20.0)]], [1.0])
    boundary = PlasmaBoundaryState((ion,), reference_plane_m=1e-6)
    transport = trace_boundary_state_field_3d(
        boundary, {"Ar+": "energetic_bombardment"}, vertices, faces, areas,
        source_bounds=(0.0, 0.1, 0.0, 0.1), source_z=1.0,
        nodal_potential_v=potential, potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.1, mesh_length_unit_m=1e-6,
        n_position=16, seed=41, fixed_dt=0.0025, max_steps=2000, device="cpu")

    impact_energy = transport.surface_fluxes.energetic_fluxes[0].event_energy_eV
    assert np.allclose(potential[:, :, 0], target_surface_voltage, rtol=2e-12)
    assert np.allclose(impact_energy, 10.0, atol=3e-3)
