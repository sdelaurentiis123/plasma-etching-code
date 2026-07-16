import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.boundary_transport_3d import trace_boundary_state_field_3d
from petch.charging_poisson import EPS0
from petch.charging_poisson_3d import (
    CompatibleQ1SurfaceChargeProjector3D,
    NodalPoissonSystem3D,
    assemble_q1_stiffness_3d,
    lump_mixed_surface_density_3d,
    lump_triangle_sheet_charge_3d,
    triangle_sheet_face_charge_coupling_3d,
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


def test_floating_conductor_projection_preserves_component_charge_and_is_equipotential():
    epsilon = np.ones((3, 2, 4))
    fixed = np.zeros((4, 3, 5), dtype=bool); fixed[:, :, -1] = True
    conductor = np.zeros(fixed.shape, dtype=int)
    conductor[0:2, :, 0:3] = 7
    system = NodalPoissonSystem3D(
        epsilon, 20e-9, fixed, floating_conductor_node_ids=conductor)
    charge = np.zeros(system.shape)
    charge[0, 0, 0] = 3.0e-18
    charge[1, 2, 2] = -0.7e-18
    charge[3, 1, 1] = 0.4e-18

    canonical = system.canonicalize_charge(charge)
    voltage, diagnostics = system.solve(canonical)
    selected = conductor == 7

    assert system.has_floating_conductors
    assert system.floating_conductor_ids == (7,)
    assert np.isclose(canonical[selected].sum(), charge[selected].sum(), atol=3e-31)
    assert np.array_equal(canonical[conductor == 0], charge[conductor == 0])
    assert np.ptp(voltage[selected]) < 2e-12
    assert diagnostics.maximum_floating_conductor_voltage_spread_v < 2e-12
    assert diagnostics.floating_conductor_ids == (7,)
    assert np.isclose(
        diagnostics.floating_conductor_charge_c[0], charge[selected].sum(), atol=3e-31)
    assert np.isclose(canonical.sum(), charge.sum(), atol=3e-31)


def test_distinct_floating_conductors_pool_independently():
    epsilon = np.ones((4, 1, 3))
    fixed = np.zeros((5, 2, 4), dtype=bool); fixed[:, :, -1] = True
    conductor = np.zeros(fixed.shape, dtype=int)
    conductor[0:2, :, 0:2] = 2
    conductor[3:5, :, 0:2] = 9
    system = NodalPoissonSystem3D(
        epsilon, 25e-9, fixed, floating_conductor_node_ids=conductor)
    charge = np.zeros(system.shape)
    charge[0, 0, 0] = -2e-18
    charge[4, 1, 0] = 4e-18

    voltage, diagnostics = system.solve(charge)

    assert np.ptp(voltage[conductor == 2]) < 2e-12
    assert np.ptp(voltage[conductor == 9]) < 2e-12
    assert not np.isclose(
        diagnostics.floating_conductor_voltage_v[0],
        diagnostics.floating_conductor_voltage_v[1], rtol=0.0, atol=1e-6)
    assert np.allclose(
        diagnostics.floating_conductor_charge_c, (-2e-18, 4e-18),
        rtol=0.0, atol=3e-31)


def test_mixed_surface_coupling_pools_conductor_faces_without_losing_charge():
    vertices = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    fixed = np.zeros((2, 2, 3), dtype=bool); fixed[:, :, -1] = True
    conductor = np.zeros(fixed.shape, dtype=int); conductor[:, :, 0] = 4
    system = NodalPoissonSystem3D(
        np.ones((1, 1, 2)), 1.0, fixed,
        floating_conductor_node_ids=conductor)
    density = np.array([3.0e-18, -0.5e-18])

    coupled = lump_mixed_surface_density_3d(
        system, vertices, faces, density, np.array([4, 4]))
    voltage, diagnostics = system.solve(coupled)

    expected = 0.5 * density.sum()
    assert np.isclose(coupled.sum(), expected, rtol=0.0, atol=2e-31)
    assert np.isclose(
        diagnostics.floating_conductor_charge_c[0], expected,
        rtol=0.0, atol=2e-31)
    assert np.ptp(voltage[conductor == 4]) < 2e-12

    physical = lump_mixed_surface_density_3d(
        system, vertices, faces, density, np.full(2, 4),
        canonicalize=False)
    assert np.min(physical) >= 0.0
    assert np.isclose(physical.sum(), expected, rtol=0.0, atol=2e-31)
    assert np.count_nonzero(physical) == 1


def test_mixed_compatible_projector_preserves_field_and_each_conductor_inventory():
    vertices = np.array([
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0], [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0], [2.0, 1.0, 0.0],
    ])
    faces = np.array([
        [0, 2, 3], [0, 3, 1],
        [2, 4, 5], [2, 5, 3],
    ])
    fixed = np.zeros((3, 2, 3), dtype=bool)
    fixed[:, :, -1] = True
    conductor = np.zeros(fixed.shape, dtype=int)
    conductor[:2, :, :2] = 4
    system = NodalPoissonSystem3D(
        np.ones((2, 1, 2)), 1.0, fixed,
        floating_conductor_node_ids=conductor)
    face_conductor_id = np.array([4, 4, 0, 0])
    projector = CompatibleQ1SurfaceChargeProjector3D.from_mixed_poisson_system(
        system, vertices, faces, face_conductor_id)
    face_charge = np.array([4.0e-18, -1.0e-18, 2.5e-18, -0.7e-18])
    compatible = projector.project_face_charge(face_charge)
    area = projector.physical_face_area_m2
    before = lump_mixed_surface_density_3d(
        system, vertices, faces, face_charge / area, face_conductor_id)
    after = lump_mixed_surface_density_3d(
        system, vertices, faces, compatible / area, face_conductor_id)
    potential_before, diagnostics_before = system.solve(before)
    potential_after, diagnostics_after = system.solve(after)

    assert projector.nullity >= 1
    assert projector.unresolved_fraction(face_charge) > 0.1
    assert np.allclose(after, before, rtol=0.0, atol=5e-32)
    assert np.allclose(
        potential_after, potential_before, rtol=2e-14, atol=2e-20)
    assert np.isclose(
        compatible[face_conductor_id == 4].sum(),
        face_charge[face_conductor_id == 4].sum(),
        rtol=0.0, atol=5e-32)
    assert np.isclose(compatible.sum(), face_charge.sum(), rtol=0.0, atol=5e-32)
    conductor_density = (
        compatible[face_conductor_id == 4] / area[face_conductor_id == 4])
    assert np.allclose(
        conductor_density, np.mean(conductor_density),
        rtol=2e-14, atol=2e-31)
    assert np.allclose(
        diagnostics_after.floating_conductor_charge_c,
        diagnostics_before.floating_conductor_charge_c,
        rtol=0.0, atol=5e-32)
    assert np.allclose(
        projector.project_face_charge(compatible), compatible,
        rtol=2e-14, atol=5e-32)
    assert projector.unresolved_linear_functional_fraction(
        np.ones(len(faces))) < 3e-14
    one_conductor_face = np.zeros(len(faces))
    one_conductor_face[0] = 1.0
    assert projector.unresolved_linear_functional_fraction(
        one_conductor_face) > 0.1


def test_periodic_q1_poisson_identifies_endpoint_voltage_and_charge_exactly():
    epsilon = np.ones((4, 3, 5))
    fixed = np.zeros((5, 4, 6), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(
        epsilon, (11e-9, 17e-9, 23e-9), fixed, periodic_axes=(0, 1))
    charge = np.random.default_rng(431).normal(scale=2e-18, size=system.shape)
    charge[:, :, -1] = 0.0
    canonical = system.canonicalize_charge(charge)
    voltage, diagnostics = system.solve(charge)

    assert system.reduced_shape == (4, 3, 6)
    assert np.array_equal(voltage[0, :, :], voltage[-1, :, :])
    assert np.array_equal(voltage[:, 0, :], voltage[:, -1, :])
    assert np.isclose(canonical.sum(), charge.sum(), rtol=0.0, atol=2e-31)
    assert np.allclose(
        system.reduce_charge(canonical), system.reduce_charge(charge),
        rtol=0.0, atol=2e-31)
    assert abs(diagnostics.charge_balance_c) < 2e-27


def test_periodic_q1_poisson_reproduces_uniform_parallel_plate_exactly():
    cell_shape = (3, 4, 7)
    spacing = np.array([13e-9, 19e-9, 23e-9])
    epsilon = np.full(cell_shape, 3.9)
    fixed = np.zeros(tuple(np.asarray(cell_shape) + 1), dtype=bool)
    fixed[:, :, -1] = True
    periodic = NodalPoissonSystem3D(
        epsilon, spacing, fixed, periodic_axes=(0, 1))
    nonperiodic = NodalPoissonSystem3D(epsilon, spacing, fixed)
    sigma = 6.2e-4
    charge = _uniform_surface_charge(cell_shape, sigma, spacing)

    periodic_voltage, periodic_diagnostics = periodic.solve(charge)
    nonperiodic_voltage, _ = nonperiodic.solve(charge)
    expected_surface = sigma * cell_shape[2] * spacing[2] / (EPS0 * 3.9)
    expected_line = np.linspace(expected_surface, 0.0, cell_shape[2] + 1)

    assert np.allclose(
        periodic_voltage, expected_line[None, None, :], rtol=2e-12, atol=2e-12)
    assert np.allclose(periodic_voltage, nonperiodic_voltage, rtol=2e-12, atol=2e-12)
    assert periodic_diagnostics.max_abs_free_charge_residual_c < 1e-25
    assert abs(periodic_diagnostics.charge_balance_c) < 1e-24


def test_periodic_q1_poisson_is_invariant_to_endpoint_representative():
    epsilon = np.ones((3, 2, 4))
    fixed = np.zeros((4, 3, 5), dtype=bool); fixed[:, :, -1] = True
    system = NodalPoissonSystem3D(epsilon, 19e-9, fixed, periodic_axes=(0, 1))
    left = np.zeros(system.shape); left[0, 1, 2] = 7e-18
    right = np.zeros(system.shape); right[-1, 1, 2] = 7e-18

    voltage_left, _ = system.solve(left)
    voltage_right, _ = system.solve(right)

    assert np.array_equal(system.reduce_charge(left), system.reduce_charge(right))
    assert np.array_equal(voltage_left, voltage_right)


def test_periodic_q1_poisson_refuses_invalid_axes_and_conflicting_voltages():
    epsilon = np.ones((2, 2, 2))
    fixed = np.zeros((3, 3, 3), dtype=bool); fixed[:, :, -1] = True
    with pytest.raises(ValueError, match="integer axes"):
        NodalPoissonSystem3D(epsilon, 1.0, fixed, periodic_axes=(0.5,))
    voltage = np.zeros_like(epsilon, shape=fixed.shape)
    voltage[0, 0, -1] = 1.0
    with pytest.raises(ValueError, match="identified Dirichlet voltages"):
        NodalPoissonSystem3D(
            epsilon, 1.0, fixed, voltage, periodic_axes=(0, 1))

    partial = np.zeros((3, 3, 3), dtype=bool)
    partial[0, :, -1] = True
    with pytest.raises(ValueError, match="Dirichlet masks must agree"):
        NodalPoissonSystem3D(
            epsilon, 1.0, partial, periodic_axes=(0, 1))


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


def test_triangle_sheet_projection_retains_float32_endpoint_tolerance_after_driver_cast():
    vertices = np.array([
        [0.30, 0.10, 0.10],
        [0.30, 0.11, 0.10],
        [0.30, 0.10, 0.11],
    ], dtype=np.float32).astype(np.float64)
    faces = np.array([[0, 1, 2]])

    charge = lump_triangle_sheet_charge_3d(
        (31, 31, 31), vertices, faces, np.array([1.0]), grid_spacing=0.01)
    coupling, _ = triangle_sheet_face_charge_coupling_3d(
        (31, 31, 31), vertices, faces, grid_spacing=0.01)

    assert np.isclose(charge.sum(), 0.00005, rtol=2e-6)
    assert np.isclose(np.asarray(coupling.sum(axis=0)).item(), 1.0, atol=2e-14)


def test_triangle_sheet_projection_still_refuses_a_resolved_out_of_grid_vertex():
    vertices = np.array([
        [0.3001, 0.10, 0.10],
        [0.3001, 0.11, 0.10],
        [0.3001, 0.10, 0.11],
    ])
    with pytest.raises(ValueError, match="outside the nodal grid"):
        lump_triangle_sheet_charge_3d(
            (31, 31, 31), vertices, np.array([[0, 1, 2]]),
            np.array([1.0]), grid_spacing=0.01)


def _triangulated_plane(nx=3, ny=3):
    vertices = np.asarray([
        [float(i), float(j), 0.0]
        for i in range(nx + 1) for j in range(ny + 1)])

    def vertex(i, j):
        return i * (ny + 1) + j

    faces = []
    for i in range(nx):
        for j in range(ny):
            lower_left = vertex(i, j)
            lower_right = vertex(i + 1, j)
            upper_right = vertex(i + 1, j + 1)
            upper_left = vertex(i, j + 1)
            faces.extend((
                [lower_left, lower_right, upper_right],
                [lower_left, upper_right, upper_left]))
    return vertices, np.asarray(faces, dtype=int)


def test_triangle_face_charge_coupling_matches_direct_q1_projection():
    vertices, faces = _triangulated_plane()
    coupling, physical_area = triangle_sheet_face_charge_coupling_3d(
        (4, 4, 2), vertices, faces)
    face_charge = np.random.default_rng(812).normal(scale=2e-18, size=len(faces))
    direct = lump_triangle_sheet_charge_3d(
        (4, 4, 2), vertices, faces, face_charge / physical_area)

    assert np.allclose(
        np.asarray(coupling @ face_charge).reshape(direct.shape), direct,
        rtol=2e-15, atol=2e-33)
    assert np.allclose(np.asarray(coupling.sum(axis=0)).ravel(), 1.0, atol=2e-14)
    assert np.min(coupling.data) >= 0.0


def test_compatible_q1_surface_charge_removes_only_field_invisible_modes():
    vertices, faces = _triangulated_plane()
    projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        (4, 4, 2), vertices, faces)
    face_charge = np.random.default_rng(913).normal(scale=2e-18, size=len(faces))
    compatible = projector.project_face_charge(face_charge)
    node_before = projector.node_charge_from_face_charge(face_charge)
    node_after = projector.node_charge_from_face_charge(compatible)

    fixed = np.zeros((4, 4, 2), dtype=bool)
    fixed[:, :, -1] = True
    poisson = NodalPoissonSystem3D(np.ones((3, 3, 1)), 1.0, fixed)
    potential_before, _ = poisson.solve(node_before)
    potential_after, _ = poisson.solve(node_after)

    assert projector.nullity == 4
    assert projector.rank == 14
    assert projector.unresolved_fraction(face_charge) > 0.1
    assert np.allclose(node_after, node_before, rtol=0.0, atol=6e-33)
    assert np.allclose(potential_after, potential_before, rtol=2e-14, atol=2e-21)
    assert np.isclose(compatible.sum(), face_charge.sum(), rtol=0.0, atol=6e-33)
    assert np.allclose(
        projector.project_face_charge(compatible), compatible,
        rtol=2e-15, atol=6e-33)
    assert projector.unresolved_fraction(compatible) < 3e-15


def test_compatible_surface_projector_uses_periodic_poisson_charge_space():
    vertices, faces = _triangulated_plane()
    fixed = np.zeros((4, 4, 2), dtype=bool); fixed[:, :, -1] = True
    poisson = NodalPoissonSystem3D(
        np.ones((3, 3, 1)), 1.0, fixed, periodic_axes=(0, 1))
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson, vertices, faces)
    face_charge = np.random.default_rng(915).normal(scale=2e-18, size=len(faces))
    compatible = projector.project_face_charge(face_charge)
    reduced_before = projector.node_charge_from_face_charge(face_charge)
    reduced_after = projector.node_charge_from_face_charge(compatible)
    full_before, _ = triangle_sheet_face_charge_coupling_3d(
        poisson.shape, vertices, faces)
    full_before = np.asarray(full_before @ face_charge).reshape(poisson.shape)

    assert projector._node_shape == poisson.reduced_shape
    assert np.allclose(
        reduced_before, poisson.reduce_charge(full_before), rtol=0.0, atol=5e-33)
    assert np.allclose(reduced_after, reduced_before, rtol=0.0, atol=5e-33)
    assert projector.nullity > 4
    assert projector.unresolved_fraction(compatible) < 4e-15
    potential_before, _ = poisson.solve(full_before)
    potential_after, _ = poisson.solve(
        poisson.canonicalize_reduced_charge(reduced_after))
    assert np.allclose(potential_after, potential_before, rtol=2e-14, atol=2e-21)


def test_compatible_q1_surface_charge_reconstructs_in_range_nodal_load():
    vertices, faces = _triangulated_plane()
    projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        (4, 4, 2), vertices, faces)
    original = np.random.default_rng(914).normal(scale=3e-18, size=len(faces))
    node = projector.node_charge_from_face_charge(original)
    reconstructed = projector.face_charge_from_node_charge(node)

    assert np.allclose(
        projector.node_charge_from_face_charge(reconstructed), node,
        rtol=0.0, atol=4e-33)
    assert np.isclose(reconstructed.sum(), original.sum(), rtol=0.0, atol=4e-33)
    assert projector.unresolved_fraction(reconstructed) < 3e-15


def test_compatible_q1_projector_identifies_statistics_that_see_null_modes():
    vertices, faces = _triangulated_plane()
    projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        (4, 4, 2), vertices, faces)
    global_charge = np.ones(len(faces))
    one_triangle = np.zeros(len(faces)); one_triangle[0] = 1.0

    assert projector.unresolved_linear_functional_fraction(global_charge) < 3e-15
    assert projector.unresolved_linear_functional_fraction(one_triangle) > 0.1
    with pytest.raises(ValueError, match="one finite value"):
        projector.unresolved_linear_functional_fraction(np.ones(len(faces) - 1))


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
