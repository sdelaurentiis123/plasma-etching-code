import numpy as np

from petch.feature_step_3d import make_rectangular_trench_geometry_3d
from petch.threed import advect_3d, reinit_cr2


def _cut_edge_displacements(before, after, dx):
    displacement = []
    for axis in range(before.ndim):
        left = [slice(None)] * before.ndim; right = [slice(None)] * before.ndim
        left[axis] = slice(None, -1); right[axis] = slice(1, None)
        left = tuple(left); right = tuple(right)
        cut = (before[left] >= 0.0) != (before[right] >= 0.0)
        before_fraction = before[left][cut] / (before[left][cut] - before[right][cut])
        after_fraction = after[left][cut] / (after[left][cut] - after[right][cut])
        displacement.append(np.abs(after_fraction - before_fraction) * dx)
    return np.concatenate(displacement)


def test_cr2_anchors_distorted_slanted_interface_and_restores_gradient():
    dx = 0.1
    coordinate = np.arange(13) * dx
    x, y, z = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
    exact_distance = (0.63719 + 0.15 * x - 0.10 * y - z) / np.sqrt(
        1.0 + 0.15**2 + 0.10**2)
    distorted = exact_distance * (0.55 + 0.25 * np.cos(2.0 * x) * np.cos(1.5 * y))

    redistanced = reinit_cr2(distorted, dx, 4.0 * dx)

    displacement = _cut_edge_displacements(distorted, redistanced, dx)
    gradient = np.gradient(redistanced, dx)
    gradient_norm = np.sqrt(sum(component * component for component in gradient))
    near_interface = np.abs(redistanced) < 2.0 * dx
    assert np.array_equal(redistanced >= 0.0, distorted >= 0.0)
    assert np.max(displacement) < 0.08 * dx
    assert np.mean(displacement) < 0.02 * dx
    assert abs(np.mean(gradient_norm[near_interface]) - 1.0) < 0.08


def _vertical_crossing(phi, dx, i, j):
    line = phi[i, j]
    crossing = np.flatnonzero((line[:-1] >= 0.0) & (line[1:] < 0.0))
    assert crossing.size == 1
    lower = int(crossing[0])
    fraction = line[lower] / (line[lower] - line[lower + 1])
    return dx * (lower + fraction)


def test_cr2_accumulates_repeated_subcell_motion_instead_of_erasing_it():
    dx = 0.1
    coordinate = np.arange(13) * dx
    x, y, z = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
    normalizer = np.sqrt(1.0 + 0.15**2 + 0.10**2)
    phi = (0.63719 + 0.15 * x - 0.10 * y - z) / normalizer
    phi = reinit_cr2(phi, dx, 4.0 * dx)
    initial = _vertical_crossing(phi, dx, 6, 6)
    shift = 0.01 * dx

    for _ in range(40):
        phi = reinit_cr2(phi - shift, dx, 4.0 * dx)

    final = _vertical_crossing(phi, dx, 6, 6)
    expected = 40.0 * shift * normalizer
    assert np.isclose(initial - final, expected, rtol=0.03, atol=0.002 * dx)


def test_cr2_accumulates_floor_motion_beneath_a_pinned_narrow_mask():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.5, cell_length=0.1, domain_height=2.35, dx=0.02,
        opening_width=0.08, mask_thickness=0.7,
        substrate_top=1.4, etched_depth=0.0)
    phi = geometry.phi.copy()
    pinned_mask = geometry.material_id == 2
    speed = np.full(phi.shape, 1.0e-4)
    initial = _vertical_crossing(
        phi, geometry.dx, phi.shape[0] // 2, phi.shape[1] // 2)

    duration_per_step = 12.5
    unreinitialized = geometry.phi.copy()
    for _ in range(16):
        unreinitialized = advect_3d(
            unreinitialized, speed, geometry.dx, duration_per_step)
        unreinitialized[pinned_mask] = geometry.phi[pinned_mask]
        phi = advect_3d(phi, speed, geometry.dx, duration_per_step)
        phi[pinned_mask] = geometry.phi[pinned_mask]
        phi = reinit_cr2(phi, geometry.dx, 4.0 * geometry.dx)
        phi[pinned_mask] = geometry.phi[pinned_mask]

    final = _vertical_crossing(
        phi, geometry.dx, phi.shape[0] // 2, phi.shape[1] // 2)
    reference = _vertical_crossing(
        unreinitialized, geometry.dx,
        unreinitialized.shape[0] // 2, unreinitialized.shape[1] // 2)
    expected = 16.0 * duration_per_step * speed.flat[0]
    assert initial - final > 0.8 * expected
    assert np.isclose(initial - final, initial - reference, rtol=0.03, atol=0.002 * geometry.dx)
