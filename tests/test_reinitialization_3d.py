import numpy as np

from petch.threed import reinit_cr2


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

