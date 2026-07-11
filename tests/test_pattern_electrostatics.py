"""Frozen-field prerequisites for the Huang et al. pattern-distortion problem.

These gates validate only the nonlocal electrostatic kernel. They do not claim stochastic charging or
profile twisting, which require self-consistent arrival histories and 3-D evolving geometry.
"""
import numpy as np

from petch.charging_general import GROUND, INSULATOR, poisson_field


def _charged_pattern(nx, left_charge, right_charge, sweeps=1500):
    nz = 61; center = nx // 2
    material = np.zeros((nx, nz), dtype=np.int8)
    material[:, 45:60] = INSULATOR
    material[:, -1] = GROUND
    sigma = np.zeros((nx, nz))
    sigma[center - 35:center - 24, 45] = left_charge
    sigma[center + 25:center + 36, 45] = right_charge
    potential = poisson_field(
        material, sigma, cell_size_m=2e-9, eps_insulator=3.9,
        sweeps=sweeps, omega=1.7)
    lateral_field = -np.gradient(potential, axis=0)
    return potential, lateral_field, center


def test_symmetric_pattern_has_no_preferred_lateral_field():
    potential, lateral_field, center = _charged_pattern(121, 1e-5, 1e-5)
    assert np.array_equal(potential, potential[::-1])
    assert lateral_field[center, 20] == 0.0


def test_dense_to_sparse_pattern_field_direction_and_padding_converge():
    # Huang et al. (JVST A 38, 023001): the more-positive dense region is on the left and the
    # less-positive sparse region on the right, so Ex must point left-to-right.
    fields = []
    for nx in (121, 161):
        _, lateral_field, center = _charged_pattern(nx, 1.5e-5, 0.5e-5)
        fields.append(float(lateral_field[center, 20]))
    assert fields[0] > 0.0 and fields[1] > 0.0
    assert np.isclose(fields[0], fields[1], rtol=0.05)

