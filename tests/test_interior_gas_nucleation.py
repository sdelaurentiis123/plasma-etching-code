"""Gates for the interior-gas-nucleation guard (material-seam tear fix).

The physical invariant: gas cannot nucleate in the interior of solid — every
legitimate new gas cell (including cavity pinch-off) is 6-adjacent to
pre-step gas. Cells violating it are restored to their pre-step values.
"""

import numpy as np

from petch.feature_step_3d import _suppress_interior_gas_nucleation


def _solid_block_with_top_gas(shape=(6, 6, 8), gas_above=5):
    phi = np.ones(shape)
    phi[:, :, gas_above:] = -1.0
    return phi


def test_isolated_interior_gas_is_restored_with_count():
    previous = _solid_block_with_top_gas()
    new = previous.copy()
    new[3, 3, 1] = -0.25  # deep interior, nowhere near the top gas
    fixed, _, count = _suppress_interior_gas_nucleation(
        previous, None, new, None)
    assert count == 1
    assert fixed[3, 3, 1] == previous[3, 3, 1]


def test_gas_adjacent_to_prior_gas_is_preserved():
    previous = _solid_block_with_top_gas()
    new = previous.copy()
    new[3, 3, 4] = -0.1  # directly below the old gas front: legitimate etch
    fixed, _, count = _suppress_interior_gas_nucleation(
        previous, None, new, None)
    assert count == 0
    assert fixed[3, 3, 4] == -0.1


def test_seam_sheet_between_two_materials_is_restored():
    """The exact failure mode: a one-cell gas sheet at an on-grid seam."""
    shape = (6, 6, 10)
    previous = np.ones(shape)
    previous[:, :, 8:] = -1.0
    seam_z = 4
    m_lower = np.ones(shape)
    m_lower[:, :, seam_z:] = -1.0
    m_lower[:, :, seam_z] = 0.0
    m_upper = -np.ones(shape)
    m_upper[:, :, seam_z:8] = 1.0
    m_upper[:, :, seam_z] = 0.0
    m_upper[:, :, 8:] = -1.0
    new_lower = m_lower.copy()
    new_lower[:, :, seam_z] = -1e-12   # the signed-zero drift
    new_upper = m_upper.copy()
    new_upper[:, :, seam_z] = -1e-12
    new_phi = np.maximum(new_lower, new_upper)
    fixed_phi, fixed_materials, count = _suppress_interior_gas_nucleation(
        previous, {1: m_lower, 2: m_upper}, new_phi,
        {1: new_lower, 2: new_upper})
    assert count == shape[0] * shape[1]
    assert np.all(fixed_phi[:, :, seam_z] == previous[:, :, seam_z])
    assert np.all(fixed_materials[1][:, :, seam_z] == m_lower[:, :, seam_z])
    assert np.all(fixed_materials[2][:, :, seam_z] == m_upper[:, :, seam_z])


def test_periodic_wrap_adjacency_counts():
    """Old gas at the far periodic plane is a neighbor across the seam."""
    shape = (6, 6, 8)
    previous = np.ones(shape)
    previous[4, :, :] = -1.0        # a gas column in the core
    previous[-1, :, :] = previous[0, :, :]  # duplicate-endpoint convention
    new = previous.copy()
    new[3, 2, 3] = -0.2             # adjacent to the gas column: legitimate
    _, _, count = _suppress_interior_gas_nucleation(
        previous, None, new, None, periodic_lateral=True)
    assert count == 0
    # Same cell moved away from any gas: suppressed.
    new2 = previous.copy()
    new2[1, 2, 3] = -0.2
    _, _, count2 = _suppress_interior_gas_nucleation(
        previous, None, new2, None, periodic_lateral=True)
    assert count2 == 1


def test_cavity_pinch_off_untouched():
    """An enclosed cavity keeps prior-gas adjacency; the guard must not
    interfere with legitimate topology events."""
    shape = (6, 6, 10)
    previous = np.ones(shape)
    previous[2:4, 2:4, 2:6] = -1.0  # open gas finger
    new = previous.copy()
    new[2:4, 2:4, 6] = -0.5         # finger grows upward: every cell adjacent
    _, _, count = _suppress_interior_gas_nucleation(
        previous, None, new, None)
    assert count == 0
