import numpy as np
import pytest

from petch.feature_step_3d import make_rectangular_trench_geometry_3d
from petch.krueger_replay_3d import make_krueger_2024_poisson_system_3d


def _geometry():
    return make_rectangular_trench_geometry_3d(
        cell_width=0.13,
        cell_length=0.02,
        domain_height=2.8,
        dx=0.05,
        opening_width=0.09,
        mask_thickness=0.85,
        substrate_top=1.8,
        etched_depth=0.2,
    )


def test_krueger_poisson_system_matches_published_electrostatic_boundaries():
    geometry = _geometry()
    system = make_krueger_2024_poisson_system_3d(
        geometry, mask_relative_permittivity=3.3)

    assert system.shape == geometry.phi.shape
    assert system.periodic_axes == (0, 1)
    assert np.all(system.dirichlet_mask[:, :, 0])
    assert not np.any(system.dirichlet_mask[:, :, 1:])
    assert set(np.unique(system.epsilon_r)) == {1.0, 3.3, 3.9}

    potential, diagnostics = system.solve(np.zeros(system.shape))
    np.testing.assert_array_equal(potential, np.zeros(system.shape))
    assert diagnostics.charge_balance_c == 0.0


def test_krueger_poisson_system_refuses_hidden_mask_material_assumption():
    with pytest.raises(ValueError, match="positive and finite"):
        make_krueger_2024_poisson_system_3d(
            _geometry(), mask_relative_permittivity=0.0)
