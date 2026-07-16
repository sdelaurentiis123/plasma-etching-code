import numpy as np
import pytest

from petch.charging_poisson_3d import NodalPoissonSystem3D
from petch.conductor_terminal_3d import RemotePadElectronCollector3D


def _two_conductor_system():
    fixed = np.zeros((3, 2, 3), dtype=bool)
    fixed[:, :, -1] = True
    conductor = np.zeros_like(fixed, dtype=int)
    conductor[0, :, 0] = 1
    conductor[2, :, 0] = 2
    return NodalPoissonSystem3D(
        np.ones((2, 1, 2)),
        np.full(3, 1.0e-7),
        fixed,
        floating_conductor_node_ids=conductor,
    )


def _collector(**overrides):
    values = dict(
        collector_perimeter_m_by_conductor={1: 2.0e-3, 2: 5.0e-3},
        electron_current_per_perimeter_a_m=3.0e-6,
        coefficient_bounds_a_m=(1.0e-7, 1.0e-5),
        source="Nozawa et al., JJAP 34, 2107 (1995), Fig. 5 limit",
        coefficient_evidence="manufactured value for conservation test only",
        topology_evidence="manufactured two-pad topology",
    )
    values.update(overrides)
    return RemotePadElectronCollector3D(**values)


def test_remote_pad_current_routes_exactly_to_floating_conductor_inventories():
    system = _two_conductor_system()
    result = _collector().current_contribution(system)

    assert result.signed_current_a_by_conductor[1] == pytest.approx(-6.0e-9)
    assert result.signed_current_a_by_conductor[2] == pytest.approx(-15.0e-9)
    assert result.positive_node_current_a.sum() == 0.0
    assert result.negative_node_current_a.sum() == pytest.approx(21.0e-9)
    assert result.signed_total_current_a == pytest.approx(-21.0e-9)
    assert result.absolute_total_current_a == pytest.approx(21.0e-9)
    for conductor_id, expected in ((1, 6.0e-9), (2, 15.0e-9)):
        node = system.floating_conductor_representative_node(conductor_id)
        assert result.negative_node_current_a[node] == pytest.approx(expected)
    assert result.provenance["volume_plasma_charge_added"] is False
    assert result.provenance["resolved_line_resistance"] is False


def test_remote_pad_current_refuses_unknown_conductor_and_unbounded_coefficient():
    system = _two_conductor_system()
    with pytest.raises(ValueError, match="absent"):
        _collector(
            collector_perimeter_m_by_conductor={3: 1.0e-3},
        ).current_contribution(system)
    with pytest.raises(ValueError, match="invalid"):
        _collector(
            electron_current_per_perimeter_a_m=2.0e-5,
        )


def test_shared_pad_topology_is_one_current_inventory_not_four_local_sources():
    fixed = np.zeros((5, 2, 3), dtype=bool)
    fixed[:, :, -1] = True
    conductor = np.zeros_like(fixed, dtype=int)
    conductor[[0, 1, 3, 4], :, 0] = 1
    system = NodalPoissonSystem3D(
        np.ones((4, 1, 2)),
        np.full(3, 1.0e-7),
        fixed,
        floating_conductor_node_ids=conductor,
    )
    result = _collector(
        collector_perimeter_m_by_conductor={1: 8.0e-3},
    ).current_contribution(system)

    assert tuple(system.floating_conductor_ids) == (1,)
    assert result.signed_current_a_by_conductor[1] == pytest.approx(-24.0e-9)
    assert np.count_nonzero(result.negative_node_current_a) == 1
