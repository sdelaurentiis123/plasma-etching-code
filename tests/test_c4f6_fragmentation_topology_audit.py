import json

import pytest

from scripts.audit_c4f6_fragmentation_topology import (
    DEFAULT_OUTPUT,
    audit,
)


def test_direct_fragmentation_is_rejected_without_using_depth():
    result = audit()
    diagnostics = result["diagnostics"]

    assert diagnostics["minimum_CF2+/CF+_enhancement_over_direct_EI"] > 4.4
    assert diagnostics["maximum_CF2+/CF+_enhancement_over_direct_EI"] < 5.4
    assert diagnostics["maximum_to_minimum_reactor_CF3+/CF+"] > 2.2
    assert result["physics_decision"]["direct_fragmentation_only_is_adequate"] is False
    assert result["physics_decision"]["secondary_fragment_ionization_required"]
    assert result["physics_decision"]["ion_neutral_conversion_required"]
    assert not any(result["calibration_firewall"].values())


def test_committed_topology_board_is_exact_replay():
    result = audit()
    committed = json.loads(
        (DEFAULT_OUTPUT / "audit.json").read_text(encoding="utf-8")
    )

    assert committed == result
    assert len(result["reactor_conditions"]) == 7
    assert result["direct_70_eV_EI_reference"]["CF2+/CF+"] == pytest.approx(
        0.1044776, rel=1e-6
    )
    assert result["certification"]["supports_absolute_reactor_flux"] is False
    assert result["certification"]["supports_krueger_boundary"] is False
