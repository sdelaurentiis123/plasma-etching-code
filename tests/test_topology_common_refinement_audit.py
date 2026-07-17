from pathlib import Path
import sys

import pytest


SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import topology_common_refinement_audit as audit
from petch.feature_step_3d import _periodic_physical_volume_topology_signature


def test_keyhole_fixture_preserves_declared_physical_domain_and_open_topology():
    geometry = audit.build_keyhole_geometry(0.05)

    assert geometry.phi.shape == (13, 3, 21)
    assert _periodic_physical_volume_topology_signature(
        geometry, (1,)) == (1, 0, False, ((1, 1),))

    with pytest.raises(ValueError, match="divide every physical extent"):
        audit.build_keyhole_geometry(0.03)


def test_rounded_comparator_fixtures_separate_open_and_prescribed_sealed_states():
    open_geometry = audit.build_rounded_keyhole_geometry(0.025)
    sealed_geometry = audit.build_rounded_keyhole_geometry(0.025, sealed=True)

    assert _periodic_physical_volume_topology_signature(
        open_geometry, (1,)) == (1, 0, False, ((1, 1),))
    assert _periodic_physical_volume_topology_signature(
        sealed_geometry, (1,)) == (1, 1, False, ((1, 1),))
    assert audit.ROUNDED_PHYSICAL_GEOMETRY["sealed_cap_thickness_um"] == 0.10


def test_refinement_gate_uses_physical_cell_crossing_not_arbitrary_percentage():
    levels = [
        {"dx_um": 0.05, "closure_time_s": 4.25, "reopening_time_s": 0.75},
        {"dx_um": 0.025, "closure_time_s": 4.375, "reopening_time_s": 7.875},
        {"dx_um": 0.0125, "closure_time_s": 4.125, "reopening_time_s": 7.5},
    ]

    result = audit._refinement_summary(levels)

    assert result["passed"] is True
    assert result["all_adjacent_pairs_passed"] is False
    assert result["coarsest_level_is_diagnostic_only"] is True
    assert result["adjacent_pairs"][0]["passed"] is False
    assert result["authoritative_pair"]["passed"] is True
    assert result["authoritative_pair"][
        "reopening_one_coarse_cell_crossing_time_s"] == pytest.approx(0.5)
