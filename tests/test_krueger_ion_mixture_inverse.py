import pytest

from scripts.audit_krueger_ion_mixture_inverse import build_audit
from scripts.krueger_2024_trench_pilot import _parse_guo_aggregate_ion_mixture


def test_combined_iead_inverse_is_rank_deficient_and_depth_firewalled():
    audit = build_audit()

    assert audit["status"] == "structurally_nonidentifiable_from_published_iead"
    inverse = audit["published_inverse_problem"]
    assert inverse["aggregate_operator_rank"] == 1
    assert inverse["composition_contrast_rank"] == 0
    assert inverse["composition_nullity"] == 10
    assert audit["figure_resolution_gate"]["mass_shift_per_digitization_bin"] == 0.24
    assert not audit["calibration_firewall"][
        "feature_depth_used_to_choose_mixture"]

    boards = audit["cross_reactor_measured_composition_sensitivities"]["boards"]
    for board in boards.values():
        assert abs(board["resolved_plus_unresolved_closure"] - 1.0) < 1.0e-12
        assert board["resolved_fraction_of_total_positive_current"]["Ar+"] > 0.3

    propagated = audit[
        "cross_reactor_measured_composition_sensitivities"]["surface_propagation"]
    for cases in propagated.values():
        assert set(cases) == {
            "unresolved_as_inert", "unresolved_as_CF", "unresolved_as_CF2",
            "unresolved_as_CF3", "unresolved_as_C3F3",
        }
        for case in cases.values():
            assert case["feature_depth_used"] is False
            assert abs(case["atom_ledger_residual"]) < 1.0e-10
            assert "score_only_after_solve" in case


def test_feature_cli_mixture_parser_closes_fraction_and_refuses_bad_json():
    assert _parse_guo_aggregate_ion_mixture(
        '{"CF+":0.1,"CF2":0.2}') == {"CF+": 0.1, "CF2": 0.2}
    with pytest.raises(ValueError, match="JSON object"):
        _parse_guo_aggregate_ion_mixture("not-json")
    with pytest.raises(ValueError, match="sum to at most one"):
        _parse_guo_aggregate_ion_mixture('{"CF":0.7,"CF2":0.4}')
