import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_partial_redistribution_audit"
    / "tian_2017_partial_redistribution_audit.json"
)


def _result():
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_partial_redistribution_audit_is_held_out_and_fails_closed():
    result = _result()
    assert result["summary"]["held_out_marker_count"] == 22
    assert result["source_model_targets_used_to_fit_parameters"] is False
    assert result["model"]["fitted_parameters"] == []
    assert result["summary"]["formal_gate_pass"] is False
    assert result["formal_experimental_validation"] is False


def test_partial_redistribution_materially_improves_the_missing_physics_board():
    summary = _result()["summary"]
    assert summary["best_combined_mape_percent"] < (
        summary["previous_complete_redistribution_best_mape_percent"])
    assert summary["best_line_trapping_ordering_failure_count"] < (
        summary["previous_complete_redistribution_ordering_failure_count"])


def test_split_spatial_moment_is_accurate_but_not_mislabeled_validation():
    diagnostic = _result()["base_case_split_spatial_moment_diagnostic"]
    assert diagnostic["combined_mape_percent"] < 2.0
    assert diagnostic["target_used_during_model_selection"] is True
    assert diagnostic["formal_validation"] is False


def test_frequency_propagator_has_grid_range_and_solver_receipts():
    convergence = _result()["numerical_convergence"]
    assert len(convergence) == 2
    for row in convergence:
        assert row["coarse_to_refined_relative_change"] < 1.0e-5
        assert row["refined_to_extended_half_range_relative_change"] < 5.0e-5
        assert row["refined_linear_solver_relative_residual"] < 1.0e-7
