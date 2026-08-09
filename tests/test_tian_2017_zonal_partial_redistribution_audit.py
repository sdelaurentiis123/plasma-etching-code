import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RESULT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_zonal_partial_redistribution_audit"
    / "tian_2017_zonal_partial_redistribution_audit.json"
)


def _result():
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_zonal_mechanism_is_accurate_without_being_mislabeled_validation():
    result = _result()
    assert result["summary"]["combined_mape_percent"] < 5.0
    assert result["summary"]["maximum_absolute_line_error_percent"] < 10.0
    assert result["summary"]["line_order_reproduced"] is True
    assert result["atomic_or_depth_parameter_fitted"] is False
    assert result["formal_experimental_validation"] is False
    assert result["formal_gate_pass"] is False
    assert result["source_model_target_visible_during_model_selection"] is True


def test_zonal_mechanism_has_independent_numerical_and_conservation_receipts():
    convergence = _result()["numerical_convergence"]
    assert len(convergence) == 2
    for row in convergence:
        assert row["order_10_to_12_relative_change"] < 0.02
        assert row["order_10_to_12_wafer_escape_relative_change"] < 0.01
        assert row["frequency_refinement_relative_change"] < 1.0e-4
        assert row["transition_probability_conservation_error_maximum"] < 1.0e-12
        assert row["terminal_probability_conservation_error_maximum"] < 1.0e-6
        assert row["refined_linear_solver_relative_residual"] < 1.0e-7

    for row in _result()["rows"]:
        assert row["escape_boundary_labels"][0] == (
            "lower_endcap_wafer_plane")
        assert sum(
            row["partial_redistribution_escape_boundary_probability"]
        ) + row["partial_redistribution_quench_probability"] == pytest.approx(
            1.0, abs=1.0e-6)


def test_zonal_audit_names_the_missing_experiment_or_source_export():
    blocker = _result()["remaining_blocker"]
    assert "Ar(1s2) emitter field" in blocker
    assert "raw source-state field export" in blocker
