import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_resonance_escape_audit"
    / "tian_2017_resonance_escape_audit.json"
)


def test_tian_escape_audit_is_held_out_and_fails_closed():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    assert result["summary"]["held_out_marker_count"] == 22
    assert result["source_model_targets_used_to_fit_parameters"] is False
    assert result["model"]["fitted_parameters"] == []
    assert result["summary"]["formal_gate_pass"] is False
    assert result["formal_experimental_validation"] is False


def test_tian_escape_audit_exposes_spatial_state_not_a_scalar_multiplier():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    assert result["summary"]["best_combined_mape_percent"] > 10.0
    assert result["summary"]["best_line_trapping_ordering_failure_count"] > 0
    missing = result["discovered_missing_state"]
    assert "spatial gas-temperature and ground-state-density field" in missing
    assert "partial frequency redistribution after resonant absorption" in missing
