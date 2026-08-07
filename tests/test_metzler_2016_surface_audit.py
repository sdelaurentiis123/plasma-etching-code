import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "metzler_2016_surface_closure"
    / "legacy_mixed_layer_audit.json")


def test_metzler_surface_audit_records_falsification_not_calibration():
    report = json.loads(AUDIT.read_text(encoding="utf-8"))
    assert report["audit_id"] == "METZLER-2016-SURFACE-CLOSURE-LEGACY-R1"
    assert report["source"]["points_used"] == 6
    assert "no model parameter" in report["operation"]
    assert not report["verdict"]["legacy_delta_energy_closure_passes"]
    assert not report["verdict"]["unreported_IEDF_closes_model"]
    assert len(report["replays"]) == 36


def test_every_delta_replay_overpredicts_even_across_film_composition_bracket():
    report = json.loads(AUDIT.read_text(encoding="utf-8"))
    delta = [
        row for row in report["replays"]
        if row["spectrum_sensitivity_case"] == "delta_at_reported_maximum"
    ]
    assert {row["assumed_initial_film_F_over_C"] for row in delta} == {
        0.4, 1.0, 1.5, 2.0}
    assert min(row["prediction_over_measurement"] for row in delta) > 4.0
    assert (
        report["integration_convergence"][
            "maximum_relative_prediction_change"]
        < 1e-3)
