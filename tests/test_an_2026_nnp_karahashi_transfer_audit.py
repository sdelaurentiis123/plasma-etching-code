from pathlib import Path

import pytest

from scripts.audit_an_2026_nnp_karahashi_transfer import build_audit


ROOT = Path(__file__).parents[1]


@pytest.fixture(scope="module")
def audit():
    return build_audit(ROOT)


def test_exact_overlap_has_no_fit_or_interpolation(audit):
    assert len(audit["point_board"]) == 9
    assert {
        row["species"] for row in audit["point_board"]
    } == {"CF+", "CF2+", "CF3+"}
    firewall = audit["calibration_firewall"]
    assert firewall["petch_fitted_parameters"] == []
    assert not firewall["karahashi_yields_used_in_nnp_training_loss"]
    assert not firewall["comparison_is_blind"]
    overlap = audit["experiment_and_exact_overlap"]
    assert overlap["model_points_interpolated"] == 0
    assert overlap["model_points_extrapolated"] == 0


def test_source_receipts_are_pinned(audit):
    receipts = audit["source_receipts"]
    assert receipts["transcribed_model_data_sha256"] == (
        "c897427ffb8055ffe71f042530efc8acda30c05bcdfe52df0fb4542af1c295d4")
    assert receipts["author_figure3_yaml_sha256"] == (
        "7fdb3b72bc55eb47e8e2fbf6b218ddcacdb869136b2978c4a90bebe3a3dafb9a")
    assert receipts["author_repository_commit"] == (
        "4bcd035090b9f652cda10150c4da4b662143b34e")
    code = receipts["author_code_audit"]
    assert code["source_sha256"] == (
        "a0471f888bd885a84b7188a4111aa18a4c7c6a68d8061e0b7d5e6193204babd0")
    assert code["implementation_or_weights_copied_into_petch"] is False
    assert audit["model_evidence_boundary"][
        "product_escape_transport_solved"] is False
    assert audit["model_evidence_boundary"][
        "diagnostic_escape_depth_status"
    ] == "arbitrary_sensitivity_not_measured_parameter"
    assert receipts["karahashi_data_sha256"] == (
        receipts["loader_pinned_karahashi_data_sha256"])


def test_atomistic_model_materially_improves_same_support_error(audit):
    comparison = audit["same_support_baseline"]
    assert comparison["support_is_identical"]
    assert comparison["nn_p_zbl_metrics"][
        "mean_absolute_relative_error"] == pytest.approx(
            0.233484, rel=2e-5)
    assert comparison["nn_p_zbl_metrics"][
        "mean_absolute_relative_error"] < 0.25
    assert comparison["guo_metrics"][
        "mean_absolute_relative_error"] > 0.45
    assert comparison["mape_reduction_fraction"] > 0.50


def test_low_energy_product_escape_defect_is_not_hidden(audit):
    board = {
        (row["species"], row["energy_eV"]): row
        for row in audit["point_board"]
    }
    assert board[("CF2+", 500.0)]["signed_relative_error"] > 0.65
    assert board[("CF3+", 250.0)]["signed_relative_error"] > 0.70
    assert not audit["physics_diagnosis"][
        "atomistic_transfer_validated_over_full_overlap"]
    assert audit["model_evidence_boundary"][
        "released_sio2_energy_maximum_eV"] == 1000.0
    assert not audit["model_evidence_boundary"][
        "executable_code_or_potential_imported"]


def test_high_energy_subset_is_recorded_as_diagnostic_only(audit):
    subsets = audit["diagnostic_subsets_not_preregistered"]
    assert subsets["energy_at_least_750_eV"]["point_count"] == 6
    assert subsets["energy_at_least_750_eV"][
        "mean_absolute_relative_error"] < 0.11
    assert subsets["energy_exactly_1000_eV"]["point_count"] == 3
    assert subsets["energy_exactly_1000_eV"][
        "maximum_absolute_relative_error"] < 0.10
    assert audit["physics_diagnosis"][
        "high_energy_prompt_event_kernel_supported"
    ] == "promising_diagnostic_not_preregistered_validation"


def test_committed_audit_is_current(audit):
    path = (
        ROOT / "results" / "curated" / "an_nnp_karahashi_transfer"
        / "audit.json")
    assert path.exists()
    import json
    assert json.loads(path.read_text(encoding="utf-8")) == audit
