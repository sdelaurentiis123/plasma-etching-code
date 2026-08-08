import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_malyshev_1998_lam_transport.py"
OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"
CSV_PATH = OUTPUT / "malyshev_1998_lam_transport_diagnostic.csv"
AUDIT_PATH = OUTPUT / "malyshev_1998_lam_transport_diagnostic.json"
REPORT_PATH = OUTPUT / "MALYSHEV_1998_TRANSPORT_DIAGNOSTIC.md"


def _module():
    spec = importlib.util.spec_from_file_location(
        "malyshev_transport_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_transport_board_is_exact_script_replay():
    csv_payload, audit, report = _module().build_outputs()

    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(AUDIT_PATH.read_text(encoding="utf-8")) == audit
    assert REPORT_PATH.read_text(encoding="utf-8") == report
    assert audit["csv_sha256"] == hashlib.sha256(
        csv_payload.encode("utf-8")).hexdigest()


def test_transport_board_accounts_for_rows_and_fails_closed():
    rows = list(csv.DictReader(io.StringIO(
        CSV_PATH.read_text(encoding="utf-8"))))
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))

    assert len(rows) == 23
    accounting = audit["marker_accounting"]
    assert accounting["successful_measured_state_rows"] == 23
    assert accounting["transport_attainable_rows"] == 22
    assert accounting["transport_unattainable_rows"] == 1
    assert sum(
        row["target_is_transport_attainable"] == "true"
        for row in rows
    ) == 22
    failed = next(
        row for row in rows
        if row["target_is_transport_attainable"] == "false")
    assert float(failed["required_wall_return_frequency_s_inv"]) > (
        float(failed["absorbing_wall_limit_frequency_s_inv"]))
    assert failed["effective_wall_recombination_probability"] == ""
    assert failed["matched_wall_return_frequency_s_inv"] == ""

    for row in rows:
        assert row["supports_prediction"] == "false"
        assert row["supports_local_wall_probability_prediction"] == "false"
        assert row["supports_wafer_flux"] == "false"
        assert row["supports_feature_depth"] == "false"


def test_transport_reconstruction_retires_wrong_temperature_constant_only():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    diffusivity = audit["diffusivity"]

    assert diffusivity[
        "room_temperature_298p15K_one_atm_cm2_s"] == pytest.approx(
            0.15167523928085008)
    assert diffusivity["333K_one_atm_cm2_s"] == pytest.approx(
        0.1862154115776554)
    assert diffusivity["333K_reduced_diffusivity_m_inv_s"] == (
        pytest.approx(4.103975103414651e20))
    assert diffusivity["new_to_economou_ratio"] < 0.67
    assert not diffusivity["supports_prediction"]
    assert audit["coefficient_selection_target"] is None
    assert audit["reactor_fit_target"] is None
    assert audit["feature_depth_target"] is None


def test_published_gamma_replay_is_labeled_same_board_not_validation():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    replay = audit["published_gamma_forward_replay"]

    assert replay["gamma"] == 0.035
    assert "same Figures 7-8 board" in replay["gamma_role"]
    assert replay["formal_gate"] is None
    assert replay["mean_absolute_error_percentage_point"] > 9.0
    assert replay["root_mean_square_error_percentage_point"] > 12.0
    assert audit["gas_temperature"]["measured_at_each_power"] is False
    assert audit["cross_source_wall_comparison"]["status"] == "not_scored"
