import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_malyshev_1998_lam_eq7_inversion.py"
OUTPUT = ROOT / "results" / "curated" / "reactor_global_chlorine"
CSV_PATH = OUTPUT / "malyshev_1998_lam_eq7_inversion.csv"
AUDIT_PATH = OUTPUT / "malyshev_1998_lam_eq7_inversion.json"


def _module():
    spec = importlib.util.spec_from_file_location(
        "malyshev_eq7_inversion_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_eq7_board_is_exact_script_replay():
    csv_payload, audit = _module().build_outputs()

    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(AUDIT_PATH.read_text(encoding="utf-8")) == audit
    assert audit["csv_sha256"] == hashlib.sha256(
        csv_payload.encode("utf-8")).hexdigest()


def test_eq7_board_closes_algebra_and_accounts_for_every_marker():
    rows = list(csv.DictReader(io.StringIO(
        CSV_PATH.read_text(encoding="utf-8"))))
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))

    assert len(rows) == 23
    accounting = audit["marker_accounting"]
    assert accounting["audited_marker_total"] == 38
    assert accounting["successful_inversions"] == 23
    assert accounting["excluded"] == {
        "diagnostic_flow_check": 1,
        "nonphysical_or_zero_derived_dissociation": 4,
        "electron_temperature_support_missing": 2,
        "electron_density_support_missing": 8,
    }
    assert accounting["electron_state_method_pairs"] == {
        "Te:exact_marker|ne:exact_marker": 11,
        "Te:exact_marker|ne:linear_interpolation": 11,
        "Te:linear_interpolation|ne:linear_interpolation": 1,
    }

    for row in rows:
        measured = float(row["relative_cl2_density_percent"])
        reproduced = float(row["reproduced_relative_cl2_density_percent"])
        assert reproduced == pytest.approx(measured, abs=1.0e-7)
        assert row["supports_wall_probability_inference"] == "false"
        assert row["supports_wafer_flux"] == "false"
        assert row["supports_feature_depth"] == "false"


def test_eq7_audit_reports_ill_conditioning_instead_of_hiding_it():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    summary = audit["required_wall_return_frequency_s_inv"]

    assert summary["maximum_all"] > 1.0e4
    assert summary["median_all"] < 200.0
    assert summary["finite_reported_cl2_upper_envelope_count"] < 23
    assert "no formal pass/fail" in audit["uncertainty_boundary"]
    assert audit["coefficient_selection_target"] is None
    assert audit["feature_depth_target"] is None


def test_source_footnote_rate_is_quarantined_by_its_own_printed_law():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    consistency = audit["source_internal_consistency"]

    assert consistency["status"] == (
        "printed_nominal_values_do_not_share_one_condition")
    assert consistency["footnote_to_printed_upper_ratio"] > 6.0
    assert consistency["footnote14_destruction_frequency_s_inv"] == 700.0
    assert consistency[
        "printed_upper_destruction_frequency_s_inv"] < 120.0
    assert "quarantined" in consistency["interpretation"]
