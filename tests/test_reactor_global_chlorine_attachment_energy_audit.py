import csv
import hashlib
import io
import json
from pathlib import Path

from scripts.audit_malyshev_1998_attachment_energy_support import (
    AUDIT_PATH,
    CSV_PATH,
    REPORT_PATH,
    build_outputs,
)


def test_lam_attachment_energy_audit_is_deterministic_and_fail_closed():
    csv_payload, audit, report = build_outputs()
    rows = list(csv.DictReader(io.StringIO(csv_payload)))
    assert len(rows) == 62
    assert audit["coefficient_selection_target"] is None
    assert audit["reactor_fit_target"] is None
    assert audit["feature_depth_target"] is None
    assert audit["marker_accounting"] == {
        "lam_measured_Te_markers": 62,
        "particle_rate_support_complete": 0,
        "incident_energy_support_complete": 0,
    }
    assert audit["electron_power_closure"]["status"] == "blocked"
    assert audit["kernel_tolerance"]["value"] == 1.0e-6
    assert all(row["supports_absorbed_power_solve"] == "false" for row in rows)
    assert all(row["supports_wafer_flux"] == "false" for row in rows)
    assert all(row["supports_feature_depth"] == "false" for row in rows)
    assert hashlib.sha256(csv_payload.encode("utf-8")).hexdigest() == (
        audit["csv_sha256"])
    assert "0/62" in report
    assert "not a power or depth result" in report


def test_checked_in_lam_attachment_energy_receipt_matches_builder():
    csv_payload, audit, report = build_outputs()
    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(AUDIT_PATH.read_text(encoding="utf-8")) == audit
    assert REPORT_PATH.read_text(encoding="utf-8") == report


def test_attachment_energy_receipt_paths_stay_inside_curated_results():
    expected = Path("results/curated/reactor_global_chlorine")
    for path in (CSV_PATH, AUDIT_PATH, REPORT_PATH):
        assert expected.as_posix() in path.as_posix()
