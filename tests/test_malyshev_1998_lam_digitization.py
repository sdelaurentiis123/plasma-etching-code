import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "digitize_malyshev_1998_lam_dissociation.py"
DATA = ROOT / "data" / "experimental" / "malyshev_1998_lam"
CSV_PATH = DATA / "figures7_8_chlorine_dissociation.csv"
MANIFEST_PATH = DATA / "digitization_manifest.json"
PACKAGE_CSV_PATH = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "malyshev_1998_lam_chlorine_dissociation.csv"
)


def _module():
    spec = importlib.util.spec_from_file_location("malyshev_digitizer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_lam_dataset_is_exact_script_replay():
    module = _module()
    csv_payload = module.csv_text()
    csv_sha256 = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert PACKAGE_CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == (
        module.manifest(csv_sha256))


def test_lam_markers_preserve_observable_and_use_boundaries():
    rows = list(csv.DictReader(io.StringIO(_module().csv_text())))
    assert len(rows) == 38
    assert {float(row["window_to_wafer_gap_cm"]) for row in rows} == {6.5, 11.0}
    assert {float(row["pressure_mTorr"]) for row in rows} == {0.5, 1.0, 2.0, 10.0}
    for row in rows:
        relative_cl2 = float(row["relative_cl2_density_percent"])
        dissociation = float(row["cl2_dissociation_percent"])
        assert dissociation == pytest.approx(100.0 - relative_cl2, abs=1e-4)
        assert row["error_bar_semantics"].endswith("not_sigma")
        assert row["tcp_power_semantics"] == (
            "power_into_matching_network_not_absorbed_power")
        assert row["supports_absorbed_power"] == "false"
        assert row["supports_wafer_flux"] == "false"


def test_lam_dataset_excludes_documented_low_power_anomaly_and_model_curves():
    rows = list(csv.DictReader(io.StringIO(_module().csv_text())))
    anomalous = [
        row for row in rows
        if float(row["pressure_mTorr"]) == 10.0
        and float(row["window_to_wafer_gap_cm"]) == 6.5
        and 50.0 <= float(row["tcp_source_power_W"]) <= 120.0
    ]
    assert anomalous == []
    assert {row["validation_role"] for row in rows} == {
        "reactor_dissociation_validation_candidate",
        "diagnostic_flow_check",
    }
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert "never digitized" in manifest["exclusions"]["smooth_curves"]
    assert manifest["use_boundary"]["supports_reactor_dissociation_validation"]
    assert not manifest["use_boundary"]["supports_feature_depth_validation"]
