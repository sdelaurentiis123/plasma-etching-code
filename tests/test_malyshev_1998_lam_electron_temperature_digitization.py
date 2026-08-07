import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "scripts" /
    "digitize_malyshev_1998_lam_electron_temperature.py"
)
DATA = ROOT / "data" / "experimental" / "malyshev_1998_lam"
CSV_PATH = DATA / "figure3_electron_temperature.csv"
MANIFEST_PATH = DATA / "electron_temperature_manifest.json"


def _module():
    spec = importlib.util.spec_from_file_location(
        "malyshev_electron_temperature_digitizer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_electron_temperature_board_is_exact_script_replay():
    module = _module()
    csv_payload = module.csv_text()
    csv_sha256 = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()

    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == (
        module.manifest(csv_sha256)
    )


def test_native_axis_calibration_and_digitization_domain_are_explicit():
    module = _module()
    for panel in module.PANELS.values():
        assert panel.power_at_pixel(panel.left_px) == pytest.approx(0.0)
        assert panel.power_at_pixel(panel.right_px) == pytest.approx(950.0)
        assert panel.temperature_at_pixel(panel.top_px) == pytest.approx(
            panel.maximum_temperature_eV)
        assert panel.temperature_at_pixel(panel.bottom_px) == pytest.approx(
            panel.minimum_temperature_eV)

    rows = list(csv.DictReader(io.StringIO(module.csv_text())))
    assert len(rows) == 62
    counts = {}
    for row in rows:
        key = (
            float(row["window_to_wafer_gap_cm"]),
            float(row["pressure_mTorr"]),
        )
        counts[key] = counts.get(key, 0) + 1
        # The complete measured board lies inside the exact Hamilton
        # Maxwellian-rate support; no Te extrapolation is needed to join them.
        assert 0.3 <= float(row["electron_temperature_eV"]) <= 5.0

    assert counts == {
        (11.0, 0.5): 3,
        (11.0, 1.0): 7,
        (11.0, 2.0): 8,
        (11.0, 10.0): 10,
        (11.0, 20.0): 7,
        (6.5, 0.5): 3,
        (6.5, 1.0): 5,
        (6.5, 2.0): 9,
        (6.5, 10.0): 10,
    }


def test_temperature_board_cannot_masquerade_as_power_flux_or_depth_data():
    rows = list(csv.DictReader(CSV_PATH.open(encoding="utf-8")))
    assert {row["validation_role"] for row in rows} == {
        "measured_electron_state_conditioning_input"
    }
    assert {row["tcp_power_semantics"] for row in rows} == {
        "power_into_matching_network_not_absorbed_power"
    }
    assert {row["reported_measurement_uncertainty"] for row in rows} == {
        "not_reported_in_article"
    }
    assert all(row["supports_absorbed_power"] == "false" for row in rows)
    assert all(row["supports_wafer_flux"] == "false" for row in rows)

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["digitization"]["marker_count"] == 62
    assert "never digitized" in manifest["exclusions"]["smooth_curves"]
    assert "omitted" in manifest["exclusions"]["low_power_overlap"]
    assert manifest["use_boundary"]["supports_measured_Te_conditioning"]
    assert not manifest["use_boundary"]["supports_absorbed_power_validation"]
    assert not manifest["use_boundary"]["supports_wafer_flux_validation"]
    assert not manifest["use_boundary"]["supports_feature_depth_validation"]
