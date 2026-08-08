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
    ROOT / "scripts" / "digitize_malyshev_1998_lam_electron_density.py"
)
DATA = ROOT / "data" / "experimental" / "malyshev_1998_lam"
CSV_PATH = DATA / "figure11_electron_density.csv"
MANIFEST_PATH = DATA / "electron_density_manifest.json"
PACKAGE_CSV_PATH = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "malyshev_1998_lam_electron_density.csv"
)


def _module():
    spec = importlib.util.spec_from_file_location(
        "malyshev_electron_density_digitizer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_electron_density_board_is_exact_script_replay():
    module = _module()
    csv_payload = module.csv_text()
    csv_sha256 = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()

    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert PACKAGE_CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == (
        module.manifest(csv_sha256)
    )


def test_axes_groups_and_gap_response_reproduce_the_native_figure():
    module = _module()
    for panel in module.PANELS.values():
        assert panel.power_at_pixel(panel.left_px) == pytest.approx(0.0)
        assert panel.power_at_pixel(panel.right_px) == pytest.approx(950.0)
        assert panel.density_at_pixel(panel.top_px) == pytest.approx(
            panel.maximum_electron_density_cm3)
        assert panel.density_at_pixel(panel.bottom_px) == pytest.approx(0.0)

    rows = list(csv.DictReader(io.StringIO(module.csv_text())))
    counts = {}
    for row in rows:
        key = (
            float(row["window_to_wafer_gap_cm"]),
            float(row["pressure_mTorr"]),
        )
        counts[key] = counts.get(key, 0) + 1
        assert float(row["volume_average_electron_density_cm3"]) > 0.0
    assert counts == {
        (11.0, 0.5): 4,
        (11.0, 1.0): 5,
        (11.0, 2.0): 6,
        (11.0, 10.0): 4,
        (6.5, 2.0): 5,
        (6.5, 10.0): 3,
    }

    large_500 = next(
        row for row in rows
        if float(row["window_to_wafer_gap_cm"]) == 11.0
        and float(row["pressure_mTorr"]) == 2.0
        and 490.0 < float(row["tcp_source_power_W"]) < 510.0
    )
    small_500 = next(
        row for row in rows
        if float(row["window_to_wafer_gap_cm"]) == 6.5
        and float(row["pressure_mTorr"]) == 2.0
        and 490.0 < float(row["tcp_source_power_W"]) < 510.0
    )
    ratio = (
        float(large_500["volume_average_electron_density_cm3"])
        / float(small_500["volume_average_electron_density_cm3"])
    )
    assert ratio == pytest.approx(1.99, rel=0.01)


def test_volume_average_density_cannot_masquerade_as_local_flux_evidence():
    rows = list(csv.DictReader(CSV_PATH.open(encoding="utf-8")))
    assert len(rows) == 27
    assert {row["reported_measurement_uncertainty"] for row in rows} == {
        "not_reported_in_article"
    }
    assert {row["validation_role"] for row in rows} == {
        "measured_volume_average_electron_state_conditioning_input"
    }
    assert all(
        row["supports_local_wafer_electron_density"] == "false"
        for row in rows
    )
    assert all(row["supports_wafer_flux"] == "false" for row in rows)
    assert all(
        row["volume_average_conversion"]
        == "radial_symmetry_and_axial_sin_pi_h_over_gap"
        for row in rows
    )

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["digitization"]["marker_count"] == 27
    assert "never digitized" in manifest["exclusions"]["linear_lines"]
    assert "omitted" in manifest["exclusions"]["low_power_panel_b_overlap"]
    assert manifest["use_boundary"][
        "supports_measured_volume_average_ne_conditioning"
    ]
    assert not manifest["use_boundary"][
        "supports_local_wafer_electron_density"
    ]
    assert not manifest["use_boundary"]["supports_wafer_flux_validation"]
    assert not manifest["use_boundary"]["supports_feature_depth_validation"]
