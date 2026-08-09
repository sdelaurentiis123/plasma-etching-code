import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "base_case_spatial_moments.json"
)


def test_spatial_moment_digitization_is_reproducible():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "digitize_tian_2017_base_spatial_moments.py"),
        "--verify",
    ], check=True)


def test_spatial_moments_are_source_model_data_with_closed_claim_boundary():
    result = json.loads(DATA.read_text(encoding="utf-8"))
    assert result["digitization"]["not_measurement"] is True
    assert result["digitization"]["fitted_to_trapping_factor"] is False
    assert result["digitization"]["model_reduction"] is True
    field = result["axisymmetric_zone_field"]
    assert len(field["radial_edges_cm"]) == 4
    assert len(field["axial_edges_cm"]) == 5
    assert len(field["zones"]) == 3
    assert "feature-depth calibration" in result["claim_boundary"]["not_valid"]
