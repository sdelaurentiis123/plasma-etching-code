import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
CSV_PATH = DATA / "digitized_figures_5_10_5_12.csv"
MANIFEST_PATH = DATA / "digitization_manifest.json"


def _rows():
    with CSV_PATH.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _series(name):
    return [
        float(row["value"]) for row in _rows() if row["series"] == name
    ]


def test_tian_digitization_is_checksum_pinned_and_model_labeled():
    rows = _rows()
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert len(rows) == manifest["digitization"]["marker_count"] == 97
    assert manifest["source"]["pdf_sha256"] == (
        "4d260ab9e85240bd051ccb3ba32cd047b4ac1ddb6c309dbbca4093822e37790b")
    assert manifest["source"]["render_sha256_by_pdf_page"] == {
        "171": "c0fc7cc6f06aef081ce46f81116bf06963d8886d20754cbad85dbb0aba0da63c",
        "173": "21f084f3cd7c309db86dc656be68c2ff991481aac56f990b80847c15a7f3c863",
    }
    assert manifest["output"]["sha256"] == hashlib.sha256(
        CSV_PATH.read_bytes()).hexdigest()
    assert {row["evidence_type"] for row in rows} == {
        "source_equipment_model_digitized"}
    assert manifest["experiment_context"]["not_measurement"] is True


def test_tian_mixture_and_trapping_trends_replay_source_markers():
    photon = _series("total_photon_flux")
    beta = _series("total_photon_to_total_ion_flux_ratio_beta")
    trap_1067 = _series("Ar_106.7_nm_trapping_factor")
    trap_1048 = _series("Ar_104.8_nm_trapping_factor")
    trap_139 = _series("Cl_139_nm_trapping_factor")
    assert all(left > right for left, right in zip(photon, photon[1:]))
    assert all(left > right for left, right in zip(beta, beta[1:]))
    assert trap_1067[0] > 500.0
    assert trap_1067[-1] < 20.0
    assert trap_1048[0] > 200.0
    assert trap_1048[-1] < 60.0
    assert 3.0 < min(trap_139) < 5.0
    assert 10.0 < max(trap_139) < 12.0


def test_tian_digitization_claim_boundary_excludes_validation_and_depth():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    invalid = manifest["claim_boundary"]["not_valid"]
    assert "experimental validation" in invalid
    assert "feature-depth calibration or validation" in invalid
    assert "surface photo-etch yield" in invalid
