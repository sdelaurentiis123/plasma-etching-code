import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = (
    ROOT / "data" / "experimental" / "wang_olthoff_1999")
CSV_PATH = DATA_DIRECTORY / "figure9_chlorine_ion_flux.csv"
MANIFEST_PATH = DATA_DIRECTORY / "digitization_manifest.json"
PROTOCOL_PATH = (
    ROOT / "results" / "curated" / "reactor_global_chlorine"
    / "wang_olthoff_1999_preregistered_gate.json")


def _rows():
    with CSV_PATH.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_wang_olthoff_figure9_data_are_complete_and_pixel_pinned():
    rows = _rows()
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert len(rows) == 24
    assert manifest["source"]["pdf_sha256"] == (
        "17702f0ffb904cca42760867c693c382b174c68f02ef7b5064ec95933adea460")
    assert manifest["source"]["rendered_page_sha256"] == (
        "44a1953b90353711dfb6798b357ecb90eb0fbc10abb1c49d725a21a39a54b12a")
    assert manifest["source"]["rendered_page_size_px"] == [5096, 6716]
    assert manifest["output"]["sha256"] == hashlib.sha256(
        CSV_PATH.read_bytes()).hexdigest()
    assert manifest["digitization"]["measurement_uncertainty"].startswith(
        "not assigned")


def test_wang_olthoff_pure_chlorine_species_close_measured_total():
    rows = _rows()
    series = {}
    for row in rows:
        if row["panel"] == "a":
            series.setdefault(row["species"], []).append(
                float(row["ion_flux_mA_cm2"]))
    total = np.asarray(series["total_positive_ion"])
    atomic = np.asarray(series["Cl+"])
    molecular = np.asarray(series["Cl2+"])
    np.testing.assert_array_less(total[:-1], total[1:])
    np.testing.assert_array_less(molecular[1:], molecular[:-1])
    np.testing.assert_array_less(0.8, atomic / total)
    np.testing.assert_allclose(
        (atomic + molecular) / total,
        np.ones(3),
        rtol=0.0,
        atol=0.03,
    )


def test_wang_olthoff_gate_is_all_held_out_and_caps_claims():
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    assert not protocol["frozen_before_model"][
        "mixed_Ar_Cl2_global_solver_exists"]
    assert protocol["frozen_before_model"]["coefficient_selection_target"] is None
    assert protocol["splits"]["calibration"] == []
    assert len(protocol["splits"]["held_out"]) == 2
    assert protocol["acceptance"]["total_positive_ion_flux"][
        "maximum_mean_absolute_percentage_error"] == 0.3
    claim_cap = protocol["claim_cap"]
    assert "conditional" in claim_cap["if_all_gates_pass"]
    assert "feature-depth prediction" in claim_cap["not_earned"]
    assert "elementary Cl2 ionization branching identification" in (
        claim_cap["not_earned"])
