import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "data"
    / "experimental"
    / "humbird_graves_2004"
    / "seminar_surface_state_curves.csv"
)
MANIFEST = DATA.with_name("digitization_manifest.json")


def _rows():
    with DATA.open(newline="") as stream:
        return list(csv.DictReader(stream))


def _series(panel, quantity):
    selected = [
        row for row in _rows()
        if row["source_panel"] == panel and row["quantity"] == quantity
    ]
    selected.sort(key=lambda row: float(row["cf2_fluence_1e15_cm2"]))
    return (
        np.asarray(
            [float(row["cf2_fluence_1e15_cm2"]) for row in selected]),
        np.asarray([float(row["digitized_value"]) for row in selected]),
    )


def test_manifest_and_csv_keep_non_archival_evidence_grade_explicit():
    manifest = json.loads(MANIFEST.read_text())
    rows = _rows()
    assert manifest["observation_count"] == len(rows)
    assert manifest["source"]["evidence_grade"] == (
        "primary_author_seminar_not_peer_reviewed")
    assert not manifest["digitization"]["curve_fit_used"]
    assert not manifest["digitization"]["copyrighted_rasters_redistributed"]
    assert {
        row["evidence_grade"] for row in rows
    } == {"primary_author_seminar_not_peer_reviewed"}


def test_board_resolves_energy_and_thermal_f_surface_response():
    _, carbon_20ev = _series(
        "20_eV_surface_state", "surface_C_uptake_ML")
    _, carbon_200ev = _series(
        "200_eV_surface_state", "surface_C_uptake_ML")
    assert carbon_20ev[-1] < 6.0
    assert carbon_200ev[-1] > 15.0

    _, yield_9_to_1 = _series(
        "200_eV_etch_yield", "Si_etch_yield_per_ion")
    _, yield_10_percent_f = _series(
        "10_percent_F_etch_yield", "Si_etch_yield_per_ion")
    _, yield_20_percent_f = _series(
        "20_percent_F_etch_yield", "Si_etch_yield_per_ion")
    assert yield_9_to_1[-1] < yield_10_percent_f[-1]
    assert yield_10_percent_f[-1] < yield_20_percent_f[-1]


def test_silicon_has_a_finite_cumulative_transport_history():
    for panel in (
            "200_eV_surface_state",
            "10_percent_F_surface_state",
            "20_percent_F_surface_state"):
        _, cumulative = _series(panel, "cumulative_Si_etch_ML")
        assert np.all(np.diff(cumulative) > 0.0)
