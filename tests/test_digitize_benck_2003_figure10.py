import csv
import io
import json

import pytest

from scripts.digitize_benck_2003_figure10 import (
    CSV_PATH,
    MANIFEST_PATH,
    cross_figure_relative_differences,
    csv_text,
    current_at_pixel,
    manifest_text,
)


def _rows():
    return list(csv.DictReader(io.StringIO(csv_text())))


def test_figure10_replays_and_reconciles_the_independent_figure9_drawing():
    payload = csv_text()

    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)
    assert len(_rows()) == 15
    assert current_at_pixel(1065.5) == pytest.approx(0.1, rel=0.004)
    assert current_at_pixel(1637.5) == pytest.approx(0.01, rel=0.006)
    assert current_at_pixel(2213.5) == pytest.approx(0.001, rel=0.004)

    differences = cross_figure_relative_differences()
    assert set(differences) == {
        "Ar+", "CF+", "CF2+", "CF3+", "total_positive_ion_current"
    }
    assert max(map(abs, differences.values())) < 0.036


def test_pressure_board_keeps_physics_scope_and_uncertainties_separate():
    board = {}
    for row in _rows():
        board.setdefault(row["species"], []).append(
            float(row["ion_current_density_mA_cm2"])
        )

    assert board["total_positive_ion_current"] == sorted(
        board["total_positive_ion_current"], reverse=True
    )
    assert board["CF+"] == sorted(board["CF+"], reverse=True)
    assert board["CF2+"] == sorted(board["CF2+"], reverse=True)
    assert board["CF3+"][2] > board["CF3+"][1]

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["cross_figure_reconciliation"][
        "passed_within_individual_digitization_bound"
    ]
    assert manifest["experiment"]["corrected_transmission_relative_uncertainty"] == 0.20
    assert "Krueger" in " ".join(manifest["claim_boundary"]["not_valid"])
