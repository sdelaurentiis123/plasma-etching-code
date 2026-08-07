import csv
import io
import json

import pytest

from scripts.digitize_benck_2003_figure9 import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    current_at_pixel,
    manifest_text,
)


def _rows():
    return list(csv.DictReader(io.StringIO(csv_text())))


def test_committed_benck_board_is_reproduced_from_full_page_pixels():
    payload = csv_text()

    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)
    assert len(_rows()) == 19
    assert current_at_pixel(4182.5) == pytest.approx(0.1, rel=0.007)
    assert current_at_pixel(4753.5) == pytest.approx(0.01, rel=0.013)
    assert current_at_pixel(5333.5) == pytest.approx(0.001, rel=0.007)


def test_benck_species_trends_match_the_source_text_without_transplant():
    board = {}
    for row in _rows():
        board.setdefault(row["species"], []).append(
            float(row["ion_current_density_mA_cm2"])
        )

    assert board["total_positive_ion_current"] == sorted(
        board["total_positive_ion_current"], reverse=True
    )
    assert board["Ar+"] == sorted(board["Ar+"], reverse=True)
    assert board["CF+"] == sorted(board["CF+"])
    assert board["CF2+"] == sorted(board["CF2+"])
    assert board["CF3+"] == sorted(board["CF3+"])
    assert board["CF+"][3] / board["CF+"][2] < 1.10
    assert board["CF2+"][3] / board["CF2+"][2] < 1.10
    assert board["CF3+"][3] / board["CF3+"][2] > 1.30

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert "Krueger" in " ".join(manifest["claim_boundary"]["not_valid"])
    assert manifest["experiment"]["oxygen_added"] is False
    assert manifest["digitization"]["source_measurement_uncertainty_kept_separate"]
