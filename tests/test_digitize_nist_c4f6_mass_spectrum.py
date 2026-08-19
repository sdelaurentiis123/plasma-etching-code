import csv
import io
import json

import pytest

from scripts.digitize_nist_c4f6_mass_spectrum import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    intensity_at_pixel,
    m_over_z_at_pixel,
    manifest_text,
)


def _rows():
    return list(csv.DictReader(io.StringIO(csv_text())))


def test_committed_nist_c4f6_sticks_replay_from_source_pixels():
    payload = csv_text()

    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)
    assert len(_rows()) == 23
    assert m_over_z_at_pixel(400.5) == pytest.approx(93.0, abs=0.2)
    assert intensity_at_pixel(66.0) == 100.0
    assert intensity_at_pixel(506.0) == 0.0


def test_direct_parent_fragmentation_forbids_a_light_cfx_only_reactor_map():
    board = {
        row["assignment"]: float(row["relative_intensity_percent"])
        for row in _rows()
    }

    assert board["C3F3+"] == 100.0
    assert board["C4F6+"] == pytest.approx(43.864, abs=0.001)
    assert board["CF+"] > board["CF3+"] > board["CF2+"]
    assert board["CF2+"] / board["CF+"] == pytest.approx(0.1045, rel=0.002)

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["digitization"]["raw_source_committed"] is False
    forbidden = " ".join(manifest["model_consequence"]["forbidden"])
    assert "only to CF+/CF2+/CF3+" in forbidden
    assert "825 nm" in forbidden
