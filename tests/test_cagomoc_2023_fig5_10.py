import hashlib
import json

import pytest

from scripts.digitize_cagomoc_2023_fig5_10 import (
    CSV_PATH,
    MANIFEST_PATH,
    POINTS,
    csv_text,
    manifest,
    rows,
)


@pytest.fixture(scope="module")
def record():
    payload = csv_text()
    return payload, manifest(hashlib.sha256(payload.encode("utf-8")).hexdigest())


def test_digitization_is_exactly_replayable(record):
    payload, audit = record
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == audit


def test_all_six_source_setpoints_and_nonmonotone_response_are_retained():
    assert [point.radical_to_ion_ratio for point in POINTS] == [
        0,
        25,
        50,
        100,
        200,
        300,
    ]
    yields = [
        float(row["si_removal_yield_per_cf3_ion"]) for row in rows()
    ]
    assert yields[0] == pytest.approx(1.8094, abs=1.0e-4)
    assert yields[1] == pytest.approx(3.9449, abs=1.0e-4)
    assert yields[2] == pytest.approx(3.8207, abs=1.0e-4)
    assert yields[3] == pytest.approx(2.2659, abs=1.0e-4)
    assert yields[4:] == [0.0, 0.0]
    assert yields[1] > 2.0 * yields[0]
    assert yields[1] > yields[2] > yields[3] > yields[4]


def test_md_board_cannot_identify_krueger_boundary(record):
    _, audit = record
    scope = audit["method_scope"]
    assert scope["evidence_class"].startswith(
        "classical molecular dynamics"
    )
    assert scope["ion"] == "CF3+ represented as a fast neutral at surface impact"
    assert scope["ion_energy_eV"] == 2000
    assert scope["target"] == "flat SiO2"
    boundaries = " ".join(audit["claim_boundaries"])
    assert "Do not use this curve to identify Krueger" in boundaries
    assert "slow in-hole redeposition may be underestimated" in boundaries


def test_pixel_placement_and_uncertainty_are_explicit(record):
    _, audit = record
    pixel = audit["pixel_calibration"]
    assert pixel["maximum_setpoint_placement_offset_ratio_units"] < 0.5
    assert audit["digitization"]["yield_bound"] < 0.03
    assert "occluded" in audit["digitization"]["zero_point_policy"]
