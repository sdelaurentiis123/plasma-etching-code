import csv
import io
import json

import pytest

from scripts.audit_mahorowala_1998_cl2_fixed_time import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    manifest_text,
    rows,
)


def test_table2_2_transcription_and_dimensional_depth_are_exact():
    payload = csv_text()
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)

    parsed = list(csv.DictReader(io.StringIO(payload)))
    assert parsed == rows()
    assert len(parsed) == 13
    assert sum(row["quantitative_status"] == "usable" for row in parsed) == 11
    assert {
        int(row["run"])
        for row in parsed
        if row["quantitative_status"] != "usable"
    } == {8, 12}

    run1 = parsed[0]
    assert float(run1["derived_poly_si_removed_nm"]) == pytest.approx(
        3150.0 * 75.0 / 60.0 * 0.1
    )
    depths = [
        float(row["derived_poly_si_removed_nm"])
        for row in parsed
        if row["quantitative_status"] == "usable"
    ]
    assert min(depths) == 112.5
    assert max(depths) == 459.375


def test_board_does_not_promote_missing_reactor_boundary_to_prediction():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    identifiability = manifest["identifiability"]

    assert identifiability["absolute_time_published"]
    assert identifiability["absolute_rate_published"]
    assert not identifiability["species_resolved_wafer_flux_published"]
    assert not identifiability["measured_iead_published"]
    assert not identifiability["feature_profile_formal_pass_granted"]
    assert manifest["transcription"]["reported_measurement_uncertainty"] is None
