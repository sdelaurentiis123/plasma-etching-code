import csv
import io
import json

import pytest

from scripts.digitize_hong_2023_tio2 import (
    FIGURE2_CSV_PATH,
    FIGURE3_CSV_PATH,
    MANIFEST_PATH,
    figure2_csv_text,
    figure2_rows,
    figure3_csv_text,
    figure3_rows,
    manifest_text,
)


def test_generated_hong_tables_and_manifest_are_current():
    figure2 = figure2_csv_text()
    figure3 = figure3_csv_text()
    assert FIGURE2_CSV_PATH.read_text(encoding="utf-8") == figure2
    assert FIGURE3_CSV_PATH.read_text(encoding="utf-8") == figure3
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(
        figure2, figure3
    )
    assert len(list(csv.DictReader(io.StringIO(figure2)))) == 42
    assert len(list(csv.DictReader(io.StringIO(figure3)))) == 6


def test_digitization_replays_paper_text_and_arde_direction():
    rows = figure2_rows()
    c4_p1_selectivity = {
        row["plasma_mode"]: float(row["value"])
        for row in rows
        if row["quantity"] == "tio2_acl_selectivity"
        and row["chemistry"] == "C4F8/SF6/Ar"
        and row["specimen"] == "P1"
    }
    assert c4_p1_selectivity["S30_B70"] == pytest.approx(3.9, abs=0.11)
    assert c4_p1_selectivity["S70_B30"] == pytest.approx(5.9, abs=0.11)

    arde = {
        (row["chemistry"], row["plasma_mode"]): float(row["value"])
        for row in rows
        if row["quantity"] == "arde_p1_depth_over_p3_depth"
    }
    assert arde[("C4F8/SF6/Ar", "CW")] == pytest.approx(1.5, abs=0.03)
    assert arde[("C4F8/SF6/Ar", "S70_B30")] == pytest.approx(1.1, abs=0.03)
    assert arde[("BCl3/CF4/Ar", "CW")] > arde[("C4F8/SF6/Ar", "CW")]


def test_error_bars_and_figure3_are_kept_as_distinct_evidence():
    arde_rows = [
        row
        for row in figure2_rows()
        if row["quantity"] == "arde_p1_depth_over_p3_depth"
    ]
    assert all(row["experimental_error_low"] for row in arde_rows)
    assert all(
        float(row["experimental_error_low"])
        < float(row["value"])
        < float(row["experimental_error_high"])
        for row in arde_rows
    )
    assert {row["evidence_type"] for row in figure3_rows()} == {
        "author_annotated_sem"
    }
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    boundary = " ".join(manifest["claim_boundary"]["not_valid"])
    assert "Zhu" in boundary
    assert "radial" in boundary
