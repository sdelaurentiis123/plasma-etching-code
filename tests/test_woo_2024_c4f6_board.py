import csv
import json
from pathlib import Path

import pytest

from scripts.digitize_woo_2024_c4f6_board import (
    AUDIT_PATH,
    CSV_PATH,
    MANIFEST_PATH,
    _payloads,
)


def test_committed_woo_outputs_are_exact_replay():
    csv_payload, manifest_payload, audit_payload, _, _ = _payloads()

    assert CSV_PATH.read_text(encoding="utf-8") == csv_payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_payload
    assert AUDIT_PATH.read_text(encoding="utf-8") == audit_payload


def test_figure4_1_keeps_all_original_pixel_rate_points():
    with CSV_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 10
    assert {row["material"] for row in rows} == {"SiO2", "ACL"}
    assert {
        float(row["c4f6_fraction_of_cf4_plus_c4f6_percent"])
        for row in rows
    } == {37.5, 43.75, 50.0, 56.25, 62.5}

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    check = manifest["digitization"]["diagnostics"]
    assert check["maximum_abs_endpoint_difference_nm_min"] <= 0.2


def test_source_reporting_conflicts_cannot_silently_become_constants():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    consistency = audit["source_internal_consistency"]

    assert not consistency["ion_current_percentage_consistent"]
    assert consistency["ion_current_arithmetic_percent_increase"] == pytest.approx(
        212.79761904761907
    )
    assert not consistency["power_sweep_fraction_consistent"]
    assert consistency["power_sweep_fraction_from_printed_flows_percent"] == 62.5


def test_equal_depth_sem_sweep_is_not_a_blind_depth_pass():
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    depth = audit["feature_depth_classification"]
    boundary = audit["boundary_identifiability"]

    assert depth["exposure_time_adjusted_to_equalize_depth"]
    assert not depth["value_blind_held_out_depth"]
    assert not depth["may_calibrate_simulation_time_from_depth"]
    assert not depth["formal_absolute_feature_depth_pass"]
    assert not boundary["species_resolved_ion_flux_measured"]
    assert not boundary["iead_measured"]
    assert not boundary["absolute_neutral_flux_measured"]
    assert not boundary["knobs_to_feature_depth_identified"]
