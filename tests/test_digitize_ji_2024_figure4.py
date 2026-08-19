import csv

from scripts import digitize_ji_2024_figure4 as figure4


def _rows():
    return list(csv.DictReader(figure4.CSV_PATH.open(encoding="utf-8")))


def test_figure4_replays_twenty_checksum_bound_values():
    figure4.verify_committed_files()
    rows = _rows()

    assert len(rows) == 5
    assert {row["panel"] for row in rows} == {"4a", "4b", "4c", "4d", "4e"}
    assert {float(row["designed_gap_nm"]) for row in rows} == {
        70.0, 100.0, 350.0, 530.0, 750.0,
    }
    assert all(row["source_measurement_uncertainty_reported"] == "false" for row in rows)


def test_figure4_preserves_observed_loading_response_and_boundary_ambiguity():
    manifest = figure4.manifest("placeholder")
    checks = manifest["derived_checks"]

    assert checks["digitized_100nm_point_already_shifted"] is True
    assert checks["strict_threshold_identified"] is False
    assert checks["empirical_transition_bracket_nm"] == [100.0, 350.0]
    assert checks["loading_gap_h1_below_every_high_gap_h1"] is True
    assert checks["loading_gap_angle_above_every_high_gap_angle"] is True
    assert checks["all_widths_within_15nm"] is True


def test_figure4_forbids_sharp_threshold_or_oxford_coefficient_transfer():
    manifest = figure4.manifest("placeholder")
    implication = manifest["freddie_geometry_implication"]
    boundary = manifest["claim_boundary"]

    assert implication["threshold_transfer_allowed"] is False
    assert implication["member_at_or_below_changed_100nm_point_width_nm"] == [320.0]
    assert implication["member_in_unsampled_transition_interval_width_nm"] == [280.0]
    assert "an Oxford pattern-loading threshold" in boundary["not_valid"]
    assert "an absolute etch-depth or time calibration" in boundary["not_valid"]
