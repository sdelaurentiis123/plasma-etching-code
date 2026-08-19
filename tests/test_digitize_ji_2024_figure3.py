import csv

from scripts import digitize_ji_2024_figure3 as figure3


def _series(observable):
    rows = list(csv.DictReader(figure3.CSV_PATH.open(encoding="utf-8")))
    selected = [row for row in rows if row["observable"] == observable]
    selected.sort(key=lambda row: float(row["rf_power_W"]))
    return selected


def test_figure3_replays_fifteen_checksum_bound_points():
    figure3.verify_committed_files()
    rows = list(csv.DictReader(figure3.CSV_PATH.open(encoding="utf-8")))

    assert len(rows) == 15
    assert {row["panel"] for row in rows} == {"3b1", "3b2", "3b3"}
    assert {float(row["rf_power_W"]) for row in rows} == {
        90.0, 120.0, 150.0, 180.0, 210.0,
    }


def test_figure3_preserves_the_distinct_physical_responses():
    height = [float(row["value"]) for row in _series("upper_triangle_height")]
    radius = [float(row["value"]) for row in _series("upper_tip_corner_radius")]
    gap = [float(row["value"]) for row in _series("interfeature_gap")]

    assert height.index(max(height)) == 3
    assert all(left >= right for left, right in zip(radius, radius[1:]))
    assert all(left > right for left, right in zip(gap, gap[1:]))
    assert gap[0] - gap[-1] > 75.0


def test_figure3_manifest_forbids_oxford_coefficient_transfer():
    manifest = figure3.manifest("placeholder")
    physics = manifest["physics_use"]

    assert physics[
        "removal_only_model_can_reproduce_strict_gap_narrowing"
    ] is False
    assert "positive deposited or retained surface-volume channel" in physics[
        "required_model_topology"
    ]
    assert "an Oxford NPG80 TiO2 coefficient" in physics["not_valid"]
    assert manifest["experiment"]["etch_time_s"] is None
    for calibration in manifest["pixel_calibration"].values():
        assert max(
            abs(observed - expected)
            for observed, expected in zip(
                calibration["tick_replay"], calibration["tick_values"]
            )
        ) < 0.11
