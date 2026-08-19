import csv

from scripts import digitize_benck_2003_figure13 as figure13


def test_figure13_replays_committed_voltage_board():
    figure13.verify_committed_files()
    rows = list(csv.DictReader(figure13.CSV_PATH.open(encoding="utf-8")))

    assert len(rows) == 12
    assert {row["panel"] for row in rows} == {"13a", "13b"}
    assert {row["power_W"] for row in rows} == {"200.0"}
    assert all(float(row["plasma_potential_peak_to_peak_V"]) > 0.0 for row in rows)
    assert all(float(row["digitization_absolute_bound_V"]) == 0.6 for row in rows)


def test_figure13_preserves_feed_and_pressure_trends():
    rows = list(csv.DictReader(figure13.CSV_PATH.open(encoding="utf-8")))
    feed = [
        row for row in rows
        if row["panel"] == "13a" and row["series"] == "5_sccm"
    ]
    feed.sort(key=lambda row: float(row["c4f6_feed_percent"]))
    assert [int(row["c4f6_feed_percent"]) for row in feed] == [25, 50, 75, 100]
    assert all(
        float(left["plasma_potential_peak_to_peak_V"])
        < float(right["plasma_potential_peak_to_peak_V"])
        for left, right in zip(feed, feed[1:])
    )

    pure = [
        row for row in rows
        if row["panel"] == "13b" and row["series"] == "100_percent_C4F6"
    ]
    pure.sort(key=lambda row: float(row["pressure_Pa"]))
    pure_values = [float(row["plasma_potential_peak_to_peak_V"]) for row in pure]
    assert pure_values[1] > pure_values[0]
    assert pure_values[1] > pure_values[2]

    mixture = [
        row for row in rows
        if row["panel"] == "13b" and row["series"] == "50_percent_C4F6"
    ]
    mixture.sort(key=lambda row: float(row["pressure_Pa"]))
    assert all(
        float(left["plasma_potential_peak_to_peak_V"])
        < float(right["plasma_potential_peak_to_peak_V"])
        for left, right in zip(mixture, mixture[1:])
    )


def test_figure13_affine_tick_replay_is_subpixel_consistent():
    manifest = figure13.manifest("placeholder")
    for panel in ("panel_13a", "panel_13b"):
        calibration = manifest["pixel_calibration"][panel]
        for replay in (
            calibration["left_tick_replay_V"],
            calibration["right_tick_replay_V"],
        ):
            assert max(
                abs(observed - expected)
                for observed, expected in zip(replay, figure13.VPP_AT_TICKS)
            ) < 0.08
