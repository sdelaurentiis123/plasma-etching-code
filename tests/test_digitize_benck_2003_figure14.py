import csv
from io import StringIO

import pytest

from scripts.digitize_benck_2003_figure14 import (
    CF2_CF_AT_TICKS,
    LEFT_X_PX,
    LEFT_Y_TICKS_PX,
    RIGHT_X_PX,
    RIGHT_Y_TICKS_PX,
    csv_text,
    ratio_at_pixel,
    verify_committed_files,
)


def test_committed_benck_figure14_board_replays_pixels():
    verify_committed_files()
    for pixel, expected in zip(LEFT_Y_TICKS_PX, CF2_CF_AT_TICKS):
        assert ratio_at_pixel(LEFT_X_PX, pixel) == pytest.approx(
            expected, abs=0.05
        )
    for pixel, expected in zip(RIGHT_Y_TICKS_PX, CF2_CF_AT_TICKS):
        assert ratio_at_pixel(RIGHT_X_PX, pixel) == pytest.approx(
            expected, abs=0.05
        )


def test_benck_figure14_retains_only_co_conditioned_five_sccm_series():
    rows = list(csv.DictReader(StringIO(csv_text())))

    assert [int(row["c4f6_feed_percent"]) for row in rows] == [25, 50, 75, 100]
    assert {float(row["flow_sccm"]) for row in rows} == {5.0}
    assert {float(row["pressure_Pa"]) for row in rows} == {1.33}
    ratios = [float(row["neutral_CF2_to_CF_density_ratio"]) for row in rows]
    assert ratios == sorted(ratios)
    assert 10.0 < ratios[0] < 12.0
    assert 16.0 < ratios[-1] < 18.0
