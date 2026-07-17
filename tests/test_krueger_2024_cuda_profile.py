import cProfile
from pathlib import Path
import pstats
import sys

import pytest


SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import krueger_2024_cuda_profile as profile


def _profile_target():
    return sum(index * index for index in range(25))


def test_top_profile_rows_are_sorted_and_json_ready():
    profiler = cProfile.Profile()
    profiler.runcall(_profile_target)

    rows = profile.top_profile_rows(pstats.Stats(profiler), limit=6)

    assert rows
    assert all(
        left["cumulative_time_s"] >= right["cumulative_time_s"]
        for left, right in zip(rows[:-1], rows[1:]))
    assert all(isinstance(row["path"], str) for row in rows)


def test_profiler_cli_refuses_an_unbounded_physical_run():
    with pytest.raises(SystemExit):
        profile.parse_args((
            "--positive-warmup-steps", "1",
            "--profile-steps", "4", "--step-duration-s", "0.025"))

    accepted = profile.parse_args((
        "--positive-warmup-steps", "1",
        "--profile-steps", "2", "--step-duration-s", "0.025"))
    assert (
        (accepted.positive_warmup_steps + accepted.profile_steps)
        * accepted.step_duration_s) == pytest.approx(0.075)
