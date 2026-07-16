from pathlib import Path

import numpy as np
import pytest

from petch.hwang_giapis_notch_validation_2d import (
    HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256,
    load_hwang_giapis_1997_fig13_profile,
    score_hwang_giapis_1997_fig13_profile,
)


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (
    ROOT / "data" / "experimental" / "hwang_giapis_1997"
    / "fig13_notch_profile.csv")


def test_fig13_profile_is_checksum_verified_and_pixel_replayable():
    observations = load_hwang_giapis_1997_fig13_profile(PROFILE)

    assert len(observations) == 22
    assert max(item.notch_depth_um for item in observations) == pytest.approx(
        0.215688259)
    assert max(abs(item.replayed_notch_depth_um - item.notch_depth_um)
               for item in observations) < 5e-6
    assert max(abs(
        item.replayed_height_above_oxide_um - item.height_above_oxide_um)
        for item in observations) < 5e-6


def test_fig13_profile_refuses_the_wrong_checksum():
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_hwang_giapis_1997_fig13_profile(
            PROFILE, expected_sha256="0" * len(
                HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256))


def test_fig13_score_never_claims_strict_validation_without_measurement_error():
    observations = load_hwang_giapis_1997_fig13_profile(PROFILE)
    height = (np.arange(60) + 0.5) * 0.005
    experimental_height = np.asarray([
        item.height_above_oxide_um for item in observations])
    experimental_depth = np.asarray([
        item.notch_depth_um for item in observations])
    prediction = np.interp(
        height, np.sort(experimental_height),
        experimental_depth[np.argsort(experimental_height)])

    score = score_hwang_giapis_1997_fig13_profile(
        observations, prediction)

    assert score.rmse_um < 0.02
    assert score.strict_validation_pass is False
    assert "DEVELOPMENT_REPLAY_ONLY" in score.claim_status
