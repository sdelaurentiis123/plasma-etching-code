import pytest

from scripts.closeout_zhu_npg80_profile_caches import summarize_profiles


def _profile(depth, width, *, exhausted=False):
    return {
        "profile": {
            "etched_depth_nm": depth,
            "top_cd_nm": width - 10.0,
            "middle_cd_nm": width,
            "bottom_cd_nm": width + 5.0,
            "sidewall_angle_from_wafer_deg": 87.0,
            "bow_nm": 3.0,
        },
        "cr_mask": {
            "center_remaining_thickness_nm": 0.0 if exhausted else 5.0,
            "mask_exhausted_at_center": exhausted,
        },
        "terminal_reason": "requested_duration",
        "maximum_transport_relative_particle_balance_error": 1e-15,
        "maximum_state_remap_relative_conservation_residual": 2e-15,
    }


def test_complete_slice_summary_preserves_physical_ranges():
    result = summarize_profiles(
        105.0,
        [_profile(680.0, 100.0), _profile(685.0, 110.0, exhausted=True)],
    )

    assert result["width_nm"] == 105.0
    assert result["profile_count"] == 2
    assert result["etch_depth_nm"] == [680.0, 685.0]
    assert result["middle_cd_nm"] == [100.0, 110.0]
    assert result["cr_center_exhausted_fraction"] == pytest.approx(0.5)
