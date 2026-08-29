from scripts.summarize_zhu_npg80_gds_square_profiles import _csv_text, build


def _profile(width, depth, bottom, *, exhausted=False):
    return {
        "width_nm": width,
        "profile": {
            "etched_depth_nm": depth,
            "top_cd_nm": width - 10.0,
            "middle_cd_nm": width,
            "bottom_cd_nm": bottom,
            "sidewall_angle_from_wafer_deg": 87.0,
            "bow_nm": 5.0,
        },
        "cr_mask": {
            "center_remaining_thickness_nm": 0.0 if exhausted else 10.0,
            "mask_exhausted_at_center": exhausted,
            "material_layer_retired": exhausted,
        },
        "terminal_reason": "requested_duration",
        "maximum_transport_relative_particle_balance_error": 1e-15,
        "maximum_state_remap_relative_conservation_residual": 2e-15,
    }


def test_exact_gds_summary_reduces_profile_envelopes_without_target_data():
    audit = {
        "schema": "petch.zhu-npg80-exact-gds-square-profile-board.v1",
        "condition_id": "condition",
        "smoke_only": False,
        "target_sem_used": False,
        "target_depth_used": False,
        "geometry": {
            "gds_sha256": "abc",
            "pitch_nm": 350.0,
            "square_width_nm": [105.0],
            "mask_polarity_assumption": "conditional islands",
            "mask_polarity_confirmed_by_operator": False,
        },
        "claim_boundary": {"unique_absolute_oxford_profile_predicted": False},
        "profiles": [
            _profile(105.0, 600.0, 110.0),
            _profile(105.0, 680.0, 120.0, exhausted=True),
        ],
    }

    summary = build(audit)
    row = summary["rows"][0]

    assert summary["target_sem_used"] is False
    assert row["etch_depth_nm"] == [600.0, 680.0]
    assert row["bottom_cd_nm"] == [110.0, 120.0]
    assert row["cr_center_exhausted_fraction"] == 0.5
    assert row["profile_count"] == 2
    assert "etch_depth_min_nm" in _csv_text(summary)
