from scripts.hwang_giapis_1997_fig13_validation import (
    _continuation_compatibility_hash,
    _potential_guard_transition_allowed,
)


def test_fig13_continuation_compatibility_ignores_sampling_controls():
    first = {
        "AR": 2.6,
        "W": 25,
        "mouth": 120,
        "source_launch_plane": "sheath_lower_boundary",
        "n_per_iter": 1000,
        "n_iter": 2000,
        "seed": 1,
        "relax": 8.0,
        "trace_step_cap_factor": 40.0,
        "trace_adaptive_horizon": True,
        "trace_emergency_step_cap_factor": 1280.0,
        "trace_relative_tail_tolerance": 0.0,
        "allow_trajectory_truncation": False,
        "potential_guard_policy": "legacy_clip",
        "potential_emergency_abs_v": 250.0,
    }
    second = {
        **first,
        "n_per_iter": 2000,
        "n_iter": 500,
        "seed": 2,
        "relax": 4.0,
        "trace_step_cap_factor": 80.0,
        "trace_adaptive_horizon": False,
        "trace_emergency_step_cap_factor": 2560.0,
        "trace_relative_tail_tolerance": 1.0e-4,
        "allow_trajectory_truncation": True,
        "potential_guard_policy": "source_faithful_refuse",
        "potential_emergency_abs_v": 500.0,
    }
    assert (
        _continuation_compatibility_hash(first)
        == _continuation_compatibility_hash(second)
    )


def test_fig13_continuation_allows_only_legacy_to_source_faithful_guard_upgrade():
    assert _potential_guard_transition_allowed(
        "legacy_clip", "source_faithful_refuse")
    assert _potential_guard_transition_allowed(
        "source_faithful_refuse", "source_faithful_refuse")
    assert not _potential_guard_transition_allowed(
        "source_faithful_refuse", "legacy_clip")
    assert not _potential_guard_transition_allowed(
        "unknown", "source_faithful_refuse")


def test_fig13_continuation_compatibility_separates_launch_operators():
    corrected = {
        "AR": 2.6,
        "W": 25,
        "mouth": 120,
        "source_launch_plane": "sheath_lower_boundary",
    }
    legacy_explicit = {
        **corrected,
        "source_launch_plane": "feature_mouth_legacy",
    }
    legacy_implicit = {
        key: value for key, value in corrected.items()
        if key != "source_launch_plane"
    }
    assert (
        _continuation_compatibility_hash(legacy_implicit)
        == _continuation_compatibility_hash(legacy_explicit)
    )
    assert (
        _continuation_compatibility_hash(corrected)
        != _continuation_compatibility_hash(legacy_explicit)
    )
