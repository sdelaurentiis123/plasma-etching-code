from scripts.hwang_giapis_1997_fig13_validation import (
    _continuation_compatibility_hash,
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
    }
    second = {
        **first,
        "n_per_iter": 2000,
        "n_iter": 500,
        "seed": 2,
        "relax": 4.0,
    }
    assert (
        _continuation_compatibility_hash(first)
        == _continuation_compatibility_hash(second)
    )


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
