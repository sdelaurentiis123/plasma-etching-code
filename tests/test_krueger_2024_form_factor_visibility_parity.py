from types import SimpleNamespace

import numpy as np

from scripts.krueger_2024_form_factor_visibility_parity import compare_visibility_events


def _reference(face, termination, wraps=None):
    face = np.asarray(face, dtype=int)
    return SimpleNamespace(
        hit_face=face,
        termination=np.asarray(termination, dtype=np.int8),
        wrap_count=np.zeros(len(face), dtype=int) if wraps is None else np.asarray(wraps, dtype=int))


def test_visibility_comparison_passes_only_exact_events():
    result = compare_visibility_events(
        [0, 0, 1, 1], [1, -1, 0, -1],
        _reference([1, -1, 0, -1], [1, 2, 1, 2]), [1.0, 2.0])

    assert result["all_gates_pass"]
    assert result["any_event_mismatch_count"] == 0
    assert result["maximum_source_row_total_variation"] == 0.0


def test_visibility_comparison_separates_false_escape_target_and_exhaustion():
    result = compare_visibility_events(
        [0, 0, 1, 1], [-1, 0, 0, -1],
        _reference([1, 1, 0, -1], [1, 1, 1, 3], [0, 0, 0, 4]), [1.0, 2.0])

    assert not result["all_gates_pass"]
    assert result["false_fast_escape_count"] == 1
    assert result["target_face_mismatch_count"] == 1
    assert result["reference_wrap_exhaustion_count"] == 1
    assert result["any_event_mismatch_count"] == 3
    assert result["maximum_source_row_total_variation"] == 1.0
