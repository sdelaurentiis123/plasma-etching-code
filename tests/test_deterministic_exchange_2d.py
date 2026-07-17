import numpy as np
import pytest

from petch.deterministic_exchange_2d import (
    build_deterministic_line_exchange_2d,
    unobstructed_crossed_string_exchange_2d,
)
from petch.neutral_radiosity_3d import solve_diffuse_neutral_radiosity_3d


def test_crossed_string_parallel_equal_segments_matches_analytic_result():
    lower = np.array([[0.0, 0.0], [1.0, 0.0]])
    upper = np.array([[0.0, 1.0], [1.0, 1.0]])
    exchange = unobstructed_crossed_string_exchange_2d(lower, upper)
    assert exchange == pytest.approx(np.sqrt(2.0) - 1.0, abs=2e-15)
    # Endpoint ordering is not a physical input.
    assert unobstructed_crossed_string_exchange_2d(
        lower[::-1], upper) == pytest.approx(exchange, abs=2e-15)


def test_open_parallel_segments_are_reciprocal_and_row_closing():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[-0.5, 1.0], [1.5, 1.0]],
    ])
    normals = np.array([[0.0, 1.0], [0.0, -1.0]])
    result = build_deterministic_line_exchange_2d(segments, normals)
    length = result.segment_length
    assert result.exchange_length[0, 1] == pytest.approx(
        result.exchange_length[1, 0], abs=1e-15)
    assert length[0] * result.transfer_fraction[0, 1] == pytest.approx(
        length[1] * result.transfer_fraction[1, 0], abs=1e-15)
    assert np.allclose(
        result.transfer_fraction.sum(axis=1) + result.escape_fraction, 1.0,
        rtol=0.0, atol=1e-15)


def test_back_to_back_segments_exchange_only_with_open_boundary():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 1.0], [1.0, 1.0]],
    ])
    result = build_deterministic_line_exchange_2d(
        segments, np.array([[0.0, -1.0], [0.0, 1.0]]))
    assert np.array_equal(result.transfer_fraction, np.zeros((2, 2)))
    assert np.array_equal(result.escape_fraction, np.ones(2))


def test_opaque_middle_surface_blocks_parallel_plate_exchange():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 2.0], [1.0, 2.0]],
        [[-0.2, 1.0], [1.2, 1.0]],
    ])
    normals = np.array([
        [0.0, 1.0], [0.0, -1.0], [0.0, -1.0],
    ])
    result = build_deterministic_line_exchange_2d(
        segments, normals, relative_tolerance=1e-7,
        minimum_refinement_level=3, maximum_refinement_level=7)
    assert result.transfer_fraction[0, 1] == pytest.approx(0.0, abs=1e-14)
    assert result.refinement_level[0, 1] >= 3


def test_partial_blocking_converges_under_tighter_refinement():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 2.0], [1.0, 2.0]],
        [[0.0, 1.0], [0.45, 1.0]],
    ])
    normals = np.array([
        [0.0, 1.0], [0.0, -1.0], [0.0, -1.0],
    ])
    coarse = build_deterministic_line_exchange_2d(
        segments, normals, relative_tolerance=2e-3,
        minimum_refinement_level=2, maximum_refinement_level=12)
    fine = build_deterministic_line_exchange_2d(
        segments, normals, relative_tolerance=2e-4,
        minimum_refinement_level=3, maximum_refinement_level=15)
    assert 0.0 < fine.transfer_fraction[0, 1]
    assert fine.transfer_fraction[0, 1] < (
        unobstructed_crossed_string_exchange_2d(segments[0], segments[1])
        / fine.segment_length[0])
    assert coarse.transfer_fraction[0, 1] == pytest.approx(
        fine.transfer_fraction[0, 1], rel=8e-3, abs=2e-5)


def test_common_radiosity_adapter_preserves_particle_balance():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 1.0], [1.0, 1.0]],
    ])
    normals = np.array([[0.0, 1.0], [0.0, -1.0]])
    exchange = build_deterministic_line_exchange_2d(segments, normals)
    factors = exchange.as_diffuse_form_factors_3d()
    solution = solve_diffuse_neutral_radiosity_3d(
        np.array([1.0, 0.0]), exchange.segment_length,
        factors.source_face, factors.target_face, factors.transfer_fraction,
        factors.escape_fraction, np.array([0.25, 0.25]))
    assert solution.iterations_converged
    assert solution.relative_balance_error < 1e-13
    assert solution.source_rate_s == pytest.approx(
        solution.reacted_rate_s + solution.escaped_rate_s, rel=2e-14)
