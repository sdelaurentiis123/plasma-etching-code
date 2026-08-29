import numpy as np

from petch.mask_footprints import (
    centered_cross_footprint_levelset,
    centered_inverse_square_hole_footprint_levelset,
    centered_rectangle_footprint_levelset,
    centered_square_footprint_levelset,
)


def test_square_and_rectangle_fields_match_exact_gds_zero_sets():
    square = centered_square_footprint_levelset(
        pitch=350.0, dx=2.5, square_width=105.0)
    rectangle = centered_rectangle_footprint_levelset(
        cell_width=350.0,
        cell_length=350.0,
        dx=2.5,
        rectangle_width=250.0,
        rectangle_length=105.0,
    )
    center = 70

    assert square[center, center] == 52.5
    assert square[center + 21, center] == 0.0
    assert rectangle[center + 50, center] == 0.0
    assert rectangle[center, center + 21] == 0.0
    assert np.array_equal(square[0], square[-1])
    assert np.array_equal(rectangle[:, 0], rectangle[:, -1])


def test_cross_is_union_of_250_by_105_rectangles():
    cross = centered_cross_footprint_levelset(
        pitch=350.0, dx=2.5, outer_width=250.0, arm_width=105.0)
    center = 70

    assert cross[center, center] == 52.5
    assert cross[center + 50, center] == 0.0
    assert cross[center, center + 50] == 0.0
    assert cross[center + 30, center + 30] < 0.0
    assert np.array_equal(cross[0], cross[-1])
    assert np.array_equal(cross[:, 0], cross[:, -1])


def test_inverse_hole_is_positive_blanket_outside_exact_square_opening():
    inverse = centered_inverse_square_hole_footprint_levelset(
        pitch=350.0, dx=2.5, opening_width=105.0)
    center = 70

    assert inverse[center, center] == -52.5
    assert inverse[center + 21, center] == 0.0
    assert inverse[0, 0] > 0.0
    assert np.array_equal(inverse[0], inverse[-1])
    assert np.array_equal(inverse[:, 0], inverse[:, -1])
