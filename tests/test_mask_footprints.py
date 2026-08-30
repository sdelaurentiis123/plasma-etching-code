import numpy as np
import pytest

from petch.mask_footprints import (
    centered_cross_footprint_levelset,
    centered_inverse_square_hole_footprint_levelset,
    centered_rectangle_footprint_levelset,
    centered_square_footprint_levelset,
    polygon_union_footprint_levelset,
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


def test_polygon_union_levelset_has_correct_sign_and_periodic_distance():
    rectangle = np.asarray([
        [0.0, 0.75],
        [0.5, 0.75],
        [0.5, 1.25],
        [0.0, 1.25],
    ])
    field = polygon_union_footprint_levelset(
        cell_width=2.0,
        cell_length=2.0,
        dx=0.25,
        polygons=[rectangle],
    )

    assert field.shape == (9, 9)
    assert np.allclose(field[0, :], field[-1, :])
    assert np.allclose(field[:, 0], field[:, -1])
    assert field[1, 4] > 0.0
    assert field[4, 4] < 0.0
    assert field[7, 4] < 0.0
    assert field[7, 4] == pytest.approx(-0.25)
    assert field[0, 3] == pytest.approx(0.0)


def test_polygon_union_combines_disconnected_mask_solids():
    first = np.asarray([[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]])
    second = np.asarray([[1.25, 1.25], [1.75, 1.25], [1.75, 1.75], [1.25, 1.75]])
    field = polygon_union_footprint_levelset(
        cell_width=2.0,
        cell_length=2.0,
        dx=0.25,
        polygons=[first, second],
    )

    assert field[2, 2] > 0.0
    assert field[6, 6] > 0.0
    assert field[4, 4] < 0.0


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"cell_width": 2.1, "cell_length": 2.0, "dx": 0.25, "polygons": [[[0, 0], [1, 0], [1, 1]]]}, "integer multiples"),
        ({"cell_width": 2.0, "cell_length": 2.0, "dx": 0.25, "polygons": []}, "at least one"),
        ({"cell_width": 2.0, "cell_length": 2.0, "dx": 0.25, "polygons": [[[0, 0], [3, 0], [1, 1]]]}, "outside"),
        ({"cell_width": 2.0, "cell_length": 2.0, "dx": 0.25, "polygons": [[[0, 0], [1, 0], [2, 0]]]}, "zero-area"),
    ],
)
def test_polygon_union_rejects_ambiguous_geometry(kwargs, match):
    with pytest.raises(ValueError, match=match):
        polygon_union_footprint_levelset(**kwargs)
