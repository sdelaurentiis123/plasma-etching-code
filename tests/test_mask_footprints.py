import numpy as np
import pytest

from petch.mask_footprints import polygon_union_footprint_levelset


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
