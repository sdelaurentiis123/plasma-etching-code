import numpy as np
import pytest

from petch.extruded_mask_transport import (
    AngularOrdinate,
    cosine_flux_hemisphere_ordinates,
    direct_extruded_mask_transmission,
    gaussian_transverse_angle_ordinates,
)


def test_angular_quadratures_are_normalized_and_symmetric():
    ion = gaussian_transverse_angle_ordinates(
        np.deg2rad(1.5), order_per_component=7)
    neutral = cosine_flux_hemisphere_ordinates(
        polar_cosine_order=8, azimuth_count=16)
    for rule in (ion, neutral):
        assert sum(item.weight for item in rule) == pytest.approx(1.0, abs=2e-14)
        assert sum(item.weight * item.tangent_x for item in rule) == pytest.approx(
            0.0, abs=2e-13)
        assert sum(item.weight * item.tangent_y for item in rule) == pytest.approx(
            0.0, abs=2e-13)


def test_open_cell_transmits_every_ordinate_exactly():
    opening = np.ones((17, 13), dtype=bool)
    result = direct_extruded_mask_transmission(
        opening,
        mask_height=10.0,
        grid_spacing=0.5,
        ordinates=gaussian_transverse_angle_ordinates(
            np.deg2rad(3.0), order_per_component=5),
    )
    assert np.allclose(
        result.transmission, np.ones_like(opening, dtype=float),
        rtol=0.0, atol=2e-14)
    assert not result.transmission.flags.writeable


def test_normal_incidence_reproduces_the_opening_mask():
    opening = np.ones((19, 15), dtype=bool)
    opening[4:7, 2:13] = False
    result = direct_extruded_mask_transmission(
        opening,
        mask_height=30.0,
        grid_spacing=0.5,
        ordinates=(AngularOrdinate(0.0, 0.0, 1.0),),
    )
    assert np.array_equal(result.transmission, opening.astype(float))


def test_single_slit_matches_geometric_overlap_fraction():
    # Periodic x cell with a 20-cell open slit.  A characteristic displaced by
    # four cells over the mask survives on 16 of the 20 floor cells.
    opening = np.zeros((40, 5), dtype=bool)
    opening[10:30, :] = True
    result = direct_extruded_mask_transmission(
        opening,
        mask_height=8.0,
        grid_spacing=1.0,
        ordinates=(AngularOrdinate(0.5, 0.0, 1.0),),
        subdivisions_per_crossed_cell=2.0,
    )
    assert np.count_nonzero(result.transmission[:, 2]) == 16
    assert np.mean(result.transmission[opening]) == pytest.approx(0.8)


def test_nonperiodic_exterior_is_solid():
    opening = np.ones((11, 11), dtype=bool)
    result = direct_extruded_mask_transmission(
        opening,
        mask_height=5.0,
        grid_spacing=1.0,
        ordinates=(AngularOrdinate(1.0, 0.0, 1.0),),
        periodic_lateral=False,
    )
    assert np.all(result.transmission[:6] == 1.0)
    assert np.all(result.transmission[6:] == 0.0)


@pytest.mark.parametrize("sigma", [-1.0, np.nan])
def test_invalid_angular_width_refuses(sigma):
    with pytest.raises(ValueError):
        gaussian_transverse_angle_ordinates(sigma)
