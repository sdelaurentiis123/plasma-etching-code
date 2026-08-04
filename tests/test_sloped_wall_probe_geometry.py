"""Gates for the straight-taper probe geometry (check (c) of the wall-slope study).

The Gaussian constriction in ``make_necked_trench_geometry_3d`` pins the wall
vertical at its own minimum, so an aperture sweep cannot vary the wall angle
independently.  ``make_sloped_wall_geometry_3d`` prescribes the angle instead
and pins the aperture at the measurement-band midpoint, which is what makes
net-velocity-versus-angle a controlled measurement
(RESULTS_WALL_SLOPE_FALSIFICATION_2026-08-04.md).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from mouth_equilibrium_probe import (  # noqa: E402
    ANGLE_BAND_HI_UM,
    ANGLE_BAND_LO_UM,
    ANGLE_BAND_MID_UM,
    MASK_THICKNESS_UM,
    NECK_DEPTH_UM,
    SUBSTRATE_TOP_UM,
    make_sloped_wall_geometry_3d,
    measure_aperture_profile,
)


def _band_fit(geometry):
    """Realised (aperture at band midpoint, wall angle) from the level set."""
    depths, apertures = measure_aperture_profile(geometry)
    apex_depth_nm = NECK_DEPTH_UM * 1e3
    band = ((depths >= apex_depth_nm - ANGLE_BAND_HI_UM * 1e3)
            & (depths <= apex_depth_nm - ANGLE_BAND_LO_UM * 1e3))
    assert band.sum() >= 2, "measurement band must contain at least two rows"
    slope, intercept = np.polyfit(depths[band], apertures[band] / 2.0, 1)
    mid_depth = apex_depth_nm - ANGLE_BAND_MID_UM * 1e3
    aperture_mid = 2.0 * (slope * mid_depth + intercept)
    return aperture_mid, float(np.degrees(np.arctan(abs(slope))))


@pytest.mark.parametrize("angle_deg", [1.0, 3.0, 8.0, 12.0])
def test_prescribed_wall_angle_is_realised(angle_deg):
    """The wedge wall must carry the angle it was asked for."""
    geometry = make_sloped_wall_geometry_3d(
        wall_angle_deg=angle_deg, aperture_nm=45.0, dx=0.01)
    _, realised = _band_fit(geometry)
    assert realised == pytest.approx(angle_deg, abs=0.1)


@pytest.mark.parametrize("angle_deg", [1.0, 3.0, 8.0, 12.0])
def test_band_aperture_is_pinned_across_angles(angle_deg):
    """Angle must vary at FIXED local aperture, else the sweep confounds the two."""
    geometry = make_sloped_wall_geometry_3d(
        wall_angle_deg=angle_deg, aperture_nm=45.0, dx=0.01)
    aperture_mid, _ = _band_fit(geometry)
    assert aperture_mid == pytest.approx(45.0, abs=0.5)


def test_wall_is_straight_not_stationary_at_the_band():
    """Contrast with the Gaussian: a straight wedge has non-zero slope in-band.

    A smooth minimum has zero slope by construction, which is why the aperture
    sweep read alpha ~ 0 at every prescribed neck.
    """
    from mouth_equilibrium_probe import make_necked_trench_geometry_3d

    wedge = make_sloped_wall_geometry_3d(
        wall_angle_deg=8.0, aperture_nm=45.0, dx=0.01)
    _, wedge_angle = _band_fit(wedge)
    gaussian = make_necked_trench_geometry_3d(neck_width_um=0.045, dx=0.01)
    depths, apertures = measure_aperture_profile(gaussian)
    apex_depth_nm = NECK_DEPTH_UM * 1e3
    near_apex = np.abs(depths - apex_depth_nm) <= 10.0
    slope = np.polyfit(depths[near_apex], apertures[near_apex] / 2.0, 1)[0]
    gaussian_angle = float(np.degrees(np.arctan(abs(slope))))
    assert wedge_angle > 5.0
    assert gaussian_angle < 1.0


def test_steep_angle_for_narrow_aperture_is_refused():
    """The apex would invert; the builder must refuse rather than emit garbage."""
    with pytest.raises(ValueError):
        make_sloped_wall_geometry_3d(
            wall_angle_deg=45.0, aperture_nm=20.0, dx=0.01)


def test_mask_material_spans_the_measurement_band():
    """Band faces must be mask, not substrate, or the probe reads the wrong film."""
    geometry = make_sloped_wall_geometry_3d(
        wall_angle_deg=8.0, aperture_nm=45.0, dx=0.01)
    _, _, z = geometry.coordinate_arrays
    mask_top = SUBSTRATE_TOP_UM + MASK_THICKNESS_UM
    z_apex = mask_top - NECK_DEPTH_UM
    in_band = (z >= z_apex + ANGLE_BAND_LO_UM) & (z <= z_apex + ANGLE_BAND_HI_UM)
    assert in_band.any()
    band_materials = geometry.material_id[:, :, in_band]
    assert np.any(band_materials == 2)
    assert not np.any(band_materials == 1)
