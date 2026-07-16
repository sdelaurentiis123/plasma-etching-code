import numpy as np
import pytest

from petch.feature_step_3d import FeatureGeometry3D
from petch.profile_observables_3d import (
    measure_feature_centerline_3d,
    measure_feature_centerline_ensemble_3d,
    measure_trench_profile_ensemble_3d,
    measure_trench_profile_observables_3d,
)


def _tilted_hole(slope):
    dx = 0.05
    shape = (41, 41, 21)
    x, y, z = (np.arange(size) * dx for size in shape)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    center_x = 1.0 + float(slope) * (1.0 - Z)
    radius = 0.3
    phi = np.sqrt((X - center_x) ** 2 + (Y - 1.0) ** 2) - radius
    material = np.where(phi > 0.0, 1, 0)
    return FeatureGeometry3D(phi, material, dx, 1e-6)


def _notched_bowed_trench(dx, *, bow=0.1, left_notch=0.12, right_notch=0.04):
    shape = tuple(int(round(length / dx)) + 1 for length in (2.0, 0.4, 1.2))
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    depth = 1.0 - Z
    bow_shape = np.maximum(0.0, 1.0 - ((depth - 0.53) / 0.17) ** 2)
    notch_shape = np.maximum(0.0, 1.0 - ((depth - 0.87) / 0.09) ** 2)
    left = 0.7 - bow * bow_shape - left_notch * notch_shape
    right = 1.3 + bow * bow_shape + right_notch * notch_shape
    phi = np.maximum(left - X, X - right)
    return FeatureGeometry3D(phi, np.where(phi >= 0.0, 1, 0), dx, 1e-6)


def _trench_contract():
    return dict(
        lateral_bounds_m=(0.3e-6, 1.7e-6, 0.0, 0.4e-6),
        opening_center_x_m=1.0e-6, feature_top_z_m=1.0e-6,
        reference_depth_interval_m=(0.0, 0.1e-6),
        bow_depth_interval_m=(0.35e-6, 0.7e-6),
        notch_depth_interval_m=(0.78e-6, 0.98e-6),
        minimum_longitudinal_rows=3)


def test_centerline_recovers_manufactured_tilt_and_onset_ar():
    result = measure_feature_centerline_3d(
        _tilted_hole(0.2), lateral_bounds_m=(0.5e-6, 1.5e-6, 0.5e-6, 1.5e-6),
        opening_width_m=0.6e-6, onset_displacement_m=0.1e-6,
        reference_slice_count=1)

    deepest = np.flatnonzero(result.valid)[-1]
    assert result.displacement_xy_m[deepest, 0] == pytest.approx(0.2e-6, abs=0.03e-6)
    assert abs(result.displacement_xy_m[deepest, 1]) < 0.01e-6
    assert result.maximum_lateral_displacement_m == pytest.approx(0.2e-6, abs=0.03e-6)
    assert result.onset_aspect_ratio == pytest.approx(0.5 / 0.6, abs=0.12)


def test_symmetric_tilt_ensemble_has_zero_mean_and_nonzero_variance():
    result = measure_feature_centerline_ensemble_3d(
        (_tilted_hole(0.2), _tilted_hole(-0.2)),
        lateral_bounds_m=(0.5e-6, 1.5e-6, 0.5e-6, 1.5e-6),
        opening_width_m=0.6e-6, onset_displacement_m=0.1e-6,
        reference_slice_count=1)

    valid = result.valid_realization_count == 2
    assert np.max(np.abs(result.mean_displacement_xy_m[valid])) < 0.01e-6
    assert result.has_nonzero_twist_variance
    assert result.maximum_systematic_z_score < 0.5
    assert not result.statistical_claim_ready

    low, high = result.confidence_interval_displacement_xy_m(1.96)
    assert low.shape == result.mean_displacement_xy_m.shape
    assert high.shape == result.mean_displacement_xy_m.shape


def test_trench_observables_recover_manufactured_notch_bow_and_refine():
    coarse = measure_trench_profile_observables_3d(
        _notched_bowed_trench(0.1), **_trench_contract())
    fine = measure_trench_profile_observables_3d(
        _notched_bowed_trench(0.05), **_trench_contract())

    assert fine.maximum_left_notch_depth_m == pytest.approx(0.12e-6, abs=0.01e-6)
    assert fine.maximum_right_notch_depth_m == pytest.approx(0.04e-6, abs=0.01e-6)
    assert fine.notch_asymmetry_m == pytest.approx(0.08e-6, abs=0.015e-6)
    assert fine.maximum_bow_expansion_m == pytest.approx(0.2e-6, abs=0.01e-6)
    fine_error = abs(fine.maximum_bow_expansion_m - 0.2e-6)
    coarse_error = abs(coarse.maximum_bow_expansion_m - 0.2e-6)
    assert fine_error <= coarse_error + 1e-15


def test_trench_observable_ensemble_reports_declared_uncertainty():
    result = measure_trench_profile_ensemble_3d(
        (_notched_bowed_trench(0.05, left_notch=0.10),
         _notched_bowed_trench(0.05, left_notch=0.14)),
        confidence_multiplier=1.96, **_trench_contract())

    assert result.left_notch_depth.mean_m == pytest.approx(0.12e-6, abs=0.01e-6)
    assert result.left_notch_depth.realization_count == 2
    assert result.left_notch_depth.confidence_halfwidth_m > 0.0
    assert result.maximum_bow_expansion.mean_m == pytest.approx(0.2e-6, abs=0.01e-6)
