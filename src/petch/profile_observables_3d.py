"""Geometry-native centerline and stochastic twist observables for 3-D features."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .feature_step_3d import FeatureGeometry3D


@dataclass(frozen=True)
class EnsembleScalarEstimate3D:
    """A scalar ensemble estimate with an explicitly declared confidence multiplier."""

    mean_m: float
    sample_standard_deviation_m: float
    standard_error_m: float
    confidence_halfwidth_m: float
    realization_count: int
    confidence_multiplier: float

    def __post_init__(self):
        values = np.asarray([
            self.mean_m, self.sample_standard_deviation_m, self.standard_error_m,
            self.confidence_halfwidth_m, self.confidence_multiplier], dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values[1:] < 0.0)
                or int(self.realization_count) != self.realization_count
                or self.realization_count < 2 or self.confidence_multiplier <= 0.0):
            raise ValueError("invalid scalar ensemble estimate")
        object.__setattr__(self, "realization_count", int(self.realization_count))

    @property
    def confidence_interval_m(self):
        return (
            self.mean_m - self.confidence_halfwidth_m,
            self.mean_m + self.confidence_halfwidth_m,
        )


@dataclass(frozen=True)
class TrenchProfileObservables3D:
    """Geometry-native sidewall, bow, and notch measurements for one line opening.

    Left and right are coordinate labels only.  A benchmark must explicitly map them to inner and
    outer line sidewalls; this object never infers pattern connectivity from a cropped geometry.
    """

    depth_m: np.ndarray
    left_boundary_m: np.ndarray
    right_boundary_m: np.ndarray
    valid: np.ndarray
    reference_depth_interval_m: tuple[float, float]
    bow_depth_interval_m: tuple[float, float]
    notch_depth_interval_m: tuple[float, float]
    reference_left_boundary_m: float
    reference_right_boundary_m: float

    def __post_init__(self):
        depth = np.asarray(self.depth_m, dtype=float).copy()
        left = np.asarray(self.left_boundary_m, dtype=float).copy()
        right = np.asarray(self.right_boundary_m, dtype=float).copy()
        valid = np.asarray(self.valid, dtype=bool).copy()
        intervals = tuple(tuple(float(value) for value in item) for item in (
            self.reference_depth_interval_m,
            self.bow_depth_interval_m,
            self.notch_depth_interval_m,
        ))
        references = np.asarray([
            self.reference_left_boundary_m, self.reference_right_boundary_m], dtype=float)
        if (depth.ndim != 1 or left.shape != depth.shape or right.shape != depth.shape
                or valid.shape != depth.shape or not np.any(valid)
                or np.any(~np.isfinite(depth)) or np.any(depth < 0.0)
                or np.any(~np.isfinite(left[valid])) or np.any(~np.isfinite(right[valid]))
                or np.any(right[valid] <= left[valid]) or np.any(~np.isfinite(references))
                or references[1] <= references[0]
                or any(len(item) != 2 or item[0] < 0.0 or item[1] < item[0]
                       or not np.all(np.isfinite(item)) for item in intervals)):
            raise ValueError("invalid trench-profile observables")
        for interval in intervals:
            if not np.any(valid & (depth >= interval[0]) & (depth <= interval[1])):
                raise ValueError("a declared profile-observable depth interval is unresolved")
        for value in (depth, left, right, valid):
            value.setflags(write=False)
        object.__setattr__(self, "depth_m", depth)
        object.__setattr__(self, "left_boundary_m", left)
        object.__setattr__(self, "right_boundary_m", right)
        object.__setattr__(self, "valid", valid)
        object.__setattr__(self, "reference_depth_interval_m", intervals[0])
        object.__setattr__(self, "bow_depth_interval_m", intervals[1])
        object.__setattr__(self, "notch_depth_interval_m", intervals[2])

    def _depth_mask(self, interval):
        return self.valid & (self.depth_m >= interval[0]) & (self.depth_m <= interval[1])

    @property
    def opening_width_m(self):
        return self.right_boundary_m - self.left_boundary_m

    @property
    def reference_opening_width_m(self):
        return self.reference_right_boundary_m - self.reference_left_boundary_m

    @property
    def left_outward_recession_m(self):
        return self.reference_left_boundary_m - self.left_boundary_m

    @property
    def right_outward_recession_m(self):
        return self.right_boundary_m - self.reference_right_boundary_m

    @property
    def maximum_left_notch_depth_m(self):
        mask = self._depth_mask(self.notch_depth_interval_m)
        return float(max(0.0, np.max(self.left_outward_recession_m[mask])))

    @property
    def maximum_right_notch_depth_m(self):
        mask = self._depth_mask(self.notch_depth_interval_m)
        return float(max(0.0, np.max(self.right_outward_recession_m[mask])))

    @property
    def maximum_notch_depth_m(self):
        return max(self.maximum_left_notch_depth_m, self.maximum_right_notch_depth_m)

    @property
    def notch_asymmetry_m(self):
        return self.maximum_left_notch_depth_m - self.maximum_right_notch_depth_m

    @property
    def maximum_bow_width_m(self):
        mask = self._depth_mask(self.bow_depth_interval_m)
        return float(np.max(self.opening_width_m[mask]))

    @property
    def maximum_bow_expansion_m(self):
        return max(0.0, self.maximum_bow_width_m - self.reference_opening_width_m)


@dataclass(frozen=True)
class TrenchProfileEnsemble3D:
    """Ensemble uncertainty for the profile observables used by notch/bow campaigns."""

    members: tuple[TrenchProfileObservables3D, ...]
    maximum_notch_depth: EnsembleScalarEstimate3D
    left_notch_depth: EnsembleScalarEstimate3D
    right_notch_depth: EnsembleScalarEstimate3D
    notch_asymmetry: EnsembleScalarEstimate3D
    maximum_bow_width: EnsembleScalarEstimate3D
    maximum_bow_expansion: EnsembleScalarEstimate3D

    def __post_init__(self):
        members = tuple(self.members)
        estimates = (
            self.maximum_notch_depth, self.left_notch_depth, self.right_notch_depth,
            self.notch_asymmetry, self.maximum_bow_width, self.maximum_bow_expansion,
        )
        if (len(members) < 2
                or any(not isinstance(item, TrenchProfileObservables3D) for item in members)
                or any(not isinstance(item, EnsembleScalarEstimate3D) for item in estimates)
                or any(item.realization_count != len(members) for item in estimates)):
            raise ValueError("invalid trench-profile ensemble")
        object.__setattr__(self, "members", members)


@dataclass(frozen=True)
class FeatureCenterline3D:
    """Open-feature gas centroid and lateral displacement versus physical depth."""

    depth_m: np.ndarray
    centroid_xy_m: np.ndarray
    displacement_xy_m: np.ndarray
    equivalent_diameter_m: np.ndarray
    valid: np.ndarray
    opening_width_m: float
    onset_displacement_m: float

    def __post_init__(self):
        depth = np.asarray(self.depth_m, dtype=float).copy()
        centroid = np.asarray(self.centroid_xy_m, dtype=float).copy()
        displacement = np.asarray(self.displacement_xy_m, dtype=float).copy()
        diameter = np.asarray(self.equivalent_diameter_m, dtype=float).copy()
        valid = np.asarray(self.valid, dtype=bool).copy()
        if (depth.ndim != 1 or centroid.shape != (len(depth), 2)
                or displacement.shape != centroid.shape or diameter.shape != depth.shape
                or valid.shape != depth.shape or not np.any(valid)
                or np.any(~np.isfinite(depth)) or np.any(depth < 0.0)
                or np.any(~np.isfinite(centroid[valid]))
                or np.any(~np.isfinite(displacement[valid]))
                or np.any(~np.isfinite(diameter[valid])) or np.any(diameter[valid] <= 0.0)
                or not np.isfinite(self.opening_width_m) or self.opening_width_m <= 0.0
                or not np.isfinite(self.onset_displacement_m)
                or self.onset_displacement_m <= 0.0):
            raise ValueError("invalid feature-centerline observable")
        for value in (depth, centroid, displacement, diameter, valid):
            value.setflags(write=False)
        object.__setattr__(self, "depth_m", depth)
        object.__setattr__(self, "centroid_xy_m", centroid)
        object.__setattr__(self, "displacement_xy_m", displacement)
        object.__setattr__(self, "equivalent_diameter_m", diameter)
        object.__setattr__(self, "valid", valid)

    @property
    def lateral_displacement_m(self):
        return np.linalg.norm(self.displacement_xy_m, axis=1)

    @property
    def maximum_lateral_displacement_m(self):
        return float(np.max(self.lateral_displacement_m[self.valid]))

    @property
    def onset_depth_m(self):
        selected = np.flatnonzero(
            self.valid & (self.lateral_displacement_m >= self.onset_displacement_m))
        return None if not len(selected) else float(self.depth_m[selected[0]])

    @property
    def onset_aspect_ratio(self):
        depth = self.onset_depth_m
        return None if depth is None else depth / self.opening_width_m


@dataclass(frozen=True)
class FeatureCenterlineEnsemble3D:
    """Ensemble centerline moments on one common physical depth grid."""

    members: tuple[FeatureCenterline3D, ...]
    mean_displacement_xy_m: np.ndarray
    standard_deviation_displacement_xy_m: np.ndarray
    valid_realization_count: np.ndarray
    maximum_systematic_z_score: float
    onset_aspect_ratio: np.ndarray

    def __post_init__(self):
        members = tuple(self.members)
        mean = np.asarray(self.mean_displacement_xy_m, dtype=float).copy()
        deviation = np.asarray(self.standard_deviation_displacement_xy_m, dtype=float).copy()
        count = np.asarray(self.valid_realization_count, dtype=int).copy()
        onset = np.asarray(self.onset_aspect_ratio, dtype=float).copy()
        if (len(members) < 2 or any(not isinstance(item, FeatureCenterline3D) for item in members)
                or mean.ndim != 2 or mean.shape[1] != 2 or deviation.shape != mean.shape
                or count.shape != (len(mean),) or onset.shape != (len(members),)
                or np.any(count < 0) or np.any(count > len(members))
                or np.any(~np.isfinite(mean[count > 0]))
                or np.any(~np.isfinite(deviation[count > 1]))
                or np.any(deviation[count > 1] < 0.0)
                or np.any(~np.isfinite(onset) & ~np.isnan(onset))
                or not np.isfinite(self.maximum_systematic_z_score)
                or self.maximum_systematic_z_score < 0.0):
            raise ValueError("invalid feature-centerline ensemble")
        for value in (mean, deviation, count, onset):
            value.setflags(write=False)
        object.__setattr__(self, "members", members)
        object.__setattr__(self, "mean_displacement_xy_m", mean)
        object.__setattr__(self, "standard_deviation_displacement_xy_m", deviation)
        object.__setattr__(self, "valid_realization_count", count)
        object.__setattr__(self, "onset_aspect_ratio", onset)

    @property
    def has_nonzero_twist_variance(self):
        return bool(np.any(self.standard_deviation_displacement_xy_m > 0.0))

    @property
    def statistical_claim_ready(self):
        # N/sample doubling and an independent isotropy threshold remain campaign gates.
        return False

    @property
    def standard_error_displacement_xy_m(self):
        return np.divide(
            self.standard_deviation_displacement_xy_m,
            np.sqrt(self.valid_realization_count)[:, None],
            out=np.full_like(self.standard_deviation_displacement_xy_m, np.nan),
            where=self.valid_realization_count[:, None] > 1)

    def confidence_interval_displacement_xy_m(self, confidence_multiplier=1.96):
        """Return componentwise intervals using an explicitly declared multiplier."""
        if not np.isfinite(confidence_multiplier) or confidence_multiplier <= 0.0:
            raise ValueError("confidence_multiplier must be finite and positive")
        halfwidth = float(confidence_multiplier) * self.standard_error_displacement_xy_m
        return (
            self.mean_displacement_xy_m - halfwidth,
            self.mean_displacement_xy_m + halfwidth,
        )


def measure_feature_centerline_3d(
        geometry: FeatureGeometry3D, *, lateral_bounds_m, opening_width_m,
        onset_displacement_m, minimum_slice_cells=3, reference_slice_count=2):
    """Measure the connected-opening proxy inside a declared lateral region of interest.

    The caller supplies the mask-opening bounds; this prevents exterior vacuum from being mistaken
    for a hole. For each z plane, the centroid of negative-level-set nodes inside that region is
    measured. Depth is referenced to the highest valid plane and displacement to the mean centroid
    of the highest ``reference_slice_count`` valid planes.
    """
    if not isinstance(geometry, FeatureGeometry3D):
        raise TypeError("geometry must be FeatureGeometry3D")
    bounds = np.asarray(lateral_bounds_m, dtype=float)
    if (bounds.shape != (4,) or np.any(~np.isfinite(bounds))
            or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]
            or not np.isfinite(opening_width_m) or opening_width_m <= 0.0
            or not np.isfinite(onset_displacement_m) or onset_displacement_m <= 0.0
            or int(minimum_slice_cells) != minimum_slice_cells or minimum_slice_cells < 1
            or int(reference_slice_count) != reference_slice_count or reference_slice_count < 1):
        raise ValueError("invalid centerline measurement contract")
    unit = geometry.mesh_length_unit_m
    origin = np.asarray(geometry.mesh_origin_m, dtype=float)
    coordinates = tuple(
        origin[axis] + np.arange(size) * geometry.dx * unit
        for axis, size in enumerate(geometry.phi.shape))
    x, y, z = coordinates
    lateral = ((x[:, None] >= bounds[0]) & (x[:, None] <= bounds[1])
               & (y[None, :] >= bounds[2]) & (y[None, :] <= bounds[3]))
    order = np.argsort(z)[::-1]
    centroid = np.full((len(z), 2), np.nan)
    diameter = np.full(len(z), np.nan)
    valid = np.zeros(len(z), dtype=bool)
    cell_area = (geometry.dx * unit) ** 2
    for output_index, k in enumerate(order):
        gas = (geometry.phi[:, :, k] < 0.0) & lateral
        selected = np.argwhere(gas)
        if len(selected) < int(minimum_slice_cells):
            continue
        centroid[output_index] = [
            float(np.mean(x[selected[:, 0]])),
            float(np.mean(y[selected[:, 1]]))]
        diameter[output_index] = 2.0 * np.sqrt(len(selected) * cell_area / np.pi)
        valid[output_index] = True
    if not np.any(valid):
        raise ValueError("declared lateral region contains no resolved open feature")
    valid_index = np.flatnonzero(valid)
    reference_index = valid_index[:min(int(reference_slice_count), len(valid_index))]
    reference = np.mean(centroid[reference_index], axis=0)
    displacement = centroid - reference
    top_z = float(z[order[valid_index[0]]])
    depth = top_z - z[order]
    return FeatureCenterline3D(
        depth, centroid, displacement, diameter, valid,
        float(opening_width_m), float(onset_displacement_m))


def measure_feature_centerline_ensemble_3d(
        geometries, *, lateral_bounds_m, opening_width_m,
        onset_displacement_m, minimum_slice_cells=3, reference_slice_count=2):
    """Measure common-grid centerline statistics for independent finite-arrival geometries."""
    geometries = tuple(geometries)
    if len(geometries) < 2:
        raise ValueError("a centerline ensemble requires at least two geometries")
    members = tuple(measure_feature_centerline_3d(
        geometry, lateral_bounds_m=lateral_bounds_m,
        opening_width_m=opening_width_m,
        onset_displacement_m=onset_displacement_m,
        minimum_slice_cells=minimum_slice_cells,
        reference_slice_count=reference_slice_count) for geometry in geometries)
    reference_depth = members[0].depth_m
    if any(not np.array_equal(item.depth_m, reference_depth) for item in members[1:]):
        raise ValueError("ensemble geometries must share one physical depth grid")
    stack = np.stack([
        np.where(item.valid[:, None], item.displacement_xy_m, np.nan)
        for item in members])
    count = np.sum(np.isfinite(stack[:, :, 0]), axis=0)
    mean = np.full(stack.shape[1:], np.nan)
    np.divide(
        np.nansum(stack, axis=0), count[:, None], out=mean,
        where=count[:, None] > 0)
    centered = stack - mean[None, :, :]
    deviation = np.full_like(mean, np.nan)
    np.sqrt(np.divide(
        np.nansum(centered ** 2, axis=0), (count - 1)[:, None],
        out=np.zeros_like(mean), where=(count - 1)[:, None] > 0), out=deviation,
        where=(count - 1)[:, None] > 0)
    radial_standard_error = np.sqrt(np.divide(
        np.nansum(deviation ** 2, axis=1), count,
        out=np.full(len(count), np.nan), where=count > 0))
    mean_magnitude = np.linalg.norm(mean, axis=1)
    zscore = np.divide(
        mean_magnitude, radial_standard_error,
        out=np.zeros_like(mean_magnitude), where=radial_standard_error > 0.0)
    onset = np.asarray([
        np.nan if item.onset_aspect_ratio is None else item.onset_aspect_ratio
        for item in members])
    return FeatureCenterlineEnsemble3D(
        members, mean, deviation, count,
        float(np.nanmax(zscore, initial=0.0)), onset)


def _linear_zero_crossing(coordinate_a, value_a, coordinate_b, value_b):
    denominator = value_b - value_a
    if denominator == 0.0:
        return 0.5 * (coordinate_a + coordinate_b)
    fraction = float(np.clip(-value_a / denominator, 0.0, 1.0))
    return float(coordinate_a + fraction * (coordinate_b - coordinate_a))


def measure_trench_profile_observables_3d(
        geometry: FeatureGeometry3D, *, lateral_bounds_m, opening_center_x_m,
        feature_top_z_m, reference_depth_interval_m, bow_depth_interval_m,
        notch_depth_interval_m, minimum_longitudinal_rows=3):
    """Measure one bounded line opening without assuming a benchmark-specific sidewall label.

    At each depth and longitudinal row, this follows the connected gas interval containing the
    declared opening center and linearly interpolates its two hard level-set crossings.  Rows whose
    interval touches the declared ROI are rejected rather than silently treating an exterior-vacuum
    boundary as a material sidewall.
    """
    if not isinstance(geometry, FeatureGeometry3D):
        raise TypeError("geometry must be FeatureGeometry3D")
    bounds = np.asarray(lateral_bounds_m, dtype=float)
    intervals = tuple(np.asarray(item, dtype=float) for item in (
        reference_depth_interval_m, bow_depth_interval_m, notch_depth_interval_m))
    if (bounds.shape != (4,) or np.any(~np.isfinite(bounds))
            or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]
            or not np.isfinite(opening_center_x_m)
            or not bounds[0] < opening_center_x_m < bounds[1]
            or not np.isfinite(feature_top_z_m)
            or int(minimum_longitudinal_rows) != minimum_longitudinal_rows
            or minimum_longitudinal_rows < 1
            or any(item.shape != (2,) or np.any(~np.isfinite(item))
                   or item[0] < 0.0 or item[1] < item[0] for item in intervals)):
        raise ValueError("invalid trench-profile measurement contract")
    unit = geometry.mesh_length_unit_m
    origin = np.asarray(geometry.mesh_origin_m, dtype=float)
    x, y, z = tuple(
        origin[axis] + np.arange(size) * geometry.dx * unit
        for axis, size in enumerate(geometry.phi.shape))
    x_index = np.flatnonzero((x >= bounds[0]) & (x <= bounds[1]))
    y_index = np.flatnonzero((y >= bounds[2]) & (y <= bounds[3]))
    if len(x_index) < 3 or len(y_index) < int(minimum_longitudinal_rows):
        raise ValueError("declared trench-profile ROI is not spatially resolved")
    center_index = int(x_index[np.argmin(np.abs(x[x_index] - opening_center_x_m))])
    z_index = np.flatnonzero(z <= feature_top_z_m + 0.5 * geometry.dx * unit)
    if not len(z_index):
        raise ValueError("feature_top_z_m lies below the geometry")
    z_index = z_index[np.argsort(z[z_index])[::-1]]
    depth = np.maximum(0.0, float(feature_top_z_m) - z[z_index])
    left = np.full(len(z_index), np.nan)
    right = np.full(len(z_index), np.nan)
    valid = np.zeros(len(z_index), dtype=bool)
    first_x = int(x_index[0])
    last_x = int(x_index[-1])
    for output_index, k in enumerate(z_index):
        row_left = []
        row_right = []
        for j in y_index:
            values = geometry.phi[:, j, k]
            if values[center_index] >= 0.0:
                continue
            inner_left = center_index
            while inner_left > first_x and values[inner_left - 1] < 0.0:
                inner_left -= 1
            inner_right = center_index
            while inner_right < last_x and values[inner_right + 1] < 0.0:
                inner_right += 1
            if (inner_left <= first_x or inner_right >= last_x
                    or values[inner_left - 1] < 0.0
                    or values[inner_right + 1] < 0.0):
                continue
            row_left.append(_linear_zero_crossing(
                x[inner_left - 1], values[inner_left - 1],
                x[inner_left], values[inner_left]))
            row_right.append(_linear_zero_crossing(
                x[inner_right], values[inner_right],
                x[inner_right + 1], values[inner_right + 1]))
        if len(row_left) < int(minimum_longitudinal_rows):
            continue
        left[output_index] = float(np.mean(row_left))
        right[output_index] = float(np.mean(row_right))
        valid[output_index] = True
    reference = valid & (depth >= intervals[0][0]) & (depth <= intervals[0][1])
    if not np.any(reference):
        raise ValueError("reference sidewall interval has no resolved bounded opening")
    return TrenchProfileObservables3D(
        depth, left, right, valid,
        tuple(float(value) for value in intervals[0]),
        tuple(float(value) for value in intervals[1]),
        tuple(float(value) for value in intervals[2]),
        float(np.mean(left[reference])), float(np.mean(right[reference])))


def _ensemble_scalar_estimate(values, confidence_multiplier):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError("ensemble observable values must be finite")
    deviation = float(np.std(values, ddof=1))
    standard_error = deviation / np.sqrt(len(values))
    return EnsembleScalarEstimate3D(
        float(np.mean(values)), deviation, standard_error,
        float(confidence_multiplier) * standard_error,
        len(values), float(confidence_multiplier))


def measure_trench_profile_ensemble_3d(
        geometries, *, confidence_multiplier=1.96, **measurement_contract):
    """Measure notch/bow ensemble intervals from independent realized geometries."""
    geometries = tuple(geometries)
    if len(geometries) < 2:
        raise ValueError("a trench-profile ensemble requires at least two geometries")
    if not np.isfinite(confidence_multiplier) or confidence_multiplier <= 0.0:
        raise ValueError("confidence_multiplier must be finite and positive")
    members = tuple(measure_trench_profile_observables_3d(
        geometry, **measurement_contract) for geometry in geometries)
    fields = (
        "maximum_notch_depth_m", "maximum_left_notch_depth_m",
        "maximum_right_notch_depth_m", "notch_asymmetry_m",
        "maximum_bow_width_m", "maximum_bow_expansion_m",
    )
    estimates = tuple(_ensemble_scalar_estimate(
        [getattr(member, field) for member in members], confidence_multiplier)
        for field in fields)
    return TrenchProfileEnsemble3D(members, *estimates)
