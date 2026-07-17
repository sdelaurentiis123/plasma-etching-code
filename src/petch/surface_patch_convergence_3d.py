"""Fixed-physical-patch convergence instruments for triangle surface fields.

The triangle operator remains authoritative.  This module only integrates its outputs over
Cartesian surface patches whose dimensions and origin are declared in SI units.  Consequently a
mesh refinement, marching-cubes retriangulation, or future AMR hierarchy cannot silently redefine
the spatial scale at which a numerical or experimental claim is scored.

Patches are split by material and dominant oriented gas-normal class so a wall and a floor that
share one Cartesian box never become one balance equation.  The stable integer patch keys permit
comparisons between independently triangulated representations of the same physical surface.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256

import numpy as np
from scipy.stats import t as student_t


# A patch mean is a claim about a spatially resolved region, not an arbitrarily small triangle
# sliver left by Cartesian clipping.  The default is deliberately a geometry threshold rather
# than a physics-error tolerance: every patch remains in the conservative integrated gate, while
# means are scored only where at least ten percent of the represented tangential footprint exists.
# This is a geometry-independent round fraction, not a value tuned to a benchmark result.  Formal
# campaigns must pass it explicitly and report sensitivity; the default exists for API convenience.
# Integrated inventory on every patch remains authoritative regardless of this mean-support rule.
DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION = 0.1


def _readonly(value, dtype):
    array = np.asarray(value, dtype=dtype).copy()
    array.setflags(write=False)
    return array


def _digest_array(digest, name, value):
    array = np.ascontiguousarray(np.asarray(value))
    digest.update(name.encode("utf-8") + b"\0")
    digest.update(array.dtype.str.encode("ascii") + b"\0")
    digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    digest.update(array.tobytes(order="C"))


def _validated_surface_inputs(
        face_centroids, face_gas_normals, face_material_id, *, mesh_length_unit_m,
        mesh_origin_m, patch_origin_m):
    centroid = np.asarray(face_centroids, dtype=float)
    normal = np.asarray(face_gas_normals, dtype=float)
    material = np.asarray(face_material_id, dtype=int)
    origin = np.asarray(mesh_origin_m, dtype=float)
    patch_origin = origin if patch_origin_m is None else np.asarray(patch_origin_m, dtype=float)
    if (centroid.ndim != 2 or centroid.shape[1] != 3
            or normal.shape != centroid.shape or material.shape != (len(centroid),)
            or len(centroid) == 0 or np.any(~np.isfinite(centroid))
            or np.any(~np.isfinite(normal)) or np.any(material <= 0)
            or not np.allclose(np.linalg.norm(normal, axis=1), 1.0, rtol=0.0, atol=2e-6)
            or origin.shape != (3,) or patch_origin.shape != (3,)
            or np.any(~np.isfinite(origin)) or np.any(~np.isfinite(patch_origin))
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("invalid physical surface-patch inputs")
    return centroid, normal, material, origin, patch_origin


def physical_surface_patch_keys_3d(
        face_centroids, face_gas_normals, face_material_id, patch_scale_m, *,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None):
    """Return stable patch keys and one group index per face.

    Each key is ``(material, normal_axis, normal_sign, ix, iy, iz)``.  Indices refer to fixed
    physical Cartesian boxes anchored at ``patch_origin_m``; they do not refer to mesh cells.
    """
    centroid, normal, material, origin, patch_origin = _validated_surface_inputs(
        face_centroids, face_gas_normals, face_material_id,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        patch_origin_m=patch_origin_m)
    if not np.isfinite(patch_scale_m) or patch_scale_m <= 0.0:
        raise ValueError("patch_scale_m must be finite and positive")
    physical = origin + centroid * float(mesh_length_unit_m)
    cell = np.floor(
        (physical - patch_origin) / float(patch_scale_m) + 1e-12).astype(np.int64)
    dominant_axis = np.argmax(np.abs(normal), axis=1)
    dominant_sign = (normal[np.arange(len(normal)), dominant_axis] >= 0.0).astype(int)
    face_key = np.column_stack((material, dominant_axis, dominant_sign, cell))
    patch_key, group = np.unique(face_key, axis=0, return_inverse=True)
    return _readonly(patch_key, np.int64), _readonly(group, np.int64)


def physical_surface_patch_groups_3d(
        face_centroids, face_gas_normals, face_material_id, patch_scale_m, *,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None):
    """Compatibility wrapper returning only fixed physical group indices."""
    return physical_surface_patch_keys_3d(
        face_centroids, face_gas_normals, face_material_id, patch_scale_m,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        patch_origin_m=patch_origin_m)[1]


@dataclass(frozen=True)
class PhysicalPatchField3D:
    """Area-conservative aggregation of one scalar triangle field."""

    patch_scale_m: float
    patch_origin_m: np.ndarray
    patch_key: np.ndarray
    contribution_face_index: np.ndarray
    contribution_patch_index: np.ndarray
    contribution_area_m2: np.ndarray
    patch_area_m2: np.ndarray
    patch_projected_support_area_m2: np.ndarray
    patch_nominal_projected_area_m2: np.ndarray
    patch_projected_support_fraction: np.ndarray
    integrated_field_area: np.ndarray
    mean_field: np.ndarray
    maximum_absolute_face_value: float
    periodic_domain_origin_m: np.ndarray
    periodic_domain_lengths_m: np.ndarray
    scheme_sha256: str

    def __post_init__(self):
        key = _readonly(self.patch_key, np.int64)
        face_index = _readonly(self.contribution_face_index, np.int64)
        patch_index = _readonly(self.contribution_patch_index, np.int64)
        contribution_area = _readonly(self.contribution_area_m2, float)
        origin = _readonly(self.patch_origin_m, float)
        area = _readonly(self.patch_area_m2, float)
        projected = _readonly(self.patch_projected_support_area_m2, float)
        nominal = _readonly(self.patch_nominal_projected_area_m2, float)
        support_fraction = _readonly(self.patch_projected_support_fraction, float)
        integrated = _readonly(self.integrated_field_area, float)
        mean = _readonly(self.mean_field, float)
        periodic_origin = _readonly(self.periodic_domain_origin_m, float)
        periodic_lengths = _readonly(self.periodic_domain_lengths_m, float)
        count = len(key)
        values = np.asarray([
            self.patch_scale_m, self.maximum_absolute_face_value], dtype=float)
        if (key.ndim != 2 or key.shape[1] != 6 or face_index.ndim != 1
                or patch_index.shape != face_index.shape
                or contribution_area.shape != face_index.shape
                or origin.shape != (3,) or area.shape != (count,)
                or projected.shape != (count,) or nominal.shape != (count,)
                or support_fraction.shape != (count,)
                or integrated.shape != (count,) or mean.shape != (count,)
                or periodic_origin.shape != (3,) or periodic_lengths.shape != (3,)
                or np.any(face_index < 0) or np.any(patch_index < 0)
                or np.any(patch_index >= count)
                or np.any(~np.isfinite(contribution_area))
                or np.any(contribution_area <= 0.0)
                or np.any(~np.isfinite(area)) or np.any(area <= 0.0)
                or np.any(~np.isfinite(projected)) or np.any(projected <= 0.0)
                or np.any(~np.isfinite(nominal)) or np.any(nominal <= 0.0)
                or np.any(~np.isfinite(support_fraction))
                or np.any(support_fraction <= 0.0)
                or np.any(~np.isfinite(integrated)) or np.any(~np.isfinite(mean))
                or np.any(~np.isfinite(periodic_origin))
                or np.any(~np.isfinite(periodic_lengths))
                or np.any(periodic_lengths < 0.0)
                or not np.allclose(
                    support_fraction, projected / nominal, rtol=2e-14, atol=0.0)
                or np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or not isinstance(self.scheme_sha256, str) or len(self.scheme_sha256) != 64):
            raise ValueError("invalid physical patch-field receipt")
        object.__setattr__(self, "patch_key", key)
        object.__setattr__(self, "contribution_face_index", face_index)
        object.__setattr__(self, "contribution_patch_index", patch_index)
        object.__setattr__(self, "contribution_area_m2", contribution_area)
        object.__setattr__(self, "patch_origin_m", origin)
        object.__setattr__(self, "patch_area_m2", area)
        object.__setattr__(self, "patch_projected_support_area_m2", projected)
        object.__setattr__(self, "patch_nominal_projected_area_m2", nominal)
        object.__setattr__(self, "patch_projected_support_fraction", support_fraction)
        object.__setattr__(self, "integrated_field_area", integrated)
        object.__setattr__(self, "mean_field", mean)
        object.__setattr__(self, "periodic_domain_origin_m", periodic_origin)
        object.__setattr__(self, "periodic_domain_lengths_m", periodic_lengths)


@dataclass(frozen=True)
class ReplicatedPhysicalPatchScore3D:
    """Independent-replicate uncertainty for one field on one physical patch scale.

    The authority field is normally produced by the mean transport operator.  Because a nonlinear
    downstream solve need not commute with averaging, the score keeps two errors separate: the
    Student-t interval of the replicate mean and the authority-to-replicate-mean discrepancy.  Their
    sum is the conservative envelope used by the prospective stopping gate.
    """

    authority: PhysicalPatchField3D
    replicate_count: int
    confidence_level: float
    student_t_critical: float
    absolute_tolerance: float
    relative_tolerance: float
    minimum_mean_support_fraction: float
    mean_eligible_patch_mask: np.ndarray
    eligible_mean_patch_count: int
    excluded_mean_patch_count: int
    excluded_mean_surface_area_m2: float
    excluded_mean_surface_area_fraction: float
    excluded_mean_projected_support_area_m2: float
    excluded_mean_projected_support_fraction: float
    replicate_mean_integrated_field_area: np.ndarray
    replicate_mean_field: np.ndarray
    integrated_confidence_half_width: np.ndarray
    mean_confidence_half_width: np.ndarray
    maximum_integrated_confidence_mixed_normalized: float
    maximum_mean_confidence_mixed_normalized: float
    maximum_integrated_authority_bias_mixed_normalized: float
    maximum_mean_authority_bias_mixed_normalized: float
    maximum_integrated_combined_mixed_normalized: float
    maximum_mean_combined_mixed_normalized: float
    maximum_mean_combined_all_patches_mixed_normalized: float

    def __post_init__(self):
        if not isinstance(self.authority, PhysicalPatchField3D):
            raise TypeError("authority must be a PhysicalPatchField3D")
        patch_count = len(self.authority.patch_key)
        arrays = {}
        for name in (
                "replicate_mean_integrated_field_area", "replicate_mean_field",
                "integrated_confidence_half_width", "mean_confidence_half_width"):
            value = _readonly(getattr(self, name), float)
            if (value.shape != (patch_count,) or np.any(~np.isfinite(value))
                    or ("half_width" in name and np.any(value < 0.0))):
                raise ValueError("invalid replicated physical-patch array")
            arrays[name] = value
        eligible = _readonly(self.mean_eligible_patch_mask, bool)
        if eligible.shape != (patch_count,) or not np.any(eligible):
            raise ValueError("at least one physical patch must support a mean gate")
        scalars = np.asarray([
            self.confidence_level, self.student_t_critical,
            self.absolute_tolerance, self.relative_tolerance,
            self.minimum_mean_support_fraction,
            self.excluded_mean_surface_area_m2,
            self.excluded_mean_surface_area_fraction,
            self.excluded_mean_projected_support_area_m2,
            self.excluded_mean_projected_support_fraction,
            self.maximum_integrated_confidence_mixed_normalized,
            self.maximum_mean_confidence_mixed_normalized,
            self.maximum_integrated_authority_bias_mixed_normalized,
            self.maximum_mean_authority_bias_mixed_normalized,
            self.maximum_integrated_combined_mixed_normalized,
            self.maximum_mean_combined_mixed_normalized,
            self.maximum_mean_combined_all_patches_mixed_normalized,
        ], dtype=float)
        if (int(self.replicate_count) != self.replicate_count
                or self.replicate_count < 4 or np.any(~np.isfinite(scalars))
                or not 0.0 < self.confidence_level < 1.0
                or self.student_t_critical <= 0.0
                or self.absolute_tolerance <= 0.0 or self.relative_tolerance < 0.0
                or not 0.0 < self.minimum_mean_support_fraction <= 1.0
                or int(self.eligible_mean_patch_count) != self.eligible_mean_patch_count
                or int(self.excluded_mean_patch_count) != self.excluded_mean_patch_count
                or self.eligible_mean_patch_count != int(np.count_nonzero(eligible))
                or self.excluded_mean_patch_count != int(np.count_nonzero(~eligible))
                or self.eligible_mean_patch_count + self.excluded_mean_patch_count != patch_count
                or np.any(scalars[5:] < 0.0)
                or self.excluded_mean_surface_area_fraction > 1.0 + 1e-14
                or self.excluded_mean_projected_support_fraction > 1.0 + 1e-14):
            raise ValueError("invalid replicated physical-patch score")
        object.__setattr__(self, "replicate_count", int(self.replicate_count))
        object.__setattr__(self, "eligible_mean_patch_count", int(
            self.eligible_mean_patch_count))
        object.__setattr__(self, "excluded_mean_patch_count", int(
            self.excluded_mean_patch_count))
        object.__setattr__(self, "mean_eligible_patch_mask", eligible)
        for name, value in arrays.items():
            object.__setattr__(self, name, value)

    @property
    def all_mixed_tolerances_pass(self):
        return bool(
            self.maximum_integrated_combined_mixed_normalized <= 1.0
            and self.maximum_mean_combined_mixed_normalized <= 1.0)


def _clip_polygon_axis_3d(polygon, axis, boundary, keep_greater):
    """Clip one planar polygon by one axis-aligned half-space."""
    if len(polygon) == 0:
        return polygon
    output = []
    previous = polygon[-1]
    previous_signed = previous[axis] - boundary
    previous_inside = previous_signed >= 0.0 if keep_greater else previous_signed <= 0.0
    for current in polygon:
        current_signed = current[axis] - boundary
        current_inside = current_signed >= 0.0 if keep_greater else current_signed <= 0.0
        if current_inside != previous_inside:
            denominator = current[axis] - previous[axis]
            if denominator != 0.0:
                fraction = (boundary - previous[axis]) / denominator
                output.append(previous + fraction * (current - previous))
        if current_inside:
            output.append(current)
        previous = current
        previous_inside = current_inside
    return np.asarray(output, dtype=float).reshape(-1, 3)


def _polygon_area_3d(polygon):
    if len(polygon) < 3:
        return 0.0
    anchor = polygon[0]
    return 0.5 * sum(
        np.linalg.norm(np.cross(polygon[index] - anchor, polygon[index + 1] - anchor))
        for index in range(1, len(polygon) - 1))


def _triangle_patch_overlap_contributions_3d(
        physical_triangle, patch_scale_m, patch_origin_m, *, geometry_epsilon=None):
    """Return exact constant-P0 triangle area in every intersected Cartesian patch."""
    scale = float(patch_scale_m)
    epsilon = (
        np.finfo(float).eps if geometry_epsilon is None
        else float(geometry_epsilon))
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("geometry epsilon must be finite and positive")
    # Arithmetic in mesh-unit -> metre conversion can place a coplanar triangle one ULP to either
    # side of an exact physical patch plane.  Snap only a whole, numerically coplanar coordinate;
    # snapping individual vertices can collapse legitimate thin marching-surface triangles.
    physical_triangle = np.asarray(physical_triangle, dtype=float).copy()
    normalized = (physical_triangle - patch_origin_m) / scale
    nearest_plane = np.rint(normalized)
    axis_scale = np.maximum(1.0, np.max(np.abs(normalized), axis=0))
    plane_tolerance = 2.0 * epsilon * axis_scale
    for axis in range(3):
        candidate = nearest_plane[0, axis]
        if (np.all(nearest_plane[:, axis] == candidate)
                and np.all(np.abs(normalized[:, axis] - candidate)
                           <= plane_tolerance[axis])):
            normalized[:, axis] = candidate
            physical_triangle[:, axis] = patch_origin_m[axis] + scale * candidate
    # Enumerate both cells touched by a non-coplanar boundary vertex.  Degenerate roundoff pieces
    # are removed after exact clipping, rather than using a coordinate tolerance that can erase a
    # real fraction of an extremely thin triangle.
    lower_index = np.floor(np.min(normalized, axis=0)).astype(int)
    upper_index = np.floor(np.max(normalized, axis=0)).astype(int)
    upper_index = np.maximum(upper_index, lower_index)
    output = []
    for ix in range(lower_index[0], upper_index[0] + 1):
        for iy in range(lower_index[1], upper_index[1] + 1):
            for iz in range(lower_index[2], upper_index[2] + 1):
                cell = np.asarray([ix, iy, iz], dtype=int)
                lower = patch_origin_m + scale * cell
                upper = lower + scale
                polygon = np.asarray(physical_triangle, dtype=float)
                for axis in range(3):
                    polygon = _clip_polygon_axis_3d(
                        polygon, axis, lower[axis], True)
                    polygon = _clip_polygon_axis_3d(
                        polygon, axis, upper[axis], False)
                    if len(polygon) < 3:
                        break
                area = _polygon_area_3d(polygon)
                if area > 0.0:
                    output.append((tuple(int(value) for value in cell), area))
    return output


def _validated_periodic_domain_3d(
        periodic_domain_lengths_m, periodic_domain_origin_m, mesh_origin_m):
    """Return an explicit fundamental-domain contract.

    Zero length marks a nonperiodic coordinate.  A positive length marks the single represented
    periodic cell; physical-patch support is normalized by its actual overlap with that cell rather
    than by an unavailable repeated width.  This matters when, for example, a 40 nm scoring patch
    spans a geometry represented by one 20 nm periodic cell.
    """
    lengths = (
        np.zeros(3, dtype=float) if periodic_domain_lengths_m is None
        else np.asarray(periodic_domain_lengths_m, dtype=float))
    origin = (
        np.asarray(mesh_origin_m, dtype=float) if periodic_domain_origin_m is None
        else np.asarray(periodic_domain_origin_m, dtype=float))
    if (lengths.shape != (3,) or origin.shape != (3,)
            or np.any(~np.isfinite(lengths)) or np.any(lengths < 0.0)
            or np.any(~np.isfinite(origin))):
        raise ValueError("periodic domain requires finite nonnegative lengths and an origin")
    return origin, lengths


def _nominal_projected_patch_areas_3d(
        patch_key, patch_scale_m, patch_origin_m, periodic_origin_m,
        periodic_lengths_m):
    """Nominal represented tangential footprint for every oriented surface patch."""
    scale = float(patch_scale_m)
    output = np.empty(len(patch_key), dtype=float)
    for patch_index, key in enumerate(np.asarray(patch_key, dtype=np.int64)):
        normal_axis = int(key[1])
        cell = key[3:]
        footprint = 1.0
        for axis in range(3):
            if axis == normal_axis:
                continue
            periodic_length = float(periodic_lengths_m[axis])
            if periodic_length == 0.0:
                represented_extent = scale
            else:
                patch_lower = float(patch_origin_m[axis] + scale * cell[axis])
                patch_upper = patch_lower + scale
                domain_lower = float(periodic_origin_m[axis])
                domain_upper = domain_lower + periodic_length
                represented_extent = max(
                    0.0, min(patch_upper, domain_upper) - max(patch_lower, domain_lower))
                # A surface contribution outside the declared fundamental domain means the domain
                # provenance is wrong.  Refuse instead of manufacturing a denominator.
                if represented_extent <= 0.0:
                    raise ValueError(
                        "surface patch lies outside the declared periodic fundamental domain")
            footprint *= represented_extent
        output[patch_index] = footprint
    return output


def aggregate_surface_field_on_physical_patches_3d(
        face_field, face_area_m2, verts, faces, face_gas_normals,
        face_material_id, patch_scale_m, *, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None,
        periodic_domain_lengths_m=None, periodic_domain_origin_m=None):
    """Conservatively integrate a P0 face field over fixed physical patches.

    Triangles are clipped against every intersected physical patch.  A triangle crossing a patch
    boundary therefore contributes to both patches in proportion to actual overlap area; assigning
    the whole triangle by its centroid would not be invariant under retriangulation.
    """
    field = np.asarray(face_field, dtype=float)
    area = np.asarray(face_area_m2, dtype=float)
    supplied_vertices = np.asarray(verts)
    geometry_epsilon = (
        np.finfo(supplied_vertices.dtype).eps
        if np.issubdtype(supplied_vertices.dtype, np.floating)
        else np.finfo(float).eps)
    vertices = np.asarray(supplied_vertices, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    normals = np.asarray(face_gas_normals, dtype=float)
    material = np.asarray(face_material_id, dtype=int)
    if (field.ndim != 1 or area.shape != field.shape or len(field) == 0
            or np.any(~np.isfinite(field)) or np.any(~np.isfinite(area))
            or np.any(area <= 0.0) or vertices.ndim != 2 or vertices.shape[1] != 3
            or triangles.shape != (len(field), 3) or np.any(triangles < 0)
            or np.any(triangles >= len(vertices))
            or normals.shape != (len(field), 3) or material.shape != (len(field),)
            or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(normals))
            or np.any(material <= 0)
            or not np.allclose(
                np.linalg.norm(normals, axis=1), 1.0, rtol=0.0, atol=2e-6)):
        raise ValueError("invalid surface field, area, or triangle geometry")
    if (not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0
            or not np.isfinite(patch_scale_m) or patch_scale_m <= 0.0):
        raise ValueError("physical length scales must be finite and positive")
    mesh_origin = np.asarray(mesh_origin_m, dtype=float)
    patch_origin = (
        mesh_origin if patch_origin_m is None
        else np.asarray(patch_origin_m, dtype=float))
    if (mesh_origin.shape != (3,) or patch_origin.shape != (3,)
            or np.any(~np.isfinite(mesh_origin)) or np.any(~np.isfinite(patch_origin))):
        raise ValueError("mesh and patch origins must be finite three-vectors")
    periodic_origin, periodic_lengths = _validated_periodic_domain_3d(
        periodic_domain_lengths_m, periodic_domain_origin_m, mesh_origin)
    physical_vertices = mesh_origin + vertices * float(mesh_length_unit_m)
    geometric_area = 0.5 * np.linalg.norm(np.cross(
        physical_vertices[triangles[:, 1]] - physical_vertices[triangles[:, 0]],
        physical_vertices[triangles[:, 2]] - physical_vertices[triangles[:, 0]]), axis=1)
    if not np.allclose(area, geometric_area, rtol=2e-10, atol=1e-30):
        raise ValueError("declared face areas disagree with physical triangle geometry")

    contribution = []
    for face_index, triangle in enumerate(physical_vertices[triangles]):
        dominant_axis = int(np.argmax(np.abs(normals[face_index])))
        dominant_sign = int(normals[face_index, dominant_axis] >= 0.0)
        raw_local = _triangle_patch_overlap_contributions_3d(
            triangle, patch_scale_m, patch_origin,
            geometry_epsilon=geometry_epsilon)
        # Inclusive plane clipping can leave a machine-epsilon polygon on the opposite side of an
        # exact patch boundary.  Such a piece is neither resolved area nor a physical patch; if kept,
        # a max norm over patches becomes a disguised face/roundoff statistic.  Remove only pieces
        # below a triangle-relative roundoff bound, then restore the exact authoritative face area
        # by conservative renormalization of the retained pieces.
        overlap_floor = 64.0 * geometry_epsilon * area[face_index]
        local = [item for item in raw_local if item[1] > overlap_floor]
        if not local:
            raise RuntimeError("triangle-to-patch overlap retained no resolved area")
        local_total = sum(value for _, value in local)
        raw_total = sum(value for _, value in raw_local)
        # Plane snapping can change the clipped geometric area by at most the precision of the
        # supplied mesh.  Accept only that bounded change; the scale below still returns the exact
        # authoritative face area and therefore closes inventory exactly.
        area_rtol = max(2e-10, 128.0 * geometry_epsilon)
        if not np.isclose(
                raw_total, area[face_index], rtol=area_rtol, atol=1e-30):
            raise RuntimeError("triangle-to-patch overlap failed area conservation")
        closure_scale = area[face_index] / local_total
        for cell, overlap_area in local:
            key = (
                int(material[face_index]), dominant_axis, dominant_sign,
                int(cell[0]), int(cell[1]), int(cell[2]))
            contribution.append((key, face_index, overlap_area * closure_scale))
    patch_keys = sorted({item[0] for item in contribution})
    key_to_patch = {key: index for index, key in enumerate(patch_keys)}
    face_index = np.asarray([item[1] for item in contribution], dtype=int)
    patch_index = np.asarray([key_to_patch[item[0]] for item in contribution], dtype=int)
    contribution_area = np.asarray([item[2] for item in contribution], dtype=float)
    key = np.asarray(patch_keys, dtype=np.int64)
    patch_area = np.bincount(
        patch_index, weights=contribution_area, minlength=len(key))
    projected_support_area = np.bincount(
        patch_index,
        weights=(contribution_area * np.abs(
            normals[face_index, key[patch_index, 1]])),
        minlength=len(key))
    nominal_projected_area = _nominal_projected_patch_areas_3d(
        key, patch_scale_m, patch_origin, periodic_origin, periodic_lengths)
    projected_support_fraction = projected_support_area / nominal_projected_area
    integrated = np.bincount(
        patch_index, weights=contribution_area * field[face_index], minlength=len(key))
    mean = integrated / patch_area
    digest = sha256(b"petch.physical-surface-patch-scheme-3d.v3\0")
    digest.update(np.float64(patch_scale_m).tobytes())
    digest.update(np.float64(mesh_length_unit_m).tobytes())
    digest.update(np.float64(geometry_epsilon).tobytes())
    _digest_array(digest, "patch_origin_m", patch_origin)
    _digest_array(digest, "periodic_domain_origin_m", periodic_origin)
    _digest_array(digest, "periodic_domain_lengths_m", periodic_lengths)
    _digest_array(digest, "patch_key", key)
    _digest_array(digest, "verts", vertices)
    _digest_array(digest, "faces", triangles)
    _digest_array(digest, "contribution_face_index", face_index)
    _digest_array(digest, "contribution_patch_index", patch_index)
    _digest_array(digest, "contribution_area_m2", contribution_area)
    _digest_array(digest, "face_area_m2", area)
    return PhysicalPatchField3D(
        patch_scale_m=float(patch_scale_m), patch_origin_m=patch_origin,
        patch_key=key, contribution_face_index=face_index,
        contribution_patch_index=patch_index,
        contribution_area_m2=contribution_area, patch_area_m2=patch_area,
        patch_projected_support_area_m2=projected_support_area,
        patch_nominal_projected_area_m2=nominal_projected_area,
        patch_projected_support_fraction=projected_support_fraction,
        integrated_field_area=integrated, mean_field=mean,
        maximum_absolute_face_value=float(np.max(np.abs(field), initial=0.0)),
        periodic_domain_origin_m=periodic_origin,
        periodic_domain_lengths_m=periodic_lengths,
        scheme_sha256=digest.hexdigest())


def _mean_support_partition_3d(authority, minimum_support_fraction):
    threshold = float(minimum_support_fraction)
    if not np.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise ValueError("minimum mean support fraction must be in (0, 1]")
    eligible = authority.patch_projected_support_fraction >= threshold
    if not np.any(eligible):
        raise ValueError("no physical patch meets the declared mean-support fraction")
    excluded = ~eligible
    excluded_surface_area = float(np.sum(authority.patch_area_m2[excluded]))
    total_surface_area = float(np.sum(authority.patch_area_m2))
    excluded_projected_area = float(np.sum(
        authority.patch_projected_support_area_m2[excluded]))
    total_projected_area = float(np.sum(authority.patch_projected_support_area_m2))
    return {
        "threshold": threshold,
        "eligible": eligible,
        "eligible_count": int(np.count_nonzero(eligible)),
        "excluded_count": int(np.count_nonzero(excluded)),
        "excluded_surface_area_m2": excluded_surface_area,
        "excluded_surface_area_fraction": excluded_surface_area / total_surface_area,
        "excluded_projected_support_area_m2": excluded_projected_area,
        "excluded_projected_support_fraction": excluded_projected_area / total_projected_area,
    }


def score_replicated_surface_field_on_physical_patches_3d(
        authority_face_field, replicate_face_fields, face_area_m2, verts, faces,
        face_gas_normals, face_material_id, patch_scale_m, *, absolute_tolerance,
        relative_tolerance, confidence_level=0.95, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None,
        periodic_domain_lengths_m=None, periodic_domain_origin_m=None,
        minimum_mean_support_fraction=DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION):
    """Score independent face-field replicates on one exact-overlap physical patch grid.

    Replicates must be independent estimator scrambles, not nested levels from one scramble.  The
    triangle-to-patch geometry is constructed once from the authority field and then reused exactly
    for every replicate, so only the sampled physics varies.
    """
    authority_field = np.asarray(authority_face_field, dtype=float)
    replicates = np.asarray(replicate_face_fields, dtype=float)
    confidence = float(confidence_level)
    absolute = float(absolute_tolerance)
    relative = float(relative_tolerance)
    if (authority_field.ndim != 1 or replicates.ndim != 2
            or replicates.shape[1:] != authority_field.shape
            or replicates.shape[0] < 4 or np.any(~np.isfinite(authority_field))
            or np.any(~np.isfinite(replicates))
            or not np.isfinite(confidence) or not 0.0 < confidence < 1.0
            or not np.isfinite(absolute) or absolute <= 0.0
            or not np.isfinite(relative) or relative < 0.0):
        raise ValueError(
            "replicated patch scoring requires at least four compatible finite fields")
    authority = aggregate_surface_field_on_physical_patches_3d(
        face_field=authority_field, face_area_m2=face_area_m2, verts=verts,
        faces=faces, face_gas_normals=face_gas_normals,
        face_material_id=face_material_id, patch_scale_m=patch_scale_m,
        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
        patch_origin_m=patch_origin_m,
        periodic_domain_lengths_m=periodic_domain_lengths_m,
        periodic_domain_origin_m=periodic_domain_origin_m)
    support = _mean_support_partition_3d(
        authority, minimum_mean_support_fraction)
    integrated = np.empty((replicates.shape[0], len(authority.patch_key)))
    for replicate_index, field in enumerate(replicates):
        integrated[replicate_index] = np.bincount(
            authority.contribution_patch_index,
            weights=(
                authority.contribution_area_m2
                * field[authority.contribution_face_index]),
            minlength=len(authority.patch_key))
    mean_field = integrated / authority.patch_area_m2[None, :]
    critical = float(student_t.ppf(
        0.5 + 0.5 * confidence, replicates.shape[0] - 1))
    integrated_mean = np.mean(integrated, axis=0)
    field_mean = np.mean(mean_field, axis=0)
    integrated_half = critical * np.std(
        integrated, axis=0, ddof=1) / np.sqrt(replicates.shape[0])
    field_half = critical * np.std(
        mean_field, axis=0, ddof=1) / np.sqrt(replicates.shape[0])
    integrated_scale = (
        absolute * authority.patch_nominal_projected_area_m2
        + relative * np.maximum(
            np.abs(authority.integrated_field_area), np.abs(integrated_mean)))
    field_scale = (
        absolute
        + relative * np.maximum(np.abs(authority.mean_field), np.abs(field_mean)))
    integrated_bias = np.abs(
        authority.integrated_field_area - integrated_mean)
    field_bias = np.abs(authority.mean_field - field_mean)

    def maximum_normalized(numerator, denominator, mask=None):
        normalized = numerator / denominator
        if mask is not None:
            normalized = normalized[np.asarray(mask, dtype=bool)]
        return float(np.max(normalized, initial=0.0))

    return ReplicatedPhysicalPatchScore3D(
        authority=authority, replicate_count=replicates.shape[0],
        confidence_level=confidence, student_t_critical=critical,
        absolute_tolerance=absolute, relative_tolerance=relative,
        minimum_mean_support_fraction=support["threshold"],
        mean_eligible_patch_mask=support["eligible"],
        eligible_mean_patch_count=support["eligible_count"],
        excluded_mean_patch_count=support["excluded_count"],
        excluded_mean_surface_area_m2=support["excluded_surface_area_m2"],
        excluded_mean_surface_area_fraction=support[
            "excluded_surface_area_fraction"],
        excluded_mean_projected_support_area_m2=support[
            "excluded_projected_support_area_m2"],
        excluded_mean_projected_support_fraction=support[
            "excluded_projected_support_fraction"],
        replicate_mean_integrated_field_area=integrated_mean,
        replicate_mean_field=field_mean,
        integrated_confidence_half_width=integrated_half,
        mean_confidence_half_width=field_half,
        maximum_integrated_confidence_mixed_normalized=maximum_normalized(
            integrated_half, integrated_scale),
        maximum_mean_confidence_mixed_normalized=maximum_normalized(
            field_half, field_scale, support["eligible"]),
        maximum_integrated_authority_bias_mixed_normalized=maximum_normalized(
            integrated_bias, integrated_scale),
        maximum_mean_authority_bias_mixed_normalized=maximum_normalized(
            field_bias, field_scale, support["eligible"]),
        maximum_integrated_combined_mixed_normalized=maximum_normalized(
            integrated_half + integrated_bias, integrated_scale),
        maximum_mean_combined_mixed_normalized=maximum_normalized(
            field_half + field_bias, field_scale, support["eligible"]),
        maximum_mean_combined_all_patches_mixed_normalized=maximum_normalized(
            field_half + field_bias, field_scale))


@dataclass(frozen=True)
class PhysicalPatchRefinementScore3D:
    """Mixed absolute/relative discrepancy at one declared physical scale."""

    patch_scale_m: float
    patch_count: int
    maximum_patch_area_relative_error: float
    integrated_absolute_linf: float
    integrated_mixed_normalized_linf: float
    mean_absolute_linf: float
    mean_mixed_normalized_linf: float
    mean_all_patch_mixed_normalized_linf: float
    minimum_mean_support_fraction: float
    eligible_mean_patch_count: int
    excluded_mean_patch_count: int
    excluded_mean_surface_area_fraction: float
    excluded_mean_projected_support_fraction: float
    common_patch_surface_area_m2: np.ndarray
    common_patch_projected_support_area_m2: np.ndarray
    common_patch_projected_support_fraction: np.ndarray
    patch_nominal_projected_area_m2: np.ndarray
    mean_mixed_normalized_by_patch: np.ndarray
    reference_maximum_absolute_face_value: float
    candidate_maximum_absolute_face_value: float
    absolute_tolerance: float
    relative_tolerance: float

    def __post_init__(self):
        arrays = {}
        for name in (
                "common_patch_surface_area_m2",
                "common_patch_projected_support_area_m2",
                "common_patch_projected_support_fraction",
                "patch_nominal_projected_area_m2",
                "mean_mixed_normalized_by_patch"):
            array = _readonly(getattr(self, name), float)
            if (array.shape != (int(self.patch_count),)
                    or np.any(~np.isfinite(array)) or np.any(array < 0.0)):
                raise ValueError("invalid physical patch-refinement array")
            arrays[name] = array
        values = np.asarray([
            self.patch_scale_m, self.maximum_patch_area_relative_error,
            self.integrated_absolute_linf, self.integrated_mixed_normalized_linf,
            self.mean_absolute_linf, self.mean_mixed_normalized_linf,
            self.mean_all_patch_mixed_normalized_linf,
            self.minimum_mean_support_fraction,
            self.excluded_mean_surface_area_fraction,
            self.excluded_mean_projected_support_fraction,
            self.reference_maximum_absolute_face_value,
            self.candidate_maximum_absolute_face_value,
            self.absolute_tolerance, self.relative_tolerance], dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or int(self.patch_count) != self.patch_count or self.patch_count <= 0
                or int(self.eligible_mean_patch_count) != self.eligible_mean_patch_count
                or int(self.excluded_mean_patch_count) != self.excluded_mean_patch_count
                or self.eligible_mean_patch_count <= 0
                or self.excluded_mean_patch_count < 0
                or self.eligible_mean_patch_count + self.excluded_mean_patch_count
                    != self.patch_count
                or not 0.0 < self.minimum_mean_support_fraction <= 1.0
                or self.excluded_mean_surface_area_fraction > 1.0 + 1e-14
                or self.excluded_mean_projected_support_fraction > 1.0 + 1e-14):
            raise ValueError("invalid physical patch-refinement score")
        if (np.any(arrays["common_patch_surface_area_m2"] <= 0.0)
                or np.any(arrays["common_patch_projected_support_area_m2"] <= 0.0)
                or np.any(arrays["common_patch_projected_support_fraction"] <= 0.0)
                or np.any(arrays["patch_nominal_projected_area_m2"] <= 0.0)):
            raise ValueError("physical patch-refinement support must be positive")
        for name, array in arrays.items():
            object.__setattr__(self, name, array)


def compare_physical_patch_fields_3d(
        reference, candidate, *, absolute_tolerance, relative_tolerance,
        minimum_mean_support_fraction=DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION):
    """Compare two independently triangulated patch receipts without tiny-denominator blowup.

    The mixed normalized error is ``abs(a-b)/(atol + rtol*max(abs(a),abs(b)))``.  A declared
    tolerance is met when this value is at most one.  Patch-key disagreement is a topology or
    coverage change and is refused rather than hidden by intersecting only convenient patches.
    """
    if not isinstance(reference, PhysicalPatchField3D) or not isinstance(
            candidate, PhysicalPatchField3D):
        raise TypeError("reference and candidate must be PhysicalPatchField3D")
    if (not np.isfinite(absolute_tolerance) or absolute_tolerance <= 0.0
            or not np.isfinite(relative_tolerance) or relative_tolerance < 0.0
            ):
        raise ValueError("a positive absolute tolerance and nonnegative relative tolerance are required")
    if not np.isclose(
            reference.patch_scale_m, candidate.patch_scale_m, rtol=0.0, atol=0.0):
        raise ValueError("physical patch scales differ")
    if not np.array_equal(reference.patch_key, candidate.patch_key):
        raise ValueError("physical patch keys differ; topology or patch coverage changed")
    if (not np.array_equal(
            reference.periodic_domain_lengths_m, candidate.periodic_domain_lengths_m)
            or not np.array_equal(
                reference.periodic_domain_origin_m, candidate.periodic_domain_origin_m)):
        raise ValueError("physical patch periodic-domain contracts differ")

    area_difference = np.abs(reference.patch_area_m2 - candidate.patch_area_m2)
    area_scale = np.maximum(reference.patch_area_m2, candidate.patch_area_m2)
    area_relative = area_difference / area_scale

    integrated_difference = np.abs(
        reference.integrated_field_area - candidate.integrated_field_area)
    integrated_scale = (
        float(absolute_tolerance) * np.maximum(
            reference.patch_nominal_projected_area_m2,
            candidate.patch_nominal_projected_area_m2)
        + float(relative_tolerance) * np.maximum(
            np.abs(reference.integrated_field_area),
            np.abs(candidate.integrated_field_area)))
    mean_difference = np.abs(reference.mean_field - candidate.mean_field)
    mean_scale = (
        float(absolute_tolerance)
        + float(relative_tolerance) * np.maximum(
            np.abs(reference.mean_field), np.abs(candidate.mean_field)))
    threshold = float(minimum_mean_support_fraction)
    if not np.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise ValueError("minimum mean support fraction must be in (0, 1]")
    # A mean is comparable only where both representations resolve the declared footprint.  The
    # integrated discrepancy above still includes every patch, including excluded corner slivers.
    eligible = np.minimum(
        reference.patch_projected_support_fraction,
        candidate.patch_projected_support_fraction) >= threshold
    if not np.any(eligible):
        raise ValueError("no common physical patch meets the declared mean-support fraction")
    excluded = ~eligible
    common_surface_area = np.maximum(
        reference.patch_area_m2, candidate.patch_area_m2)
    common_projected_area = np.maximum(
        reference.patch_projected_support_area_m2,
        candidate.patch_projected_support_area_m2)
    mean_normalized = mean_difference / mean_scale
    return PhysicalPatchRefinementScore3D(
        patch_scale_m=reference.patch_scale_m, patch_count=len(reference.patch_key),
        maximum_patch_area_relative_error=float(np.max(area_relative, initial=0.0)),
        integrated_absolute_linf=float(np.max(integrated_difference, initial=0.0)),
        integrated_mixed_normalized_linf=float(np.max(
            integrated_difference / integrated_scale, initial=0.0)),
        mean_absolute_linf=float(np.max(mean_difference[eligible], initial=0.0)),
        mean_mixed_normalized_linf=float(np.max(
            mean_normalized[eligible], initial=0.0)),
        mean_all_patch_mixed_normalized_linf=float(np.max(
            mean_normalized, initial=0.0)),
        minimum_mean_support_fraction=threshold,
        eligible_mean_patch_count=int(np.count_nonzero(eligible)),
        excluded_mean_patch_count=int(np.count_nonzero(excluded)),
        excluded_mean_surface_area_fraction=float(
            np.sum(common_surface_area[excluded]) / np.sum(common_surface_area)),
        excluded_mean_projected_support_fraction=float(
            np.sum(common_projected_area[excluded]) / np.sum(common_projected_area)),
        common_patch_surface_area_m2=common_surface_area,
        common_patch_projected_support_area_m2=common_projected_area,
        common_patch_projected_support_fraction=np.minimum(
            reference.patch_projected_support_fraction,
            candidate.patch_projected_support_fraction),
        patch_nominal_projected_area_m2=np.maximum(
            reference.patch_nominal_projected_area_m2,
            candidate.patch_nominal_projected_area_m2),
        mean_mixed_normalized_by_patch=mean_normalized,
        reference_maximum_absolute_face_value=(
            reference.maximum_absolute_face_value),
        candidate_maximum_absolute_face_value=(
            candidate.maximum_absolute_face_value),
        absolute_tolerance=float(absolute_tolerance),
        relative_tolerance=float(relative_tolerance))


def score_surface_field_refinement_at_physical_scales_3d(
        reference_surface, candidate_surface, patch_scales_m, *,
        absolute_tolerance, relative_tolerance,
        minimum_mean_support_fraction=DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION):
    """Score one scalar field at no fewer than two fixed physical patch scales.

    Each surface mapping supplies the keyword arguments accepted by
    :func:`aggregate_surface_field_on_physical_patches_3d`, including ``face_field``.
    """
    scales = tuple(float(value) for value in patch_scales_m)
    if (len(scales) < 2 or len(set(scales)) != len(scales)
            or any(not np.isfinite(value) or value <= 0.0 for value in scales)):
        raise ValueError("at least two distinct positive physical patch scales are required")
    output = []
    for scale in scales:
        reference = aggregate_surface_field_on_physical_patches_3d(
            patch_scale_m=scale, **dict(reference_surface))
        candidate = aggregate_surface_field_on_physical_patches_3d(
            patch_scale_m=scale, **dict(candidate_surface))
        output.append(compare_physical_patch_fields_3d(
            reference, candidate, absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
            minimum_mean_support_fraction=minimum_mean_support_fraction))
    return tuple(output)


def score_replicated_surface_field_at_physical_scales_3d(
        authority_face_field, replicate_face_fields, face_area_m2, verts, faces,
        face_gas_normals, face_material_id, patch_scales_m, *, absolute_tolerance,
        relative_tolerance, confidence_level=0.95, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), patch_origin_m=None,
        periodic_domain_lengths_m=None, periodic_domain_origin_m=None,
        minimum_mean_support_fraction=DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION):
    """Apply the independent-replicate uncertainty gate at two or more physical scales."""
    scales = tuple(float(value) for value in patch_scales_m)
    if (len(scales) < 2 or len(set(scales)) != len(scales)
            or any(not np.isfinite(value) or value <= 0.0 for value in scales)):
        raise ValueError("at least two distinct positive physical patch scales are required")
    return tuple(
        score_replicated_surface_field_on_physical_patches_3d(
            authority_face_field, replicate_face_fields, face_area_m2, verts, faces,
            face_gas_normals, face_material_id, scale,
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
            confidence_level=confidence_level,
            mesh_length_unit_m=mesh_length_unit_m,
            mesh_origin_m=mesh_origin_m, patch_origin_m=patch_origin_m,
            periodic_domain_lengths_m=periodic_domain_lengths_m,
            periodic_domain_origin_m=periodic_domain_origin_m,
            minimum_mean_support_fraction=minimum_mean_support_fraction)
        for scale in scales)
