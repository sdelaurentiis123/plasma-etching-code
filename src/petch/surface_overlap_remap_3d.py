"""Exact planar triangle-overlap transfer for immutable 3-D surfaces.

This is a deliberately narrow geometry authority.  It supports one coplanar,
consistently oriented patch per old/new surface pair, with the two patches
parallel and separated by at most half a characteristic face length.  It does
not approximate curved, nonparallel, folded, or large-normal-motion geometry;
those cases refuse and remain on the certified KNN reference path.

The sparse matrix stores physical overlap area from old faces (columns) to new
faces (rows).  Candidate pairs come from ``TriangleSurface3D``'s periodic
spatial index, but exact convex-polygon clipping owns every nonzero.  Periodic
images of one physical old/new face pair are combined into one matrix entry,
so fundamental-cell inventory can never be counted twice silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from types import MappingProxyType
from typing import Mapping

import numpy as np
from numba import njit

from .surface_mesh_3d import TriangleSurface3D


_OVERLAP_SCHEMA = b"petch-planar-surface-overlap-3d-v1"
_APPLICATION_SCHEMA = b"petch-planar-surface-overlap-application-3d-v1"
_MAXIMUM_NORMAL_MOTION_TO_FACE_SCALE = 0.5


def _readonly(value, dtype):
    array = np.array(value, dtype=dtype, copy=True, order="C")
    array.setflags(write=False)
    return array


def _digest_array(digest, name, value, dtype):
    array = np.ascontiguousarray(value, dtype=dtype)
    encoded = str(name).encode("utf-8")
    digest.update(np.asarray([len(encoded)], dtype="<u8").tobytes())
    digest.update(encoded)
    digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
    digest.update(array.tobytes())


def _digest_text(digest, name, value):
    encoded = str(value).encode("utf-8")
    _digest_array(digest, name, np.frombuffer(encoded, dtype=np.uint8), "u1")


@njit(cache=True)
def _triangle_overlap_area_2d_batch(subject_triangle, clip_triangle):
    """Sutherland-Hodgman overlap area for corresponding triangle pairs."""
    pair_count = subject_triangle.shape[0]
    output_area = np.zeros(pair_count, dtype=np.float64)
    for pair in range(pair_count):
        polygon = np.empty((8, 2), dtype=np.float64)
        workspace = np.empty((8, 2), dtype=np.float64)
        polygon[:3] = subject_triangle[pair]
        polygon_count = 3
        clip = clip_triangle[pair]
        clip_signed_area = (
            (clip[1, 0] - clip[0, 0]) * (clip[2, 1] - clip[0, 1])
            - (clip[1, 1] - clip[0, 1]) * (clip[2, 0] - clip[0, 0]))
        orientation = 1.0 if clip_signed_area >= 0.0 else -1.0
        for edge_index in range(3):
            if polygon_count == 0:
                break
            edge_start = clip[edge_index]
            edge_stop = clip[(edge_index + 1) % 3]
            edge_x = edge_stop[0] - edge_start[0]
            edge_y = edge_stop[1] - edge_start[1]
            previous = polygon[polygon_count - 1]
            previous_signed = orientation * (
                edge_x * (previous[1] - edge_start[1])
                - edge_y * (previous[0] - edge_start[0]))
            previous_inside = previous_signed >= 0.0
            output_count = 0
            for vertex_index in range(polygon_count):
                current = polygon[vertex_index]
                current_signed = orientation * (
                    edge_x * (current[1] - edge_start[1])
                    - edge_y * (current[0] - edge_start[0]))
                current_inside = current_signed >= 0.0
                if current_inside != previous_inside:
                    denominator = previous_signed - current_signed
                    if denominator != 0.0:
                        fraction = previous_signed / denominator
                        workspace[output_count, 0] = (
                            previous[0] + fraction * (current[0] - previous[0]))
                        workspace[output_count, 1] = (
                            previous[1] + fraction * (current[1] - previous[1]))
                        output_count += 1
                if current_inside:
                    workspace[output_count] = current
                    output_count += 1
                previous = current
                previous_signed = current_signed
                previous_inside = current_inside
            for vertex_index in range(output_count):
                polygon[vertex_index] = workspace[vertex_index]
            polygon_count = output_count
        if polygon_count >= 3:
            twice_area = 0.0
            for vertex_index in range(polygon_count):
                following = (vertex_index + 1) % polygon_count
                twice_area += (
                    polygon[vertex_index, 0] * polygon[following, 1]
                    - polygon[vertex_index, 1] * polygon[following, 0])
            output_area[pair] = 0.5 * abs(twice_area)
    return output_area


def _surface_unit_normals(surface):
    triangle = surface.triangles
    normal = np.cross(
        triangle[:, 1] - triangle[:, 0],
        triangle[:, 2] - triangle[:, 0])
    return normal / np.linalg.norm(normal, axis=1)[:, None]


def _expanded_new_value(value, size, *, name, nonnegative):
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        array = np.full(int(size), float(array))
    if (array.shape != (int(size),) or np.any(~np.isfinite(array))
            or (nonnegative and np.any(array < 0.0))):
        qualifier = " finite nonnegative" if nonnegative else " finite"
        raise ValueError(f"{name} must be scalar or one{qualifier} value per new face")
    return array


@dataclass(frozen=True, eq=False)
class SurfaceOverlapApplication3D:
    """Immutable transferred field and its coverage/inventory receipt."""

    values: np.ndarray
    semantics: str
    transfer_fingerprint: str
    application_fingerprint: str
    covered_fraction: np.ndarray
    material_ledger: Mapping[int, Mapping[str, float]]
    maximum_relative_balance_error: float
    metadata: Mapping[str, object] = field(repr=False)

    def __post_init__(self):
        values = _readonly(self.values, "<f8")
        covered = _readonly(self.covered_fraction, "<f8")
        if (values.ndim != 1 or covered.shape != values.shape
                or np.any(~np.isfinite(values))
                or np.any(~np.isfinite(covered)) or np.any(covered < 0.0)
                or np.any(covered > 1.0 + 5e-12)
                or self.semantics not in {"extensive", "intensive"}
                or not isinstance(self.transfer_fingerprint, str)
                or len(self.transfer_fingerprint) != 64
                or not isinstance(self.application_fingerprint, str)
                or len(self.application_fingerprint) != 64
                or not np.isfinite(self.maximum_relative_balance_error)
                or self.maximum_relative_balance_error < 0.0):
            raise ValueError("invalid surface-overlap application")
        ledger = MappingProxyType({
            int(material): MappingProxyType({
                str(name): float(value) for name, value in item.items()
            })
            for material, item in self.material_ledger.items()
        })
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "covered_fraction", covered)
        object.__setattr__(self, "material_ledger", ledger)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, eq=False)
class SurfaceOverlapTransfer3D:
    """Sparse exact-area old-face to new-face planar overlap matrix."""

    old_surface: TriangleSurface3D
    new_surface: TriangleSurface3D
    row_offsets: np.ndarray
    old_face_index: np.ndarray
    overlap_area: np.ndarray
    old_covered_area: np.ndarray
    new_covered_area: np.ndarray
    old_uncovered_area: np.ndarray
    new_uncovered_area: np.ndarray
    projection_axis: int
    orientation_sign: int
    reference_unit_normal: np.ndarray
    signed_normal_offset: float
    maximum_normal_distance: float
    candidate_image_count: int
    positive_image_overlap_count: int
    combined_periodic_image_count: int
    fingerprint: str
    geometry_receipt: Mapping[str, object] = field(repr=False)

    def __post_init__(self):
        if (not isinstance(self.old_surface, TriangleSurface3D)
                or not isinstance(self.new_surface, TriangleSurface3D)):
            raise TypeError("overlap transfer requires TriangleSurface3D inputs")
        offsets = _readonly(self.row_offsets, "<i8")
        old_face = _readonly(self.old_face_index, "<i8")
        overlap = _readonly(self.overlap_area, "<f8")
        old_covered = _readonly(self.old_covered_area, "<f8")
        new_covered = _readonly(self.new_covered_area, "<f8")
        old_uncovered = _readonly(self.old_uncovered_area, "<f8")
        new_uncovered = _readonly(self.new_uncovered_area, "<f8")
        normal = _readonly(self.reference_unit_normal, "<f8")
        old_count = len(self.old_surface.faces)
        new_count = len(self.new_surface.faces)
        integer_values = (
            self.projection_axis, self.orientation_sign,
            self.candidate_image_count, self.positive_image_overlap_count,
            self.combined_periodic_image_count)
        if (offsets.shape != (new_count + 1,) or offsets[0] != 0
                or np.any(np.diff(offsets) < 0) or offsets[-1] != len(old_face)
                or overlap.shape != old_face.shape
                or old_covered.shape != (old_count,)
                or old_uncovered.shape != (old_count,)
                or new_covered.shape != (new_count,)
                or new_uncovered.shape != (new_count,)
                or normal.shape != (3,) or np.any(~np.isfinite(normal))
                or not np.isclose(np.linalg.norm(normal), 1.0, rtol=0.0, atol=2e-14)
                or np.any(old_face < 0) or np.any(old_face >= old_count)
                or np.any(~np.isfinite(overlap)) or np.any(overlap <= 0.0)
                or np.any(~np.isfinite(old_covered)) or np.any(old_covered < 0.0)
                or np.any(~np.isfinite(new_covered)) or np.any(new_covered < 0.0)
                or np.any(~np.isfinite(old_uncovered)) or np.any(old_uncovered < 0.0)
                or np.any(~np.isfinite(new_uncovered)) or np.any(new_uncovered < 0.0)
                or any(int(value) != value for value in integer_values)
                or int(self.projection_axis) not in (0, 1, 2)
                or int(self.orientation_sign) not in (-1, 1)
                or min(integer_values[2:]) < 0
                or not np.isfinite(self.signed_normal_offset)
                or not np.isfinite(self.maximum_normal_distance)
                or self.maximum_normal_distance < 0.0
                or not isinstance(self.fingerprint, str) or len(self.fingerprint) != 64):
            raise ValueError("invalid planar surface-overlap transfer")
        row_face = np.repeat(np.arange(new_count), np.diff(offsets))
        if np.any(
                self.old_surface.face_material_id[old_face]
                != self.new_surface.face_material_id[row_face]):
            raise ValueError("surface overlap crosses a material boundary")
        for row in range(new_count):
            start, stop = offsets[row:row + 2]
            if len(np.unique(old_face[start:stop])) != stop - start:
                raise ValueError("periodic images were not combined by physical face")
        object.__setattr__(self, "row_offsets", offsets)
        object.__setattr__(self, "old_face_index", old_face)
        object.__setattr__(self, "overlap_area", overlap)
        object.__setattr__(self, "old_covered_area", old_covered)
        object.__setattr__(self, "new_covered_area", new_covered)
        object.__setattr__(self, "old_uncovered_area", old_uncovered)
        object.__setattr__(self, "new_uncovered_area", new_uncovered)
        object.__setattr__(self, "reference_unit_normal", normal)
        object.__setattr__(self, "projection_axis", int(self.projection_axis))
        object.__setattr__(self, "orientation_sign", int(self.orientation_sign))
        object.__setattr__(self, "candidate_image_count", int(self.candidate_image_count))
        object.__setattr__(
            self, "positive_image_overlap_count", int(self.positive_image_overlap_count))
        object.__setattr__(
            self, "combined_periodic_image_count", int(self.combined_periodic_image_count))
        object.__setattr__(
            self, "geometry_receipt", MappingProxyType(dict(self.geometry_receipt)))

    @property
    def row_face_index(self):
        value = np.repeat(
            np.arange(len(self.new_surface.faces)), np.diff(self.row_offsets))
        value.setflags(write=False)
        return value

    def _old_values(self, supplied, *, nonnegative):
        values = np.asarray(supplied, dtype=float)
        if (values.shape != (len(self.old_surface.faces),)
                or np.any(~np.isfinite(values))
                or (nonnegative and np.any(values < 0.0))):
            qualifier = " finite nonnegative" if nonnegative else " finite"
            raise ValueError(f"old values require one{qualifier} value per old face")
        return values

    def _application_fingerprint(self, semantics, old_values, fill, closure):
        digest = sha256()
        digest.update(_APPLICATION_SCHEMA)
        _digest_text(digest, "transfer", self.fingerprint)
        _digest_text(digest, "semantics", semantics)
        _digest_text(digest, "closure", closure)
        _digest_array(digest, "old_values", old_values, "<f8")
        _digest_array(digest, "fill", fill, "<f8")
        return digest.hexdigest()

    def apply_extensive(self, old_density, *, newly_exposed_density):
        """Transfer nonnegative inventory and close old/new ledgers separately."""
        old = self._old_values(old_density, nonnegative=True)
        fill = _expanded_new_value(
            newly_exposed_density, len(self.new_surface.faces),
            name="newly_exposed_density", nonnegative=True)
        row = self.row_face_index
        retained_by_new = np.zeros(len(self.new_surface.faces), dtype=float)
        np.add.at(
            retained_by_new, row,
            self.overlap_area * old[self.old_face_index])
        newly_exposed_by_new = self.new_uncovered_area * fill
        new_inventory = retained_by_new + newly_exposed_by_new
        if self.old_surface.fingerprint == self.new_surface.fingerprint:
            values = np.array(old, copy=True)
        else:
            values = new_inventory / self.new_surface.face_area
        removed_by_old = self.old_uncovered_area * old
        material_ledger = {}
        maximum_error = 0.0
        for material in sorted(set(self.old_surface.face_material_id.tolist())):
            old_selected = self.old_surface.face_material_id == material
            new_selected = self.new_surface.face_material_id == material
            entry_selected = (
                self.new_surface.face_material_id[row] == material)
            old_total = float(np.dot(
                old[old_selected], self.old_surface.face_area[old_selected]))
            retained = float(np.sum(
                self.overlap_area[entry_selected]
                * old[self.old_face_index[entry_selected]]))
            removed = float(np.sum(removed_by_old[old_selected]))
            exposed = float(np.sum(newly_exposed_by_new[new_selected]))
            new_total = float(np.dot(
                values[new_selected], self.new_surface.face_area[new_selected]))
            old_residual = old_total - retained - removed
            new_residual = new_total - retained - exposed
            scale = max(
                abs(old_total), abs(new_total), abs(retained), abs(removed),
                abs(exposed), np.finfo(float).tiny)
            relative = max(abs(old_residual), abs(new_residual)) / scale
            maximum_error = max(maximum_error, relative)
            material_ledger[int(material)] = {
                "old_inventory": old_total,
                "retained_inventory": retained,
                "removed_inventory": removed,
                "newly_exposed_inventory": exposed,
                "new_inventory": new_total,
                "old_balance_residual": old_residual,
                "new_balance_residual": new_residual,
                "relative_balance_error": relative,
                "removed_area": float(np.sum(self.old_uncovered_area[old_selected])),
                "newly_exposed_area": float(np.sum(
                    self.new_uncovered_area[new_selected])),
            }
        return SurfaceOverlapApplication3D(
            values=values,
            semantics="extensive",
            transfer_fingerprint=self.fingerprint,
            application_fingerprint=self._application_fingerprint(
                "extensive", old, fill, "caller_declared_newly_exposed_density"),
            covered_fraction=(
                self.new_covered_area / self.new_surface.face_area),
            material_ledger=material_ledger,
            maximum_relative_balance_error=maximum_error,
            metadata={
                "newly_exposed_closure": "caller_declared_density",
                "old_uncovered_closure": "removed_inventory",
                "matrix_semantics": "exact_physical_overlap_area",
            })

    def apply_intensive(self, old_values, *, uncovered_fill=None):
        """Apply convex overlap averages, refusing undeclared exposed area."""
        old = self._old_values(old_values, nonnegative=False)
        area_scale = np.maximum(self.new_surface.face_area, np.finfo(float).tiny)
        uncovered_significant = self.new_uncovered_area > 5e-12 * area_scale
        if uncovered_fill is None and np.any(uncovered_significant):
            raise ValueError(
                "intensive overlap transfer requires uncovered_fill for newly exposed area")
        if uncovered_fill is None:
            fill = np.zeros(len(self.new_surface.faces), dtype=float)
        else:
            fill = _expanded_new_value(
                uncovered_fill, len(self.new_surface.faces),
                name="uncovered_fill", nonnegative=False)
        row = self.row_face_index
        numerator = np.zeros(len(self.new_surface.faces), dtype=float)
        np.add.at(
            numerator, row,
            self.overlap_area * old[self.old_face_index])
        if uncovered_fill is None:
            denominator = self.new_covered_area
        else:
            numerator += self.new_uncovered_area * fill
            denominator = self.new_covered_area + self.new_uncovered_area
        if np.any(denominator <= 0.0):
            raise RuntimeError("intensive overlap transfer has an unresolved empty row")
        if self.old_surface.fingerprint == self.new_surface.fingerprint:
            values = np.array(old, copy=True)
        else:
            values = numerator / denominator

        tolerance = 128.0 * np.finfo(float).eps
        for new_face in range(len(self.new_surface.faces)):
            start, stop = self.row_offsets[new_face:new_face + 2]
            contributors = old[self.old_face_index[start:stop]]
            if uncovered_fill is not None and self.new_uncovered_area[new_face] > 0.0:
                contributors = np.r_[contributors, fill[new_face]]
            if len(contributors) == 0:
                raise RuntimeError("intensive overlap transfer has no declared contributor")
            scale = max(float(np.max(np.abs(contributors))), 1.0)
            if (values[new_face] < float(np.min(contributors)) - tolerance * scale
                    or values[new_face] > float(np.max(contributors)) + tolerance * scale):
                raise RuntimeError("intensive overlap transfer violated convex monotonicity")
        material_ledger = {}
        for material in sorted(set(self.old_surface.face_material_id.tolist())):
            old_selected = self.old_surface.face_material_id == material
            new_selected = self.new_surface.face_material_id == material
            material_ledger[int(material)] = {
                "old_minimum": float(np.min(old[old_selected])),
                "old_maximum": float(np.max(old[old_selected])),
                "new_minimum": float(np.min(values[new_selected])),
                "new_maximum": float(np.max(values[new_selected])),
                "covered_area": float(np.sum(self.new_covered_area[new_selected])),
                "newly_exposed_area": float(np.sum(
                    self.new_uncovered_area[new_selected])),
            }
        return SurfaceOverlapApplication3D(
            values=values,
            semantics="intensive",
            transfer_fingerprint=self.fingerprint,
            application_fingerprint=self._application_fingerprint(
                "intensive", old, fill,
                "covered_rows_only" if uncovered_fill is None
                else "caller_declared_uncovered_fill"),
            covered_fraction=(
                self.new_covered_area / self.new_surface.face_area),
            material_ledger=material_ledger,
            maximum_relative_balance_error=0.0,
            metadata={
                "uncovered_fill_closure": (
                    "not_required" if uncovered_fill is None
                    else "caller_declared_value"),
                "monotonicity": "convex_overlap_average",
                "matrix_semantics": "exact_physical_overlap_area",
            })


def _validate_planar_pair(
        old_surface, new_surface, *, projection_axis, orientation_sign,
        maximum_normal_distance, parallel_tolerance):
    if (not isinstance(old_surface, TriangleSurface3D)
            or not isinstance(new_surface, TriangleSurface3D)):
        raise TypeError("overlap transfer requires TriangleSurface3D inputs")
    if (old_surface.periodic_lengths != new_surface.periodic_lengths
            or old_surface.periodic_origin != new_surface.periodic_origin):
        raise ValueError("old and new overlap surfaces require one periodic cell")
    if (isinstance(projection_axis, (bool, np.bool_))
            or int(projection_axis) not in (0, 1, 2)):
        raise ValueError("projection_axis must be 0, 1, or 2")
    projection_axis = int(projection_axis)
    if (isinstance(orientation_sign, (bool, np.bool_))
            or int(orientation_sign) not in (-1, 1)):
        raise ValueError("orientation_sign must be -1 or +1")
    orientation_sign = int(orientation_sign)
    maximum_normal_distance = float(maximum_normal_distance)
    parallel_tolerance = float(parallel_tolerance)
    if (not np.isfinite(maximum_normal_distance) or maximum_normal_distance < 0.0
            or not np.isfinite(parallel_tolerance)
            or not 0.0 < parallel_tolerance <= 1e-4):
        raise ValueError("invalid planar-overlap geometry tolerances")
    if (set(old_surface.face_material_id.tolist())
            != set(new_surface.face_material_id.tolist())):
        raise ValueError(
            "material surface appeared or disappeared; initialize/retire explicitly")

    old_normal = _surface_unit_normals(old_surface)
    new_normal = _surface_unit_normals(new_surface)
    reference = old_normal[0]
    for name, normal in (("old", old_normal), ("new", new_normal)):
        if (np.any(np.argmax(np.abs(normal), axis=1) != projection_axis)
                or np.any(np.sign(normal[:, projection_axis]) != orientation_sign)
                or np.any(1.0 - normal @ reference > parallel_tolerance)):
            raise ValueError(
                f"{name} surface is not one declared oriented parallel patch")
    for axis, length in enumerate(old_surface.periodic_lengths):
        if length is not None and abs(reference[axis]) > parallel_tolerance:
            raise ValueError(
                "periodicity with a normal component is unsupported by planar overlap")

    coordinate_scale = max(
        float(np.max(np.abs(old_surface.vertices), initial=0.0)),
        float(np.max(np.abs(new_surface.vertices), initial=0.0)), 1.0)
    planarity_tolerance = 512.0 * np.finfo(float).eps * coordinate_scale
    old_offset = float(np.mean(old_surface.vertices @ reference))
    new_offset = float(np.mean(new_surface.vertices @ reference))
    if (np.max(np.abs(old_surface.vertices @ reference - old_offset))
            > planarity_tolerance
            or np.max(np.abs(new_surface.vertices @ reference - new_offset))
            > planarity_tolerance):
        raise ValueError(
            "general noncoplanar/curved overlap geometry is not supported")
    signed_offset = new_offset - old_offset
    separation = abs(signed_offset)
    if separation > maximum_normal_distance + planarity_tolerance:
        raise ValueError("parallel-surface normal motion exceeds the declared bound")
    characteristic_face_length = float(np.sqrt(min(
        np.median(old_surface.face_area),
        np.median(new_surface.face_area))))
    if (separation
            > _MAXIMUM_NORMAL_MOTION_TO_FACE_SCALE * characteristic_face_length
            + planarity_tolerance):
        raise ValueError(
            "large normal motion is unsupported by the planar overlap operator")
    return (
        projection_axis, orientation_sign, maximum_normal_distance,
        parallel_tolerance, reference, signed_offset, separation,
        planarity_tolerance)


def build_surface_overlap_transfer_3d(
        old_surface, new_surface, *, projection_axis, orientation_sign,
        maximum_normal_distance, parallel_tolerance=1e-10):
    """Build an exact planar overlap matrix or refuse unsupported geometry."""
    (projection_axis, orientation_sign, maximum_normal_distance,
     parallel_tolerance, reference, signed_offset, separation,
     planarity_tolerance) = _validate_planar_pair(
        old_surface, new_surface,
        projection_axis=projection_axis,
        orientation_sign=orientation_sign,
        maximum_normal_distance=maximum_normal_distance,
        parallel_tolerance=parallel_tolerance)
    tangent_axis = tuple(axis for axis in range(3) if axis != projection_axis)
    projection_scale = 1.0 / abs(reference[projection_axis])
    old_count = len(old_surface.faces)
    new_count = len(new_surface.faces)

    if old_surface.fingerprint == new_surface.fingerprint:
        offsets = np.arange(new_count + 1, dtype=int)
        old_face = np.arange(old_count, dtype=int)
        overlap = np.array(old_surface.face_area, copy=True)
        candidate_image_count = old_count
        positive_image_overlap_count = old_count
        combined_periodic_image_count = 0
    else:
        pair_new = []
        pair_old = []
        pair_shift = []
        candidate_image_count = 0
        for material in sorted(set(old_surface.face_material_id.tolist())):
            new_index = np.flatnonzero(
                new_surface.face_material_id == material)
            candidate = old_surface.certified_triangle_image_candidates(
                new_surface.face_centroid[new_index],
                query_radius=new_surface.face_radius[new_index] + separation,
                material_id=int(material))
            candidate_image_count += len(candidate.face_index)
            for local_row, new_face in enumerate(new_index):
                start, stop = candidate.row_offsets[local_row:local_row + 2]
                count = stop - start
                if count:
                    pair_new.extend([int(new_face)] * count)
                    pair_old.extend(candidate.face_index[start:stop].tolist())
                    pair_shift.extend(candidate.periodic_shift[start:stop].tolist())
        if pair_new:
            pair_new = np.asarray(pair_new, dtype=int)
            pair_old = np.asarray(pair_old, dtype=int)
            pair_shift = np.asarray(pair_shift, dtype=float)
            old_triangle = (
                old_surface.triangles[pair_old] + pair_shift[:, None, :])
            new_triangle = new_surface.triangles[pair_new]
            projected_overlap = _triangle_overlap_area_2d_batch(
                old_triangle[:, :, tangent_axis],
                new_triangle[:, :, tangent_axis])
            exact_overlap = projected_overlap * projection_scale
            area_floor = (
                64.0 * np.finfo(float).eps
                * np.minimum(
                    old_surface.face_area[pair_old],
                    new_surface.face_area[pair_new]))
            positive = exact_overlap > area_floor
            pair_new = pair_new[positive]
            pair_old = pair_old[positive]
            pair_shift = pair_shift[positive]
            exact_overlap = exact_overlap[positive]
        else:
            pair_new = np.empty(0, dtype=int)
            pair_old = np.empty(0, dtype=int)
            pair_shift = np.empty((0, 3), dtype=float)
            exact_overlap = np.empty(0, dtype=float)
        positive_image_overlap_count = len(exact_overlap)

        # Sort image contributions deterministically, then combine all images
        # of the same physical old/new face pair into one fundamental-cell entry.
        order = np.lexsort((
            pair_shift[:, 2] if len(pair_shift) else np.empty(0),
            pair_shift[:, 1] if len(pair_shift) else np.empty(0),
            pair_shift[:, 0] if len(pair_shift) else np.empty(0),
            pair_old, pair_new))
        pair_new = pair_new[order]
        pair_old = pair_old[order]
        exact_overlap = exact_overlap[order]
        if len(exact_overlap):
            group_start = np.r_[
                True,
                (np.diff(pair_new) != 0) | (np.diff(pair_old) != 0)]
            start = np.flatnonzero(group_start)
            combined_new = pair_new[start]
            combined_old = pair_old[start]
            combined_overlap = np.add.reduceat(exact_overlap, start)
            combined_periodic_image_count = int(
                len(exact_overlap) - len(combined_overlap))
        else:
            combined_new = np.empty(0, dtype=int)
            combined_old = np.empty(0, dtype=int)
            combined_overlap = np.empty(0, dtype=float)
            combined_periodic_image_count = 0
        offsets = np.zeros(new_count + 1, dtype=int)
        if len(combined_new):
            offsets[1:] = np.cumsum(np.bincount(
                combined_new, minlength=new_count))
        old_face = combined_old
        overlap = combined_overlap

    row_face = np.repeat(np.arange(new_count), np.diff(offsets))
    old_covered = np.bincount(
        old_face, weights=overlap, minlength=old_count).astype(float)
    new_covered = np.bincount(
        row_face, weights=overlap, minlength=new_count).astype(float)
    geometry_scale = max(
        float(np.max(old_surface.face_area)),
        float(np.max(new_surface.face_area)), np.finfo(float).tiny)
    roundoff_area = 512.0 * np.finfo(float).eps * geometry_scale
    old_area_tolerance = 5e-12 * old_surface.face_area + roundoff_area
    new_area_tolerance = 5e-12 * new_surface.face_area + roundoff_area
    if (np.any(old_covered > old_surface.face_area + old_area_tolerance)
            or np.any(new_covered
                      > new_surface.face_area + new_area_tolerance)):
        raise RuntimeError(
            "overlap area exceeds a physical face; overlapping facets or periodic "
            "double counting are unsupported")
    old_uncovered = np.maximum(old_surface.face_area - old_covered, 0.0)
    new_uncovered = np.maximum(new_surface.face_area - new_covered, 0.0)
    row_column_residual = abs(
        float(np.sum(old_covered)) - float(np.sum(new_covered)))

    digest = sha256()
    digest.update(_OVERLAP_SCHEMA)
    _digest_text(digest, "old_surface", old_surface.fingerprint)
    _digest_text(digest, "new_surface", new_surface.fingerprint)
    _digest_array(digest, "row_offsets", offsets, "<i8")
    _digest_array(digest, "old_face_index", old_face, "<i8")
    _digest_array(digest, "overlap_area", overlap, "<f8")
    _digest_array(digest, "old_covered_area", old_covered, "<f8")
    _digest_array(digest, "new_covered_area", new_covered, "<f8")
    _digest_array(digest, "projection", [projection_axis, orientation_sign], "<i8")
    _digest_array(
        digest, "geometry_settings",
        [maximum_normal_distance, parallel_tolerance,
         _MAXIMUM_NORMAL_MOTION_TO_FACE_SCALE], "<f8")
    fingerprint = digest.hexdigest()
    receipt = {
        "operator": "exact_planar_triangle_overlap",
        "candidate_authority": "TriangleSurface3D_periodic_spatial_index",
        "nonzero_authority": "float64_convex_polygon_clipping",
        "periodic_policy": "combine_images_by_physical_old_new_face_pair",
        "total_overlap_area": float(np.sum(overlap)),
        "total_old_uncovered_area": float(np.sum(old_uncovered)),
        "total_new_uncovered_area": float(np.sum(new_uncovered)),
        "row_column_area_residual": row_column_residual,
        "maximum_old_area_excess": float(np.max(
            old_covered - old_surface.face_area, initial=0.0)),
        "maximum_new_area_excess": float(np.max(
            new_covered - new_surface.face_area, initial=0.0)),
        "normal_motion_to_characteristic_face_scale_limit": (
            _MAXIMUM_NORMAL_MOTION_TO_FACE_SCALE),
        "limitation": (
            "one coplanar oriented patch; nonparallel, curved, folded, and "
            "large-normal-motion geometry refuses"),
    }
    return SurfaceOverlapTransfer3D(
        old_surface=old_surface,
        new_surface=new_surface,
        row_offsets=offsets,
        old_face_index=old_face,
        overlap_area=overlap,
        old_covered_area=old_covered,
        new_covered_area=new_covered,
        old_uncovered_area=old_uncovered,
        new_uncovered_area=new_uncovered,
        projection_axis=projection_axis,
        orientation_sign=orientation_sign,
        reference_unit_normal=reference,
        signed_normal_offset=signed_offset,
        maximum_normal_distance=maximum_normal_distance,
        candidate_image_count=candidate_image_count,
        positive_image_overlap_count=positive_image_overlap_count,
        combined_periodic_image_count=combined_periodic_image_count,
        fingerprint=fingerprint,
        geometry_receipt=receipt)
