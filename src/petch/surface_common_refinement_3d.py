"""Conservative common-refinement transfer for nearby moving triangle surfaces.

The old and new surfaces need not share connectivity or exactly coincident tangent planes.  An
indexed candidate search finds material-local, similarly oriented triangle pairs.  Each pair is
projected onto a well-conditioned shared coordinate plane and float64 convex clipping computes its
common-refinement area.  The physical coupling is the smaller of the source- and target-equivalent
areas, so no matrix entry can create surface capacity.

This is a first-order moving-surface operator, not a general mesh Boolean.  Large motion, opposite
or nearly orthogonal normals, ambiguous multiply covered sheets, and material appearance or
disappearance refuse.  Unmatched old area is explicitly removed; unmatched new area is explicitly
initialized by the caller.  The planar exact operator remains the certification reference.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from itertools import product
from types import MappingProxyType

import numpy as np
from scipy.spatial import cKDTree

from .surface_mesh_3d import TriangleSurface3D
from .surface_overlap_remap_3d import (
    SurfaceOverlapTransfer3D, _triangle_overlap_area_2d_batch,
)


_SCHEMA = b"petch-moving-surface-common-refinement-3d-v1"


def _digest_array(digest, name, value, dtype):
    array = np.ascontiguousarray(value, dtype=dtype)
    digest.update(str(name).encode("utf-8") + b"\0")
    digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
    digest.update(array.tobytes())


def _unit_normals(surface):
    triangle = surface.triangles
    normal = np.cross(
        triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    return normal / np.linalg.norm(normal, axis=1)[:, None]


def _periodic_shifts(periodic_lengths):
    choices = [
        (0.0,) if length is None else (-float(length), 0.0, float(length))
        for length in periodic_lengths]
    return np.asarray(tuple(product(*choices)), dtype=float)


@dataclass(frozen=True, eq=False)
class SurfaceCommonRefinementTransfer3D:
    """Sparse conservative coupling plus immutable moving-surface diagnostics."""

    _sparse: SurfaceOverlapTransfer3D = field(repr=False)
    candidate_pair_count: int
    aligned_pair_count: int
    positive_pair_image_count: int
    combined_pair_count: int
    minimum_normal_dot: float
    fingerprint: str
    geometry_receipt: dict = field(repr=False)

    def __post_init__(self):
        integer = (
            self.candidate_pair_count, self.aligned_pair_count,
            self.positive_pair_image_count, self.combined_pair_count)
        if (not isinstance(self._sparse, SurfaceOverlapTransfer3D)
                or any(int(value) != value or int(value) < 0 for value in integer)
                or self.aligned_pair_count > self.candidate_pair_count
                or self.positive_pair_image_count > self.aligned_pair_count
                or self.combined_pair_count > self.positive_pair_image_count
                or not np.isfinite(self.minimum_normal_dot)
                or not 0.0 < self.minimum_normal_dot <= 1.0
                or not isinstance(self.fingerprint, str) or len(self.fingerprint) != 64
                or self._sparse.fingerprint != self.fingerprint):
            raise ValueError("invalid common-refinement surface transfer")
        object.__setattr__(
            self, "geometry_receipt", MappingProxyType(dict(self.geometry_receipt)))

    @property
    def old_surface(self):
        return self._sparse.old_surface

    @property
    def new_surface(self):
        return self._sparse.new_surface

    @property
    def row_offsets(self):
        return self._sparse.row_offsets

    @property
    def old_face_index(self):
        return self._sparse.old_face_index

    @property
    def overlap_area(self):
        return self._sparse.overlap_area

    @property
    def old_covered_area(self):
        return self._sparse.old_covered_area

    @property
    def new_covered_area(self):
        return self._sparse.new_covered_area

    @property
    def old_uncovered_area(self):
        return self._sparse.old_uncovered_area

    @property
    def new_uncovered_area(self):
        return self._sparse.new_uncovered_area

    @property
    def row_face_index(self):
        return self._sparse.row_face_index

    def apply_extensive(self, old_density, *, newly_exposed_density):
        return self._sparse.apply_extensive(
            old_density, newly_exposed_density=newly_exposed_density)

    def apply_intensive(self, old_values, *, uncovered_fill=None):
        return self._sparse.apply_intensive(old_values, uncovered_fill=uncovered_fill)


def _candidate_pairs(old_surface, new_surface, *, maximum_normal_distance):
    """Return periodic old-image/new-face candidates using conservative bounding spheres."""
    old_triangle = old_surface.triangles
    new_triangle = new_surface.triangles
    old_radius = np.max(
        np.linalg.norm(old_triangle - old_surface.face_centroid[:, None, :], axis=2), axis=1)
    new_radius = np.max(
        np.linalg.norm(new_triangle - new_surface.face_centroid[:, None, :], axis=2), axis=1)
    shifts = _periodic_shifts(old_surface.periodic_lengths)
    old_image_centroid = np.concatenate([
        old_surface.face_centroid + shift for shift in shifts], axis=0)
    old_image_face = np.tile(np.arange(len(old_surface.faces), dtype=int), len(shifts))
    old_image_shift = np.repeat(shifts, len(old_surface.faces), axis=0)
    tree = cKDTree(old_image_centroid)
    maximum_old_radius = float(np.max(old_radius))
    pair_old = []
    pair_new = []
    pair_shift = []
    for new_face, (centroid, radius) in enumerate(zip(
            new_surface.face_centroid, new_radius)):
        image_index = tree.query_ball_point(
            centroid,
            float(maximum_normal_distance) + float(radius) + maximum_old_radius)
        for image in sorted(image_index):
            old_face = int(old_image_face[image])
            if (old_surface.face_material_id[old_face]
                    != new_surface.face_material_id[new_face]):
                continue
            pair_old.append(old_face)
            pair_new.append(new_face)
            pair_shift.append(old_image_shift[image])
    if not pair_old:
        return (
            np.empty(0, dtype=int), np.empty(0, dtype=int),
            np.empty((0, 3), dtype=float))
    return (
        np.asarray(pair_old, dtype=int), np.asarray(pair_new, dtype=int),
        np.asarray(pair_shift, dtype=float))


def build_surface_common_refinement_transfer_3d(
        old_surface, new_surface, *, maximum_normal_distance,
        minimum_normal_dot=0.5):
    """Build a sparse first-order common-refinement coupling for moving surfaces."""
    if (not isinstance(old_surface, TriangleSurface3D)
            or not isinstance(new_surface, TriangleSurface3D)):
        raise TypeError("common refinement requires TriangleSurface3D inputs")
    if (old_surface.periodic_lengths != new_surface.periodic_lengths
            or old_surface.periodic_origin != new_surface.periodic_origin):
        raise ValueError("common-refinement surfaces require one periodic cell")
    if set(old_surface.face_material_id.tolist()) != set(new_surface.face_material_id.tolist()):
        raise ValueError(
            "material surface appeared or disappeared; initialize/retire explicitly")
    maximum_normal_distance = float(maximum_normal_distance)
    minimum_normal_dot = float(minimum_normal_dot)
    if (not np.isfinite(maximum_normal_distance) or maximum_normal_distance < 0.0
            or not np.isfinite(minimum_normal_dot)
            or not 0.0 < minimum_normal_dot <= 1.0):
        raise ValueError("invalid common-refinement geometry controls")

    old_normal = _unit_normals(old_surface)
    new_normal = _unit_normals(new_surface)
    if old_surface.fingerprint == new_surface.fingerprint:
        combined_new = np.arange(len(new_surface.faces), dtype=int)
        combined_old = np.arange(len(old_surface.faces), dtype=int)
        combined_area = np.array(old_surface.face_area, dtype=float, copy=True)
        candidate_count = len(combined_area)
        aligned_count = len(combined_area)
        positive_image_count = len(combined_area)
    else:
        pair_old, pair_new, pair_shift = _candidate_pairs(
            old_surface, new_surface,
            maximum_normal_distance=maximum_normal_distance)
        candidate_count = len(pair_old)
        if candidate_count:
            normal_dot = np.einsum(
                "ij,ij->i", old_normal[pair_old], new_normal[pair_new])
            aligned = normal_dot >= minimum_normal_dot
            pair_old = pair_old[aligned]
            pair_new = pair_new[aligned]
            pair_shift = pair_shift[aligned]
        aligned_count = len(pair_old)

        entry_new = []
        entry_old = []
        entry_area = []
        for projection_axis in range(3):
            if not len(pair_old):
                break
            # Use one chart for every child of a source face.  Pair-specific charts can make two
            # target triangles cover overlapping images of the same old triangle and manufacture
            # source area.  The source-dominant chart makes its projected intersections a true
            # partition whenever the nearby target sheet is single-valued.
            selected_axis = np.argmax(np.abs(old_normal[pair_old]), axis=1)
            selected = np.flatnonzero(selected_axis == projection_axis)
            if not selected.size:
                continue
            old_index = pair_old[selected]
            new_index = pair_new[selected]
            shifted_old = (
                old_surface.triangles[old_index] + pair_shift[selected, None, :])
            target = new_surface.triangles[new_index]
            average_normal = old_normal[old_index] + new_normal[new_index]
            average_normal /= np.linalg.norm(average_normal, axis=1)[:, None]
            old_projection = np.einsum("fvc,fc->fv", shifted_old, average_normal)
            new_projection = np.einsum("fvc,fc->fv", target, average_normal)
            interval_gap = np.maximum.reduce([
                np.min(new_projection, axis=1) - np.max(old_projection, axis=1),
                np.min(old_projection, axis=1) - np.max(new_projection, axis=1),
                np.zeros(len(selected)),
            ])
            within_distance = interval_gap <= maximum_normal_distance
            if not np.any(within_distance):
                continue
            old_index = old_index[within_distance]
            new_index = new_index[within_distance]
            shifted_old = shifted_old[within_distance]
            target = target[within_distance]
            tangent_axis = tuple(axis for axis in range(3) if axis != projection_axis)
            projected = _triangle_overlap_area_2d_batch(
                shifted_old[:, :, tangent_axis], target[:, :, tangent_axis])
            old_factor = np.abs(old_normal[old_index, projection_axis])
            new_factor = np.abs(new_normal[new_index, projection_axis])
            positive = projected > 0.0
            if not np.any(positive):
                continue
            old_equivalent = projected[positive] / old_factor[positive]
            new_equivalent = projected[positive] / new_factor[positive]
            entry_new.append(new_index[positive])
            entry_old.append(old_index[positive])
            entry_area.append(np.minimum(old_equivalent, new_equivalent))

        if entry_area:
            global_new = np.concatenate(entry_new)
            global_old = np.concatenate(entry_old)
            global_area = np.concatenate(entry_area)
            positive_image_count = len(global_area)
            order = np.lexsort((global_old, global_new))
            global_new = global_new[order]
            global_old = global_old[order]
            global_area = global_area[order]
            start = np.flatnonzero(np.r_[
                True,
                (np.diff(global_new) != 0) | (np.diff(global_old) != 0),
            ])
            combined_new = global_new[start]
            combined_old = global_old[start]
            combined_area = np.add.reduceat(global_area, start)
        else:
            positive_image_count = 0
            combined_new = np.empty(0, dtype=int)
            combined_old = np.empty(0, dtype=int)
            combined_area = np.empty(0, dtype=float)

    offsets = np.zeros(len(new_surface.faces) + 1, dtype=int)
    if combined_new.size:
        offsets[1:] = np.cumsum(np.bincount(
            combined_new, minlength=len(new_surface.faces)))
    old_covered = np.bincount(
        combined_old, weights=combined_area,
        minlength=len(old_surface.faces)).astype(float)
    new_covered = np.bincount(
        combined_new, weights=combined_area,
        minlength=len(new_surface.faces)).astype(float)
    old_ratio = old_covered / old_surface.face_area
    new_ratio = new_covered / new_surface.face_area
    maximum_old_ratio = float(np.max(old_ratio, initial=0.0))
    maximum_new_ratio = float(np.max(new_ratio, initial=0.0))
    raw_overlap_area = float(np.sum(combined_area))
    # Curvature and periodic seam splitting can make either directed chart cover a face more than
    # once.  Project the geometric first guess symmetrically onto both capacity constraints.  Hard
    # support was already fixed by distance, material, and normal alignment; this projection only
    # decreases those local weights and therefore cannot create or teleport inventory.
    projection_iterations = 0
    if combined_area.size:
        for projection_iterations in range(65):
            old_covered = np.bincount(
                combined_old, weights=combined_area,
                minlength=len(old_surface.faces)).astype(float)
            new_covered = np.bincount(
                combined_new, weights=combined_area,
                minlength=len(new_surface.faces)).astype(float)
            maximum_capacity_ratio = max(
                float(np.max(old_covered / old_surface.face_area, initial=0.0)),
                float(np.max(new_covered / new_surface.face_area, initial=0.0)))
            if maximum_capacity_ratio <= 1.0 + 2e-14:
                break
            old_scale = np.minimum(
                1.0,
                np.sqrt(old_surface.face_area / np.maximum(
                    old_covered, np.finfo(float).tiny)))
            new_scale = np.minimum(
                1.0,
                np.sqrt(new_surface.face_area / np.maximum(
                    new_covered, np.finfo(float).tiny)))
            combined_area *= old_scale[combined_old] * new_scale[combined_new]
        else:  # pragma: no cover - monotone positive scaling must contract
            raise RuntimeError("common-refinement capacity projection did not contract")
        # Close any remaining representational excess exactly; each pass can only decrease the
        # opposite partition's sum, so one column/row pass preserves both constraints.
        old_covered = np.bincount(
            combined_old, weights=combined_area,
            minlength=len(old_surface.faces)).astype(float)
        combined_area *= np.minimum(
            1.0,
            old_surface.face_area / np.maximum(
                old_covered, np.finfo(float).tiny))[combined_old]
        new_covered = np.bincount(
            combined_new, weights=combined_area,
            minlength=len(new_surface.faces)).astype(float)
        combined_area *= np.minimum(
            1.0,
            new_surface.face_area / np.maximum(
                new_covered, np.finfo(float).tiny))[combined_new]
        old_covered = np.bincount(
            combined_old, weights=combined_area,
            minlength=len(old_surface.faces)).astype(float)
        new_covered = np.bincount(
            combined_new, weights=combined_area,
            minlength=len(new_surface.faces)).astype(float)
    old_uncovered = np.maximum(old_surface.face_area - old_covered, 0.0)
    new_uncovered = np.maximum(new_surface.face_area - new_covered, 0.0)

    digest = sha256()
    digest.update(_SCHEMA)
    digest.update(old_surface.fingerprint.encode("ascii"))
    digest.update(new_surface.fingerprint.encode("ascii"))
    _digest_array(digest, "row_offsets", offsets, "<i8")
    _digest_array(digest, "old_face_index", combined_old, "<i8")
    _digest_array(digest, "overlap_area", combined_area, "<f8")
    _digest_array(digest, "settings", [
        maximum_normal_distance, minimum_normal_dot], "<f8")
    fingerprint = digest.hexdigest()
    total_old_area = float(np.sum(old_surface.face_area))
    total_new_area = float(np.sum(new_surface.face_area))
    receipt = {
        "operator": "indexed_pairwise_tangent_common_refinement",
        "nonzero_authority": "float64_convex_polygon_clipping",
        "candidate_pair_count": candidate_count,
        "aligned_pair_count": aligned_count,
        "positive_pair_image_count": positive_image_count,
        "combined_pair_count": len(combined_area),
        "minimum_normal_dot": minimum_normal_dot,
        "maximum_normal_distance": maximum_normal_distance,
        "maximum_raw_old_coverage_ratio": maximum_old_ratio,
        "maximum_raw_new_coverage_ratio": maximum_new_ratio,
        "raw_overlap_area": raw_overlap_area,
        "capacity_projection_area_reduction": (
            raw_overlap_area - float(np.sum(combined_area))),
        "capacity_projection_iterations": projection_iterations,
        "old_matched_area_fraction": (
            float(np.sum(old_covered)) / total_old_area if total_old_area else 1.0),
        "new_matched_area_fraction": (
            float(np.sum(new_covered)) / total_new_area if total_new_area else 1.0),
        "total_overlap_area": float(np.sum(combined_area)),
        "total_old_uncovered_area": float(np.sum(old_uncovered)),
        "total_new_uncovered_area": float(np.sum(new_uncovered)),
        "limitation": (
            "first-order nearby similarly oriented surfaces; large motion, opposite sheets, "
            "and multiply covered correspondence refuse"),
    }
    sparse = SurfaceOverlapTransfer3D(
        old_surface=old_surface, new_surface=new_surface,
        row_offsets=offsets, old_face_index=combined_old,
        overlap_area=combined_area, old_covered_area=old_covered,
        new_covered_area=new_covered, old_uncovered_area=old_uncovered,
        new_uncovered_area=new_uncovered, projection_axis=0,
        orientation_sign=1, reference_unit_normal=np.asarray([1.0, 0.0, 0.0]),
        signed_normal_offset=0.0,
        maximum_normal_distance=maximum_normal_distance,
        candidate_image_count=candidate_count,
        positive_image_overlap_count=positive_image_count,
        combined_periodic_image_count=positive_image_count - len(combined_area),
        fingerprint=fingerprint, geometry_receipt=receipt)
    return SurfaceCommonRefinementTransfer3D(
        sparse, candidate_count, aligned_count, positive_image_count,
        len(combined_area), minimum_normal_dot, fingerprint, receipt)
