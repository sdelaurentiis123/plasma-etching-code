"""Orientation-local exact overlap transfer for piecewise-planar surfaces.

The single-patch authority in :mod:`petch.surface_overlap_remap_3d` is exact but deliberately
refuses a trench containing a floor, walls, and corner facets.  This module partitions each
immutable surface into material-local coplanar orientation groups, evaluates every geometrically
compatible old/new patch pair with that authority, and assembles one sparse global overlap matrix.

It is still a bounded operator: curved or rotating nearby faces refuse, large normal motion
refuses, and genuinely new or disappearing area must be handled by the caller's explicit closure.
No nearest-neighbor value is substituted when overlap geometry is unresolved.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from types import MappingProxyType

import numpy as np

from .surface_mesh_3d import TriangleSurface3D
from .surface_overlap_remap_3d import (
    SurfaceOverlapTransfer3D, build_surface_overlap_transfer_3d,
)


_SCHEMA = b"petch-partitioned-surface-overlap-3d-v1"


def _readonly(value, dtype):
    array = np.asarray(value, dtype=dtype).copy()
    array.setflags(write=False)
    return array


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


@dataclass(frozen=True)
class _PlanarPatch3D:
    face_index: np.ndarray
    material_id: int
    projection_axis: int
    orientation_sign: int
    unit_normal: np.ndarray
    plane_offset: float

    def __post_init__(self):
        face = _readonly(self.face_index, int)
        normal = _readonly(self.unit_normal, float)
        if (face.ndim != 1 or face.size == 0 or np.any(face < 0)
                or normal.shape != (3,)
                or not np.isclose(np.linalg.norm(normal), 1.0, atol=2e-14, rtol=0.0)
                or int(self.material_id) <= 0
                or int(self.projection_axis) not in (0, 1, 2)
                or int(self.orientation_sign) not in (-1, 1)
                or not np.isfinite(self.plane_offset)):
            raise ValueError("invalid planar surface patch")
        object.__setattr__(self, "face_index", face)
        object.__setattr__(self, "unit_normal", normal)


def _partition_planar_patches(surface, *, parallel_tolerance, planarity_tolerance):
    normal = _unit_normals(surface)
    triangle = surface.triangles
    projection = np.argmax(np.abs(normal), axis=1)
    orientation = np.sign(normal[np.arange(len(normal)), projection]).astype(int)
    unassigned = np.ones(len(surface.faces), dtype=bool)
    patches = []
    while np.any(unassigned):
        seed = int(np.flatnonzero(unassigned)[0])
        reference = normal[seed]
        offset = float(np.mean(triangle[seed] @ reference))
        material = int(surface.face_material_id[seed])
        candidate = (
            unassigned
            & (surface.face_material_id == material)
            & (projection == projection[seed])
            & (orientation == orientation[seed])
            & (1.0 - normal @ reference <= parallel_tolerance))
        plane_error = np.max(
            np.abs(np.einsum("fvc,c->fv", triangle, reference) - offset), axis=1)
        selected = np.flatnonzero(candidate & (plane_error <= planarity_tolerance))
        if selected.size == 0:  # pragma: no cover - the seed always satisfies its plane
            raise RuntimeError("planar partition lost its seed face")
        unassigned[selected] = False
        weighted = np.sum(normal[selected] * surface.face_area[selected, None], axis=0)
        weighted /= np.linalg.norm(weighted)
        patches.append(_PlanarPatch3D(
            selected, material, int(projection[seed]), int(orientation[seed]),
            weighted, float(np.mean(triangle[selected].reshape(-1, 3) @ weighted))))
    return tuple(patches)


def _compact_patch_surface(surface, patch, *, normal_axis_shift=0.0):
    global_faces = surface.faces[patch.face_index]
    used_vertex, inverse = np.unique(global_faces, return_inverse=True)
    vertices = np.array(surface.vertices[used_vertex], copy=True)
    vertices[:, patch.projection_axis] += float(normal_axis_shift)
    faces = inverse.reshape(global_faces.shape)
    lengths = list(surface.periodic_lengths)
    # Periodic images along the patch normal are alternative cells, not tangential overlap
    # contributors.  A nearest normal image is selected explicitly by the caller.
    lengths[patch.projection_axis] = None
    compact = TriangleSurface3D(
        vertices, faces, surface.face_material_id[patch.face_index],
        periodic_lengths=tuple(lengths), periodic_origin=surface.periodic_origin)
    return compact


@dataclass(frozen=True, eq=False)
class PartitionedSurfaceOverlapTransfer3D:
    """Global sparse overlap plus immutable orientation-patch receipts."""

    _sparse: SurfaceOverlapTransfer3D = field(repr=False)
    old_patch_count: int
    new_patch_count: int
    candidate_patch_pair_count: int
    positive_patch_pair_count: int
    patch_receipts: tuple
    fingerprint: str
    geometry_receipt: dict = field(repr=False)

    def __post_init__(self):
        if (not isinstance(self._sparse, SurfaceOverlapTransfer3D)
                or any(int(value) != value or int(value) < 0 for value in (
                    self.old_patch_count, self.new_patch_count,
                    self.candidate_patch_pair_count, self.positive_patch_pair_count))
                or self.positive_patch_pair_count > self.candidate_patch_pair_count
                or not isinstance(self.fingerprint, str) or len(self.fingerprint) != 64
                or self._sparse.fingerprint != self.fingerprint):
            raise ValueError("invalid partitioned surface-overlap transfer")
        object.__setattr__(self, "patch_receipts", tuple(
            MappingProxyType(dict(item)) for item in self.patch_receipts))
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


def build_partitioned_surface_overlap_transfer_3d(
        old_surface, new_surface, *, maximum_normal_distance,
        parallel_tolerance=1e-10, planarity_tolerance=None):
    """Build exact overlap on every compatible material/orientation-local planar patch."""
    if (not isinstance(old_surface, TriangleSurface3D)
            or not isinstance(new_surface, TriangleSurface3D)):
        raise TypeError("partitioned overlap requires TriangleSurface3D inputs")
    if (old_surface.periodic_lengths != new_surface.periodic_lengths
            or old_surface.periodic_origin != new_surface.periodic_origin):
        raise ValueError("partitioned overlap surfaces require one periodic cell")
    if set(old_surface.face_material_id.tolist()) != set(new_surface.face_material_id.tolist()):
        raise ValueError(
            "material surface appeared or disappeared; initialize/retire explicitly")
    maximum_normal_distance = float(maximum_normal_distance)
    parallel_tolerance = float(parallel_tolerance)
    scale = max(
        float(np.max(np.abs(old_surface.vertices), initial=0.0)),
        float(np.max(np.abs(new_surface.vertices), initial=0.0)), 1.0)
    if planarity_tolerance is None:
        planarity_tolerance = 1024.0 * np.finfo(float).eps * scale
    planarity_tolerance = float(planarity_tolerance)
    if (not np.isfinite(maximum_normal_distance) or maximum_normal_distance < 0.0
            or not np.isfinite(parallel_tolerance)
            or not 0.0 < parallel_tolerance <= 1e-4
            or not np.isfinite(planarity_tolerance) or planarity_tolerance <= 0.0):
        raise ValueError("invalid partitioned overlap tolerances")

    old_patch = _partition_planar_patches(
        old_surface, parallel_tolerance=parallel_tolerance,
        planarity_tolerance=planarity_tolerance)
    new_patch = _partition_planar_patches(
        new_surface, parallel_tolerance=parallel_tolerance,
        planarity_tolerance=planarity_tolerance)
    entry_new = []
    entry_old = []
    entry_area = []
    receipts = []
    candidate_count = 0
    positive_count = 0
    image_count = 0
    combined_image_count = 0
    identical_surface = old_surface.fingerprint == new_surface.fingerprint
    if identical_surface:
        entry_new.append(np.arange(len(new_surface.faces), dtype=int))
        entry_old.append(np.arange(len(old_surface.faces), dtype=int))
        entry_area.append(np.asarray(old_surface.face_area))
        candidate_count = len(old_patch)
        positive_count = len(old_patch)
        image_count = len(old_surface.faces)
        for patch_index, patch in enumerate(old_patch):
            receipts.append({
                "old_patch": patch_index,
                "new_patch": patch_index,
                "material_id": patch.material_id,
                "projection_axis": patch.projection_axis,
                "orientation_sign": patch.orientation_sign,
                "old_face_count": len(patch.face_index),
                "new_face_count": len(patch.face_index),
                "normal_image_shift": 0.0,
                "signed_normal_offset": 0.0,
                "overlap_area": float(np.sum(old_surface.face_area[patch.face_index])),
                "operator_fingerprint": old_surface.fingerprint,
            })
    for old_index, old in (() if identical_surface else enumerate(old_patch)):
        for new_index, new in enumerate(new_patch):
            if (old.material_id != new.material_id
                    or old.projection_axis != new.projection_axis
                    or old.orientation_sign != new.orientation_sign
                    or 1.0 - float(np.dot(old.unit_normal, new.unit_normal))
                    > parallel_tolerance):
                continue
            axis = old.projection_axis
            normal_period = old_surface.periodic_lengths[axis]
            normal_shift = 0.0
            raw_separation = float(
                np.mean(new_surface.triangles[new.face_index].reshape(-1, 3)
                        @ old.unit_normal)
                - np.mean(old_surface.triangles[old.face_index].reshape(-1, 3)
                          @ old.unit_normal))
            if normal_period is not None:
                axis_projection = old.unit_normal[axis]
                if abs(axis_projection) <= parallel_tolerance:
                    raise RuntimeError("dominant patch normal lost its projection axis")
                image = int(np.rint(raw_separation / (normal_period * axis_projection)))
                normal_shift = -image * normal_period
            shifted_separation = raw_separation + normal_shift * old.unit_normal[axis]
            if abs(shifted_separation) > maximum_normal_distance + planarity_tolerance:
                continue
            candidate_count += 1
            old_compact = _compact_patch_surface(old_surface, old)
            new_compact = _compact_patch_surface(
                new_surface, new, normal_axis_shift=normal_shift)
            local = build_surface_overlap_transfer_3d(
                old_compact, new_compact,
                projection_axis=axis, orientation_sign=old.orientation_sign,
                maximum_normal_distance=maximum_normal_distance,
                parallel_tolerance=parallel_tolerance)
            image_count += local.candidate_image_count
            combined_image_count += local.combined_periodic_image_count
            local_new = local.row_face_index
            if local.overlap_area.size:
                positive_count += 1
                entry_new.append(new.face_index[local_new])
                entry_old.append(old.face_index[local.old_face_index])
                entry_area.append(local.overlap_area)
            receipts.append({
                "old_patch": old_index,
                "new_patch": new_index,
                "material_id": old.material_id,
                "projection_axis": axis,
                "orientation_sign": old.orientation_sign,
                "old_face_count": len(old.face_index),
                "new_face_count": len(new.face_index),
                "normal_image_shift": normal_shift,
                "signed_normal_offset": local.signed_normal_offset,
                "overlap_area": float(np.sum(local.overlap_area)),
                "operator_fingerprint": local.fingerprint,
            })

    if entry_area:
        global_new = np.concatenate(entry_new)
        global_old = np.concatenate(entry_old)
        global_area = np.concatenate(entry_area)
        order = np.lexsort((global_old, global_new))
        global_new = global_new[order]
        global_old = global_old[order]
        global_area = global_area[order]
        start = np.flatnonzero(np.r_[
            True, (np.diff(global_new) != 0) | (np.diff(global_old) != 0)])
        combined_new = global_new[start]
        combined_old = global_old[start]
        combined_area = np.add.reduceat(global_area, start)
        combined_image_count += len(global_area) - len(combined_area)
    else:
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
    area_roundoff = 1024.0 * np.finfo(float).eps * max(
        float(np.max(old_surface.face_area)), float(np.max(new_surface.face_area)))
    if (np.any(old_covered > old_surface.face_area + area_roundoff)
            or np.any(new_covered > new_surface.face_area + area_roundoff)):
        raise RuntimeError(
            "partitioned overlap double-counted physical area across orientation patches")
    old_uncovered = np.maximum(old_surface.face_area - old_covered, 0.0)
    new_uncovered = np.maximum(new_surface.face_area - new_covered, 0.0)

    # A nearby face with a different normal is not "new surface"; it is unresolved curved or
    # rotating correspondence.  Refuse instead of laundering that numerical gap through the
    # caller's physically meaningful removal/exposure closure.
    old_normal = _unit_normals(old_surface)
    new_normal = _unit_normals(new_surface)
    for new_face in np.flatnonzero(new_covered <= area_roundoff):
        material = int(new_surface.face_material_id[new_face])
        nearest = old_surface.nearest(
            new_surface.face_centroid[new_face:new_face + 1],
            material_id=material,
            maximum_distance=max(maximum_normal_distance, planarity_tolerance))
        if (nearest.found[0]
                and 1.0 - float(np.dot(
                    old_normal[nearest.face_index[0]], new_normal[new_face]))
                > parallel_tolerance):
            raise ValueError(
                "nearby nonparallel surface requires a refinement-convergent curved "
                f"overlap operator; new_face={new_face}, old_face={nearest.face_index[0]}")
    for old_face in np.flatnonzero(old_covered <= area_roundoff):
        material = int(old_surface.face_material_id[old_face])
        nearest = new_surface.nearest(
            old_surface.face_centroid[old_face:old_face + 1],
            material_id=material,
            maximum_distance=max(maximum_normal_distance, planarity_tolerance))
        if (nearest.found[0]
                and 1.0 - float(np.dot(
                    old_normal[old_face], new_normal[nearest.face_index[0]]))
                > parallel_tolerance):
            raise ValueError(
                "nearby nonparallel surface requires a refinement-convergent curved "
                f"overlap operator; old_face={old_face}, new_face={nearest.face_index[0]}")

    digest = sha256()
    digest.update(_SCHEMA)
    digest.update(old_surface.fingerprint.encode("ascii"))
    digest.update(new_surface.fingerprint.encode("ascii"))
    _digest_array(digest, "row_offsets", offsets, "<i8")
    _digest_array(digest, "old_face_index", combined_old, "<i8")
    _digest_array(digest, "overlap_area", combined_area, "<f8")
    _digest_array(digest, "settings", [
        maximum_normal_distance, parallel_tolerance, planarity_tolerance], "<f8")
    fingerprint = digest.hexdigest()
    receipt = {
        "operator": "partitioned_exact_planar_triangle_overlap",
        "old_patch_count": len(old_patch),
        "new_patch_count": len(new_patch),
        "candidate_patch_pair_count": candidate_count,
        "positive_patch_pair_count": positive_count,
        "total_overlap_area": float(np.sum(combined_area)),
        "total_old_uncovered_area": float(np.sum(old_uncovered)),
        "total_new_uncovered_area": float(np.sum(new_uncovered)),
        "row_column_area_residual": abs(float(np.sum(old_covered) - np.sum(new_covered))),
        "parallel_tolerance": parallel_tolerance,
        "planarity_tolerance": planarity_tolerance,
        "limitation": (
            "piecewise-planar small-normal-motion patches; nearby rotation/curvature requires "
            "a refinement-convergent tangent-chart operator"),
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
        candidate_image_count=image_count,
        positive_image_overlap_count=len(combined_area),
        combined_periodic_image_count=combined_image_count,
        fingerprint=fingerprint, geometry_receipt=receipt)
    return PartitionedSurfaceOverlapTransfer3D(
        sparse, len(old_patch), len(new_patch), candidate_count, positive_count,
        tuple(receipts), fingerprint, receipt)
