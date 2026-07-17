"""Immutable deterministic triangle-surface geometry for 3-D feature models.

The surface is the immutable geometry authority shared by feature evolution and
surface-state remapping.  It deliberately contains no process physics itself.
All coordinates use one caller-declared unit; returned distances and closest
points use that same unit.

Periodic queries use whole-triangle nearest images.  Therefore each stored
triangle must be geometrically continuous in the primary periodic cell; a
triangle whose vertices themselves straddle a seam must be unwrapped before it
is supplied here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from itertools import product
from types import MappingProxyType

import numpy as np
from scipy.spatial import cKDTree


_SURFACE_SCHEMA = b"petch-triangle-surface-3d-v1"


def _readonly_copy(value, dtype):
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


def _points_3d(value, *, name):
    points = np.asarray(value, dtype=float)
    if points.ndim == 1:
        points = points[None, :]
    if points.ndim != 2 or points.shape[1] != 3 or np.any(~np.isfinite(points)):
        raise ValueError(f"{name} requires finite points with shape (n, 3)")
    return np.array(points, dtype=float, copy=True, order="C")


def _closest_points_on_triangle_pairs(points, triangles):
    """Return exact distances for corresponding point/triangle pairs."""
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3 or triangles.ndim != 3
            or triangles.shape[1:] != (3, 3) or len(triangles) == 0
            or len(points) != len(triangles)
            or np.any(~np.isfinite(points)) or np.any(~np.isfinite(triangles))):
        raise ValueError("triangle distance requires finite nonempty 3-D geometry")

    a = triangles[:, 0]
    b = triangles[:, 1]
    c = triangles[:, 2]
    ab = b - a
    ac = c - a
    normal = np.cross(ab, ac)
    normal_squared = np.einsum("ij,ij->i", normal, normal)

    ap = points - a
    plane_parameter = np.einsum("ij,ij->i", ap, normal) / normal_squared
    projection = points - plane_parameter[:, None] * normal
    projected = projection - a
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", projected, ab)
    d21 = np.einsum("ij,ij->i", projected, ac)
    denominator = d00 * d11 - d01 * d01
    barycentric_b = (d11 * d20 - d01 * d21) / denominator
    barycentric_c = (d00 * d21 - d01 * d20) / denominator
    barycentric_a = 1.0 - barycentric_b - barycentric_c
    tolerance = 256.0 * np.finfo(float).eps
    inside = ((barycentric_a >= -tolerance)
              & (barycentric_b >= -tolerance)
              & (barycentric_c >= -tolerance))

    def segment_closest(start, end):
        edge = end - start
        edge_squared = np.einsum("ij,ij->i", edge, edge)
        parameter = np.einsum(
            "ij,ij->i", points - start, edge) / edge_squared
        parameter = np.clip(parameter, 0.0, 1.0)
        closest = start + parameter[:, None] * edge
        delta = points - closest
        return closest, np.einsum("ij,ij->i", delta, delta)

    closest_ab, distance_ab = segment_closest(a, b)
    closest_bc, distance_bc = segment_closest(b, c)
    closest_ca, distance_ca = segment_closest(c, a)
    edge_distances = np.stack((distance_ab, distance_bc, distance_ca), axis=1)
    edge_choice = np.argmin(edge_distances, axis=1)
    edge_points = np.stack((closest_ab, closest_bc, closest_ca), axis=1)
    selected_edge_point = edge_points[np.arange(len(triangles)), edge_choice]
    closest = np.where(inside[:, None], projection, selected_edge_point)
    delta = points - closest
    distance = np.sqrt(np.maximum(np.einsum("ij,ij->i", delta, delta), 0.0))
    return distance, closest


def _closest_points_on_triangles(point, triangles):
    """Return exact distances and closest points for one point and many triangles."""
    point = np.asarray(point, dtype=float)
    triangles = np.asarray(triangles, dtype=float)
    if point.shape != (3,) or triangles.ndim != 3:
        raise ValueError("triangle distance requires one 3-D point")
    return _closest_points_on_triangle_pairs(
        np.broadcast_to(point, (len(triangles), 3)), triangles)


@dataclass(frozen=True, eq=False)
class TriangleQueryResult3D:
    """Immutable nearest-surface result, preserving input-point order."""

    distance: np.ndarray
    face_index: np.ndarray
    closest_point: np.ndarray
    periodic_shift: np.ndarray
    candidate_count: np.ndarray

    def __post_init__(self):
        distance = _readonly_copy(self.distance, "<f8")
        face = _readonly_copy(self.face_index, "<i8")
        closest = _readonly_copy(self.closest_point, "<f8")
        shift = _readonly_copy(self.periodic_shift, "<f8")
        count = _readonly_copy(self.candidate_count, "<i8")
        size = len(distance)
        if (distance.shape != (size,) or face.shape != (size,)
                or closest.shape != (size, 3) or shift.shape != (size, 3)
                or count.shape != (size,) or np.any(count < 0)
                or np.any((face < 0) != ~np.isfinite(distance))
                or np.any(face < -1)):
            raise ValueError("invalid triangle-query result")
        found = face >= 0
        if (np.any(~np.isfinite(closest[found])) or np.any(~np.isfinite(shift[found]))
                or np.any(distance[found] < 0.0)
                or np.any(~np.isnan(closest[~found]))):
            raise ValueError("invalid triangle-query result values")
        object.__setattr__(self, "distance", distance)
        object.__setattr__(self, "face_index", face)
        object.__setattr__(self, "closest_point", closest)
        object.__setattr__(self, "periodic_shift", shift)
        object.__setattr__(self, "candidate_count", count)

    @property
    def found(self):
        value = self.face_index >= 0
        value.setflags(write=False)
        return value


@dataclass(frozen=True, eq=False)
class TriangleCentroidQueryResult3D:
    """K nearest physical face centroids with one periodic image per face."""

    distance: np.ndarray
    face_index: np.ndarray
    periodic_shift: np.ndarray
    candidate_count: np.ndarray

    def __post_init__(self):
        distance = _readonly_copy(self.distance, "<f8")
        face = _readonly_copy(self.face_index, "<i8")
        shift = _readonly_copy(self.periodic_shift, "<f8")
        count = _readonly_copy(self.candidate_count, "<i8")
        if (distance.ndim != 2 or face.shape != distance.shape
                or shift.shape != distance.shape + (3,)
                or count.shape != (distance.shape[0],)
                or np.any(~np.isfinite(distance)) or np.any(distance < 0.0)
                or np.any(face < 0) or np.any(~np.isfinite(shift))
                or np.any(count < distance.shape[1])):
            raise ValueError("invalid triangle-centroid query result")
        for row in face:
            if len(np.unique(row)) != len(row):
                raise ValueError("centroid query returned a physical face twice")
        object.__setattr__(self, "distance", distance)
        object.__setattr__(self, "face_index", face)
        object.__setattr__(self, "periodic_shift", shift)
        object.__setattr__(self, "candidate_count", count)


@dataclass(frozen=True, eq=False)
class TriangleImageCandidateResult3D:
    """Immutable CSR-like certified triangle-image candidate rows."""

    row_offsets: np.ndarray
    face_index: np.ndarray
    periodic_shift: np.ndarray
    centroid_distance: np.ndarray

    def __post_init__(self):
        offsets = _readonly_copy(self.row_offsets, "<i8")
        face = _readonly_copy(self.face_index, "<i8")
        shift = _readonly_copy(self.periodic_shift, "<f8")
        distance = _readonly_copy(self.centroid_distance, "<f8")
        if (offsets.ndim != 1 or len(offsets) == 0 or offsets[0] != 0
                or np.any(np.diff(offsets) < 0) or offsets[-1] != len(face)
                or shift.shape != (len(face), 3) or distance.shape != face.shape
                or np.any(face < 0) or np.any(~np.isfinite(shift))
                or np.any(~np.isfinite(distance)) or np.any(distance < 0.0)):
            raise ValueError("invalid triangle-image candidate result")
        object.__setattr__(self, "row_offsets", offsets)
        object.__setattr__(self, "face_index", face)
        object.__setattr__(self, "periodic_shift", shift)
        object.__setattr__(self, "centroid_distance", distance)


@dataclass(frozen=True, eq=False)
class _TriangleImageSpatialIndex3D:
    """Non-authoritative candidate index over explicit periodic triangle images."""

    triangles: np.ndarray
    centroids: np.ndarray
    face_index: np.ndarray
    periodic_shift: np.ndarray
    radius: np.ndarray
    tree: cKDTree
    maximum_radius: float


@dataclass(frozen=True, eq=False)
class TriangleSurface3D:
    """Exact immutable triangle surface with deterministic nearest queries.

    The fingerprint is intentionally order-sensitive: per-face physical state
    uses face order as an identity contract.  Numerically identical arrays have
    the same fingerprint independent of input dtype or memory layout.
    """

    vertices: np.ndarray
    faces: np.ndarray
    face_material_id: np.ndarray
    periodic_lengths: tuple[float | None, float | None, float | None] = (
        None, None, None)
    periodic_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    triangles: np.ndarray = field(init=False, repr=False)
    face_area: np.ndarray = field(init=False)
    face_centroid: np.ndarray = field(init=False)
    face_radius: np.ndarray = field(init=False, repr=False)
    fingerprint: str = field(init=False)
    _periodic_shifts: np.ndarray = field(init=False, repr=False)
    _spatial_indices: MappingProxyType = field(init=False, repr=False)

    def __post_init__(self):
        vertices = np.array(self.vertices, dtype="<f8", copy=True, order="C")
        faces = np.array(self.faces, dtype="<i8", copy=True, order="C")
        material = np.array(
            self.face_material_id, dtype="<i8", copy=True, order="C")
        if (vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 3
                or faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0
                or material.shape != (len(faces),)
                or np.any(~np.isfinite(vertices)) or np.any(faces < 0)
                or np.any(faces >= len(vertices)) or np.any(material <= 0)):
            raise ValueError("invalid triangle surface arrays")

        raw_lengths = tuple(self.periodic_lengths)
        raw_origin = tuple(self.periodic_origin)
        if len(raw_lengths) != 3 or len(raw_origin) != 3:
            raise ValueError("periodic metadata must have three axes")
        origin = tuple(float(value) for value in raw_origin)
        if np.any(~np.isfinite(origin)):
            raise ValueError("periodic origin must be finite")
        lengths = []
        choices = []
        for axis, supplied in enumerate(raw_lengths):
            if supplied is None:
                lengths.append(None)
                choices.append((0.0,))
                continue
            if isinstance(supplied, (bool, np.bool_)):
                raise ValueError("periodic lengths must be positive finite values")
            length = float(supplied)
            if not np.isfinite(length) or length <= 0.0:
                raise ValueError("periodic lengths must be positive finite values")
            scale = max(abs(origin[axis]), abs(origin[axis] + length),
                        float(np.max(np.abs(vertices[:, axis]))), 1.0)
            # Marching-cubes and Warp geometry are float32 producers even though this immutable
            # authority stores float64.  Accept only their roundoff-sized boundary overshoot, then
            # canonicalize it to the exact primary-cell endpoint.  Interior vertices are never
            # snapped, and a geometrically meaningful excursion still refuses.
            tolerance = max(
                256.0 * np.finfo(float).eps * scale,
                16.0 * np.finfo(np.float32).eps * scale)
            if (np.any(vertices[:, axis] < origin[axis] - tolerance)
                    or np.any(vertices[:, axis] > origin[axis] + length + tolerance)):
                raise ValueError("periodic surface vertices must lie in the primary cell")
            below = vertices[:, axis] < origin[axis]
            above = vertices[:, axis] > origin[axis] + length
            vertices[below, axis] = origin[axis]
            vertices[above, axis] = origin[axis] + length
            lengths.append(length)
            choices.append((-length, 0.0, length))
        lengths = tuple(lengths)
        shifts = np.asarray(tuple(product(*choices)), dtype="<f8")

        triangles = vertices[faces]
        edge_ab = triangles[:, 1] - triangles[:, 0]
        edge_ac = triangles[:, 2] - triangles[:, 0]
        edge_bc = triangles[:, 2] - triangles[:, 1]
        normal = np.cross(edge_ab, edge_ac)
        double_area = np.linalg.norm(normal, axis=1)
        maximum_edge_squared = np.maximum.reduce((
            np.einsum("ij,ij->i", edge_ab, edge_ab),
            np.einsum("ij,ij->i", edge_ac, edge_ac),
            np.einsum("ij,ij->i", edge_bc, edge_bc),
        ))
        degeneracy_floor = 64.0 * np.finfo(float).eps * np.maximum(
            maximum_edge_squared, np.finfo(float).tiny)
        if np.any(double_area <= degeneracy_floor):
            raise ValueError("triangle surface requires nondegenerate faces")
        area = 0.5 * double_area
        centroid = np.mean(triangles, axis=1)
        radius = np.max(
            np.linalg.norm(triangles - centroid[:, None, :], axis=2), axis=1)

        digest = sha256()
        digest.update(_SURFACE_SCHEMA)
        _digest_array(digest, "vertices", vertices, "<f8")
        _digest_array(digest, "faces", faces, "<i8")
        _digest_array(digest, "face_material_id", material, "<i8")
        _digest_array(
            digest, "periodic_axis_mask",
            [value is not None for value in lengths], "<i8")
        _digest_array(
            digest, "periodic_lengths",
            [0.0 if value is None else value for value in lengths], "<f8")
        _digest_array(digest, "periodic_origin", origin, "<f8")

        for array in (
                vertices, faces, material, triangles, area, centroid, radius, shifts):
            array.setflags(write=False)
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)
        object.__setattr__(self, "face_material_id", material)
        object.__setattr__(self, "periodic_lengths", lengths)
        object.__setattr__(self, "periodic_origin", origin)
        object.__setattr__(self, "triangles", triangles)
        object.__setattr__(self, "face_area", area)
        object.__setattr__(self, "face_centroid", centroid)
        object.__setattr__(self, "face_radius", radius)
        object.__setattr__(self, "_periodic_shifts", shifts)
        object.__setattr__(self, "fingerprint", digest.hexdigest())

        def make_index(face_selection):
            selected = np.asarray(face_selection, dtype=int)
            image_triangles = (
                triangles[selected, None, :, :] + shifts[None, :, None, :]
            ).reshape(-1, 3, 3)
            image_centroids = (
                centroid[selected, None, :] + shifts[None, :, :]
            ).reshape(-1, 3)
            image_faces = np.repeat(selected, len(shifts))
            image_shifts = np.tile(shifts, (len(selected), 1))
            image_radius = np.repeat(radius[selected], len(shifts))
            for array in (
                    image_triangles, image_centroids, image_faces,
                    image_shifts, image_radius):
                array.setflags(write=False)
            return _TriangleImageSpatialIndex3D(
                triangles=image_triangles,
                centroids=image_centroids,
                face_index=image_faces,
                periodic_shift=image_shifts,
                radius=image_radius,
                tree=cKDTree(image_centroids),
                maximum_radius=float(np.max(image_radius)),
            )

        all_faces = np.arange(len(faces), dtype=int)
        all_index = make_index(all_faces)
        material_indices = {None: all_index}
        unique_material = np.unique(material)
        for material_id in unique_material:
            selected = np.flatnonzero(material == material_id)
            material_indices[int(material_id)] = (
                all_index if len(selected) == len(faces) else make_index(selected))
        object.__setattr__(
            self, "_spatial_indices", MappingProxyType(material_indices))

    def _canonical_points(self, supplied):
        points = _points_3d(supplied, name="surface query")
        for axis, length in enumerate(self.periodic_lengths):
            if length is not None:
                points[:, axis] = (
                    self.periodic_origin[axis]
                    + np.mod(points[:, axis] - self.periodic_origin[axis], length))
        return points

    def _material_faces(self, material_id):
        if material_id is None:
            return np.arange(len(self.faces), dtype=int)
        if (isinstance(material_id, (bool, np.bool_))
                or int(material_id) != material_id or int(material_id) <= 0):
            raise ValueError("material_id must be a positive integer")
        selected = np.flatnonzero(self.face_material_id == int(material_id))
        if selected.size == 0:
            raise ValueError(f"triangle surface has no material {int(material_id)}")
        return selected

    def _query(self, supplied, *, material_id, maximum_distance, bounded, indexed):
        points = self._canonical_points(supplied)
        face_selection = self._material_faces(material_id)
        if maximum_distance is None:
            if bounded:
                raise ValueError("bounded query requires maximum_distance")
            limit = None
        else:
            limit = float(maximum_distance)
            if not np.isfinite(limit) or limit <= 0.0:
                raise ValueError("maximum_distance must be positive and finite")

        spatial = self._spatial_indices[
            None if material_id is None else int(material_id)]
        image_triangles = spatial.triangles
        image_centroids = spatial.centroids
        image_faces = spatial.face_index
        image_shifts = spatial.periodic_shift
        image_radius = spatial.radius

        distance_out = np.full(len(points), np.inf, dtype=float)
        face_out = np.full(len(points), -1, dtype=int)
        closest_out = np.full((len(points), 3), np.nan, dtype=float)
        shift_out = np.zeros((len(points), 3), dtype=float)
        candidate_count = np.zeros(len(points), dtype=int)
        coordinate_scale = max(
            float(np.max(np.abs(points), initial=0.0)),
            float(np.max(np.abs(image_triangles), initial=0.0)), 1.0)
        roundoff = 256.0 * np.finfo(float).eps * coordinate_scale

        indexed_rows = None
        indexed_offsets = None
        indexed_distance = None
        indexed_closest = None
        if indexed and bounded and len(points):
            broad_rows = spatial.tree.query_ball_point(
                points, r=limit + spatial.maximum_radius + roundoff)
            indexed_rows = []
            for point, broad in zip(points, broad_rows):
                broad = np.asarray(broad, dtype=int)
                if broad.size:
                    center_distance = np.linalg.norm(
                        image_centroids[broad] - point[None, :], axis=1)
                    broad = broad[
                        center_distance
                        <= limit + image_radius[broad] + roundoff]
                indexed_rows.append(broad)
            indexed_offsets = np.zeros(len(points) + 1, dtype=int)
            indexed_offsets[1:] = np.cumsum(
                [len(candidate) for candidate in indexed_rows])
            if indexed_offsets[-1]:
                flat_candidate = np.concatenate(indexed_rows)
                pair_points = np.repeat(
                    points, np.diff(indexed_offsets), axis=0)
                indexed_distance, indexed_closest = (
                    _closest_points_on_triangle_pairs(
                        pair_points, image_triangles[flat_candidate]))

        for row, point in enumerate(points):
            if indexed and bounded:
                candidate = indexed_rows[row]
            elif indexed:
                # One centroid-nearest image establishes a finite exact upper
                # bound.  A triangle outside d_best + its enclosing radius
                # cannot improve that bound.  The tree only chooses this
                # certified superset; exact geometry and the legacy tie rule
                # below remain authoritative.
                _, seed = spatial.tree.query(point, k=1)
                seed = int(seed)
                seed_distance = float(_closest_points_on_triangles(
                    point, image_triangles[seed:seed + 1])[0][0])
                broad = np.asarray(spatial.tree.query_ball_point(
                    point,
                    r=seed_distance + spatial.maximum_radius + roundoff),
                    dtype=int)
                center_distance = np.linalg.norm(
                    image_centroids[broad] - point[None, :], axis=1)
                candidate = broad[
                    center_distance
                    <= seed_distance + image_radius[broad] + roundoff]
            elif bounded:
                center_distance = np.linalg.norm(image_centroids - point[None, :], axis=1)
                candidate = np.flatnonzero(
                    center_distance <= limit + image_radius + roundoff)
            else:
                candidate = np.arange(len(image_triangles), dtype=int)
            candidate_count[row] = len(candidate)
            if candidate.size == 0:
                continue
            if indexed and bounded:
                start, stop = indexed_offsets[row:row + 2]
                exact_distance = indexed_distance[start:stop]
                exact_closest = indexed_closest[start:stop]
            else:
                exact_distance, exact_closest = _closest_points_on_triangles(
                    point, image_triangles[candidate])
            # Distance is primary, followed by physical face order and then a
            # lexicographic image shift.  Exact ties are therefore replayable.
            order = np.lexsort((
                image_shifts[candidate, 2], image_shifts[candidate, 1],
                image_shifts[candidate, 0], image_faces[candidate], exact_distance))
            winner_local = int(order[0])
            winner = int(candidate[winner_local])
            winner_distance = float(exact_distance[winner_local])
            if limit is not None and winner_distance > limit + roundoff:
                continue
            distance_out[row] = winner_distance
            face_out[row] = int(image_faces[winner])
            closest_out[row] = exact_closest[winner_local]
            shift_out[row] = image_shifts[winner]
        return TriangleQueryResult3D(
            distance_out, face_out, closest_out, shift_out, candidate_count)

    def nearest(self, points, *, material_id=None, maximum_distance=None):
        """Return exact nearest triangles through a certified spatial index.

        The index only produces a superset of candidates.  Exact triangle
        distance remains primary, followed by physical face order and the
        lexicographic periodic image shift.  With ``maximum_distance``, a point
        with no face inside the bound has distance ``inf`` and face index ``-1``.
        """
        return self._query(
            points, material_id=material_id,
            maximum_distance=maximum_distance,
            bounded=maximum_distance is not None, indexed=True)

    def nearest_brute_force(self, points, *, material_id=None, maximum_distance=None):
        """Reference nearest query evaluating every eligible periodic image."""
        return self._query(
            points, material_id=material_id,
            maximum_distance=maximum_distance, bounded=False, indexed=False)

    def nearest_face_centroids(self, points, *, count, material_id=None):
        """Return exact K-nearest physical centroids through the image index.

        Each physical face appears at most once.  Its distance and shift come
        from its nearest periodic image; exact image ties use lexicographic
        shift order.  Physical faces are then ordered by distance and face id,
        matching an exhaustive all-image calculation.
        """
        points = self._canonical_points(points)
        face_selection = self._material_faces(material_id)
        if (isinstance(count, (bool, np.bool_))
                or int(count) != count or int(count) <= 0
                or int(count) > len(face_selection)):
            raise ValueError(
                "count must be a positive integer no larger than eligible faces")
        count = int(count)
        spatial = self._spatial_indices[
            None if material_id is None else int(material_id)]
        raw_count = min(
            len(spatial.centroids), count * len(self._periodic_shifts))
        _, initial = spatial.tree.query(points, k=raw_count)
        initial = np.asarray(initial, dtype=int)
        if raw_count == 1:
            initial = initial[:, None]

        output_distance = np.empty((len(points), count), dtype=float)
        output_face = np.empty((len(points), count), dtype=int)
        output_shift = np.empty((len(points), count, 3), dtype=float)
        candidate_count = np.empty(len(points), dtype=int)
        coordinate_scale = max(
            float(np.max(np.abs(points), initial=0.0)),
            float(np.max(np.abs(spatial.centroids), initial=0.0)), 1.0)
        roundoff = 256.0 * np.finfo(float).eps * coordinate_scale

        def reduce_images(point, candidate):
            candidate = np.asarray(candidate, dtype=int)
            distance = np.linalg.norm(
                spatial.centroids[candidate] - point[None, :], axis=1)
            # Group by physical face.  Within a face, choose exact distance and
            # then lexicographic shift, reproducing np.argmin on the declared
            # itertools.product image order even when the tree reorders ties.
            grouped = np.lexsort((
                spatial.periodic_shift[candidate, 2],
                spatial.periodic_shift[candidate, 1],
                spatial.periodic_shift[candidate, 0],
                distance,
                spatial.face_index[candidate],
            ))
            grouped_candidate = candidate[grouped]
            grouped_distance = distance[grouped]
            keep = np.r_[
                True,
                np.diff(spatial.face_index[grouped_candidate]) != 0,
            ]
            physical_candidate = grouped_candidate[keep]
            physical_distance = grouped_distance[keep]
            order = np.lexsort((
                spatial.face_index[physical_candidate], physical_distance))
            return physical_candidate[order], physical_distance[order]

        for row, point in enumerate(points):
            preliminary_image, preliminary_distance = reduce_images(
                point, initial[row])
            if len(preliminary_image) < count:
                raise RuntimeError(
                    "periodic centroid query lost physical source faces")
            kth_distance = float(preliminary_distance[count - 1])
            candidate = np.asarray(spatial.tree.query_ball_point(
                point, r=kth_distance + roundoff), dtype=int)
            physical_image, physical_distance = reduce_images(point, candidate)
            if len(physical_image) < count:
                raise RuntimeError(
                    "certified centroid query lost physical source faces")
            selected = physical_image[:count]
            output_distance[row] = physical_distance[:count]
            output_face[row] = spatial.face_index[selected]
            output_shift[row] = spatial.periodic_shift[selected]
            candidate_count[row] = len(candidate)
        return TriangleCentroidQueryResult3D(
            output_distance, output_face, output_shift, candidate_count)

    def certified_triangle_image_candidates(
            self, points, *, query_radius, material_id=None):
        """Return every image whose enclosing sphere can meet a query sphere.

        This is a candidate-only primitive.  Rows are deterministically sorted
        by physical face and lexicographic image shift.  Downstream geometry
        must still evaluate its exact predicate; the tree never owns a hit or
        overlap decision.
        """
        points = self._canonical_points(points)
        radius = np.asarray(query_radius, dtype=float)
        if radius.ndim == 0:
            radius = np.full(len(points), float(radius))
        if (radius.shape != (len(points),) or np.any(~np.isfinite(radius))
                or np.any(radius < 0.0)):
            raise ValueError(
                "query_radius must be finite nonnegative scalar or one per point")
        self._material_faces(material_id)
        spatial = self._spatial_indices[
            None if material_id is None else int(material_id)]
        coordinate_scale = max(
            float(np.max(np.abs(points), initial=0.0)),
            float(np.max(np.abs(spatial.centroids), initial=0.0)), 1.0)
        roundoff = 256.0 * np.finfo(float).eps * coordinate_scale
        broad_rows = spatial.tree.query_ball_point(
            points, r=radius + spatial.maximum_radius + roundoff)
        rows = []
        row_distance = []
        offsets = np.zeros(len(points) + 1, dtype=int)
        for row, (point, broad) in enumerate(zip(points, broad_rows)):
            broad = np.asarray(broad, dtype=int)
            if broad.size:
                distance = np.linalg.norm(
                    spatial.centroids[broad] - point[None, :], axis=1)
                keep = (
                    distance
                    <= radius[row] + spatial.radius[broad] + roundoff)
                broad = broad[keep]
                distance = distance[keep]
                order = np.lexsort((
                    spatial.periodic_shift[broad, 2],
                    spatial.periodic_shift[broad, 1],
                    spatial.periodic_shift[broad, 0],
                    spatial.face_index[broad],
                ))
                broad = broad[order]
                distance = distance[order]
            else:
                distance = np.empty(0, dtype=float)
            rows.append(broad)
            row_distance.append(distance)
            offsets[row + 1] = offsets[row] + len(broad)
        if offsets[-1]:
            image = np.concatenate(rows)
            face = spatial.face_index[image]
            shift = spatial.periodic_shift[image]
            distance = np.concatenate(row_distance)
        else:
            face = np.empty(0, dtype=int)
            shift = np.empty((0, 3), dtype=float)
            distance = np.empty(0, dtype=float)
        return TriangleImageCandidateResult3D(
            offsets, face, shift, distance)
