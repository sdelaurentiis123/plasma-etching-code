"""Read-only geometry-backend contract for the common 3-D feature engine.

This module is intentionally not wired into :func:`petch.feature_step_3d.advance_feature_step_3d`
yet.  It defines the narrow storage/geometry seam needed by future sparse and AMR backends, and a
uniform wrapper that delegates every geometric operation to the current dense reference operator.

Coordinates carrying the ``_mesh`` suffix use the declared mesh coordinate unit.  Coordinates
carrying the ``_m`` suffix are SI positions and include ``FeatureGeometry3D.mesh_origin_m``.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Protocol, runtime_checkable

import numpy as np
from scipy.ndimage import map_coordinates

from .feature_geometry_state_3d import FeatureGeometry3D, face_material_ids_3d
from .threed import extract_mesh_3d


_BACKEND_SCHEMA = "petch-feature-geometry-backend-3d-v1"


def _canonical_digest_array(digest, name, value, dtype):
    """Hash one named array without depending on native dtype or memory layout."""
    array = np.ascontiguousarray(value, dtype=dtype)
    encoded_name = str(name).encode("utf-8")
    digest.update(np.asarray([len(encoded_name)], dtype="<u8").tobytes())
    digest.update(encoded_name)
    digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
    digest.update(array.tobytes())


@dataclass(frozen=True)
class PeriodicGeometryMetadata3D:
    """Canonical duplicate-endpoint periodic metadata for one geometry backend."""

    axes: tuple[int, ...]
    lengths_mesh_units: tuple[float | None, float | None, float | None]
    lengths_m: tuple[float | None, float | None, float | None]
    duplicate_endpoint_planes: tuple[bool, bool, bool]

    def __post_init__(self):
        axes = tuple(int(axis) for axis in self.axes)
        mesh = tuple(
            None if value is None else float(value)
            for value in self.lengths_mesh_units)
        physical = tuple(None if value is None else float(value) for value in self.lengths_m)
        duplicate = tuple(bool(value) for value in self.duplicate_endpoint_planes)
        if (axes != tuple(sorted(set(axes))) or any(axis < 0 or axis >= 3 for axis in axes)
                or len(mesh) != 3 or len(physical) != 3 or len(duplicate) != 3
                or any((axis in axes) != duplicate[axis] for axis in range(3))
                or any((mesh[axis] is None) != (axis not in axes) for axis in range(3))
                or any((physical[axis] is None) != (axis not in axes) for axis in range(3))
                or any(value is not None and (not np.isfinite(value) or value <= 0.0)
                       for value in mesh + physical)):
            raise ValueError("invalid periodic geometry metadata")
        object.__setattr__(self, "axes", axes)
        object.__setattr__(self, "lengths_mesh_units", mesh)
        object.__setattr__(self, "lengths_m", physical)
        object.__setattr__(self, "duplicate_endpoint_planes", duplicate)


@dataclass(frozen=True)
class FeatureSurfaceMesh3D:
    """Immutable surface extracted by a geometry backend in mesh coordinates."""

    vertices_mesh: np.ndarray
    faces: np.ndarray
    centroids_mesh: np.ndarray
    areas_mesh2: np.ndarray
    face_material_id: np.ndarray
    mesh_length_unit_m: float
    mesh_origin_m: tuple[float, float, float]

    def __post_init__(self):
        vertices = np.asarray(self.vertices_mesh).copy()
        faces = np.asarray(self.faces, dtype=np.int32).copy()
        centroids = np.asarray(self.centroids_mesh, dtype=float).copy()
        areas = np.asarray(self.areas_mesh2, dtype=float).copy()
        material = np.asarray(self.face_material_id, dtype=int).copy()
        origin = tuple(float(value) for value in self.mesh_origin_m)
        if (vertices.ndim != 2 or vertices.shape[1] != 3
                or faces.ndim != 2 or faces.shape[1] != 3
                or centroids.shape != (len(faces), 3)
                or areas.shape != (len(faces),)
                or material.shape != (len(faces),)
                or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(centroids))
                or np.any(~np.isfinite(areas)) or np.any(areas <= 0.0)
                or np.any(faces < 0) or np.any(faces >= len(vertices))
                or np.any(material <= 0)
                or not np.isfinite(self.mesh_length_unit_m)
                or self.mesh_length_unit_m <= 0.0
                or len(origin) != 3 or np.any(~np.isfinite(origin))):
            raise ValueError("invalid extracted feature surface")
        for array in (vertices, faces, centroids, areas, material):
            array.setflags(write=False)
        object.__setattr__(self, "vertices_mesh", vertices)
        object.__setattr__(self, "faces", faces)
        object.__setattr__(self, "centroids_mesh", centroids)
        object.__setattr__(self, "areas_mesh2", areas)
        object.__setattr__(self, "face_material_id", material)
        object.__setattr__(self, "mesh_length_unit_m", float(self.mesh_length_unit_m))
        object.__setattr__(self, "mesh_origin_m", origin)

    @property
    def vertices_m(self):
        return (np.asarray(self.mesh_origin_m)[None, :]
                + self.vertices_mesh * self.mesh_length_unit_m)

    @property
    def centroids_m(self):
        return (np.asarray(self.mesh_origin_m)[None, :]
                + self.centroids_mesh * self.mesh_length_unit_m)

    @property
    def areas_m2(self):
        return self.areas_mesh2 * self.mesh_length_unit_m ** 2

    @property
    def fingerprint(self):
        digest = sha256()
        digest.update(b"petch-feature-surface-mesh-3d-v1")
        _canonical_digest_array(digest, "vertices", self.vertices_mesh, "<f8")
        _canonical_digest_array(digest, "faces", self.faces, "<i8")
        _canonical_digest_array(digest, "centroids", self.centroids_mesh, "<f8")
        _canonical_digest_array(digest, "areas", self.areas_mesh2, "<f8")
        _canonical_digest_array(digest, "material", self.face_material_id, "<i8")
        _canonical_digest_array(
            digest, "physical_metadata",
            [self.mesh_length_unit_m, *self.mesh_origin_m], "<f8")
        return digest.hexdigest()


@runtime_checkable
class FeatureGeometryBackend3D(Protocol):
    """Storage-neutral read-only contract required before sparse evolution is introduced."""

    @property
    def backend_kind(self) -> str: ...

    @property
    def shape(self) -> tuple[int, int, int]: ...

    @property
    def finest_spacing_mesh_units(self) -> float: ...

    @property
    def finest_spacing_m(self) -> float: ...

    @property
    def domain_extent_mesh_units(self) -> tuple[float, float, float]: ...

    @property
    def domain_bounds_m(self) -> tuple[
            tuple[float, float, float], tuple[float, float, float]]: ...

    @property
    def periodic_metadata(self) -> PeriodicGeometryMetadata3D: ...

    @property
    def fingerprint(self) -> str: ...

    def sample_signed_distance_mesh(
            self, points_mesh, *, material_id: int | None = None) -> np.ndarray: ...

    def sample_signed_distance_m(
            self, points_m, *, material_id: int | None = None) -> np.ndarray: ...

    def sample_material_owner_mesh(self, points_mesh) -> np.ndarray: ...

    def sample_material_owner_m(self, points_m) -> np.ndarray: ...

    def extract_surface(self) -> FeatureSurfaceMesh3D: ...


@dataclass(frozen=True)
class UniformFeatureGeometryBackend3D:
    """Read-only adapter around the current dense :class:`FeatureGeometry3D` authority."""

    geometry: FeatureGeometry3D
    periodic_axes: tuple[int, ...] = ()

    def __post_init__(self):
        if not isinstance(self.geometry, FeatureGeometry3D):
            raise TypeError("uniform backend requires FeatureGeometry3D")
        raw = tuple(self.periodic_axes)
        if any(isinstance(axis, (bool, np.bool_)) or int(axis) != axis for axis in raw):
            raise ValueError("periodic axes must contain unique integer axes")
        axes = tuple(sorted(int(axis) for axis in raw))
        if len(set(axes)) != len(axes) or any(axis < 0 or axis >= 3 for axis in axes):
            raise ValueError("periodic axes must contain unique axes from 0, 1, 2")
        object.__setattr__(self, "periodic_axes", axes)

    @property
    def backend_kind(self):
        return "uniform_dense_reference_v1"

    @property
    def shape(self):
        return tuple(int(value) for value in self.geometry.phi.shape)

    @property
    def finest_spacing_mesh_units(self):
        return float(self.geometry.dx)

    @property
    def finest_spacing_m(self):
        return float(self.geometry.dx * self.geometry.mesh_length_unit_m)

    @property
    def domain_extent_mesh_units(self):
        return tuple(float(value) for value in (
            (np.asarray(self.shape) - 1) * self.geometry.dx))

    @property
    def domain_extent_m(self):
        return tuple(
            float(value * self.geometry.mesh_length_unit_m)
            for value in self.domain_extent_mesh_units)

    @property
    def domain_bounds_m(self):
        lower = tuple(float(value) for value in self.geometry.mesh_origin_m)
        upper = tuple(
            lower[axis] + self.domain_extent_m[axis]
            for axis in range(3))
        return lower, upper

    @property
    def periodic_metadata(self):
        extent_mesh = self.domain_extent_mesh_units
        extent_m = self.domain_extent_m
        return PeriodicGeometryMetadata3D(
            self.periodic_axes,
            tuple(extent_mesh[axis] if axis in self.periodic_axes else None
                  for axis in range(3)),
            tuple(extent_m[axis] if axis in self.periodic_axes else None
                  for axis in range(3)),
            tuple(axis in self.periodic_axes for axis in range(3)))

    def _canonical_mesh_points(self, supplied):
        points = np.asarray(supplied, dtype=float)
        if points.ndim == 1:
            points = points[None, :]
        if points.ndim != 2 or points.shape[1] != 3 or np.any(~np.isfinite(points)):
            raise ValueError("geometry samples require finite points with shape (n,3)")
        points = points.copy()
        extent = np.asarray(self.domain_extent_mesh_units)
        tolerance = 256.0 * np.finfo(float).eps * max(
            float(np.max(extent)), float(np.max(np.abs(points), initial=0.0)), 1.0)
        for axis in range(3):
            if axis in self.periodic_axes:
                points[:, axis] = np.mod(points[:, axis], extent[axis])
            else:
                if np.any(points[:, axis] < -tolerance) or np.any(
                        points[:, axis] > extent[axis] + tolerance):
                    raise ValueError("nonperiodic geometry sample lies outside the domain")
                points[:, axis] = np.clip(points[:, axis], 0.0, extent[axis])
        return points

    def _field_for_material(self, material_id):
        if material_id is None:
            return self.geometry.phi
        if isinstance(material_id, (bool, np.bool_)) or int(material_id) != material_id:
            raise ValueError("material_id must be a positive integer")
        material_id = int(material_id)
        if material_id <= 0:
            raise ValueError("material_id must be a positive integer")
        layers = self.geometry.material_levelsets
        if layers is not None:
            if material_id not in layers:
                raise ValueError(f"geometry has no material level set {material_id}")
            return layers[material_id]
        solid_materials = tuple(
            int(value) for value in np.unique(self.geometry.material_id) if int(value) > 0)
        if solid_materials != (material_id,):
            raise ValueError("geometry has no authoritative level set for the requested material")
        return self.geometry.phi

    def sample_signed_distance_mesh(self, points_mesh, *, material_id=None):
        points = self._canonical_mesh_points(points_mesh)
        field = self._field_for_material(material_id)
        coordinate = (points / self.geometry.dx).T
        return np.asarray(map_coordinates(
            field, coordinate, order=1, mode="nearest", prefilter=False), dtype=float)

    def sample_signed_distance_m(self, points_m, *, material_id=None):
        points = np.asarray(points_m, dtype=float)
        if points.ndim == 1:
            points = points[None, :]
        if points.ndim != 2 or points.shape[1] != 3 or np.any(~np.isfinite(points)):
            raise ValueError("geometry samples require finite points with shape (n,3)")
        mesh = ((points - np.asarray(self.geometry.mesh_origin_m)[None, :])
                / self.geometry.mesh_length_unit_m)
        return (self.sample_signed_distance_mesh(mesh, material_id=material_id)
                * self.geometry.mesh_length_unit_m)

    def sample_material_owner_mesh(self, points_mesh):
        points = self._canonical_mesh_points(points_mesh)
        combined = self.sample_signed_distance_mesh(points)
        owner = np.zeros(len(points), dtype=int)
        solid = combined >= 0.0
        if not np.any(solid):
            return owner
        layers = self.geometry.material_levelsets
        if layers is not None:
            material_ids = np.asarray(tuple(layers), dtype=int)
            coordinate = (points[solid] / self.geometry.dx).T
            values = np.vstack([
                map_coordinates(
                    layers[int(material_id)], coordinate,
                    order=1, mode="nearest", prefilter=False)
                for material_id in material_ids])
            owner[solid] = material_ids[np.argmax(values, axis=0)]
            return owner
        index = np.rint(points[solid] / self.geometry.dx).astype(int)
        for axis in range(3):
            index[:, axis] = np.clip(index[:, axis], 0, self.shape[axis] - 1)
        selected = self.geometry.material_id[tuple(index.T)]
        if np.any(selected <= 0):
            raise RuntimeError("solid sample has no material owner")
        owner[solid] = selected
        return owner

    def sample_material_owner_m(self, points_m):
        points = np.asarray(points_m, dtype=float)
        if points.ndim == 1:
            points = points[None, :]
        if points.ndim != 2 or points.shape[1] != 3 or np.any(~np.isfinite(points)):
            raise ValueError("geometry samples require finite points with shape (n,3)")
        mesh = ((points - np.asarray(self.geometry.mesh_origin_m)[None, :])
                / self.geometry.mesh_length_unit_m)
        return self.sample_material_owner_mesh(mesh)

    def extract_surface(self):
        vertices, faces, centroids, areas = extract_mesh_3d(
            self.geometry.phi, self.geometry.dx)
        material = face_material_ids_3d(centroids, self.geometry)
        return FeatureSurfaceMesh3D(
            vertices, faces, centroids, areas, material,
            self.geometry.mesh_length_unit_m, self.geometry.mesh_origin_m)

    @property
    def fingerprint(self):
        digest = sha256()
        digest.update(_BACKEND_SCHEMA.encode("ascii"))
        digest.update(self.backend_kind.encode("ascii"))
        _canonical_digest_array(digest, "phi", self.geometry.phi, "<f8")
        _canonical_digest_array(
            digest, "material_id", self.geometry.material_id, "<i8")
        layers = self.geometry.material_levelsets
        material_ids = () if layers is None else tuple(sorted(layers))
        _canonical_digest_array(digest, "material_levelset_ids", material_ids, "<i8")
        for material_id in material_ids:
            _canonical_digest_array(
                digest, f"material_levelset_{material_id}",
                layers[material_id], "<f8")
        _canonical_digest_array(
            digest, "uniform_metadata",
            [self.geometry.dx, self.geometry.mesh_length_unit_m,
             *self.geometry.mesh_origin_m], "<f8")
        _canonical_digest_array(digest, "periodic_axes", self.periodic_axes, "<i8")
        return digest.hexdigest()
