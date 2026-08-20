"""Dependency-neutral dense geometry state for the common 3-D feature engine.

The class and face-ownership operator in this module are the exact uniform-grid authorities formerly
defined in :mod:`petch.feature_step_3d`.  Keeping them below both the feature-step orchestrator and
the geometry-backend adapter prevents an import cycle while preserving the legacy re-exports.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class FeatureGeometry3D:
    """Eulerian material geometry in declared mesh units; material id zero is gas."""

    phi: np.ndarray
    material_id: np.ndarray
    dx: float
    mesh_length_unit_m: float
    mesh_origin_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    material_levelsets: Mapping[int, np.ndarray] | None = None

    def __post_init__(self):
        phi = np.asarray(self.phi, dtype=float).copy()
        material = np.asarray(self.material_id, dtype=int).copy()
        origin = tuple(float(value) for value in self.mesh_origin_m)
        layers = (None if self.material_levelsets is None
                  else {int(key): np.asarray(value, dtype=float).copy()
                        for key, value in self.material_levelsets.items()})
        if (phi.ndim != 3 or min(phi.shape) < 2 or material.shape != phi.shape
                or np.any(~np.isfinite(phi)) or np.any(material < 0)
                or not np.isfinite(self.dx) or self.dx <= 0.0
                or not np.isfinite(self.mesh_length_unit_m) or self.mesh_length_unit_m <= 0.0
                or len(origin) != 3 or np.any(~np.isfinite(origin))):
            raise ValueError("invalid 3-D feature geometry")
        if not np.any(phi < 0.0) or not np.any(phi > 0.0):
            raise ValueError("phi must contain both gas and solid")
        if layers is not None:
            material_ids = set(np.unique(material)) - {0}
            if (not layers or any(key <= 0 for key in layers)
                    or any(value.shape != phi.shape or np.any(~np.isfinite(value))
                           for value in layers.values())):
                raise ValueError("invalid material level-set fields")
            if set(layers) != material_ids:
                ranges = {
                    int(key): (
                        float(np.min(value)), float(np.max(value)),
                        int(np.count_nonzero(value >= 0.0)),
                    )
                    for key, value in sorted(layers.items())
                }
                raise ValueError(
                    "invalid material level-set fields: registered ids "
                    f"{sorted(layers)} do not match labeled solid ids "
                    f"{sorted(material_ids)}; "
                    "ranges=(minimum, maximum, nonnegative_node_count) "
                    f"{ranges}"
                )
            union = np.maximum.reduce(tuple(layers.values()))
            if np.any((union >= 0.0) != (phi >= 0.0)):
                raise ValueError("material level sets do not reconstruct the combined solid")
            for value in layers.values():
                value.setflags(write=False)
        phi.setflags(write=False); material.setflags(write=False)
        object.__setattr__(self, "phi", phi)
        object.__setattr__(self, "material_id", material)
        object.__setattr__(self, "dx", float(self.dx))
        object.__setattr__(self, "mesh_length_unit_m", float(self.mesh_length_unit_m))
        object.__setattr__(self, "mesh_origin_m", origin)
        object.__setattr__(
            self, "material_levelsets",
            None if layers is None else MappingProxyType(layers))

    @property
    def coordinate_arrays(self):
        return tuple(np.arange(size) * self.dx for size in self.phi.shape)


def face_material_ids_3d(centroids, geometry):
    """Assign each interface triangle by probing locally into its positive-phi solid.

    A global nearest-solid lookup is ambiguous at a material junction: an unetched substrate face
    inside a mask opening can be closer to the mask corner than to the next substrate grid node.
    The signed-distance gradient gives the physical solid-side normal and therefore the local owner.
    A nearest-solid search remains only as a fallback at degenerate zero-gradient CSG corners.
    """
    centroids = np.asarray(centroids, dtype=float)
    if geometry.material_levelsets is not None:
        material_ids = np.asarray(tuple(geometry.material_levelsets), dtype=int)
        coordinates = (centroids / geometry.dx).T
        values = np.vstack([
            map_coordinates(
                geometry.material_levelsets[int(material_id)], coordinates,
                order=1, mode="nearest", prefilter=False)
            for material_id in material_ids])
        return material_ids[np.argmax(values, axis=0)]
    solid = (geometry.phi > 0.0) & (geometry.material_id > 0)
    index = np.column_stack(np.where(solid))
    if index.size == 0:
        raise ValueError("geometry contains no labeled solid material")
    gradient = np.gradient(geometry.phi, geometry.dx)
    nearest_grid = np.rint(centroids / geometry.dx).astype(int)
    for axis in range(3):
        nearest_grid[:, axis] = np.clip(
            nearest_grid[:, axis], 0, geometry.phi.shape[axis] - 1)
    solid_normal = np.column_stack([
        component[tuple(nearest_grid.T)] for component in gradient])
    magnitude = np.linalg.norm(solid_normal, axis=1)
    valid_normal = magnitude > 1e-12
    solid_normal[valid_normal] /= magnitude[valid_normal, None]

    material = np.zeros(centroids.shape[0], dtype=int)
    unresolved = np.ones(centroids.shape[0], dtype=bool)
    for distance in (0.35, 0.75, 1.25):
        selected = np.where(unresolved & valid_normal)[0]
        if not selected.size:
            break
        probe = centroids[selected] + distance * geometry.dx * solid_normal[selected]
        probe_index = np.rint(probe / geometry.dx).astype(int)
        for axis in range(3):
            probe_index[:, axis] = np.clip(
                probe_index[:, axis], 0, geometry.phi.shape[axis] - 1)
        local_solid = geometry.phi[tuple(probe_index.T)] > 0.0
        local_material = geometry.material_id[tuple(probe_index.T)]
        accepted = local_solid & (local_material > 0)
        material[selected[accepted]] = local_material[accepted]
        unresolved[selected[accepted]] = False

    if np.any(unresolved):
        points = index * geometry.dx
        _, nearest = cKDTree(points).query(centroids[unresolved])
        chosen = index[np.asarray(nearest, dtype=int)]
        material[unresolved] = geometry.material_id[tuple(chosen.T)]
    return material
