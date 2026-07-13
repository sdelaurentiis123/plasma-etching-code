"""Conservative signed surface-charge transfer across one material-interface motion.

Declared v1 closure
-------------------
Charge resides on material surface elements.  ``material_normal_displacement > 0`` means the
material surface advances into vacuum/deposition-like growth, so charge rides with the retained
material surface.  ``material_normal_displacement < 0`` means the charged surface layer is etched
away; its charge is removed with that material and itemized in the ledger.  Newly exposed surface is
uncharged.  Charge-left-behind and partial-retention closures are intentionally not implemented.

The operator is geometric and solver-independent.  It maps a piecewise-constant signed sheet charge
between old and new triangle meshes, preserves retained positive and negative charge separately, and
never uses net-charge cancellation as its conservation scale.  The compatible Q1 projection remains
the separate authoritative field-coupling operator.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.spatial import cKDTree


def _triangle_geometry(vertices, faces, length_unit_m):
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    if (vertices.ndim != 2 or vertices.shape[1] != 3
            or faces.ndim != 2 or faces.shape[1] != 3
            or np.any(~np.isfinite(vertices)) or np.any(faces < 0)
            or np.any(faces >= len(vertices))
            or not np.isfinite(length_unit_m) or length_unit_m <= 0.0):
        raise ValueError("invalid triangulated surface")
    triangle = vertices[faces]
    cross = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    area_mesh2 = 0.5 * np.linalg.norm(cross, axis=1)
    if np.any(~np.isfinite(area_mesh2)) or np.any(area_mesh2 <= 0.0):
        raise ValueError("surface-charge triangles must have positive area")
    return triangle.mean(axis=1), area_mesh2 * float(length_unit_m) ** 2


@dataclass(frozen=True)
class SurfaceChargeRemap3DResult:
    """New sheet charge and an explicit retained/removed charge ledger."""

    sigma_c_per_m2: np.ndarray
    face_charge_c: np.ndarray
    removed_charge_by_old_face_c: np.ndarray
    retained_positive_charge_c: float
    retained_negative_charge_c: float
    removed_positive_charge_c: float
    removed_negative_charge_c: float
    relative_charge_balance_error: float
    material_ledger: Mapping[int, Mapping[str, float]]
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        face_charge = np.asarray(self.face_charge_c, dtype=float).copy()
        removed = np.asarray(self.removed_charge_by_old_face_c, dtype=float).copy()
        if (sigma.ndim != 1 or face_charge.shape != sigma.shape or removed.ndim != 1
                or np.any(~np.isfinite(sigma)) or np.any(~np.isfinite(face_charge))
                or np.any(~np.isfinite(removed))):
            raise ValueError("invalid surface-charge remap arrays")
        values = np.asarray([
            self.retained_positive_charge_c, self.retained_negative_charge_c,
            self.removed_positive_charge_c, self.removed_negative_charge_c,
            self.relative_charge_balance_error], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("surface-charge ledger magnitudes must be finite and nonnegative")
        for value in (sigma, face_charge, removed):
            value.setflags(write=False)
        ledger = {
            int(material): MappingProxyType(dict(item))
            for material, item in self.material_ledger.items()}
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "face_charge_c", face_charge)
        object.__setattr__(self, "removed_charge_by_old_face_c", removed)
        object.__setattr__(self, "material_ledger", MappingProxyType(ledger))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))

    @property
    def retained_net_charge_c(self):
        return self.retained_positive_charge_c - self.retained_negative_charge_c

    @property
    def removed_net_charge_c(self):
        return self.removed_positive_charge_c - self.removed_negative_charge_c


def _conserve_nonnegative_density(raw, target_charge_c, area_m2):
    raw = np.maximum(np.asarray(raw, dtype=float), 0.0)
    area = np.asarray(area_m2, dtype=float)
    if target_charge_c == 0.0:
        return np.zeros_like(raw)
    achieved = float(np.dot(raw, area))
    if achieved <= 0.0:
        return np.full(raw.shape, target_charge_c / float(np.sum(area)))
    return raw * (target_charge_c / achieved)


def remap_surface_charge_3d(
        old_vertices, old_faces, old_sigma_c_per_m2, old_material_id,
        material_normal_displacement, new_vertices, new_faces, new_material_id, *,
        mesh_length_unit_m=1.0, neighbor_count=4, maximum_distance=None,
        displacement_tolerance=0.0):
    """Transfer retained signed charge and ledger charge removed with etched material.

    Distances and ``material_normal_displacement`` use mesh coordinate units.  A displacement within
    ``displacement_tolerance`` is retained; negative motion beyond that tolerance is removed.  The
    no-motion, identical-mesh path returns the original sheet-charge array bitwise.
    """
    old_centroid, old_area_m2 = _triangle_geometry(
        old_vertices, old_faces, mesh_length_unit_m)
    new_centroid, new_area_m2 = _triangle_geometry(
        new_vertices, new_faces, mesh_length_unit_m)
    old_sigma_input = np.asarray(old_sigma_c_per_m2)
    old_sigma = np.asarray(old_sigma_c_per_m2, dtype=float)
    old_material = np.asarray(old_material_id, dtype=int)
    new_material = np.asarray(new_material_id, dtype=int)
    displacement = np.asarray(material_normal_displacement, dtype=float)
    if (old_sigma.shape != old_area_m2.shape or old_material.shape != old_area_m2.shape
            or displacement.shape != old_area_m2.shape
            or new_material.shape != new_area_m2.shape
            or np.any(~np.isfinite(old_sigma)) or np.any(~np.isfinite(displacement))
            or np.any(old_material <= 0) or np.any(new_material <= 0)
            or int(neighbor_count) != neighbor_count or neighbor_count <= 0
            or not np.isfinite(displacement_tolerance) or displacement_tolerance < 0.0):
        raise ValueError("invalid surface-charge remap state")
    if maximum_distance is None:
        span = np.ptp(np.vstack((old_centroid, new_centroid)), axis=0)
        maximum_distance = max(float(np.linalg.norm(span)), np.finfo(float).tiny)
    if not np.isfinite(maximum_distance) or maximum_distance <= 0.0:
        raise ValueError("maximum remap distance must be finite and positive")

    identical = (
        np.array_equal(np.asarray(old_vertices), np.asarray(new_vertices))
        and np.array_equal(np.asarray(old_faces), np.asarray(new_faces))
        and np.array_equal(old_material, new_material)
        and np.all(displacement == 0.0))
    old_charge = old_sigma * old_area_m2
    removed_mask = displacement < -float(displacement_tolerance)
    removed_charge = np.where(removed_mask, old_charge, 0.0)
    retained_mask = ~removed_mask
    output_sigma = np.zeros(new_area_m2.shape)
    maximum_nearest = 0.0
    material_ledger = {}

    if identical:
        output_sigma = old_sigma_input.copy()
    elif (np.array_equal(np.asarray(old_faces), np.asarray(new_faces))
          and np.array_equal(old_material, new_material)
          and old_area_m2.shape == new_area_m2.shape
          and not np.any(removed_mask)):
        # Unchanged connectivity supplies an exact Lagrangian face correspondence.  Carry each
        # face's charge directly; density changes only with its new physical area.
        output_sigma = old_charge / new_area_m2
    else:
        for material in sorted(set(old_material.tolist()) | set(new_material.tolist())):
            source = np.where((old_material == material) & retained_mask)[0]
            target = np.where(new_material == material)[0]
            if source.size == 0:
                continue
            if target.size == 0:
                raise ValueError(
                    f"retained charged material {material} has no new surface")
            count = min(int(neighbor_count), source.size)
            distance, local = cKDTree(old_centroid[source]).query(
                new_centroid[target], k=count)
            if count == 1:
                distance = np.asarray(distance)[:, None]
                local = np.asarray(local)[:, None]
            nearest = float(np.max(distance[:, 0]))
            maximum_nearest = max(maximum_nearest, nearest)
            if nearest > float(maximum_distance):
                raise ValueError(
                    f"surface-charge remap distance {nearest:g} exceeds {maximum_distance:g}")
            source_index = source[np.asarray(local, dtype=int)]
            coordinate_scale = max(
                float(np.max(np.abs(old_centroid[source]))),
                float(np.max(np.abs(new_centroid[target]))), 1.0)
            distance_floor = 64.0 * np.finfo(float).eps * coordinate_scale
            exact = distance[:, 0] <= distance_floor
            weight = old_area_m2[source_index] / np.maximum(
                distance * distance, distance_floor ** 2)
            if np.any(exact):
                weight[exact] = 0.0
                weight[exact, 0] = 1.0
            weight /= weight.sum(axis=1, keepdims=True)
            positive_raw = np.sum(weight * np.maximum(old_sigma[source_index], 0.0), axis=1)
            negative_raw = np.sum(weight * np.maximum(-old_sigma[source_index], 0.0), axis=1)
            positive_target = float(np.sum(np.maximum(old_charge[source], 0.0)))
            negative_target = float(np.sum(np.maximum(-old_charge[source], 0.0)))
            output_sigma[target] = (
                _conserve_nonnegative_density(positive_raw, positive_target, new_area_m2[target])
                - _conserve_nonnegative_density(negative_raw, negative_target, new_area_m2[target]))

    new_charge = output_sigma * new_area_m2
    retained_positive = float(np.sum(np.maximum(old_charge[retained_mask], 0.0)))
    retained_negative = float(np.sum(np.maximum(-old_charge[retained_mask], 0.0)))
    removed_positive = float(np.sum(np.maximum(removed_charge, 0.0)))
    removed_negative = float(np.sum(np.maximum(-removed_charge, 0.0)))
    old_net = float(np.sum(old_charge))
    ledger_net = float(np.sum(new_charge) + np.sum(removed_charge))
    throughput = float(np.sum(np.abs(old_charge)))
    relative_error = abs(old_net - ledger_net) / max(throughput, np.finfo(float).tiny)

    for material in sorted(set(old_material.tolist()) | set(new_material.tolist())):
        old_selected = old_material == material
        new_selected = new_material == material
        material_removed = removed_charge[old_selected]
        material_ledger[int(material)] = dict(
            old_charge_c=float(np.sum(old_charge[old_selected])),
            retained_charge_c=float(np.sum(new_charge[new_selected])),
            removed_charge_c=float(np.sum(material_removed)),
            removed_positive_charge_c=float(np.sum(np.maximum(material_removed, 0.0))),
            removed_negative_charge_c=float(np.sum(np.maximum(-material_removed, 0.0))))

    return SurfaceChargeRemap3DResult(
        sigma_c_per_m2=output_sigma, face_charge_c=new_charge,
        removed_charge_by_old_face_c=removed_charge,
        retained_positive_charge_c=retained_positive,
        retained_negative_charge_c=retained_negative,
        removed_positive_charge_c=removed_positive,
        removed_negative_charge_c=removed_negative,
        relative_charge_balance_error=relative_error,
        material_ledger=material_ledger,
        diagnostics=dict(
            closure="charge rides advancing material; etched charged layer is removed; "
                    "newly exposed surface is uncharged",
            neighbor_count=int(neighbor_count),
            maximum_nearest_distance=float(maximum_nearest),
            maximum_allowed_distance=float(maximum_distance),
            removed_face_count=int(np.count_nonzero(removed_mask)),
            retained_face_count=int(np.count_nonzero(retained_mask)),
            old_total_charge_c=old_net,
            new_total_charge_c=float(np.sum(new_charge)),
            removed_total_charge_c=float(np.sum(removed_charge))))
