"""Bridge an exactly extruded triangle surface to deterministic 2-D exchange.

Marching cubes represents an extruded line element as rectangular strips split into triangles.
This module groups those triangles by their physical cross-section segment, constructs the
deterministic crossed-string operator on the unique segments, and expands the result back to the
existing face-level :class:`DiffuseFormFactors3D` contract.

The bridge is deliberately strict. It refuses non-extruded triangles, incomplete strips,
inconsistent gas normals, and face fields that vary along the declared symmetry direction. A
genuinely three-dimensional surface must use a three-dimensional transport backend.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256

import numpy as np

from .deterministic_exchange_2d import (
    DeterministicLineExchange2D, build_deterministic_line_exchange_2d,
)
from .neutral_radiosity_3d import DiffuseFormFactors3D


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


def _canonical_segment_key(segment, tolerance):
    quantized = np.rint(np.asarray(segment, dtype=float) / tolerance).astype(np.int64)
    first, second = map(tuple, quantized)
    if second < first:
        first, second = second, first
    return first + second


def _triangle_section_segment(triangle, section_axes, tolerance):
    projected = np.asarray(triangle[:, section_axes], dtype=float)
    pair = ((0, 1), (0, 2), (1, 2))
    distance = np.asarray([
        np.linalg.norm(projected[left] - projected[right]) for left, right in pair])
    left, right = pair[int(np.argmax(distance))]
    if distance.max(initial=0.0) <= tolerance:
        raise ValueError("extruded triangle collapses to a point in the section")
    segment = np.stack((projected[left], projected[right]))
    remaining = next(index for index in range(3) if index not in (left, right))
    endpoint_distance = min(
        np.linalg.norm(projected[remaining] - segment[0]),
        np.linalg.norm(projected[remaining] - segment[1]))
    if endpoint_distance > tolerance:
        raise ValueError(
            "triangle is not an extrusion of one cross-section line element")
    if tuple(segment[1]) < tuple(segment[0]):
        segment = segment[::-1]
    return segment


@dataclass(frozen=True)
class ExtrudedTriangleExchange3D:
    """Face-level deterministic exchange plus its symmetry-reduction receipt."""

    extrusion_axis: int
    extrusion_length: float
    face_area: np.ndarray
    face_group_index: np.ndarray
    group_area: np.ndarray
    group_segments_2d: np.ndarray
    group_gas_normals_2d: np.ndarray
    line_exchange: DeterministicLineExchange2D
    form_factors: DiffuseFormFactors3D
    maximum_group_area_relative_error: float
    maximum_area_reciprocity_error: float
    fingerprint: str

    def __post_init__(self):
        face_area = _readonly(self.face_area)
        face_group = _readonly(self.face_group_index, int)
        group_area = _readonly(self.group_area)
        group_segment = _readonly(self.group_segments_2d)
        group_normal = _readonly(self.group_gas_normals_2d)
        face_count = len(face_area)
        group_count = len(group_area)
        if (self.extrusion_axis not in (0, 1, 2)
                or not np.isfinite(self.extrusion_length) or self.extrusion_length <= 0.0
                or face_area.shape != (face_count,) or np.any(face_area <= 0.0)
                or face_group.shape != (face_count,)
                or np.any(face_group < 0) or np.any(face_group >= group_count)
                or group_segment.shape != (group_count, 2, 2)
                or group_normal.shape != (group_count, 2)
                or self.form_factors.face_count != face_count
                or self.line_exchange.segments.shape != group_segment.shape):
            raise ValueError("invalid extruded triangle-exchange result")
        for name, value in (
                ("face_area", face_area), ("face_group_index", face_group),
                ("group_area", group_area), ("group_segments_2d", group_segment),
                ("group_gas_normals_2d", group_normal)):
            object.__setattr__(self, name, value)

    @property
    def face_count(self):
        return len(self.face_area)

    @property
    def group_count(self):
        return len(self.group_area)

    def area_weighted_group_mean(self, face_field):
        value = np.asarray(face_field, dtype=float)
        if value.shape != (self.face_count,):
            raise ValueError("extruded face field must contain one scalar per triangle")
        weighted = np.bincount(
            self.face_group_index, weights=self.face_area * value,
            minlength=self.group_count)
        return weighted / self.group_area

    def certify_face_field(self, face_field, *, relative_tolerance, absolute_tolerance):
        """Return group means or refuse variation along the declared extrusion axis."""
        value = np.asarray(face_field, dtype=float)
        mean = self.area_weighted_group_mean(value)
        difference = np.abs(value - mean[self.face_group_index])
        tolerance = (
            float(absolute_tolerance)
            + float(relative_tolerance) * np.abs(mean[self.face_group_index]))
        if (not np.all(np.isfinite(value)) or relative_tolerance < 0.0
                or absolute_tolerance < 0.0):
            raise ValueError("invalid extrusion-invariance certification inputs")
        if np.any(difference > tolerance):
            worst = int(np.argmax(difference - tolerance))
            raise ValueError(
                "face field violates the declared extrusion invariance; "
                f"face={worst}, group={int(self.face_group_index[worst])}, "
                f"difference={difference[worst]:.6g}, tolerance={tolerance[worst]:.6g}")
        return mean


def build_extruded_triangle_exchange_3d(
        vertices, faces, gas_normals, *, extrusion_axis=1, extrusion_length=None,
        geometry_tolerance=None, normal_tolerance=2.0e-6,
        area_relative_tolerance=2.0e-6, exchange_relative_tolerance=1.0e-5,
        exchange_absolute_tolerance=1.0e-12,
        minimum_refinement_level=2, maximum_refinement_level=18):
    """Construct deterministic face-level exchange for an exactly extruded triangle mesh."""
    vertex = np.asarray(vertices, dtype=float)
    face = np.asarray(faces, dtype=int)
    normal = np.asarray(gas_normals, dtype=float)
    axis = int(extrusion_axis)
    if (vertex.ndim != 2 or vertex.shape[1] != 3
            or face.ndim != 2 or face.shape[1] != 3 or len(face) == 0
            or normal.shape != (len(face), 3)
            or np.any(face < 0) or np.any(face >= len(vertex))
            or axis not in (0, 1, 2)
            or np.any(~np.isfinite(vertex)) or np.any(~np.isfinite(normal))):
        raise ValueError("invalid extruded triangle mesh")
    span = np.ptp(vertex, axis=0)
    scale = max(float(np.max(span)), 1.0)
    tolerance = (
        128.0 * np.finfo(np.float32).eps * scale
        if geometry_tolerance is None else float(geometry_tolerance))
    if tolerance <= 0.0 or normal_tolerance <= 0.0 or area_relative_tolerance <= 0.0:
        raise ValueError("invalid extrusion certification tolerances")
    length = float(span[axis] if extrusion_length is None else extrusion_length)
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("extrusion length must be finite and positive")
    section_axes = tuple(index for index in range(3) if index != axis)
    triangle = vertex[face]
    cross = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    face_area = 0.5 * np.linalg.norm(cross, axis=1)
    normal_length = np.linalg.norm(normal, axis=1)
    if np.any(face_area <= 0.0) or np.any(normal_length <= 0.0):
        raise ValueError("extruded mesh contains a zero-measure face or normal")
    unit_normal = normal / normal_length[:, None]
    if np.any(np.abs(unit_normal[:, axis]) > normal_tolerance):
        raise ValueError("surface normal has a component along the declared extrusion axis")

    segment_by_face = np.stack([
        _triangle_section_segment(item, section_axes, tolerance)
        for item in triangle])
    key_by_face = [
        _canonical_segment_key(item, tolerance) for item in segment_by_face]
    keys = sorted(set(key_by_face))
    key_to_group = {key: index for index, key in enumerate(keys)}
    face_group = np.asarray([key_to_group[key] for key in key_by_face], dtype=int)
    group_count = len(keys)
    group_segment = np.empty((group_count, 2, 2), dtype=float)
    group_normal = np.empty((group_count, 2), dtype=float)
    group_area = np.bincount(
        face_group, weights=face_area, minlength=group_count).astype(float)
    maximum_area_error = 0.0
    for group in range(group_count):
        member = np.flatnonzero(face_group == group)
        representative = segment_by_face[member[0]]
        group_segment[group] = representative
        projected_normal = unit_normal[member][:, section_axes]
        weighted_normal = np.sum(projected_normal * face_area[member, None], axis=0)
        weighted_length = float(np.linalg.norm(weighted_normal))
        if weighted_length <= 0.0:
            raise ValueError("extruded group has cancelling gas normals")
        group_normal[group] = weighted_normal / weighted_length
        alignment = projected_normal @ group_normal[group]
        if np.any(alignment < 1.0 - normal_tolerance):
            raise ValueError("extruded group contains inconsistent gas normals")
        expected_area = float(
            np.linalg.norm(representative[1] - representative[0]) * length)
        relative_error = abs(group_area[group] - expected_area) / expected_area
        maximum_area_error = max(maximum_area_error, relative_error)
        if relative_error > area_relative_tolerance:
            raise ValueError(
                "triangle strip does not cover exactly one extrusion period; "
                f"group={group}, relative_area_error={relative_error:.6g}")

    line_exchange = build_deterministic_line_exchange_2d(
        group_segment, group_normal,
        relative_tolerance=exchange_relative_tolerance,
        absolute_tolerance=exchange_absolute_tolerance,
        minimum_refinement_level=minimum_refinement_level,
        maximum_refinement_level=maximum_refinement_level,
        geometry_tolerance=tolerance)

    source_parts = []
    target_parts = []
    fraction_parts = []
    for source_group in range(group_count):
        source_face = np.flatnonzero(face_group == source_group)
        for target_group in np.flatnonzero(
                line_exchange.transfer_fraction[source_group] > 0.0):
            target_face = np.flatnonzero(face_group == target_group)
            source_parts.append(np.repeat(source_face, len(target_face)))
            target_parts.append(np.tile(target_face, len(source_face)))
            target_weight = face_area[target_face] / group_area[target_group]
            fraction_parts.append(np.tile(
                line_exchange.transfer_fraction[source_group, target_group]
                * target_weight, len(source_face)))
    source_index = (
        np.concatenate(source_parts) if source_parts else np.empty(0, dtype=int))
    target_index = (
        np.concatenate(target_parts) if target_parts else np.empty(0, dtype=int))
    transfer_fraction = (
        np.concatenate(fraction_parts) if fraction_parts else np.empty(0, dtype=float))
    escape_fraction = line_exchange.escape_fraction[face_group]
    factors = DiffuseFormFactors3D(
        len(face), source_index, target_index, transfer_fraction,
        escape_fraction, rays_per_face=1)

    reciprocity = np.abs(
        face_area[source_index] * transfer_fraction
        - face_area[target_index] * np.asarray([
            line_exchange.transfer_fraction[face_group[target], face_group[source]]
            * face_area[source] / group_area[face_group[source]]
            for source, target in zip(source_index, target_index)]))
    maximum_reciprocity = float(np.max(reciprocity, initial=0.0))

    digest = sha256()
    digest.update(b"petch.extruded-triangle-exchange-3d.v1\0")
    for name, value, dtype in (
            ("vertices", vertex, "<f8"), ("faces", face, "<i8"),
            ("gas_normals", unit_normal, "<f8"),
            ("face_group", face_group, "<i8")):
        array = np.ascontiguousarray(value, dtype=dtype)
        digest.update(name.encode("ascii") + b"\0")
        digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
        digest.update(array.tobytes())
    digest.update(line_exchange.fingerprint.encode("ascii"))
    digest.update(np.asarray([axis, length, tolerance], dtype="<f8").tobytes())
    return ExtrudedTriangleExchange3D(
        extrusion_axis=axis, extrusion_length=length,
        face_area=face_area, face_group_index=face_group,
        group_area=group_area, group_segments_2d=group_segment,
        group_gas_normals_2d=group_normal, line_exchange=line_exchange,
        form_factors=factors,
        maximum_group_area_relative_error=maximum_area_error,
        maximum_area_reciprocity_error=maximum_reciprocity,
        fingerprint=digest.hexdigest())
