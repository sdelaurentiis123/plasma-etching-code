"""Sparse material-local transfer weights between immutable triangle surfaces.

This module packages the inverse-distance/K-nearest-predecessor rule as an
explicit, auditable contract.  Feature evolution can select it explicitly as
the indexed reference backend; the historical inline implementation remains
the default until overlap transfer passes the moving-surface promotion gates.
It does not implement triangle-overlap remapping or signed charge transfer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_mesh_3d import TriangleSurface3D


_TRANSFER_SCHEMA = b"petch-surface-transfer-weights-3d-v1"


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


def _digest_text(digest, name, value):
    encoded = str(value).encode("utf-8")
    _digest_array(digest, name, np.frombuffer(encoded, dtype=np.uint8), "u1")


def _normalize_nonnegative_weights(raw):
    """Normalize one sparse row without sacrificing a tiny tail to closure.

    Closing a probability row by replacing its final entry with
    ``1 - sum(previous)`` is unsafe when that final entry is tiny: ordinary
    binary64 normalization can make the preceding sum one ulp larger than
    unity. The resulting negative coefficient is numerical bookkeeping, but
    accepting it would violate the interpolation contract. Instead, put the
    closure residual into the largest coefficient. Its pre-closure value is at
    least ``1 / len(row)``, so roundoff cannot change its sign.
    """
    weight = np.array(raw, dtype=float, copy=True, order="C")
    if weight.ndim != 1 or len(weight) == 0:
        raise ValueError(
            "surface-transfer weight row must be nonempty and one dimensional")
    if np.any(~np.isfinite(weight)) or np.any(weight < 0.0):
        raise ValueError(
            "raw surface-transfer weights must be finite and nonnegative")
    total = float(np.sum(weight))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError(
            "raw surface-transfer weight sum must be positive and finite")

    weight /= total
    closure_index = int(np.argmax(weight))
    weight[closure_index] = 0.0
    complement = float(np.sum(weight))
    weight[closure_index] = 1.0 - complement
    if weight[closure_index] < 0.0:
        raise RuntimeError("nonnegative surface-transfer row closure failed")
    if not np.allclose(
            float(np.sum(weight)), 1.0, rtol=0.0, atol=8e-16):
        raise RuntimeError("surface-transfer row closure did not sum to unity")
    return weight


def _conserve_nonnegative_density(raw, target_integral, new_area, *, upper_bound):
    """Match the legacy nonnegative per-material area-integral correction."""
    raw = np.maximum(np.asarray(raw, dtype=float), 0.0)
    area = np.asarray(new_area, dtype=float)
    target = float(target_integral)
    scale = max(abs(target), 1.0)
    if target < -1e-13 * scale:
        raise ValueError("negative extensive-transfer target")
    if target <= 1e-15 * scale:
        return np.zeros_like(raw)
    if upper_bound is None:
        raw_integral = float(np.dot(raw, area))
        if raw_integral <= 0.0:
            raw = np.ones_like(raw)
            raw_integral = float(area.sum())
        return raw * (target / raw_integral)

    upper = float(upper_bound)
    if not np.isfinite(upper) or upper < 0.0:
        raise ValueError("upper_bound must be finite and nonnegative")
    capacity = upper * float(area.sum())
    if target > capacity * (1.0 + 5e-13):
        raise ValueError("surface contraction exceeds bounded extensive capacity")
    seed = raw if np.any(raw > 0.0) else np.ones_like(raw)

    def integral(multiplier):
        return float(np.dot(np.minimum(multiplier * seed, upper), area))

    lower = 0.0
    upper_multiplier = 1.0
    while integral(upper_multiplier) < target:
        upper_multiplier *= 2.0
        if upper_multiplier > 1e300:
            raise RuntimeError("bounded extensive transfer failed to bracket target")
    for _ in range(80):
        midpoint = 0.5 * (lower + upper_multiplier)
        if integral(midpoint) < target:
            lower = midpoint
        else:
            upper_multiplier = midpoint
    return np.minimum(upper_multiplier * seed, upper)


@dataclass(frozen=True, eq=False)
class SurfaceTransferApplication3D:
    """Immutable values plus the semantics and material-integral audit."""

    values: np.ndarray
    semantics: str
    transfer_fingerprint: str
    material_integrals: Mapping[int, Mapping[str, float]]
    maximum_relative_integral_error: float
    metadata: Mapping[str, object]

    def __post_init__(self):
        values = _readonly_copy(self.values, "<f8")
        if values.ndim != 1 or np.any(~np.isfinite(values)):
            raise ValueError("transferred values must be one finite dimension")
        if self.semantics not in {"intensive", "extensive"}:
            raise ValueError("unknown surface-transfer semantics")
        error = float(self.maximum_relative_integral_error)
        if not np.isfinite(error) or error < 0.0:
            raise ValueError("invalid surface-transfer integral error")
        material = MappingProxyType({
            int(key): MappingProxyType({str(name): float(value)
                                        for name, value in item.items()})
            for key, item in self.material_integrals.items()
        })
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "material_integrals", material)
        object.__setattr__(self, "maximum_relative_integral_error", error)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, eq=False)
class SurfaceTransferWeights3D:
    """Immutable CSR-like old-face to new-face predecessor weights."""

    old_surface: TriangleSurface3D
    new_surface: TriangleSurface3D
    row_offsets: np.ndarray
    old_face_index: np.ndarray
    weight: np.ndarray
    source_periodic_shift: np.ndarray
    centroid_distance: np.ndarray
    row_sum: np.ndarray
    neighbor_count: int
    maximum_allowed_distance: float
    maximum_exact_surface_distance: float
    maximum_nearest_centroid_distance: float
    fingerprint: str
    application_metadata: Mapping[str, Mapping[str, object]] = field(repr=False)

    def __post_init__(self):
        def invalid(reason):
            raise ValueError(
                f"invalid sparse surface-transfer weights: {reason}")

        if (not isinstance(self.old_surface, TriangleSurface3D)
                or not isinstance(self.new_surface, TriangleSurface3D)):
            raise TypeError("surface transfer requires TriangleSurface3D inputs")
        offsets = _readonly_copy(self.row_offsets, "<i8")
        source = _readonly_copy(self.old_face_index, "<i8")
        weight = _readonly_copy(self.weight, "<f8")
        shift = _readonly_copy(self.source_periodic_shift, "<f8")
        distance = _readonly_copy(self.centroid_distance, "<f8")
        row_sum = _readonly_copy(self.row_sum, "<f8")
        row_count = len(self.new_surface.faces)
        scalar_diagnostics = np.asarray([
            self.maximum_allowed_distance,
            self.maximum_exact_surface_distance,
            self.maximum_nearest_centroid_distance,
        ], dtype=float)
        if offsets.shape != (row_count + 1,):
            invalid(
                f"row_offsets_shape={offsets.shape}, expected="
                f"{(row_count + 1,)}")
        if offsets[0] != 0:
            invalid(f"row_offsets_first={int(offsets[0])}, expected=0")
        offset_step = np.diff(offsets)
        if np.any(offset_step <= 0):
            bad = np.flatnonzero(offset_step <= 0)
            invalid(
                f"nonpositive_row_lengths count={len(bad)}, "
                f"first_row={int(bad[0])}, "
                f"first_length={int(offset_step[bad[0]])}")
        if offsets[-1] != len(source):
            invalid(
                f"row_offsets_last={int(offsets[-1])}, "
                f"source_count={len(source)}")
        if weight.shape != source.shape:
            invalid(
                f"weight_shape={weight.shape}, source_shape={source.shape}")
        if distance.shape != source.shape:
            invalid(
                f"distance_shape={distance.shape}, source_shape={source.shape}")
        if shift.shape != (len(source), 3):
            invalid(
                f"periodic_shift_shape={shift.shape}, expected="
                f"{(len(source), 3)}")
        if row_sum.shape != (row_count,):
            invalid(
                f"row_sum_shape={row_sum.shape}, expected={(row_count,)}")
        invalid_source = (source < 0) | (source >= len(self.old_surface.faces))
        if np.any(invalid_source):
            bad = np.flatnonzero(invalid_source)
            invalid(
                f"source_index_out_of_range count={len(bad)}, "
                f"first_entry={int(bad[0])}, "
                f"first_index={int(source[bad[0]])}, "
                f"old_face_count={len(self.old_surface.faces)}")
        if np.any(~np.isfinite(weight)):
            bad = np.flatnonzero(~np.isfinite(weight))
            invalid(
                f"nonfinite_weights count={len(bad)}, "
                f"first_entry={int(bad[0])}, "
                f"first_value={weight[bad[0]]!r}")
        if np.any(weight < 0.0):
            bad = np.flatnonzero(weight < 0.0)
            invalid(
                f"negative_weights count={len(bad)}, "
                f"first_entry={int(bad[0])}, "
                f"minimum={float(np.min(weight)):.17g}")
        if np.any(~np.isfinite(distance)):
            bad = np.flatnonzero(~np.isfinite(distance))
            invalid(
                f"nonfinite_centroid_distances count={len(bad)}, "
                f"first_entry={int(bad[0])}, "
                f"first_value={distance[bad[0]]!r}")
        if np.any(distance < 0.0):
            bad = np.flatnonzero(distance < 0.0)
            invalid(
                f"negative_centroid_distances count={len(bad)}, "
                f"first_entry={int(bad[0])}, "
                f"minimum={float(np.min(distance)):.17g}")
        if np.any(~np.isfinite(shift)):
            bad = np.argwhere(~np.isfinite(shift))[0]
            invalid(
                f"nonfinite_periodic_shifts first_entry="
                f"({int(bad[0])}, {int(bad[1])}), "
                f"first_value={shift[tuple(bad)]!r}")
        row_sum_error = np.abs(row_sum - 1.0)
        if not np.allclose(row_sum, 1.0, rtol=0.0, atol=8e-16):
            bad = np.flatnonzero(row_sum_error > 8e-16)
            first = int(bad[0])
            invalid(
                f"row_sum_not_unit count={len(bad)}, first_row={first}, "
                f"first_sum={float(row_sum[first]):.17g}, "
                f"maximum_absolute_error={float(np.max(row_sum_error)):.17g}, "
                "absolute_tolerance=8e-16")
        if np.any(~np.isfinite(scalar_diagnostics)):
            invalid(
                "nonfinite_scalar_diagnostics values="
                f"{scalar_diagnostics.tolist()}")
        if self.maximum_allowed_distance <= 0.0:
            invalid(
                "maximum_allowed_distance_not_positive="
                f"{float(self.maximum_allowed_distance):.17g}")
        if self.maximum_exact_surface_distance < 0.0:
            invalid(
                "maximum_exact_surface_distance_negative="
                f"{float(self.maximum_exact_surface_distance):.17g}")
        if self.maximum_nearest_centroid_distance < 0.0:
            invalid(
                "maximum_nearest_centroid_distance_negative="
                f"{float(self.maximum_nearest_centroid_distance):.17g}")
        if isinstance(self.neighbor_count, (bool, np.bool_)):
            invalid(f"neighbor_count_is_boolean={self.neighbor_count!r}")
        if int(self.neighbor_count) != self.neighbor_count:
            invalid(f"neighbor_count_not_integral={self.neighbor_count!r}")
        if int(self.neighbor_count) <= 0:
            invalid(f"neighbor_count_not_positive={self.neighbor_count!r}")
        if not isinstance(self.fingerprint, str):
            invalid(f"fingerprint_not_string={type(self.fingerprint).__name__}")
        if len(self.fingerprint) != 64:
            invalid(f"fingerprint_length={len(self.fingerprint)}, expected=64")
        metadata_keys = set(self.application_metadata)
        if metadata_keys != {"intensive", "extensive"}:
            invalid(
                f"application_metadata_keys={sorted(metadata_keys)!r}, "
                "expected=['extensive', 'intensive']")
        for row in range(row_count):
            start, stop = offsets[row:row + 2]
            selected = source[start:stop]
            if len(np.unique(selected)) != len(selected):
                raise ValueError("a transfer row contains duplicate physical old faces")
            if np.any(
                    self.old_surface.face_material_id[selected]
                    != self.new_surface.face_material_id[row]):
                raise ValueError("surface transfer crosses a material boundary")
            computed_sum = float(np.sum(weight[start:stop]))
            if (abs(computed_sum - row_sum[row]) > 8e-16
                    or abs(computed_sum - 1.0) > 8e-16):
                raise ValueError("surface-transfer row weights do not close")
        metadata = MappingProxyType({
            str(mode): MappingProxyType(dict(item))
            for mode, item in self.application_metadata.items()
        })
        object.__setattr__(self, "row_offsets", offsets)
        object.__setattr__(self, "old_face_index", source)
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "source_periodic_shift", shift)
        object.__setattr__(self, "centroid_distance", distance)
        object.__setattr__(self, "row_sum", row_sum)
        object.__setattr__(self, "application_metadata", metadata)

    def _raw_apply(self, old_values):
        values = np.asarray(old_values, dtype=float)
        if (values.shape != (len(self.old_surface.faces),)
                or np.any(~np.isfinite(values))):
            raise ValueError("old_values must be one finite value per old face")
        output = np.empty(len(self.new_surface.faces), dtype=float)
        for row in range(len(output)):
            start, stop = self.row_offsets[row:row + 2]
            output[row] = float(np.dot(
                self.weight[start:stop], values[self.old_face_index[start:stop]]))
        return values, output

    def apply_intensive(self, old_values, *, lower_bound=None, upper_bound=None):
        """Interpolate an intensive field without imposing area conservation."""
        _, output = self._raw_apply(old_values)
        if lower_bound is not None:
            lower = float(lower_bound)
            if not np.isfinite(lower):
                raise ValueError("lower_bound must be finite")
            output = np.maximum(output, lower)
        if upper_bound is not None:
            upper = float(upper_bound)
            if not np.isfinite(upper):
                raise ValueError("upper_bound must be finite")
            if lower_bound is not None and upper < float(lower_bound):
                raise ValueError("upper_bound must not be below lower_bound")
            output = np.minimum(output, upper)
        integrals = self._material_integrals(old_values, output)
        return SurfaceTransferApplication3D(
            output, "intensive", self.fingerprint, integrals,
            self._maximum_integral_error(integrals),
            {
                **dict(self.application_metadata["intensive"]),
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
            })

    def apply_extensive(self, old_density, *, upper_bound=None):
        """Transfer a nonnegative density and preserve each material area integral."""
        values, raw = self._raw_apply(old_density)
        if np.any(values < 0.0):
            raise ValueError("extensive density must be nonnegative")
        output = np.empty_like(raw)
        for material in sorted(set(self.old_surface.face_material_id.tolist())):
            old_selected = self.old_surface.face_material_id == material
            new_selected = self.new_surface.face_material_id == material
            target = float(np.dot(
                values[old_selected], self.old_surface.face_area[old_selected]))
            output[new_selected] = _conserve_nonnegative_density(
                raw[new_selected], target,
                self.new_surface.face_area[new_selected],
                upper_bound=upper_bound)
        integrals = self._material_integrals(values, output)
        return SurfaceTransferApplication3D(
            output, "extensive", self.fingerprint, integrals,
            self._maximum_integral_error(integrals),
            {
                **dict(self.application_metadata["extensive"]),
                "upper_bound": upper_bound,
            })

    def _material_integrals(self, old_values, new_values):
        old_values = np.asarray(old_values, dtype=float)
        new_values = np.asarray(new_values, dtype=float)
        output = {}
        for material in sorted(set(self.old_surface.face_material_id.tolist())):
            old_selected = self.old_surface.face_material_id == material
            new_selected = self.new_surface.face_material_id == material
            before = float(np.dot(
                old_values[old_selected], self.old_surface.face_area[old_selected]))
            after = float(np.dot(
                new_values[new_selected], self.new_surface.face_area[new_selected]))
            output[int(material)] = {
                "old_area_integral": before,
                "new_area_integral": after,
                "relative_difference": abs(after - before) / max(abs(before), 1.0),
            }
        return output

    @staticmethod
    def _maximum_integral_error(material_integrals):
        return max(
            (float(item["relative_difference"])
             for item in material_integrals.values()), default=0.0)


def build_surface_transfer_3d(
        old_surface, new_surface, *, neighbor_count=4, maximum_distance):
    """Build deterministic sparse material-local predecessor weights.

    ``maximum_distance`` gates exact new-centroid-to-old-triangle distance,
    while K-nearest predecessor weights retain the legacy area-weighted inverse
    squared centroid-distance rule.
    """
    if (not isinstance(old_surface, TriangleSurface3D)
            or not isinstance(new_surface, TriangleSurface3D)):
        raise TypeError("surface transfer requires TriangleSurface3D inputs")
    if (old_surface.periodic_lengths != new_surface.periodic_lengths
            or old_surface.periodic_origin != new_surface.periodic_origin):
        raise ValueError("old and new surfaces require identical periodic cells")
    if (isinstance(neighbor_count, (bool, np.bool_))
            or int(neighbor_count) != neighbor_count or int(neighbor_count) <= 0):
        raise ValueError("neighbor_count must be a positive integer")
    neighbor_count = int(neighbor_count)
    maximum_distance = float(maximum_distance)
    if not np.isfinite(maximum_distance) or maximum_distance <= 0.0:
        raise ValueError("maximum_distance must be positive and finite")
    old_materials = set(old_surface.face_material_id.tolist())
    new_materials = set(new_surface.face_material_id.tolist())
    if old_materials != new_materials:
        raise ValueError(
            "material surface appeared or disappeared; initialize/retire state explicitly")

    # An unchanged mesh is a true Lagrangian identity, including the unusual
    # case of coincident centroids.  Preserve that identity explicitly.
    if old_surface.fingerprint == new_surface.fingerprint:
        rows = len(new_surface.faces)
        offsets = np.arange(rows + 1, dtype=int)
        source = np.arange(rows, dtype=int)
        weight = np.ones(rows, dtype=float)
        shift = np.zeros((rows, 3), dtype=float)
        centroid_distance = np.zeros(rows, dtype=float)
        row_sum = np.ones(rows, dtype=float)
        maximum_exact = 0.0
        maximum_centroid = 0.0
    else:
        row_source = [None] * len(new_surface.faces)
        row_weight = [None] * len(new_surface.faces)
        row_shift = [None] * len(new_surface.faces)
        row_distance = [None] * len(new_surface.faces)
        maximum_exact = 0.0
        maximum_centroid = 0.0
        for material in sorted(old_materials):
            old_index = np.flatnonzero(old_surface.face_material_id == material)
            new_index = np.flatnonzero(new_surface.face_material_id == material)
            exact = old_surface.nearest(
                new_surface.face_centroid[new_index], material_id=int(material),
                maximum_distance=maximum_distance)
            if not np.all(exact.found):
                # The fast bounded query intentionally reports absence, not a
                # fabricated distance.  Pay for an unbounded exact diagnostic
                # only on this refusal path.
                rejected = old_surface.nearest(
                    new_surface.face_centroid[new_index],
                    material_id=int(material))
                material_maximum = float(np.max(rejected.distance))
                union = old_surface.nearest(
                    new_surface.face_centroid[new_index])
                raise ValueError(
                    f"material-local surface transfer distance {material_maximum:g} "
                    f"exceeds {maximum_distance:g} for material {int(material)}; "
                    f"union-surface distance={float(np.max(union.distance)):g}; "
                    "borrowing across materials is forbidden")
            material_maximum = float(np.max(exact.distance))
            maximum_exact = max(maximum_exact, material_maximum)

            count = min(neighbor_count, len(old_index))
            centroid = old_surface.nearest_face_centroids(
                new_surface.face_centroid[new_index], count=count,
                material_id=int(material))
            for local_row, new_face in enumerate(new_index):
                point = new_surface.face_centroid[new_face]
                selected_source = centroid.face_index[local_row]
                selected_distance = centroid.distance[local_row]
                selected_shift = centroid.periodic_shift[local_row]
                maximum_centroid = max(
                    maximum_centroid, float(selected_distance[0]))
                scale = max(
                    maximum_distance,
                    float(np.max(np.abs(
                        old_surface.face_centroid[old_index]), initial=0.0)),
                    float(np.max(np.abs(point), initial=0.0)), 1.0)
                distance_floor = 64.0 * np.finfo(float).eps * scale
                if selected_distance[0] <= distance_floor:
                    selected_source = selected_source[:1]
                    selected_distance = selected_distance[:1]
                    selected_shift = selected_shift[:1]
                    selected_weight = np.ones(1, dtype=float)
                else:
                    raw_weight = (
                        old_surface.face_area[selected_source]
                        / np.maximum(selected_distance ** 2, distance_floor ** 2))
                    # Close through the largest positive coefficient. The
                    # former final-entry closure could turn an O(1e-27) tail
                    # into -2.22e-16 when preceding weights rounded above one.
                    selected_weight = _normalize_nonnegative_weights(raw_weight)
                row_source[new_face] = selected_source
                row_weight[new_face] = selected_weight
                row_shift[new_face] = selected_shift
                row_distance[new_face] = selected_distance

        offsets = np.zeros(len(new_surface.faces) + 1, dtype=int)
        for row, selected in enumerate(row_source):
            if selected is None or len(selected) == 0:
                raise RuntimeError("surface-transfer row was not constructed")
            offsets[row + 1] = offsets[row] + len(selected)
        source = np.concatenate(row_source)
        weight = np.concatenate(row_weight)
        shift = np.concatenate(row_shift, axis=0)
        centroid_distance = np.concatenate(row_distance)
        row_sum = np.asarray([
            float(np.sum(weight[offsets[row]:offsets[row + 1]]))
            for row in range(len(new_surface.faces))
        ])

    digest = sha256()
    digest.update(_TRANSFER_SCHEMA)
    _digest_text(digest, "old_surface_fingerprint", old_surface.fingerprint)
    _digest_text(digest, "new_surface_fingerprint", new_surface.fingerprint)
    _digest_array(digest, "row_offsets", offsets, "<i8")
    _digest_array(digest, "old_face_index", source, "<i8")
    _digest_array(digest, "weight", weight, "<f8")
    _digest_array(digest, "source_periodic_shift", shift, "<f8")
    _digest_array(digest, "centroid_distance", centroid_distance, "<f8")
    _digest_array(
        digest, "settings", [neighbor_count, maximum_distance], "<f8")
    application_metadata = {
        "intensive": {
            "input_representation": "value per face",
            "row_operation": "material-local area-weighted inverse-distance interpolation",
            "area_integral_preserved": False,
        },
        "extensive": {
            "input_representation": "nonnegative density per surface area",
            "row_operation": "material-local area-weighted inverse-distance interpolation",
            "conservation_correction": "per-material new-area integral rescaling",
            "area_integral_preserved": True,
        },
    }
    return SurfaceTransferWeights3D(
        old_surface=old_surface,
        new_surface=new_surface,
        row_offsets=offsets,
        old_face_index=source,
        weight=weight,
        source_periodic_shift=shift,
        centroid_distance=centroid_distance,
        row_sum=row_sum,
        neighbor_count=neighbor_count,
        maximum_allowed_distance=maximum_distance,
        maximum_exact_surface_distance=maximum_exact,
        maximum_nearest_centroid_distance=maximum_centroid,
        fingerprint=digest.hexdigest(),
        application_metadata=application_metadata,
    )
