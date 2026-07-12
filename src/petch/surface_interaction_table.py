"""Versionable physical surface-interaction tables for MD, experiment, and calibrated data."""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class SurfaceInteractionDomainError(ValueError):
    """Raised when an interaction table is asked to extrapolate without explicit authorization."""


@dataclass(frozen=True)
class InteractionAxis:
    name: str
    values: np.ndarray
    unit: str
    interpolation: str = "linear"

    def __post_init__(self):
        values = np.asarray(self.values, dtype=float).copy()
        if (not self.name or not self.unit or values.ndim != 1 or values.size < 2
                or np.any(~np.isfinite(values)) or np.any(np.diff(values) <= 0.0)
                or self.interpolation not in {"linear", "log"}
                or (self.interpolation == "log" and np.any(values <= 0.0))):
            raise ValueError("interaction axes require named, unit-bearing, increasing finite values")
        values.setflags(write=False)
        object.__setattr__(self, "values", values)

    @property
    def interpolation_values(self):
        return np.log(self.values) if self.interpolation == "log" else self.values

    def transform(self, value):
        value = np.asarray(value, dtype=float)
        if np.any(~np.isfinite(value)) or (self.interpolation == "log" and np.any(value <= 0.0)):
            raise ValueError(f"invalid coordinate for axis {self.name!r}")
        return np.log(value) if self.interpolation == "log" else value


@dataclass(frozen=True)
class SurfaceInteractionEvaluation:
    values: Mapping[str, np.ndarray]
    standard_uncertainty: Mapping[str, np.ndarray]
    extrapolated_fraction: float
    outside_axes: tuple[str, ...]
    table_fingerprint: str

    def __post_init__(self):
        def freeze(mapping):
            output = {}
            for name, value in mapping.items():
                array = np.asarray(value, dtype=float).copy(); array.setflags(write=False)
                output[name] = array
            return MappingProxyType(output)
        object.__setattr__(self, "values", freeze(self.values))
        object.__setattr__(self, "standard_uncertainty", freeze(self.standard_uncertainty))
        object.__setattr__(self, "outside_axes", tuple(self.outside_axes))


@dataclass(frozen=True)
class SurfaceInteractionInterpolationAudit:
    """Interior-node leave-one-out error, kept separate from source uncertainty."""

    axis_name: str
    output_name: str
    coordinates: np.ndarray
    observed: np.ndarray
    predicted: np.ndarray
    absolute_error: np.ndarray
    table_fingerprint: str

    def __post_init__(self):
        for name in ("coordinates", "observed", "predicted", "absolute_error"):
            array = np.asarray(getattr(self, name), dtype=float).copy()
            array.setflags(write=False); object.__setattr__(self, name, array)


@dataclass(frozen=True)
class SurfaceInteractionTable:
    """Regular-grid physical response table with explicit validity and conservation contracts.

    Axes can represent energy, incidence cosine, flux ratio, coverage, temperature, or other physical
    state. Outputs can represent yields, branching fractions, damage depth, or mixed-layer observables.
    Interpolation is multilinear in each axis' declared linear or logarithmic coordinate. Extrapolation
    is refused by default and, when explicitly enabled, its affected point fraction is returned.
    """

    material: str
    incident_species: tuple[str, ...]
    axes: tuple[InteractionAxis, ...]
    outputs: Mapping[str, np.ndarray]
    output_units: Mapping[str, str]
    provenance: Mapping[str, object]
    standard_uncertainty: Mapping[str, np.ndarray | float] = field(default_factory=dict)
    bounds: Mapping[str, tuple[float | None, float | None]] = field(default_factory=dict)
    conservation_groups: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    schema_version: int = 1

    def __post_init__(self):
        axes = tuple(self.axes); species = tuple(self.incident_species)
        if (not self.material or not species or any(not item for item in species)
                or not axes or len({axis.name for axis in axes}) != len(axes)
                or int(self.schema_version) != self.schema_version or self.schema_version <= 0):
            raise ValueError("interaction table requires material, species, unique axes, and schema")
        shape = tuple(axis.values.size for axis in axes)
        outputs = {}
        for name, value in self.outputs.items():
            array = np.asarray(value, dtype=float).copy()
            if not name or array.shape != shape or np.any(~np.isfinite(array)):
                raise ValueError("every interaction output must be a finite full-grid array")
            array.setflags(write=False); outputs[name] = array
        units = dict(self.output_units)
        if not outputs or set(units) != set(outputs) or any(not unit for unit in units.values()):
            raise ValueError("every interaction output requires one explicit unit")
        bounds = dict(self.bounds)
        if not set(bounds).issubset(outputs):
            raise ValueError("bounds may name only table outputs")
        for name, limits in bounds.items():
            if len(limits) != 2:
                raise ValueError("output bounds require (lower, upper)")
            lower, upper = limits
            if (lower is not None and (not np.isfinite(lower)
                                       or np.any(outputs[name] < float(lower)))):
                raise ValueError(f"output {name!r} violates its lower physical bound")
            if (upper is not None and (not np.isfinite(upper)
                                       or np.any(outputs[name] > float(upper)))):
                raise ValueError(f"output {name!r} violates its upper physical bound")
            if lower is not None and upper is not None and lower > upper:
                raise ValueError("output lower bound exceeds upper bound")
        uncertainty = {}
        for name, value in self.standard_uncertainty.items():
            if name not in outputs:
                raise ValueError("uncertainty may name only table outputs")
            array = np.broadcast_to(np.asarray(value, dtype=float), shape).copy()
            if np.any(~np.isfinite(array)) or np.any(array < 0.0):
                raise ValueError("standard uncertainties must be finite and nonnegative")
            array.setflags(write=False); uncertainty[name] = array
        groups = {name: tuple(items) for name, items in self.conservation_groups.items()}
        for group, names in groups.items():
            if not group or not names or not set(names).issubset(outputs):
                raise ValueError("conservation groups require named table outputs")
            total = sum(outputs[name] for name in names)
            if not np.allclose(total, 1.0, rtol=1e-12, atol=1e-12):
                raise ValueError(f"conservation group {group!r} must sum to one at every node")
        provenance = dict(self.provenance)
        if not provenance.get("source") or not provenance.get("evidence_type"):
            raise ValueError("table provenance requires source and evidence_type")
        # Fingerprinting/serialization must never silently stringify nonportable private objects.
        try:
            json.dumps(provenance, sort_keys=True, separators=(",", ":"))
        except TypeError as error:
            raise ValueError("table provenance must be JSON serializable") from error
        object.__setattr__(self, "axes", axes)
        object.__setattr__(self, "incident_species", species)
        object.__setattr__(self, "outputs", MappingProxyType(outputs))
        object.__setattr__(self, "output_units", MappingProxyType(units))
        object.__setattr__(self, "provenance", MappingProxyType(provenance))
        object.__setattr__(self, "standard_uncertainty", MappingProxyType(uncertainty))
        object.__setattr__(self, "bounds", MappingProxyType(bounds))
        object.__setattr__(self, "conservation_groups", MappingProxyType(groups))
        object.__setattr__(self, "schema_version", int(self.schema_version))

    def to_payload(self):
        return dict(
            schema_version=self.schema_version,
            material=self.material,
            incident_species=list(self.incident_species),
            axes=[dict(
                name=axis.name, values=axis.values.tolist(), unit=axis.unit,
                interpolation=axis.interpolation) for axis in self.axes],
            outputs={name: value.tolist() for name, value in self.outputs.items()},
            output_units=dict(self.output_units),
            provenance=dict(self.provenance),
            standard_uncertainty={
                name: value.tolist() for name, value in self.standard_uncertainty.items()},
            bounds={name: list(value) for name, value in self.bounds.items()},
            conservation_groups={
                name: list(value) for name, value in self.conservation_groups.items()})

    @classmethod
    def from_payload(cls, payload):
        payload = dict(payload)
        return cls(
            material=payload["material"],
            incident_species=tuple(payload["incident_species"]),
            axes=tuple(InteractionAxis(**axis) for axis in payload["axes"]),
            outputs=payload["outputs"], output_units=payload["output_units"],
            provenance=payload["provenance"],
            standard_uncertainty=payload.get("standard_uncertainty", {}),
            bounds={name: tuple(value) for name, value in payload.get("bounds", {}).items()},
            conservation_groups={
                name: tuple(value)
                for name, value in payload.get("conservation_groups", {}).items()},
            schema_version=payload.get("schema_version", 1))

    @property
    def fingerprint(self):
        encoded = json.dumps(
            self.to_payload(), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        return sha256(encoded).hexdigest()

    def evaluate(self, coordinates: Mapping[str, object], *, extrapolation="refuse"):
        coordinates = dict(coordinates)
        if set(coordinates) != {axis.name for axis in self.axes}:
            raise ValueError("coordinates must provide every and only interaction-table axis")
        if extrapolation not in {"refuse", "linear"}:
            raise ValueError("extrapolation must be 'refuse' or explicitly 'linear'")
        broadcast = np.broadcast_arrays(*[
            np.asarray(coordinates[axis.name], dtype=float) for axis in self.axes])
        transformed = []; outside = np.zeros(broadcast[0].shape, dtype=bool); outside_axes = []
        for axis, value in zip(self.axes, broadcast):
            transformed_value = axis.transform(value)
            axis_outside = (value < axis.values[0]) | (value > axis.values[-1])
            if np.any(axis_outside):
                outside_axes.append(axis.name); outside |= axis_outside
            transformed.append(transformed_value)
        if np.any(outside) and extrapolation == "refuse":
            raise SurfaceInteractionDomainError(
                f"interaction coordinates leave validated axes {tuple(outside_axes)}")
        points = np.column_stack([value.ravel() for value in transformed])
        grid = tuple(axis.interpolation_values for axis in self.axes)

        def interpolate(array):
            result = RegularGridInterpolator(
                grid, array, bounds_error=False, fill_value=None)(points)
            return result.reshape(broadcast[0].shape)

        values = {name: interpolate(value) for name, value in self.outputs.items()}
        uncertainty = {
            name: np.maximum(interpolate(value), 0.0)
            for name, value in self.standard_uncertainty.items()}
        for name, value in values.items():
            lower, upper = self.bounds.get(name, (None, None))
            if ((lower is not None and np.any(value < lower - 1e-12))
                    or (upper is not None and np.any(value > upper + 1e-12))):
                raise SurfaceInteractionDomainError(
                    f"interpolated output {name!r} violates its physical bound")
        for group, names in self.conservation_groups.items():
            total = sum(values[name] for name in names)
            if not np.allclose(total, 1.0, rtol=1e-11, atol=1e-11):
                raise RuntimeError(f"interpolated conservation group {group!r} does not close")
        return SurfaceInteractionEvaluation(
            values=values, standard_uncertainty=uncertainty,
            extrapolated_fraction=float(np.mean(outside)) if outside.size else float(outside),
            outside_axes=tuple(outside_axes), table_fingerprint=self.fingerprint)

    def leave_one_out_interpolation_audit(self, output_name):
        """Withhold each interior node of a one-axis table and predict it from its neighbors.

        Endpoint extrapolation is intentionally excluded. The result measures interpolation model
        error at released coordinates; it does not replace or combine with the source's MD/experimental
        standard uncertainty.
        """
        if len(self.axes) != 1:
            raise ValueError("leave-one-out audit currently requires exactly one interaction axis")
        if output_name not in self.outputs:
            raise ValueError(f"unknown interaction output: {output_name!r}")
        axis = self.axes[0]
        if axis.values.size < 3:
            raise ValueError("leave-one-out audit requires at least three axis nodes")
        coordinate = axis.interpolation_values
        observed = self.outputs[output_name]
        held_index = np.arange(1, axis.values.size - 1)
        predicted = np.empty(held_index.size)
        for output_index, withheld in enumerate(held_index):
            retained = np.arange(axis.values.size) != withheld
            predicted[output_index] = np.interp(
                coordinate[withheld], coordinate[retained], observed[retained])
        held_observed = observed[held_index]
        return SurfaceInteractionInterpolationAudit(
            axis_name=axis.name, output_name=output_name,
            coordinates=axis.values[held_index], observed=held_observed,
            predicted=predicted, absolute_error=np.abs(predicted - held_observed),
            table_fingerprint=self.fingerprint)
