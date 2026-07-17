"""Versioned reactor-to-feature coupling contracts.

The feature engine consumes :class:`~petch.boundary_state.PlasmaBoundaryState`.
This module defines how a reactor calculation, diagnostic reconstruction, or
reduced-order surrogate is allowed to produce that state.  Providers are bound
to a complete operating-point query and must report applicability, uncertainty,
and provenance; a boundary for one recipe cannot silently be reused for another.

The reverse contract is intentionally compact.  A feature calculation may
return area-averaged wall-loss probabilities and product fluxes to a reactor
model without exposing its internal mesh or surface-state representation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np

from .boundary_state import PlasmaBoundaryState


def _is_sha256(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _immutable_scalar_mapping(values, *, field_name, nonnegative=False):
    output = {}
    for name, value in dict(values).items():
        number = float(value)
        if (
            not isinstance(name, str)
            or not name
            or not np.isfinite(number)
            or (nonnegative and number < 0.0)
        ):
            raise ValueError(f"invalid {field_name}")
        output[name] = number
    return MappingProxyType(output)


def _immutable_string_mapping(values, *, field_name):
    output = {}
    for name, value in dict(values).items():
        if (
            not isinstance(name, str)
            or not name
            or not isinstance(value, str)
            or not value
        ):
            raise ValueError(f"invalid {field_name}")
        output[name] = value
    return MappingProxyType(output)


def _json_value(value):
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError("manifest values must be finite")
        return number
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, Mapping):
        return {
            str(key): _json_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"unsupported manifest value: {type(value).__name__}")


def _digest(payload):
    encoded = json.dumps(
        _json_value(payload),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def _array_manifest(value):
    array = np.ascontiguousarray(np.asarray(value))
    return {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "sha256": sha256(array.tobytes(order="C")).hexdigest(),
    }


def _boundary_sha256(boundary):
    provenance = dict(boundary.provenance)
    provenance.pop("reactor_coupling", None)
    return _digest({
        "reference_plane_m": float(boundary.reference_plane_m),
        "provenance": provenance,
        "species": [{
            "name": species.name,
            "charge_number": int(species.charge_number),
            "mass_amu": float(species.mass_amu),
            "flux_m2_s": float(species.flux_m2_s),
            "velocity_sqrt_eV": _array_manifest(species.velocity_sqrt_eV),
            "weight": _array_manifest(species.weight),
            "phase_rad": (
                None if species.phase_rad is None
                else _array_manifest(species.phase_rad)
            ),
            "position_m": (
                None if species.position_m is None
                else _array_manifest(species.position_m)
            ),
            "density_model": (
                None if species.density_model is None
                else type(species.density_model).__name__
            ),
            "density_model_2d": (
                None if species.density_model_2d is None
                else type(species.density_model_2d).__name__
            ),
            "provenance": dict(species.provenance),
        } for species in boundary.species],
    })


@dataclass(frozen=True)
class ReactorBoundaryQuery:
    """One complete operating point requested from a reactor provider.

    ``recipe`` contains dimensional scalar controls and ``recipe_units`` must
    cover the same names.  The schema is deliberately equipment-agnostic:
    pressure, gas flows, source power, bias harmonics, magnetic fields, and
    temperatures can all be represented without teaching the feature engine
    reactor-specific vocabulary.
    """

    condition_id: str
    tool_id: str
    recipe: Mapping[str, float]
    recipe_units: Mapping[str, str]
    wafer_position_m: tuple[float, float]
    process_time_s: float
    substrate_temperature_K: float
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        recipe = _immutable_scalar_mapping(self.recipe, field_name="reactor recipe")
        units = _immutable_string_mapping(
            self.recipe_units, field_name="reactor recipe units")
        position = tuple(float(value) for value in self.wafer_position_m)
        if (
            not str(self.condition_id).strip()
            or not str(self.tool_id).strip()
            or not recipe
            or set(recipe) != set(units)
            or len(position) != 2
            or np.any(~np.isfinite(position))
            or not np.isfinite(self.process_time_s)
            or self.process_time_s < 0.0
            or not np.isfinite(self.substrate_temperature_K)
            or self.substrate_temperature_K <= 0.0
        ):
            raise ValueError("invalid reactor boundary query")
        provenance = MappingProxyType(dict(self.provenance))
        _json_value(provenance)
        object.__setattr__(self, "recipe", recipe)
        object.__setattr__(self, "recipe_units", units)
        object.__setattr__(self, "wafer_position_m", position)
        object.__setattr__(self, "process_time_s", float(self.process_time_s))
        object.__setattr__(
            self, "substrate_temperature_K", float(self.substrate_temperature_K))
        object.__setattr__(self, "provenance", provenance)

    @property
    def manifest(self):
        return {
            "condition_id": self.condition_id,
            "tool_id": self.tool_id,
            "recipe": dict(self.recipe),
            "recipe_units": dict(self.recipe_units),
            "wafer_position_m": list(self.wafer_position_m),
            "process_time_s": self.process_time_s,
            "substrate_temperature_K": self.substrate_temperature_K,
            "provenance": dict(self.provenance),
        }

    @property
    def query_sha256(self):
        return _digest(self.manifest)


@dataclass(frozen=True)
class SurfaceFeedbackState:
    """Area-averaged feature-to-reactor feedback at one operating point.

    ``species_loss_probability`` is the net probability that an incident gas
    species is removed from the gas phase after all feature-scale re-emission
    and reflection.  ``product_flux_m2_s`` contains products returned to the
    reactor.  Both quantities are homogenized over ``reference_area_m2``.
    """

    species_loss_probability: Mapping[str, float]
    product_flux_m2_s: Mapping[str, float]
    net_current_density_A_m2: float
    surface_temperature_K: float
    reference_area_m2: float
    provenance: Mapping[str, object] = field(default_factory=dict)
    supports_prediction: bool = False

    def __post_init__(self):
        loss = _immutable_scalar_mapping(
            self.species_loss_probability,
            field_name="surface species-loss probability",
            nonnegative=True,
        )
        products = _immutable_scalar_mapping(
            self.product_flux_m2_s,
            field_name="surface product flux",
            nonnegative=True,
        )
        if (
            any(value > 1.0 for value in loss.values())
            or not np.isfinite(self.net_current_density_A_m2)
            or not np.isfinite(self.surface_temperature_K)
            or self.surface_temperature_K <= 0.0
            or not np.isfinite(self.reference_area_m2)
            or self.reference_area_m2 <= 0.0
            or not isinstance(self.supports_prediction, (bool, np.bool_))
        ):
            raise ValueError("invalid surface feedback state")
        provenance = MappingProxyType(dict(self.provenance))
        _json_value(provenance)
        object.__setattr__(self, "species_loss_probability", loss)
        object.__setattr__(self, "product_flux_m2_s", products)
        object.__setattr__(
            self, "net_current_density_A_m2", float(self.net_current_density_A_m2))
        object.__setattr__(
            self, "surface_temperature_K", float(self.surface_temperature_K))
        object.__setattr__(self, "reference_area_m2", float(self.reference_area_m2))
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "supports_prediction", bool(self.supports_prediction))

    @property
    def manifest(self):
        return {
            "species_loss_probability": dict(self.species_loss_probability),
            "product_flux_m2_s": dict(self.product_flux_m2_s),
            "net_current_density_A_m2": self.net_current_density_A_m2,
            "surface_temperature_K": self.surface_temperature_K,
            "reference_area_m2": self.reference_area_m2,
            "provenance": dict(self.provenance),
            "supports_prediction": self.supports_prediction,
        }

    @property
    def feedback_sha256(self):
        return _digest(self.manifest)


@dataclass(frozen=True)
class ReactorBoundaryPrediction:
    """A provider result carrying the kinetic boundary and its evidence."""

    boundary: PlasmaBoundaryState
    query_sha256: str
    provider_name: str
    provider_version: str
    supports_prediction: bool
    applicability_reasons: tuple[str, ...] = ()
    relative_standard_uncertainty: Mapping[str, float] = field(default_factory=dict)
    provenance: Mapping[str, object] = field(default_factory=dict)
    feedback_sha256: str | None = None

    def __post_init__(self):
        if (
            not isinstance(self.boundary, PlasmaBoundaryState)
            or not _is_sha256(self.query_sha256)
            or not str(self.provider_name).strip()
            or not str(self.provider_version).strip()
            or not isinstance(self.supports_prediction, (bool, np.bool_))
            or (
                self.feedback_sha256 is not None
                and not _is_sha256(self.feedback_sha256)
            )
        ):
            raise ValueError("invalid reactor boundary prediction")
        uncertainty = _immutable_scalar_mapping(
            self.relative_standard_uncertainty,
            field_name="reactor boundary uncertainty",
            nonnegative=True,
        )
        provenance = MappingProxyType(dict(self.provenance))
        _json_value(provenance)
        boundary_support = bool(
            dict(self.boundary.provenance).get("supports_prediction", False))
        if self.supports_prediction and not boundary_support:
            raise ValueError(
                "provider cannot promote a boundary whose own provenance is nonpredictive")
        object.__setattr__(
            self, "supports_prediction", bool(self.supports_prediction))
        object.__setattr__(
            self, "applicability_reasons", tuple(self.applicability_reasons))
        object.__setattr__(
            self, "relative_standard_uncertainty", uncertainty)
        object.__setattr__(self, "provenance", provenance)

    @property
    def manifest(self):
        return {
            "query_sha256": self.query_sha256,
            "feedback_sha256": self.feedback_sha256,
            "boundary_sha256": _boundary_sha256(self.boundary),
            "provider_name": self.provider_name,
            "provider_version": self.provider_version,
            "supports_prediction": self.supports_prediction,
            "applicability_reasons": list(self.applicability_reasons),
            "relative_standard_uncertainty": dict(
                self.relative_standard_uncertainty),
            "provenance": dict(self.provenance),
        }


class ReactorBoundaryProvider(Protocol):
    """Protocol implemented by reactor solvers, diagnostics, and surrogates."""

    def predict_boundary(
        self,
        query: ReactorBoundaryQuery,
        feedback: SurfaceFeedbackState | None = None,
    ) -> ReactorBoundaryPrediction:
        ...


@dataclass(frozen=True)
class BoundReactorBoundaryProvider:
    """Bind an existing boundary to exactly one query and optional feedback.

    This adapter is the safe entry point for published HPEM decks and measured
    distributions.  It deliberately refuses interpolation; a future surrogate
    provider must declare and test its own applicability domain.
    """

    query_sha256: str
    boundary: PlasmaBoundaryState
    provider_name: str
    provider_version: str
    supports_prediction: bool
    relative_standard_uncertainty: Mapping[str, float] = field(default_factory=dict)
    provenance: Mapping[str, object] = field(default_factory=dict)
    feedback_sha256: str | None = None

    def __post_init__(self):
        prediction = ReactorBoundaryPrediction(
            boundary=self.boundary,
            query_sha256=self.query_sha256,
            provider_name=self.provider_name,
            provider_version=self.provider_version,
            supports_prediction=self.supports_prediction,
            relative_standard_uncertainty=self.relative_standard_uncertainty,
            provenance=self.provenance,
            feedback_sha256=self.feedback_sha256,
        )
        object.__setattr__(
            self, "relative_standard_uncertainty",
            prediction.relative_standard_uncertainty,
        )
        object.__setattr__(self, "provenance", prediction.provenance)
        object.__setattr__(
            self, "supports_prediction", prediction.supports_prediction)

    def predict_boundary(self, query, feedback=None):
        if not isinstance(query, ReactorBoundaryQuery):
            raise TypeError("reactor provider requires a ReactorBoundaryQuery")
        supplied_feedback_sha256 = (
            None if feedback is None else feedback.feedback_sha256)
        if query.query_sha256 != self.query_sha256:
            raise ValueError(
                "bound reactor boundary does not match the requested operating point")
        if supplied_feedback_sha256 != self.feedback_sha256:
            raise ValueError(
                "bound reactor boundary does not match the requested surface feedback")
        return ReactorBoundaryPrediction(
            boundary=self.boundary,
            query_sha256=self.query_sha256,
            provider_name=self.provider_name,
            provider_version=self.provider_version,
            supports_prediction=self.supports_prediction,
            relative_standard_uncertainty=self.relative_standard_uncertainty,
            provenance=self.provenance,
            feedback_sha256=self.feedback_sha256,
        )


def resolve_reactor_boundary(
    provider: ReactorBoundaryProvider,
    query: ReactorBoundaryQuery,
    feedback: SurfaceFeedbackState | None = None,
    *,
    claim_mode: str = "development",
):
    """Resolve a provider while enforcing the requested scientific claim mode."""
    if claim_mode not in {"development", "predictive"}:
        raise ValueError("claim_mode must be development or predictive")
    prediction = provider.predict_boundary(query, feedback)
    if not isinstance(prediction, ReactorBoundaryPrediction):
        raise TypeError("reactor provider returned the wrong result type")
    if prediction.query_sha256 != query.query_sha256:
        raise ValueError("reactor provider returned a boundary for a different query")
    expected_feedback = None if feedback is None else feedback.feedback_sha256
    if prediction.feedback_sha256 != expected_feedback:
        raise ValueError("reactor provider returned a boundary for different feedback")
    if claim_mode == "predictive" and not prediction.supports_prediction:
        raise ValueError(
            "predictive feature runs require a predictive reactor boundary")
    coupling_manifest = prediction.manifest
    existing = dict(prediction.boundary.provenance).get("reactor_coupling")
    if existing is not None and existing != coupling_manifest:
        raise ValueError("boundary already carries a different reactor-coupling manifest")
    boundary = PlasmaBoundaryState(
        species=prediction.boundary.species,
        reference_plane_m=prediction.boundary.reference_plane_m,
        provenance=dict(
            prediction.boundary.provenance,
            reactor_coupling=coupling_manifest,
        ),
    )
    return ReactorBoundaryPrediction(
        boundary=boundary,
        query_sha256=prediction.query_sha256,
        provider_name=prediction.provider_name,
        provider_version=prediction.provider_version,
        supports_prediction=prediction.supports_prediction,
        applicability_reasons=prediction.applicability_reasons,
        relative_standard_uncertainty=prediction.relative_standard_uncertainty,
        provenance=prediction.provenance,
        feedback_sha256=prediction.feedback_sha256,
    )
