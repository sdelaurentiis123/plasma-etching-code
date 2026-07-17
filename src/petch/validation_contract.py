"""Generic calibration and held-out experimental validation contracts.

Campaign-specific observation loaders remain responsible for source checksums,
digitization replay, and physical meaning.  This module supplies the common
claim discipline: parameter roles are declared before calibration, numerical
controls cannot be fitted, held-out rows cannot enter a calibration reveal, and
formal validation requires uncertainty and predictive boundary/mechanism
evidence.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping

import numpy as np


_PARAMETER_ROLES = frozenset({
    "physical_mechanism",
    "reactor_boundary",
    "numerical_control",
    "model_discrepancy",
})
_OBSERVATION_SPLITS = frozenset({
    "boundary_input",
    "calibration",
    "held_out_transfer",
})


def _is_sha256(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _digest(payload):
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


@dataclass(frozen=True)
class ValidationParameter:
    """One declared input with its scientific role and admissible range."""

    name: str
    role: str
    unit: str
    lower: float
    upper: float
    source: str
    calibration_allowed: bool
    supports_prediction: bool

    def __post_init__(self):
        if (
            not str(self.name).strip()
            or self.role not in _PARAMETER_ROLES
            or not str(self.unit).strip()
            or not np.isfinite(self.lower)
            or not np.isfinite(self.upper)
            or self.upper <= self.lower
            or not str(self.source).strip()
            or not isinstance(self.calibration_allowed, (bool, np.bool_))
            or not isinstance(self.supports_prediction, (bool, np.bool_))
            or (
                self.calibration_allowed
                and self.role in {"numerical_control", "model_discrepancy"}
            )
        ):
            raise ValueError("invalid validation parameter")
        object.__setattr__(
            self, "calibration_allowed", bool(self.calibration_allowed))
        object.__setattr__(
            self, "supports_prediction", bool(self.supports_prediction))


@dataclass(frozen=True)
class ValidationObservation:
    """One source-bound experimental observable."""

    observation_id: str
    condition_id: str
    observable: str
    value: float
    unit: str
    split: str
    source: str
    source_sha256: str
    digitization_uncertainty: float = 0.0
    measurement_uncertainty: float | None = None

    def __post_init__(self):
        if (
            not str(self.observation_id).strip()
            or not str(self.condition_id).strip()
            or not str(self.observable).strip()
            or not np.isfinite(self.value)
            or not str(self.unit).strip()
            or self.split not in _OBSERVATION_SPLITS
            or not str(self.source).strip()
            or not _is_sha256(self.source_sha256)
            or not np.isfinite(self.digitization_uncertainty)
            or self.digitization_uncertainty < 0.0
            or (
                self.measurement_uncertainty is not None
                and (
                    not np.isfinite(self.measurement_uncertainty)
                    or self.measurement_uncertainty < 0.0
                )
            )
        ):
            raise ValueError("invalid validation observation")


@dataclass(frozen=True)
class ValidationProtocol:
    """A preregistered split, parameter set, and intended-use statement."""

    protocol_id: str
    revision: str
    intended_use: str
    observations: tuple[ValidationObservation, ...]
    parameters: tuple[ValidationParameter, ...]
    maximum_calibrated_parameters: int
    mechanism_id: str
    boundary_provider_id: str

    def __post_init__(self):
        observations = tuple(self.observations)
        parameters = tuple(self.parameters)
        calibratable = tuple(
            item for item in parameters if item.calibration_allowed)
        if (
            not str(self.protocol_id).strip()
            or not str(self.revision).strip()
            or not str(self.intended_use).strip()
            or not observations
            or any(
                not isinstance(item, ValidationObservation)
                for item in observations
            )
            or len({item.observation_id for item in observations})
            != len(observations)
            or not any(item.split == "calibration" for item in observations)
            or not any(
                item.split == "held_out_transfer" for item in observations)
            or not parameters
            or any(
                not isinstance(item, ValidationParameter)
                for item in parameters
            )
            or len({item.name for item in parameters}) != len(parameters)
            or int(self.maximum_calibrated_parameters)
            != self.maximum_calibrated_parameters
            or self.maximum_calibrated_parameters < 0
            or len(calibratable) > self.maximum_calibrated_parameters
            or not str(self.mechanism_id).strip()
            or not str(self.boundary_provider_id).strip()
        ):
            raise ValueError("invalid validation protocol")
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(
            self,
            "maximum_calibrated_parameters",
            int(self.maximum_calibrated_parameters),
        )

    @property
    def calibration_observations(self):
        return tuple(
            item for item in self.observations if item.split == "calibration")

    @property
    def held_out_observations(self):
        return tuple(
            item
            for item in self.observations
            if item.split == "held_out_transfer"
        )

    @property
    def calibratable_parameters(self):
        return tuple(
            item for item in self.parameters if item.calibration_allowed)

    @property
    def commit_sha256(self):
        return _digest({
            "protocol_id": self.protocol_id,
            "revision": self.revision,
            "intended_use": self.intended_use,
            "observations": [{
                "observation_id": item.observation_id,
                "condition_id": item.condition_id,
                "observable": item.observable,
                "value": item.value,
                "unit": item.unit,
                "split": item.split,
                "source": item.source,
                "source_sha256": item.source_sha256,
                "digitization_uncertainty": item.digitization_uncertainty,
                "measurement_uncertainty": item.measurement_uncertainty,
            } for item in sorted(
                self.observations, key=lambda value: value.observation_id)],
            "parameters": [{
                "name": item.name,
                "role": item.role,
                "unit": item.unit,
                "lower": item.lower,
                "upper": item.upper,
                "source": item.source,
                "calibration_allowed": item.calibration_allowed,
                "supports_prediction": item.supports_prediction,
            } for item in sorted(self.parameters, key=lambda value: value.name)],
            "maximum_calibrated_parameters": self.maximum_calibrated_parameters,
            "mechanism_id": self.mechanism_id,
            "boundary_provider_id": self.boundary_provider_id,
        })


@dataclass(frozen=True)
class CalibrationReveal:
    """Values disclosed after a protocol has frozen its data split and bounds.

    ``operator_epoch_sha256`` identifies the canonical executable-source and
    numerical-operator epoch used for calibration.  It is included in the
    reveal checksum so held-out predictions cannot silently reuse calibrated
    values under a different operator.
    """

    protocol_commit_sha256: str
    operator_epoch_sha256: str
    parameter_values: Mapping[str, float]
    reveal_sha256: str

    @classmethod
    def from_protocol(
        cls, protocol, parameter_values, *, operator_epoch_sha256
    ):
        if not isinstance(protocol, ValidationProtocol):
            raise TypeError("calibration reveal requires a ValidationProtocol")
        if not _is_sha256(operator_epoch_sha256):
            raise ValueError("operator epoch must be a sha256 digest")
        values = {
            str(name): float(value)
            for name, value in dict(parameter_values).items()
        }
        expected = {
            item.name: item for item in protocol.calibratable_parameters}
        if set(values) != set(expected):
            raise ValueError(
                "calibration reveal must contain exactly the preregistered "
                "calibratable parameters")
        if any(
            not np.isfinite(value)
            or value < expected[name].lower
            or value > expected[name].upper
            for name, value in values.items()
        ):
            raise ValueError("calibration reveal lies outside declared bounds")
        reveal_sha256 = _digest({
            "protocol_commit_sha256": protocol.commit_sha256,
            "operator_epoch_sha256": operator_epoch_sha256,
            "parameter_values": {
                name: values[name] for name in sorted(values)},
        })
        return cls(
            protocol_commit_sha256=protocol.commit_sha256,
            operator_epoch_sha256=operator_epoch_sha256,
            parameter_values=values,
            reveal_sha256=reveal_sha256,
        )

    def __post_init__(self):
        values = {
            str(name): float(value)
            for name, value in dict(self.parameter_values).items()
        }
        if (
            not _is_sha256(self.protocol_commit_sha256)
            or not _is_sha256(self.operator_epoch_sha256)
            or any(not name or not np.isfinite(value)
                   for name, value in values.items())
            or not _is_sha256(self.reveal_sha256)
        ):
            raise ValueError("invalid calibration reveal")
        expected_sha256 = _digest({
            "protocol_commit_sha256": self.protocol_commit_sha256,
            "operator_epoch_sha256": self.operator_epoch_sha256,
            "parameter_values": {
                name: values[name] for name in sorted(values)},
        })
        if self.reveal_sha256 != expected_sha256:
            raise ValueError("calibration reveal checksum does not match its content")
        object.__setattr__(
            self, "parameter_values", MappingProxyType(values))


@dataclass(frozen=True)
class HeldOutPrediction:
    """One prediction for an observation that was never used in calibration."""

    observation_id: str
    predicted_value: float
    numerical_uncertainty: float
    parameter_uncertainty: float
    model_discrepancy_uncertainty: float
    within_declared_scope: bool
    boundary_supports_prediction: bool
    mechanism_supports_prediction: bool
    run_manifest_sha256: str
    calibration_reveal_sha256: str
    operator_epoch_sha256: str

    def __post_init__(self):
        uncertainty = np.asarray([
            self.numerical_uncertainty,
            self.parameter_uncertainty,
            self.model_discrepancy_uncertainty,
        ], dtype=float)
        if (
            not str(self.observation_id).strip()
            or not np.isfinite(self.predicted_value)
            or np.any(~np.isfinite(uncertainty))
            or np.any(uncertainty < 0.0)
            or not isinstance(self.within_declared_scope, (bool, np.bool_))
            or not isinstance(
                self.boundary_supports_prediction, (bool, np.bool_))
            or not isinstance(
                self.mechanism_supports_prediction, (bool, np.bool_))
            or not _is_sha256(self.run_manifest_sha256)
            or not _is_sha256(self.calibration_reveal_sha256)
            or not _is_sha256(self.operator_epoch_sha256)
        ):
            raise ValueError("invalid held-out prediction")

    @property
    def prediction_uncertainty_bound(self):
        return (
            self.numerical_uncertainty
            + self.parameter_uncertainty
            + self.model_discrepancy_uncertainty
        )


@dataclass(frozen=True)
class ValidationScore:
    """Conservative held-out score with decomposed failure reasons."""

    formal_validation_passed: bool
    uncertainty_complete: bool
    absolute_error: Mapping[str, float]
    combined_uncertainty_bound: Mapping[str, float]
    reasons: tuple[str, ...]
    protocol_commit_sha256: str
    calibration_reveal_sha256: str
    operator_epoch_sha256: str

    def __post_init__(self):
        error = MappingProxyType({
            str(name): float(value)
            for name, value in dict(self.absolute_error).items()
        })
        bound = MappingProxyType({
            str(name): float(value)
            for name, value in dict(self.combined_uncertainty_bound).items()
        })
        if (
            not error
            or set(error) != set(bound)
            or any(
                not np.isfinite(value) or value < 0.0
                for value in (*error.values(), *bound.values())
            )
            or not _is_sha256(self.protocol_commit_sha256)
            or not _is_sha256(self.calibration_reveal_sha256)
            or not _is_sha256(self.operator_epoch_sha256)
        ):
            raise ValueError("invalid validation score")
        object.__setattr__(
            self, "formal_validation_passed",
            bool(self.formal_validation_passed),
        )
        object.__setattr__(
            self, "uncertainty_complete", bool(self.uncertainty_complete))
        object.__setattr__(self, "absolute_error", error)
        object.__setattr__(self, "combined_uncertainty_bound", bound)
        object.__setattr__(self, "reasons", tuple(self.reasons))


def score_held_out_predictions(protocol, reveal, predictions):
    """Score every and only held-out observation from the committed protocol."""
    if (
        not isinstance(protocol, ValidationProtocol)
        or not isinstance(reveal, CalibrationReveal)
        or reveal.protocol_commit_sha256 != protocol.commit_sha256
    ):
        raise ValueError("score requires a matching protocol and calibration reveal")
    predictions = tuple(predictions)
    if any(not isinstance(item, HeldOutPrediction) for item in predictions):
        raise TypeError("predictions must be HeldOutPrediction values")
    by_id = {item.observation_id: item for item in predictions}
    expected = {
        item.observation_id for item in protocol.held_out_observations}
    if (
        len(by_id) != len(predictions)
        or set(by_id) != expected
        or any(
            item.calibration_reveal_sha256 != reveal.reveal_sha256
            for item in predictions
        )
    ):
        raise ValueError(
            "predictions must cover every and only held-out observations")
    if any(
        item.operator_epoch_sha256 != reveal.operator_epoch_sha256
        for item in predictions
    ):
        raise ValueError(
            "held-out prediction operator epoch does not match calibration"
        )

    error = {}
    bound = {}
    reasons = []
    uncertainty_complete = True
    for observation in protocol.held_out_observations:
        prediction = by_id[observation.observation_id]
        measurement = observation.measurement_uncertainty
        if measurement is None:
            uncertainty_complete = False
            measurement = 0.0
            reasons.append(
                f"{observation.observation_id}: measurement uncertainty is missing")
        total_bound = (
            observation.digitization_uncertainty
            + measurement
            + prediction.prediction_uncertainty_bound
        )
        absolute_error = abs(prediction.predicted_value - observation.value)
        error[observation.observation_id] = absolute_error
        bound[observation.observation_id] = total_bound
        if absolute_error > total_bound:
            reasons.append(
                f"{observation.observation_id}: error exceeds combined uncertainty")
        if not prediction.within_declared_scope:
            reasons.append(
                f"{observation.observation_id}: prediction is outside declared scope")
        if not prediction.boundary_supports_prediction:
            reasons.append(
                f"{observation.observation_id}: reactor boundary is nonpredictive")
        if not prediction.mechanism_supports_prediction:
            reasons.append(
                f"{observation.observation_id}: surface mechanism is nonpredictive")

    nonpredictive_calibrated = tuple(
        item.name
        for item in protocol.calibratable_parameters
        if not item.supports_prediction
    )
    if nonpredictive_calibrated:
        reasons.append(
            "calibrated parameters lack predictive evidence: "
            + ", ".join(nonpredictive_calibrated)
        )
    return ValidationScore(
        formal_validation_passed=not reasons,
        uncertainty_complete=uncertainty_complete,
        absolute_error=error,
        combined_uncertainty_bound=bound,
        reasons=tuple(reasons),
        protocol_commit_sha256=protocol.commit_sha256,
        calibration_reveal_sha256=reveal.reveal_sha256,
        operator_epoch_sha256=reveal.operator_epoch_sha256,
    )
