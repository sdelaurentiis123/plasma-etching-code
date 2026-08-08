"""Evidence-aware RF-to-plasma absorbed-power boundaries.

A generator setpoint, forward-minus-reflected power, and plasma-absorbed
power are different observables.  This module keeps those nodes separate and
requires downstream hardware losses to be subtracted before an RF
measurement can become a predictive absorbed-power input.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Sequence

import numpy as np

PREDICTIVE_EVIDENCE_KINDS = frozenset({"measured", "validated_model"})
POWER_EVIDENCE_KINDS = frozenset({
    "measured",
    "validated_model",
    "published_range_member",
    "assumed",
    "sensitivity",
    "unresolved",
})


def _validate_evidence(kind: str) -> str:
    value = str(kind)
    if value not in POWER_EVIDENCE_KINDS:
        raise ValueError(f"invalid power evidence kind: {value!r}")
    return value


def _combined_evidence(measurement: str, loss: str) -> str:
    measurement = _validate_evidence(measurement)
    loss = _validate_evidence(loss)
    if measurement == "measured" and loss == "measured":
        return "measured"
    if (
        measurement in PREDICTIVE_EVIDENCE_KINDS
        and loss in PREDICTIVE_EVIDENCE_KINDS
    ):
        return "validated_model"
    for weaker in (
        "unresolved", "assumed", "sensitivity", "published_range_member"
    ):
        if weaker in {measurement, loss}:
            return weaker
    raise RuntimeError("unreachable power-evidence combination")


@dataclass(frozen=True)
class AbsorbedPowerEstimate:
    """A bounded absorbed-power result with its measurement chain exposed."""

    lower_W: float
    upper_W: float
    point_W: float | None
    boundary_kind: str
    measurement_source: str
    loss_source: str
    measurement_evidence: str
    loss_evidence: str
    provenance: Mapping[str, object] = None

    def __post_init__(self):
        bounds = np.asarray([self.lower_W, self.upper_W], dtype=float)
        if (
            np.any(~np.isfinite(bounds))
            or self.lower_W < 0.0
            or self.upper_W <= 0.0
            or self.lower_W > self.upper_W
            or not str(self.boundary_kind).strip()
            or not str(self.measurement_source).strip()
            or not str(self.loss_source).strip()
        ):
            raise ValueError("invalid absorbed-power estimate")
        point = self.point_W
        if point is not None:
            point = float(point)
            if (
                not np.isfinite(point)
                or point <= 0.0
                or point < self.lower_W
                or point > self.upper_W
            ):
                raise ValueError("absorbed-power point must lie within bounds")
        object.__setattr__(self, "lower_W", float(self.lower_W))
        object.__setattr__(self, "upper_W", float(self.upper_W))
        object.__setattr__(self, "point_W", point)
        object.__setattr__(
            self, "measurement_evidence",
            _validate_evidence(self.measurement_evidence))
        object.__setattr__(
            self, "loss_evidence", _validate_evidence(self.loss_evidence))
        object.__setattr__(
            self, "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)))

    @property
    def evidence_kind(self) -> str:
        return _combined_evidence(
            self.measurement_evidence, self.loss_evidence)

    @property
    def supports_prediction(self) -> bool:
        return (
            self.point_W is not None
            and self.measurement_evidence in PREDICTIVE_EVIDENCE_KINDS
            and self.loss_evidence in PREDICTIVE_EVIDENCE_KINDS
        )

    def require_point_W(self) -> float:
        """Return the declared point or fail instead of choosing from a range."""
        if self.point_W is None:
            raise ValueError(
                "power boundary has no point estimate; run the bounds as "
                "explicit sensitivities")
        return self.point_W


@dataclass(frozen=True)
class MeasuredAbsorbedPowerBoundary:
    """A direct calorimetric or otherwise plasma-specific measurement."""

    absorbed_power_W: float
    absolute_uncertainty_W: float
    source: str
    method: str
    evidence_kind: str = "measured"

    def estimate(self) -> AbsorbedPowerEstimate:
        values = np.asarray([
            self.absorbed_power_W, self.absolute_uncertainty_W,
        ], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.absorbed_power_W <= 0.0
            or self.absolute_uncertainty_W < 0.0
            or self.absolute_uncertainty_W >= self.absorbed_power_W
            or not str(self.source).strip()
            or not str(self.method).strip()
        ):
            raise ValueError("invalid measured absorbed-power boundary")
        evidence = _validate_evidence(self.evidence_kind)
        return AbsorbedPowerEstimate(
            lower_W=self.absorbed_power_W - self.absolute_uncertainty_W,
            upper_W=self.absorbed_power_W + self.absolute_uncertainty_W,
            point_W=self.absorbed_power_W,
            boundary_kind="direct_absorbed_power_measurement",
            measurement_source=self.source,
            loss_source=(
                "measurement is specific to plasma absorption; no RF "
                "hardware-loss inference"),
            measurement_evidence=evidence,
            loss_evidence=evidence,
            provenance={"method": self.method},
        )


@dataclass(frozen=True)
class ReactorDiagnosticConditionedPowerFraction:
    """Constant source-to-plasma fraction inferred from a reactor diagnostic.

    This boundary is useful for an out-of-sample equipment transfer, but it is
    deliberately weaker than an absorbed-power measurement. A density or
    temperature diagnostic can identify an effective fraction only through a
    declared reactor model; it cannot turn the source setpoint into measured
    plasma absorption.
    """

    absorbed_fraction: float
    calibration_condition_id: str
    calibration_observable: str
    calibration_source: str
    feature_depth_used: bool = False

    def __post_init__(self):
        fraction = float(self.absorbed_fraction)
        if (
            not np.isfinite(fraction)
            or not 0.0 < fraction <= 1.0
            or not str(self.calibration_condition_id).strip()
            or not str(self.calibration_observable).strip()
            or not str(self.calibration_source).strip()
            or self.feature_depth_used is not False
        ):
            raise ValueError(
                "invalid reactor-diagnostic-conditioned power fraction")
        object.__setattr__(self, "absorbed_fraction", fraction)

    def estimate(
        self,
        source_power_W: float,
        *,
        source: str,
        source_evidence: str = "measured",
    ) -> AbsorbedPowerEstimate:
        power = float(source_power_W)
        if not np.isfinite(power) or power <= 0.0 or not str(source).strip():
            raise ValueError("invalid source-power setpoint")
        absorbed = self.absorbed_fraction * power
        return AbsorbedPowerEstimate(
            lower_W=absorbed,
            upper_W=absorbed,
            point_W=absorbed,
            boundary_kind=(
                "reactor_diagnostic_conditioned_constant_power_fraction"),
            measurement_source=source,
            loss_source=(
                "effective constant source-to-plasma fraction conditioned "
                f"on {self.calibration_observable} at "
                f"{self.calibration_condition_id}; {self.calibration_source}"
            ),
            measurement_evidence=_validate_evidence(source_evidence),
            loss_evidence="sensitivity",
            provenance={
                "absorbed_fraction": self.absorbed_fraction,
                "calibration_condition_id": self.calibration_condition_id,
                "calibration_observable": self.calibration_observable,
                "calibration_source": self.calibration_source,
                "feature_depth_used": False,
            },
        )


def _rf_delivery_estimate(
        *, delivered_power_W: float, hardware_loss_lower_W: float,
        hardware_loss_upper_W: float, hardware_loss_point_W: float | None,
        boundary_kind: str, measurement_source: str, loss_source: str,
        measurement_evidence: str, loss_evidence: str,
        provenance: Mapping[str, object],
        ) -> AbsorbedPowerEstimate:
    values = np.asarray([
        delivered_power_W,
        hardware_loss_lower_W,
        hardware_loss_upper_W,
    ], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or delivered_power_W <= 0.0
        or hardware_loss_lower_W < 0.0
        or hardware_loss_upper_W < hardware_loss_lower_W
        or hardware_loss_upper_W >= delivered_power_W
    ):
        raise ValueError("invalid RF delivery or hardware-loss interval")
    point = None
    if hardware_loss_point_W is not None:
        hardware_loss_point_W = float(hardware_loss_point_W)
        if (
            not np.isfinite(hardware_loss_point_W)
            or hardware_loss_point_W < hardware_loss_lower_W
            or hardware_loss_point_W > hardware_loss_upper_W
        ):
            raise ValueError("hardware-loss point must lie within bounds")
        point = delivered_power_W - hardware_loss_point_W
    return AbsorbedPowerEstimate(
        lower_W=delivered_power_W - hardware_loss_upper_W,
        upper_W=delivered_power_W - hardware_loss_lower_W,
        point_W=point,
        boundary_kind=boundary_kind,
        measurement_source=measurement_source,
        loss_source=loss_source,
        measurement_evidence=measurement_evidence,
        loss_evidence=loss_evidence,
        provenance=provenance,
    )


@dataclass(frozen=True)
class MatchedRFPowerBoundary:
    """Forward/reflected RF power followed by explicit match/coil losses."""

    forward_power_W: float
    reflected_power_W: float
    hardware_loss_lower_W: float
    hardware_loss_upper_W: float
    measurement_source: str
    loss_source: str
    measurement_evidence: str = "measured"
    loss_evidence: str = "unresolved"
    hardware_loss_point_W: float | None = None

    def estimate(self) -> AbsorbedPowerEstimate:
        forward = float(self.forward_power_W)
        reflected = float(self.reflected_power_W)
        if (
            not np.isfinite(forward)
            or not np.isfinite(reflected)
            or forward <= 0.0
            or reflected < 0.0
            or reflected >= forward
        ):
            raise ValueError("reflected power must be below forward power")
        net_power = forward - reflected
        return _rf_delivery_estimate(
            delivered_power_W=net_power,
            hardware_loss_lower_W=self.hardware_loss_lower_W,
            hardware_loss_upper_W=self.hardware_loss_upper_W,
            hardware_loss_point_W=self.hardware_loss_point_W,
            boundary_kind="matched_rf_forward_minus_reflected",
            measurement_source=self.measurement_source,
            loss_source=self.loss_source,
            measurement_evidence=self.measurement_evidence,
            loss_evidence=self.loss_evidence,
            provenance={
                "forward_power_W": forward,
                "reflected_power_W": reflected,
                "net_rf_power_W": net_power,
            },
        )


@dataclass(frozen=True)
class DirectDriveRFPowerBoundary:
    """Calibrated real power at a direct-drive output, before hardware loss."""

    output_real_power_W: float
    hardware_loss_lower_W: float
    hardware_loss_upper_W: float
    measurement_source: str
    loss_source: str
    measurement_evidence: str = "measured"
    loss_evidence: str = "unresolved"
    hardware_loss_point_W: float | None = None

    def estimate(self) -> AbsorbedPowerEstimate:
        return _rf_delivery_estimate(
            delivered_power_W=self.output_real_power_W,
            hardware_loss_lower_W=self.hardware_loss_lower_W,
            hardware_loss_upper_W=self.hardware_loss_upper_W,
            hardware_loss_point_W=self.hardware_loss_point_W,
            boundary_kind="direct_drive_output_real_power",
            measurement_source=self.measurement_source,
            loss_source=self.loss_source,
            measurement_evidence=self.measurement_evidence,
            loss_evidence=self.loss_evidence,
            provenance={"output_real_power_W": self.output_real_power_W},
        )


def time_average_real_power_W(
        voltage_V: Sequence[float], current_A: Sequence[float]) -> float:
    """Return ``mean(v(t) i(t))`` for simultaneous, uniform waveform samples."""
    voltage = np.asarray(voltage_V, dtype=float)
    current = np.asarray(current_A, dtype=float)
    if (
        voltage.ndim != 1
        or current.ndim != 1
        or voltage.size < 2
        or voltage.shape != current.shape
        or np.any(~np.isfinite(voltage))
        or np.any(~np.isfinite(current))
    ):
        raise ValueError(
            "voltage/current waveforms must be equal finite 1-D samples")
    power = float(np.mean(voltage * current))
    if not np.isfinite(power) or power <= 0.0:
        raise ValueError("waveform real power must be positive")
    return power
