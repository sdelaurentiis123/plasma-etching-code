"""Evidence-gated DC self-bias transfer for the Oxford RIE-80 family.

Generator forward power is not an electrode-voltage boundary condition.  This
module keeps the sparse measurements that do exist as interval/censored data
and constructs deterministic *sensitivity histories* for an unmeasured target
condition.  It intentionally does not turn equipment-family measurements into
a probability distribution or an absolute-depth certification.

All voltages use the magnitude of the negative powered-electrode self-bias.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import math
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Sequence

import numpy as np


_BIAS_RELATIONS = frozenset({
    "measured_point",
    "measured_range",
    "greater_than",
    "less_than_approx",
})
_TOOL_RELATIONS = frozenset({"exact_ngp80", "oxford_rie80_family"})


def _optional_float(value: str | float | None) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("non-finite numeric evidence value")
    return result


@dataclass(frozen=True)
class Oxford80RIECondition:
    """Machine settings that are actually known at recipe level."""

    tool_model: str
    rf_power_W: float
    pressure_mTorr: float
    gas_flows_sccm: Mapping[str, float]
    electrode_temperature_C: float | None = None
    duration_s: float | None = None

    def __post_init__(self):
        flows = {
            str(species): float(flow)
            for species, flow in self.gas_flows_sccm.items()
            if float(flow) > 0.0
        }
        scalars = (self.rf_power_W, self.pressure_mTorr)
        if (
            not str(self.tool_model).strip()
            or any(not math.isfinite(value) or value <= 0.0 for value in scalars)
            or not flows
            or any(not species.strip() or not math.isfinite(flow)
                   for species, flow in flows.items())
            or (
                self.electrode_temperature_C is not None
                and not math.isfinite(self.electrode_temperature_C)
            )
            or (
                self.duration_s is not None
                and (not math.isfinite(self.duration_s) or self.duration_s <= 0.0)
            )
        ):
            raise ValueError("invalid Oxford RIE condition")
        object.__setattr__(self, "gas_flows_sccm", MappingProxyType(flows))

    @property
    def reduced_drive_W_per_mTorr(self) -> float:
        """Reported forward-power setpoint divided by chamber pressure."""

        return self.rf_power_W / self.pressure_mTorr

    @property
    def active_gases(self) -> frozenset[str]:
        return frozenset(self.gas_flows_sccm)


@dataclass(frozen=True)
class Oxford80SelfBiasObservation:
    """One primary-source voltage statement, including censor semantics."""

    source_id: str
    tool_model: str
    tool_relation: str
    rf_power_W: float
    pressure_mTorr: float
    gas_flows_sccm: Mapping[str, float | None]
    bias_relation: str
    bias_lower_V: float | None
    bias_upper_V: float | None
    run_phase: str
    electrode_temperature_C: float | None
    source_page: str
    source_pdf_sha256: str
    source_locator: str
    loading_context: str

    def __post_init__(self):
        flows = {
            str(species): (
                None if flow is None else float(flow)
            )
            for species, flow in self.gas_flows_sccm.items()
        }
        bounds = (self.bias_lower_V, self.bias_upper_V)
        if (
            not str(self.source_id).strip()
            or not str(self.tool_model).strip()
            or self.tool_relation not in _TOOL_RELATIONS
            or self.bias_relation not in _BIAS_RELATIONS
            or not math.isfinite(self.rf_power_W)
            or self.rf_power_W <= 0.0
            or not math.isfinite(self.pressure_mTorr)
            or self.pressure_mTorr <= 0.0
            or not flows
            or any(
                not species.strip()
                or (flow is not None and (not math.isfinite(flow) or flow < 0.0))
                for species, flow in flows.items()
            )
            or any(value is not None and
                   (not math.isfinite(value) or value <= 0.0)
                   for value in bounds)
            or not str(self.run_phase).strip()
            or not str(self.source_page).strip()
            or len(str(self.source_pdf_sha256)) != 64
            or not str(self.source_locator).strip()
        ):
            raise ValueError("invalid Oxford RIE-80 self-bias observation")
        if self.bias_relation == "measured_point":
            valid_bounds = (
                self.bias_lower_V is not None
                and self.bias_upper_V is not None
                and self.bias_lower_V == self.bias_upper_V
            )
        elif self.bias_relation == "measured_range":
            valid_bounds = (
                self.bias_lower_V is not None
                and self.bias_upper_V is not None
                and self.bias_lower_V <= self.bias_upper_V
            )
        elif self.bias_relation == "greater_than":
            valid_bounds = (
                self.bias_lower_V is not None and self.bias_upper_V is None
            )
        else:
            valid_bounds = (
                self.bias_lower_V is None and self.bias_upper_V is not None
            )
        if not valid_bounds:
            raise ValueError("bias bounds contradict their censor relation")
        object.__setattr__(self, "gas_flows_sccm", MappingProxyType(flows))

    @property
    def reduced_drive_W_per_mTorr(self) -> float:
        return self.rf_power_W / self.pressure_mTorr

    @property
    def active_gases(self) -> frozenset[str]:
        return frozenset(
            species for species, flow in self.gas_flows_sccm.items()
            if flow is None or flow > 0.0
        )

    @property
    def midpoint_V(self) -> float:
        if self.bias_relation not in {"measured_point", "measured_range"}:
            raise ValueError("a censored observation has no measured midpoint")
        assert self.bias_lower_V is not None
        assert self.bias_upper_V is not None
        return 0.5 * (self.bias_lower_V + self.bias_upper_V)


@dataclass(frozen=True)
class SelfBiasSensitivityHistory:
    """Deterministic voltage history used only for propagated sensitivity."""

    name: str
    time_s: np.ndarray
    bias_magnitude_V: np.ndarray
    source_ids: tuple[str, ...]
    interpretation: str
    endpoints_are_censor_thresholds: bool = False
    measured_on_target_condition: bool = False
    supports_absolute_depth_prediction: bool = False

    def __post_init__(self):
        time = np.asarray(self.time_s, dtype=float).copy()
        bias = np.asarray(self.bias_magnitude_V, dtype=float).copy()
        if (
            not str(self.name).strip()
            or time.ndim != 1
            or time.size < 2
            or bias.shape != time.shape
            or np.any(~np.isfinite(time))
            or np.any(~np.isfinite(bias))
            or time[0] != 0.0
            or np.any(np.diff(time) <= 0.0)
            or np.any(bias <= 0.0)
            or not self.source_ids
            or any(not str(source_id).strip() for source_id in self.source_ids)
            or not str(self.interpretation).strip()
            or self.measured_on_target_condition
            or self.supports_absolute_depth_prediction
        ):
            raise ValueError("invalid self-bias sensitivity history")
        time.setflags(write=False)
        bias.setflags(write=False)
        object.__setattr__(self, "time_s", time)
        object.__setattr__(self, "bias_magnitude_V", bias)

    def at(self, time_s: float | np.ndarray) -> float | np.ndarray:
        query = np.asarray(time_s, dtype=float)
        if np.any(~np.isfinite(query)):
            raise ValueError("history query time must be finite")
        if np.any(query < self.time_s[0]) or np.any(query > self.time_s[-1]):
            raise ValueError("history query lies outside the recipe duration")
        value = np.interp(query, self.time_s, self.bias_magnitude_V)
        return float(value) if value.ndim == 0 else value


@dataclass(frozen=True)
class Oxford80SelfBiasTransfer:
    """Target-free transfer result with explicit non-prediction gates."""

    target: Oxford80RIECondition
    observations: tuple[Oxford80SelfBiasObservation, ...]
    matched_chemistry_reduced_drive_source_id: str
    matched_chemistry_reduced_drive_anchor_V: float
    printed_reference_window_V: tuple[float, float]
    histories: tuple[SelfBiasSensitivityHistory, ...]
    sem_target_used: bool = False
    measured_target_bias_used: bool = False
    printed_window_is_probability_interval: bool = False
    censored_data_extend_outside_printed_window: bool = True
    supports_unique_target_bias: bool = False
    supports_absolute_depth_prediction: bool = False

    def __post_init__(self):
        source_ids = {item.source_id for item in self.observations}
        lo, hi = self.printed_reference_window_V
        if (
            not self.observations
            or len(source_ids) != len(self.observations)
            or self.matched_chemistry_reduced_drive_source_id not in source_ids
            or not math.isfinite(self.matched_chemistry_reduced_drive_anchor_V)
            or self.matched_chemistry_reduced_drive_anchor_V <= 0.0
            or not math.isfinite(lo)
            or not math.isfinite(hi)
            or lo <= 0.0
            or hi <= lo
            or not self.histories
            or self.sem_target_used
            or self.measured_target_bias_used
            or self.printed_window_is_probability_interval
            or not self.censored_data_extend_outside_printed_window
            or self.supports_unique_target_bias
            or self.supports_absolute_depth_prediction
        ):
            raise ValueError("invalid evidence-gated self-bias transfer")


def load_oxford80_self_bias_evidence(
    path: str | Path,
) -> tuple[Oxford80SelfBiasObservation, ...]:
    """Load checksum-bound primary observations from the committed CSV."""

    observations = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            flows = {
                gas: _optional_float(row.get(f"{gas}_sccm"))
                for gas in ("CHF3", "SF6", "O2", "Ar")
                if str(row.get(f"{gas}_sccm", "")).strip()
            }
            observations.append(Oxford80SelfBiasObservation(
                source_id=row["source_id"],
                tool_model=row["tool_model"],
                tool_relation=row["tool_relation"],
                rf_power_W=float(row["rf_power_W"]),
                pressure_mTorr=float(row["pressure_mTorr"]),
                gas_flows_sccm=flows,
                bias_relation=row["bias_relation"],
                bias_lower_V=_optional_float(row["bias_lower_V"]),
                bias_upper_V=_optional_float(row["bias_upper_V"]),
                run_phase=row["run_phase"],
                electrode_temperature_C=_optional_float(
                    row["electrode_temperature_C"]),
                source_page=row["source_page"],
                source_pdf_sha256=row["source_pdf_sha256"],
                source_locator=row["source_locator"],
                loading_context=row["loading_context"],
            ))
    if not observations:
        raise ValueError("Oxford RIE-80 evidence table is empty")
    return tuple(observations)


def build_oxford80_self_bias_transfer(
    target: Oxford80RIECondition,
    observations: Sequence[Oxford80SelfBiasObservation],
) -> Oxford80SelfBiasTransfer:
    """Build target-free voltage witnesses without fitting an etch outcome.

    The central anchor is selected mechanically: same active-gas set, measured
    (not censored) voltage, then minimum logarithmic separation in ``P/p``.
    Censored exact-tool statements remain censored and are represented only by
    a threshold history; they are not silently converted into endpoint data.
    """

    items = tuple(observations)
    candidates = [
        item for item in items
        if item.bias_relation in {"measured_point", "measured_range"}
        and item.active_gases == target.active_gases
    ]
    if not candidates:
        raise ValueError("no measured same-chemistry Oxford RIE-80 bias anchor")
    anchor = min(
        candidates,
        key=lambda item: abs(math.log(
            item.reduced_drive_W_per_mTorr
            / target.reduced_drive_W_per_mTorr
        )),
    )
    duration = target.duration_s
    if duration is None:
        raise ValueError("target duration is required for bias histories")

    exact_start = next(
        item for item in items
        if item.tool_relation == "exact_ngp80"
        and item.run_phase == "start"
        and item.bias_relation == "greater_than"
    )
    exact_end = next(
        item for item in items
        if item.tool_relation == "exact_ngp80"
        and item.run_phase == "end"
        and item.bias_relation == "less_than_approx"
    )
    assert exact_start.bias_lower_V is not None
    assert exact_end.bias_upper_V is not None

    same_pressure = next(
        item for item in items
        if item.source_id == "zou-2004-plasmalab80-sf6-o2-30mtorr"
    )
    high_family = next(
        item for item in items
        if item.source_id == "plattner-2003-plasmalab80-chf3-ar"
    )
    assert same_pressure.bias_lower_V is not None
    assert same_pressure.bias_upper_V is not None
    assert high_family.bias_upper_V is not None

    histories = (
        SelfBiasSensitivityHistory(
            name="same-chemistry-reduced-drive anchor",
            time_s=np.array([0.0, duration]),
            bias_magnitude_V=np.array([anchor.midpoint_V, anchor.midpoint_V]),
            source_ids=(anchor.source_id,),
            interpretation=(
                "Exact active-gas set and exact P/p match in the Oxford RIE-80 "
                "family; mixture fractions and target material differ."
            ),
        ),
        SelfBiasSensitivityHistory(
            name="exact-NGP80 conditioning thresholds",
            time_s=np.array([0.0, duration]),
            bias_magnitude_V=np.array([
                exact_start.bias_lower_V,
                exact_end.bias_upper_V,
            ]),
            source_ids=(exact_start.source_id, exact_end.source_id),
            interpretation=(
                "Linear connection between printed censor thresholds only; "
                "the source states the true start is higher and end lower."
            ),
            endpoints_are_censor_thresholds=True,
        ),
        SelfBiasSensitivityHistory(
            name="same-pressure SF6/O2 low witness",
            time_s=np.array([0.0, duration]),
            bias_magnitude_V=np.array([
                same_pressure.bias_lower_V,
                same_pressure.bias_lower_V,
            ]),
            source_ids=(same_pressure.source_id,),
            interpretation=(
                "Low end of a measured 30 mTorr O2-sweep range on a "
                "PlasmaLab 80+; chemistry and loading differ."
            ),
        ),
        SelfBiasSensitivityHistory(
            name="same-pressure SF6/O2 high witness",
            time_s=np.array([0.0, duration]),
            bias_magnitude_V=np.array([
                same_pressure.bias_upper_V,
                same_pressure.bias_upper_V,
            ]),
            source_ids=(same_pressure.source_id,),
            interpretation=(
                "High end of a measured 30 mTorr O2-sweep range on a "
                "PlasmaLab 80+; chemistry and loading differ."
            ),
        ),
        SelfBiasSensitivityHistory(
            name="same-family CHF3/Ar high witness",
            time_s=np.array([0.0, duration]),
            bias_magnitude_V=np.array([
                high_family.bias_upper_V,
                high_family.bias_upper_V,
            ]),
            source_ids=(high_family.source_id,),
            interpretation=(
                "Measured CHF3/Ar value on a PlasmaLab RIE 80; pressure, "
                "composition, loading, and exact machine differ."
            ),
        ),
    )

    finite_bounds = [
        value
        for item in items
        for value in (item.bias_lower_V, item.bias_upper_V)
        if value is not None
    ]
    return Oxford80SelfBiasTransfer(
        target=target,
        observations=items,
        matched_chemistry_reduced_drive_source_id=anchor.source_id,
        matched_chemistry_reduced_drive_anchor_V=anchor.midpoint_V,
        printed_reference_window_V=(min(finite_bounds), max(finite_bounds)),
        histories=histories,
    )
