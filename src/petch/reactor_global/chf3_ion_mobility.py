"""Measured CHF2+ mobility in CHF3 from Basurto and de Urquijo (2002).

The source reports a mass-resolved swarm observable, reduced mobility, rather
than an elastic differential cross section.  This module preserves that
boundary: it provides a differentiable measured mobility closure over the
digitized support and derived drift/momentum-relaxation scales, but it does
not manufacture an angular collision kernel or a sheath IEAD.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.interpolate import PchipInterpolator

from .network import E_CHARGE_C
from .transport import ATOMIC_MASS_UNIT_KG


CM2_TO_M2 = 1.0e-4
TOWNSEND_V_M2 = 1.0e-21
LOSCHMIDT_NUMBER_DENSITY_M3 = 2.686780111e25
CHF2_PLUS_MASS_AMU = 51.0
BASURTO_MOBILITY_CSV_SHA256 = (
    "96dc38fc862888ae2252be123a2209166ac3a6d35c2e2c0718ab33aead227775"
)
_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BASURTO_MOBILITY_CSV = (
    _ROOT / "data" / "experimental" / "basurto_2002_chf3"
    / "figure1_chf2_chf3_reduced_mobility.csv"
)


def _scalar_or_array(value: np.ndarray) -> float | np.ndarray:
    return float(value) if value.ndim == 0 else value


@dataclass(frozen=True)
class CHF2MobilityState:
    """Measured-swarm closure evaluated at one reduced field and density."""

    reduced_field_Td: float
    reduced_mobility_cm2_V_s: float
    actual_mobility_m2_V_s: float
    drift_speed_m_s: float
    effective_momentum_relaxation_frequency_s_inv: float
    drift_relaxation_length_m: float
    total_neutral_density_m3: float
    source_relative_uncertainty: float
    digitization_reduced_field_relative_bound: float
    digitization_reduced_mobility_relative_bound: float
    supports_measured_swarm_transport: bool = True
    supports_elastic_differential_cross_section: bool = False
    supports_target_sheath_iead: bool = False
    supports_absolute_depth_prediction: bool = False

    def __post_init__(self):
        values = np.asarray([
            self.reduced_field_Td,
            self.reduced_mobility_cm2_V_s,
            self.actual_mobility_m2_V_s,
            self.drift_speed_m_s,
            self.effective_momentum_relaxation_frequency_s_inv,
            self.drift_relaxation_length_m,
            self.total_neutral_density_m3,
            self.source_relative_uncertainty,
            self.digitization_reduced_field_relative_bound,
            self.digitization_reduced_mobility_relative_bound,
        ])
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or self.source_relative_uncertainty >= 1.0
            or self.digitization_reduced_field_relative_bound >= 1.0
            or self.digitization_reduced_mobility_relative_bound >= 1.0
            or not self.supports_measured_swarm_transport
            or self.supports_elastic_differential_cross_section
            or self.supports_target_sheath_iead
            or self.supports_absolute_depth_prediction
        ):
            raise ValueError("invalid CHF2+ mobility state")


@dataclass(frozen=True)
class Basurto2002CHF2CHF3MobilityModel:
    """C1 interpolation of the measured CHF2+--CHF3 reduced mobility."""

    reduced_field_Td: np.ndarray
    reduced_mobility_cm2_V_s: np.ndarray
    digitization_reduced_field_relative_bound: float
    digitization_reduced_mobility_relative_bound: float
    source_relative_uncertainty: float = 0.04
    ion_mass_amu: float = CHF2_PLUS_MASS_AMU
    source: str = "basurto-2002-chf3-ion-mobility"
    provenance: Mapping[str, object] = field(default_factory=lambda: {
        "observable": "mass-resolved CHF2+ reduced mobility in CHF3",
        "measurement_pressure_mTorr": [5.0, 100.0],
        "measurement_temperature_K": [293.0, 310.0],
        "interpolation": "PCHIP in log(E/N), no extrapolation",
        "missing": (
            "elastic differential cross section, angular scattering kernel, "
            "other molecular-ion/neutral pairs, target ion fractions"
        ),
        "coefficient_selected_from_depth_target": None,
    })

    def __post_init__(self):
        field_values = np.asarray(self.reduced_field_Td, dtype=float).copy()
        mobility = np.asarray(
            self.reduced_mobility_cm2_V_s, dtype=float).copy()
        scalars = np.asarray([
            self.digitization_reduced_field_relative_bound,
            self.digitization_reduced_mobility_relative_bound,
            self.source_relative_uncertainty,
            self.ion_mass_amu,
        ])
        if (
            field_values.ndim != 1
            or field_values.size < 8
            or mobility.shape != field_values.shape
            or np.any(~np.isfinite(field_values))
            or np.any(~np.isfinite(mobility))
            or np.any(field_values <= 0.0)
            or np.any(mobility <= 0.0)
            or np.any(np.diff(field_values) <= 0.0)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars <= 0.0)
            or np.any(scalars[:3] >= 1.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid Basurto CHF2+ mobility model")
        field_values.setflags(write=False)
        mobility.setflags(write=False)
        object.__setattr__(self, "reduced_field_Td", field_values)
        object.__setattr__(self, "reduced_mobility_cm2_V_s", mobility)
        object.__setattr__(self, "provenance", MappingProxyType(
            dict(self.provenance)))
        object.__setattr__(self, "_mobility_pchip", PchipInterpolator(
            np.log(field_values), mobility, extrapolate=False))

    @property
    def reduced_field_support_Td(self) -> tuple[float, float]:
        return float(self.reduced_field_Td[0]), float(
            self.reduced_field_Td[-1])

    def _checked_field(self, reduced_field_Td) -> np.ndarray:
        field_values = np.asarray(reduced_field_Td, dtype=float)
        lower, upper = self.reduced_field_support_Td
        if (
            np.any(~np.isfinite(field_values))
            or np.any(field_values < lower)
            or np.any(field_values > upper)
        ):
            raise ValueError("reduced field is outside digitized support")
        return field_values

    def reduced_mobility(self, reduced_field_Td) -> float | np.ndarray:
        field_values = self._checked_field(reduced_field_Td)
        result = np.asarray(self._mobility_pchip(np.log(field_values)))
        return _scalar_or_array(result)

    def reduced_mobility_derivative_cm2_V_s_per_Td(
        self, reduced_field_Td,
    ) -> float | np.ndarray:
        field_values = self._checked_field(reduced_field_Td)
        result = (
            np.asarray(self._mobility_pchip.derivative()(
                np.log(field_values))) / field_values
        )
        return _scalar_or_array(result)

    def evaluate(
        self,
        *,
        reduced_field_Td: float,
        total_neutral_density_m3: float,
    ) -> CHF2MobilityState:
        field_value = float(reduced_field_Td)
        density = float(total_neutral_density_m3)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("total neutral density must be positive")
        reduced_mobility = float(self.reduced_mobility(field_value))
        actual_mobility = (
            reduced_mobility * CM2_TO_M2
            * LOSCHMIDT_NUMBER_DENSITY_M3 / density
        )
        electric_field_V_m = field_value * TOWNSEND_V_M2 * density
        drift_speed = actual_mobility * electric_field_V_m
        ion_mass_kg = self.ion_mass_amu * ATOMIC_MASS_UNIT_KG
        relaxation_frequency = E_CHARGE_C / (
            ion_mass_kg * actual_mobility)
        return CHF2MobilityState(
            reduced_field_Td=field_value,
            reduced_mobility_cm2_V_s=reduced_mobility,
            actual_mobility_m2_V_s=actual_mobility,
            drift_speed_m_s=drift_speed,
            effective_momentum_relaxation_frequency_s_inv=(
                relaxation_frequency),
            drift_relaxation_length_m=drift_speed / relaxation_frequency,
            total_neutral_density_m3=density,
            source_relative_uncertainty=self.source_relative_uncertainty,
            digitization_reduced_field_relative_bound=(
                self.digitization_reduced_field_relative_bound),
            digitization_reduced_mobility_relative_bound=(
                self.digitization_reduced_mobility_relative_bound),
        )


def load_basurto_2002_chf2_chf3_mobility_model(
    path: Path = DEFAULT_BASURTO_MOBILITY_CSV,
) -> Basurto2002CHF2CHF3MobilityModel:
    """Load the checksum-pinned Figure-1 open-circle series."""
    payload = Path(path).read_bytes()
    if sha256(payload).hexdigest() != BASURTO_MOBILITY_CSV_SHA256:
        raise ValueError("Basurto mobility CSV checksum changed")
    records = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    if len(records) < 8:
        raise ValueError("Basurto mobility CSV is incomplete")
    field_values = np.asarray([
        float(row["reduced_field_Td"]) for row in records], dtype=float)
    mobility = np.asarray([
        float(row["reduced_mobility_cm2_V_s"]) for row in records],
        dtype=float,
    )
    field_bounds = {
        float(row["digitization_reduced_field_relative_bound"])
        for row in records
    }
    mobility_bounds = {
        float(row["digitization_reduced_mobility_relative_bound"])
        for row in records
    }
    source_upper = {
        float(row["source_measurement_relative_uncertainty_upper"])
        for row in records
    }
    if len(field_bounds) != 1 or len(mobility_bounds) != 1 or len(
        source_upper) != 1:
        raise ValueError("inconsistent Basurto uncertainty columns")
    return Basurto2002CHF2CHF3MobilityModel(
        reduced_field_Td=field_values,
        reduced_mobility_cm2_V_s=mobility,
        digitization_reduced_field_relative_bound=field_bounds.pop(),
        digitization_reduced_mobility_relative_bound=mobility_bounds.pop(),
        source_relative_uncertainty=source_upper.pop(),
    )
