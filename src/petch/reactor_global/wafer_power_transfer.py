"""Diagnostic-conditioned chlorine transfer from a global state to a wafer.

The first implementation follows the dimensional estimates Mahorowala used
for the Lam TCP 9400SE: an independently supplied center ion current fixes the
global-model axial-flux normalization, an independently supplied center
Cl/ion ratio fixes the radical transmission, and bottom-electrode power closes
the mean ion sheath-energy gain through ``P = A e Gamma_i DeltaV``.

No etched rate or depth appears in this transfer.  It is deterministic and
analytic, but it remains a facility-conditioned source-estimate mode until the
reference current, radical flux, bias absorption, and plasma potential are
measured with uncertainties on the actual reactor.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

from .network import E_CHARGE_C


BOLTZMANN_J_K = 1.380649e-23
ATOMIC_MASS_KG = 1.66053906660e-27
CHLORINE_ATOM_MASS_AMU = 35.45


def isotropic_thermal_particle_flux_m2_s(
    number_density_m3: float,
    temperature_K: float,
    *,
    mass_amu: float = CHLORINE_ATOM_MASS_AMU,
) -> float:
    """Return the one-sided Maxwellian wall flux ``n vbar / 4``."""
    density = float(number_density_m3)
    temperature = float(temperature_K)
    mass = float(mass_amu)
    if (
        not math.isfinite(density)
        or density < 0.0
        or not math.isfinite(temperature)
        or temperature <= 0.0
        or not math.isfinite(mass)
        or mass <= 0.0
    ):
        raise ValueError("invalid thermal particle-flux state")
    mean_speed = math.sqrt(
        8.0 * BOLTZMANN_J_K * temperature
        / (math.pi * mass * ATOMIC_MASS_KG)
    )
    return float(0.25 * density * mean_speed)


@dataclass(frozen=True)
class ChlorineWaferTransferEvidence:
    reference_total_ion_flux: str
    reference_neutral_to_ion_ratio: str
    electrode_area: str
    bias_power_coupling: str
    plasma_potential: str
    equipment_transfer: str
    reference_facts_measured: bool = False
    bias_coupling_measured: bool = False
    plasma_potential_measured: bool = False

    def __post_init__(self):
        for name in (
            "reference_total_ion_flux",
            "reference_neutral_to_ion_ratio",
            "electrode_area",
            "bias_power_coupling",
            "plasma_potential",
            "equipment_transfer",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError("wafer-transfer evidence strings cannot be empty")

    @property
    def supports_prediction(self) -> bool:
        return bool(
            self.reference_facts_measured
            and self.bias_coupling_measured
            and self.plasma_potential_measured
        )


@dataclass(frozen=True)
class ChlorineWaferBoundaryProjection:
    positive_ion_flux_m2_s: Mapping[str, float]
    atomic_chlorine_flux_m2_s: float
    mean_sheath_energy_gain_eV: float
    mean_impact_energy_eV: float
    ion_flux_scale: float
    radical_flux_scale: float
    bias_power_delivered_to_ions_W: float
    reconstructed_ion_power_W: float
    power_closure_relative_residual: float
    evidence_supports_prediction: bool
    supports_feature_depth: bool = False

    def __post_init__(self):
        ions = {
            str(name): float(value)
            for name, value in self.positive_ion_flux_m2_s.items()
        }
        scalars = (
            self.atomic_chlorine_flux_m2_s,
            self.mean_sheath_energy_gain_eV,
            self.mean_impact_energy_eV,
            self.ion_flux_scale,
            self.radical_flux_scale,
            self.bias_power_delivered_to_ions_W,
            self.reconstructed_ion_power_W,
            abs(self.power_closure_relative_residual),
        )
        if (
            not ions
            or any(not name or not math.isfinite(value) or value < 0.0
                   for name, value in ions.items())
            or any(not math.isfinite(value) or value < 0.0 for value in scalars)
            or sum(ions.values()) <= 0.0
            or self.mean_impact_energy_eV < self.mean_sheath_energy_gain_eV
            or abs(self.power_closure_relative_residual) > 1.0e-12
            or self.supports_feature_depth
        ):
            raise ValueError("invalid chlorine wafer-boundary projection")
        object.__setattr__(
            self, "positive_ion_flux_m2_s", MappingProxyType(ions))

    @property
    def total_positive_ion_flux_m2_s(self) -> float:
        return float(sum(self.positive_ion_flux_m2_s.values()))

    @property
    def neutral_to_total_ion_flux_ratio(self) -> float:
        return float(
            self.atomic_chlorine_flux_m2_s
            / self.total_positive_ion_flux_m2_s)


@dataclass(frozen=True)
class DiagnosticConditionedChlorineWaferTransfer:
    """Freeze two non-etch diagnostics, then transfer other reactor states."""

    reference_model_positive_ion_flux_m2_s: Mapping[str, float]
    reference_model_atomic_chlorine_density_m3: float
    reference_gas_temperature_K: float
    reference_wafer_total_ion_flux_m2_s: float
    reference_wafer_neutral_to_ion_flux_ratio: float
    electrode_area_m2: float
    bias_power_to_ion_fraction: float
    plasma_potential_eV: float
    evidence: ChlorineWaferTransferEvidence

    def __post_init__(self):
        ions = {
            str(name): float(value)
            for name, value in self.reference_model_positive_ion_flux_m2_s.items()
        }
        values = (
            self.reference_model_atomic_chlorine_density_m3,
            self.reference_gas_temperature_K,
            self.reference_wafer_total_ion_flux_m2_s,
            self.reference_wafer_neutral_to_ion_flux_ratio,
            self.electrode_area_m2,
            self.bias_power_to_ion_fraction,
            self.plasma_potential_eV,
        )
        if (
            not ions
            or any(not name or not math.isfinite(value) or value < 0.0
                   for name, value in ions.items())
            or sum(ions.values()) <= 0.0
            or any(not math.isfinite(value) for value in values)
            or self.reference_model_atomic_chlorine_density_m3 < 0.0
            or self.reference_gas_temperature_K <= 0.0
            or self.reference_wafer_total_ion_flux_m2_s <= 0.0
            or self.reference_wafer_neutral_to_ion_flux_ratio < 0.0
            or self.electrode_area_m2 <= 0.0
            or not 0.0 <= self.bias_power_to_ion_fraction <= 1.0
            or self.plasma_potential_eV < 0.0
            or not isinstance(self.evidence, ChlorineWaferTransferEvidence)
        ):
            raise ValueError("invalid diagnostic-conditioned wafer transfer")
        raw_neutral = isotropic_thermal_particle_flux_m2_s(
            self.reference_model_atomic_chlorine_density_m3,
            self.reference_gas_temperature_K,
        )
        if raw_neutral <= 0.0:
            raise ValueError("reference model atomic-chlorine flux is zero")
        object.__setattr__(
            self,
            "reference_model_positive_ion_flux_m2_s",
            MappingProxyType(ions),
        )

    @property
    def ion_flux_scale(self) -> float:
        return float(
            self.reference_wafer_total_ion_flux_m2_s
            / sum(self.reference_model_positive_ion_flux_m2_s.values())
        )

    @property
    def radical_flux_scale(self) -> float:
        target = (
            self.reference_wafer_neutral_to_ion_flux_ratio
            * self.reference_wafer_total_ion_flux_m2_s
        )
        raw = isotropic_thermal_particle_flux_m2_s(
            self.reference_model_atomic_chlorine_density_m3,
            self.reference_gas_temperature_K,
        )
        return float(target / raw)

    def predict(
        self,
        *,
        model_positive_ion_flux_m2_s: Mapping[str, float],
        model_atomic_chlorine_density_m3: float,
        gas_temperature_K: float,
        applied_bias_power_W: float,
    ) -> ChlorineWaferBoundaryProjection:
        model_ions = {
            str(name): float(value)
            for name, value in model_positive_ion_flux_m2_s.items()
        }
        bias_power = float(applied_bias_power_W)
        if (
            set(model_ions) != set(self.reference_model_positive_ion_flux_m2_s)
            or any(not math.isfinite(value) or value < 0.0
                   for value in model_ions.values())
            or not math.isfinite(bias_power)
            or bias_power < 0.0
        ):
            raise ValueError("invalid state supplied to wafer transfer")
        ions = {
            name: self.ion_flux_scale * value
            for name, value in model_ions.items()
        }
        total_ion_flux = sum(ions.values())
        if total_ion_flux <= 0.0:
            raise ValueError("wafer ion flux is zero")
        neutral_flux = (
            self.radical_flux_scale
            * isotropic_thermal_particle_flux_m2_s(
                model_atomic_chlorine_density_m3,
                gas_temperature_K,
            )
        )
        delivered_power = self.bias_power_to_ion_fraction * bias_power
        sheath_gain = (
            delivered_power
            / (self.electrode_area_m2 * E_CHARGE_C * total_ion_flux)
        )
        reconstructed_power = (
            self.electrode_area_m2
            * E_CHARGE_C
            * total_ion_flux
            * sheath_gain
        )
        scale = max(abs(delivered_power), abs(reconstructed_power), 1.0)
        residual = (reconstructed_power - delivered_power) / scale
        return ChlorineWaferBoundaryProjection(
            positive_ion_flux_m2_s=ions,
            atomic_chlorine_flux_m2_s=neutral_flux,
            mean_sheath_energy_gain_eV=sheath_gain,
            mean_impact_energy_eV=sheath_gain + self.plasma_potential_eV,
            ion_flux_scale=self.ion_flux_scale,
            radical_flux_scale=self.radical_flux_scale,
            bias_power_delivered_to_ions_W=delivered_power,
            reconstructed_ion_power_W=reconstructed_power,
            power_closure_relative_residual=residual,
            evidence_supports_prediction=self.evidence.supports_prediction,
            supports_feature_depth=False,
        )
