"""Deterministic power-closed RF-sheath transfer for species-resolved IEADs.

The upstream wafer transfer closes the *mean* bias-energy gain exactly from
``P_bias = A e Gamma_i <Delta E>``.  This module preserves that measured power
closure while replacing the monoenergetic boundary by finite-transit ion
trajectories through a periodic Child-profile sheath.  Entry phase is a fixed
periodic quadrature, never Monte Carlo, and each ion species retains its own
mass-dependent energy distribution.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.optimize import brentq

from ..sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage
from .network import E_CHARGE_C


@dataclass(frozen=True)
class SpeciesResolvedIonEnergyDistribution:
    species: str
    ion_mass_amu: float
    flux_m2_s: float
    energy_eV: np.ndarray
    weight: np.ndarray

    def __post_init__(self):
        energy = np.asarray(self.energy_eV, dtype=float).copy()
        weight = np.asarray(self.weight, dtype=float).copy()
        if (
            not str(self.species).strip()
            or not math.isfinite(self.ion_mass_amu)
            or self.ion_mass_amu <= 0.0
            or not math.isfinite(self.flux_m2_s)
            or self.flux_m2_s < 0.0
            or energy.ndim != 1
            or energy.size < 8
            or weight.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(energy < 0.0)
            or np.any(~np.isfinite(weight))
            or np.any(weight < 0.0)
            or not np.isclose(np.sum(weight), 1.0, rtol=0.0, atol=1.0e-12)
        ):
            raise ValueError("invalid species-resolved ion energy distribution")
        energy.setflags(write=False)
        weight.setflags(write=False)
        object.__setattr__(self, "energy_eV", energy)
        object.__setattr__(self, "weight", weight)

    @property
    def mean_energy_eV(self) -> float:
        return float(np.sum(self.weight * self.energy_eV))

    @property
    def standard_deviation_eV(self) -> float:
        mean = self.mean_energy_eV
        return float(np.sqrt(np.sum(self.weight * (self.energy_eV - mean) ** 2)))

    def probability_inside(self, lower_eV: float, upper_eV: float) -> float:
        lower = float(lower_eV)
        upper = float(upper_eV)
        if (
            not math.isfinite(lower)
            or not math.isfinite(upper)
            or lower < 0.0
            or upper < lower
        ):
            raise ValueError("invalid IEAD support interval")
        return float(np.sum(
            self.weight
            * ((self.energy_eV >= lower) & (self.energy_eV <= upper))))


@dataclass(frozen=True)
class PowerClosedRFSheathProjection:
    distributions: Mapping[str, SpeciesResolvedIonEnergyDistribution]
    target_mean_bias_energy_gain_eV: float
    realized_mean_bias_energy_gain_eV: float
    plasma_potential_eV: float
    electron_temperature_eV: float
    electron_density_m3: float
    bias_dc_component_v: float
    sheath_dc_v: float
    sheath_rf_amplitude_v: float
    frequency_hz: float
    sheath_thickness_m_by_species: Mapping[str, float]
    delivered_bias_power_W: float
    reconstructed_bias_power_W: float
    power_closure_relative_residual: float
    source: str
    evidence_supports_prediction: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        distributions = dict(self.distributions)
        thickness = {
            str(name): float(value)
            for name, value in self.sheath_thickness_m_by_species.items()
        }
        scalars = (
            self.target_mean_bias_energy_gain_eV,
            self.realized_mean_bias_energy_gain_eV,
            self.plasma_potential_eV,
            self.electron_temperature_eV,
            self.electron_density_m3,
            self.bias_dc_component_v,
            self.sheath_dc_v,
            self.sheath_rf_amplitude_v,
            self.frequency_hz,
            self.delivered_bias_power_W,
            self.reconstructed_bias_power_W,
            abs(self.power_closure_relative_residual),
        )
        if (
            not distributions
            or any(
                name != item.species
                or not isinstance(item, SpeciesResolvedIonEnergyDistribution)
                for name, item in distributions.items()
            )
            or set(distributions) != set(thickness)
            or any(not math.isfinite(value) or value <= 0.0
                   for value in thickness.values())
            or any(not math.isfinite(value) or value < 0.0 for value in scalars)
            or self.electron_temperature_eV <= 0.0
            or self.electron_density_m3 <= 0.0
            or self.frequency_hz <= 0.0
            or abs(self.power_closure_relative_residual) > 2.0e-8
            or not str(self.source).strip()
            or self.evidence_supports_prediction
            or self.supports_feature_depth
        ):
            raise ValueError("invalid power-closed RF-sheath projection")
        object.__setattr__(
            self, "distributions", MappingProxyType(distributions))
        object.__setattr__(
            self, "sheath_thickness_m_by_species", MappingProxyType(thickness))


@dataclass(frozen=True)
class DiagnosticConditionedRFSheathTransfer:
    """Turn bias power and a global plasma state into deterministic IEADs.

    ``collapse_fraction=1`` applies full sinusoidal modulation to the inferred
    *bias* component, ``V_s(t)=V_p+V_bias[1+sin(omega t)]``.  Keeping the
    static plasma-potential sheath out of the RF amplitude is required by the
    zero-bias limit: with no delivered bias power the sheath is static rather
    than spuriously oscillating between zero and ``2 V_p``.  The closure is
    intentionally evidence-gated: a frequency or waveform transferred from
    an equipment-class source does not support a formal predictive boundary
    until the actual tool waveform is measured.
    """

    ion_mass_amu: Mapping[str, float]
    electrode_area_m2: float
    plasma_potential_eV: float
    frequency_hz: float
    collapse_fraction: float = 1.0
    phase_count: int = 192
    steps_per_period: int = 256
    steps_per_transit: int = 256
    source: str = "equipment-class RF-sheath sensitivity"

    def __post_init__(self):
        masses = {
            str(name): float(value) for name, value in self.ion_mass_amu.items()
        }
        if (
            not masses
            or any(not name or not math.isfinite(value) or value <= 0.0
                   for name, value in masses.items())
            or not math.isfinite(self.electrode_area_m2)
            or self.electrode_area_m2 <= 0.0
            or not math.isfinite(self.plasma_potential_eV)
            or self.plasma_potential_eV < 0.0
            or not math.isfinite(self.frequency_hz)
            or self.frequency_hz <= 0.0
            or not math.isfinite(self.collapse_fraction)
            or not 0.0 <= self.collapse_fraction <= 1.0
            or int(self.phase_count) < 16
            or int(self.steps_per_period) < 32
            or int(self.steps_per_transit) < 32
            or not str(self.source).strip()
        ):
            raise ValueError("invalid diagnostic-conditioned RF-sheath transfer")
        object.__setattr__(self, "ion_mass_amu", MappingProxyType(masses))

    def _evaluate(
        self,
        bias_dc_v: float,
        *,
        positive_ion_flux_m2_s: Mapping[str, float],
        electron_temperature_eV: float,
        electron_density_m3: float,
    ):
        sheath_dc = self.plasma_potential_eV + float(bias_dc_v)
        amplitude = self.collapse_fraction * float(bias_dc_v)
        waveform = PeriodicSheathVoltage.sinusoidal(
            dc_v=sheath_dc,
            amplitude_v=amplitude,
            frequency_hz=self.frequency_hz,
            source=self.source,
            evidence_kind="assumed",
        )
        phases = (
            2.0 * np.pi
            * (np.arange(int(self.phase_count), dtype=float) + 0.5)
            / int(self.phase_count)
        )
        weight = np.full(phases.size, 1.0 / phases.size)
        distributions = {}
        thickness = {}
        total_flux = float(sum(positive_ion_flux_m2_s.values()))
        mean_impact = 0.0
        for species in sorted(positive_ion_flux_m2_s):
            flux = float(positive_ion_flux_m2_s[species])
            sheath = CollisionlessWaveformSheath(
                waveform=waveform,
                Te_eV=float(electron_temperature_eV),
                ion_mass_amu=self.ion_mass_amu[species],
                density_m3=float(electron_density_m3),
            )
            energy = sheath.ion_impact_energies(
                phases,
                steps_per_period=int(self.steps_per_period),
                steps_per_transit=int(self.steps_per_transit),
            )
            distribution = SpeciesResolvedIonEnergyDistribution(
                species=species,
                ion_mass_amu=self.ion_mass_amu[species],
                flux_m2_s=flux,
                energy_eV=energy,
                weight=weight,
            )
            distributions[species] = distribution
            thickness[species] = sheath.thickness
            mean_impact += flux * distribution.mean_energy_eV / total_flux
        bohm_entry_energy = 0.5 * float(electron_temperature_eV)
        mean_bias_gain = (
            mean_impact - bohm_entry_energy - self.plasma_potential_eV)
        return mean_bias_gain, distributions, thickness, sheath_dc, amplitude

    def predict(
        self,
        *,
        positive_ion_flux_m2_s: Mapping[str, float],
        electron_temperature_eV: float,
        electron_density_m3: float,
        delivered_bias_power_W: float,
    ) -> PowerClosedRFSheathProjection:
        ions = {
            str(name): float(value)
            for name, value in positive_ion_flux_m2_s.items()
        }
        if (
            set(ions) != set(self.ion_mass_amu)
            or any(not math.isfinite(value) or value < 0.0
                   for value in ions.values())
            or sum(ions.values()) <= 0.0
            or not math.isfinite(electron_temperature_eV)
            or electron_temperature_eV <= 0.0
            or not math.isfinite(electron_density_m3)
            or electron_density_m3 <= 0.0
            or not math.isfinite(delivered_bias_power_W)
            or delivered_bias_power_W < 0.0
        ):
            raise ValueError("invalid plasma state supplied to RF-sheath transfer")
        total_flux = float(sum(ions.values()))
        target = (
            float(delivered_bias_power_W)
            / (self.electrode_area_m2 * E_CHARGE_C * total_flux)
        )

        def residual(bias_dc_v):
            realized, *_ = self._evaluate(
                bias_dc_v,
                positive_ion_flux_m2_s=ions,
                electron_temperature_eV=electron_temperature_eV,
                electron_density_m3=electron_density_m3,
            )
            return realized - target

        lower = 0.0
        residual_lower = residual(lower)
        if residual_lower > 1.0e-10:
            raise ValueError(
                "source-driven sheath modulation already exceeds the requested bias power")
        upper = max(2.0 * target + self.plasma_potential_eV + 10.0, 10.0)
        residual_upper = residual(upper)
        expansion = 0
        while residual_upper < 0.0 and expansion < 12:
            upper *= 2.0
            residual_upper = residual(upper)
            expansion += 1
        if residual_upper < 0.0:
            raise RuntimeError("could not bracket RF-sheath power closure")
        bias_dc = float(brentq(
            residual,
            lower,
            upper,
            xtol=1.0e-10,
            rtol=4.0 * np.finfo(float).eps,
            maxiter=48,
        ))
        realized, distributions, thickness, sheath_dc, amplitude = self._evaluate(
            bias_dc,
            positive_ion_flux_m2_s=ions,
            electron_temperature_eV=electron_temperature_eV,
            electron_density_m3=electron_density_m3,
        )
        reconstructed = (
            self.electrode_area_m2 * E_CHARGE_C * total_flux * realized)
        # The trajectory integrator resolves the static Child sheath to finite
        # time-step accuracy.  At the exact zero-power root its residual can be
        # a few 1e-12 eV on either side of zero; do not turn that roundoff into
        # a negative physical power in the immutable result contract.
        if float(delivered_bias_power_W) == 0.0 and abs(realized) < 1.0e-8:
            realized = 0.0
            reconstructed = 0.0
        scale = max(
            abs(float(delivered_bias_power_W)), abs(reconstructed), 1.0)
        power_residual = (
            reconstructed - float(delivered_bias_power_W)) / scale
        return PowerClosedRFSheathProjection(
            distributions=distributions,
            target_mean_bias_energy_gain_eV=target,
            realized_mean_bias_energy_gain_eV=realized,
            plasma_potential_eV=self.plasma_potential_eV,
            electron_temperature_eV=float(electron_temperature_eV),
            electron_density_m3=float(electron_density_m3),
            bias_dc_component_v=bias_dc,
            sheath_dc_v=sheath_dc,
            sheath_rf_amplitude_v=amplitude,
            frequency_hz=self.frequency_hz,
            sheath_thickness_m_by_species=thickness,
            delivered_bias_power_W=float(delivered_bias_power_W),
            reconstructed_bias_power_W=reconstructed,
            power_closure_relative_residual=power_residual,
            source=self.source,
            evidence_supports_prediction=False,
            supports_feature_depth=False,
        )

    def project_from_bias_dc_component(
        self,
        *,
        positive_ion_flux_m2_s: Mapping[str, float],
        electron_temperature_eV: float,
        electron_density_m3: float,
        bias_dc_component_v: float,
    ) -> PowerClosedRFSheathProjection:
        """Project an independently identified bias voltage into an IEAD.

        This is the inverse-facing companion to :meth:`predict`.  It is useful
        when a measured self-bias, voltage waveform, or a separate material
        diagnostic identifies the voltage boundary directly.  The implied
        delivered ion power is reconstructed from exactly the same
        ``A e Gamma_i <Delta E>`` ledger; no etch depth enters this operator.
        """
        ions = {
            str(name): float(value)
            for name, value in positive_ion_flux_m2_s.items()
        }
        bias = float(bias_dc_component_v)
        if (
            set(ions) != set(self.ion_mass_amu)
            or any(not math.isfinite(value) or value < 0.0
                   for value in ions.values())
            or sum(ions.values()) <= 0.0
            or not math.isfinite(electron_temperature_eV)
            or electron_temperature_eV <= 0.0
            or not math.isfinite(electron_density_m3)
            or electron_density_m3 <= 0.0
            or not math.isfinite(bias)
            or bias < 0.0
        ):
            raise ValueError("invalid plasma state supplied to RF-sheath transfer")
        realized, distributions, thickness, sheath_dc, amplitude = self._evaluate(
            bias,
            positive_ion_flux_m2_s=ions,
            electron_temperature_eV=electron_temperature_eV,
            electron_density_m3=electron_density_m3,
        )
        if bias == 0.0 and abs(realized) < 5.0e-2:
            # Static trajectory time-stepping misses the analytic Bohm-plus-
            # plasma-potential energy by only millielectronvolts.
            realized = 0.0
        if realized < -1.0e-8:
            raise RuntimeError("bias projection produced negative delivered ion power")
        realized = max(float(realized), 0.0)
        delivered = (
            self.electrode_area_m2 * E_CHARGE_C * sum(ions.values()) * realized
        )
        return PowerClosedRFSheathProjection(
            distributions=distributions,
            target_mean_bias_energy_gain_eV=realized,
            realized_mean_bias_energy_gain_eV=realized,
            plasma_potential_eV=self.plasma_potential_eV,
            electron_temperature_eV=float(electron_temperature_eV),
            electron_density_m3=float(electron_density_m3),
            bias_dc_component_v=bias,
            sheath_dc_v=sheath_dc,
            sheath_rf_amplitude_v=amplitude,
            frequency_hz=self.frequency_hz,
            sheath_thickness_m_by_species=thickness,
            delivered_bias_power_W=delivered,
            reconstructed_bias_power_W=delivered,
            power_closure_relative_residual=0.0,
            source=self.source,
            evidence_supports_prediction=False,
            supports_feature_depth=False,
        )
