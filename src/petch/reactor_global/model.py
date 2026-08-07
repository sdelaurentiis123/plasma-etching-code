"""Steady zero-dimensional argon particle and power balances."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np
from scipy.optimize import least_squares

from petch.sheath import bohm_speed

from .argon import ARGON_MASS_AMU, build_lee_lieberman_argon_volume_network
from .geometry import CylindricalReactor, ElectropositiveEdgeFactors
from .network import E_CHARGE_C, RateContext

BOLTZMANN_J_K = 1.380649e-23
PASCAL_PER_MTORR = 0.13332236842105263
_PREDICTIVE_EVIDENCE = {"measured", "validated_model"}


@dataclass(frozen=True)
class ArgonGlobalCondition:
    """One pure-argon global-model operating point."""

    condition_id: str
    absorbed_power_W: float
    pressure_Pa: float
    gas_temperature_K: float
    geometry: CylindricalReactor
    ion_wall_energy_factor_Te: float
    ion_wall_energy_source: str
    ion_wall_energy_evidence: str

    def __post_init__(self):
        values = np.asarray([
            self.absorbed_power_W,
            self.pressure_Pa,
            self.gas_temperature_K,
            self.ion_wall_energy_factor_Te,
        ], dtype=float)
        if (
            not str(self.condition_id).strip()
            or np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not isinstance(self.geometry, CylindricalReactor)
            or not str(self.ion_wall_energy_source).strip()
            or self.ion_wall_energy_evidence not in {
                "measured", "validated_model", "published_range_member",
                "assumed", "sensitivity",
            }
        ):
            raise ValueError("invalid argon global-model condition")
        object.__setattr__(self, "absorbed_power_W", float(self.absorbed_power_W))
        object.__setattr__(self, "pressure_Pa", float(self.pressure_Pa))
        object.__setattr__(
            self, "gas_temperature_K", float(self.gas_temperature_K))
        object.__setattr__(
            self, "ion_wall_energy_factor_Te",
            float(self.ion_wall_energy_factor_Te))

    @property
    def neutral_ground_density_m3(self) -> float:
        return float(
            self.pressure_Pa
            / (BOLTZMANN_J_K * self.gas_temperature_K))


@dataclass(frozen=True)
class ArgonTransportState:
    """Transport quantities needed by the Lee--Lieberman wall closures."""

    ion_mean_free_path_m: float
    ambipolar_diffusion_m2_s: float
    metastable_effective_diffusion_m2_s: float
    source: str
    evidence_kind: str
    provenance: Mapping[str, object] = None

    def __post_init__(self):
        values = np.asarray([
            self.ion_mean_free_path_m,
            self.ambipolar_diffusion_m2_s,
            self.metastable_effective_diffusion_m2_s,
        ], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
            or self.evidence_kind not in {
                "measured", "validated_model", "published_model",
                "assumed", "sensitivity",
            }
        ):
            raise ValueError("invalid argon transport state")
        object.__setattr__(
            self, "ion_mean_free_path_m", float(self.ion_mean_free_path_m))
        object.__setattr__(
            self, "ambipolar_diffusion_m2_s",
            float(self.ambipolar_diffusion_m2_s))
        object.__setattr__(
            self, "metastable_effective_diffusion_m2_s",
            float(self.metastable_effective_diffusion_m2_s))
        object.__setattr__(
            self, "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)))

    @property
    def supports_prediction(self) -> bool:
        return self.evidence_kind in _PREDICTIVE_EVIDENCE


class ArgonTransportProvider(Protocol):
    """Pressure/temperature-dependent transport provider."""

    name: str
    version: str

    def predict(
            self, condition: ArgonGlobalCondition,
            electron_temperature_eV: float) -> ArgonTransportState:
        ...


@dataclass(frozen=True)
class FixedArgonTransportProvider:
    """A fixed transport state for tests or declared sensitivity studies."""

    state: ArgonTransportState
    name: str = "fixed_argon_transport"
    version: str = "1"

    def __post_init__(self):
        if (
            not isinstance(self.state, ArgonTransportState)
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid fixed argon transport provider")

    def predict(
            self, condition: ArgonGlobalCondition,
            electron_temperature_eV: float) -> ArgonTransportState:
        if not isinstance(condition, ArgonGlobalCondition):
            raise TypeError("argon condition is required")
        if (
            not np.isfinite(electron_temperature_eV)
            or electron_temperature_eV <= 0.0
        ):
            raise ValueError("electron temperature must be positive")
        return self.state


@dataclass(frozen=True)
class ArgonGlobalSolution:
    """A solved operating point and its conservation residuals."""

    condition_id: str
    electron_temperature_eV: float
    electron_density_m3: float
    metastable_density_m3: float
    ion_density_m3: float
    ground_density_m3: float
    axial_ion_flux_m2_s: float
    radial_ion_flux_m2_s: float
    edge_factors: ElectropositiveEdgeFactors
    transport: ArgonTransportState
    ion_balance_residual: float
    metastable_balance_residual: float
    power_balance_residual: float
    absorbed_power_W: float
    modeled_power_loss_W: float
    solver_evaluations: int
    transport_provider: str
    transport_provider_version: str
    ion_wall_energy_evidence: str

    @property
    def maximum_normalized_residual(self) -> float:
        return float(max(
            abs(self.ion_balance_residual),
            abs(self.metastable_balance_residual),
            abs(self.power_balance_residual),
        ))

    @property
    def supports_prediction(self) -> bool:
        return (
            self.transport.supports_prediction
            and self.ion_wall_energy_evidence in _PREDICTIVE_EVIDENCE
            and self.maximum_normalized_residual <= 1.0e-8
        )


class LeeLiebermanArgonGlobalModel:
    """Pure-Ar steady global model using the sourced five-reaction deck."""

    def __init__(self, transport_provider: ArgonTransportProvider):
        if (
            not hasattr(transport_provider, "predict")
            or not str(getattr(transport_provider, "name", "")).strip()
            or not str(getattr(transport_provider, "version", "")).strip()
        ):
            raise ValueError("a versioned argon transport provider is required")
        self.transport_provider = transport_provider
        self.network = build_lee_lieberman_argon_volume_network()
        self._reaction_index = {
            reaction.name: index
            for index, reaction in enumerate(self.network.reactions)
        }

    def _ledger(
            self, log_state: np.ndarray,
            condition: ArgonGlobalCondition) -> dict[str, object]:
        electron_density, metastable_density, electron_temperature = np.exp(
            np.asarray(log_state, dtype=float))
        ground_density = condition.neutral_ground_density_m3
        ion_density = electron_density
        context = RateContext(
            electron_temperature_eV=electron_temperature,
            gas_temperature_K=condition.gas_temperature_K,
        )
        densities = {
            "e": electron_density,
            "Ar": ground_density,
            "Ar*": metastable_density,
            "Ar+": ion_density,
        }
        rates = self.network.event_rates_m3_s(densities, context)
        rate = {
            name: rates[index]
            for name, index in self._reaction_index.items()
        }
        transport = self.transport_provider.predict(
            condition, electron_temperature)
        speed = bohm_speed(electron_temperature, ARGON_MASS_AMU)
        edge = condition.geometry.electropositive_edge_factors(
            ion_mean_free_path_m=transport.ion_mean_free_path_m,
            bohm_speed_m_s=speed,
            ambipolar_diffusion_m2_s=transport.ambipolar_diffusion_m2_s,
        )
        effective_area = condition.geometry.effective_loss_area_m2(edge)
        wall_loss_frequency = (
            speed * effective_area / condition.geometry.volume_m3)
        ion_wall_loss = ion_density * wall_loss_frequency
        metastable_wall_frequency = (
            transport.metastable_effective_diffusion_m2_s
            / condition.geometry.diffusion_length_m ** 2)
        metastable_wall_loss = (
            metastable_density * metastable_wall_frequency)

        ion_production = (
            rate["e_Ar_ground_ionization"]
            + rate["e_Ar_metastable_step_ionization"]
            + rate["Ar_metastable_pooling_associative_ionization"]
        )
        metastable_production = rate["e_Ar_metastable_excitation"]
        metastable_volume_loss = (
            rate["e_Ar_metastable_step_ionization"]
            + rate["e_Ar_metastable_superelastic_quench"]
            + 2.0 * rate["Ar_metastable_pooling_associative_ionization"]
        )
        volume_power_density = (
            self.network.electron_power_loss_density_W_m3(
                densities, context))
        wall_energy_eV = (
            condition.ion_wall_energy_factor_Te + 2.0
        ) * electron_temperature
        wall_power_density = (
            E_CHARGE_C * ion_wall_loss * wall_energy_eV)
        modeled_power_density = volume_power_density + wall_power_density
        absorbed_power_density = (
            condition.absorbed_power_W / condition.geometry.volume_m3)
        residual = np.array([
            (ion_production - ion_wall_loss)
            / max(ion_production, ion_wall_loss, 1.0),
            (metastable_production
             - metastable_volume_loss
             - metastable_wall_loss)
            / max(
                metastable_production,
                metastable_volume_loss + metastable_wall_loss,
                1.0,
            ),
            (absorbed_power_density - modeled_power_density)
            / max(absorbed_power_density, abs(modeled_power_density), 1.0),
        ])
        return {
            "residual": residual,
            "electron_density_m3": electron_density,
            "metastable_density_m3": metastable_density,
            "electron_temperature_eV": electron_temperature,
            "ground_density_m3": ground_density,
            "transport": transport,
            "edge_factors": edge,
            "bohm_speed_m_s": speed,
            "modeled_power_density_W_m3": modeled_power_density,
        }

    def residuals(
            self, *, electron_temperature_eV: float,
            electron_density_m3: float, metastable_density_m3: float,
            condition: ArgonGlobalCondition) -> np.ndarray:
        values = np.asarray([
            electron_density_m3,
            metastable_density_m3,
            electron_temperature_eV,
        ], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("candidate densities and temperature must be positive")
        return np.asarray(
            self._ledger(np.log(values), condition)["residual"]).copy()

    def solve(
            self, condition: ArgonGlobalCondition, *,
            initial_electron_temperature_eV: float = 3.0,
            initial_electron_density_m3: float = 1.0e17,
            initial_metastable_density_m3: float = 1.0e15,
            residual_tolerance: float = 1.0e-8,
            maximum_evaluations: int = 1000) -> ArgonGlobalSolution:
        if not isinstance(condition, ArgonGlobalCondition):
            raise TypeError("argon condition is required")
        initial = np.asarray([
            initial_electron_density_m3,
            initial_metastable_density_m3,
            initial_electron_temperature_eV,
        ], dtype=float)
        if np.any(~np.isfinite(initial)) or np.any(initial <= 0.0):
            raise ValueError("initial state must be positive and finite")
        ground = condition.neutral_ground_density_m3
        lower = np.log([1.0e6, 1.0e6, 0.1])
        upper = np.log([
            max(ground, 1.0e7),
            max(ground, 1.0e7),
            100.0,
        ])
        result = least_squares(
            lambda state: self._ledger(state, condition)["residual"],
            x0=np.clip(np.log(initial), lower, upper),
            bounds=(lower, upper),
            xtol=1.0e-13,
            ftol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=int(maximum_evaluations),
        )
        ledger = self._ledger(result.x, condition)
        residual = np.asarray(ledger["residual"])
        if (
            not result.success
            or np.any(~np.isfinite(residual))
            or np.max(np.abs(residual)) > float(residual_tolerance)
        ):
            raise RuntimeError(
                "argon global solve failed conservation gate: "
                f"success={result.success}, residual={residual.tolist()}, "
                f"message={result.message}")
        electron_density = float(ledger["electron_density_m3"])
        edge = ledger["edge_factors"]
        speed = float(ledger["bohm_speed_m_s"])
        modeled_power_W = float(
            ledger["modeled_power_density_W_m3"]
            * condition.geometry.volume_m3)
        return ArgonGlobalSolution(
            condition_id=condition.condition_id,
            electron_temperature_eV=float(
                ledger["electron_temperature_eV"]),
            electron_density_m3=electron_density,
            metastable_density_m3=float(ledger["metastable_density_m3"]),
            ion_density_m3=electron_density,
            ground_density_m3=float(ledger["ground_density_m3"]),
            axial_ion_flux_m2_s=float(edge.axial * electron_density * speed),
            radial_ion_flux_m2_s=float(edge.radial * electron_density * speed),
            edge_factors=edge,
            transport=ledger["transport"],
            ion_balance_residual=float(residual[0]),
            metastable_balance_residual=float(residual[1]),
            power_balance_residual=float(residual[2]),
            absorbed_power_W=condition.absorbed_power_W,
            modeled_power_loss_W=modeled_power_W,
            solver_evaluations=int(result.nfev),
            transport_provider=str(self.transport_provider.name),
            transport_provider_version=str(self.transport_provider.version),
            ion_wall_energy_evidence=condition.ion_wall_energy_evidence,
        )
