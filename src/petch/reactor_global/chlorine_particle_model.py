"""Open, fixed-pressure chlorine particle balances.

This is a deliberately bounded rung toward a recipe-to-wafer reactor model.
It solves the six-species chlorine particle system at a supplied electron
temperature, with explicit feed, pressure control, exact neutral wall
transport, positive-ion wall loss/return, and quasineutrality. Active-plasma
and neutral-control volumes are distinct: chemistry and plasma-wall exchange
occur only in the active region, while feed, exhaust, and pressure inventory
close over the control volume. The exhaust frequency is solved as the
throttle response needed to maintain the pressure setpoint; it is not inferred
from flow and pressure before chemistry is solved.

No electron power balance is present.  Consequently the solution can
reproduce a source model or condition on a measured electron state, but it
cannot claim predictive knobs-to-flux closure.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np
from scipy.optimize import least_squares

from .chlorine_transport import ChlorineNeutralWallTransport
from .geometry import CylindricalReactor
from .network import RateContext, ReactionNetwork

BOLTZMANN_J_K = 1.380649e-23
STANDARD_CUBIC_CENTIMETER_M3 = 1.0e-6
SECONDS_PER_MINUTE = 60.0

REACTOR_SCALAR_EVIDENCE_KINDS = frozenset({
    "measured",
    "interpolated_measurement",
    "reported_equipment",
    "validated_model",
    "published_model",
    "assumed",
    "sensitivity",
})
_PREDICTIVE_EVIDENCE_KINDS = frozenset({"measured", "validated_model"})
_CHLORINE_SPECIES = frozenset({"e", "Cl2", "Cl", "Cl2+", "Cl+", "Cl-"})
_POSITIVE_ION_SPECIES = frozenset({"Cl2+", "Cl+"})
_SOLVED_DENSITY_ORDER = ("Cl2", "Cl", "Cl2+", "Cl+", "Cl-", "e")
_BALANCE_AND_AUDIT_NAMES = frozenset({
    "Cl2",
    "Cl",
    "Cl2+",
    "Cl+",
    "Cl-",
    "quasineutrality",
    "neutral_pressure",
    "chlorine_atom_inventory",
    "electron_wall_current",
})


@dataclass(frozen=True)
class ReactorScalarInput:
    """One scalar reactor input with units and an evidence chain."""

    value: float
    unit: str
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None

    def __post_init__(self):
        value = float(self.value)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not np.isfinite(value)
            or value <= 0.0
            or not str(self.unit).strip()
            or not str(self.source).strip()
            or self.evidence_kind not in REACTOR_SCALAR_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid sourced reactor scalar")
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "relative_uncertainty", uncertainty)

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_uncertainty is not None
        )


def standard_volume_flow_molecules_s(
    flow_sccm: float,
    *,
    standard_temperature_K: float,
    standard_pressure_Pa: float,
) -> float:
    """Convert SCCM to molecules/s at explicitly supplied standard conditions."""
    values = np.asarray([
        flow_sccm,
        standard_temperature_K,
        standard_pressure_Pa,
    ], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("flow and standard conditions must be positive")
    standard_number_density = (
        standard_pressure_Pa
        / (BOLTZMANN_J_K * standard_temperature_K)
    )
    return float(
        flow_sccm
        * STANDARD_CUBIC_CENTIMETER_M3
        / SECONDS_PER_MINUTE
        * standard_number_density
    )


@dataclass(frozen=True)
class PositiveIonWallTransport:
    """Flux velocities from volume density to the axial and radial walls."""

    axial_flux_velocity_m_s: float
    radial_flux_velocity_m_s: float
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        values = np.asarray([
            self.axial_flux_velocity_m_s,
            self.radial_flux_velocity_m_s,
        ], dtype=float)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
            or self.evidence_kind not in REACTOR_SCALAR_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid positive-ion wall transport")
        object.__setattr__(
            self, "axial_flux_velocity_m_s", float(values[0]))
        object.__setattr__(
            self, "radial_flux_velocity_m_s", float(values[1]))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_uncertainty is not None
        )


@dataclass(frozen=True)
class ChlorineChargedTransportState:
    """Species-resolved positive-ion transport and negative-ion confinement."""

    geometry: CylindricalReactor
    positive_ion_transport: Mapping[str, PositiveIonWallTransport]
    negative_ion_confinement_source: str
    negative_ion_confinement_evidence: str
    negative_ion_confinement_relative_uncertainty: float | None = None

    def __post_init__(self):
        transport = dict(self.positive_ion_transport)
        uncertainty = self.negative_ion_confinement_relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not isinstance(self.geometry, CylindricalReactor)
            or set(transport) != _POSITIVE_ION_SPECIES
            or any(
                not isinstance(item, PositiveIonWallTransport)
                for item in transport.values()
            )
            or not str(self.negative_ion_confinement_source).strip()
            or self.negative_ion_confinement_evidence
            not in REACTOR_SCALAR_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid chlorine charged-transport state")
        object.__setattr__(
            self, "positive_ion_transport", MappingProxyType(transport))
        object.__setattr__(
            self,
            "negative_ion_confinement_relative_uncertainty",
            uncertainty,
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            all(
                item.supports_prediction
                for item in self.positive_ion_transport.values()
            )
            and self.negative_ion_confinement_evidence
            in _PREDICTIVE_EVIDENCE_KINDS
            and self.negative_ion_confinement_relative_uncertainty is not None
        )

    def wall_loss_frequency_s_inv(self, species_name: str) -> float:
        if species_name not in self.positive_ion_transport:
            raise KeyError(f"no charged transport for {species_name}")
        transport = self.positive_ion_transport[species_name]
        axial_area = 2.0 * np.pi * self.geometry.radius_m ** 2
        radial_area = (
            2.0 * np.pi * self.geometry.radius_m * self.geometry.length_m)
        return float(
            (
                transport.axial_flux_velocity_m_s * axial_area
                + transport.radial_flux_velocity_m_s * radial_area
            )
            / self.geometry.volume_m3
        )


@dataclass(frozen=True)
class ChlorineFixedPressureCondition:
    """A pressure-controlled chlorine operating point with supplied ``Te``."""

    condition_id: str
    geometry: CylindricalReactor
    neutral_control_volume: ReactorScalarInput
    pressure: ReactorScalarInput
    gas_temperature: ReactorScalarInput
    electron_temperature: ReactorScalarInput
    chlorine_molecule_feed: ReactorScalarInput
    source_power: ReactorScalarInput

    def __post_init__(self):
        required_units = {
            "neutral control volume": (self.neutral_control_volume, "m3"),
            "pressure": (self.pressure, "Pa"),
            "gas temperature": (self.gas_temperature, "K"),
            "electron temperature": (self.electron_temperature, "eV"),
            "chlorine molecule feed": (
                self.chlorine_molecule_feed, "molecule s^-1"),
            "source power": (self.source_power, "W"),
        }
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
        ):
            raise ValueError("invalid chlorine fixed-pressure condition")
        for name, (state, expected_unit) in required_units.items():
            if not isinstance(state, ReactorScalarInput):
                raise TypeError(f"{name} must be a sourced reactor scalar")
            if state.unit != expected_unit:
                raise ValueError(
                    f"{name} must use the explicit unit {expected_unit}")
        if self.neutral_control_volume.value < self.geometry.volume_m3:
            raise ValueError(
                "neutral control volume cannot be smaller than active plasma"
            )

    @property
    def target_neutral_density_m3(self) -> float:
        return float(
            self.pressure.value
            / (BOLTZMANN_J_K * self.gas_temperature.value)
        )

    @property
    def active_plasma_volume_m3(self) -> float:
        return self.geometry.volume_m3

    @property
    def active_volume_fraction(self) -> float:
        return float(
            self.active_plasma_volume_m3
            / self.neutral_control_volume.value
        )


class ChlorineChargedTransportProvider(Protocol):
    """Condition- and state-dependent charged-wall transport provider."""

    name: str
    version: str

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineChargedTransportState:
        ...


@dataclass(frozen=True)
class FixedChlorineChargedTransportProvider:
    """A declared fixed state for reproduction and sensitivity checks."""

    state: ChlorineChargedTransportState
    name: str = "fixed_chlorine_charged_transport"
    version: str = "1"

    def __post_init__(self):
        if (
            not isinstance(self.state, ChlorineChargedTransportState)
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid fixed chlorine transport provider")

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineChargedTransportState:
        if not isinstance(condition, ChlorineFixedPressureCondition):
            raise TypeError("chlorine fixed-pressure condition is required")
        _positive_mapping(densities_m3, expected_keys=_CHLORINE_SPECIES)
        if self.state.geometry != condition.geometry:
            raise ValueError(
                "charged transport geometry does not match condition")
        return self.state


class ChlorineNeutralWallTransportProvider(Protocol):
    """Condition- and state-dependent neutral-wall transport provider."""

    name: str
    version: str
    neutral_density_basis: str

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineNeutralWallTransport:
        ...


@dataclass(frozen=True)
class FixedChlorineNeutralWallTransportProvider:
    """A declared fixed neutral state for reproduction/sensitivity only."""

    state: ChlorineNeutralWallTransport
    name: str = "fixed_chlorine_neutral_wall_transport"
    version: str = "1"
    neutral_density_basis: str = "condition_target"

    def __post_init__(self):
        if (
            not isinstance(self.state, ChlorineNeutralWallTransport)
            or not str(self.name).strip()
            or not str(self.version).strip()
            or self.neutral_density_basis != "condition_target"
        ):
            raise ValueError("invalid fixed chlorine neutral provider")

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
    ) -> ChlorineNeutralWallTransport:
        if not isinstance(condition, ChlorineFixedPressureCondition):
            raise TypeError("chlorine fixed-pressure condition is required")
        _positive_mapping(densities_m3, expected_keys=_CHLORINE_SPECIES)
        _require_neutral_transport_matches_condition(self.state, condition)
        return self.state


@dataclass(frozen=True)
class ChlorineParticleSolution:
    """Solved particle state and its independent conservation ledgers.

    Volumetric source/loss fields are averaged over the neutral control
    volume; densities and wafer-facing fluxes refer to the active plasma.
    """

    condition_id: str
    active_plasma_volume_m3: float
    neutral_control_volume_m3: float
    active_volume_fraction: float
    densities_m3: Mapping[str, float]
    axial_positive_ion_flux_m2_s: Mapping[str, float]
    positive_ion_wall_loss_m3_s: Mapping[str, float]
    exhaust_loss_frequency_s_inv: float
    chlorine_wall_atom_loss_m3_s: float
    chlorine_wall_molecule_return_m3_s: float
    normalized_balance_residuals: Mapping[str, float]
    chlorine_atom_inventory_residual_m3_s: float
    electron_current_balance_residual_m3_s: float
    solver_evaluations: int
    neutral_transport_supports_prediction: bool
    charged_transport_supports_prediction: bool
    condition_inputs_support_prediction: bool

    def __post_init__(self):
        densities = _positive_mapping(
            self.densities_m3, expected_keys=_CHLORINE_SPECIES)
        fluxes = _nonnegative_mapping(
            self.axial_positive_ion_flux_m2_s,
            expected_keys=_POSITIVE_ION_SPECIES,
        )
        wall_loss = _nonnegative_mapping(
            self.positive_ion_wall_loss_m3_s,
            expected_keys=_POSITIVE_ION_SPECIES,
        )
        residuals = {
            str(name): float(value)
            for name, value in self.normalized_balance_residuals.items()
        }
        scalar_values = np.asarray([
            self.active_plasma_volume_m3,
            self.neutral_control_volume_m3,
            self.active_volume_fraction,
            self.exhaust_loss_frequency_s_inv,
            self.chlorine_wall_atom_loss_m3_s,
            self.chlorine_wall_molecule_return_m3_s,
            self.chlorine_atom_inventory_residual_m3_s,
            self.electron_current_balance_residual_m3_s,
        ], dtype=float)
        if (
            not str(self.condition_id).strip()
            or np.any(~np.isfinite(scalar_values))
            or np.any(scalar_values[:3] <= 0.0)
            or np.any(scalar_values[3:6] < 0.0)
            or self.neutral_control_volume_m3 < self.active_plasma_volume_m3
            or not np.isclose(
                self.active_volume_fraction,
                self.active_plasma_volume_m3
                / self.neutral_control_volume_m3,
                rtol=1.0e-14,
                atol=0.0,
            )
            or self.active_volume_fraction > 1.0
            or not residuals
            or set(residuals) != _BALANCE_AND_AUDIT_NAMES
            or any(not np.isfinite(value) for value in residuals.values())
            or int(self.solver_evaluations) <= 0
            or any(
                not isinstance(value, bool)
                for value in (
                    self.neutral_transport_supports_prediction,
                    self.charged_transport_supports_prediction,
                    self.condition_inputs_support_prediction,
                )
            )
        ):
            raise ValueError("invalid chlorine particle solution")
        object.__setattr__(self, "densities_m3", MappingProxyType(densities))
        object.__setattr__(
            self,
            "axial_positive_ion_flux_m2_s",
            MappingProxyType(fluxes),
        )
        object.__setattr__(
            self,
            "positive_ion_wall_loss_m3_s",
            MappingProxyType(wall_loss),
        )
        object.__setattr__(
            self,
            "normalized_balance_residuals",
            MappingProxyType(residuals),
        )
        object.__setattr__(
            self,
            "exhaust_loss_frequency_s_inv",
            float(self.exhaust_loss_frequency_s_inv),
        )
        for name in (
            "active_plasma_volume_m3",
            "neutral_control_volume_m3",
            "active_volume_fraction",
            "chlorine_wall_atom_loss_m3_s",
            "chlorine_wall_molecule_return_m3_s",
            "chlorine_atom_inventory_residual_m3_s",
            "electron_current_balance_residual_m3_s",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))
        object.__setattr__(
            self, "solver_evaluations", int(self.solver_evaluations))
        if not np.isclose(
            self.chlorine_wall_atom_loss_m3_s,
            2.0 * self.chlorine_wall_molecule_return_m3_s,
            rtol=1.0e-14,
            atol=0.0,
        ):
            raise ValueError("chlorine wall return does not conserve atoms")

    @property
    def maximum_normalized_residual(self) -> float:
        return float(max(
            abs(value) for value in self.normalized_balance_residuals.values()
        ))

    @property
    def chlorine_atom_dissociation_fraction(self) -> float:
        atoms = self.densities_m3["Cl"]
        molecules = self.densities_m3["Cl2"]
        return float(atoms / (atoms + 2.0 * molecules))

    @property
    def total_axial_positive_ion_flux_m2_s(self) -> float:
        return float(sum(self.axial_positive_ion_flux_m2_s.values()))

    @property
    def supports_particle_reproduction(self) -> bool:
        return self.maximum_normalized_residual <= 1.0e-8

    @property
    def supports_prediction(self) -> bool:
        """The fixed-Te particle model has no electron power closure."""
        return False

    @property
    def missing_prediction_closures(self) -> tuple[str, ...]:
        missing = ["electron_power_balance"]
        if not self.condition_inputs_support_prediction:
            missing.append("condition_input_evidence")
        if not self.neutral_transport_supports_prediction:
            missing.append("neutral_transport_evidence")
        if not self.charged_transport_supports_prediction:
            missing.append("charged_transport_evidence")
        return tuple(missing)


class FixedElectronTemperatureChlorineParticleModel:
    """Solve the open chlorine particle ledger at an externally supplied Te."""

    def __init__(self, network: ReactionNetwork):
        if not isinstance(network, ReactionNetwork):
            raise TypeError("chlorine reaction network is required")
        if set(network.species_names) != _CHLORINE_SPECIES:
            raise ValueError("reaction network is not the six-species chlorine deck")
        self.network = network
        self._species_index = {
            name: index for index, name in enumerate(network.species_names)
        }

    def _ledger(
        self,
        log_state: np.ndarray,
        *,
        condition: ChlorineFixedPressureCondition,
        charged_transport_provider: ChlorineChargedTransportProvider,
        neutral_wall_transport_provider: (
            ChlorineNeutralWallTransportProvider),
    ) -> dict[str, object]:
        state = np.exp(np.asarray(log_state, dtype=float))
        densities = dict(zip(_SOLVED_DENSITY_ORDER, state[:6]))
        exhaust_frequency = float(state[6])
        neutral_wall_transport = neutral_wall_transport_provider.predict(
            condition, densities)
        if not isinstance(
            neutral_wall_transport, ChlorineNeutralWallTransport
        ):
            raise TypeError(
                "neutral transport provider returned an invalid state")
        _require_neutral_transport_matches_condition(
            neutral_wall_transport,
            condition,
            expected_neutral_density_m3=(
                densities["Cl2"] + densities["Cl"]
                if neutral_wall_transport_provider.neutral_density_basis
                == "state_total_neutral_particles"
                else condition.target_neutral_density_m3
            ),
        )
        charged_transport = charged_transport_provider.predict(
            condition, densities)
        if not isinstance(charged_transport, ChlorineChargedTransportState):
            raise TypeError(
                "charged transport provider returned an invalid state")
        context = RateContext(
            electron_temperature_eV=condition.electron_temperature.value,
            gas_temperature_K=condition.gas_temperature.value,
        )
        event_rates = self.network.event_rates_m3_s(densities, context)
        active_fraction = condition.active_volume_fraction
        network_source = (
            self.network.stoichiometric_matrix @ event_rates
            * active_fraction
        )
        network_turnover = (
            np.abs(self.network.stoichiometric_matrix) @ event_rates)
        network_turnover *= active_fraction
        source = {
            name: float(network_source[index])
            for name, index in self._species_index.items()
        }
        turnover = {
            name: float(network_turnover[index])
            for name, index in self._species_index.items()
        }
        feed_rate_density = (
            condition.chlorine_molecule_feed.value
            / condition.neutral_control_volume.value
        )
        active_neutral_wall = neutral_wall_transport.evaluate_volume_rates(
            densities["Cl"])
        neutral_wall_atom_loss = (
            active_fraction * active_neutral_wall.chlorine_atom_loss_m3_s
        )
        neutral_wall_molecule_return = (
            active_fraction
            * active_neutral_wall.chlorine_molecule_return_m3_s
        )
        positive_wall_loss = {
            species: (
                active_fraction
                * densities[species]
                * charged_transport.wall_loss_frequency_s_inv(species)
            )
            for species in _POSITIVE_ION_SPECIES
        }

        external_terms = {
            "Cl2": (
                feed_rate_density,
                -exhaust_frequency * densities["Cl2"],
                neutral_wall_molecule_return,
                positive_wall_loss["Cl2+"],
            ),
            "Cl": (
                -exhaust_frequency * densities["Cl"],
                -neutral_wall_atom_loss,
                positive_wall_loss["Cl+"],
            ),
            "Cl2+": (
                -positive_wall_loss["Cl2+"],
            ),
            "Cl+": (
                -positive_wall_loss["Cl+"],
            ),
            "Cl-": (),
        }
        normalized_balances = {
            name: (
                source[name] + sum(values)
            ) / max(
                turnover[name] + sum(abs(value) for value in values),
                1.0,
            )
            for name, values in external_terms.items()
        }
        quasineutrality = (
            densities["Cl2+"]
            + densities["Cl+"]
            - densities["Cl-"]
            - densities["e"]
        )
        normalized_balances["quasineutrality"] = (
            quasineutrality
            / max(
                densities["Cl2+"] + densities["Cl+"],
                densities["Cl-"] + densities["e"],
                1.0,
            )
        )
        neutral_density = densities["Cl2"] + densities["Cl"]
        normalized_balances["neutral_pressure"] = (
            neutral_density - condition.target_neutral_density_m3
        ) / condition.target_neutral_density_m3
        residual = np.asarray([
            normalized_balances[name]
            for name in (
                "Cl2", "Cl", "Cl2+", "Cl+", "Cl-",
                "quasineutrality", "neutral_pressure",
            )
        ])
        chlorine_atom_residual = (
            2.0 * feed_rate_density
            - exhaust_frequency
            * (2.0 * densities["Cl2"] + densities["Cl"])
        )
        electron_wall_loss = sum(positive_wall_loss.values())
        electron_current_residual = source["e"] - electron_wall_loss
        normalized_balances["chlorine_atom_inventory"] = (
            chlorine_atom_residual
            / max(
                2.0 * feed_rate_density,
                exhaust_frequency
                * (2.0 * densities["Cl2"] + densities["Cl"]),
                1.0,
            )
        )
        normalized_balances["electron_wall_current"] = (
            electron_current_residual
            / max(abs(source["e"]), electron_wall_loss, 1.0)
        )
        return {
            "residual": residual,
            "normalized_balances": normalized_balances,
            "densities": densities,
            "exhaust_frequency_s_inv": exhaust_frequency,
            "active_neutral_wall": active_neutral_wall,
            "neutral_wall_atom_loss_m3_s": neutral_wall_atom_loss,
            "neutral_wall_molecule_return_m3_s": (
                neutral_wall_molecule_return),
            "neutral_wall_transport": neutral_wall_transport,
            "positive_wall_loss": positive_wall_loss,
            "charged_transport": charged_transport,
            "chlorine_atom_residual_m3_s": chlorine_atom_residual,
            "electron_current_residual_m3_s": electron_current_residual,
        }

    def solve(
        self,
        condition: ChlorineFixedPressureCondition,
        *,
        charged_transport_provider: ChlorineChargedTransportProvider,
        neutral_wall_transport_provider: (
            ChlorineNeutralWallTransportProvider),
        initial_densities_m3: Mapping[str, float] | None = None,
        initial_exhaust_loss_frequency_s_inv: float | None = None,
        residual_tolerance: float = 1.0e-8,
        maximum_evaluations: int = 3000,
    ) -> ChlorineParticleSolution:
        if not isinstance(condition, ChlorineFixedPressureCondition):
            raise TypeError("chlorine fixed-pressure condition is required")
        if (
            not hasattr(charged_transport_provider, "predict")
            or not str(
                getattr(charged_transport_provider, "name", "")
            ).strip()
            or not str(
                getattr(charged_transport_provider, "version", "")
            ).strip()
        ):
            raise TypeError(
                "versioned chlorine charged-transport provider is required")
        if (
            not hasattr(neutral_wall_transport_provider, "predict")
            or not str(
                getattr(neutral_wall_transport_provider, "name", "")
            ).strip()
            or not str(
                getattr(neutral_wall_transport_provider, "version", "")
            ).strip()
            or getattr(
                neutral_wall_transport_provider,
                "neutral_density_basis",
                None,
            ) not in {
                "condition_target", "state_total_neutral_particles"
            }
        ):
            raise TypeError(
                "versioned chlorine neutral-transport provider is required")

        target_density = condition.target_neutral_density_m3
        if initial_densities_m3 is None:
            charge_seed = min(1.0e16, 1.0e-3 * target_density)
            initial = {
                "Cl2": 0.55 * target_density,
                "Cl": 0.45 * target_density,
                "Cl2+": charge_seed,
                "Cl+": charge_seed,
                "Cl-": charge_seed,
                "e": charge_seed,
            }
        else:
            initial = _positive_mapping(
                initial_densities_m3,
                expected_keys=_CHLORINE_SPECIES,
            )
        if initial_exhaust_loss_frequency_s_inv is None:
            initial_exhaust = (
                condition.chlorine_molecule_feed.value
                / (
                    target_density
                    * condition.neutral_control_volume.value
                )
            )
        else:
            initial_exhaust = float(initial_exhaust_loss_frequency_s_inv)
        if not np.isfinite(initial_exhaust) or initial_exhaust <= 0.0:
            raise ValueError("initial exhaust frequency must be positive")

        initial_vector = np.asarray([
            initial[name] for name in _SOLVED_DENSITY_ORDER
        ] + [initial_exhaust])
        density_floor = max(1.0, 1.0e-18 * target_density)
        lower = np.log([density_floor] * 6 + [1.0e-9])
        upper = np.log([10.0 * target_density] * 6 + [1.0e6])
        result = least_squares(
            lambda state: self._ledger(
                state,
                condition=condition,
                charged_transport_provider=charged_transport_provider,
                neutral_wall_transport_provider=(
                    neutral_wall_transport_provider),
            )["residual"],
            x0=np.clip(np.log(initial_vector), lower, upper),
            bounds=(lower, upper),
            xtol=1.0e-13,
            ftol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=int(maximum_evaluations),
        )
        ledger = self._ledger(
            result.x,
            condition=condition,
            charged_transport_provider=charged_transport_provider,
            neutral_wall_transport_provider=neutral_wall_transport_provider,
        )
        residual = np.asarray(ledger["residual"])
        if (
            not result.success
            or np.any(~np.isfinite(residual))
            or np.max(np.abs(residual)) > float(residual_tolerance)
        ):
            raise RuntimeError(
                "chlorine particle solve failed conservation gate: "
                f"success={result.success}, residual={residual.tolist()}, "
                f"message={result.message}"
            )
        densities = ledger["densities"]
        charged_transport = ledger["charged_transport"]
        neutral_wall_transport = ledger["neutral_wall_transport"]
        ratio = densities["Cl"] / densities["Cl2"]
        neutral_wall_transport.wall_boundary.require_applicable(
            cl_to_cl2_ratio=ratio,
            pressure_Pa=condition.pressure.value,
            icp_power_W=condition.source_power.value,
            gas_temperature_K=condition.gas_temperature.value,
        )
        axial_flux = {
            species: (
                densities[species]
                * charged_transport.positive_ion_transport[
                    species
                ].axial_flux_velocity_m_s
            )
            for species in _POSITIVE_ION_SPECIES
        }
        condition_support = all(
            state.supports_prediction
            for state in (
                condition.pressure,
                condition.gas_temperature,
                condition.electron_temperature,
                condition.chlorine_molecule_feed,
                condition.source_power,
                condition.neutral_control_volume,
            )
        )
        return ChlorineParticleSolution(
            condition_id=condition.condition_id,
            active_plasma_volume_m3=condition.active_plasma_volume_m3,
            neutral_control_volume_m3=condition.neutral_control_volume.value,
            active_volume_fraction=condition.active_volume_fraction,
            densities_m3=densities,
            axial_positive_ion_flux_m2_s=axial_flux,
            positive_ion_wall_loss_m3_s=ledger["positive_wall_loss"],
            exhaust_loss_frequency_s_inv=ledger[
                "exhaust_frequency_s_inv"],
            chlorine_wall_atom_loss_m3_s=(
                ledger["neutral_wall_atom_loss_m3_s"]),
            chlorine_wall_molecule_return_m3_s=(
                ledger["neutral_wall_molecule_return_m3_s"]),
            normalized_balance_residuals=ledger["normalized_balances"],
            chlorine_atom_inventory_residual_m3_s=ledger[
                "chlorine_atom_residual_m3_s"],
            electron_current_balance_residual_m3_s=ledger[
                "electron_current_residual_m3_s"],
            solver_evaluations=int(result.nfev),
            neutral_transport_supports_prediction=(
                neutral_wall_transport.supports_prediction),
            charged_transport_supports_prediction=(
                charged_transport.supports_prediction),
            condition_inputs_support_prediction=condition_support,
        )


def _positive_mapping(
    values: Mapping[str, float],
    *,
    expected_keys: frozenset[str],
) -> dict[str, float]:
    converted = {str(name): float(value) for name, value in values.items()}
    if (
        set(converted) != expected_keys
        or any(not np.isfinite(value) or value <= 0.0
               for value in converted.values())
    ):
        raise ValueError("mapping keys or positive values are invalid")
    return converted


def _nonnegative_mapping(
    values: Mapping[str, float],
    *,
    expected_keys: frozenset[str],
) -> dict[str, float]:
    converted = {str(name): float(value) for name, value in values.items()}
    if (
        set(converted) != expected_keys
        or any(not np.isfinite(value) or value < 0.0
               for value in converted.values())
    ):
        raise ValueError("mapping keys or nonnegative values are invalid")
    return converted


def _require_neutral_transport_matches_condition(
    transport: ChlorineNeutralWallTransport,
    condition: ChlorineFixedPressureCondition,
    *,
    expected_neutral_density_m3: float | None = None,
) -> None:
    if transport.geometry != condition.geometry:
        raise ValueError("neutral transport geometry does not match condition")
    expected_density = (
        condition.target_neutral_density_m3
        if expected_neutral_density_m3 is None
        else float(expected_neutral_density_m3)
    )
    if not np.isfinite(expected_density) or expected_density <= 0.0:
        raise ValueError("expected neutral transport density is invalid")
    if not np.isclose(
        transport.diffusivity.total_neutral_density_m3,
        expected_density,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "neutral transport density does not match current neutral state")
    if not np.isclose(
        transport.diffusivity.gas_temperature_K,
        condition.gas_temperature.value,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "neutral transport temperature does not match condition")
