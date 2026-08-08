"""Absorbed-power chlorine global model driven by a deterministic EEPF.

This is the reference coupling rung between collision kinetics and the
existing pressure, wall, and charged-transport ledgers.  It intentionally
requires absorbed plasma power as a separate evidence object and keeps the
equipment setpoint outside the electron power balance.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping, Protocol

import numpy as np
from scipy.optimize import least_squares

from .chlorine_particle_model import (
    BOLTZMANN_J_K,
    ChlorineChargedTransportProvider,
    ChlorineChargedTransportState,
    ChlorineFixedPressureCondition,
    ChlorineNeutralWallTransportProvider,
    ReactorScalarInput,
)
from .electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionChemistryState,
)
from .electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    TwoTermBoltzmannCondition,
    TwoTermBoltzmannSolution,
)
from .network import E_CHARGE_C, RateContext, ReactionNetwork
from .power import AbsorbedPowerEstimate


_SPECIES = frozenset({"e", "Cl2", "Cl", "Cl2+", "Cl+", "Cl-"})
_POSITIVE_IONS = frozenset({"Cl2+", "Cl+"})
_STATE_ORDER = ("Cl2", "Cl", "Cl2+", "Cl+", "Cl-", "e")


def _mapping(values, expected, *, positive: bool) -> MappingProxyType:
    converted = {str(name): float(value) for name, value in values.items()}
    invalid_value = (
        (lambda value: value <= 0.0)
        if positive
        else (lambda value: value < 0.0)
    )
    if (
        set(converted) != set(expected)
        or any(not math.isfinite(value) or invalid_value(value)
               for value in converted.values())
    ):
        raise ValueError("invalid chlorine EEPF model mapping")
    return MappingProxyType(converted)


@dataclass(frozen=True)
class EEDFChlorineCondition:
    """Pressure-controlled knobs plus an explicit absorbed-power boundary."""

    condition_id: str
    geometry: object
    neutral_control_volume: ReactorScalarInput
    pressure: ReactorScalarInput
    gas_temperature: ReactorScalarInput
    chlorine_molecule_feed: ReactorScalarInput
    source_power: ReactorScalarInput
    absorbed_power: AbsorbedPowerEstimate
    reduced_field_bounds_Td: tuple[float, float]
    source_frequency: ReactorScalarInput | None = None

    def __post_init__(self):
        from .geometry import CylindricalReactor

        required_units = {
            "neutral control volume": (self.neutral_control_volume, "m3"),
            "pressure": (self.pressure, "Pa"),
            "gas temperature": (self.gas_temperature, "K"),
            "chlorine molecule feed": (
                self.chlorine_molecule_feed, "molecule s^-1"),
            "source power": (self.source_power, "W"),
        }
        bounds = tuple(float(value) for value in self.reduced_field_bounds_Td)
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
            or not isinstance(self.absorbed_power, AbsorbedPowerEstimate)
            or self.absorbed_power.point_W is None
            or len(bounds) != 2
            or not all(math.isfinite(value) and value > 0.0 for value in bounds)
            or bounds[1] <= bounds[0]
        ):
            raise ValueError("invalid EEPF chlorine condition")
        if self.source_frequency is not None and (
            not isinstance(self.source_frequency, ReactorScalarInput)
            or self.source_frequency.unit != "Hz"
        ):
            raise ValueError("source frequency must use the explicit unit Hz")
        for name, (state, unit) in required_units.items():
            if not isinstance(state, ReactorScalarInput):
                raise TypeError(f"{name} must be a sourced reactor scalar")
            if state.unit != unit:
                raise ValueError(f"{name} must use the explicit unit {unit}")
        if self.neutral_control_volume.value < self.geometry.volume_m3:
            raise ValueError("neutral control volume is smaller than plasma")
        object.__setattr__(self, "reduced_field_bounds_Td", bounds)

    @property
    def target_neutral_density_m3(self) -> float:
        return float(
            self.pressure.value
            / (BOLTZMANN_J_K * self.gas_temperature.value)
        )

    @property
    def active_volume_fraction(self) -> float:
        return float(
            self.geometry.volume_m3 / self.neutral_control_volume.value)

    @property
    def absorbed_power_density_W_m3(self) -> float:
        return float(
            self.absorbed_power.require_point_W()
            / self.neutral_control_volume.value
        )

    @property
    def angular_field_frequency_over_density_m3_s(self) -> float:
        if self.source_frequency is None:
            return 0.0
        return float(
            2.0 * math.pi * self.source_frequency.value
            / self.target_neutral_density_m3
        )

    def transport_condition(
        self,
        mean_energy_equivalent_temperature_eV: float,
    ) -> ChlorineFixedPressureCondition:
        return ChlorineFixedPressureCondition(
            condition_id=self.condition_id,
            geometry=self.geometry,
            neutral_control_volume=self.neutral_control_volume,
            pressure=self.pressure,
            gas_temperature=self.gas_temperature,
            electron_temperature=ReactorScalarInput(
                value=mean_energy_equivalent_temperature_eV,
                unit="eV",
                source=(
                    "two-term EEPF mean-energy equivalent used only by "
                    "legacy charged-transport closure"
                ),
                evidence_kind="published_model",
                relative_uncertainty=None,
            ),
            chlorine_molecule_feed=self.chlorine_molecule_feed,
            source_power=self.source_power,
        )


@dataclass(frozen=True)
class PositiveIonWallEnergyState:
    """Positive-ion kinetic/electrostatic energy charged to electron power."""

    energy_eV_per_lost_ion: Mapping[str, float]
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None

    def __post_init__(self):
        energies = _mapping(
            self.energy_eV_per_lost_ion, _POSITIVE_IONS, positive=False)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not str(self.source).strip()
            or not str(self.evidence_kind).strip()
            or (
                uncertainty is not None
                and (not math.isfinite(uncertainty) or not 0.0 <= uncertainty < 1.0)
            )
        ):
            raise ValueError("invalid positive-ion wall energy state")
        object.__setattr__(self, "energy_eV_per_lost_ion", energies)
        object.__setattr__(self, "relative_uncertainty", uncertainty)

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in {"measured", "validated_model"}
            and self.relative_uncertainty is not None
        )


class PositiveIonWallEnergyProvider(Protocol):
    name: str
    version: str

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
        charged_transport: ChlorineChargedTransportState,
        electron_solution: TwoTermBoltzmannSolution,
    ) -> PositiveIonWallEnergyState:
        ...


@dataclass(frozen=True)
class FixedPositiveIonWallEnergyProvider:
    state: PositiveIonWallEnergyState
    name: str = "fixed_positive_ion_wall_energy"
    version: str = "1"

    def predict(
        self,
        condition: ChlorineFixedPressureCondition,
        densities_m3: Mapping[str, float],
        charged_transport: ChlorineChargedTransportState,
        electron_solution: TwoTermBoltzmannSolution,
    ) -> PositiveIonWallEnergyState:
        if not isinstance(self.state, PositiveIonWallEnergyState):
            raise TypeError("a positive-ion wall energy state is required")
        return self.state


@dataclass(frozen=True)
class EEDFChlorineSolution:
    condition_id: str
    densities_m3: Mapping[str, float]
    reduced_electric_field_Td: float
    mean_electron_energy_eV: float
    axial_positive_ion_flux_m2_s: Mapping[str, float]
    exhaust_loss_frequency_s_inv: float
    absorbed_power_density_W_m3: float
    collisional_power_density_W_m3: float
    charged_wall_power_density_W_m3: float
    modeled_power_density_W_m3: float
    normalized_residuals: Mapping[str, float]
    electron_solution: TwoTermBoltzmannSolution
    collision_chemistry_state: ElectronCollisionChemistryState
    solver_evaluations: int
    absorbed_power_supports_prediction: bool
    wall_energy_supports_prediction: bool
    supports_implicit_differentiation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        densities = _mapping(self.densities_m3, _SPECIES, positive=True)
        fluxes = _mapping(
            self.axial_positive_ion_flux_m2_s,
            _POSITIVE_IONS,
            positive=False,
        )
        residuals = {
            str(name): float(value)
            for name, value in self.normalized_residuals.items()
        }
        scalars = np.asarray((
            self.reduced_electric_field_Td,
            self.mean_electron_energy_eV,
            self.exhaust_loss_frequency_s_inv,
            self.absorbed_power_density_W_m3,
            self.collisional_power_density_W_m3,
            self.charged_wall_power_density_W_m3,
            self.modeled_power_density_W_m3,
        ))
        if (
            np.any(~np.isfinite(scalars))
            or np.any(scalars[:4] <= 0.0)
            or np.any(scalars[4:] < 0.0)
            or not residuals
            or any(not math.isfinite(value) for value in residuals.values())
            or int(self.solver_evaluations) <= 0
        ):
            raise ValueError("invalid EEPF chlorine solution")
        object.__setattr__(self, "densities_m3", densities)
        object.__setattr__(self, "axial_positive_ion_flux_m2_s", fluxes)
        object.__setattr__(self, "normalized_residuals", MappingProxyType(residuals))
        object.__setattr__(self, "solver_evaluations", int(self.solver_evaluations))
        for name, value in zip((
            "reduced_electric_field_Td",
            "mean_electron_energy_eV",
            "exhaust_loss_frequency_s_inv",
            "absorbed_power_density_W_m3",
            "collisional_power_density_W_m3",
            "charged_wall_power_density_W_m3",
            "modeled_power_density_W_m3",
        ), scalars):
            object.__setattr__(self, name, float(value))

    @property
    def maximum_normalized_residual(self) -> float:
        return max(abs(value) for value in self.normalized_residuals.values())


class EEDFChlorineAbsorbedPowerModel:
    """Solve chlorine particles, EEPF, reduced field, and power together."""

    def __init__(
        self,
        electron_solver: DeterministicTwoTermBoltzmannSolver,
        collision_chemistry: ElectronCollisionChemistry,
        heavy_reaction_network: ReactionNetwork,
    ):
        if not isinstance(
            electron_solver, DeterministicTwoTermBoltzmannSolver
        ):
            raise TypeError("a deterministic electron solver is required")
        if not isinstance(collision_chemistry, ElectronCollisionChemistry):
            raise TypeError("electron collision chemistry is required")
        if not isinstance(heavy_reaction_network, ReactionNetwork):
            raise TypeError("a heavy reaction network is required")
        if set(heavy_reaction_network.species_names) != _SPECIES:
            raise ValueError("heavy reaction network species are not chlorine")
        for reaction in heavy_reaction_network.reactions:
            if (
                "e" in reaction.reactants
                or "e" in reaction.products
                or "e" in reaction.kinetic_orders
            ):
                raise ValueError(
                    "supplemental network must contain heavy reactions only"
                )
        if electron_solver.collision_deck is not collision_chemistry.collision_deck:
            raise ValueError("electron solver and chemistry must share one deck")
        self.electron_solver = electron_solver
        self.collision_chemistry = collision_chemistry
        self.heavy_reaction_network = heavy_reaction_network
        self._species_index = {
            name: index
            for index, name in enumerate(heavy_reaction_network.species_names)
        }

    def _ledger(
        self,
        log_state: np.ndarray,
        *,
        condition: EEDFChlorineCondition,
        charged_transport_provider: ChlorineChargedTransportProvider,
        neutral_wall_transport_provider: ChlorineNeutralWallTransportProvider,
        wall_energy_provider: PositiveIonWallEnergyProvider,
        maximum_tail_population_fraction: float,
        electron_cache: dict[tuple[float, ...], TwoTermBoltzmannSolution] | None = None,
    ) -> dict[str, object]:
        state = np.exp(np.asarray(log_state, dtype=float))
        densities = dict(zip(_STATE_ORDER, state[:6]))
        exhaust_frequency = float(state[6])
        reduced_field_Td = float(state[7])
        targets = self.electron_solver.collision_deck.targets
        total_target_density = sum(densities[name] for name in targets)
        target_fractions = {
            name: densities[name] / total_target_density for name in targets
        }
        electron_condition = TwoTermBoltzmannCondition(
            reduced_electric_field_Td=reduced_field_Td,
            gas_temperature_K=condition.gas_temperature.value,
            target_mole_fractions=target_fractions,
            growth_model="temporal_growth",
            angular_field_frequency_over_density_m3_s=(
                condition.angular_field_frequency_over_density_m3_s),
        )
        electron_key = (
            reduced_field_Td,
            condition.gas_temperature.value,
            *(target_fractions[name] for name in targets),
            maximum_tail_population_fraction,
            condition.angular_field_frequency_over_density_m3_s,
        )
        electron_solution = (
            None if electron_cache is None else electron_cache.get(electron_key)
        )
        if electron_solution is None:
            electron_solution = self.electron_solver.solve(
                electron_condition,
                relative_tolerance=1.0e-8,
                maximum_iterations=200,
                maximum_tail_population_fraction=(
                    maximum_tail_population_fraction),
            )
            if electron_cache is not None:
                electron_cache[electron_key] = electron_solution
        chemistry = self.collision_chemistry.evaluate(
            electron_solution, electron_condition, densities)
        equivalent_temperature = (
            2.0 / 3.0 * electron_solution.distribution.mean_energy_eV)
        transport_condition = condition.transport_condition(
            equivalent_temperature)
        neutral_transport = neutral_wall_transport_provider.predict(
            transport_condition, densities)
        charged_transport = charged_transport_provider.predict(
            transport_condition, densities)
        wall_energy = wall_energy_provider.predict(
            transport_condition,
            densities,
            charged_transport,
            electron_solution,
        )

        heavy_context = RateContext(
            electron_temperature_eV=equivalent_temperature,
            gas_temperature_K=condition.gas_temperature.value,
        )
        heavy_rates = self.heavy_reaction_network.event_rates_m3_s(
            densities, heavy_context)
        heavy_source_vector = (
            self.heavy_reaction_network.stoichiometric_matrix @ heavy_rates)
        heavy_turnover_vector = (
            np.abs(self.heavy_reaction_network.stoichiometric_matrix)
            @ heavy_rates)
        collision_turnover = {name: 0.0 for name in _SPECIES}
        for mapping in self.collision_chemistry.mappings:
            rate = chemistry.event_rates_m3_s[mapping.reaction_name]
            for name, coefficient in mapping.heavy_reactants.items():
                collision_turnover[name] += coefficient * rate
            for name, coefficient in mapping.heavy_products.items():
                collision_turnover[name] += coefficient * rate
            process = self.electron_solver.collision_deck.processes[
                mapping.process_index]
            if process.electron_number_change:
                collision_turnover["e"] += (
                    abs(process.electron_number_change) * rate)

        active_fraction = condition.active_volume_fraction
        sources = {
            name: active_fraction * (
                chemistry.species_sources_m3_s[name]
                + heavy_source_vector[self._species_index[name]]
            )
            for name in _SPECIES
        }
        turnover = {
            name: active_fraction * (
                collision_turnover[name]
                + heavy_turnover_vector[self._species_index[name]]
            )
            for name in _SPECIES
        }
        feed_density = (
            condition.chlorine_molecule_feed.value
            / condition.neutral_control_volume.value)
        neutral_wall = neutral_transport.evaluate_volume_rates(densities["Cl"])
        atom_wall_loss = (
            active_fraction * neutral_wall.chlorine_atom_loss_m3_s)
        molecule_wall_return = (
            active_fraction * neutral_wall.chlorine_molecule_return_m3_s)
        positive_wall_loss = {
            name: active_fraction * densities[name]
            * charged_transport.wall_loss_frequency_s_inv(name)
            for name in _POSITIVE_IONS
        }
        external = {
            "Cl2": (
                feed_density,
                -exhaust_frequency * densities["Cl2"],
                molecule_wall_return,
                positive_wall_loss["Cl2+"],
            ),
            "Cl": (
                -exhaust_frequency * densities["Cl"],
                -atom_wall_loss,
                positive_wall_loss["Cl+"],
            ),
            "Cl2+": (-positive_wall_loss["Cl2+"],),
            "Cl+": (-positive_wall_loss["Cl+"],),
            "Cl-": (),
        }
        balances = {
            name: (sources[name] + sum(terms)) / max(
                turnover[name] + sum(abs(value) for value in terms), 1.0)
            for name, terms in external.items()
        }
        balances["quasineutrality"] = (
            densities["Cl2+"] + densities["Cl+"]
            - densities["Cl-"] - densities["e"]
        ) / max(
            densities["Cl2+"] + densities["Cl+"],
            densities["Cl-"] + densities["e"],
            1.0,
        )
        balances["neutral_pressure"] = (
            densities["Cl2"] + densities["Cl"]
            - condition.target_neutral_density_m3
        ) / condition.target_neutral_density_m3

        collisional_power = (
            active_fraction * chemistry.collisional_field_power_gain_W_m3)
        electron_wall_energy = (
            electron_solution.transport_moments
            .mean_wall_loss_electron_energy_eV)
        charged_wall_power = E_CHARGE_C * sum(
            positive_wall_loss[name]
            * (electron_wall_energy + wall_energy.energy_eV_per_lost_ion[name])
            for name in _POSITIVE_IONS
        )
        modeled_power = collisional_power + charged_wall_power
        balances["electron_power"] = (
            condition.absorbed_power_density_W_m3 - modeled_power
        ) / max(
            condition.absorbed_power_density_W_m3, modeled_power, 1.0)
        residual_order = (
            "Cl2", "Cl", "Cl2+", "Cl+", "Cl-",
            "quasineutrality", "neutral_pressure", "electron_power",
        )
        return {
            "residual": np.asarray([balances[name] for name in residual_order]),
            "balances": balances,
            "densities": densities,
            "exhaust_frequency": exhaust_frequency,
            "reduced_field_Td": reduced_field_Td,
            "electron_solution": electron_solution,
            "electron_condition": electron_condition,
            "chemistry": chemistry,
            "charged_transport": charged_transport,
            "neutral_transport": neutral_transport,
            "wall_energy": wall_energy,
            "positive_wall_loss": positive_wall_loss,
            "collisional_power": collisional_power,
            "charged_wall_power": charged_wall_power,
            "modeled_power": modeled_power,
        }

    def solve(
        self,
        condition: EEDFChlorineCondition,
        *,
        charged_transport_provider: ChlorineChargedTransportProvider,
        neutral_wall_transport_provider: ChlorineNeutralWallTransportProvider,
        wall_energy_provider: PositiveIonWallEnergyProvider,
        initial_densities_m3: Mapping[str, float] | None = None,
        initial_exhaust_loss_frequency_s_inv: float | None = None,
        initial_reduced_electric_field_Td: float | None = None,
        residual_tolerance: float = 1.0e-7,
        maximum_evaluations: int = 1200,
        maximum_tail_population_fraction: float = 1.0e-6,
    ) -> EEDFChlorineSolution:
        if not isinstance(condition, EEDFChlorineCondition):
            raise TypeError("an EEPF chlorine condition is required")
        target_density = condition.target_neutral_density_m3
        if initial_densities_m3 is None:
            seed = min(1.0e16, 1.0e-3 * target_density)
            initial = {
                "Cl2": 0.7 * target_density,
                "Cl": 0.3 * target_density,
                "Cl2+": seed,
                "Cl+": seed,
                "Cl-": seed,
                "e": seed,
            }
        else:
            initial = dict(_mapping(
                initial_densities_m3, _SPECIES, positive=True))
        initial_exhaust = (
            condition.chlorine_molecule_feed.value
            / (target_density * condition.neutral_control_volume.value)
            if initial_exhaust_loss_frequency_s_inv is None
            else float(initial_exhaust_loss_frequency_s_inv)
        )
        field_bounds = condition.reduced_field_bounds_Td
        initial_field = (
            math.sqrt(field_bounds[0] * field_bounds[1])
            if initial_reduced_electric_field_Td is None
            else float(initial_reduced_electric_field_Td)
        )
        if (
            not math.isfinite(initial_exhaust)
            or initial_exhaust <= 0.0
            or not field_bounds[0] <= initial_field <= field_bounds[1]
        ):
            raise ValueError("invalid EEPF chlorine initial state")
        initial_vector = np.asarray(
            [initial[name] for name in _STATE_ORDER]
            + [initial_exhaust, initial_field])
        density_floor = max(1.0, 1.0e-18 * target_density)
        lower = np.log(
            [density_floor] * 6 + [1.0e-9, field_bounds[0]])
        upper = np.log(
            [10.0 * target_density] * 6 + [1.0e6, field_bounds[1]])
        # This cache is solve-local: independent conditions share no mutable
        # state, while finite-difference probes that leave E/N and composition
        # unchanged reuse the exact same immutable EEPF solution.
        electron_cache: dict[
            tuple[float, ...], TwoTermBoltzmannSolution
        ] = {}

        def residual(log_state):
            return self._ledger(
                log_state,
                condition=condition,
                charged_transport_provider=charged_transport_provider,
                neutral_wall_transport_provider=neutral_wall_transport_provider,
                wall_energy_provider=wall_energy_provider,
                maximum_tail_population_fraction=(
                    maximum_tail_population_fraction),
                electron_cache=electron_cache,
            )["residual"]

        result = least_squares(
            residual,
            x0=np.clip(np.log(initial_vector), lower, upper),
            bounds=(lower, upper),
            xtol=1.0e-11,
            ftol=1.0e-11,
            gtol=1.0e-11,
            max_nfev=int(maximum_evaluations),
        )
        ledger = self._ledger(
            result.x,
            condition=condition,
            charged_transport_provider=charged_transport_provider,
            neutral_wall_transport_provider=neutral_wall_transport_provider,
            wall_energy_provider=wall_energy_provider,
            maximum_tail_population_fraction=maximum_tail_population_fraction,
            electron_cache=electron_cache,
        )
        max_residual = float(np.max(np.abs(ledger["residual"])))
        if (
            not result.success
            or not math.isfinite(max_residual)
            or max_residual > float(residual_tolerance)
        ):
            dominant_balance = max(
                ledger["balances"],
                key=lambda name: abs(ledger["balances"][name]),
            )
            raise RuntimeError(
                "EEPF chlorine power solve failed conservation gate: "
                f"success={result.success}, residual={max_residual}, "
                f"dominant={dominant_balance}="
                f"{ledger['balances'][dominant_balance]}, "
                f"E/N={ledger['reduced_field_Td']} Td, "
                f"message={result.message}"
            )
        densities = ledger["densities"]
        charged_transport = ledger["charged_transport"]
        axial_flux = {
            name: densities[name]
            * charged_transport.positive_ion_transport[
                name].axial_flux_velocity_m_s
            for name in _POSITIVE_IONS
        }
        return EEDFChlorineSolution(
            condition_id=condition.condition_id,
            densities_m3=densities,
            reduced_electric_field_Td=ledger["reduced_field_Td"],
            mean_electron_energy_eV=(
                ledger["electron_solution"].distribution.mean_energy_eV),
            axial_positive_ion_flux_m2_s=axial_flux,
            exhaust_loss_frequency_s_inv=ledger["exhaust_frequency"],
            absorbed_power_density_W_m3=condition.absorbed_power_density_W_m3,
            collisional_power_density_W_m3=ledger["collisional_power"],
            charged_wall_power_density_W_m3=ledger["charged_wall_power"],
            modeled_power_density_W_m3=ledger["modeled_power"],
            normalized_residuals=ledger["balances"],
            electron_solution=ledger["electron_solution"],
            collision_chemistry_state=ledger["chemistry"],
            solver_evaluations=result.nfev,
            absorbed_power_supports_prediction=(
                condition.absorbed_power.supports_prediction),
            wall_energy_supports_prediction=(
                ledger["wall_energy"].supports_prediction),
        )
