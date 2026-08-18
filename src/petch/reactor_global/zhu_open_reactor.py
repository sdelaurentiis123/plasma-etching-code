"""Open fixed-pressure CCP global model for the Zhu NPG80 chemistry.

This module joins the measured parent-feed EEPF operator to the conserved
daughter network.  It solves all heavy-particle balances, quasineutrality,
neutral pressure, throttle/exhaust frequency, reduced field, and absorbed
power in one deterministic log-density system.

Machine transfer inputs remain explicit.  In particular, generator forward
power is not silently treated as absorbed power, and the published NPG80
recipe does not determine reactor height, ion momentum mean free path, wall
state, or mean all-wall sheath energy.  Those quantities define an ensemble
outside this solver; no etched-profile target enters the equations.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from petch.sheath import bohm_speed

from .chlorine_particle_model import BOLTZMANN_J_K
from .electron_collision_chemistry import ElectronCollisionChemistryState
from .electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    TwoTermBoltzmannCondition,
    TwoTermBoltzmannSolution,
)
from .geometry import CylindricalReactor
from .network import E_CHARGE_C, RateContext
from .neutral_transport import solve_cylindrical_neutral_wall_loss
from .zhu_parent_collision_chemistry import ZhuParentCollisionChemistry
from .zhu_supplemental_chemistry import ZhuSupplementalChemistry


_FEED_SPECIES = frozenset({"CHF3", "SF6", "O2"})


def _finite_mapping(
    values: Mapping[str, float], *, nonnegative: bool,
) -> MappingProxyType:
    converted = {str(name): float(value) for name, value in values.items()}
    if (
        not converted
        or any(not name.strip() for name in converted)
        or any(
            not math.isfinite(value)
            or (value < 0.0 if nonnegative else value <= 0.0)
            for value in converted.values()
        )
    ):
        raise ValueError("invalid reactor mapping")
    return MappingProxyType(converted)


@dataclass(frozen=True)
class ZhuOpenReactorCondition:
    """One pressure-controlled condition with explicit machine closures."""

    condition_id: str
    geometry: CylindricalReactor
    neutral_control_volume_m3: float
    pressure_Pa: float
    gas_temperature_K: float
    feed_molecules_s: Mapping[str, float]
    absorbed_power_W: float
    source_frequency_hz: float
    reduced_field_bounds_Td: tuple[float, float]
    ion_temperature_eV: float
    ion_momentum_mean_free_path_m: float
    mean_positive_ion_wall_energy_eV: float
    neutral_reduced_diffusivity_m_inv_s: float
    neutral_wall_probabilities: Mapping[str, float]
    source: str
    absorbed_power_source: str
    machine_closure_source: str

    def __post_init__(self):
        scalars = np.asarray((
            self.neutral_control_volume_m3,
            self.pressure_Pa,
            self.gas_temperature_K,
            self.absorbed_power_W,
            self.source_frequency_hz,
            self.ion_temperature_eV,
            self.ion_momentum_mean_free_path_m,
            self.mean_positive_ion_wall_energy_eV,
            self.neutral_reduced_diffusivity_m_inv_s,
        ), dtype=float)
        bounds = tuple(float(value) for value in self.reduced_field_bounds_Td)
        feed = _finite_mapping(self.feed_molecules_s, nonnegative=False)
        probabilities = _finite_mapping(
            self.neutral_wall_probabilities, nonnegative=True)
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars <= 0.0)
            or self.neutral_control_volume_m3 < self.geometry.volume_m3
            or set(feed) != _FEED_SPECIES
            or len(bounds) != 2
            or not 0.0 < bounds[0] < bounds[1]
            or any(value > 1.0 for value in probabilities.values())
            or not str(self.source).strip()
            or not str(self.absorbed_power_source).strip()
            or not str(self.machine_closure_source).strip()
        ):
            raise ValueError("invalid Zhu open-reactor condition")
        object.__setattr__(self, "feed_molecules_s", feed)
        object.__setattr__(self, "neutral_wall_probabilities", probabilities)
        object.__setattr__(self, "reduced_field_bounds_Td", bounds)
        for name, value in zip((
            "neutral_control_volume_m3",
            "pressure_Pa",
            "gas_temperature_K",
            "absorbed_power_W",
            "source_frequency_hz",
            "ion_temperature_eV",
            "ion_momentum_mean_free_path_m",
            "mean_positive_ion_wall_energy_eV",
            "neutral_reduced_diffusivity_m_inv_s",
        ), scalars):
            object.__setattr__(self, name, float(value))

    @property
    def target_neutral_density_m3(self) -> float:
        return float(
            self.pressure_Pa / (BOLTZMANN_J_K * self.gas_temperature_K))

    @property
    def active_volume_fraction(self) -> float:
        return float(self.geometry.volume_m3 / self.neutral_control_volume_m3)

    @property
    def absorbed_power_density_W_m3(self) -> float:
        return float(self.absorbed_power_W / self.neutral_control_volume_m3)

    @property
    def supports_unique_machine_state(self) -> bool:
        return False


@dataclass(frozen=True)
class ZhuOpenReactorSolution:
    condition_id: str
    densities_m3: Mapping[str, float]
    reduced_electric_field_Td: float
    mean_electron_energy_eV: float
    exhaust_loss_frequency_s_inv: float
    axial_positive_ion_flux_m2_s: Mapping[str, float]
    neutral_thermal_flux_m2_s: Mapping[str, float]
    positive_ion_wall_loss_m3_s: Mapping[str, float]
    neutral_wall_loss_m3_s: Mapping[str, float]
    absorbed_power_density_W_m3: float
    parent_collision_power_density_W_m3: float
    supplemental_collision_power_density_W_m3: float
    charged_wall_power_density_W_m3: float
    normalized_residuals: Mapping[str, float]
    electron_solution: TwoTermBoltzmannSolution
    parent_collision_state: ElectronCollisionChemistryState
    solver_evaluations: int
    electron_collision_basis_neutral_fraction: float
    machine_closure_supports_prediction: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux_prediction: bool = False
    supports_feature_depth_prediction: bool = False

    def __post_init__(self):
        densities = _finite_mapping(self.densities_m3, nonnegative=False)
        axial = _finite_mapping(
            self.axial_positive_ion_flux_m2_s, nonnegative=True)
        neutral_flux = _finite_mapping(
            self.neutral_thermal_flux_m2_s, nonnegative=True)
        positive_loss = _finite_mapping(
            self.positive_ion_wall_loss_m3_s, nonnegative=True)
        neutral_loss = _finite_mapping(
            self.neutral_wall_loss_m3_s, nonnegative=True)
        residuals = {
            str(name): float(value)
            for name, value in self.normalized_residuals.items()
        }
        values = np.asarray((
            self.reduced_electric_field_Td,
            self.mean_electron_energy_eV,
            self.exhaust_loss_frequency_s_inv,
            self.absorbed_power_density_W_m3,
            self.parent_collision_power_density_W_m3,
            self.supplemental_collision_power_density_W_m3,
            self.charged_wall_power_density_W_m3,
            self.electron_collision_basis_neutral_fraction,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values[:4] <= 0.0)
            or np.any(values[4:-1] < 0.0)
            or not 0.0 < values[-1] <= 1.0
            or set(axial) != set(positive_loss)
            or set(neutral_flux) != set(neutral_loss)
            or not residuals
            or any(not math.isfinite(value) for value in residuals.values())
            or int(self.solver_evaluations) <= 0
        ):
            raise ValueError("invalid Zhu open-reactor solution")
        object.__setattr__(self, "densities_m3", densities)
        object.__setattr__(self, "axial_positive_ion_flux_m2_s", axial)
        object.__setattr__(self, "neutral_thermal_flux_m2_s", neutral_flux)
        object.__setattr__(self, "positive_ion_wall_loss_m3_s", positive_loss)
        object.__setattr__(self, "neutral_wall_loss_m3_s", neutral_loss)
        object.__setattr__(self, "normalized_residuals", MappingProxyType(residuals))
        object.__setattr__(self, "solver_evaluations", int(self.solver_evaluations))
        for name, value in zip((
            "reduced_electric_field_Td",
            "mean_electron_energy_eV",
            "exhaust_loss_frequency_s_inv",
            "absorbed_power_density_W_m3",
            "parent_collision_power_density_W_m3",
            "supplemental_collision_power_density_W_m3",
            "charged_wall_power_density_W_m3",
            "electron_collision_basis_neutral_fraction",
        ), values):
            object.__setattr__(self, name, float(value))

    @property
    def maximum_normalized_residual(self) -> float:
        return float(max(abs(value) for value in self.normalized_residuals.values()))

    @property
    def total_axial_positive_ion_flux_m2_s(self) -> float:
        return float(sum(self.axial_positive_ion_flux_m2_s.values()))


def positive_ion_wall_return(species_name: str) -> Mapping[str, float]:
    """Return the atom-conserving neutralization product at a material wall."""
    direct = {
        "CF3+": "CF3", "CHF2+": "CHF2", "CF2+": "CF2",
        "CHF+": "CHF", "CF+": "CF", "CH+": "CH",
        "F+": "F", "F2+": "F2",
        "SF5+": "SF5", "SF4+": "SF4", "SF3+": "SF3",
        "SF2+": "SF2", "SF+": "SF", "S+": "S",
        "SF4++": "SF4", "SF2++": "SF2",
        "O2+": "O2", "O+": "O",
        "H2+": "H2", "H+": "H",
    }
    try:
        return MappingProxyType({direct[species_name]: 1.0})
    except KeyError as exc:
        raise KeyError(f"no wall-return closure for {species_name}") from exc


class ZhuOpenReactorModel:
    """Deterministic coupled particle/EEPF/power solver."""

    def __init__(
        self,
        electron_solver: DeterministicTwoTermBoltzmannSolver,
        parent_chemistry: ZhuParentCollisionChemistry,
        supplemental_chemistry: ZhuSupplementalChemistry,
    ):
        if not isinstance(electron_solver, DeterministicTwoTermBoltzmannSolver):
            raise TypeError("a deterministic electron solver is required")
        if not isinstance(parent_chemistry, ZhuParentCollisionChemistry):
            raise TypeError("Zhu parent chemistry is required")
        if not isinstance(supplemental_chemistry, ZhuSupplementalChemistry):
            raise TypeError("Zhu supplemental chemistry is required")
        if electron_solver.collision_deck is not parent_chemistry.mixed_deck:
            raise ValueError("electron solver and parent chemistry must share a deck")
        reactor_species = supplemental_chemistry.network.species
        reactor_names = {species.name for species in reactor_species}
        parent_names = {species.name for species in parent_chemistry.species}
        if not parent_names <= reactor_names:
            raise ValueError("reactor species do not contain every parent product")
        self.electron_solver = electron_solver
        self.parent_chemistry = parent_chemistry
        self.supplemental_chemistry = supplemental_chemistry
        self.species = reactor_species
        self.species_by_name = {species.name: species for species in self.species}
        self.species_order = tuple(species.name for species in self.species)
        self.heavy_order = tuple(name for name in self.species_order if name != "e")
        self.neutral_names = tuple(
            species.name for species in self.species
            if species.role in {"neutral", "excited_neutral"})
        self.positive_names = tuple(
            species.name for species in self.species
            if species.role == "positive_ion")
        self.negative_names = tuple(
            species.name for species in self.species
            if species.role == "negative_ion")
        self.parent_species_names = tuple(
            species.name for species in self.parent_chemistry.species)
        self._supplemental_index = {
            name: index for index, name in enumerate(
                self.supplemental_chemistry.network.species_names)
        }
        for name in self.positive_names:
            products = positive_ion_wall_return(name)
            self._assert_wall_return_conserves(name, products)

    def _assert_wall_return_conserves(
        self, ion_name: str, products: Mapping[str, float],
    ) -> None:
        ion = self.species_by_name[ion_name]
        elements = set(ion.composition)
        for element in elements:
            product_count = sum(
                amount * self.species_by_name[name].composition.get(element, 0)
                for name, amount in products.items()
            )
            if product_count != ion.composition.get(element, 0):
                raise ValueError(f"wall return for {ion_name} loses {element}")

    def _neutral_wall_frequencies(
        self, condition: ZhuOpenReactorCondition,
    ) -> dict[str, float]:
        frequencies = {name: 0.0 for name in self.neutral_names}
        reference_mass = self.species_by_name["CF3"].mass_amu
        for name, probability in condition.neutral_wall_probabilities.items():
            if name not in frequencies:
                raise ValueError(f"wall probability names non-neutral {name}")
            species = self.species_by_name[name]
            diffusivity = (
                condition.neutral_reduced_diffusivity_m_inv_s
                / condition.target_neutral_density_m3
                * math.sqrt(reference_mass / species.mass_amu)
            )
            mean_speed = math.sqrt(
                8.0 * BOLTZMANN_J_K * condition.gas_temperature_K
                / (math.pi * species.mass_amu * 1.66053906660e-27)
            )
            frequencies[name] = solve_cylindrical_neutral_wall_loss(
                geometry=condition.geometry,
                diffusivity_m2_s=diffusivity,
                mean_thermal_speed_m_s=mean_speed,
                wall_reaction_probability=probability,
            ).exact_loss_frequency_s_inv
        return frequencies

    def _charged_wall_transport(
        self,
        densities: Mapping[str, float],
        equivalent_temperature_eV: float,
        condition: ZhuOpenReactorCondition,
    ) -> tuple[dict[str, float], dict[str, float]]:
        negative_charge = sum(
            -self.species_by_name[name].charge_number * densities[name]
            for name in self.negative_names)
        electronegativity = negative_charge / densities["e"]
        wall_frequencies = {}
        axial_velocities = {}
        for name in self.positive_names:
            species = self.species_by_name[name]
            charge = species.charge_number
            speed = (
                math.sqrt(charge)
                * bohm_speed(equivalent_temperature_eV, species.mass_amu)
            )
            edge = condition.geometry.electronegative_edge_factors(
                electronegativity=electronegativity,
                electron_to_ion_temperature_ratio=(
                    equivalent_temperature_eV / condition.ion_temperature_eV),
                ion_mean_free_path_m=(
                    condition.ion_momentum_mean_free_path_m),
                include_high_pressure_diffusion=False,
            )
            wall_frequencies[name] = (
                condition.geometry.effective_loss_area_m2(edge)
                / condition.geometry.volume_m3
                * speed
            )
            axial_velocities[name] = edge.axial * speed
        return wall_frequencies, axial_velocities

    def jacobian_sparsity(self):
        """Return the exact structural dependency graph for the log solve.

        This graph changes no physics.  It lets SciPy color independent
        density columns together, avoiding redundant Boltzmann solves during
        finite-difference Jacobian construction.
        """
        residual_order = (
            *self.heavy_order,
            "quasineutrality", "neutral_pressure", "electron_power",
        )
        column_order = (*self.species_order, "exhaust", "reduced_field")
        row_index = {name: index for index, name in enumerate(residual_order)}
        column_index = {name: index for index, name in enumerate(column_order)}
        matrix = lil_matrix(
            (len(residual_order), len(column_order)), dtype=bool)
        eedf_columns = {"e", "CHF3", "SF6", "O2", "reduced_field"}

        network = self.supplemental_chemistry.network
        for reaction in network.reactions:
            affected = (
                set(reaction.reactants) | set(reaction.products)
            ) - {"e"}
            dependencies = set(reaction.kinetic_orders)
            if "e" in reaction.kinetic_orders:
                dependencies |= eedf_columns
            for affected_name in affected:
                for dependency in dependencies:
                    matrix[
                        row_index[affected_name], column_index[dependency]
                    ] = True
        for mapping in self.parent_chemistry.collision_chemistry.mappings:
            affected = (
                set(mapping.heavy_reactants)
                | set(mapping.heavy_products)
            )
            for affected_name in affected:
                for dependency in eedf_columns:
                    matrix[
                        row_index[affected_name], column_index[dependency]
                    ] = True
        for name in self.neutral_names:
            matrix[row_index[name], column_index[name]] = True
            matrix[row_index[name], column_index["exhaust"]] = True
        charged_transport_dependencies = {
            "e", *self.negative_names, "CHF3", "SF6", "O2", "reduced_field",
        }
        for ion_name in self.positive_names:
            dependencies = charged_transport_dependencies | {ion_name}
            for dependency in dependencies:
                matrix[row_index[ion_name], column_index[dependency]] = True
            for product in positive_ion_wall_return(ion_name):
                for dependency in dependencies:
                    matrix[row_index[product], column_index[dependency]] = True
        for name in ("F", "H", "O", "O(1d)"):
            if name in self.neutral_names:
                matrix[row_index[name], column_index[name]] = True
        matrix[row_index["F2"], column_index["F"]] = True
        matrix[row_index["H2"], column_index["H"]] = True
        matrix[row_index["O2"], column_index["O"]] = True
        matrix[row_index["O"], column_index["O(1d)"]] = True
        for name in ("e", *self.positive_names, *self.negative_names):
            matrix[
                row_index["quasineutrality"], column_index[name]
            ] = True
        for name in self.neutral_names:
            matrix[row_index["neutral_pressure"], column_index[name]] = True
        # Power sees every electron-reaction reactant, all wall-current
        # carriers, and the EEPF state.
        power_dependencies = set(eedf_columns)
        for reaction in network.reactions:
            if "e" in reaction.reactants or "e" in reaction.products:
                power_dependencies |= set(reaction.kinetic_orders)
        power_dependencies |= set(self.positive_names)
        power_dependencies |= set(self.negative_names)
        for dependency in power_dependencies:
            matrix[
                row_index["electron_power"], column_index[dependency]
            ] = True
        return matrix.tocsr()

    def _ledger(
        self,
        log_state: np.ndarray,
        *,
        condition: ZhuOpenReactorCondition,
        maximum_tail_population_fraction: float,
        electron_cache: dict[tuple[float, ...], TwoTermBoltzmannSolution],
        electron_continuation: dict[str, TwoTermBoltzmannSolution],
    ) -> dict[str, object]:
        values = np.exp(np.asarray(log_state, dtype=float))
        densities = dict(zip(self.species_order, values[:len(self.species_order)]))
        exhaust_frequency = float(values[-2])
        reduced_field = float(values[-1])
        target_density = sum(densities[name] for name in _FEED_SPECIES)
        target_fractions = {
            name: densities[name] / target_density for name in _FEED_SPECIES
        }
        electron_condition = TwoTermBoltzmannCondition(
            reduced_electric_field_Td=reduced_field,
            gas_temperature_K=condition.gas_temperature_K,
            target_mole_fractions=target_fractions,
            growth_model="temporal_growth",
            angular_field_frequency_over_density_m3_s=(
                2.0 * math.pi * condition.source_frequency_hz / target_density),
        )
        electron_key = (
            reduced_field,
            *(target_fractions[name] for name in sorted(_FEED_SPECIES)),
            electron_condition.angular_field_frequency_over_density_m3_s,
            maximum_tail_population_fraction,
        )
        electron_solution = electron_cache.get(electron_key)
        if electron_solution is None:
            electron_solution = self.electron_solver.solve(
                electron_condition,
                initial_solution=electron_continuation.get("latest"),
                relative_tolerance=1.0e-8,
                maximum_iterations=220,
                maximum_tail_population_fraction=maximum_tail_population_fraction,
            )
            electron_cache[electron_key] = electron_solution
        electron_continuation["latest"] = electron_solution
        parent_densities = {
            name: densities[name] for name in self.parent_species_names
        }
        parent_state = self.parent_chemistry.collision_chemistry.evaluate(
            electron_solution, electron_condition, parent_densities,
            closure_relative_tolerance=2.0e-7,
        )
        equivalent_temperature = (
            2.0 / 3.0 * electron_solution.distribution.mean_energy_eV)
        context = RateContext(
            electron_temperature_eV=equivalent_temperature,
            gas_temperature_K=condition.gas_temperature_K,
        )
        network = self.supplemental_chemistry.network
        supplemental_rates = network.event_rates_m3_s(densities, context)
        supplemental_source = network.stoichiometric_matrix @ supplemental_rates
        supplemental_turnover = (
            np.abs(network.stoichiometric_matrix) @ supplemental_rates)
        parent_turnover = {name: 0.0 for name in self.species_order}
        for mapping in self.parent_chemistry.collision_chemistry.mappings:
            rate = parent_state.event_rates_m3_s[mapping.reaction_name]
            for name, coefficient in mapping.heavy_reactants.items():
                parent_turnover[name] += coefficient * rate
            for name, coefficient in mapping.heavy_products.items():
                parent_turnover[name] += coefficient * rate
            process = self.parent_chemistry.mixed_deck.processes[
                mapping.process_index]
            if process.electron_number_change:
                parent_turnover["e"] += abs(process.electron_number_change) * rate
        active_fraction = condition.active_volume_fraction
        internal_sources = {
            name: active_fraction * (
                parent_state.species_sources_m3_s.get(name, 0.0)
                + supplemental_source[self._supplemental_index[name]]
            )
            for name in self.species_order
        }
        turnover = {
            name: active_fraction * (
                parent_turnover[name]
                + supplemental_turnover[self._supplemental_index[name]]
            )
            for name in self.species_order
        }
        neutral_wall_frequency = self._neutral_wall_frequencies(condition)
        positive_wall_frequency, axial_velocity = self._charged_wall_transport(
            densities, equivalent_temperature, condition)
        external_terms = {name: [] for name in self.heavy_order}
        for name in self.neutral_names:
            external_terms[name].append(-exhaust_frequency * densities[name])
        for name, feed in condition.feed_molecules_s.items():
            external_terms[name].append(
                feed / condition.neutral_control_volume_m3)
        neutral_wall_loss = {}
        neutral_thermal_flux = {}
        for name, frequency in neutral_wall_frequency.items():
            species = self.species_by_name[name]
            mean_speed = math.sqrt(
                8.0 * BOLTZMANN_J_K * condition.gas_temperature_K
                / (math.pi * species.mass_amu * 1.66053906660e-27)
            )
            neutral_thermal_flux[name] = 0.25 * densities[name] * mean_speed
            loss = active_fraction * frequency * densities[name]
            neutral_wall_loss[name] = loss
            external_terms[name].append(-loss)
            if name == "F":
                external_terms["F2"].append(0.5 * loss)
            elif name == "H":
                external_terms["H2"].append(0.5 * loss)
            elif name == "O":
                external_terms["O2"].append(0.5 * loss)
            elif name == "O(1d)":
                external_terms["O"].append(loss)
        positive_wall_loss = {}
        axial_flux = {}
        for name, frequency in positive_wall_frequency.items():
            loss = active_fraction * frequency * densities[name]
            positive_wall_loss[name] = loss
            axial_flux[name] = densities[name] * axial_velocity[name]
            external_terms[name].append(-loss)
            for product, amount in positive_ion_wall_return(name).items():
                external_terms[product].append(amount * loss)
        balances = {}
        for name in self.heavy_order:
            terms = external_terms[name]
            numerator = internal_sources[name] + sum(terms)
            denominator = max(
                turnover[name] + sum(abs(value) for value in terms), 1.0)
            balances[name] = numerator / denominator
        positive_charge = sum(
            self.species_by_name[name].charge_number * densities[name]
            for name in self.positive_names)
        negative_charge = densities["e"] + sum(
            -self.species_by_name[name].charge_number * densities[name]
            for name in self.negative_names)
        balances["quasineutrality"] = (
            positive_charge - negative_charge
        ) / max(positive_charge, negative_charge, 1.0)
        total_neutral = sum(densities[name] for name in self.neutral_names)
        balances["neutral_pressure"] = (
            total_neutral - condition.target_neutral_density_m3
        ) / condition.target_neutral_density_m3
        parent_power = (
            active_fraction * parent_state.collisional_field_power_gain_W_m3)
        supplemental_power = active_fraction * E_CHARGE_C * sum(
            reaction.electron_energy_loss_rate_eV_m3_s(densities, context)
            for reaction in network.reactions
        )
        electron_wall_energy = (
            electron_solution.transport_moments.mean_wall_loss_electron_energy_eV)
        charged_wall_power = E_CHARGE_C * sum(
            positive_wall_loss[name]
            * self.species_by_name[name].charge_number
            * (electron_wall_energy + condition.mean_positive_ion_wall_energy_eV)
            for name in self.positive_names
        )
        modeled_power = parent_power + supplemental_power + charged_wall_power
        balances["electron_power"] = (
            condition.absorbed_power_density_W_m3 - modeled_power
        ) / max(condition.absorbed_power_density_W_m3, modeled_power, 1.0)
        residual_order = (
            *self.heavy_order,
            "quasineutrality", "neutral_pressure", "electron_power",
        )
        return {
            "residual": np.asarray([balances[name] for name in residual_order]),
            "balances": balances,
            "densities": densities,
            "exhaust_frequency": exhaust_frequency,
            "reduced_field": reduced_field,
            "electron_solution": electron_solution,
            "parent_state": parent_state,
            "equivalent_temperature": equivalent_temperature,
            "positive_wall_loss": positive_wall_loss,
            "neutral_wall_loss": neutral_wall_loss,
            "neutral_thermal_flux": neutral_thermal_flux,
            "axial_flux": axial_flux,
            "parent_power": parent_power,
            "supplemental_power": supplemental_power,
            "charged_wall_power": charged_wall_power,
        }

    def solve(
        self,
        condition: ZhuOpenReactorCondition,
        *,
        initial_densities_m3: Mapping[str, float] | None = None,
        initial_exhaust_loss_frequency_s_inv: float | None = None,
        initial_reduced_electric_field_Td: float | None = None,
        residual_tolerance: float = 2.0e-6,
        maximum_evaluations: int = 1000,
        maximum_tail_population_fraction: float = 2.0e-6,
    ) -> ZhuOpenReactorSolution:
        if not isinstance(condition, ZhuOpenReactorCondition):
            raise TypeError("a Zhu open-reactor condition is required")
        target_density = condition.target_neutral_density_m3
        if initial_densities_m3 is None:
            total_feed = sum(condition.feed_molecules_s.values())
            feed_fraction = {
                name: value / total_feed
                for name, value in condition.feed_molecules_s.items()
            }
            density_floor = max(1.0e4, 1.0e-16 * target_density)
            densities = {name: density_floor for name in self.species_order}
            # A log-density Newton step cannot discover a strongly produced
            # daughter that starts tens of decades below its balance scale:
            # d(nu*n)/d(log n)=nu*n then vanishes numerically.  Reserve a
            # composition-neutral fraction for every daughter.  This is only
            # an interior numerical seed; pressure and all species balances
            # still determine the converged inventory.
            neutral_seed_fraction = 0.80
            for name, fraction in feed_fraction.items():
                densities[name] = neutral_seed_fraction * fraction * target_density
            daughter_neutrals = tuple(
                name for name in self.neutral_names if name not in _FEED_SPECIES)
            daughter_seed = (
                (1.0 - neutral_seed_fraction) * target_density
                / len(daughter_neutrals)
            )
            for name in daughter_neutrals:
                densities[name] = daughter_seed
            charged_seed = min(1.0e15, 1.0e-5 * target_density)
            for name in self.positive_names:
                densities[name] = charged_seed / len(self.positive_names)
            for name in self.negative_names:
                densities[name] = 0.5 * charged_seed / len(self.negative_names)
            densities["e"] = 0.5 * charged_seed
        else:
            densities = {
                str(name): float(value)
                for name, value in initial_densities_m3.items()
            }
            if (
                set(densities) != set(self.species_order)
                or any(not math.isfinite(value) or value <= 0.0
                       for value in densities.values())
            ):
                raise ValueError("invalid initial density state")
        initial_exhaust = (
            sum(condition.feed_molecules_s.values())
            / (target_density * condition.neutral_control_volume_m3)
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
            raise ValueError("invalid initial exhaust or reduced field")
        initial = np.asarray([
            *(densities[name] for name in self.species_order),
            initial_exhaust,
            initial_field,
        ])
        density_floor = max(1.0, 1.0e-24 * target_density)
        lower = np.log([
            *([density_floor] * len(self.species_order)),
            1.0e-9,
            field_bounds[0],
        ])
        upper = np.log([
            *([10.0 * target_density] * len(self.species_order)),
            1.0e6,
            field_bounds[1],
        ])
        electron_cache: dict[tuple[float, ...], TwoTermBoltzmannSolution] = {}
        electron_continuation: dict[str, TwoTermBoltzmannSolution] = {}

        def residual(log_state):
            return self._ledger(
                log_state,
                condition=condition,
                maximum_tail_population_fraction=(
                    maximum_tail_population_fraction),
                electron_cache=electron_cache,
                electron_continuation=electron_continuation,
            )["residual"]

        result = least_squares(
            residual,
            x0=np.clip(np.log(initial), lower, upper),
            bounds=(lower, upper),
            xtol=2.0e-10,
            ftol=2.0e-10,
            gtol=2.0e-10,
            max_nfev=int(maximum_evaluations),
            jac_sparsity=self.jacobian_sparsity(),
        )
        ledger = self._ledger(
            result.x,
            condition=condition,
            maximum_tail_population_fraction=maximum_tail_population_fraction,
            electron_cache=electron_cache,
            electron_continuation=electron_continuation,
        )
        maximum_residual = float(np.max(np.abs(ledger["residual"])))
        if (
            not result.success
            or not math.isfinite(maximum_residual)
            or maximum_residual > residual_tolerance
        ):
            dominant = max(
                ledger["balances"],
                key=lambda name: abs(ledger["balances"][name]),
            )
            raise RuntimeError(
                "Zhu open-reactor solve failed conservation gate: "
                f"success={result.success}, residual={maximum_residual}, "
                f"dominant={dominant}={ledger['balances'][dominant]}, "
                f"E/N={ledger['reduced_field']} Td, message={result.message}"
            )
        return ZhuOpenReactorSolution(
            condition_id=condition.condition_id,
            densities_m3=ledger["densities"],
            reduced_electric_field_Td=ledger["reduced_field"],
            mean_electron_energy_eV=(
                ledger["electron_solution"].distribution.mean_energy_eV),
            exhaust_loss_frequency_s_inv=ledger["exhaust_frequency"],
            axial_positive_ion_flux_m2_s=ledger["axial_flux"],
            neutral_thermal_flux_m2_s=ledger["neutral_thermal_flux"],
            positive_ion_wall_loss_m3_s=ledger["positive_wall_loss"],
            neutral_wall_loss_m3_s=ledger["neutral_wall_loss"],
            absorbed_power_density_W_m3=condition.absorbed_power_density_W_m3,
            parent_collision_power_density_W_m3=ledger["parent_power"],
            supplemental_collision_power_density_W_m3=(
                ledger["supplemental_power"]),
            charged_wall_power_density_W_m3=ledger["charged_wall_power"],
            normalized_residuals=ledger["balances"],
            electron_solution=ledger["electron_solution"],
            parent_collision_state=ledger["parent_state"],
            solver_evaluations=result.nfev,
            electron_collision_basis_neutral_fraction=(
                sum(ledger["densities"][name] for name in _FEED_SPECIES)
                / sum(
                    ledger["densities"][name]
                    for name in self.neutral_names
                )
            ),
        )
