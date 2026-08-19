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
from .network import E_CHARGE_C, RateContext, ReactionNetwork
from .neutral_transport import solve_cylindrical_neutral_wall_loss
from .zhu_daughter_electron_collisions import ZhuAugmentedCollisionChemistry
from .zhu_parent_collision_chemistry import ZhuParentCollisionChemistry
from .zhu_supplemental_chemistry import ZhuSupplementalChemistry


_FEED_SPECIES = ("CHF3", "SF6", "O2")
_FEED_SPECIES_SET = frozenset(_FEED_SPECIES)


def compile_bimolecular_kinetic_pairs(
    network: ReactionNetwork,
) -> np.ndarray:
    """Compile one exact pair of density indices per bimolecular event."""
    index = {
        name: position for position, name in enumerate(network.species_names)
    }
    pairs = []
    for reaction in network.reactions:
        expanded = []
        for name, order in reaction.kinetic_orders.items():
            integer_order = int(order)
            if order != integer_order or integer_order < 0:
                raise ValueError(
                    "compiled mass action requires integer kinetic orders")
            expanded.extend([index[name]] * integer_order)
        if len(expanded) != 2:
            raise ValueError(
                "compiled mass action requires bimolecular reactions")
        pairs.append(expanded)
    result = np.asarray(pairs, dtype=np.intp)
    result.setflags(write=False)
    return result


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


def wall_resolved_charged_power_density_W_m3(
    *,
    powered_positive_ion_loss_m3_s: Mapping[str, float],
    grounded_positive_ion_loss_m3_s: Mapping[str, float],
    ion_charge_numbers: Mapping[str, int],
    electron_wall_energy_eV: float,
    powered_electrode_sheath_drop_V: float,
    grounded_surface_sheath_drop_V: float,
) -> tuple[float, float]:
    """Return powered and grounded charged-particle wall-power densities.

    Sheath inputs are potential drops.  Multiplication by the ion charge
    number therefore gives each ion's directed sheath energy in eV.  The
    electron wall-loss energy is allocated over the same ambipolar loss
    channels, so equal sheath drops recover the legacy all-wall expression.
    """
    powered = _finite_mapping(
        powered_positive_ion_loss_m3_s, nonnegative=True)
    grounded = _finite_mapping(
        grounded_positive_ion_loss_m3_s, nonnegative=True)
    charges = {str(name): int(value) for name, value in ion_charge_numbers.items()}
    electron_energy = float(electron_wall_energy_eV)
    powered_drop = float(powered_electrode_sheath_drop_V)
    grounded_drop = float(grounded_surface_sheath_drop_V)
    if (
        set(powered) != set(grounded)
        or set(powered) != set(charges)
        or any(value <= 0 for value in charges.values())
        or any(
            not math.isfinite(value) or value <= 0.0
            for value in (electron_energy, powered_drop, grounded_drop)
        )
    ):
        raise ValueError("invalid wall-resolved charged-power inputs")

    def channel_power(losses, sheath_drop):
        return E_CHARGE_C * sum(
            losses[name]
            * charges[name]
            * (electron_energy + sheath_drop)
            for name in charges
        )

    return (
        float(channel_power(powered, powered_drop)),
        float(channel_power(grounded, grounded_drop)),
    )


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
    powered_electrode_sheath_drop_V: float | None = None
    grounded_surface_sheath_drop_V: float | None = None

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
        resolved_drops = (
            self.powered_electrode_sheath_drop_V,
            self.grounded_surface_sheath_drop_V,
        )
        if any(value is None for value in resolved_drops):
            if not all(value is None for value in resolved_drops):
                raise ValueError(
                    "powered and grounded sheath drops must be supplied together")
        else:
            resolved_drops = tuple(float(value) for value in resolved_drops)
            if any(
                not math.isfinite(value) or value <= 0.0
                for value in resolved_drops
            ):
                raise ValueError("sheath drops must be positive and finite")
        if (
            not str(self.condition_id).strip()
            or not isinstance(self.geometry, CylindricalReactor)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars <= 0.0)
            or self.neutral_control_volume_m3 < self.geometry.volume_m3
            or set(feed) != _FEED_SPECIES_SET
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
        if resolved_drops[0] is not None:
            object.__setattr__(
                self, "powered_electrode_sheath_drop_V", resolved_drops[0])
            object.__setattr__(
                self, "grounded_surface_sheath_drop_V", resolved_drops[1])
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

    @property
    def uses_wall_resolved_sheath_power(self) -> bool:
        return self.powered_electrode_sheath_drop_V is not None


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
    powered_electrode_positive_ion_wall_loss_m3_s: Mapping[str, float]
    grounded_positive_ion_wall_loss_m3_s: Mapping[str, float]
    neutral_wall_loss_m3_s: Mapping[str, float]
    absorbed_power_density_W_m3: float
    parent_collision_power_density_W_m3: float
    supplemental_collision_power_density_W_m3: float
    charged_wall_power_density_W_m3: float
    powered_electrode_charged_wall_power_density_W_m3: float
    grounded_charged_wall_power_density_W_m3: float
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
        powered_positive_loss = _finite_mapping(
            self.powered_electrode_positive_ion_wall_loss_m3_s,
            nonnegative=True,
        )
        grounded_positive_loss = _finite_mapping(
            self.grounded_positive_ion_wall_loss_m3_s, nonnegative=True)
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
            self.powered_electrode_charged_wall_power_density_W_m3,
            self.grounded_charged_wall_power_density_W_m3,
            self.electron_collision_basis_neutral_fraction,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values[:4] <= 0.0)
            or np.any(values[4:-1] < 0.0)
            or not 0.0 < values[-1] <= 1.0
            or set(axial) != set(positive_loss)
            or set(powered_positive_loss) != set(positive_loss)
            or set(grounded_positive_loss) != set(positive_loss)
            or any(
                not math.isclose(
                    positive_loss[name],
                    powered_positive_loss[name] + grounded_positive_loss[name],
                    rel_tol=2.0e-14,
                    abs_tol=0.0,
                )
                for name in positive_loss
            )
            or not math.isclose(
                self.charged_wall_power_density_W_m3,
                self.powered_electrode_charged_wall_power_density_W_m3
                + self.grounded_charged_wall_power_density_W_m3,
                rel_tol=2.0e-14,
                abs_tol=0.0,
            )
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
        object.__setattr__(
            self,
            "powered_electrode_positive_ion_wall_loss_m3_s",
            powered_positive_loss,
        )
        object.__setattr__(
            self, "grounded_positive_ion_wall_loss_m3_s", grounded_positive_loss)
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
            "powered_electrode_charged_wall_power_density_W_m3",
            "grounded_charged_wall_power_density_W_m3",
            "electron_collision_basis_neutral_fraction",
        ), values):
            object.__setattr__(self, name, float(value))

    @property
    def maximum_normalized_residual(self) -> float:
        return float(max(abs(value) for value in self.normalized_residuals.values()))

    @property
    def total_axial_positive_ion_flux_m2_s(self) -> float:
        return float(sum(self.axial_positive_ion_flux_m2_s.values()))

    @property
    def implied_total_neutral_reduced_electric_field_Td(self) -> float:
        """Field divided by total neutral density under transparent daughters.

        The solved electron operator currently normalizes collision frequency
        to the represented CHF3/SF6/O2 basis.  Multiplying by that basis'
        neutral fraction preserves the dimensional electric field while
        expressing it against the full pressure density.  This is a
        diagnostic conversion, not a daughter-collision closure.
        """
        return float(
            self.reduced_electric_field_Td
            * self.electron_collision_basis_neutral_fraction
        )


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
        "H2+": "H2", "H+": "H", "HF+": "HF",
    }
    try:
        return MappingProxyType({direct[species_name]: 1.0})
    except KeyError as exc:
        raise KeyError(f"no wall-return closure for {species_name}") from exc


def log_flux_ratio_residual(conservation_residual: np.ndarray) -> np.ndarray:
    """Desaturate bounded production/loss residuals for nonlinear solves.

    A species balance written as ``(P-L)/(P+L)`` is an excellent bounded
    conservation grade but a poor optimization coordinate: its derivative
    vanishes when a continuation state is orders of magnitude away from the
    new steady inventory.  ``2*atanh(r)`` is exactly ``log(P/L)`` for that
    balance, has the same and only root, and retains useful log-density
    derivatives far from closure.  The physical bounded residual remains the
    post-solve acceptance gate.
    """

    values = np.asarray(conservation_residual, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("conservation residual must be one finite vector")
    limit = 1.0 - 8.0 * np.finfo(float).eps
    return 2.0 * np.arctanh(np.clip(values, -limit, limit))


class ZhuOpenReactorModel:
    """Deterministic coupled particle/EEPF/power solver."""

    def __init__(
        self,
        electron_solver: DeterministicTwoTermBoltzmannSolver,
        parent_chemistry: (
            ZhuParentCollisionChemistry | ZhuAugmentedCollisionChemistry
        ),
        supplemental_chemistry: ZhuSupplementalChemistry,
    ):
        if not isinstance(electron_solver, DeterministicTwoTermBoltzmannSolver):
            raise TypeError("a deterministic electron solver is required")
        if not isinstance(parent_chemistry, (
            ZhuParentCollisionChemistry, ZhuAugmentedCollisionChemistry,
        )):
            raise TypeError("Zhu collision chemistry is required")
        if not isinstance(supplemental_chemistry, ZhuSupplementalChemistry):
            raise TypeError("Zhu supplemental chemistry is required")
        if electron_solver.collision_deck is not parent_chemistry.mixed_deck:
            raise ValueError("electron solver and parent chemistry must share a deck")
        if (
            isinstance(parent_chemistry, ZhuAugmentedCollisionChemistry)
            and set(parent_chemistry.supplemental_reactions_replaced)
            != set(supplemental_chemistry.electron_collision_rows_replaced)
        ):
            raise ValueError(
                "augmented collision rows must be removed from supplemental "
                "chemistry"
            )
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
        self.electron_collision_targets = tuple(
            self.electron_solver.collision_deck.targets
        )
        self._supplemental_index = {
            name: index for index, name in enumerate(
                self.supplemental_chemistry.network.species_names)
        }
        self._supplemental_kinetic_pairs = compile_bimolecular_kinetic_pairs(
            self.supplemental_chemistry.network)
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

    def _supplemental_event_rates_m3_s(
        self,
        densities: Mapping[str, float],
        context: RateContext,
        coefficient_cache: dict[tuple[float, float], np.ndarray],
    ) -> np.ndarray:
        """Evaluate this fixed bimolecular mechanism as a dense vector plan."""
        key = (
            context.electron_temperature_eV,
            context.gas_temperature_K,
        )
        coefficients = coefficient_cache.get(key)
        if coefficients is None:
            coefficients = np.asarray([
                reaction.rate_coefficient.coefficient_si(context)
                for reaction in self.supplemental_chemistry.network.reactions
            ])
            coefficients.setflags(write=False)
            coefficient_cache[key] = coefficients
        density = np.asarray([
            densities[name]
            for name in self.supplemental_chemistry.network.species_names
        ])
        pairs = self._supplemental_kinetic_pairs
        rates = (
            coefficients
            * density[pairs[:, 0]]
            * density[pairs[:, 1]]
        )
        rates.setflags(write=False)
        return rates

    def _charged_wall_transport(
        self,
        densities: Mapping[str, float],
        equivalent_temperature_eV: float,
        condition: ZhuOpenReactorCondition,
    ) -> tuple[
        dict[str, float],
        dict[str, float],
        dict[str, tuple[float, float]],
    ]:
        negative_charge = sum(
            -self.species_by_name[name].charge_number * densities[name]
            for name in self.negative_names)
        electronegativity = negative_charge / densities["e"]
        wall_frequencies = {}
        axial_velocities = {}
        edge = condition.geometry.electronegative_edge_factors(
            electronegativity=electronegativity,
            electron_to_ion_temperature_ratio=(
                equivalent_temperature_eV / condition.ion_temperature_eV),
            ion_mean_free_path_m=condition.ion_momentum_mean_free_path_m,
            include_high_pressure_diffusion=False,
        )
        powered_area_over_volume = (
            edge.axial * math.pi * condition.geometry.radius_m ** 2
            / condition.geometry.volume_m3
        )
        grounded_area_over_volume = (
            (
                edge.axial * math.pi * condition.geometry.radius_m ** 2
                + edge.radial
                * 2.0
                * math.pi
                * condition.geometry.radius_m
                * condition.geometry.length_m
            )
            / condition.geometry.volume_m3
        )
        loss_area_over_volume = powered_area_over_volume + grounded_area_over_volume
        resolved_frequencies = {}
        for name in self.positive_names:
            species = self.species_by_name[name]
            charge = species.charge_number
            speed = (
                math.sqrt(charge)
                * bohm_speed(equivalent_temperature_eV, species.mass_amu)
            )
            wall_frequencies[name] = loss_area_over_volume * speed
            resolved_frequencies[name] = (
                powered_area_over_volume * speed,
                grounded_area_over_volume * speed,
            )
            axial_velocities[name] = edge.axial * speed
        return wall_frequencies, axial_velocities, resolved_frequencies

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
        eedf_columns = {
            "e", *self.electron_collision_targets, "reduced_field",
        }

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
            *eedf_columns, *self.negative_names,
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
        neutral_wall_frequency: Mapping[str, float],
        supplemental_coefficient_cache: dict[tuple[float, float], np.ndarray],
        electron_cache: dict[tuple[float, ...], TwoTermBoltzmannSolution],
        electron_continuation: dict[str, TwoTermBoltzmannSolution],
    ) -> dict[str, object]:
        values = np.exp(np.asarray(log_state, dtype=float))
        densities = dict(zip(self.species_order, values[:len(self.species_order)]))
        exhaust_frequency = float(values[-2])
        reduced_field = float(values[-1])
        target_density = sum(
            densities[name] for name in self.electron_collision_targets
        )
        target_fractions = {
            name: densities[name] / target_density
            for name in self.electron_collision_targets
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
            *(target_fractions[name] for name in self.electron_collision_targets),
            electron_condition.angular_field_frequency_over_density_m3_s,
            maximum_tail_population_fraction,
        )
        electron_solution = electron_cache.get(electron_key)
        if electron_solution is None:
            electron_solution = self.electron_solver.solve(
                electron_condition,
                initial_solution=electron_continuation.get("reference"),
                relative_tolerance=1.0e-8,
                maximum_iterations=220,
                maximum_tail_population_fraction=maximum_tail_population_fraction,
            )
            electron_cache[electron_key] = electron_solution
        # A moving "latest" warm start makes the nonlinear residual depend on
        # finite-difference column order.  Freeze the first solved EEPF as the
        # common reference so identical reactor states are bitwise replayable.
        electron_continuation.setdefault("reference", electron_solution)
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
        supplemental_rates = self._supplemental_event_rates_m3_s(
            densities, context, supplemental_coefficient_cache)
        # This network is small enough that dispatching two 66x259 GEMVs to a
        # multithreaded BLAS costs far more than the arithmetic.  The explicit
        # reductions are deterministic, avoid thread oversubscription inside
        # finite-difference Jacobians, and are algebraically identical.
        weighted_events = (
            network.stoichiometric_matrix
            * supplemental_rates[np.newaxis, :]
        )
        supplemental_source = np.sum(weighted_events, axis=1)
        supplemental_turnover = np.sum(np.abs(weighted_events), axis=1)
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
        (
            positive_wall_frequency,
            axial_velocity,
            resolved_wall_frequency,
        ) = self._charged_wall_transport(
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
        powered_positive_wall_loss = {}
        grounded_positive_wall_loss = {}
        axial_flux = {}
        for name, frequency in positive_wall_frequency.items():
            loss = active_fraction * frequency * densities[name]
            powered_frequency, grounded_frequency = resolved_wall_frequency[name]
            powered_loss = active_fraction * powered_frequency * densities[name]
            grounded_loss = active_fraction * grounded_frequency * densities[name]
            positive_wall_loss[name] = loss
            powered_positive_wall_loss[name] = powered_loss
            grounded_positive_wall_loss[name] = grounded_loss
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
            (
                rate * reaction.electron_energy_loss_eV
                if reaction.electron_energy_loss_eV is not None
                else reaction.electron_energy_loss_rate_eV_m3_s(
                    densities, context)
            )
            for rate, reaction in zip(supplemental_rates, network.reactions)
        )
        electron_wall_energy = (
            electron_solution.transport_moments.mean_wall_loss_electron_energy_eV)
        if condition.uses_wall_resolved_sheath_power:
            powered_sheath_drop = condition.powered_electrode_sheath_drop_V
            grounded_sheath_drop = condition.grounded_surface_sheath_drop_V
        else:
            powered_sheath_drop = condition.mean_positive_ion_wall_energy_eV
            grounded_sheath_drop = condition.mean_positive_ion_wall_energy_eV
        (
            powered_charged_wall_power,
            grounded_charged_wall_power,
        ) = wall_resolved_charged_power_density_W_m3(
            powered_positive_ion_loss_m3_s=powered_positive_wall_loss,
            grounded_positive_ion_loss_m3_s=grounded_positive_wall_loss,
            ion_charge_numbers={
                name: self.species_by_name[name].charge_number
                for name in self.positive_names
            },
            electron_wall_energy_eV=electron_wall_energy,
            powered_electrode_sheath_drop_V=powered_sheath_drop,
            grounded_surface_sheath_drop_V=grounded_sheath_drop,
        )
        charged_wall_power = (
            powered_charged_wall_power + grounded_charged_wall_power)
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
            "internal_sources": internal_sources,
            "turnover": turnover,
            "external_terms": external_terms,
            "positive_wall_loss": positive_wall_loss,
            "powered_positive_wall_loss": powered_positive_wall_loss,
            "grounded_positive_wall_loss": grounded_positive_wall_loss,
            "neutral_wall_loss": neutral_wall_loss,
            "neutral_thermal_flux": neutral_thermal_flux,
            "axial_flux": axial_flux,
            "parent_power": parent_power,
            "supplemental_power": supplemental_power,
            "charged_wall_power": charged_wall_power,
            "powered_charged_wall_power": powered_charged_wall_power,
            "grounded_charged_wall_power": grounded_charged_wall_power,
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
        nonlinear_verbose: int = 0,
    ) -> ZhuOpenReactorSolution:
        if not isinstance(condition, ZhuOpenReactorCondition):
            raise TypeError("a Zhu open-reactor condition is required")
        if int(nonlinear_verbose) not in {0, 1, 2}:
            raise ValueError("nonlinear_verbose must be 0, 1, or 2")
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
            # A dissociating fluorocarbon discharge does not retain a nearly
            # frozen feed inventory.  Reserve most of the neutral pressure for
            # daughters so the cold start lies inside their reaction basin;
            # this is a composition-neutral numerical seed, not a fitted
            # steady-state fraction.
            neutral_seed_fraction = 0.40
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
            # The CHF3/SF6 feed is strongly electronegative.  Start the log
            # solve on a charge-neutral interior state at a few parts per
            # thousand ionization instead of a nearly unionized boundary,
            # where derivatives of rare charged daughters vanish.
            charged_seed = 5.0e-3 * target_density
            for name in self.positive_names:
                densities[name] = charged_seed / len(self.positive_names)
            for name in self.negative_names:
                densities[name] = 0.99 * charged_seed / len(self.negative_names)
            densities["e"] = 0.01 * charged_seed
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
        supplemental_coefficient_cache: dict[
            tuple[float, float], np.ndarray
        ] = {}
        # These are exact condition-only closures.  Solving the cylindrical
        # diffusion eigenproblem inside every nonlinear residual evaluation
        # made a 0-D solve needlessly scale with the Jacobian color count.
        neutral_wall_frequency = self._neutral_wall_frequencies(condition)

        def residual(log_state):
            physical = self._ledger(
                log_state,
                condition=condition,
                maximum_tail_population_fraction=(
                    maximum_tail_population_fraction),
                neutral_wall_frequency=neutral_wall_frequency,
                supplemental_coefficient_cache=(
                    supplemental_coefficient_cache),
                electron_cache=electron_cache,
                electron_continuation=electron_continuation,
            )["residual"]
            transformed = physical.copy()
            transformed[:len(self.heavy_order)] = log_flux_ratio_residual(
                physical[:len(self.heavy_order)]
            )
            return transformed

        result = least_squares(
            residual,
            x0=np.clip(np.log(initial), lower, upper),
            bounds=(lower, upper),
            xtol=2.0e-10,
            ftol=2.0e-10,
            gtol=2.0e-10,
            max_nfev=int(maximum_evaluations),
            jac_sparsity=self.jacobian_sparsity(),
            verbose=int(nonlinear_verbose),
        )
        ledger = self._ledger(
            result.x,
            condition=condition,
            maximum_tail_population_fraction=maximum_tail_population_fraction,
            neutral_wall_frequency=neutral_wall_frequency,
            supplemental_coefficient_cache=supplemental_coefficient_cache,
            electron_cache=electron_cache,
            electron_continuation=electron_continuation,
        )
        maximum_residual = float(np.max(np.abs(ledger["residual"])))
        if (
            result.status < 0
            or not math.isfinite(maximum_residual)
            or maximum_residual > residual_tolerance
        ):
            dominant = max(
                ledger["balances"],
                key=lambda name: abs(ledger["balances"][name]),
            )
            dominant_diagnostic = ""
            if dominant in ledger["densities"]:
                dominant_diagnostic = (
                    f", density={ledger['densities'][dominant]:.12g}, "
                    f"internal={ledger['internal_sources'][dominant]:.12g}, "
                    f"turnover={ledger['turnover'][dominant]:.12g}, "
                    "external="
                    f"{sum(ledger['external_terms'][dominant]):.12g}"
                )
            raise RuntimeError(
                "Zhu open-reactor solve failed conservation gate: "
                f"success={result.success}, residual={maximum_residual}, "
                f"dominant={dominant}={ledger['balances'][dominant]}, "
                f"E/N={ledger['reduced_field']} Td"
                f"{dominant_diagnostic}, message={result.message}"
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
            powered_electrode_positive_ion_wall_loss_m3_s=(
                ledger["powered_positive_wall_loss"]),
            grounded_positive_ion_wall_loss_m3_s=(
                ledger["grounded_positive_wall_loss"]),
            neutral_wall_loss_m3_s=ledger["neutral_wall_loss"],
            absorbed_power_density_W_m3=condition.absorbed_power_density_W_m3,
            parent_collision_power_density_W_m3=ledger["parent_power"],
            supplemental_collision_power_density_W_m3=(
                ledger["supplemental_power"]),
            charged_wall_power_density_W_m3=ledger["charged_wall_power"],
            powered_electrode_charged_wall_power_density_W_m3=(
                ledger["powered_charged_wall_power"]),
            grounded_charged_wall_power_density_W_m3=(
                ledger["grounded_charged_wall_power"]),
            normalized_residuals=ledger["balances"],
            electron_solution=ledger["electron_solution"],
            parent_collision_state=ledger["parent_state"],
            solver_evaluations=result.nfev,
            electron_collision_basis_neutral_fraction=(
                sum(
                    ledger["densities"][name]
                    for name in self.electron_collision_targets
                )
                / sum(
                    ledger["densities"][name]
                    for name in self.neutral_names
                )
            ),
        )
