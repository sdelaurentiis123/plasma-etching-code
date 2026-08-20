"""Conserved Bosch silicon etching with an explicit fluorocarbon-film clock.

Bosch cycling alternates an ``SF6`` silicon-removal phase with a ``C4F8``
passivation phase.  This module composes two existing common-engine laws rather
than introducing an empirical depth multiplier:

* :class:`BelenSiliconSF6O2Mechanism` supplies the bare-silicon reaction rate;
* :class:`LaMagnaGarozzoFluorocarbonMechanism` supplies finite-film deposition,
  ion-assisted film removal, and the exact fraction of a step for which the
  substrate is exposed after that film is depleted.

The La Magna substrate-removal rate is deliberately discarded.  It is a
SiO2 law and is used here only as a conserved fluorocarbon-film clock.  Bare Si
removal is the unchanged Belen rate multiplied by the film model's exposed
time fraction.  Film and silicon inventories close independently in the
material ledger, and both use the common feature engine's conservative remap.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .fluorocarbon_lamagna import (
    LaMagnaFluorocarbonState, LaMagnaGarozzoFluorocarbonMechanism,
)
from .silicon_sf6o2 import BelenSiliconSF6O2Mechanism, BelenSiliconState
from .surface_exchange import SurfaceMaterialExchange, unresolved_surface_exchange
from .surface_kinetics import MechanismValidity, SurfaceFluxes


@dataclass(frozen=True)
class BoschSiliconFluorocarbonState:
    """Intensive surface states plus conserved Si and film inventories."""

    available_site_fraction: np.ndarray | float = 1.0
    removed_si_atoms_m2: np.ndarray | float = 0.0
    etchant_coverage: np.ndarray | float = 0.0
    polymer_coverage: np.ndarray | float = 0.0
    etchant_on_polymer_coverage: np.ndarray | float = 0.0
    polymer_film_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        values = np.broadcast_arrays(*(
            np.asarray(getattr(self, name), dtype=float)
            for name in self.__dataclass_fields__))
        output = [np.array(value, copy=True) for value in values]
        available, removed, etchant, polymer, etchant_polymer, film = output
        bounded = (available, etchant, polymer, etchant_polymer)
        if (any(np.any(~np.isfinite(value)) for value in output)
                or any(np.any((value < 0.0) | (value > 1.0)) for value in bounded)
                or np.any(removed < 0.0) or np.any(film < 0.0)):
            raise ValueError("invalid Bosch silicon/fluorocarbon surface state")
        for name, value in zip(self.__dataclass_fields__, output):
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape)
        return cls(np.ones(shape), zero, zero, zero, zero, zero)

    def conservative_surface_fields(self):
        return {
            name: getattr(self, name) for name in self.__dataclass_fields__
        }

    def conservative_surface_upper_bounds(self):
        return {
            "available_site_fraction": 1.0,
            "removed_si_atoms_m2": None,
            "etchant_coverage": 1.0,
            "polymer_coverage": 1.0,
            "etchant_on_polymer_coverage": 1.0,
            "polymer_film_units_m2": None,
        }

    def surface_field_remap_modes(self):
        return {
            "available_site_fraction": "intensive",
            "removed_si_atoms_m2": "conservative",
            "etchant_coverage": "intensive",
            "polymer_coverage": "intensive",
            "etchant_on_polymer_coverage": "intensive",
            "polymer_film_units_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(self.conservative_surface_fields()):
            raise ValueError("Bosch surface remap fields do not match its state contract")
        return type(self)(**fields)

    def silicon_state(self):
        return BelenSiliconState(
            self.available_site_fraction, self.removed_si_atoms_m2)

    def film_state(self):
        return LaMagnaFluorocarbonState(
            self.etchant_coverage, self.polymer_coverage,
            self.etchant_on_polymer_coverage, self.polymer_film_units_m2,
            np.zeros_like(self.polymer_film_units_m2))


@dataclass(frozen=True)
class BoschSiliconFluorocarbonStepResult:
    state: BoschSiliconFluorocarbonState
    etch_velocity_m_s: np.ndarray
    normal_growth_velocity_m_s: np.ndarray
    substrate_exposure_fraction: np.ndarray
    silicon_removal_rate_m2_s: np.ndarray
    polymer_deposition_rate_m2_s: np.ndarray
    polymer_removal_rate_m2_s: np.ndarray
    removed_si_atoms_m2: np.ndarray
    deposited_polymer_units_m2: np.ndarray
    removed_polymer_units_m2: np.ndarray
    fluorine_coverage: np.ndarray
    oxygen_coverage: np.ndarray
    available_site_fraction: np.ndarray
    etchant_coverage: np.ndarray
    polymer_coverage: np.ndarray
    etchant_on_polymer_coverage: np.ndarray
    transport_fixed_point_change: np.ndarray
    material_exchange: SurfaceMaterialExchange
    validity: MechanismValidity
    product_populations: tuple = ()

    def __post_init__(self):
        if not isinstance(self.state, BoschSiliconFluorocarbonState):
            raise TypeError("invalid Bosch surface result state")
        nonnegative = set(self.__dataclass_fields__) - {
            "state", "material_exchange", "validity", "product_populations",
            "transport_fixed_point_change",
        }
        for name in nonnegative | {"transport_fixed_point_change"}:
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if (np.any(~np.isfinite(value))
                    or (name in nonnegative and np.any(value < 0.0))):
                raise ValueError(f"invalid Bosch result field: {name}")
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if np.any(self.substrate_exposure_fraction > 1.0):
            raise ValueError("substrate exposure fraction must lie in [0, 1]")
        if not isinstance(self.material_exchange, SurfaceMaterialExchange):
            raise TypeError("Bosch result requires a material-exchange ledger")
        if not isinstance(self.validity, MechanismValidity):
            raise TypeError("Bosch result requires mechanism validity")
        object.__setattr__(self, "product_populations", tuple(self.product_populations))
        if self.product_populations:
            raise ValueError("Bosch v1 leaves volatile products explicitly unresolved")


def _filtered_fluxes(fluxes, *, neutral_species, projectile_species):
    neutral_species = set(neutral_species)
    projectile_species = set(projectile_species)
    return SurfaceFluxes(
        {name: value for name, value in fluxes.neutral_flux_m2_s.items()
         if name in neutral_species},
        tuple(population for population in fluxes.energetic_fluxes
              if population.name in projectile_species))


class BoschSiliconFluorocarbonMechanism:
    """Belen bare-Si removal gated by a finite La Magna fluorocarbon film."""

    quasi_steady_surface_state = True

    def __init__(
            self, silicon_mechanism: BelenSiliconSF6O2Mechanism,
            film_mechanism: LaMagnaGarozzoFluorocarbonMechanism):
        if not isinstance(silicon_mechanism, BelenSiliconSF6O2Mechanism):
            raise TypeError("silicon_mechanism must be BelenSiliconSF6O2Mechanism")
        if not isinstance(film_mechanism, LaMagnaGarozzoFluorocarbonMechanism):
            raise TypeError(
                "film_mechanism must be LaMagnaGarozzoFluorocarbonMechanism")
        silicon = silicon_mechanism.parameters
        film = film_mechanism.parameters
        if silicon.fluorine_species not in film.etchant_species:
            raise ValueError("the Belen fluorine species must etch the film")
        if set(silicon.projectile_species) != set(film.projectile_species):
            raise ValueError("silicon and film laws must share projectile species")
        if set(film.polymer_species) & {
                silicon.fluorine_species, silicon.oxygen_species}:
            raise ValueError("film precursor species must be distinct from Belen F/O")
        self.silicon_mechanism = silicon_mechanism
        self.film_mechanism = film_mechanism
        self.provenance = MappingProxyType({
            "model": "bosch-silicon-fluorocarbon-composite-v1",
            "composition_rule": (
                "unchanged Belen bare-Si rate times La Magna finite-film "
                "substrate exposure fraction"),
            "discarded_quantity": (
                "La Magna SiO2 substrate-removal rate; retained only as the "
                "film-depletion exposure clock"),
            "silicon_mechanism": dict(silicon_mechanism.provenance),
            "film_mechanism": dict(film_mechanism.provenance),
            "known_omissions": [
                "Bosch-cycle transients enter through the caller's measured waveform",
                "fluorocarbon film composition, crosslinking, and carbonization are reduced",
                "oxygen attack on an existing fluorocarbon film is not represented",
                "volatile SiFx and fluorocarbon product return transport are unresolved",
            ],
        })

    @staticmethod
    def initial_state(shape=()):
        return BoschSiliconFluorocarbonState.bare(shape)

    @property
    def _supported_neutrals(self):
        silicon = self.silicon_mechanism.parameters
        film = self.film_mechanism.parameters
        return ({silicon.fluorine_species, silicon.oxygen_species}
                | set(film.etchant_species) | set(film.polymer_species))

    @property
    def _supported_projectiles(self):
        return set(self.silicon_mechanism.parameters.projectile_species)

    def validity(self, fluxes: SurfaceFluxes):
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in self._supported_neutrals
            and np.any(np.asarray(value) > 0.0)))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in self._supported_projectiles
            and np.any(np.asarray(population.flux_m2_s) > 0.0)}))
        reasons = []
        if unsupported_neutral:
            reasons.append(
                "positive neutral flux has no declared Bosch reaction channel")
        if unsupported_energetic:
            reasons.append(
                "positive energetic flux has no declared Bosch channel: "
                + ", ".join(unsupported_energetic))

        silicon = self.silicon_mechanism.parameters
        film = self.film_mechanism.parameters
        silicon_validity = self.silicon_mechanism.validity(_filtered_fluxes(
            fluxes,
            neutral_species={silicon.fluorine_species, silicon.oxygen_species},
            projectile_species=silicon.projectile_species))
        film_validity = self.film_mechanism.validity(_filtered_fluxes(
            fluxes,
            neutral_species=set(film.etchant_species) | set(film.polymer_species),
            projectile_species=film.projectile_species))
        reasons.extend(silicon_validity.reasons)
        reasons.extend(film_validity.reasons)
        nonpredictive = tuple(
            f"silicon.{name}" for name in silicon_validity.nonpredictive_parameters
        ) + tuple(
            f"film.{name}" for name in film_validity.nonpredictive_parameters)
        omissions = tuple(dict.fromkeys(
            (*silicon_validity.known_model_form_omissions,
             *film_validity.known_model_form_omissions,
             *self.provenance["known_omissions"])))
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=omissions,
            parameter_evidence_supports_prediction=not nonpredictive,
            nonpredictive_parameters=nonpredictive)

    def neutral_reaction_probability(
            self, state: BoschSiliconFluorocarbonState):
        if not isinstance(state, BoschSiliconFluorocarbonState):
            raise TypeError("Bosch neutral probabilities require Bosch state")
        silicon_probability = dict(
            self.silicon_mechanism.neutral_reaction_probability(
                state.silicon_state()))
        film_probability = dict(
            self.film_mechanism.neutral_reaction_probability(state.film_state()))
        film_present = state.polymer_film_units_m2 > 0.0
        output = {}
        for name in sorted(set(silicon_probability) | set(film_probability)):
            if name in silicon_probability and name in film_probability:
                output[name] = np.where(
                    film_present, film_probability[name], silicon_probability[name])
            elif name in silicon_probability:
                output[name] = np.where(
                    film_present, 0.0, silicon_probability[name])
            else:
                output[name] = film_probability[name]
        return MappingProxyType(output)

    def advance(
            self, state: BoschSiliconFluorocarbonState, fluxes: SurfaceFluxes,
            duration_s: float, *, strict=True):
        if not isinstance(state, BoschSiliconFluorocarbonState):
            raise TypeError("Bosch advance requires Bosch state")
        if not np.isfinite(duration_s) or duration_s < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(fluxes)
        if strict and not validity.within_declared_scope:
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))

        silicon = self.silicon_mechanism.parameters
        film = self.film_mechanism.parameters
        silicon_fluxes = _filtered_fluxes(
            fluxes,
            neutral_species={silicon.fluorine_species, silicon.oxygen_species},
            projectile_species=silicon.projectile_species)
        film_fluxes = _filtered_fluxes(
            fluxes,
            neutral_species=set(film.etchant_species) | set(film.polymer_species),
            projectile_species=film.projectile_species)
        film_result = self.film_mechanism.advance(
            state.film_state(), film_fluxes, duration_s, strict=False)
        silicon_result = self.silicon_mechanism.advance(
            state.silicon_state(), silicon_fluxes, duration_s, strict=False)

        exposure = film_result.substrate_exposure_fraction
        bare_si_rate = (
            silicon_result.chemical_removal_rate_m2_s
            + silicon_result.physical_sputter_rate_m2_s
            + silicon_result.ion_enhanced_removal_rate_m2_s)
        silicon_rate = bare_si_rate * exposure
        dt = float(duration_s)
        removed_si = silicon_rate * dt
        removed_polymer = film_result.removed_polymer_units_m2
        deposited_polymer = film_result.deposited_polymer_units_m2

        end_exposed = film_result.state.polymer_film_units_m2 <= 0.0
        available = np.where(
            end_exposed, silicon_result.available_site_fraction,
            state.available_site_fraction)
        updated = BoschSiliconFluorocarbonState(
            available,
            state.removed_si_atoms_m2 + removed_si,
            film_result.state.etchant_coverage,
            film_result.state.polymer_coverage,
            film_result.state.etchant_on_polymer_coverage,
            film_result.state.polymer_film_units_m2)

        if dt > 0.0:
            etch_velocity = (
                removed_si / silicon.bulk_si_atom_density_m3
                + removed_polymer / film.polymer_unit_density_m3) / dt
            growth_velocity = (
                deposited_polymer / film.polymer_unit_density_m3) / dt
        else:
            etch_velocity = (
                silicon_rate / silicon.bulk_si_atom_density_m3
                + film_result.polymer_removal_rate_m2_s
                / film.polymer_unit_density_m3)
            growth_velocity = (
                film_result.polymer_deposition_rate_m2_s
                / film.polymer_unit_density_m3)

        fixed_point_change = np.maximum.reduce((
            np.abs(updated.available_site_fraction - state.available_site_fraction),
            np.abs(updated.etchant_coverage - state.etchant_coverage),
            np.abs(updated.polymer_coverage - state.polymer_coverage),
            np.abs(updated.etchant_on_polymer_coverage
                   - state.etchant_on_polymer_coverage),
        ))
        exchange = unresolved_surface_exchange(
            removed_units_m2={
                silicon.material_inventory_name: removed_si,
                "fluorocarbon_film_unit": removed_polymer,
            },
            deposited_units_m2={"fluorocarbon_film_unit": deposited_polymer},
            limitations=(
                "volatile SiFx and fluorocarbon product branching is unresolved",
                "incident-site balances are reported by the composed mechanisms",
            ))
        return BoschSiliconFluorocarbonStepResult(
            state=updated,
            etch_velocity_m_s=etch_velocity,
            normal_growth_velocity_m_s=growth_velocity,
            substrate_exposure_fraction=exposure,
            silicon_removal_rate_m2_s=silicon_rate,
            polymer_deposition_rate_m2_s=film_result.polymer_deposition_rate_m2_s,
            polymer_removal_rate_m2_s=film_result.polymer_removal_rate_m2_s,
            removed_si_atoms_m2=removed_si,
            deposited_polymer_units_m2=deposited_polymer,
            removed_polymer_units_m2=removed_polymer,
            fluorine_coverage=silicon_result.fluorine_coverage,
            oxygen_coverage=silicon_result.oxygen_coverage,
            available_site_fraction=available,
            etchant_coverage=film_result.etchant_coverage,
            polymer_coverage=film_result.polymer_coverage,
            etchant_on_polymer_coverage=film_result.etchant_on_polymer_coverage,
            transport_fixed_point_change=fixed_point_change,
            material_exchange=exchange,
            validity=validity)
