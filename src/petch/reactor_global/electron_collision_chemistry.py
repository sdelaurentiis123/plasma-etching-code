"""Map deterministic electron-collision moments into reactor particle ledgers.

Collision-deck product labels are descriptive strings, not a reaction
mechanism.  This module therefore requires an explicit, deck-indexed mapping
for every non-momentum row and verifies atoms, charge, electron-number change,
and the deck hash before any volume source is emitted.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .electron_collision_deck import ElectronCollisionDeck
from .electron_kinetics import (
    TwoTermBoltzmannCondition,
    TwoTermBoltzmannSolution,
)
from .network import E_CHARGE_C, Species


_MOMENTUM_KINDS = frozenset({"EFFECTIVE", "ELASTIC", "MOMENTUM"})
def _integer_stoichiometry(
    values: Mapping[str, int],
    *,
    name: str,
) -> MappingProxyType:
    converted = {}
    for species, raw_value in values.items():
        value = int(raw_value)
        if raw_value != value:
            raise ValueError(f"invalid {name} stoichiometry")
        converted[str(species)] = value
    if (
        not converted
        or any(not species.strip() for species in converted)
        or any(value <= 0 for value in converted.values())
    ):
        raise ValueError(f"invalid {name} stoichiometry")
    return MappingProxyType(converted)


@dataclass(frozen=True)
class ElectronCollisionHeavyMapping:
    """One explicit heavy-species interpretation of a collision-deck row."""

    process_index: int
    reaction_name: str
    heavy_reactants: Mapping[str, int]
    heavy_products: Mapping[str, int]
    source: str
    evidence_kind: str

    def __post_init__(self):
        if (
            int(self.process_index) != self.process_index
            or self.process_index < 0
            or not str(self.reaction_name).strip()
            or not str(self.source).strip()
            or not str(self.evidence_kind).strip()
        ):
            raise ValueError("invalid electron-collision heavy mapping")
        object.__setattr__(self, "process_index", int(self.process_index))
        object.__setattr__(
            self,
            "heavy_reactants",
            _integer_stoichiometry(
                self.heavy_reactants, name="heavy-reactant"),
        )
        object.__setattr__(
            self,
            "heavy_products",
            _integer_stoichiometry(
                self.heavy_products, name="heavy-product"),
        )


@dataclass(frozen=True)
class ElectronCollisionChemistryState:
    """Volume rates and exact particle/power closure for one EEPF state."""

    event_rates_m3_s: Mapping[str, float]
    species_sources_m3_s: Mapping[str, float]
    collision_growth_source_m3_s: float
    eepf_growth_source_m3_s: float
    electron_growth_closure_error_m3_s: float
    collisional_field_power_gain_W_m3: float
    collision_deck_sha256: str
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        rates = {str(name): float(value) for name, value in self.event_rates_m3_s.items()}
        sources = {
            str(name): float(value)
            for name, value in self.species_sources_m3_s.items()
        }
        scalars = np.asarray((
            self.collision_growth_source_m3_s,
            self.eepf_growth_source_m3_s,
            self.electron_growth_closure_error_m3_s,
            self.collisional_field_power_gain_W_m3,
        ))
        if (
            not rates
            or not sources
            or any(not math.isfinite(value) or value < 0.0 for value in rates.values())
            or any(not math.isfinite(value) for value in sources.values())
            or np.any(~np.isfinite(scalars))
            or self.collisional_field_power_gain_W_m3 < 0.0
            or len(self.collision_deck_sha256) != 64
        ):
            raise ValueError("invalid electron-collision chemistry state")
        object.__setattr__(self, "event_rates_m3_s", MappingProxyType(rates))
        object.__setattr__(self, "species_sources_m3_s", MappingProxyType(sources))
        for name, value in zip((
            "collision_growth_source_m3_s",
            "eepf_growth_source_m3_s",
            "electron_growth_closure_error_m3_s",
            "collisional_field_power_gain_W_m3",
        ), scalars):
            object.__setattr__(self, name, float(value))

    @property
    def relative_electron_growth_closure(self) -> float:
        return abs(self.electron_growth_closure_error_m3_s) / max(
            abs(self.collision_growth_source_m3_s),
            abs(self.eepf_growth_source_m3_s),
            1.0,
        )


class ElectronCollisionChemistry:
    """A fully covered, conservation-checked deck-to-reactor mapping."""

    def __init__(
        self,
        collision_deck: ElectronCollisionDeck,
        species: tuple[Species, ...],
        mappings: tuple[ElectronCollisionHeavyMapping, ...],
    ):
        if not isinstance(collision_deck, ElectronCollisionDeck):
            raise TypeError("an electron collision deck is required")
        species = tuple(species)
        mappings = tuple(mappings)
        if (
            not species
            or any(not isinstance(item, Species) for item in species)
            or len({item.name for item in species}) != len(species)
            or any(
                not isinstance(item, ElectronCollisionHeavyMapping)
                for item in mappings
            )
        ):
            raise ValueError("invalid collision-chemistry species or mappings")
        required = {
            index
            for index, process in enumerate(collision_deck.processes)
            if process.kind not in _MOMENTUM_KINDS
        }
        supplied = {item.process_index for item in mappings}
        if len(supplied) != len(mappings) or supplied != required:
            raise ValueError(
                "every non-momentum collision process requires one mapping"
            )
        self.collision_deck = collision_deck
        self.species = species
        self.mappings = tuple(sorted(
            mappings, key=lambda item: item.process_index))
        self._species = {item.name: item for item in species}
        for mapping in self.mappings:
            self._validate_mapping(mapping)

    def _validate_mapping(self, mapping: ElectronCollisionHeavyMapping) -> None:
        process = self.collision_deck.processes[mapping.process_index]
        names = set(mapping.heavy_reactants) | set(mapping.heavy_products)
        if names - set(self._species):
            raise ValueError(
                f"mapping {mapping.reaction_name} uses unknown heavy species"
            )
        if mapping.heavy_reactants.get(process.target) != 1:
            raise ValueError(
                f"mapping {mapping.reaction_name} must consume its collision target"
            )
        elements = sorted({
            element
            for name in names
            for element in self._species[name].composition
        })
        for element in elements:
            reactant_count = sum(
                coefficient * self._species[name].composition.get(element, 0)
                for name, coefficient in mapping.heavy_reactants.items()
            )
            product_count = sum(
                coefficient * self._species[name].composition.get(element, 0)
                for name, coefficient in mapping.heavy_products.items()
            )
            if reactant_count != product_count:
                raise ValueError(
                    f"mapping {mapping.reaction_name} does not conserve {element}"
                )
        heavy_charge_change = sum(
            coefficient * self._species[name].charge_number
            for name, coefficient in mapping.heavy_products.items()
        ) - sum(
            coefficient * self._species[name].charge_number
            for name, coefficient in mapping.heavy_reactants.items()
        )
        electron_change = process.electron_number_change
        if heavy_charge_change != electron_change:
            raise ValueError(
                f"mapping {mapping.reaction_name} does not conserve charge"
            )

    def evaluate(
        self,
        solution: TwoTermBoltzmannSolution,
        condition: TwoTermBoltzmannCondition,
        densities_m3: Mapping[str, float],
        *,
        closure_relative_tolerance: float = 1.0e-8,
    ) -> ElectronCollisionChemistryState:
        if not isinstance(solution, TwoTermBoltzmannSolution):
            raise TypeError("a two-term Boltzmann solution is required")
        if not isinstance(condition, TwoTermBoltzmannCondition):
            raise TypeError("a two-term Boltzmann condition is required")
        if solution.collision_deck_sha256 != self.collision_deck.payload_sha256:
            raise ValueError("Boltzmann solution and chemistry deck hashes differ")
        if (
            solution.reduced_electric_field_Td
            != condition.reduced_electric_field_Td
            or solution.gas_temperature_K != condition.gas_temperature_K
        ):
            raise ValueError("Boltzmann solution and condition differ")
        densities = {str(name): float(value) for name, value in densities_m3.items()}
        if (
            set(densities) != set(self._species)
            or any(not math.isfinite(value) or value < 0.0 for value in densities.values())
            or densities.get("e", 0.0) <= 0.0
        ):
            raise ValueError("invalid collision-chemistry densities")
        target_density = sum(densities[name] for name in self.collision_deck.targets)
        if target_density <= 0.0:
            raise ValueError("collision targets have zero total density")
        expected_fractions = {
            name: densities[name] / target_density
            for name in self.collision_deck.targets
        }
        if any(
            not math.isclose(
                condition.target_mole_fractions[name], expected_fractions[name],
                rel_tol=2.0e-13, abs_tol=2.0e-15,
            )
            for name in expected_fractions
        ):
            raise ValueError(
                "Boltzmann target fractions do not match reactor densities"
            )

        event_rates = {}
        sources = {name: 0.0 for name in self._species}
        for mapping in self.mappings:
            moment = solution.collision_moments[mapping.process_index]
            process = self.collision_deck.processes[mapping.process_index]
            rate = (
                densities["e"]
                * densities[process.target]
                * float(moment.rate_coefficient_m3_s)
            )
            event_rates[mapping.reaction_name] = rate
            for name, coefficient in mapping.heavy_reactants.items():
                sources[name] -= coefficient * rate
            for name, coefficient in mapping.heavy_products.items():
                sources[name] += coefficient * rate
            sources["e"] += process.electron_number_change * rate

        eepf_growth = (
            densities["e"]
            * target_density
            * solution.net_growth_rate_coefficient_m3_s
        )
        collision_growth = sources["e"]
        closure_error = collision_growth - eepf_growth
        scale = max(abs(collision_growth), abs(eepf_growth), 1.0)
        tolerance = float(closure_relative_tolerance)
        if (
            not math.isfinite(tolerance)
            or tolerance <= 0.0
            or abs(closure_error) > tolerance * scale
        ):
            raise RuntimeError("collision chemistry failed electron-growth closure")
        return ElectronCollisionChemistryState(
            event_rates_m3_s=event_rates,
            species_sources_m3_s=sources,
            collision_growth_source_m3_s=collision_growth,
            eepf_growth_source_m3_s=eepf_growth,
            electron_growth_closure_error_m3_s=closure_error,
            collisional_field_power_gain_W_m3=(
                E_CHARGE_C
                * densities["e"]
                * target_density
                * solution.transport_moments.reduced_field_power_gain_eV_m3_s
            ),
            collision_deck_sha256=self.collision_deck.payload_sha256,
        )
