"""Conservative material ledger at the feature/reactor boundary.

The ledger conserves declared material-origin units without pretending that an unknown reactive product
branch is known. A SiO2 mechanism can report every removed formula unit while leaving its allocation
among volatile and condensed products unresolved. A physical-sputtering mechanism may instead route the
removed solid units directly into a transported redeposit population.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np


def _immutable_inventory(values, *, field_name):
    inventory = {}
    for name, value in dict(values).items():
        array = np.asarray(value, dtype=float).copy()
        if not isinstance(name, str) or not name or np.any(~np.isfinite(array)) or np.any(array < 0.0):
            raise ValueError(f"invalid {field_name} inventory")
        array.setflags(write=False); inventory[name] = array
    return MappingProxyType(inventory)


@dataclass(frozen=True)
class SurfaceMaterialExchange:
    """Integrated material exchange over one surface-kinetics step.

    Values are nonnegative material-origin units per square metre. ``removed`` must equal ``outgoing``
    plus ``unresolved`` for every inventory. ``outgoing`` becomes transportable only when a separate
    population supplies physical species identity, mass, energy/angle distribution, and interaction law.
    ``deposited`` is material added from incident species and is outside the removal identity.
    """

    removed_units_m2: Mapping[str, np.ndarray]
    outgoing_units_m2: Mapping[str, np.ndarray]
    unresolved_units_m2: Mapping[str, np.ndarray]
    deposited_units_m2: Mapping[str, np.ndarray]
    known_limitations: tuple[str, ...] = ()

    def __post_init__(self):
        removed = _immutable_inventory(self.removed_units_m2, field_name="removed")
        outgoing = _immutable_inventory(self.outgoing_units_m2, field_name="outgoing")
        unresolved = _immutable_inventory(self.unresolved_units_m2, field_name="unresolved")
        deposited = _immutable_inventory(self.deposited_units_m2, field_name="deposited")
        if set(outgoing) - set(removed) or set(unresolved) - set(removed):
            raise ValueError("outgoing and unresolved inventories must originate in removed material")
        for name, source in removed.items():
            emitted = np.asarray(outgoing.get(name, 0.0), dtype=float)
            unknown = np.asarray(unresolved.get(name, 0.0), dtype=float)
            try:
                source_view, emitted, unknown = np.broadcast_arrays(source, emitted, unknown)
            except ValueError as error:
                raise ValueError(f"material exchange shape mismatch for {name}") from error
            scale = np.maximum(source_view, 1.0)
            if np.any(np.abs(source_view - emitted - unknown)
                      > 64.0 * np.finfo(float).eps * scale):
                raise ValueError(f"removed material does not close for {name}")
        object.__setattr__(self, "removed_units_m2", removed)
        object.__setattr__(self, "outgoing_units_m2", outgoing)
        object.__setattr__(self, "unresolved_units_m2", unresolved)
        object.__setattr__(self, "deposited_units_m2", deposited)
        object.__setattr__(self, "known_limitations", tuple(self.known_limitations))

    @property
    def product_routing_complete(self):
        return all(not np.any(value > 0.0) for value in self.unresolved_units_m2.values())

    def residual_units_m2(self, name):
        if name not in self.removed_units_m2:
            raise KeyError(name)
        return (self.removed_units_m2[name]
                - np.asarray(self.outgoing_units_m2.get(name, 0.0))
                - np.asarray(self.unresolved_units_m2.get(name, 0.0)))


@dataclass(frozen=True)
class SurfaceProductPopulation:
    """One explicit population emitted from a surface step.

    ``integrated_particle_count_m2`` is particle count per emitting face area over the step.
    ``material_units_per_particle`` maps those particles back to the conserved material-origin ledger.
    Material allocation may be known before the emission energy/angular law is. Such a population closes
    the ledger but is not transport-ready; transport must refuse it rather than silently substitute a
    cosine or Maxwellian distribution.
    """

    name: str
    source_inventory: str
    integrated_particle_count_m2: np.ndarray
    material_units_per_particle: float
    mass_amu: float
    angular_model: str | None = None
    energy_model: str | None = None
    energy_parameters: Mapping[str, object] = field(default_factory=dict)
    provenance: Mapping[str, object] = field(default_factory=dict)
    relative_standard_uncertainty: float | None = None

    def __post_init__(self):
        count = np.asarray(self.integrated_particle_count_m2, dtype=float).copy()
        energy_parameters = dict(self.energy_parameters)
        provenance = dict(self.provenance)
        if (not self.name or not self.source_inventory
                or np.any(~np.isfinite(count)) or np.any(count < 0.0)
                or not np.isfinite(self.material_units_per_particle)
                or self.material_units_per_particle <= 0.0
                or not np.isfinite(self.mass_amu) or self.mass_amu <= 0.0
                or (self.angular_model is not None and not self.angular_model)
                or (self.energy_model is not None and not self.energy_model)):
            raise ValueError("invalid emitted surface-product population")
        if self.energy_model is None and energy_parameters:
            raise ValueError("energy parameters require an energy model")
        if self.energy_model == "monoenergetic":
            if (set(energy_parameters) != {"energy_eV"}
                    or not np.isfinite(energy_parameters["energy_eV"])
                    or energy_parameters["energy_eV"] < 0.0):
                raise ValueError("monoenergetic products require nonnegative energy_eV")
        elif self.energy_model == "thompson":
            if set(energy_parameters) != {
                    "surface_binding_energy_eV", "maximum_energy_eV"}:
                raise ValueError("Thompson products require binding and maximum energies")
            binding = energy_parameters["surface_binding_energy_eV"]
            maximum = energy_parameters["maximum_energy_eV"]
            if (not np.isfinite(binding) or binding <= 0.0
                    or not np.isfinite(maximum) or maximum <= binding):
                raise ValueError("invalid Thompson emission energies")
        elif self.energy_model is not None:
            raise ValueError(f"unsupported surface-product energy model: {self.energy_model}")
        uncertainty = self.relative_standard_uncertainty
        if uncertainty is not None and (not np.isfinite(uncertainty) or uncertainty < 0.0):
            raise ValueError("product-population uncertainty must be finite and nonnegative")
        count.setflags(write=False)
        object.__setattr__(self, "integrated_particle_count_m2", count)
        object.__setattr__(self, "energy_parameters", MappingProxyType(energy_parameters))
        object.__setattr__(self, "provenance", MappingProxyType(provenance))

    @property
    def integrated_material_units_m2(self):
        return self.integrated_particle_count_m2 * self.material_units_per_particle

    @property
    def transport_ready(self):
        return self.angular_model is not None and self.energy_model is not None


def validate_surface_product_routing(exchange, populations):
    """Require explicit populations to reproduce every outgoing ledger inventory face-by-face."""
    if not isinstance(exchange, SurfaceMaterialExchange):
        raise TypeError("exchange must be SurfaceMaterialExchange")
    populations = tuple(populations)
    if any(not isinstance(item, SurfaceProductPopulation) for item in populations):
        raise TypeError("surface products must be SurfaceProductPopulation objects")
    if len({item.name for item in populations}) != len(populations):
        raise ValueError("surface-product names must be unique")
    routed = {}
    for population in populations:
        value = population.integrated_material_units_m2
        routed[population.source_inventory] = routed.get(population.source_inventory, 0.0) + value
    if set(routed) != set(exchange.outgoing_units_m2):
        raise ValueError("surface-product inventories do not match outgoing material ledger")
    for name, expected in exchange.outgoing_units_m2.items():
        try:
            actual, expected = np.broadcast_arrays(routed[name], expected)
        except ValueError as error:
            raise ValueError(f"surface-product routing shape mismatch for {name}") from error
        scale = np.maximum(expected, 1.0)
        if np.any(np.abs(actual - expected) > 64.0 * np.finfo(float).eps * scale):
            raise ValueError(f"surface-product routing does not close for {name}")
    return populations


def unresolved_surface_exchange(*, removed_units_m2, deposited_units_m2=(), limitations=()):
    """Construct a closed ledger when removal is known but product routing is not."""
    removed = dict(removed_units_m2)
    return SurfaceMaterialExchange(
        removed_units_m2=removed, outgoing_units_m2={}, unresolved_units_m2=removed,
        deposited_units_m2=dict(deposited_units_m2), known_limitations=tuple(limitations))
