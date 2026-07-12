"""Conservative material ledger at the feature/reactor boundary.

The ledger conserves declared material-origin units without pretending that an unknown reactive product
branch is known. A SiO2 mechanism can report every removed formula unit while leaving its allocation
among volatile and condensed products unresolved. A physical-sputtering mechanism may instead route the
removed solid units directly into a transported redeposit population.
"""
from __future__ import annotations

from dataclasses import dataclass
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


def unresolved_surface_exchange(*, removed_units_m2, deposited_units_m2=(), limitations=()):
    """Construct a closed ledger when removal is known but product routing is not."""
    removed = dict(removed_units_m2)
    return SurfaceMaterialExchange(
        removed_units_m2=removed, outgoing_units_m2={}, unresolved_units_m2=removed,
        deposited_units_m2=dict(deposited_units_m2), known_limitations=tuple(limitations))
