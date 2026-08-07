"""Conserved reactive-ion events from beam-resolved atomic evidence.

This module separates three things that a scalar etch yield collapses:

1. the incident molecular ion's C/F atom inventory;
2. mass-partitioned fragment energies and finite-range material transmission;
3. the measured identities of silicon-bearing desorption products.

Karahashi states that CFx+ dissociates into constituent atoms at impact and
that kinetic energy is distributed according to atomic mass.  Figure 10 then
gives SiF/SiF2/SiF4 fractions for CF3+ at three exact energies.  The Figure-10
ion incidence angle is unreported, however, so its product fractions cannot be
silently multiplied by the explicitly normal-incidence Figure-4 total yields.
This module refuses that join by default and exposes it only as a
non-production conditional projection.  Oxygen coproducts, carbon fate,
surface-fluorine exchange, delayed-product transport, and the reactor ion
mixture remain explicitly unresolved.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .ion_energy_deposition import (
    Target,
    residual_energy_after_layer_eV,
)
from .reactive_ion_beam import Karahashi2007ReactiveIonYieldTable


_ATOMS = MappingProxyType({
    "C": (6, 12.011),
    "F": (9, 18.998),
})
REACTIVE_ION_FORMULAE = MappingProxyType({
    "F+": MappingProxyType({"F": 1}),
    "CF+": MappingProxyType({"C": 1, "F": 1}),
    "CF2+": MappingProxyType({"C": 1, "F": 2}),
    "CF3+": MappingProxyType({"C": 1, "F": 3}),
})
_PRODUCT_F_COUNT = MappingProxyType({
    "SiF": 1,
    "SiF2": 2,
    "SiF4": 4,
})


def _immutable(value, name):
    array = np.asarray(value, dtype=float).copy()
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class ProjectileFragment:
    """One constituent atom type after molecular-ion dissociation."""

    element: str
    atomic_number: int
    mass_amu: float
    multiplicity: int
    energy_per_fragment_eV: np.ndarray

    def __post_init__(self):
        if (self.element not in _ATOMS
                or int(self.atomic_number) != self.atomic_number
                or self.atomic_number <= 0
                or not np.isfinite(self.mass_amu)
                or self.mass_amu <= 0.0
                or int(self.multiplicity) != self.multiplicity
                or self.multiplicity <= 0):
            raise ValueError("invalid projectile fragment")
        object.__setattr__(self, "atomic_number", int(self.atomic_number))
        object.__setattr__(self, "multiplicity", int(self.multiplicity))
        object.__setattr__(
            self, "energy_per_fragment_eV",
            _immutable(self.energy_per_fragment_eV, "fragment energy"))


def mass_partitioned_projectile_fragments(species, energy_eV):
    """Dissociate a supported CFx+ ion and conserve kinetic energy by mass."""
    formula = REACTIVE_ION_FORMULAE.get(species)
    if formula is None:
        raise ValueError(f"no resolved atomic formula for {species}")
    energy = _immutable(energy_eV, "incident ion energy")
    total_mass = sum(
        _ATOMS[element][1] * count for element, count in formula.items())
    fragments = tuple(
        ProjectileFragment(
            element=element,
            atomic_number=_ATOMS[element][0],
            mass_amu=_ATOMS[element][1],
            multiplicity=count,
            energy_per_fragment_eV=energy * _ATOMS[element][1] / total_mass,
        )
        for element, count in formula.items()
    )
    reconstructed = sum(
        item.multiplicity * item.energy_per_fragment_eV
        for item in fragments)
    tolerance = 16.0 * np.finfo(float).eps * np.maximum(energy, 1.0)
    if np.any(np.abs(reconstructed - energy) > tolerance):
        raise RuntimeError("mass-partitioned fragment energy does not close")
    return fragments


@dataclass(frozen=True)
class FragmentedProjectileTransmission:
    """CSDA transmission of every mass-partitioned constituent atom."""

    species: str
    incident_energy_eV: np.ndarray
    fragment_incident: tuple[ProjectileFragment, ...]
    residual_energy_per_fragment_eV: Mapping[str, np.ndarray]
    total_residual_energy_eV: np.ndarray
    total_deposited_energy_eV: np.ndarray

    def __post_init__(self):
        incident = _immutable(self.incident_energy_eV, "incident energy")
        residual = {
            name: _immutable(value, f"{name} residual energy")
            for name, value in dict(
                self.residual_energy_per_fragment_eV).items()
        }
        if set(residual) != {item.element for item in self.fragment_incident}:
            raise ValueError("fragment residual inventory is incomplete")
        total_residual = _immutable(
            self.total_residual_energy_eV, "total residual energy")
        total_deposited = _immutable(
            self.total_deposited_energy_eV, "total deposited energy")
        scale = np.maximum(incident, 1.0)
        if np.any(
            np.abs(incident - total_residual - total_deposited)
            > 32.0 * np.finfo(float).eps * scale
        ):
            raise ValueError("fragment transmission energy does not close")
        if np.any(
            total_residual > incident + 8.0 * np.finfo(float).eps * scale
        ):
            raise ValueError("fragment transmission creates energy")
        object.__setattr__(self, "incident_energy_eV", incident)
        object.__setattr__(
            self, "residual_energy_per_fragment_eV",
            MappingProxyType(residual))
        object.__setattr__(
            self, "total_residual_energy_eV", total_residual)
        object.__setattr__(
            self, "total_deposited_energy_eV", total_deposited)


def transmit_fragmented_projectile_through_layer(
        species, energy_eV, cosine_incidence, layer_depth_nm, target: Target):
    """Transmit a dissociated molecular ion through one material layer."""
    if not isinstance(target, Target):
        raise TypeError("target must be an ion-stopping Target")
    energy = np.asarray(energy_eV, dtype=float)
    cosine = np.asarray(cosine_incidence, dtype=float)
    depth = np.asarray(layer_depth_nm, dtype=float)
    energy, cosine, depth = np.broadcast_arrays(energy, cosine, depth)
    fragments = mass_partitioned_projectile_fragments(species, energy)
    residual = {}
    total = np.zeros(energy.shape, dtype=float)
    for fragment in fragments:
        value = residual_energy_after_layer_eV(
            fragment.energy_per_fragment_eV,
            cosine,
            depth,
            fragment.atomic_number,
            fragment.mass_amu,
            target,
        )
        residual[fragment.element] = value
        total = total + fragment.multiplicity * value
    return FragmentedProjectileTransmission(
        species=species,
        incident_energy_eV=energy,
        fragment_incident=fragments,
        residual_energy_per_fragment_eV=residual,
        total_residual_energy_eV=total,
        total_deposited_energy_eV=np.maximum(energy - total, 0.0),
    )


class Karahashi2007CF3ProductBranchTable:
    """Exact-setpoint Figure-10 SiFx branch fractions.

    No interpolation is offered: three product measurements are insufficient
    to assert a branch law between energies.
    """

    def __init__(self, observations, *, source_table_sha256):
        by_energy = {}
        row_count = 0
        for row in observations:
            row_count += 1
            if row.incident_species != "CF3+" or row.target != "SiO2":
                raise ValueError("Figure-10 product table is CF3+ on SiO2")
            if (row.ion_incidence_angle_deg is not None
                    or row.ion_incidence_angle_status
                    != "unreported_in_source"):
                raise ValueError(
                    "Figure-10 incidence angle must remain unresolved")
            energy = float(row.energy_eV)
            product = str(row.desorbed_product)
            by_energy.setdefault(energy, {})[product] = (
                float(row.product_fraction_percent_of_detected_sifx) / 100.0,
                float(row.digitization_uncertainty_percentage_points) / 100.0,
            )
        required_products = set(_PRODUCT_F_COUNT)
        if (row_count != 9
                or set(by_energy) != {500.0, 1000.0, 2000.0}
                or any(set(value) != required_products
                       for value in by_energy.values())):
            raise ValueError("incomplete Karahashi Figure-10 product board")
        if (len(source_table_sha256) != 64
                or any(char not in "0123456789abcdef"
                       for char in source_table_sha256)):
            raise ValueError("invalid Figure-10 source checksum")
        self.by_energy = MappingProxyType({
            energy: MappingProxyType(dict(value))
            for energy, value in sorted(by_energy.items())
        })
        self.source_table_sha256 = source_table_sha256
        digest = hashlib.sha256()
        digest.update(b"karahashi-2007-cf3-products-exact-v1")
        digest.update(source_table_sha256.encode("ascii"))
        for energy, products in self.by_energy.items():
            digest.update(np.float64(energy).tobytes())
            for product in sorted(products):
                digest.update(product.encode("ascii"))
                digest.update(np.asarray(products[product]).tobytes())
        self.fingerprint = digest.hexdigest()

    def evaluate_exact(self, species, energy_eV):
        energy = float(energy_eV)
        if species != "CF3+":
            raise ValueError("resolved Figure-10 products exist only for CF3+")
        if energy not in self.by_energy:
            raise ValueError(
                "Figure-10 products are available only at exact "
                "500, 1000, and 2000 eV setpoints")
        return self.by_energy[energy]


@dataclass(frozen=True)
class ReactiveIonEventOutcome:
    """A non-production conditional CF3+ Si/F product projection."""

    species: str
    energy_eV: float
    removed_sio2_formula_per_ion: float
    removed_yield_digitization_uncertainty: float
    incident_atoms_per_ion: Mapping[str, float]
    normalized_product_fraction: Mapping[str, float]
    conditional_sifx_particles_per_ion: Mapping[str, float]
    branch_normalization_factor: float
    required_f_atoms_per_ion: float
    required_f_atoms_lower: float
    required_f_atoms_upper: float
    unresolved_f_balance_per_ion: float
    unresolved_f_balance_lower: float
    unresolved_f_balance_upper: float
    condition_match_status: str
    production_eligible: bool
    unresolved_inventories: tuple[str, ...]

    def __post_init__(self):
        incident = MappingProxyType({
            name: float(value)
            for name, value in dict(self.incident_atoms_per_ion).items()
        })
        fraction = MappingProxyType({
            name: float(value)
            for name, value in dict(
                self.normalized_product_fraction).items()
        })
        particles = MappingProxyType({
            name: float(value)
            for name, value in dict(
                self.conditional_sifx_particles_per_ion).items()
        })
        scalars = (
            self.energy_eV,
            self.removed_sio2_formula_per_ion,
            self.removed_yield_digitization_uncertainty,
            self.branch_normalization_factor,
            self.required_f_atoms_per_ion,
            self.required_f_atoms_lower,
            self.required_f_atoms_upper,
            self.unresolved_f_balance_per_ion,
            self.unresolved_f_balance_lower,
            self.unresolved_f_balance_upper,
        )
        if (self.species != "CF3+"
                or any(not np.isfinite(value) for value in scalars)
                or any(value < 0.0 for value in scalars[:7])
                or set(incident) != {"C", "F"}
                or any(not np.isfinite(value) or value < 0.0
                       for value in incident.values())
                or set(fraction) != set(_PRODUCT_F_COUNT)
                or any(not np.isfinite(value) or value < 0.0
                       for value in fraction.values())
                or set(particles) != set(_PRODUCT_F_COUNT)
                or any(not np.isfinite(value) or value < 0.0
                       for value in particles.values())
                or abs(sum(fraction.values()) - 1.0)
                > 64.0 * np.finfo(float).eps
                or abs(
                    sum(particles.values())
                    - self.removed_sio2_formula_per_ion)
                > 64.0 * np.finfo(float).eps
                * max(self.removed_sio2_formula_per_ion, 1.0)
                or not (
                    self.required_f_atoms_lower
                    <= self.required_f_atoms_per_ion
                    <= self.required_f_atoms_upper)
                or not (
                    self.unresolved_f_balance_lower
                    <= self.unresolved_f_balance_per_ion
                    <= self.unresolved_f_balance_upper)
                or self.condition_match_status
                != "unresolved_figure10_incidence_angle"
                or self.production_eligible):
            raise ValueError("invalid reactive-ion event outcome")
        object.__setattr__(self, "incident_atoms_per_ion", incident)
        object.__setattr__(self, "normalized_product_fraction", fraction)
        object.__setattr__(
            self, "conditional_sifx_particles_per_ion", particles)
        object.__setattr__(
            self, "unresolved_inventories",
            tuple(self.unresolved_inventories))


class Karahashi2007ReactiveIonEventKernel:
    """Audit a hypothetical Figure-4/Figure-10 join without fitting."""

    def __init__(
            self,
            yield_table: Karahashi2007ReactiveIonYieldTable,
            product_table: Karahashi2007CF3ProductBranchTable):
        if not isinstance(yield_table, Karahashi2007ReactiveIonYieldTable):
            raise TypeError("yield_table must be Karahashi2007ReactiveIonYieldTable")
        if not isinstance(
                product_table, Karahashi2007CF3ProductBranchTable):
            raise TypeError(
                "product_table must be Karahashi2007CF3ProductBranchTable")
        self.yield_table = yield_table
        self.product_table = product_table
        digest = hashlib.sha256()
        digest.update(b"karahashi-2007-reactive-event-v1")
        digest.update(yield_table.fingerprint.encode("ascii"))
        digest.update(product_table.fingerprint.encode("ascii"))
        self.fingerprint = digest.hexdigest()

    def evaluate_conditional_cf3(
            self, energy_eV, *,
            acknowledge_unresolved_incidence_angle=False):
        if not acknowledge_unresolved_incidence_angle:
            raise ValueError(
                "Figure-10 incidence angle is unreported; the product "
                "fractions cannot be condition-matched to normal-incidence "
                "Figure-4 yields")
        energy = float(energy_eV)
        removal = float(self.yield_table.evaluate(
            "CF3+", energy, 1.0))
        removal_uncertainty = float(self.yield_table.evaluate_uncertainty(
            "CF3+", energy, 1.0))
        raw = self.product_table.evaluate_exact(
            "CF3+", energy)
        raw_fraction = {
            product: value[0] for product, value in raw.items()}
        raw_total = sum(raw_fraction.values())
        normalized = {
            product: value / raw_total
            for product, value in raw_fraction.items()
        }
        particles = {
            product: removal * value
            for product, value in normalized.items()
        }
        mean_f_per_si = sum(
            _PRODUCT_F_COUNT[product] * value
            for product, value in normalized.items())
        # Conservative L1 digitization envelope. The raw fraction errors are
        # propagated before normalization rather than treated as statistics.
        f_fraction_allowance = (
            sum(
                _PRODUCT_F_COUNT[product] * raw[product][1]
                for product in _PRODUCT_F_COUNT)
            / raw_total
        )
        required_f = removal * mean_f_per_si
        lower = max(removal - removal_uncertainty, 0.0) * max(
            mean_f_per_si - f_fraction_allowance, 0.0)
        upper = (removal + removal_uncertainty) * (
            mean_f_per_si + f_fraction_allowance)
        incident_f = float(REACTIVE_ION_FORMULAE["CF3+"]["F"])
        unresolved_f_balance = incident_f - required_f
        unresolved_f_balance_lower = incident_f - upper
        unresolved_f_balance_upper = incident_f - lower
        return ReactiveIonEventOutcome(
            species="CF3+",
            energy_eV=energy,
            removed_sio2_formula_per_ion=removal,
            removed_yield_digitization_uncertainty=removal_uncertainty,
            incident_atoms_per_ion={"C": 1.0, "F": incident_f},
            normalized_product_fraction=normalized,
            conditional_sifx_particles_per_ion=particles,
            branch_normalization_factor=1.0 / raw_total,
            required_f_atoms_per_ion=required_f,
            required_f_atoms_lower=lower,
            required_f_atoms_upper=upper,
            unresolved_f_balance_per_ion=unresolved_f_balance,
            unresolved_f_balance_lower=unresolved_f_balance_lower,
            unresolved_f_balance_upper=unresolved_f_balance_upper,
            condition_match_status="unresolved_figure10_incidence_angle",
            production_eligible=False,
            unresolved_inventories=(
                "Figure-10 ion incidence angle",
                "absolute SiFx product yield versus conditional fraction",
                "resident surface-fluorine supply or retention",
                "incident carbon: emitted CO versus retained film",
                "lattice oxygen coproduct identity",
                "prompt versus delayed product fraction",
                "formation-depth-dependent escape versus reincorporation",
            ),
        )
