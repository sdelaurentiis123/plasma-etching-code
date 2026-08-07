"""Conservative stratified fluorocarbon/silicon reaction state.

This module is the topology core for a replacement of ``mixed_layer.py``.  It
does not contain an empirical etch-yield law.  Instead it separates the
physical reservoirs demanded by the Humbird--Graves molecular-dynamics
mechanism:

1. a fluorocarbon film with explicit C, F, C--F bonds, and C--C crosslinks;
2. a finite-capacity Si--C transport layer containing carbon and silicon;
3. a subsurface Si--F reaction front; and
4. Si--F-bearing silicon in transit through the Si--C layer.

An :class:`StratifiedSiEvents` object is an integrated, mechanistically named
event budget.  :func:`advance_stratified_si` applies it without clipping and
rejects impossible budgets.  Element and bond ledgers close algebraically.
Silicon removed in a step is limited to the transport inventory present at
the *start* of that step, so newly promoted silicon cannot be volatilized with
zero residence time.

Ion transmission through the live film and transport layer is calculated by
successive inversion of the ZBL/Lindhard CSDA path.  This replaces the legacy
``E exp(-d/lambda)`` attenuation with finite range and the correct
``depth/cos(theta)`` material path.  The stopping calculation is still a
straight-path BCA-level approximation; straggling and reactive-potential
uncertainty are explicit omissions.

The topology and ledgers are chemistry mechanics, not a claim of quantum
accuracy.  A predictive kinetic provider must supply the event budget from
direct beam/MD evidence and is intentionally outside this module.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .ion_energy_deposition import (
    FLUOROCARBON_FILM,
    Target,
    TargetComponent,
    residual_energy_after_layer_eV,
)


# Beta-SiC at 3.21 g/cm3 gives about 4.83e28 formula units/m3.
SILICON_CARBIDE = Target(
    components=(
        TargetComponent(14, 28.085, 1.0),
        TargetComponent(6, 12.011, 1.0),
    ),
    atom_density_m3=9.66e28,
)


def _immutable_nonnegative(value, name):
    array = np.asarray(value, dtype=float).copy()
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    array.setflags(write=False)
    return array


def _broadcast_named(values):
    names = tuple(values)
    try:
        arrays = np.broadcast_arrays(*[
            np.asarray(values[name], dtype=float) for name in names])
    except ValueError as error:
        raise ValueError("stratified surface fields do not broadcast") from error
    return {
        name: np.array(value, copy=True)
        for name, value in zip(names, arrays)
    }


def _tolerance(*values):
    scale = np.ones(())
    for value in values:
        scale = np.maximum(scale, np.abs(np.asarray(value, dtype=float)))
    return 128.0 * np.finfo(float).eps * scale


@dataclass(frozen=True)
class StratifiedSiParameters:
    """Physical capacities and projectile identity.

    Layer capacities derive from number density times physical depth.  The
    default 3 nm Si--C depth is the 200 eV endpoint reported by Humbird et al.;
    it is a source constraint, not a monolayer saturation knob.  Kinetic
    probabilities do not belong in this parameter object.
    """

    film_atom_density_m3: float = FLUOROCARBON_FILM.atom_density_m3
    transport_atom_density_m3: float = SILICON_CARBIDE.atom_density_m3
    transport_layer_depth_nm: float = 3.0
    silicon_atom_density_m3: float = 5.0e28
    reaction_front_depth_nm: float = 1.5
    maximum_f_bonds_per_si: int = 4
    maximum_valence_per_film_c: int = 4
    projectile_atomic_number: int = 18
    projectile_mass_amu: float = 39.948
    evidence: Mapping[str, str] | None = None

    def __post_init__(self):
        positive = (
            self.film_atom_density_m3,
            self.transport_atom_density_m3,
            self.transport_layer_depth_nm,
            self.silicon_atom_density_m3,
            self.reaction_front_depth_nm,
            self.projectile_mass_amu,
        )
        if (any(not np.isfinite(value) or value <= 0.0 for value in positive)
                or int(self.maximum_f_bonds_per_si)
                != self.maximum_f_bonds_per_si
                or self.maximum_f_bonds_per_si <= 0
                or int(self.maximum_valence_per_film_c)
                != self.maximum_valence_per_film_c
                or self.maximum_valence_per_film_c <= 0
                or int(self.projectile_atomic_number)
                != self.projectile_atomic_number
                or self.projectile_atomic_number <= 0):
            raise ValueError("invalid stratified silicon parameters")
        evidence = dict(self.evidence or {
            "film_atom_density_m3": (
                "fluorocarbon stopping target used by petch; "
                "composition-density approximation"),
            "transport_atom_density_m3": (
                "beta-SiC density 3.21 g/cm3 and molar mass"),
            "transport_layer_depth_nm": (
                "Humbird et al. APL 84, 1073 (2004): about 30 A "
                "Si-C-F region at 200 eV"),
            "silicon_atom_density_m3": (
                "crystalline-Si density and molar mass"),
            "reaction_front_depth_nm": (
                "Humbird et al. APL 84, 1073 (2004): about 15 A "
                "low-energy mixed-region endpoint"),
            "maximum_f_bonds_per_si": "silicon valence",
            "maximum_valence_per_film_c": "carbon valence",
            "projectile": "mass-selected Ar+ boundary",
        })
        required = {
            "film_atom_density_m3",
            "transport_atom_density_m3",
            "transport_layer_depth_nm",
            "silicon_atom_density_m3",
            "reaction_front_depth_nm",
            "maximum_f_bonds_per_si",
            "maximum_valence_per_film_c",
            "projectile",
        }
        if set(evidence) != required or any(
                not isinstance(value, str) or not value
                for value in evidence.values()):
            raise ValueError("incomplete stratified-parameter provenance")
        object.__setattr__(self, "maximum_f_bonds_per_si",
                           int(self.maximum_f_bonds_per_si))
        object.__setattr__(self, "maximum_valence_per_film_c",
                           int(self.maximum_valence_per_film_c))
        object.__setattr__(self, "projectile_atomic_number",
                           int(self.projectile_atomic_number))
        object.__setattr__(self, "evidence", MappingProxyType(evidence))

    @property
    def transport_capacity_atoms_m2(self):
        return (
            self.transport_atom_density_m3
            * self.transport_layer_depth_nm
            * 1.0e-9
        )

    @property
    def reaction_front_si_sites_m2(self):
        return (
            self.silicon_atom_density_m3
            * self.reaction_front_depth_nm
            * 1.0e-9
        )

    @property
    def reaction_front_f_bond_capacity_m2(self):
        return (
            self.maximum_f_bonds_per_si
            * self.reaction_front_si_sites_m2
        )


@dataclass(frozen=True)
class StratifiedSiState:
    """Element and bond inventories per square metre."""

    film_c_atoms_m2: np.ndarray | float = 0.0
    film_f_atoms_m2: np.ndarray | float = 0.0
    film_cf_bonds_m2: np.ndarray | float = 0.0
    film_cc_crosslinks_m2: np.ndarray | float = 0.0
    transport_c_atoms_m2: np.ndarray | float = 0.0
    transport_si_atoms_m2: np.ndarray | float = 0.0
    transport_si_f_bonds_m2: np.ndarray | float = 0.0
    reaction_front_f_bonds_m2: np.ndarray | float = 0.0
    cumulative_removed_si_atoms_m2: np.ndarray | float = 0.0
    cumulative_drawn_si_atoms_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        values = _broadcast_named({
            item.name: getattr(self, item.name) for item in fields(self)
        })
        for name, value in values.items():
            value = _immutable_nonnegative(value, name)
            object.__setattr__(self, name, value)

    @classmethod
    def bare(cls, shape=()):
        zero = np.zeros(shape, dtype=float)
        return cls(*(zero for _ in fields(cls)))

    @property
    def shape(self):
        return self.film_c_atoms_m2.shape

    def validate(self, parameters: StratifiedSiParameters):
        if not isinstance(parameters, StratifiedSiParameters):
            raise TypeError("parameters must be StratifiedSiParameters")
        tolerance = _tolerance(
            self.film_c_atoms_m2,
            self.film_f_atoms_m2,
            self.transport_c_atoms_m2,
            self.transport_si_atoms_m2,
        )
        if np.any(
            self.film_cf_bonds_m2
            > np.minimum(
                self.film_f_atoms_m2,
                parameters.maximum_valence_per_film_c
                * self.film_c_atoms_m2,
            )
            + tolerance
        ):
            raise ValueError("C-F bond inventory exceeds film atoms or C valence")
        if np.any(
            self.film_cf_bonds_m2
            + 2.0 * self.film_cc_crosslinks_m2
            > parameters.maximum_valence_per_film_c
            * self.film_c_atoms_m2
            + tolerance
        ):
            raise ValueError("film bond inventory exceeds carbon valence")
        if np.any(
            self.transport_c_atoms_m2 + self.transport_si_atoms_m2
            > parameters.transport_capacity_atoms_m2 + tolerance
        ):
            raise ValueError("Si-C transport-layer capacity exceeded")
        if np.any(
            self.transport_si_f_bonds_m2
            > parameters.maximum_f_bonds_per_si
            * self.transport_si_atoms_m2
            + tolerance
        ):
            raise ValueError("transport Si-F bonds exceed silicon valence")
        if np.any(
            self.reaction_front_f_bonds_m2
            > parameters.reaction_front_f_bond_capacity_m2 + tolerance
        ):
            raise ValueError("reaction-front Si-F bond capacity exceeded")
        if np.any(
            np.abs(
                self.cumulative_drawn_si_atoms_m2
                - self.transport_si_atoms_m2
                - self.cumulative_removed_si_atoms_m2
            ) > _tolerance(
                self.cumulative_drawn_si_atoms_m2,
                self.transport_si_atoms_m2,
                self.cumulative_removed_si_atoms_m2,
            )
        ):
            raise ValueError("cumulative silicon inventory does not close")
        return self

    def film_thickness_nm(self, parameters: StratifiedSiParameters):
        return (
            (self.film_c_atoms_m2 + self.film_f_atoms_m2)
            / parameters.film_atom_density_m3
            * 1.0e9
        )

    def transport_thickness_nm(self, parameters: StratifiedSiParameters):
        return (
            (self.transport_c_atoms_m2 + self.transport_si_atoms_m2)
            / parameters.transport_atom_density_m3
            * 1.0e9
        )

    def conservative_surface_fields(self):
        return {
            item.name: getattr(self, item.name) for item in fields(self)
        }

    def conservative_surface_upper_bounds(self):
        return {item.name: None for item in fields(self)}

    def surface_field_remap_modes(self):
        return {item.name: "conservative" for item in fields(self)}

    def with_conservative_surface_fields(self, supplied):
        supplied = dict(supplied)
        expected = {item.name for item in fields(self)}
        if set(supplied) != expected:
            raise ValueError("stratified Si remap fields do not match")
        return type(self)(**supplied)


@dataclass(frozen=True)
class StratifiedSiEvents:
    """Integrated event counts per square metre for one accepted step.

    Every value is nonnegative.  Internal transfers appear explicitly; they
    cancel in the element ledger but change bond and reservoir topology.
    """

    deposited_film_c_atoms_m2: np.ndarray | float = 0.0
    deposited_film_f_atoms_m2: np.ndarray | float = 0.0
    formed_film_cf_bonds_m2: np.ndarray | float = 0.0
    formed_film_cc_crosslinks_m2: np.ndarray | float = 0.0
    adsorbed_atomic_f_to_front_m2: np.ndarray | float = 0.0

    broken_film_cf_bonds_m2: np.ndarray | float = 0.0
    broken_film_cc_crosslinks_m2: np.ndarray | float = 0.0
    transferred_film_f_to_front_m2: np.ndarray | float = 0.0
    transferred_film_c_to_transport_m2: np.ndarray | float = 0.0

    promoted_bulk_si_to_transport_m2: np.ndarray | float = 0.0
    promoted_front_f_bonds_to_transport_m2: np.ndarray | float = 0.0
    recycled_transport_f_bonds_to_front_m2: np.ndarray | float = 0.0

    removed_film_c_atoms_m2: np.ndarray | float = 0.0
    removed_film_f_atoms_m2: np.ndarray | float = 0.0
    exported_film_cf_bonds_m2: np.ndarray | float = 0.0
    removed_transport_c_atoms_m2: np.ndarray | float = 0.0
    removed_transport_si_atoms_m2: np.ndarray | float = 0.0
    exported_transport_si_f_bonds_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        values = _broadcast_named({
            item.name: getattr(self, item.name) for item in fields(self)
        })
        for name, value in values.items():
            value = _immutable_nonnegative(value, name)
            object.__setattr__(self, name, value)

    @property
    def shape(self):
        return self.deposited_film_c_atoms_m2.shape

    @classmethod
    def zero(cls, shape=()):
        zero = np.zeros(shape, dtype=float)
        return cls(*(zero for _ in fields(cls)))


@dataclass(frozen=True)
class StratifiedSiStepResult:
    state: StratifiedSiState
    element_ledger_residual_atoms_m2: Mapping[str, np.ndarray]
    bond_ledger_residual_bonds_m2: Mapping[str, np.ndarray]
    removed_si_atoms_m2: np.ndarray
    removed_carbon_atoms_m2: np.ndarray
    removed_f_atoms_m2: np.ndarray
    substrate_si_drawn_m2: np.ndarray

    def __post_init__(self):
        element = {
            name: _immutable_nonnegative(np.abs(value), name)
            for name, value in dict(
                self.element_ledger_residual_atoms_m2).items()
        }
        bond = {
            name: _immutable_nonnegative(np.abs(value), name)
            for name, value in dict(
                self.bond_ledger_residual_bonds_m2).items()
        }
        if set(element) != {"C", "F", "Si"}:
            raise ValueError("element ledger must contain C, F, and Si")
        if set(bond) != {"C-F", "Si-F", "C-C_crosslink"}:
            raise ValueError("bond ledger is incomplete")
        object.__setattr__(
            self, "element_ledger_residual_atoms_m2",
            MappingProxyType(element))
        object.__setattr__(
            self, "bond_ledger_residual_bonds_m2",
            MappingProxyType(bond))
        for name in (
            "removed_si_atoms_m2",
            "removed_carbon_atoms_m2",
            "removed_f_atoms_m2",
            "substrate_si_drawn_m2",
        ):
            object.__setattr__(
                self, name, _immutable_nonnegative(getattr(self, name), name))

    @property
    def maximum_absolute_ledger_residual(self):
        values = tuple(self.element_ledger_residual_atoms_m2.values())
        values += tuple(self.bond_ledger_residual_bonds_m2.values())
        return max((float(np.max(value)) for value in values), default=0.0)


@dataclass(frozen=True)
class StratifiedIonEnergies:
    incident_energy_eV: np.ndarray
    after_film_energy_eV: np.ndarray
    at_reaction_front_energy_eV: np.ndarray
    film_thickness_nm: np.ndarray
    transport_thickness_nm: np.ndarray

    def __post_init__(self):
        for item in fields(self):
            object.__setattr__(
                self, item.name,
                _immutable_nonnegative(getattr(self, item.name), item.name))


def stratified_ion_energies(
        state: StratifiedSiState, energy_eV, cosine_incidence,
        parameters: StratifiedSiParameters | None = None):
    """Transmit ions through the live film and Si--C layer by CSDA inversion."""
    parameters = parameters or StratifiedSiParameters()
    state.validate(parameters)
    energy = np.asarray(energy_eV, dtype=float)
    cosine = np.asarray(cosine_incidence, dtype=float)
    film_depth = np.asarray(state.film_thickness_nm(parameters), dtype=float)
    transport_depth = np.asarray(
        state.transport_thickness_nm(parameters), dtype=float)
    try:
        energy, cosine, film_depth, transport_depth = np.broadcast_arrays(
            energy, cosine, film_depth, transport_depth)
    except ValueError as error:
        raise ValueError("ion fields and stratified state do not broadcast") from error
    after_film = residual_energy_after_layer_eV(
        energy,
        cosine,
        film_depth,
        parameters.projectile_atomic_number,
        parameters.projectile_mass_amu,
        FLUOROCARBON_FILM,
    )
    at_front = residual_energy_after_layer_eV(
        after_film,
        cosine,
        transport_depth,
        parameters.projectile_atomic_number,
        parameters.projectile_mass_amu,
        SILICON_CARBIDE,
    )
    return StratifiedIonEnergies(
        energy,
        after_film,
        at_front,
        film_depth,
        transport_depth,
    )


def _require_not_greater(name, value, available):
    if np.any(value > available + _tolerance(value, available)):
        raise ValueError(f"{name} exceeds the start-of-step inventory")


def advance_stratified_si(
        state: StratifiedSiState,
        events: StratifiedSiEvents,
        parameters: StratifiedSiParameters | None = None,
):
    """Apply one exact event budget and return closed element/bond ledgers."""
    parameters = parameters or StratifiedSiParameters()
    if not isinstance(state, StratifiedSiState):
        raise TypeError("state must be StratifiedSiState")
    if not isinstance(events, StratifiedSiEvents):
        raise TypeError("events must be StratifiedSiEvents")
    state.validate(parameters)
    try:
        shape = np.broadcast_shapes(state.shape, events.shape)
    except ValueError as error:
        raise ValueError("state and event fields do not broadcast") from error

    def s(name):
        return np.broadcast_to(np.asarray(getattr(state, name)), shape)

    def e(name):
        return np.broadcast_to(np.asarray(getattr(events, name)), shape)

    # Removal is deliberately checked against the initial state.  Internal
    # promotion in this same budget cannot be used to satisfy removal.
    _require_not_greater(
        "removed transport silicon",
        e("removed_transport_si_atoms_m2"),
        s("transport_si_atoms_m2"),
    )
    _require_not_greater(
        "exported transport Si-F bonds",
        e("exported_transport_si_f_bonds_m2"),
        s("transport_si_f_bonds_m2"),
    )
    if np.any(
        e("exported_transport_si_f_bonds_m2")
        > parameters.maximum_f_bonds_per_si
        * e("removed_transport_si_atoms_m2")
        + _tolerance(
            e("exported_transport_si_f_bonds_m2"),
            e("removed_transport_si_atoms_m2"),
        )
    ):
        raise ValueError("exported Si-F bonds exceed removed-Si valence")
    _require_not_greater(
        "removed transport carbon",
        e("removed_transport_c_atoms_m2"),
        s("transport_c_atoms_m2"),
    )
    _require_not_greater(
        "removed film carbon",
        e("removed_film_c_atoms_m2")
        + e("transferred_film_c_to_transport_m2"),
        s("film_c_atoms_m2"),
    )
    _require_not_greater(
        "removed film fluorine",
        e("removed_film_f_atoms_m2")
        + e("transferred_film_f_to_front_m2"),
        s("film_f_atoms_m2"),
    )
    _require_not_greater(
        "destroyed or exported C-F bonds",
        e("broken_film_cf_bonds_m2") + e("exported_film_cf_bonds_m2"),
        s("film_cf_bonds_m2"),
    )
    _require_not_greater(
        "destroyed C-C crosslinks",
        e("broken_film_cc_crosslinks_m2"),
        s("film_cc_crosslinks_m2"),
    )
    if np.any(
        e("transferred_film_f_to_front_m2")
        > e("broken_film_cf_bonds_m2")
        + _tolerance(
            e("transferred_film_f_to_front_m2"),
            e("broken_film_cf_bonds_m2"),
        )
    ):
        raise ValueError("F transfer to the Si-F front requires C-F scission")
    if np.any(
        e("formed_film_cf_bonds_m2")
        > e("deposited_film_f_atoms_m2")
        + _tolerance(
            e("formed_film_cf_bonds_m2"),
            e("deposited_film_f_atoms_m2"),
        )
    ):
        raise ValueError("new C-F bonds exceed deposited fluorine")
    if np.any(
        e("formed_film_cf_bonds_m2")
        > parameters.maximum_valence_per_film_c
        * e("deposited_film_c_atoms_m2")
        + _tolerance(
            e("formed_film_cf_bonds_m2"),
            e("deposited_film_c_atoms_m2"),
        )
    ):
        raise ValueError("new C-F bonds exceed deposited-carbon valence")
    if np.any(
        e("formed_film_cf_bonds_m2")
        + 2.0 * e("formed_film_cc_crosslinks_m2")
        > parameters.maximum_valence_per_film_c
        * e("deposited_film_c_atoms_m2")
        + _tolerance(
            e("formed_film_cf_bonds_m2"),
            e("formed_film_cc_crosslinks_m2"),
            e("deposited_film_c_atoms_m2"),
        )
    ):
        raise ValueError("new film bonds exceed deposited-carbon valence")
    if np.any(
        e("promoted_front_f_bonds_to_transport_m2")
        > parameters.maximum_f_bonds_per_si
        * e("promoted_bulk_si_to_transport_m2")
        + _tolerance(
            e("promoted_front_f_bonds_to_transport_m2"),
            e("promoted_bulk_si_to_transport_m2"),
        )
    ):
        raise ValueError("promoted Si-F bonds exceed promoted-Si valence")

    # Existing transport Si-F bonds after product export are then eligible for
    # recycling.  This ordering prevents a single bond from both leaving and
    # returning to the reaction front.
    _require_not_greater(
        "recycled transport Si-F bonds",
        e("recycled_transport_f_bonds_to_front_m2"),
        s("transport_si_f_bonds_m2")
        - e("exported_transport_si_f_bonds_m2"),
    )

    film_c = (
        s("film_c_atoms_m2")
        + e("deposited_film_c_atoms_m2")
        - e("removed_film_c_atoms_m2")
        - e("transferred_film_c_to_transport_m2")
    )
    film_f = (
        s("film_f_atoms_m2")
        + e("deposited_film_f_atoms_m2")
        - e("removed_film_f_atoms_m2")
        - e("transferred_film_f_to_front_m2")
    )
    film_cf = (
        s("film_cf_bonds_m2")
        + e("formed_film_cf_bonds_m2")
        - e("broken_film_cf_bonds_m2")
        - e("exported_film_cf_bonds_m2")
    )
    film_cc = (
        s("film_cc_crosslinks_m2")
        + e("formed_film_cc_crosslinks_m2")
        - e("broken_film_cc_crosslinks_m2")
    )
    transport_c = (
        s("transport_c_atoms_m2")
        + e("transferred_film_c_to_transport_m2")
        - e("removed_transport_c_atoms_m2")
    )
    transport_si = (
        s("transport_si_atoms_m2")
        + e("promoted_bulk_si_to_transport_m2")
        - e("removed_transport_si_atoms_m2")
    )
    transport_sif = (
        s("transport_si_f_bonds_m2")
        + e("promoted_front_f_bonds_to_transport_m2")
        - e("recycled_transport_f_bonds_to_front_m2")
        - e("exported_transport_si_f_bonds_m2")
    )
    front_f = (
        s("reaction_front_f_bonds_m2")
        + e("adsorbed_atomic_f_to_front_m2")
        + e("transferred_film_f_to_front_m2")
        + e("recycled_transport_f_bonds_to_front_m2")
        - e("promoted_front_f_bonds_to_transport_m2")
    )
    removed_si = e("removed_transport_si_atoms_m2")
    drawn_si = e("promoted_bulk_si_to_transport_m2")
    updated = StratifiedSiState(
        film_c,
        film_f,
        film_cf,
        film_cc,
        transport_c,
        transport_si,
        transport_sif,
        front_f,
        s("cumulative_removed_si_atoms_m2") + removed_si,
        s("cumulative_drawn_si_atoms_m2") + drawn_si,
    ).validate(parameters)

    initial_c = s("film_c_atoms_m2") + s("transport_c_atoms_m2")
    final_c = updated.film_c_atoms_m2 + updated.transport_c_atoms_m2
    input_c = e("deposited_film_c_atoms_m2")
    output_c = (
        e("removed_film_c_atoms_m2")
        + e("removed_transport_c_atoms_m2")
    )
    initial_f = (
        s("film_f_atoms_m2")
        + s("transport_si_f_bonds_m2")
        + s("reaction_front_f_bonds_m2")
    )
    final_f = (
        updated.film_f_atoms_m2
        + updated.transport_si_f_bonds_m2
        + updated.reaction_front_f_bonds_m2
    )
    input_f = (
        e("deposited_film_f_atoms_m2")
        + e("adsorbed_atomic_f_to_front_m2")
    )
    output_f = (
        e("removed_film_f_atoms_m2")
        + e("exported_transport_si_f_bonds_m2")
    )
    initial_si = s("transport_si_atoms_m2")
    final_si = updated.transport_si_atoms_m2
    element_residual = {
        "C": initial_c + input_c - output_c - final_c,
        "F": initial_f + input_f - output_f - final_f,
        "Si": initial_si + drawn_si - removed_si - final_si,
    }

    initial_sif = (
        s("transport_si_f_bonds_m2")
        + s("reaction_front_f_bonds_m2")
    )
    final_sif = (
        updated.transport_si_f_bonds_m2
        + updated.reaction_front_f_bonds_m2
    )
    bond_residual = {
        "C-F": (
            s("film_cf_bonds_m2")
            + e("formed_film_cf_bonds_m2")
            - e("broken_film_cf_bonds_m2")
            - e("exported_film_cf_bonds_m2")
            - updated.film_cf_bonds_m2
        ),
        "Si-F": (
            initial_sif
            + e("adsorbed_atomic_f_to_front_m2")
            + e("transferred_film_f_to_front_m2")
            - e("exported_transport_si_f_bonds_m2")
            - final_sif
        ),
        "C-C_crosslink": (
            s("film_cc_crosslinks_m2")
            + e("formed_film_cc_crosslinks_m2")
            - e("broken_film_cc_crosslinks_m2")
            - updated.film_cc_crosslinks_m2
        ),
    }
    maximum = max(
        float(np.max(np.abs(value)))
        for value in (*element_residual.values(), *bond_residual.values())
    )
    throughput = max(
        1.0,
        *(float(np.max(np.abs(value))) for value in (
            input_c, input_f, drawn_si, output_c, output_f, removed_si,
        )),
    )
    if maximum > 256.0 * np.finfo(float).eps * throughput:
        raise RuntimeError("stratified surface ledger failed to close")
    return StratifiedSiStepResult(
        state=updated,
        element_ledger_residual_atoms_m2=element_residual,
        bond_ledger_residual_bonds_m2=bond_residual,
        removed_si_atoms_m2=removed_si,
        removed_carbon_atoms_m2=output_c,
        removed_f_atoms_m2=output_f,
        substrate_si_drawn_m2=drawn_si,
    )
