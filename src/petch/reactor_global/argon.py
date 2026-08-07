"""Lee--Lieberman argon volume-chemistry deck.

The rate fits are imported exactly from Table 3 of
``lee-lieberman-1994-global``.  Physical event energies are kept separate from
fit exponents and come from ``nist-asd-argon``.
"""
from __future__ import annotations

from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)

ARGON_MASS_AMU = 39.948
ELECTRON_MASS_AMU = 5.48579909065e-4
ARGON_IONIZATION_ENERGY_EV = 15.7596119
ARGON_4S_METASTABLE_ENERGY_EV = 11.54835442
ARGON_METASTABLE_IONIZATION_ENERGY_EV = (
    ARGON_IONIZATION_ENERGY_EV - ARGON_4S_METASTABLE_ENERGY_EV)

_LEE_TABLE3 = "lee-lieberman-1994-global Table 3"
_NIST_ASD = "nist-asd-argon"


def lee_lieberman_argon_species() -> tuple[Species, ...]:
    """Return the four species required by the Table-3 volume reactions."""
    return (
        Species(
            name="e",
            mass_amu=ELECTRON_MASS_AMU,
            charge_number=-1,
            composition={},
            role="electron",
            source=_NIST_ASD,
            evidence_kind="measured",
        ),
        Species(
            name="Ar",
            mass_amu=ARGON_MASS_AMU,
            charge_number=0,
            composition={"Ar": 1},
            role="neutral",
            source=_NIST_ASD,
            evidence_kind="measured",
        ),
        Species(
            name="Ar*",
            mass_amu=ARGON_MASS_AMU,
            charge_number=0,
            composition={"Ar": 1},
            role="excited_neutral",
            source=_NIST_ASD,
            evidence_kind="measured",
        ),
        Species(
            name="Ar+",
            mass_amu=ARGON_MASS_AMU,
            charge_number=1,
            composition={"Ar": 1},
            role="positive_ion",
            source=_NIST_ASD,
            evidence_kind="measured",
        ),
    )


def build_lee_lieberman_argon_volume_network() -> ReactionNetwork:
    """Build the five closed Table-3 volume reactions without wall closures.

    Table 3 omits the emitted electron from metastable pooling. It is explicit
    here because associative ionization is otherwise charge-nonconserving.
    """
    evidence = "published_compilation"
    reactions = (
        Reaction(
            name="e_Ar_ground_ionization",
            reactants={"e": 1, "Ar": 1},
            products={"Ar+": 1, "e": 2},
            kinetic_orders={"e": 1, "Ar": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                1.23e-7,
                activation_eV=18.68,
                source=f"{_LEE_TABLE3} k1",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=ARGON_IONIZATION_ENERGY_EV,
            source=f"{_LEE_TABLE3} reaction 1; energy {_NIST_ASD}",
        ),
        Reaction(
            name="e_Ar_metastable_excitation",
            reactants={"e": 1, "Ar": 1},
            products={"Ar*": 1, "e": 1},
            kinetic_orders={"e": 1, "Ar": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                3.71e-8,
                activation_eV=15.06,
                source=f"{_LEE_TABLE3} k2",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=ARGON_4S_METASTABLE_ENERGY_EV,
            source=f"{_LEE_TABLE3} reaction 2; energy {_NIST_ASD}",
        ),
        Reaction(
            name="e_Ar_metastable_step_ionization",
            reactants={"e": 1, "Ar*": 1},
            products={"Ar+": 1, "e": 2},
            kinetic_orders={"e": 1, "Ar*": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                2.05e-7,
                activation_eV=4.95,
                source=f"{_LEE_TABLE3} k3",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=ARGON_METASTABLE_IONIZATION_ENERGY_EV,
            source=f"{_LEE_TABLE3} reaction 3; energy difference {_NIST_ASD}",
        ),
        Reaction(
            name="e_Ar_metastable_superelastic_quench",
            reactants={"e": 1, "Ar*": 1},
            products={"Ar": 1, "e": 1},
            kinetic_orders={"e": 1, "Ar*": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                2.0e-7,
                source=f"{_LEE_TABLE3} k4",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=-ARGON_4S_METASTABLE_ENERGY_EV,
            source=f"{_LEE_TABLE3} reaction 4; energy {_NIST_ASD}",
        ),
        Reaction(
            name="Ar_metastable_pooling_associative_ionization",
            reactants={"Ar*": 2},
            products={"Ar": 1, "Ar+": 1, "e": 1},
            kinetic_orders={"Ar*": 2},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                6.2e-10,
                source=f"{_LEE_TABLE3} k5",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=0.0,
            source=(
                f"{_LEE_TABLE3} reaction 5 with emitted electron restored "
                "for closed charge conservation"
            ),
        ),
    )
    return ReactionNetwork(
        species=lee_lieberman_argon_species(),
        reactions=reactions,
    )
