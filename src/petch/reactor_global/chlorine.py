"""Lee--Lieberman chlorine particle-chemistry deck.

The particle reactions and rate fits are transcribed from Table 2 of
``lee-lieberman-1994-global``.  Wall loss/recombination and the separate
Tables 4--5 electron-energy channels are equipment/model closures and are not
silently folded into this volume deck.

Physical threshold energies for several composite fits have not yet been
recovered from the cited primary cross-section sources.  Their
``electron_energy_loss_eV`` fields are therefore ``None``.  The generic
network will solve particle sources and conservation, but will refuse an
electron-power calculation until that ledger is complete.
"""
from __future__ import annotations

from .argon import ELECTRON_MASS_AMU
from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    ElectronLogPolynomialRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)

CHLORINE_ATOM_MASS_AMU = 35.453
CHLORINE_MOLECULE_MASS_AMU = 2.0 * CHLORINE_ATOM_MASS_AMU
_LEE_TABLE2 = "lee-lieberman-1994-global Table 2"


def lee_lieberman_chlorine_species() -> tuple[Species, ...]:
    """Return the six tracked species in the printed chlorine mechanism."""
    measured = "measured"
    return (
        Species(
            name="e",
            mass_amu=ELECTRON_MASS_AMU,
            charge_number=-1,
            composition={},
            role="electron",
            source="CODATA electron mass",
            evidence_kind=measured,
        ),
        Species(
            name="Cl2",
            mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=0,
            composition={"Cl": 2},
            role="neutral",
            source="standard chlorine atomic weight",
            evidence_kind=measured,
        ),
        Species(
            name="Cl",
            mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=0,
            composition={"Cl": 1},
            role="neutral",
            source="standard chlorine atomic weight",
            evidence_kind=measured,
        ),
        Species(
            name="Cl2+",
            mass_amu=CHLORINE_MOLECULE_MASS_AMU,
            charge_number=1,
            composition={"Cl": 2},
            role="positive_ion",
            source=_LEE_TABLE2,
            evidence_kind="published_compilation",
        ),
        Species(
            name="Cl+",
            mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=1,
            composition={"Cl": 1},
            role="positive_ion",
            source=_LEE_TABLE2,
            evidence_kind="published_compilation",
        ),
        Species(
            name="Cl-",
            mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=-1,
            composition={"Cl": 1},
            role="negative_ion",
            source=_LEE_TABLE2,
            evidence_kind="published_compilation",
        ),
    )


def build_lee_lieberman_chlorine_particle_network() -> ReactionNetwork:
    """Build the closed Table-2 chlorine volume-particle mechanism.

    The source lumps the fine-structure neutral product of Cl-/Cl+
    neutralization.  Both neutral products are represented as tracked ground
    ``Cl`` because the source explicitly omits chlorine metastable balances.
    """
    evidence = "published_compilation"
    reactions = (
        Reaction(
            name="e_Cl2_nondissociative_ionization",
            reactants={"e": 1, "Cl2": 1},
            products={"Cl2+": 1, "e": 2},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                9.21e-8,
                activation_eV=12.9,
                source=f"{_LEE_TABLE2} first k1",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} first row",
        ),
        Reaction(
            name="e_Cl2_dissociative_ionization",
            reactants={"e": 1, "Cl2": 1},
            products={"Cl+": 1, "Cl": 1, "e": 2},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                3.88e-9,
                activation_eV=15.5,
                source=f"{_LEE_TABLE2} second k1",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} second row",
        ),
        Reaction(
            name="e_Cl2_ion_pair_production",
            reactants={"e": 1, "Cl2": 1},
            products={"Cl+": 1, "Cl-": 1, "e": 1},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                8.55e-10,
                activation_eV=12.65,
                source=f"{_LEE_TABLE2} third k1",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} third row",
        ),
        Reaction(
            name="e_Cl2_dissociation",
            reactants={"e": 1, "Cl2": 1},
            products={"Cl": 2, "e": 1},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                3.80e-8,
                activation_eV=3.824,
                source=f"{_LEE_TABLE2} k2",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} k2",
        ),
        Reaction(
            name="e_Cl2_dissociative_attachment",
            reactants={"e": 1, "Cl2": 1},
            products={"Cl": 1, "Cl-": 1},
            kinetic_orders={"e": 1, "Cl2": 1},
            rate_coefficient=(
                ElectronInverseTemperaturePolynomialRateCoefficient
                .from_cm3_per_s(
                    3.69e-10,
                    inverse_temperature_coefficients=(
                        -1.68,
                        1.457,
                        -0.44,
                        0.0572,
                        -0.0026,
                    ),
                    source=f"{_LEE_TABLE2} k3",
                    evidence_kind=evidence,
                )
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} k3",
        ),
        Reaction(
            name="e_Cl_ionization",
            reactants={"e": 1, "Cl": 1},
            products={"Cl+": 1, "e": 2},
            kinetic_orders={"e": 1, "Cl": 1},
            rate_coefficient=ElectronLogPolynomialRateCoefficient.from_cm3_per_s(
                (
                    1.419e-7,
                    -1.864e-8,
                    -5.439e-8,
                    3.306e-8,
                    -3.54e-9,
                    -2.915e-8,
                ),
                reference_temperature_eV=12.96,
                activation_eV=12.96,
                temperature_power=0.5,
                source=f"{_LEE_TABLE2} k4",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} k4",
        ),
        Reaction(
            name="Clminus_Cl2plus_neutralization",
            reactants={"Cl-": 1, "Cl2+": 1},
            products={"Cl": 1, "Cl2": 1},
            kinetic_orders={"Cl-": 1, "Cl2+": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                5.0e-8,
                source=f"{_LEE_TABLE2} k5",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=0.0,
            source=f"{_LEE_TABLE2} k5",
        ),
        Reaction(
            name="Clminus_Clplus_neutralization",
            reactants={"Cl-": 1, "Cl+": 1},
            products={"Cl": 2},
            kinetic_orders={"Cl-": 1, "Cl+": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                5.0e-8,
                source=f"{_LEE_TABLE2} k6",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=0.0,
            source=(
                f"{_LEE_TABLE2} k6; fine-structure neutral lumped to Cl"
            ),
        ),
        Reaction(
            name="e_Clminus_detachment",
            reactants={"e": 1, "Cl-": 1},
            products={"Cl": 1, "e": 2},
            kinetic_orders={"e": 1, "Cl-": 1},
            rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
                2.63e-8,
                activation_eV=5.37,
                source=f"{_LEE_TABLE2} k7",
                evidence_kind=evidence,
            ),
            electron_energy_loss_eV=None,
            source=f"{_LEE_TABLE2} k7",
        ),
    )
    return ReactionNetwork(
        species=lee_lieberman_chlorine_species(),
        reactions=reactions,
    )
