"""Lee--Graves--Lieberman SiClx gas-phase product chemistry.

The mechanism is a literal, unit-preserving transcription of Table IV in
UCB/ERL M95/9.  It expands a chlorine plasma with neutral and singly charged
``SiClx`` species for ``x = 0..4``.  The table's generalized mutual-
neutralization row is expanded into five separately conserved reactions.

The printed first-order Bohm wall-loss row is deliberately excluded: wall
loss is geometry-, state-, and species-dependent transport, not a volume
reaction coefficient.  Likewise, the fitted Arrhenius exponents are not
silently interpreted as physical energy loss per event.  The resulting deck
is therefore particle-complete but intentionally fails closed for an
electron-power calculation until an independently sourced energy ledger is
attached.
"""
from __future__ import annotations

from .argon import ELECTRON_MASS_AMU
from .chlorine import CHLORINE_ATOM_MASS_AMU
from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)


SILICON_ATOM_MASS_AMU = 28.0855
_LEE_PRODUCTS_TABLE4 = "lee-graves-lieberman-1995-etch-products Table IV"


def lee_graves_lieberman_etch_product_species() -> tuple[Species, ...]:
    """Return the electron plus ten SiClx neutrals/positive ions."""
    species = [Species(
        name="e",
        mass_amu=ELECTRON_MASS_AMU,
        charge_number=-1,
        composition={},
        role="electron",
        source="CODATA electron mass",
        evidence_kind="measured",
    )]
    for chlorine_count in range(5):
        neutral_name = (
            "Si" if chlorine_count == 0
            else "SiCl" if chlorine_count == 1
            else f"SiCl{chlorine_count}"
        )
        ion_name = f"{neutral_name}+"
        composition = {"Si": 1}
        if chlorine_count:
            composition["Cl"] = chlorine_count
        mass = (
            SILICON_ATOM_MASS_AMU
            + chlorine_count * CHLORINE_ATOM_MASS_AMU
        )
        species.extend((
            Species(
                name=neutral_name,
                mass_amu=mass,
                charge_number=0,
                composition=composition,
                role="neutral",
                source=_LEE_PRODUCTS_TABLE4,
                evidence_kind="published_compilation",
            ),
            Species(
                name=ion_name,
                mass_amu=mass,
                charge_number=1,
                composition=composition,
                role="positive_ion",
                source=_LEE_PRODUCTS_TABLE4,
                evidence_kind="published_compilation",
            ),
        ))
    return tuple(species)


def _arrhenius(
    prefactor_cm3_s: float,
    activation_eV: float,
):
    return ElectronArrheniusRateCoefficient.from_cm3_per_s(
        prefactor_cm3_s,
        activation_eV=activation_eV,
        source=_LEE_PRODUCTS_TABLE4,
        evidence_kind="regressed",
    )


def _inverse_polynomial(
    prefactor_cm3_s: float,
    coefficients: tuple[float, ...],
):
    return ElectronInverseTemperaturePolynomialRateCoefficient.from_cm3_per_s(
        prefactor_cm3_s,
        inverse_temperature_coefficients=coefficients,
        source=_LEE_PRODUCTS_TABLE4,
        evidence_kind="regressed",
    )


def lee_graves_lieberman_etch_product_reactions() -> tuple[Reaction, ...]:
    """Return Table-IV gas reactions with exact printed rate fits."""
    rows = (
        ("k1_sicl4_ionization", {"e": 1, "SiCl4": 1},
         {"SiCl4+": 1, "e": 2}, _arrhenius(7.03e-8, 12.44)),
        ("k2_sicl4_dissociation", {"e": 1, "SiCl4": 1},
         {"SiCl3": 1, "Cl": 1, "e": 1}, _arrhenius(7.27e-9, 4.73)),
        ("k3_sicl4_dissociative_ionization", {"e": 1, "SiCl4": 1},
         {"SiCl3+": 1, "Cl": 1, "e": 2}, _arrhenius(2.00e-7, 12.44)),
        ("k4_sicl3_ionization", {"e": 1, "SiCl3": 1},
         {"SiCl3+": 1, "e": 2}, _arrhenius(1.68e-8, 7.65)),
        ("k5_sicl3_dissociation", {"e": 1, "SiCl3": 1},
         {"SiCl2": 1, "Cl": 1, "e": 1}, _arrhenius(7.27e-9, 2.91)),
        ("k6_sicl3_dissociative_ionization_sicl2", {"e": 1, "SiCl3": 1},
         {"SiCl2+": 1, "Cl": 1, "e": 2},
         _inverse_polynomial(4.90e-8, (-13.9, 6.89, -1.45))),
        ("k7_sicl3_dissociative_ionization_sicl", {"e": 1, "SiCl3": 1},
         {"SiCl+": 1, "Cl": 2, "e": 2},
         _inverse_polynomial(2.41e-8, (-14.75, 2.504))),
        ("k8_sicl2_dissociation", {"e": 1, "SiCl2": 1},
         {"SiCl": 1, "Cl": 1, "e": 1}, _arrhenius(7.27e-9, 4.99)),
        ("k9_sicl2_ionization", {"e": 1, "SiCl2": 1},
         {"SiCl2+": 1, "e": 2}, _arrhenius(2.98e-8, 9.81)),
        ("k10_sicl2_dissociative_ionization", {"e": 1, "SiCl2": 1},
         {"SiCl+": 1, "Cl": 1, "e": 2}, _arrhenius(8.93e-8, 9.81)),
        ("k11_sicl_ionization", {"e": 1, "SiCl": 1},
         {"SiCl+": 1, "e": 2}, _arrhenius(7.54e-8, 6.79)),
        ("k12_sicl_dissociation", {"e": 1, "SiCl": 1},
         {"Si": 1, "Cl": 1, "e": 1}, _arrhenius(7.27e-9, 3.95)),
        ("k13_sicl_dissociative_ionization", {"e": 1, "SiCl": 1},
         {"Si+": 1, "Cl": 1, "e": 2}, _arrhenius(8.85e-8, 12.1)),
        ("k14_si_ionization", {"e": 1, "Si": 1},
         {"Si+": 1, "e": 2}, _arrhenius(7.85e-8, 7.41)),
    )
    reactions = [Reaction(
        name=name,
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=coefficient,
        electron_energy_loss_eV=None,
        source=_LEE_PRODUCTS_TABLE4,
    ) for name, reactants, products, coefficient in rows]
    neutralization = ConstantRateCoefficient.from_cm3_per_s(
        5.0e-8,
        source=_LEE_PRODUCTS_TABLE4,
        evidence_kind="estimated",
    )
    for chlorine_count in range(5):
        neutral = (
            "Si" if chlorine_count == 0
            else "SiCl" if chlorine_count == 1
            else f"SiCl{chlorine_count}"
        )
        ion = f"{neutral}+"
        reactions.append(Reaction(
            name=f"k15_{ion.lower()}_mutual_neutralization",
            reactants={ion: 1, "Cl-": 1},
            products={neutral: 1, "Cl": 1},
            kinetic_orders={ion: 1, "Cl-": 1},
            rate_coefficient=neutralization,
            electron_energy_loss_eV=0.0,
            source=_LEE_PRODUCTS_TABLE4,
        ))
    return tuple(reactions)


def build_lee_graves_lieberman_etch_product_network(
) -> ReactionNetwork:
    """Build the conserved particle mechanism from Table IV."""
    chlorine_species = (
        Species(
            name="Cl",
            mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=0,
            composition={"Cl": 1},
            role="neutral",
            source=_LEE_PRODUCTS_TABLE4,
            evidence_kind="published_compilation",
        ),
        Species(
            name="Cl-",
            mass_amu=CHLORINE_ATOM_MASS_AMU,
            charge_number=-1,
            composition={"Cl": 1},
            role="negative_ion",
            source=_LEE_PRODUCTS_TABLE4,
            evidence_kind="published_compilation",
        ),
    )
    return ReactionNetwork(
        species=lee_graves_lieberman_etch_product_species() + chlorine_species,
        reactions=lee_graves_lieberman_etch_product_reactions(),
    )
