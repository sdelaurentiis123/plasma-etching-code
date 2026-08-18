"""Conserved CHF3/O2 daughter chemistry from Lim et al. (2014).

Table I of Lim et al., DOI 10.1166/jnn.2014.10171, prints a compact
daughter-electron and neutral-neutral mechanism used with measured plasma
diagnostics.  Parent-feed electron rows and rows already represented by the
measured Zhu collision decks are deliberately omitted here.

The printed neutral rates are constants evaluated for an assumed 700 K gas.
They are retained as published-compilation values, not relabeled as target-
machine measurements.  The especially consequential CHF3 + F rate disagrees
with the 350 K Voloshin mechanism by 11.5x, so both source branches are
explicit and are never averaged or silently fitted.
"""
from __future__ import annotations

from dataclasses import dataclass

from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    Reaction,
    Species,
)


LIM_2014_DOI = "10.1166/jnn.2014.10171"
VOLOSHIN_2007_DOI = "10.1109/TPS.2007.906780"


@dataclass(frozen=True)
class Lim2014DaughterChemistry:
    species: tuple[Species, ...]
    reactions: tuple[Reaction, ...]
    chf3_f_rate_branch: str
    electron_rows_selected: tuple[int, ...]
    neutral_rows_selected: tuple[int, ...]
    source_gas_temperature_K: float = 700.0
    supports_target_temperature_transfer: bool = False


def lim2014_daughter_species() -> tuple[Species, ...]:
    masses = {"C": 12.011, "F": 18.998403163, "O": 15.9994}
    specifications = (
        ("C", {"C": 1}),
        ("CF4", {"C": 1, "F": 4}),
        ("CO", {"C": 1, "O": 1}),
        ("CO2", {"C": 1, "O": 2}),
        ("FO", {"F": 1, "O": 1}),
        ("CFO", {"C": 1, "F": 1, "O": 1}),
        ("CF2O", {"C": 1, "F": 2, "O": 1}),
    )
    return tuple(Species(
        name=name,
        mass_amu=sum(masses[element] * count for element, count in composition.items()),
        charge_number=0,
        composition=composition,
        role="neutral",
        source="Lim et al. 2014 Table I",
        evidence_kind="published_compilation",
    ) for name, composition in specifications)


def _electron_reactions() -> tuple[Reaction, ...]:
    # Lim Table I uses Te[eV] for these rows.  Rows 1--5 use kelvin and are
    # excluded along with all parent-feed duplicates.
    rows = (
        (6, {"e": 1, "CHF": 1}, {"e": 1, "H": 1, "CF": 1},
         9.31e-9, .204, 11.42, 3.30),
        (10, {"e": 1, "CF3": 1}, {"e": 1, "CF2": 1, "F": 1},
         6.48e-8, -.959, 11.25, 3.80),
        (11, {"e": 1, "CF2": 1}, {"e": 1, "CF": 1, "F": 1},
         8.11e-9, .386, 8.739, 5.40),
        (12, {"e": 1, "CF2": 1}, {"e": 1, "C": 1, "F": 2},
         1.39e-8, -1.164, 49.87, 11.00),
        (13, {"e": 1, "CF": 1}, {"e": 1, "C": 1, "F": 1},
         1.63e-8, -.002, 13.05, 5.60),
        (14, {"e": 1, "HF": 1}, {"e": 1, "H": 1, "F": 1},
         3.63e-8, -.313, 13.14, 5.85),
        (17, {"e": 1, "CH": 1}, {"e": 1, "C": 1, "H": 1},
         1.08e-8, -.296, 4.464, 4.34),
        (20, {"e": 1, "CO2": 1}, {"e": 1, "CO": 1, "O": 1},
         1.87e-8, 0.0, 13.89, 13.50),
        (21, {"e": 1, "CO": 1}, {"e": 1, "C": 1, "O": 1},
         1.87e-8, 0.0, 13.89, 13.50),
        (23, {"e": 1, "FO": 1}, {"e": 1, "F": 1, "O": 1},
         6.16e-9, 0.0, 4.30, 4.30),
        (24, {"e": 1, "CFO": 1}, {"e": 1, "CO": 1, "F": 1},
         8.11e-9, .386, 8.739, 5.40),
        (25, {"e": 1, "CF2O": 1}, {"e": 1, "CFO": 1, "F": 1},
         6.48e-8, -.959, 11.25, 3.80),
    )
    return tuple(Reaction(
        name=f"lim_2014_R{number}",
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
            prefactor,
            activation_eV=activation,
            temperature_power=power,
            source=f"Lim et al. 2014 Table I R{number}; Te in eV",
            evidence_kind="published_compilation",
        ),
        electron_energy_loss_eV=energy_loss,
        source=f"Lim et al. 2014 Table I R{number}",
    ) for (
        number, reactants, products, prefactor, power, activation, energy_loss
    ) in rows)


def _neutral_reactions(chf3_f_rate_branch: str) -> tuple[Reaction, ...]:
    rows = (
        (27, {"CHF3": 1, "H": 1}, {"H2": 1, "CF3": 1}, 1.60e-19),
        (28, {"CHF2": 1, "F": 1}, {"HF": 1, "CF2": 1}, 3.16e-11),
        (29, {"CHF2": 1, "H": 1}, {"HF": 1, "CHF": 1}, 3.22e-10),
        (30, {"CHF2": 1, "H": 1}, {"CF2": 1, "H2": 1}, 3.20e-14),
        (31, {"CHF2": 1, "CF3": 1}, {"CHF3": 1, "CF2": 1}, 1.58e-12),
        (32, {"CHF2": 1, "O": 1}, {"CF2O": 1, "H": 1}, 1.05e-11),
        (33, {"CHF": 1, "F": 1}, {"HF": 1, "CF": 1}, 3.25e-11),
        (34, {"CHF": 1, "H": 1}, {"HF": 1, "CH": 1}, 3.10e-10),
        (35, {"CHF": 1, "O": 1}, {"HF": 1, "CO": 1}, 3.25e-11),
        (36, {"F2": 1, "CF3": 1}, {"CF4": 1, "F": 1}, 6.31e-14),
        (37, {"F2": 1, "CF2": 1}, {"CF3": 1, "F": 1}, 7.94e-14),
        (38, {"F2": 1, "CF": 1}, {"CF2": 1, "F": 1}, 3.98e-12),
        (39, {"F2": 1, "H": 1}, {"HF": 1, "F": 1}, 8.20e-12),
        (40, {"F2": 1, "O(1d)": 1}, {"FO": 1, "F": 1}, 7.94e-12),
        (41, {"F2": 1, "CFO": 1}, {"CF2O": 1, "F": 1}, 5.01e-14),
        (42, {"CF3": 1, "F": 1}, {"CF4": 1}, 1.00e-12),
        (43, {"CF3": 1, "H": 1}, {"CF2": 1, "HF": 1}, 7.94e-11),
        (44, {"CF3": 1, "O": 1}, {"CF2O": 1, "F": 1}, 3.16e-11),
        (45, {"CF3": 1, "O(1d)": 1}, {"CF2O": 1, "F": 1}, 3.16e-11),
        (46, {"CF2": 1, "F": 1}, {"CF3": 1}, 4.17e-13),
        (47, {"CF2": 1, "H": 1}, {"HF": 1, "CF": 1}, 3.20e-11),
        (48, {"CF2": 1, "O": 1}, {"CFO": 1, "F": 1}, 3.16e-11),
        (49, {"CF2": 1, "O(1d)": 1}, {"CFO": 1, "F": 1}, 3.16e-11),
        (50, {"CF2": 1, "O": 1}, {"CO": 1, "F": 2}, 3.98e-12),
        (51, {"CF2": 1, "O(1d)": 1}, {"CO": 1, "F": 2}, 3.98e-12),
        (52, {"CF": 1, "F": 1}, {"CF2": 1}, 5.01e-15),
        (53, {"CF": 1, "H": 1}, {"C": 1, "HF": 1}, 1.20e-11),
        (54, {"CF": 1, "O": 1}, {"CO": 1, "F": 1}, 6.31e-11),
        (55, {"CF": 1, "O(1d)": 1}, {"CO": 1, "F": 1}, 2.00e-11),
        (56, {"CF": 1, "O2": 1}, {"CFO": 1, "O": 1}, 3.16e-11),
        (57, {"CH": 1, "HF": 1}, {"CF": 1, "H2": 1}, 3.23e-11),
        (58, {"CH": 1, "O": 1}, {"CO": 1, "H": 1}, 1.06e-10),
        (59, {"CH": 1, "F": 1}, {"C": 1, "HF": 1}, 1.02e-12),
        (60, {"H2": 1, "F": 1}, {"HF": 1, "H": 1}, 1.60e-11),
        (61, {"FO": 1, "O": 1}, {"F": 1, "O2": 1}, 2.51e-11),
        (62, {"FO": 1, "O(1d)": 1}, {"F": 1, "O2": 1}, 5.01e-11),
        (63, {"FO": 2}, {"F": 2, "O2": 1}, 2.51e-12),
        (64, {"FO": 2}, {"F2": 1, "O2": 1}, 2.51e-16),
        (65, {"CFO": 1, "CF3": 1}, {"CF4": 1, "CO": 1}, 1.00e-11),
        (66, {"CFO": 1, "CF3": 1}, {"CF2O": 1, "CF2": 1}, 1.00e-11),
        (67, {"CFO": 1, "CF2": 1}, {"CF3": 1, "CO": 1}, 3.16e-13),
        (68, {"CFO": 1, "CF2": 1}, {"CF2O": 1, "CF": 1}, 3.16e-13),
        (69, {"CFO": 1, "O": 1}, {"CO2": 1, "F": 1}, 1.00e-10),
        (70, {"CFO": 1, "O(1d)": 1}, {"CO2": 1, "F": 1}, 1.00e-10),
        (71, {"CFO": 2}, {"CF2O": 1, "CO": 1}, 1.00e-11),
        (72, {"CFO": 1, "F": 1}, {"CF2O": 1}, 7.94e-11),
        (73, {"CF2O": 1, "O(1d)": 1}, {"F2": 1, "CO2": 1}, 2.00e-11),
        (74, {"C": 1, "O2": 1}, {"CO": 1, "O": 1}, 1.58e-11),
        (75, {"CO": 1, "F": 1}, {"CFO": 1}, 1.29e-11),
    )
    if chf3_f_rate_branch == "lim_700K":
        first = Reaction(
            name="lim_2014_R26",
            reactants={"CHF3": 1, "F": 1},
            products={"HF": 1, "CF3": 1},
            kinetic_orders={"CHF3": 1, "F": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                1.58e-13,
                source="Lim et al. 2014 Table I R26; assumed Tg=700 K",
                evidence_kind="published_compilation"),
            electron_energy_loss_eV=0.0,
            source="Lim et al. 2014 Table I R26",
        )
    elif chf3_f_rate_branch == "voloshin_350K":
        first = Reaction(
            name="voloshin_2007_R13",
            reactants={"CHF3": 1, "F": 1},
            products={"HF": 1, "CF3": 1},
            kinetic_orders={"CHF3": 1, "F": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                1.82e-12,
                source="Voloshin et al. 2007 R13; evaluated at Tg=350 K",
                evidence_kind="published_compilation"),
            electron_energy_loss_eV=0.0,
            source="Voloshin et al. 2007 R13",
        )
    else:
        raise ValueError("chf3_f_rate_branch must be lim_700K or voloshin_350K")
    compiled = tuple(Reaction(
        name=f"lim_2014_R{number}",
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
            coefficient,
            source=(
                f"Lim et al. 2014 Table I R{number}; constants evaluated "
                "at assumed Tg=700 K"
            ),
            evidence_kind="published_compilation",
        ),
        electron_energy_loss_eV=0.0,
        source=f"Lim et al. 2014 Table I R{number}",
    ) for number, reactants, products, coefficient in rows)
    return (first, *compiled)


def build_lim_2014_daughter_chemistry(
    *, chf3_f_rate_branch: str = "voloshin_350K",
) -> Lim2014DaughterChemistry:
    electron = _electron_reactions()
    neutral = _neutral_reactions(chf3_f_rate_branch)
    return Lim2014DaughterChemistry(
        species=lim2014_daughter_species(),
        reactions=(*electron, *neutral),
        chf3_f_rate_branch=chf3_f_rate_branch,
        electron_rows_selected=(6, 10, 11, 12, 13, 14, 17, 20, 21, 23, 24, 25),
        neutral_rows_selected=tuple(range(26, 76)),
    )
