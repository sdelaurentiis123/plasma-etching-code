"""Conservation-checked Sandia CHF3 gas-mechanism extension.

This module transcribes Table 9 of Ho et al., SAND2001-1292.  The table is
valuable because it closes the hydrogen-bearing CHF3 branches that are absent
from CFx-only mechanisms.  It is not promoted to a modern evaluated data set:
the source itself labels the daughter-CHF2 rates as copies and the ion-ion and
dissociative-recombination rates as estimates.

The original electron-rate fit is ``k = A T**B exp(-C/T)`` with electron
temperature in kelvin and ``A`` in the corresponding molecule-cm-s units.
Petch evaluates electron rates in eV, so the conversion is performed exactly
at construction time.  No reactor or etched-profile datum enters the deck.
"""
from __future__ import annotations

from dataclasses import dataclass

from .argon import ELECTRON_MASS_AMU
from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    GasTemperatureArrheniusRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)


SANDIA_REPORT_DOI = "10.2172/782704"
SANDIA_REPORT_SHA256 = (
    "3e401d0d5c5ffb0308767bffb2d4d952b1c6110f79f95c03b781e927342a4d8f"
)
KELVIN_PER_EV = 11604.518121550082


@dataclass(frozen=True)
class SandiaTable9Row:
    number: int
    reactants: dict[str, float]
    products: dict[str, float]
    A_cm3_s_K_minus_B: float
    B: float
    C_K: float
    electron_energy_loss_eV: float
    evidence_kind: str
    note: str


def _row(
    number: int,
    reactants: dict[str, float],
    products: dict[str, float],
    A: float,
    B: float,
    C: float,
    loss: float,
    evidence_kind: str,
    note: str,
) -> SandiaTable9Row:
    return SandiaTable9Row(
        number, reactants, products, A, B, C, loss, evidence_kind, note)


# Exact Table-9 transcription. Repeated rows are retained because they
# represent separate cross-section features in the source mechanism.
SANDIA_TABLE9_ROWS = (
    _row(1, {"e": 1, "CHF3": 1}, {"e": 1, "CHF3": 1}, 1.089e-2, -1.2214, 2.645e4, .37, "regressed", "V1,4; ref 53"),
    _row(2, {"e": 1, "CHF3": 1}, {"e": 1, "CHF3": 1}, 3.019e-2, -1.3964, 1.759e4, .18, "regressed", "V2,5; ref 53"),
    _row(3, {"e": 1, "CHF3": 1}, {"e": 1, "CHF3": 1}, 5.218e-2, -1.4396, 1.406e4, .13, "regressed", "V3,6; ref 53"),
    _row(4, {"e": 1, "CHF3": 1}, {"e": 1, "CHF3": 1}, 4.122e-2, -1.3790, 1.549e5, 10.9, "estimated", "electronic excitation estimate"),
    _row(5, {"e": 1, "CHF3": 1}, {"e": 1, "CHF3": 1}, 1.584e-11, .4322, 1.415e5, 11.9, "estimated", "electronic excitation estimate"),
    _row(6, {"e": 1, "CHF3": 1}, {"CF3+": 1, "H": 1, "e": 2}, 5.361e-15, 1.3438, 1.736e5, 15.2, "regressed", "ref 53"),
    _row(7, {"e": 1, "CHF3": 1}, {"CHF2+": 1, "F": 1, "e": 2}, 2.225e-15, 1.2886, 1.906e5, 16.8, "regressed", "ref 53"),
    _row(8, {"e": 1, "CHF3": 1}, {"CF2+": 1, "HF": 1, "e": 2}, 6.533e-17, 1.4404, 1.800e5, 17.6, "regressed", "ref 53"),
    _row(9, {"e": 1, "CHF3": 1}, {"CF+": 1, "F": 2, "H": 1, "e": 2}, 6.780e-16, 1.5225, 2.335e5, 20.9, "regressed", "ref 53"),
    _row(10, {"e": 1, "CHF3": 1}, {"F+": 1, "CHF2": 1, "e": 2}, 8.120e-14, .9194, 4.273e5, 37.0, "regressed", "ref 53"),
    _row(11, {"e": 1, "CHF3": 1}, {"CHF+": 1, "F": 2, "e": 2}, 1.006e-15, 1.3223, 2.120e5, 19.8, "regressed", "ref 53"),
    _row(12, {"e": 1, "CHF3": 1}, {"CF3": 1, "H": 1, "e": 1}, 3.963e-17, 1.6416, 1.044e5, 11.0, "regressed", "ref 53"),
    _row(13, {"e": 1, "CHF3": 1}, {"CHF2": 1, "F": 1, "e": 1}, 1.187e-16, 1.3167, 1.376e5, 13.0, "regressed", "ref 53"),
    _row(14, {"e": 1, "CHF3": 1}, {"CF2": 1, "HF": 1, "e": 1}, 3.626e-14, 1.0759, 2.631e5, 23.6, "regressed", "ref 53"),
    _row(15, {"e": 1, "CHF3": 1}, {"CHF": 1, "F": 2, "e": 1}, 8.084e-13, .5725, 4.070e5, 35.0, "regressed", "ref 53"),
    _row(16, {"e": 1, "CHF3": 1}, {"CF": 1, "H": 1, "F": 2, "e": 1}, 6.752e-9, .1877, 2.358e5, 13.3, "regressed", "ref 53"),
    _row(17, {"e": 1, "CHF3": 1}, {"CF3": 1, "H": 1, "e": 1}, 2.555e-6, -.4365, 1.546e5, 11.0, "regressed", "ref 53 add-on"),
    _row(18, {"e": 1, "CHF3": 1}, {"F-": 1, "CHF2": 1}, 8.988e-5, -1.3618, 1.129e5, 1.3, "regressed", "ref 53; low-energy attachment feature"),
    _row(19, {"e": 1, "CHF3": 1}, {"F-": 1, "CHF2": 1}, 1.166e-6, -1.2306, 4.219e4, 1.3, "regressed", "ref 53; second attachment feature"),
    _row(20, {"e": 1, "CHF2": 1}, {"CF2+": 1, "H": 1, "e": 2}, 5.361e-15, 1.3438, 1.736e5, 17.2, "estimated", "copied from CHF3"),
    _row(21, {"e": 1, "CHF2": 1}, {"CHF+": 1, "F": 1, "e": 2}, 2.225e-15, 1.2886, 1.906e5, 14.3, "estimated", "copied from CHF3"),
    _row(22, {"e": 1, "CHF2": 1}, {"CF+": 1, "HF": 1, "e": 2}, 6.533e-17, 1.4404, 1.800e5, 14.6, "estimated", "copied from CHF3"),
    _row(23, {"e": 1, "CHF2": 1}, {"CF2": 1, "H": 1, "e": 1}, 3.963e-17, 1.6416, 1.044e5, 2.71, "estimated", "copied from CHF3"),
    _row(24, {"e": 1, "CHF2": 1}, {"CHF": 1, "F": 1, "e": 1}, 1.187e-16, 1.3167, 1.376e5, 4.75, "estimated", "copied from CHF3"),
    _row(25, {"e": 1, "CHF2": 1}, {"CF": 1, "H": 1, "F": 1, "e": 1}, 3.626e-14, 1.0759, 2.631e5, 8.09, "estimated", "copied from CHF3"),
    _row(26, {"e": 1, "H2": 1}, {"e": 1, "H2": 1}, 1.400e-5, -.7604, 2.264e4, .5, "regressed", "vibrational; ref 54"),
    _row(27, {"e": 1, "H2": 1}, {"e": 1, "H2": 1}, 4.185e-12, .6434, 1.407e5, 11.37, "regressed", "electronic B1; ref 54"),
    _row(28, {"e": 1, "H2": 1}, {"e": 1, "H2": 1}, 6.250e-13, .8498, 1.702e5, 11.70, "regressed", "electronic C1; ref 54"),
    _row(29, {"e": 1, "H2": 1}, {"H": 2, "e": 1}, 1.697e-8, -.0244, 1.199e5, 4.4, "regressed", "ref 54"),
    _row(30, {"e": 1, "H2": 1}, {"H2+": 1, "e": 2}, 1.329e-13, 1.0750, 1.976e5, 15.4, "regressed", "ref 54"),
    _row(31, {"e": 1, "H": 1}, {"e": 1, "H": 1}, 8.367e-10, .3014, 1.335e5, 10.2, "regressed", "electronic 2P; ref 54"),
    _row(32, {"e": 1, "H": 1}, {"H+": 1, "e": 2}, 7.332e-12, .6938, 1.694e5, 13.6, "regressed", "ref 54"),
    _row(33, {"F-": 1, "CHF2+": 1}, {"F": 1, "CF2": 1, "H": 1}, 4.0e-7, -.5, 0.0, 0.0, "estimated", "ion-ion neutralization"),
    _row(34, {"F-": 1, "CHF+": 1}, {"F": 1, "CF": 1, "H": 1}, 4.0e-7, -.5, 0.0, 0.0, "estimated", "ion-ion neutralization"),
    _row(35, {"F-": 1, "H2+": 1}, {"F": 1, "H": 2}, 4.0e-7, -.5, 0.0, 0.0, "estimated", "ion-ion neutralization"),
    _row(36, {"e": 1, "H2+": 1}, {"H": 2}, 4.0e-8, 0.0, 0.0, 0.0, "estimated", "dissociative recombination"),
    _row(37, {"e": 1, "CHF2+": 1}, {"H": 1, "CF2": 1}, 4.0e-8, 0.0, 0.0, 0.0, "estimated", "dissociative recombination"),
    _row(38, {"e": 1, "CHF+": 1}, {"H": 1, "CF": 1}, 4.0e-8, 0.0, 0.0, 0.0, "estimated", "dissociative recombination"),
)


def _species() -> tuple[Species, ...]:
    masses = {"H": 1.00784, "C": 12.011, "F": 18.998403163}

    def mass(composition: dict[str, int]) -> float:
        return sum(masses[element] * count for element, count in composition.items())

    specifications = (
        ("CHF3", {"C": 1, "H": 1, "F": 3}, 0, "neutral"),
        ("CHF2", {"C": 1, "H": 1, "F": 2}, 0, "neutral"),
        ("CHF", {"C": 1, "H": 1, "F": 1}, 0, "neutral"),
        ("CF3", {"C": 1, "F": 3}, 0, "neutral"),
        ("CF2", {"C": 1, "F": 2}, 0, "neutral"),
        ("CF", {"C": 1, "F": 1}, 0, "neutral"),
        ("H2", {"H": 2}, 0, "neutral"),
        ("HF", {"H": 1, "F": 1}, 0, "neutral"),
        ("H", {"H": 1}, 0, "neutral"),
        ("F", {"F": 1}, 0, "neutral"),
        ("F-", {"F": 1}, -1, "negative_ion"),
        ("CF3+", {"C": 1, "F": 3}, 1, "positive_ion"),
        ("CHF2+", {"C": 1, "H": 1, "F": 2}, 1, "positive_ion"),
        ("CF2+", {"C": 1, "F": 2}, 1, "positive_ion"),
        ("CHF+", {"C": 1, "H": 1, "F": 1}, 1, "positive_ion"),
        ("CF+", {"C": 1, "F": 1}, 1, "positive_ion"),
        ("F+", {"F": 1}, 1, "positive_ion"),
        ("H2+", {"H": 2}, 1, "positive_ion"),
        ("H+", {"H": 1}, 1, "positive_ion"),
    )
    return (
        Species(
            name="e", mass_amu=ELECTRON_MASS_AMU, charge_number=-1,
            composition={}, role="electron", source="CODATA electron mass",
        ),
        *(Species(
            name=name, mass_amu=mass(composition), charge_number=charge,
            composition=composition, role=role,
            source="sandia-2001-fluorocarbon-mechanisms Table 9",
            evidence_kind=("estimated" if name != "CHF3" else "measured"),
        ) for name, composition, charge, role in specifications),
    )


def _electron_rate(row: SandiaTable9Row):
    prefactor = row.A_cm3_s_K_minus_B * KELVIN_PER_EV ** row.B
    return ElectronArrheniusRateCoefficient.from_cm3_per_s(
        prefactor,
        activation_eV=row.C_K / KELVIN_PER_EV,
        temperature_power=row.B,
        source=(
            f"sandia-2001-fluorocarbon-mechanisms Table 9 row {row.number}; "
            "source kelvin fit converted exactly to Te[eV]"
        ),
        evidence_kind=row.evidence_kind,
    )


def _heavy_rate(row: SandiaTable9Row):
    return GasTemperatureArrheniusRateCoefficient(
        prefactor_si=row.A_cm3_s_K_minus_B * 1.0e-6,
        temperature_power=row.B,
        activation_temperature_K=-row.C_K,
        reference_temperature_K=1.0,
        density_order=2.0,
        source=(
            f"sandia-2001-fluorocarbon-mechanisms Table 9 row {row.number}"
        ),
        source_units="cm^3 molecule^-1 s^-1; T in K",
        evidence_kind=row.evidence_kind,
    )


def build_sandia_2001_chf3_table9_network() -> ReactionNetwork:
    """Return all 38 Table-9 reactions with exact atom/charge linting."""
    reactions = []
    for row in SANDIA_TABLE9_ROWS:
        electron_driven = "e" in row.reactants
        reactions.append(Reaction(
            name=f"sandia_chf3_table9_{row.number:02d}",
            reactants=row.reactants,
            products=row.products,
            kinetic_orders=row.reactants,
            rate_coefficient=(
                _electron_rate(row) if electron_driven else _heavy_rate(row)
            ),
            electron_energy_loss_eV=row.electron_energy_loss_eV,
            source=(
                f"sandia-2001-fluorocarbon-mechanisms Table 9 row "
                f"{row.number}: {row.note}"
            ),
        ))
    return ReactionNetwork(species=_species(), reactions=tuple(reactions))


def sandia_table9_evidence_counts() -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in SANDIA_TABLE9_ROWS:
        counts[row.evidence_kind] = counts.get(row.evidence_kind, 0) + 1
    return counts
