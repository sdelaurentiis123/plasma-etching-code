"""Daughter and heavy-particle chemistry for the Zhu NPG80 feed.

The parent CHF3, SF6, and O2 electron collisions live in
``zhu_parent_collision_chemistry`` and are intentionally absent here.  This
module installs only chemistry that begins with daughter species or closes
charged products.  Keeping the two layers separate makes double counting
machine-detectable.

SF6 daughter rates are exact transcriptions of Kokkoris et al. 2009 Table 1.
That source tabulates separate regressions for assumed Druyvesteyn and
Maxwellian EEDFs; both are exposed as sensitivity branches because neither is
identical to Petch's solved finite-frequency EEPF.  The three printed
F + SFx fall-off coefficients are excluded: they are valid at 2 Pa, whereas
the Zhu condition is 3e-2 Torr (3.9997 Pa), and the source does not print a
portable pressure law.

CHF2 daughter, hydrogen, and fluorocarbon ion-recombination rows are selected
from the conservation-checked Sandia Table-9 transcription.  The source
labels copied and estimated rows explicitly, which remains visible in every
coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .network import (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronLogTemperatureInversePolynomialRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)
from .lim_2014_chf3_oxygen_chemistry import (
    build_lim_2014_daughter_chemistry,
    lim2014_daughter_species,
)
from .sandia_chf3_mechanism import build_sandia_2001_chf3_table9_network
from .zhu_parent_collision_chemistry import zhu_parent_collision_species


KOKKORIS_2009_DOI = "10.1088/0022-3727/42/5/055209"
KOKKORIS_2009_TARGET_PRESSURE_PA = 2.0
ZHU_RECIPE_PRESSURE_PA = 3.0e-2 * 133.32236842105263
PATEAU_2014_DOI = "10.1116/1.4853675"
HUANG_2020_DOI = "10.1116/1.5125568"

_PARENT_SF6_ROWS_REPLACED_BY_MEASURED_CROSS_SECTIONS = (
    "G1", "G2", "G3", "G8", "G9", "G10", "G17", "G18",
)
_PRESSURE_SPECIFIC_ROWS_EXCLUDED = ("G35", "G36", "G37")
_SANDIA_ROWS_SELECTED = tuple(range(20, 39))


@dataclass(frozen=True)
class _KokkorisElectronRow:
    identifier: str
    reactants: dict[str, float]
    products: dict[str, float]
    energy_loss_eV: float
    druyvesteyn: tuple[float, float, float, float, float]
    maxwellian: tuple[float, float, float, float, float]


_KOKKORIS_ELECTRON_ROWS = (
    _KokkorisElectronRow(
        "G4", {"e": 1, "SF5": 1}, {"e": 1, "SF4": 1, "F": 1}, 9.6,
        (-29.36, -.2379, -14.11, -15.25, -1.204),
        (-29.57, -.2859, -13.80, 1.148, -.0781)),
    _KokkorisElectronRow(
        "G5", {"e": 1, "SF4": 1}, {"e": 1, "SF3": 1, "F": 1}, 9.6,
        (-29.36, -.2379, -14.11, -15.25, -1.204),
        (-29.57, -.2859, -13.80, 1.148, -.0781)),
    _KokkorisElectronRow(
        "G6", {"e": 1, "F2": 1}, {"e": 1, "F": 2}, 3.16,
        (-31.44, -.6986, -5.170, -1.389, -.0650),
        (-31.89, -.5549, -5.238, .4288, -.0266)),
    _KokkorisElectronRow(
        "G7", {"e": 1, "F2": 1}, {"e": 1, "F": 2}, 4.34,
        (-33.44, -.2761, -3.564, -3.946, -.0393),
        (-33.36, -.2982, -5.312, .1970, -.0124)),
    _KokkorisElectronRow(
        "G11", {"e": 1, "SF5": 1}, {"e": 2, "SF5+": 1}, 11.2,
        (-34.92, 1.487, -2.377, -29.71, -.1449),
        (-32.78, .8601, -10.76, -.0558, .0025)),
    _KokkorisElectronRow(
        "G12", {"e": 1, "SF5": 1}, {"e": 2, "SF4+": 1, "F": 1}, 14.5,
        (-36.27, 1.892, -1.387, -50.87, -.0758),
        (-33.20, 1.0177, -13.76, -.1309, .0075)),
    _KokkorisElectronRow(
        "G13", {"e": 1, "SF4": 1}, {"e": 2, "SF4+": 1}, 13.0,
        (-32.95, .8763, -10.19, -31.21, -3.989),
        (-32.01, .5939, -14.83, 2.220, -.9045)),
    _KokkorisElectronRow(
        "G14", {"e": 1, "SF4": 1}, {"e": 2, "SF3+": 1, "F": 1}, 14.5,
        (-32.75, .8222, -10.82, -40.59, -4.274),
        (-31.78, .5357, -16.26, 1.974, -.7729)),
    _KokkorisElectronRow(
        "G15", {"e": 1, "SF3": 1}, {"e": 2, "SF3+": 1}, 11.0,
        (-35.55, 1.750, -2.086, -28.70, -.1357),
        (-33.23, 1.073, -10.36, -.1016, .0055)),
    _KokkorisElectronRow(
        "G16", {"e": 1, "F2": 1}, {"e": 2, "F2+": 1}, 15.69,
        (-35.60, 1.467, -6.140, -57.14, -.4860),
        (-33.38, .8249, -15.96, .0655, -.0041)),
    _KokkorisElectronRow(
        "G19", {"e": 1, "F2": 1}, {"F": 1, "F-": 1}, 0.0,
        (-33.31, -1.487, -.2795, .0109, -.0004),
        (-32.81, -1.440, -.5283, .0558, -.0028)),
)


@dataclass(frozen=True)
class _KokkorisConstantRow:
    identifier: str
    reactants: dict[str, float]
    products: dict[str, float]
    log_rate_m3_s: float


def _constant_rows() -> tuple[_KokkorisConstantRow, ...]:
    rows: list[_KokkorisConstantRow] = []
    for partner in ("SF6", "SF5", "SF4", "SF3", "F", "F2"):
        detachment_products = {"F": 1, partner: 1, "e": 1}
        if partner == "F":
            detachment_products = {"F": 2, "e": 1}
        rows.append(_KokkorisConstantRow(
            f"G20_{partner}", {"F-": 1, partner: 1},
            detachment_products, -44.39))
        products = {"SF6": 1, partner: 1, "e": 1}
        if partner == "SF6":
            products = {"SF6": 2, "e": 1}
        rows.append(_KokkorisConstantRow(
            f"G21_{partner}", {"SF6-": 1, partner: 1},
            products, -44.98))
    rows.extend((
        _KokkorisConstantRow(
            "G38", {"F2": 1, "SF5": 1},
            {"SF6": 1, "F": 1}, -46.41),
        _KokkorisConstantRow(
            "G39", {"F2": 1, "SF4": 1},
            {"SF5": 1, "F": 1}, -46.41),
        _KokkorisConstantRow(
            "G40", {"F2": 1, "SF3": 1},
            {"SF4": 1, "F": 1}, -46.41),
        _KokkorisConstantRow(
            "G41", {"SF5": 2}, {"SF6": 1, "SF4": 1}, -41.50),
    ))
    neutral_parents = {
        "SF5+": "SF5", "SF4+": "SF4", "SF3+": "SF3", "F2+": "F2",
    }
    for positive, neutral in neutral_parents.items():
        for negative, negative_neutral in (("SF6-", "SF6"), ("F-", "F")):
            rows.append(_KokkorisConstantRow(
                f"G42_49_{positive}_{negative}",
                {positive: 1, negative: 1},
                {neutral: 1, negative_neutral: 1}, -29.93))
    rows.append(_KokkorisConstantRow(
        "G50", {"SF5+": 1, "SF6": 1},
        {"SF3+": 1, "F2": 1, "SF6": 1}, -39.65))
    return tuple(rows)


@dataclass(frozen=True)
class ZhuSupplementalChemistry:
    network: ReactionNetwork
    kokkoris_eedf_shape: str
    parent_sf6_rows_replaced: tuple[str, ...]
    pressure_specific_rows_excluded: tuple[str, ...]
    sandia_rows_selected: tuple[int, ...]
    chf3_f_rate_branch: str
    electron_collision_rows_replaced: tuple[str, ...] = ()
    supports_measured_parent_eedf: bool = True
    supports_complete_daughter_eedf: bool = False
    supports_target_pressure_falloff: bool = False
    supports_oxygen_daughter_chemistry: bool = True
    supports_sf6_o2_titration_chemistry: bool = True
    supports_chf3_neutral_chain: bool = True
    supports_complete_cross_ion_recombination: bool = False


def zhu_reactor_species() -> tuple[Species, ...]:
    """Return parent products plus the Pateau O/S/F coupling products."""
    pateau_source = "pateau-2014-sf6-o2 Table I and Tables V/VII"
    extra = (
        Species(
            name="O(1d)", mass_amu=15.9994, charge_number=0,
            composition={"O": 1}, role="excited_neutral",
            source=pateau_source, evidence_kind="published_compilation"),
        Species(
            name="O+", mass_amu=15.9994, charge_number=1,
            composition={"O": 1}, role="positive_ion",
            source=pateau_source, evidence_kind="published_compilation"),
        Species(
            name="H2", mass_amu=2.01568, charge_number=0,
            composition={"H": 2}, role="neutral",
            source="sandia-2001-fluorocarbon-mechanisms Table 9",
            evidence_kind="published_compilation"),
        Species(
            name="H2+", mass_amu=2.01568, charge_number=1,
            composition={"H": 2}, role="positive_ion",
            source="sandia-2001-fluorocarbon-mechanisms Table 9",
            evidence_kind="published_compilation"),
        Species(
            name="H+", mass_amu=1.00784, charge_number=1,
            composition={"H": 1}, role="positive_ion",
            source="sandia-2001-fluorocarbon-mechanisms Table 9",
            evidence_kind="published_compilation"),
        Species(
            name="HF+", mass_amu=20.006243163, charge_number=1,
            composition={"H": 1, "F": 1}, role="positive_ion",
            source=(
                "Huang-2020 HCl-threshold-shifted HF ionization closure"
            ),
            evidence_kind="semi_empirical"),
    )
    oxygenated = (
        ("CH", {"C": 1, "H": 1}),
        ("SO", {"S": 1, "O": 1}),
        ("SO2", {"S": 1, "O": 2}),
        ("SOF", {"S": 1, "O": 1, "F": 1}),
        ("SOF2", {"S": 1, "O": 1, "F": 2}),
        ("SOF3", {"S": 1, "O": 1, "F": 3}),
        ("SOF4", {"S": 1, "O": 1, "F": 4}),
        ("SO2F", {"S": 1, "O": 2, "F": 1}),
        ("SO2F2", {"S": 1, "O": 2, "F": 2}),
    )
    atomic_mass = {
        "C": 12.011, "H": 1.00784, "S": 32.065,
        "O": 15.9994, "F": 18.998403163,
    }
    return (
        *zhu_parent_collision_species(),
        *extra,
        *(Species(
            name=name,
            mass_amu=sum(
                atomic_mass[element] * count
                for element, count in composition.items()),
            charge_number=0,
            composition=composition,
            role="neutral",
            source=(
                "sandia-2001-fluorocarbon-mechanisms Table 9"
                if name == "CH" else pateau_source
            ),
            evidence_kind="published_compilation",
        ) for name, composition in oxygenated),
        *lim2014_daughter_species(),
    )


def _pateau_oxygen_reactions() -> tuple[Reaction, ...]:
    electron_rows = (
        (75, {"e": 1, "O": 1}, {"e": 1, "O(1d)": 1},
         4.47e-9, 0.0, 2.286, 1.96),
        (76, {"e": 1, "O": 1}, {"e": 2, "O+": 1},
         9.0e-9, .7, 13.6, 13.62),
        (77, {"e": 1, "O(1d)": 1}, {"e": 2, "O+": 1},
         9.0e-9, .7, 11.6, 11.65),
        (78, {"e": 1, "O2+": 1}, {"O": 2},
         5.2e-9, -1.0, 0.0, 0.0),
        (79, {"e": 1, "O-": 1}, {"e": 2, "O": 1},
         2.0e-7, 0.0, 5.5, 5.5),
    )
    reactions = [Reaction(
        name=f"pateau_2014_R{number}",
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
            prefactor,
            activation_eV=activation,
            temperature_power=power,
            source=f"pateau-2014-sf6-o2 Table V R{number}",
            evidence_kind="published_compilation",
        ),
        electron_energy_loss_eV=energy_loss,
        source=f"pateau-2014-sf6-o2 Table V R{number}",
    ) for (
        number, reactants, products, prefactor, power, activation, energy_loss
    ) in electron_rows]
    heavy_rows = (
        (80, {"O(1d)": 1, "O2": 1}, {"O": 1, "O2": 1}, 4.11e-11),
        (81, {"O(1d)": 1, "O": 1}, {"O": 2}, 8.1e-12),
        (82, {"O+": 1, "O2": 1}, {"O": 1, "O2+": 1}, 2.0e-11),
        (83, {"O": 1, "O-": 1}, {"O2": 1, "e": 1}, 3.0e-10),
        (84, {"O-": 1, "O2+": 1}, {"O": 1, "O2": 1}, 1.5e-7),
        (85, {"O-": 1, "O+": 1}, {"O": 2}, 2.5e-7),
    )
    reactions.extend(Reaction(
        name=f"pateau_2014_R{number}",
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
            coefficient,
            source=f"pateau-2014-sf6-o2 Table V R{number}",
            evidence_kind="published_compilation",
        ),
        electron_energy_loss_eV=0.0,
        source=f"pateau-2014-sf6-o2 Table V R{number}",
    ) for number, reactants, products, coefficient in heavy_rows)
    return tuple(reactions)


def _pateau_titration_reactions() -> tuple[Reaction, ...]:
    base_rows = (
        (119, "SF3", "SOF2", "F", 2.0e-11),
        (120, "SF2", "SOF", "F", 1.1e-10),
        (121, "SF", "SO", "F", 1.7e-10),
        (122, "SOF", "SO2", "F", 7.9e-11),
        (123, "SOF3", "SO2F2", "F", 5.0e-11),
    )
    reactions = []
    for number, reactant, product, coproduct, coefficient in base_rows:
        for oxygen in ("O", "O(1d)"):
            reactions.append(Reaction(
                name=f"pateau_2014_R{number}_{oxygen}",
                reactants={reactant: 1, oxygen: 1},
                products={product: 1, coproduct: 1},
                kinetic_orders={reactant: 1, oxygen: 1},
                rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                    coefficient,
                    source=f"pateau-2014-sf6-o2 Table VII R{number}",
                    evidence_kind="published_compilation"),
                electron_energy_loss_eV=0.0,
                source=f"pateau-2014-sf6-o2 Table VII R{number}",
            ))
    simple_rows = (
        (124, {"SOF3": 1, "F": 1}, {"SOF4": 1}, 5.2e-14),
        (125, {"SOF2": 1, "F": 1}, {"SOF3": 1}, 5.2e-14),
        (126, {"SOF": 1, "F": 1}, {"SOF2": 1}, 1.0e-13),
        (127, {"SO": 1, "F": 1}, {"SOF": 1}, 1.0e-14),
        (128, {"SO2": 1, "F": 1}, {"SO2F": 1}, 1.0e-13),
        (129, {"SO2F": 1, "F": 1}, {"SO2F2": 1}, 1.0e-11),
    )
    reactions.extend(Reaction(
        name=f"pateau_2014_R{number}",
        reactants=reactants,
        products=products,
        kinetic_orders=reactants,
        rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
            coefficient,
            source=f"pateau-2014-sf6-o2 Table VII R{number}",
            evidence_kind="published_compilation"),
        electron_energy_loss_eV=0.0,
        source=f"pateau-2014-sf6-o2 Table VII R{number}",
    ) for number, reactants, products, coefficient in simple_rows)
    continued_rows = (
        (130, "SO", {"SO2": 1}, 1.4e-13),
        (131, "SOF2", {"SO2F2": 1}, 1.0e-15),
        (132, "SF5", {"SOF4": 1, "F": 1}, 1.0e-12),
        (133, "SF4", {"SOF4": 1}, 1.0e-14),
        (134, "SF3", {"SOF3": 1}, 1.0e-10),
        (135, "SF2", {"SOF2": 1}, 1.0e-10),
        (136, "SF", {"SOF": 1}, 1.0e-10),
        (137, "SO2F2", {"SOF2": 1, "O2": 1}, 1.0e-12),
    )
    for number, reactant, products, coefficient in continued_rows:
        for oxygen in ("O", "O(1d)"):
            reactions.append(Reaction(
                name=f"pateau_2014_R{number}_{oxygen}",
                reactants={reactant: 1, oxygen: 1},
                products=products,
                kinetic_orders={reactant: 1, oxygen: 1},
                rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                    coefficient,
                    source=f"pateau-2014-sf6-o2 Table VII R{number}",
                    evidence_kind="published_compilation"),
                electron_energy_loss_eV=0.0,
                source=f"pateau-2014-sf6-o2 Table VII R{number}",
            ))
    return tuple(reactions)


def _pateau_sf6_charge_closure_reactions() -> tuple[Reaction, ...]:
    """Fill exact Table-III R38/R50 pairs absent from Kokkoris's subset."""
    reactions = []
    negative_parents = {
        "SF6-": "SF6", "SF5-": "SF5", "SF4-": "SF4",
        "SF3-": "SF3", "SF2-": "SF2", "F-": "F", "F2-": "F2",
    }
    neutral_partners = ("SF6", "SF5", "SF4", "SF3", "SF2", "SF", "S", "F", "F2")
    kokkoris_detachment = {
        (negative, partner)
        for negative in ("SF6-", "F-")
        for partner in ("SF6", "SF5", "SF4", "SF3", "F", "F2")
    }
    for negative, neutral in negative_parents.items():
        for partner in neutral_partners:
            if (negative, partner) in kokkoris_detachment:
                continue
            products = {neutral: 1, partner: 1, "e": 1}
            if neutral == partner:
                products = {neutral: 2, "e": 1}
            reactions.append(Reaction(
                name=f"pateau_2014_R38_{negative}_{partner}",
                reactants={negative: 1, partner: 1},
                products=products,
                kinetic_orders={negative: 1, partner: 1},
                rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                    5.27e-14,
                    source="pateau-2014-sf6-o2 Table III R38",
                    evidence_kind="published_compilation"),
                electron_energy_loss_eV=0.0,
                source="pateau-2014-sf6-o2 Table III R38",
            ))
    positive_parents = {
        "SF5+": "SF5", "SF4+": "SF4", "SF3+": "SF3",
        "SF2+": "SF2", "SF+": "SF", "S+": "S", "F+": "F", "F2+": "F2",
    }
    kokkoris_recombination = {
        (positive, negative)
        for positive in ("SF5+", "SF4+", "SF3+", "F2+")
        for negative in ("SF6-", "F-")
    }
    for positive, positive_neutral in positive_parents.items():
        for negative, negative_neutral in negative_parents.items():
            if (positive, negative) in kokkoris_recombination:
                continue
            products = {positive_neutral: 1, negative_neutral: 1}
            if positive_neutral == negative_neutral:
                products = {positive_neutral: 2}
            reactions.append(Reaction(
                name=f"pateau_2014_R50_{positive}_{negative}",
                reactants={positive: 1, negative: 1},
                products=products,
                kinetic_orders={positive: 1, negative: 1},
                rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                    1.0e-7,
                    source="pateau-2014-sf6-o2 Table III R50",
                    evidence_kind="published_compilation"),
                electron_energy_loss_eV=0.0,
                source="pateau-2014-sf6-o2 Table III R50",
            ))
    return tuple(reactions)


def _huang_hf_charge_closure_reactions() -> tuple[Reaction, ...]:
    """Close HF+ with the estimated rates printed by Huang et al. 2020.

    These rows are not promoted to measurements: Huang labels the
    dissociative-recombination and mutual-neutralization coefficients
    ``est.`` in its appendix.  They are still preferable to silently giving
    the newly produced HF+ ion only a material-wall sink.
    """

    source = (
        f"Huang et al. 2020 appendix ({HUANG_2020_DOI}); printed as estimated"
    )
    reactions = [Reaction(
        name="huang_2020_hfplus_e_dissociative_recombination",
        reactants={"e": 1, "HF+": 1},
        products={"H": 1, "F": 1},
        kinetic_orders={"e": 1, "HF+": 1},
        rate_coefficient=ElectronArrheniusRateCoefficient.from_cm3_per_s(
            1.0e-7,
            activation_eV=0.0,
            temperature_power=-0.5,
            source=source,
            evidence_kind="estimated",
        ),
        # The printed -10.1 eV is reaction enthalpy, not a measured electron
        # energy-loss moment, so it is not inserted into the electron ledger.
        electron_energy_loss_eV=0.0,
        source=source,
    )]
    for negative, neutral in (("F-", "F"), ("O-", "O")):
        products = {neutral: 1, "H": 1, "F": 1}
        if neutral == "F":
            products = {"F": 2, "H": 1}
        reactions.append(Reaction(
            name=f"huang_2020_hfplus_{negative}_mutual_neutralization",
            reactants={negative: 1, "HF+": 1},
            products=products,
            kinetic_orders={negative: 1, "HF+": 1},
            rate_coefficient=ConstantRateCoefficient.from_cm3_per_s(
                2.0e-7, source=source, evidence_kind="estimated"
            ),
            electron_energy_loss_eV=0.0,
            source=source,
        ))
    return tuple(reactions)


def _electron_reaction(row: _KokkorisElectronRow, shape: str) -> Reaction:
    coefficients = getattr(row, shape)
    return Reaction(
        name=f"kokkoris_2009_{row.identifier}_{shape}",
        reactants=row.reactants,
        products=row.products,
        kinetic_orders=row.reactants,
        rate_coefficient=(
            ElectronLogTemperatureInversePolynomialRateCoefficient.from_log_si(
                coefficients[0],
                temperature_power=coefficients[1],
                inverse_temperature_coefficients=coefficients[2:],
                source=(
                    f"Kokkoris et al. 2009 Table 1 {row.identifier}; "
                    f"{shape} EEDF regression"
                ),
                evidence_kind="regressed",
            )
        ),
        electron_energy_loss_eV=row.energy_loss_eV,
        source=(
            f"Kokkoris et al. 2009 Table 1 {row.identifier}; "
            "daughter-only transcription"
        ),
    )


def _constant_reaction(row: _KokkorisConstantRow) -> Reaction:
    return Reaction(
        name=f"kokkoris_2009_{row.identifier}",
        reactants=row.reactants,
        products=row.products,
        kinetic_orders=row.reactants,
        rate_coefficient=ConstantRateCoefficient(
            value_si=float(np.exp(row.log_rate_m3_s)),
            density_order=2.0,
            source=f"Kokkoris et al. 2009 Table 1 {row.identifier}",
            source_units="m^3 s^-1",
            evidence_kind="published_compilation",
        ),
        electron_energy_loss_eV=0.0,
        source=f"Kokkoris et al. 2009 Table 1 {row.identifier}",
    )


def build_zhu_supplemental_chemistry(
    *,
    kokkoris_eedf_shape: str = "druyvesteyn",
    chf3_f_rate_branch: str = "voloshin_350K",
    electron_collision_rows_replaced: tuple[str, ...] = (),
) -> ZhuSupplementalChemistry:
    """Return the non-parent volume network for a declared EEDF closure."""
    if kokkoris_eedf_shape not in {"druyvesteyn", "maxwellian"}:
        raise ValueError("kokkoris_eedf_shape must be druyvesteyn or maxwellian")
    kokkoris = [
        *(_electron_reaction(row, kokkoris_eedf_shape)
          for row in _KOKKORIS_ELECTRON_ROWS),
        *(_constant_reaction(row) for row in _constant_rows()),
    ]
    sandia = build_sandia_2001_chf3_table9_network()
    selected_names = {
        f"sandia_chf3_table9_{number:02d}" for number in _SANDIA_ROWS_SELECTED
    }
    selected_sandia = [
        reaction for reaction in sandia.reactions
        if reaction.name in selected_names
    ]
    if len(selected_sandia) != len(selected_names):
        raise RuntimeError("Sandia daughter-row selection is incomplete")
    lim = build_lim_2014_daughter_chemistry(
        chf3_f_rate_branch=chf3_f_rate_branch)
    all_reactions = tuple((
        *kokkoris,
        *selected_sandia,
        *lim.reactions,
        *_pateau_oxygen_reactions(),
        *_pateau_titration_reactions(),
        *_pateau_sf6_charge_closure_reactions(),
        *_huang_hf_charge_closure_reactions(),
    ))
    replaced = tuple(str(name) for name in electron_collision_rows_replaced)
    if len(set(replaced)) != len(replaced) or any(not name for name in replaced):
        raise ValueError("invalid replaced supplemental electron rows")
    available = {reaction.name: reaction for reaction in all_reactions}
    if set(replaced) - set(available):
        raise ValueError(
            "requested replacement names absent supplemental reactions"
        )
    if any("e" not in available[name].kinetic_orders for name in replaced):
        raise ValueError("only electron rows may be collision-deck replaced")
    network = ReactionNetwork(
        species=zhu_reactor_species(),
        reactions=tuple(
            reaction for reaction in all_reactions
            if reaction.name not in set(replaced)
        ),
    )
    return ZhuSupplementalChemistry(
        network=network,
        kokkoris_eedf_shape=kokkoris_eedf_shape,
        parent_sf6_rows_replaced=_PARENT_SF6_ROWS_REPLACED_BY_MEASURED_CROSS_SECTIONS,
        pressure_specific_rows_excluded=_PRESSURE_SPECIFIC_ROWS_EXCLUDED,
        sandia_rows_selected=_SANDIA_ROWS_SELECTED,
        chf3_f_rate_branch=lim.chf3_f_rate_branch,
        electron_collision_rows_replaced=replaced,
    )
