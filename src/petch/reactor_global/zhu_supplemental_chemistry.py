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

CHF2 daughter and fluorocarbon ion-recombination rows are selected from the
conservation-checked Sandia Table-9 transcription.  The source labels those
rows as copied or estimated, which remains visible in every coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .network import (
    ConstantRateCoefficient,
    ElectronLogTemperatureInversePolynomialRateCoefficient,
    Reaction,
    ReactionNetwork,
)
from .sandia_chf3_mechanism import build_sandia_2001_chf3_table9_network
from .zhu_parent_collision_chemistry import zhu_parent_collision_species


KOKKORIS_2009_DOI = "10.1088/0022-3727/42/5/055209"
KOKKORIS_2009_TARGET_PRESSURE_PA = 2.0
ZHU_RECIPE_PRESSURE_PA = 3.0e-2 * 133.32236842105263

_PARENT_SF6_ROWS_REPLACED_BY_MEASURED_CROSS_SECTIONS = (
    "G1", "G2", "G3", "G8", "G9", "G10", "G17", "G18",
)
_PRESSURE_SPECIFIC_ROWS_EXCLUDED = ("G35", "G36", "G37")
_SANDIA_ROWS_SELECTED = (20, 21, 22, 23, 24, 25, 33, 34, 37, 38)


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
    supports_measured_parent_eedf: bool = True
    supports_complete_daughter_eedf: bool = False
    supports_target_pressure_falloff: bool = False
    supports_complete_oxygen_heavy_chemistry: bool = False


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
    *, kokkoris_eedf_shape: str = "druyvesteyn",
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
    network = ReactionNetwork(
        species=zhu_parent_collision_species(),
        reactions=tuple((*kokkoris, *selected_sandia)),
    )
    return ZhuSupplementalChemistry(
        network=network,
        kokkoris_eedf_shape=kokkoris_eedf_shape,
        parent_sf6_rows_replaced=_PARENT_SF6_ROWS_REPLACED_BY_MEASURED_CROSS_SECTIONS,
        pressure_specific_rows_excluded=_PRESSURE_SPECIFIC_ROWS_EXCLUDED,
        sandia_rows_selected=_SANDIA_ROWS_SELECTED,
    )
