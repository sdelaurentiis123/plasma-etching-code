"""Measured secondary electron-impact ionization of CF, CF2, and CF3.

The NIST-evaluated curves in Christophorou, Olthoff, and Rao (1996) replace
the copied/Arrhenius CFx ionization rows common in older fluorocarbon global
models.  Only the six energy-resolved, atom/charge-closable channels are
executable.  The one-energy F+ anchors remain evidence, not invented curves.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
from hashlib import sha256
from pathlib import Path

from .argon import ELECTRON_MASS_AMU
from .network import (
    ElectronMaxwellianCrossSectionRateCoefficient,
    Reaction,
    ReactionNetwork,
    Species,
)


_ROOT = Path(__file__).resolve().parents[3]
DATA_DIRECTORY = (
    _ROOT / "data" / "experimental"
    / "christophorou_olthoff_rao_1996_cfx"
)
TABLE31_SHA256 = (
    "3f70148ca2cf67de4e788bb63d35d72c5aea17317b820953a7000097ef185e4f"
)
TABLE32_SHA256 = (
    "60dacd91aed11f778e4363122e79152f11e26657b3bfc3224a452dbfb7a53179"
)
TABLE33_SHA256 = (
    "a44d3a4e76eff89f352aaaf506ab0a7b2f0ef63d18413e35d316ab558cb6a1e7"
)
SOURCE = "christophorou-olthoff-rao-1996-cfx"


@dataclass(frozen=True)
class NISTCFxIonizationCurve:
    target_neutral: str
    product_ion: str
    neutral_coproducts: tuple[str, ...]
    electron_energy_eV: tuple[float, ...]
    cross_section_m2: tuple[float, ...]
    threshold_eV: float
    relative_uncertainty: float
    source_table: int

    def __post_init__(self):
        if (
            self.target_neutral not in {"CF", "CF2", "CF3"}
            or self.product_ion not in {"CF+", "CF2+", "CF3+"}
            or len(self.electron_energy_eV) < 2
            or len(self.electron_energy_eV) != len(self.cross_section_m2)
            or not 0.0 < self.relative_uncertainty < 1.0
            or self.source_table not in {31, 32, 33}
        ):
            raise ValueError("invalid NIST CFx ionization curve")

    def rate_coefficient(self) -> ElectronMaxwellianCrossSectionRateCoefficient:
        return ElectronMaxwellianCrossSectionRateCoefficient(
            electron_energy_eV=self.electron_energy_eV,
            cross_section_m2=self.cross_section_m2,
            threshold_eV=self.threshold_eV,
            relative_uncertainty=self.relative_uncertainty,
            source=(
                f"{SOURCE} Table {self.source_table}: "
                f"e+{self.target_neutral}->{self.product_ion}"
            ),
            evidence_kind="measured",
        )


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _read(path: Path, expected_sha: str) -> tuple[dict[str, str], ...]:
    if _sha(path) != expected_sha:
        raise RuntimeError(f"NIST CFx table checksum changed: {path.name}")
    with path.open(newline="", encoding="utf-8") as stream:
        return tuple(csv.DictReader(stream))


def load_nist_1996_cfx_ionization_curves() -> tuple[NISTCFxIonizationCurve, ...]:
    """Return six measured curves; deliberately exclude one-energy anchors."""
    parent = _read(
        DATA_DIRECTORY / "table31_parent_ionization.csv", TABLE31_SHA256)
    cf3_dissociation = _read(
        DATA_DIRECTORY / "table32_cf3_dissociative_ionization.csv",
        TABLE32_SHA256,
    )
    cf2_dissociation = _read(
        DATA_DIRECTORY / "table33_cf2_dissociative_ionization.csv",
        TABLE33_SHA256,
    )
    grouped: dict[tuple[int, str, str], list[dict[str, str]]] = {}
    for table, rows in (
        (31, parent), (32, cf3_dissociation), (33, cf2_dissociation),
    ):
        for row in rows:
            if row["evidence_class"].startswith("single_energy"):
                continue
            grouped.setdefault(
                (table, row["target_neutral"], row["product_ion"]), []
            ).append(row)

    thresholds = {
        ("CF3", "CF3+"): 8.9,
        ("CF2", "CF2+"): 11.4,
        ("CF", "CF+"): 9.1,
        ("CF3", "CF2+"): 17.1,
        ("CF3", "CF+"): 21.4,
        ("CF2", "CF+"): 14.3,
    }
    curves = []
    for (table, target, product), rows in grouped.items():
        coproduct = rows[0].get("neutral_coproduct", "")
        coproducts = () if not coproduct else tuple(
            (coproduct[1:] if coproduct.startswith("2") else coproduct)
            for _ in range(2 if coproduct.startswith("2") else 1)
        )
        uncertainty = {float(row["relative_uncertainty"]) for row in rows}
        if len(uncertainty) != 1 or (target, product) not in thresholds:
            raise RuntimeError("NIST CFx curve metadata changed")
        curves.append(NISTCFxIonizationCurve(
            target_neutral=target,
            product_ion=product,
            neutral_coproducts=coproducts,
            electron_energy_eV=tuple(
                float(row["electron_energy_eV"]) for row in rows),
            cross_section_m2=tuple(
                float(row["cross_section_m2"]) for row in rows),
            threshold_eV=thresholds[(target, product)],
            relative_uncertainty=uncertainty.pop(),
            source_table=table,
        ))
    return tuple(sorted(
        curves,
        key=lambda curve: (
            curve.target_neutral, curve.product_ion, curve.source_table),
    ))


def build_nist_1996_cfx_ionization_network() -> ReactionNetwork:
    """Build the measured six-channel, exactly conserved volume network."""
    carbon_mass = 12.011
    fluorine_mass = 18.998403163

    def mass(c: int, f: int) -> float:
        return c * carbon_mass + f * fluorine_mass

    species = (
        Species("e", ELECTRON_MASS_AMU, -1, {}, "electron", SOURCE),
        Species("F", fluorine_mass, 0, {"F": 1}, "neutral", SOURCE),
        Species("CF", mass(1, 1), 0, {"C": 1, "F": 1}, "neutral", SOURCE),
        Species("CF2", mass(1, 2), 0, {"C": 1, "F": 2}, "neutral", SOURCE),
        Species("CF3", mass(1, 3), 0, {"C": 1, "F": 3}, "neutral", SOURCE),
        Species("CF+", mass(1, 1), 1, {"C": 1, "F": 1}, "positive_ion", SOURCE),
        Species("CF2+", mass(1, 2), 1, {"C": 1, "F": 2}, "positive_ion", SOURCE),
        Species("CF3+", mass(1, 3), 1, {"C": 1, "F": 3}, "positive_ion", SOURCE),
    )
    reactions = []
    for curve in load_nist_1996_cfx_ionization_curves():
        products = {curve.product_ion: 1.0, "e": 2.0}
        for coproduct in curve.neutral_coproducts:
            products[coproduct] = products.get(coproduct, 0.0) + 1.0
        reactions.append(Reaction(
            name=(
                f"nist_cfx_e_{curve.target_neutral}_to_"
                f"{curve.product_ion.replace('+', 'plus')}"
            ),
            reactants={"e": 1.0, curve.target_neutral: 1.0},
            products=products,
            kinetic_orders={"e": 1.0, curve.target_neutral: 1.0},
            rate_coefficient=curve.rate_coefficient(),
            electron_energy_loss_eV=curve.threshold_eV,
            source=(
                f"{SOURCE} Table {curve.source_table}; measured curve; "
                "threshold from the same NIST review/Tarnovsky measurements"
            ),
        ))
    return ReactionNetwork(species=species, reactions=tuple(reactions))
