"""Source-backed positive-ion production for a C4F6 reactor.

The C4F6 electron-collision deck provides an aggregate parent-ionization
cross section.  It does not identify the product ion.  This module partitions
an externally computed aggregate parent-ionization event rate using the
checksum-locked NIST 70 eV mass spectrum, while retaining every measured
monoisotopic heavy product.  It then joins the independently measured NIST
CF/CF2/CF3 secondary-ionization network.

The 70 eV spectrum is a branching *prior*, not an energy-resolved branching
measurement.  Consequently the resulting source ledger is suitable for
testing reaction topology and for a declared uncertainty ensemble.  It is not
a steady reactor state, a wafer flux, or a Krueger boundary: ion-neutral
conversion, negative ions, Bohm/wall loss, residence time, and surface return
remain outside this component.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
from hashlib import sha256
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .cfx_electron_ionization import build_nist_1996_cfx_ionization_network
from .network import RateContext


_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NIST_SPECTRUM = (
    _ROOT / "data" / "experimental" / "nist_c4f6_mass_spectrum"
    / "electron_ionization_sticks.csv"
)
NIST_SPECTRUM_SHA256 = (
    "2115e1a887f53fcc99423ab81fee94fde461f00c2e0da30ac188e95748c7a665"
)
_FORMULA = re.compile(r"C(?P<carbon>[0-9]*)(?:F(?P<fluorine>[0-9]*))?$")


def _finite_nonnegative_mapping(values: Mapping[str, float]):
    converted = {str(name): float(value) for name, value in values.items()}
    if (
        any(not name.strip() for name in converted)
        or any(not math.isfinite(value) or value < 0.0
               for value in converted.values())
    ):
        raise ValueError("ion-source rates must be finite and nonnegative")
    return MappingProxyType(converted)


def _composition_from_assignment(assignment: str) -> Mapping[str, int]:
    formula = str(assignment).removesuffix("+")
    match = _FORMULA.fullmatch(formula)
    if match is None:
        raise ValueError(f"unsupported C4F6 ion assignment {assignment!r}")
    carbon_text = match.group("carbon")
    fluorine_text = match.group("fluorine")
    carbon = 1 if carbon_text == "" else int(carbon_text)
    fluorine = (
        0 if fluorine_text is None
        else (1 if fluorine_text == "" else int(fluorine_text))
    )
    if not 0 < carbon <= 4 or not 0 <= fluorine <= 6:
        raise ValueError("C4F6 fragment composition leaves the parent support")
    return MappingProxyType({"C": carbon, "F": fluorine})


@dataclass(frozen=True)
class C4F6DirectIonBranch:
    """One monoisotopic product in the NIST 70 eV EI spectrum."""

    product_ion: str
    nominal_m_over_z: int
    composition: Mapping[str, int]
    relative_intensity_percent: float
    normalized_branch_fraction: float
    digitization_absolute_intensity_bound_percent: float

    def __post_init__(self):
        composition = {
            str(name): int(value) for name, value in self.composition.items()
        }
        values = np.asarray([
            self.nominal_m_over_z,
            self.relative_intensity_percent,
            self.normalized_branch_fraction,
            self.digitization_absolute_intensity_bound_percent,
        ], dtype=float)
        if (
            not str(self.product_ion).endswith("+")
            or composition != dict(_composition_from_assignment(
                self.product_ion))
            or np.any(~np.isfinite(values))
            or self.nominal_m_over_z <= 0
            or self.relative_intensity_percent <= 0.0
            or not 0.0 < self.normalized_branch_fraction < 1.0
            or self.digitization_absolute_intensity_bound_percent <= 0.0
        ):
            raise ValueError("invalid direct C4F6 ion branch")
        object.__setattr__(self, "nominal_m_over_z", int(self.nominal_m_over_z))
        object.__setattr__(self, "composition", MappingProxyType(composition))


@dataclass(frozen=True)
class C4F6DirectIonBranchSet:
    """Complete retained monoisotopic NIST partition and its authority."""

    branches: tuple[C4F6DirectIonBranch, ...]
    source_sha256: str
    source_energy_eV: float = 70.0
    supports_energy_resolved_branching: bool = False
    supports_absolute_reactor_flux: bool = False
    supports_krueger_boundary: bool = False

    def __post_init__(self):
        names = tuple(branch.product_ion for branch in self.branches)
        fractions = np.asarray([
            branch.normalized_branch_fraction for branch in self.branches
        ])
        if (
            len(self.branches) != 15
            or len(names) != len(set(names))
            or "C3F3+" not in names
            or "C4F6+" not in names
            or not np.isclose(np.sum(fractions), 1.0, rtol=0.0, atol=2.0e-15)
            or self.source_sha256 != NIST_SPECTRUM_SHA256
            or self.source_energy_eV != 70.0
            or self.supports_energy_resolved_branching
            or self.supports_absolute_reactor_flux
            or self.supports_krueger_boundary
        ):
            raise ValueError("invalid direct C4F6 branch set")

    def partition_event_rate_m3_s(
        self, aggregate_parent_ionization_rate_m3_s: float,
    ) -> Mapping[str, float]:
        rate = float(aggregate_parent_ionization_rate_m3_s)
        if not math.isfinite(rate) or rate < 0.0:
            raise ValueError("aggregate parent ionization rate must be nonnegative")
        return MappingProxyType({
            branch.product_ion: rate * branch.normalized_branch_fraction
            for branch in self.branches
        })


def load_nist_c4f6_direct_ion_branches(
    path: str | Path = DEFAULT_NIST_SPECTRUM,
) -> C4F6DirectIonBranchSet:
    """Load the checksum-pinned spectrum and retain all monoisotopic ions."""
    source = Path(path)
    raw = source.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != NIST_SPECTRUM_SHA256:
        raise RuntimeError("NIST C4F6 EI spectrum checksum changed")
    rows = tuple(csv.DictReader(raw.decode("utf-8").splitlines()))
    retained = tuple(
        row for row in rows
        if row["assignment_class"] != "isotope_satellite"
    )
    if len(rows) != 23 or len(retained) != 15:
        raise RuntimeError("NIST C4F6 EI spectrum topology changed")
    intensities = np.asarray([
        float(row["relative_intensity_percent"]) for row in retained
    ])
    total = float(np.sum(intensities))
    branches = tuple(
        C4F6DirectIonBranch(
            product_ion=row["assignment"],
            nominal_m_over_z=int(row["nominal_m_over_z"]),
            composition=_composition_from_assignment(row["assignment"]),
            relative_intensity_percent=float(row["relative_intensity_percent"]),
            normalized_branch_fraction=float(intensity / total),
            digitization_absolute_intensity_bound_percent=float(
                row["digitization_absolute_intensity_bound_percent"]),
        )
        for row, intensity in zip(retained, intensities)
    )
    return C4F6DirectIonBranchSet(branches=branches, source_sha256=digest)


@dataclass(frozen=True)
class C4F6PositiveIonSourceLedger:
    """Direct-plus-secondary volume sources before ion chemistry or loss."""

    direct_parent_sources_m3_s: Mapping[str, float]
    secondary_cfx_sources_m3_s: Mapping[str, float]
    combined_sources_m3_s: Mapping[str, float]
    parent_branching_energy_eV: float
    known_missing_operators: tuple[str, ...]
    supports_steady_reactor_composition: bool = False
    supports_wafer_flux: bool = False
    supports_krueger_boundary: bool = False

    def __post_init__(self):
        direct = _finite_nonnegative_mapping(self.direct_parent_sources_m3_s)
        secondary = _finite_nonnegative_mapping(self.secondary_cfx_sources_m3_s)
        combined = _finite_nonnegative_mapping(self.combined_sources_m3_s)
        expected = dict(direct)
        for name, value in secondary.items():
            expected[name] = expected.get(name, 0.0) + value
        if (
            set(combined) != set(expected)
            or any(not np.isclose(combined[name], value, rtol=2.0e-15, atol=0.0)
                   for name, value in expected.items())
            or self.parent_branching_energy_eV != 70.0
            or not self.known_missing_operators
            or self.supports_steady_reactor_composition
            or self.supports_wafer_flux
            or self.supports_krueger_boundary
        ):
            raise ValueError("invalid C4F6 positive-ion source ledger")
        object.__setattr__(self, "direct_parent_sources_m3_s", direct)
        object.__setattr__(self, "secondary_cfx_sources_m3_s", secondary)
        object.__setattr__(self, "combined_sources_m3_s", combined)


class C4F6PositiveIonSourceModel:
    """Join direct parent fragmentation to measured secondary CFx ionization."""

    def __init__(self):
        self.direct = load_nist_c4f6_direct_ion_branches()
        self.secondary = build_nist_1996_cfx_ionization_network()

    def evaluate(
        self,
        *,
        aggregate_parent_ionization_rate_m3_s: float,
        electron_density_m3: float,
        neutral_cfx_densities_m3: Mapping[str, float],
        electron_temperature_eV: float,
        gas_temperature_K: float | None = None,
    ) -> C4F6PositiveIonSourceLedger:
        electron_density = float(electron_density_m3)
        temperature = float(electron_temperature_eV)
        neutral = {
            str(name): float(value)
            for name, value in neutral_cfx_densities_m3.items()
        }
        if (
            set(neutral) != {"CF", "CF2", "CF3"}
            or not math.isfinite(electron_density)
            or electron_density < 0.0
            or not math.isfinite(temperature)
            or temperature <= 0.0
            or any(not math.isfinite(value) or value < 0.0
                   for value in neutral.values())
        ):
            raise ValueError("invalid CFx secondary-ionization state")
        densities = {
            "e": electron_density,
            "F": 0.0,
            **neutral,
            "CF+": 0.0,
            "CF2+": 0.0,
            "CF3+": 0.0,
        }
        context = RateContext(
            electron_temperature_eV=temperature,
            gas_temperature_K=gas_temperature_K,
        )
        event_rates = self.secondary.event_rates_m3_s(densities, context)
        secondary_sources = {"CF+": 0.0, "CF2+": 0.0, "CF3+": 0.0}
        for reaction, event_rate in zip(self.secondary.reactions, event_rates):
            for ion in secondary_sources:
                secondary_sources[ion] += (
                    float(reaction.products.get(ion, 0.0)) * float(event_rate)
                )
        direct_sources = dict(self.direct.partition_event_rate_m3_s(
            aggregate_parent_ionization_rate_m3_s
        ))
        combined = dict(direct_sources)
        for name, value in secondary_sources.items():
            combined[name] = combined.get(name, 0.0) + value
        return C4F6PositiveIonSourceLedger(
            direct_parent_sources_m3_s=direct_sources,
            secondary_cfx_sources_m3_s=secondary_sources,
            combined_sources_m3_s=combined,
            parent_branching_energy_eV=self.direct.source_energy_eV,
            known_missing_operators=(
                "energy-resolved C4F6 primary branching",
                "neutral-fragment production and loss",
                "ion-neutral conversion and charge exchange",
                "negative-ion production and mutual neutralization",
                "Bohm, wall, exhaust, and surface-return losses",
            ),
        )
