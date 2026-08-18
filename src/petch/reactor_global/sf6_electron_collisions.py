"""NIST-evaluated SF6 electron transport and aggregate chemistry.

The source review distinguishes recommended, suggested, and deduced tables.
The constructed Boltzmann deck preserves those evidence classes and exposes
every closure needed to turn totals into one deterministic process set.  It
is a transport/primary-source input, not a unique reactor or depth result.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from .argon import ELECTRON_MASS_AMU
from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


SF6_MASS_AMU = 146.055418326
NIST_SF6_TABLE9_CSV_SHA256 = (
    "5c22dee68eec50f3bb1a0dd33c57400d6f40d599a30abbbe7cd1671da4ee75bf"
)
NIST_SF6_TABLE14_CSV_SHA256 = (
    "dc702c0a1651aba096dc4dddec8bcc16dfe97028b8ec9091e9de06f681a07e5e"
)
NIST_SF6_TABLE15_CSV_SHA256 = (
    "4d8579dbd68e2f11e386c8f38fd9e969760fcc625b6f5ae427e3bf92f024b3db"
)
NIST_SF6_TABLE16_CSV_SHA256 = (
    "3bf597105de4294cd677cd5cbe16947ca852113327eac5713914c1c9f3487063"
)
NIST_SF6_TABLE17_CSV_SHA256 = (
    "9566c88d59c4ab1c39833850dec3ef2535fcb8ccbb89c30612f3a0fc770d7663"
)
NIST_SF6_TABLE20_CSV_SHA256 = (
    "a4a7ed0c6d4227f8fbfaf5e69407d88ce03688c2b12c8a96427f90f5158c42f1"
)
NIST_SF6_TABLE24_CSV_SHA256 = (
    "000ef875fc230f5609088b52b7792a3c23aaf7288f0e21e9714d20ea69861390"
)
NIST_SF6_TABLE25_CSV_SHA256 = (
    "06a85f8cede5ce9d5e9b2a35f22aa2d9f69a77551662578702813127269ff063"
)
NIST_SF6_TABLE27_CSV_SHA256 = (
    "335d8a16ce8af84e6f3bf3bc5edbe11c73156d331ff8b1b991b993c2b131a1bb"
)
NIST_SF6_TABLE28_CSV_SHA256 = (
    "e8ee625160731eb829864f458164c4432355a1b38840eec187e737f6450e76de"
)
NIST_SF6_TABLE36_CSV_SHA256 = (
    "007006f5c347b8986f949ef6bfaf4e3902b3cedf314a5262514437e6338df5c6"
)
NIST_SF6_TABLE35_CSV_SHA256 = (
    "cfac7c0e58834a986995805616916c1c7423692fab7d2cc0a910ea05d9bc50c4"
)
NIST_SF6_TABLE37_CSV_SHA256 = (
    "fa459460ad6d98b6835b26d9d4adf22c72c6b885cc92cbd714da5a0df4e5ab03"
)
_ROOT = Path(__file__).resolve().parents[3]
_DATA = (
    _ROOT / "data" / "experimental"
    / "christophorou_olthoff_2000_sf6"
)


@dataclass(frozen=True)
class NISTSF6CrossSectionCurve:
    electron_energy_eV: np.ndarray
    cross_section_m2: np.ndarray
    evidence_class: str
    table_number: int

    def __post_init__(self):
        energy = np.asarray(self.electron_energy_eV, dtype=float).copy()
        sigma = np.asarray(self.cross_section_m2, dtype=float).copy()
        if (
            energy.ndim != 1
            or energy.size < 10
            or sigma.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(sigma))
            or np.any(energy < 0.0)
            or np.any(sigma < 0.0)
            or np.any(np.diff(energy) <= 0.0)
            or self.table_number not in {9, 14, 15, 17, 20, 28}
            or not self.evidence_class
        ):
            raise ValueError("invalid NIST SF6 cross-section curve")
        energy.setflags(write=False)
        sigma.setflags(write=False)
        object.__setattr__(self, "electron_energy_eV", energy)
        object.__setattr__(self, "cross_section_m2", sigma)


@dataclass(frozen=True)
class NISTSF6ProductCrossSectionCurve:
    """One product-resolved source curve or one-energy anchor."""

    product: str
    electron_energy_eV: np.ndarray
    cross_section_m2: np.ndarray
    evidence_class: str
    table_number: int

    def __post_init__(self):
        energy = np.asarray(self.electron_energy_eV, dtype=float).copy()
        sigma = np.asarray(self.cross_section_m2, dtype=float).copy()
        if (
            not str(self.product).strip()
            or energy.ndim != 1
            or energy.size < 1
            or sigma.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(sigma))
            or np.any(energy < 0.0)
            or np.any(sigma < 0.0)
            or np.any(np.diff(energy) <= 0.0)
            or self.table_number not in {16, 24, 25, 27}
            or not self.evidence_class
        ):
            raise ValueError("invalid NIST SF6 product cross-section curve")
        energy.setflags(write=False)
        sigma.setflags(write=False)
        object.__setattr__(self, "product", str(self.product).strip())
        object.__setattr__(self, "electron_energy_eV", energy)
        object.__setattr__(self, "cross_section_m2", sigma)


@dataclass(frozen=True)
class NISTSF6DriftCurve:
    reduced_electric_field_Td: np.ndarray
    drift_velocity_m_s: np.ndarray
    recommendation_class: tuple[str, ...]
    gas_temperature_K: str = "293-300"

    def __post_init__(self):
        field = np.asarray(self.reduced_electric_field_Td, dtype=float).copy()
        velocity = np.asarray(self.drift_velocity_m_s, dtype=float).copy()
        grade = tuple(self.recommendation_class)
        if (
            field.ndim != 1
            or field.size != 29
            or velocity.shape != field.shape
            or len(grade) != field.size
            or np.any(~np.isfinite(field))
            or np.any(~np.isfinite(velocity))
            or np.any(field < 0.0)
            or np.any(velocity < 0.0)
            or np.any(np.diff(field) <= 0.0)
            or set(grade) != {"deduced", "recommended", "suggested"}
            or self.gas_temperature_K != "293-300"
        ):
            raise ValueError("invalid NIST SF6 drift curve")
        field.setflags(write=False)
        velocity.setflags(write=False)
        object.__setattr__(self, "reduced_electric_field_Td", field)
        object.__setattr__(self, "drift_velocity_m_s", velocity)
        object.__setattr__(self, "recommendation_class", grade)

    @property
    def recommended_mask(self) -> np.ndarray:
        value = np.asarray([
            item == "recommended" for item in self.recommendation_class
        ])
        value.setflags(write=False)
        return value


@dataclass(frozen=True)
class NISTSF6EffectiveIonizationCurve:
    reduced_electric_field_Td: np.ndarray
    effective_ionization_coefficient_m2: np.ndarray

    def __post_init__(self):
        field = np.asarray(self.reduced_electric_field_Td, dtype=float).copy()
        value = np.asarray(
            self.effective_ionization_coefficient_m2, dtype=float
        ).copy()
        if (
            field.ndim != 1
            or field.size != 24
            or value.shape != field.shape
            or np.any(~np.isfinite(field))
            or np.any(~np.isfinite(value))
            or np.any(field <= 0.0)
            or np.any(np.diff(field) <= 0.0)
        ):
            raise ValueError("invalid NIST SF6 effective-ionization curve")
        field.setflags(write=False)
        value.setflags(write=False)
        object.__setattr__(self, "reduced_electric_field_Td", field)
        object.__setattr__(
            self, "effective_ionization_coefficient_m2", value
        )


@dataclass(frozen=True)
class NISTSF6AttachmentRateCurve:
    reduced_electric_field_Td: np.ndarray
    attachment_rate_coefficient_m3_s: np.ndarray

    def __post_init__(self):
        field = np.asarray(self.reduced_electric_field_Td, dtype=float).copy()
        value = np.asarray(
            self.attachment_rate_coefficient_m3_s, dtype=float
        ).copy()
        if (
            field.ndim != 1
            or field.size != 12
            or value.shape != field.shape
            or np.any(~np.isfinite(field))
            or np.any(~np.isfinite(value))
            or np.any(field <= 0.0)
            or np.any(value <= 0.0)
            or np.any(np.diff(field) <= 0.0)
        ):
            raise ValueError("invalid NIST SF6 attachment-rate curve")
        field.setflags(write=False)
        value.setflags(write=False)
        object.__setattr__(self, "reduced_electric_field_Td", field)
        object.__setattr__(
            self, "attachment_rate_coefficient_m3_s", value
        )


@dataclass(frozen=True)
class NISTEvaluatedSF6Replay:
    derived_deck: ElectronCollisionDeck
    source_curves: tuple[NISTSF6CrossSectionCurve, ...]
    vibrational_energy_loss_eV: float
    maximum_energy_eV: float
    evidence_class: str = "nist_evaluated_aggregate_sf6"
    supports_direct_transport_constraints: bool = True
    supports_resolved_primary_chemistry: bool = False
    supports_unique_reactor_state: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            len(self.source_curves) != 6
            or self.vibrational_energy_loss_eV not in {0.095, 0.117}
            or self.maximum_energy_eV != 200.0
            or self.evidence_class != "nist_evaluated_aggregate_sf6"
            or not self.supports_direct_transport_constraints
            or self.supports_resolved_primary_chemistry
            or self.supports_unique_reactor_state
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid NIST-evaluated SF6 replay")


@dataclass(frozen=True)
class NISTProductResolvedSF6Replay:
    """Evaluated totals split into auditable ion and anion products.

    Attachment is resolved directly from Tables 24, 25, and 27 subject only
    to rounding reconciliation against Table 28.  Positive-ion *curves* are
    not tabulated by the review: Table 16 supplies one-energy anchors, and
    Table 3 supplies appearance thresholds.  Those branches therefore use a
    declared threshold-onset closure and are renormalized to Table 17 at each
    energy.  The total collision probability remains exactly evaluated data.
    """

    aggregate_replay: NISTEvaluatedSF6Replay
    derived_deck: ElectronCollisionDeck
    ionization_source_anchors: tuple[NISTSF6ProductCrossSectionCurve, ...]
    attachment_source_curves: tuple[NISTSF6ProductCrossSectionCurve, ...]
    ionization_threshold_eV: tuple[tuple[str, float], ...]
    ionization_anchor_sum_m2: float
    evaluated_total_ionization_at_100eV_m2: float
    maximum_attachment_rounding_rescale_fraction: float
    maximum_sf5_source_peak_normalized_residual: float
    evidence_class: str = "nist_product_resolved_sf6_with_declared_ion_closure"
    supports_direct_attachment_products: bool = True
    supports_direct_positive_ion_curves: bool = False
    supports_unique_reactor_state: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        thresholds = dict(self.ionization_threshold_eV)
        if (
            len(self.ionization_source_anchors) != 9
            or len(self.attachment_source_curves) != 7
            or set(thresholds) != {
                "SF5+", "SF4+", "SF3+", "SF2+", "SF+", "S+", "F+",
                "SF4++", "SF2++",
            }
            or any(value <= 0.0 for value in thresholds.values())
            or self.ionization_anchor_sum_m2 <= 0.0
            or self.evaluated_total_ionization_at_100eV_m2 <= 0.0
            or not 0.0 <= self.maximum_attachment_rounding_rescale_fraction < .02
            or not 0.0 <= self.maximum_sf5_source_peak_normalized_residual < .05
            or self.evidence_class
            != "nist_product_resolved_sf6_with_declared_ion_closure"
            or not self.supports_direct_attachment_products
            or self.supports_direct_positive_ion_curves
            or self.supports_unique_reactor_state
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid NIST product-resolved SF6 replay")


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_curve(
    filename: str,
    *,
    expected_sha256: str,
    evidence_class: str,
    table_number: int,
) -> NISTSF6CrossSectionCurve:
    path = _DATA / filename
    if _digest(path) != expected_sha256:
        raise RuntimeError(f"NIST SF6 Table {table_number} CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(row["evidence_class"] != evidence_class for row in rows):
        raise RuntimeError(f"NIST SF6 Table {table_number} metadata changed")
    return NISTSF6CrossSectionCurve(
        electron_energy_eV=np.asarray([
            float(row["electron_energy_eV"]) for row in rows
        ]),
        cross_section_m2=np.asarray([
            float(row["cross_section_m2"]) for row in rows
        ]),
        evidence_class=evidence_class,
        table_number=table_number,
    )


def load_nist_2000_sf6_total_scattering():
    return _load_curve(
        "table9_total_scattering.csv",
        expected_sha256=NIST_SF6_TABLE9_CSV_SHA256,
        evidence_class="recommended_total_scattering",
        table_number=9,
    )


def load_nist_2000_sf6_momentum_transfer():
    return _load_curve(
        "table14_momentum_transfer.csv",
        expected_sha256=NIST_SF6_TABLE14_CSV_SHA256,
        evidence_class="suggested_elastic_momentum_transfer",
        table_number=14,
    )


def load_nist_2000_sf6_vibrational_excitation():
    return _load_curve(
        "table15_vibrational_excitation.csv",
        expected_sha256=NIST_SF6_TABLE15_CSV_SHA256,
        evidence_class="deduced_total_vibrational_excitation",
        table_number=15,
    )


def load_nist_2000_sf6_total_ionization():
    return _load_curve(
        "table17_total_ionization.csv",
        expected_sha256=NIST_SF6_TABLE17_CSV_SHA256,
        evidence_class="recommended_total_ionization",
        table_number=17,
    )


def load_nist_2000_sf6_total_neutral_dissociation():
    return _load_curve(
        "table20_total_neutral_dissociation.csv",
        expected_sha256=NIST_SF6_TABLE20_CSV_SHA256,
        evidence_class=(
            "deduced_total_neutral_dissociation_requires_confirmation"
        ),
        table_number=20,
    )


def load_nist_2000_sf6_total_attachment():
    return _load_curve(
        "table28_total_attachment.csv",
        expected_sha256=NIST_SF6_TABLE28_CSV_SHA256,
        evidence_class="suggested_room_temperature_total_attachment",
        table_number=28,
    )


def load_nist_2000_sf6_partial_ionization_anchors(
) -> tuple[NISTSF6ProductCrossSectionCurve, ...]:
    """Load the preferred-analysis partial-ion values at exactly 100 eV."""

    path = _DATA / "table16_partial_ionization_100eV.csv"
    if _digest(path) != NIST_SF6_TABLE16_CSV_SHA256:
        raise RuntimeError("NIST SF6 Table 16 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    evidence = "measured_case_A_preferred_analysis"
    if len(rows) != 9 or any(
        row["evidence_class"] != evidence
        or float(row["electron_energy_eV"]) != 100.0
        for row in rows
    ):
        raise RuntimeError("NIST SF6 Table 16 metadata changed")
    return tuple(
        NISTSF6ProductCrossSectionCurve(
            product=row["product_ion"],
            electron_energy_eV=np.asarray([float(row["electron_energy_eV"])]),
            cross_section_m2=np.asarray([float(row["cross_section_m2"])]),
            evidence_class=evidence,
            table_number=16,
        )
        for row in rows
    )


def _load_attachment_product_curve(
    filename: str,
    *,
    expected_sha256: str,
    product: str,
    evidence_class: str,
    table_number: int,
) -> NISTSF6ProductCrossSectionCurve:
    path = _DATA / filename
    if _digest(path) != expected_sha256:
        raise RuntimeError(f"NIST SF6 Table {table_number} CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(row["evidence_class"] != evidence_class for row in rows):
        raise RuntimeError(f"NIST SF6 Table {table_number} metadata changed")
    return NISTSF6ProductCrossSectionCurve(
        product=product,
        electron_energy_eV=np.asarray([
            float(row["electron_energy_eV"]) for row in rows
        ]),
        cross_section_m2=np.asarray([
            float(row["cross_section_m2"]) for row in rows
        ]),
        evidence_class=evidence_class,
        table_number=table_number,
    )


def load_nist_2000_sf6_partial_attachment_curves(
) -> tuple[NISTSF6ProductCrossSectionCurve, ...]:
    """Load all seven reported SF6 attachment product curves.

    Blank minor-anion cells in Table 27 are removed rather than interpreted
    as measured zero.  The reconciliation routine below supplies the explicit
    between-support closure against the evaluated total.
    """

    major = (
        _load_attachment_product_curve(
            "table24_sf6_minus_attachment.csv",
            expected_sha256=NIST_SF6_TABLE24_CSV_SHA256,
            product="SF6-",
            evidence_class="recommended_nondissociative_SF6_minus",
            table_number=24,
        ),
        _load_attachment_product_curve(
            "table25_sf5_minus_attachment.csv",
            expected_sha256=NIST_SF6_TABLE25_CSV_SHA256,
            product="SF5-",
            evidence_class="suggested_dissociative_SF5_minus",
            table_number=25,
        ),
    )
    path = _DATA / "table27_minor_anion_attachment.csv"
    if _digest(path) != NIST_SF6_TABLE27_CSV_SHA256:
        raise RuntimeError("NIST SF6 Table 27 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    evidence = "suggested_dissociative_partial_anions"
    if any(row["evidence_class"] != evidence for row in rows):
        raise RuntimeError("NIST SF6 Table 27 metadata changed")
    columns = (
        ("SF4-", "SF4_minus_cross_section_m2"),
        ("SF3-", "SF3_minus_cross_section_m2"),
        ("SF2-", "SF2_minus_cross_section_m2"),
        ("F2-", "F2_minus_cross_section_m2"),
        ("F-", "F_minus_cross_section_m2"),
    )
    minor = []
    for product, column in columns:
        reported = tuple(row for row in rows if row[column].strip())
        minor.append(NISTSF6ProductCrossSectionCurve(
            product=product,
            electron_energy_eV=np.asarray([
                float(row["electron_energy_eV"]) for row in reported
            ]),
            cross_section_m2=np.asarray([
                float(row[column]) for row in reported
            ]),
            evidence_class=evidence,
            table_number=27,
        ))
    return (*major, *minor)


def load_nist_2000_sf6_drift_curve() -> NISTSF6DriftCurve:
    path = _DATA / "table36_drift_velocity.csv"
    if _digest(path) != NIST_SF6_TABLE36_CSV_SHA256:
        raise RuntimeError("NIST SF6 Table 36 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(row["gas_temperature_K"] != "293-300" for row in rows):
        raise RuntimeError("NIST SF6 Table 36 metadata changed")
    return NISTSF6DriftCurve(
        reduced_electric_field_Td=np.asarray([
            float(row["reduced_electric_field_Td"]) for row in rows
        ]),
        drift_velocity_m_s=np.asarray([
            float(row["drift_velocity_m_s"]) for row in rows
        ]),
        recommendation_class=tuple(
            row["recommendation_class"] for row in rows
        ),
    )


def load_nist_2000_sf6_effective_ionization_curve(
) -> NISTSF6EffectiveIonizationCurve:
    path = _DATA / "table35_effective_ionization.csv"
    if _digest(path) != NIST_SF6_TABLE35_CSV_SHA256:
        raise RuntimeError("NIST SF6 Table 35 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(row["evidence_class"] != "recommended_swarm_fit" for row in rows):
        raise RuntimeError("NIST SF6 Table 35 metadata changed")
    return NISTSF6EffectiveIonizationCurve(
        reduced_electric_field_Td=np.asarray([
            float(row["reduced_electric_field_Td"]) for row in rows
        ]),
        effective_ionization_coefficient_m2=np.asarray([
            float(row["effective_ionization_coefficient_m2"])
            for row in rows
        ]),
    )


def load_nist_2000_sf6_attachment_rate_curve(
) -> NISTSF6AttachmentRateCurve:
    path = _DATA / "table37_attachment_rate.csv"
    if _digest(path) != NIST_SF6_TABLE37_CSV_SHA256:
        raise RuntimeError("NIST SF6 Table 37 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(
        row["evidence_class"]
        != "assessed_from_eta_over_n_times_drift"
        for row in rows
    ):
        raise RuntimeError("NIST SF6 Table 37 metadata changed")
    return NISTSF6AttachmentRateCurve(
        reduced_electric_field_Td=np.asarray([
            float(row["reduced_electric_field_Td"]) for row in rows
        ]),
        attachment_rate_coefficient_m3_s=np.asarray([
            float(row["attachment_rate_coefficient_m3_s"])
            for row in rows
        ]),
    )


def _with_threshold_and_tail(
    curve: NISTSF6CrossSectionCurve,
    *,
    threshold_eV: float,
    tail_energy_eV: float | None = None,
    tail_cross_section_m2: float | None = None,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    energy = curve.electron_energy_eV.copy()
    sigma = curve.cross_section_m2.copy()
    keep = energy > threshold_eV
    energy = energy[keep]
    sigma = sigma[keep]
    energy = np.concatenate(([0.0, threshold_eV], energy))
    sigma = np.concatenate(([0.0, 0.0], sigma))
    if tail_energy_eV is not None:
        if tail_energy_eV <= energy[-1]:
            raise ValueError("tail energy must exceed the source support")
        energy = np.concatenate((energy, [tail_energy_eV]))
        sigma = np.concatenate((sigma, [float(tail_cross_section_m2)]))
    return tuple(energy), tuple(sigma)


def derive_nist_evaluated_sf6_replay(
    *,
    vibrational_energy_loss_eV: float = 0.095,
    maximum_energy_eV: float = 200.0,
) -> NISTEvaluatedSF6Replay:
    """Build an aggregate SF6 set with explicit total-table closures.

    Below 2.5 eV, Table 9 total scattering is deconvolved by subtracting the
    same total-vibration and total-attachment rows that the solver re-adds.
    Table 14 supplies elastic momentum transfer from 2.75 eV upward.  This
    prevents the common double count of total scattering plus inelastic rows.
    """

    loss = float(vibrational_energy_loss_eV)
    maximum = float(maximum_energy_eV)
    if loss not in {0.095, 0.117}:
        raise ValueError("SF6 aggregate vibrational loss must be 0.095 or 0.117 eV")
    if maximum != 200.0:
        raise ValueError("the audited SF6 deck currently supports exactly 200 eV")
    total = load_nist_2000_sf6_total_scattering()
    momentum = load_nist_2000_sf6_momentum_transfer()
    vibration = load_nist_2000_sf6_vibrational_excitation()
    ionization = load_nist_2000_sf6_total_ionization()
    dissociation = load_nist_2000_sf6_total_neutral_dissociation()
    attachment = load_nist_2000_sf6_total_attachment()

    vib_energy, vib_sigma = _with_threshold_and_tail(
        vibration,
        threshold_eV=loss,
        tail_energy_eV=15.0,
        tail_cross_section_m2=0.0,
    )
    vib_energy = tuple((*vib_energy, maximum))
    vib_sigma = tuple((*vib_sigma, 0.0))
    attach_energy = tuple(np.concatenate((
        [0.0], attachment.electron_energy_eV, [20.0, maximum],
    )))
    attach_sigma = tuple(np.concatenate((
        [attachment.cross_section_m2[0]],
        attachment.cross_section_m2,
        [0.0, 0.0],
    )))

    low = total.electron_energy_eV <= 2.5
    low_energy = total.electron_energy_eV[low]
    low_elastic = total.cross_section_m2[low].copy()
    low_elastic -= np.interp(low_energy, vib_energy, vib_sigma)
    low_elastic -= np.interp(low_energy, attach_energy, attach_sigma)
    if np.any(low_elastic <= 0.0):
        raise RuntimeError("SF6 low-energy total-scattering deconvolution failed")
    momentum_energy = tuple(np.concatenate((
        [0.0], low_energy, momentum.electron_energy_eV,
    )))
    momentum_sigma = tuple(np.concatenate((
        [low_elastic[0]], low_elastic, momentum.cross_section_m2,
    )))

    diss_energy, diss_sigma = _with_threshold_and_tail(
        dissociation, threshold_eV=9.6)
    ion_energy, ion_sigma = _with_threshold_and_tail(
        ionization, threshold_eV=15.5)
    processes = (
        ElectronCollisionProcess(
            kind="MOMENTUM",
            target="SF6",
            product="SF6",
            electron_energy_eV=momentum_energy,
            cross_section_m2=momentum_sigma,
            mass_ratio=ELECTRON_MASS_AMU / SF6_MASS_AMU,
            comments=(
                "NIST Tables 9/14 evaluated momentum backbone",
                "Table 9 deconvolved through 2.5 eV; Table 14 from 2.75 eV",
                "constant low-energy elastic proxy from first Table 9 residual",
            ),
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            target="SF6",
            product="SF6(v) aggregate",
            electron_energy_eV=vib_energy,
            cross_section_m2=vib_sigma,
            energy_loss_eV=loss,
            comments=(
                "NIST Table 15 deduced total vibration",
                f"aggregate effective loss={loss} eV",
                "linear closure to zero at 15 eV",
            ),
        ),
        ElectronCollisionProcess(
            kind="EXCITATION",
            target="SF6",
            product="SFx + F aggregate",
            electron_energy_eV=diss_energy,
            cross_section_m2=diss_sigma,
            energy_loss_eV=9.6,
            comments=(
                "NIST Table 20 deduced total neutral dissociation",
                "aggregate branching unresolved; Table 21 minimum threshold",
            ),
        ),
        ElectronCollisionProcess(
            kind="IONIZATION",
            target="SF6",
            product="SFx+ aggregate",
            electron_energy_eV=ion_energy,
            cross_section_m2=ion_sigma,
            energy_loss_eV=15.5,
            comments=(
                "NIST Table 17 recommended total ionization",
                "positive-ion fragment branching unresolved",
            ),
        ),
        ElectronCollisionProcess(
            kind="ATTACHMENT",
            target="SF6",
            product="SFx- aggregate",
            electron_energy_eV=attach_energy,
            cross_section_m2=attach_sigma,
            comments=(
                "NIST Table 28 suggested total room-temperature attachment",
                "aggregate negative-ion branching unresolved",
                "linear closure from 15 eV to zero at 20 eV",
            ),
        ),
    )
    derivation = {
        "schema": "petch.nist_evaluated_sf6_replay.v1",
        "table_sha256": {
            "9": NIST_SF6_TABLE9_CSV_SHA256,
            "14": NIST_SF6_TABLE14_CSV_SHA256,
            "15": NIST_SF6_TABLE15_CSV_SHA256,
            "17": NIST_SF6_TABLE17_CSV_SHA256,
            "20": NIST_SF6_TABLE20_CSV_SHA256,
            "28": NIST_SF6_TABLE28_CSV_SHA256,
        },
        "vibrational_energy_loss_eV": loss,
        "maximum_energy_eV": maximum,
        "table9_deconvolution_upper_eV": 2.5,
        "table14_momentum_lower_eV": 2.75,
        "low_energy_momentum_closure": "constant_first_deconvolved_value_to_zero",
        "vibration_tail_closure": "linear_to_zero_at_15eV",
        "attachment_tail_closure": "linear_to_zero_at_20eV",
        "neutral_dissociation_threshold_eV": 9.6,
        "ionization_threshold_eV": 15.5,
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    deck = ElectronCollisionDeck(
        processes=processes,
        payload_sha256=digest,
        source_database="NIST-evaluated SF6 aggregate electron interactions",
        retrieved_at="2026-08-18",
        source_reference=(
            "Christophorou--Olthoff 2000 Tables 9,14,15,17,20,28; "
            f"doi:10.1063/1.1288407; derivation={json.dumps(derivation, sort_keys=True)}"
        ),
    )
    return NISTEvaluatedSF6Replay(
        derived_deck=deck,
        source_curves=(
            total, momentum, vibration, ionization, dissociation, attachment,
        ),
        vibrational_energy_loss_eV=loss,
        maximum_energy_eV=maximum,
    )


_SF6_POSITIVE_ION_APPEARANCE_EV = (
    ("SF5+", 15.5),
    ("SF4+", 18.0),
    ("SF3+", 19.0),
    ("SF2+", 26.0),
    ("SF+", 31.0),
    ("S+", 37.0),
    ("F+", 35.8),
    ("SF4++", 33.0),
    ("SF2++", 40.0),
)


def _finite_support_interpolation(
    grid_eV: np.ndarray,
    curve: NISTSF6ProductCrossSectionCurve,
    *,
    constant_below_support: bool = False,
) -> np.ndarray:
    left = curve.cross_section_m2[0] if constant_below_support else 0.0
    values = np.interp(
        grid_eV,
        curve.electron_energy_eV,
        curve.cross_section_m2,
        left=left,
        right=0.0,
    )
    if not constant_below_support:
        values[grid_eV < curve.electron_energy_eV[0]] = 0.0
    values[grid_eV > curve.electron_energy_eV[-1]] = 0.0
    return values


def derive_nist_product_resolved_sf6_replay(
    *,
    vibrational_energy_loss_eV: float = 0.095,
    maximum_energy_eV: float = 200.0,
) -> NISTProductResolvedSF6Replay:
    """Resolve primary charged products without changing evaluated totals.

    The attachment split is a reconciliation of the product tables with the
    independently rounded total table.  SF6- and the five minor anions retain
    their tabulated shapes; when their rounded sum exceeds the total they are
    rescaled together by less than one percent, and SF5- receives the exact
    non-negative residual.  The independently tabulated SF5- curve is retained
    as a source check.

    For positive ions, only the 100 eV partial values are tabulated.  Each
    branch is given a linear excess-energy onset from the earliest reported
    formation threshold in Table 3, normalized at 100 eV, and all branches are
    then renormalized to the recommended Table 17 total at every energy.  This
    closure is deterministic and conservative, but it is not promoted to a
    measured partial-ionization curve.
    """

    aggregate = derive_nist_evaluated_sf6_replay(
        vibrational_energy_loss_eV=vibrational_energy_loss_eV,
        maximum_energy_eV=maximum_energy_eV,
    )
    anchors = load_nist_2000_sf6_partial_ionization_anchors()
    attachment_curves = load_nist_2000_sf6_partial_attachment_curves()
    anchor_by_product = {curve.product: curve for curve in anchors}
    threshold_by_product = dict(_SF6_POSITIVE_ION_APPEARANCE_EV)
    if set(anchor_by_product) != set(threshold_by_product):
        raise RuntimeError("SF6 partial-ion anchors and thresholds differ")

    aggregate_ionization = next(
        process for process in aggregate.derived_deck.processes
        if process.kind == "IONIZATION"
    )
    ion_grid = np.unique(np.concatenate((
        aggregate_ionization.electron_energy_eV,
        [100.0],
        np.asarray(tuple(threshold_by_product.values())),
    )))
    ion_total = np.interp(
        ion_grid,
        aggregate_ionization.electron_energy_eV,
        aggregate_ionization.cross_section_m2,
    )
    raw_ion = []
    for product, threshold in _SF6_POSITIVE_ION_APPEARANCE_EV:
        anchor = anchor_by_product[product].cross_section_m2[0]
        onset = np.zeros(ion_grid.shape)
        above = ion_grid > threshold
        onset[above] = (
            (1.0 - threshold / ion_grid[above])
            / (1.0 - threshold / 100.0)
        )
        raw_ion.append(anchor * onset)
    raw_ion = np.asarray(raw_ion)
    raw_sum = np.sum(raw_ion, axis=0)
    ion_scale = np.divide(
        ion_total,
        raw_sum,
        out=np.zeros_like(ion_total),
        where=raw_sum > 0.0,
    )
    ion_branches = raw_ion * ion_scale[None, :]
    if not np.allclose(
        np.sum(ion_branches, axis=0), ion_total, rtol=2e-15, atol=1e-35
    ):
        raise RuntimeError("SF6 positive-ion branching failed total closure")

    aggregate_attachment = next(
        process for process in aggregate.derived_deck.processes
        if process.kind == "ATTACHMENT"
    )
    attachment_grid = np.unique(np.concatenate((
        aggregate_attachment.electron_energy_eV,
        *(curve.electron_energy_eV for curve in attachment_curves),
    )))
    attachment_total = np.interp(
        attachment_grid,
        aggregate_attachment.electron_energy_eV,
        aggregate_attachment.cross_section_m2,
    )
    attachment_by_product = {
        curve.product: curve for curve in attachment_curves
    }
    direct_products = ("SF6-", "SF4-", "SF3-", "SF2-", "F2-", "F-")
    direct_attachment = np.asarray([
        _finite_support_interpolation(
            attachment_grid,
            attachment_by_product[product],
            constant_below_support=(product == "SF6-"),
        )
        for product in direct_products
    ])
    direct_sum = np.sum(direct_attachment, axis=0)
    rounding_scale = np.minimum(
        1.0,
        np.divide(
            attachment_total,
            direct_sum,
            out=np.ones_like(attachment_total),
            where=direct_sum > 0.0,
        ),
    )
    direct_attachment *= rounding_scale[None, :]
    sf5_attachment = np.maximum(
        attachment_total - np.sum(direct_attachment, axis=0), 0.0)
    attachment_products = ("SF5-", *direct_products)
    attachment_branches = np.vstack((sf5_attachment, direct_attachment))
    if not np.allclose(
        np.sum(attachment_branches, axis=0),
        attachment_total,
        rtol=2e-15,
        atol=1e-35,
    ):
        raise RuntimeError("SF6 attachment branching failed total closure")
    sf5_source = attachment_by_product["SF5-"]
    sf5_reconciled = np.interp(
        sf5_source.electron_energy_eV,
        attachment_grid,
        sf5_attachment,
    )
    sf5_peak_residual = float(
        np.max(np.abs(sf5_reconciled - sf5_source.cross_section_m2))
        / np.max(sf5_source.cross_section_m2)
    )
    maximum_rounding_rescale = float(np.max(1.0 - rounding_scale))

    unchanged = tuple(
        process for process in aggregate.derived_deck.processes
        if process.kind not in {"IONIZATION", "ATTACHMENT"}
    )
    ion_processes = tuple(
        ElectronCollisionProcess(
            kind="IONIZATION",
            target="SF6",
            product=product,
            electron_energy_eV=tuple(ion_grid),
            cross_section_m2=tuple(ion_branches[index]),
            energy_loss_eV=threshold,
            electron_number_change=2 if product.endswith("++") else 1,
            comments=(
                "NIST Table 17 recommended total ionization preserved exactly",
                "NIST Table 16 measured preferred-analysis branch anchor at 100 eV",
                f"Table 3 earliest reported formation threshold={threshold:g} eV",
                "declared normalized linear excess-energy onset between threshold and 100 eV",
                "branch curves are closure, not direct partial-curve measurements",
            ),
        )
        for index, (product, threshold) in enumerate(
            _SF6_POSITIVE_ION_APPEARANCE_EV
        )
    )
    attachment_processes = tuple(
        ElectronCollisionProcess(
            kind="ATTACHMENT",
            target="SF6",
            product=product,
            electron_energy_eV=tuple(attachment_grid),
            cross_section_m2=tuple(attachment_branches[index]),
            comments=(
                "NIST Table 28 suggested total attachment preserved exactly",
                "NIST Tables 24/25/27 product evidence with explicit rounding reconciliation",
                (
                    "SF5- assigned exact non-negative residual after other products"
                    if product == "SF5-"
                    else "tabulated product shape retained subject to total-table rounding rescale"
                ),
                "Table 27 blank cells treated as unreported, not measured zero",
            ),
        )
        for index, product in enumerate(attachment_products)
    )
    derivation = {
        "schema": "petch.nist_product_resolved_sf6_replay.v1",
        "aggregate_payload_sha256": aggregate.derived_deck.payload_sha256,
        "product_table_sha256": {
            "16": NIST_SF6_TABLE16_CSV_SHA256,
            "24": NIST_SF6_TABLE24_CSV_SHA256,
            "25": NIST_SF6_TABLE25_CSV_SHA256,
            "27": NIST_SF6_TABLE27_CSV_SHA256,
        },
        "positive_ion_appearance_eV": dict(
            _SF6_POSITIVE_ION_APPEARANCE_EV),
        "positive_ion_shape": (
            "linear_excess_energy_onset_normalized_at_100eV_then_"
            "renormalized_to_table17"
        ),
        "attachment_shape": (
            "tables24_27_direct_shapes_rounding_scaled_if_needed;_"
            "sf5_exact_nonnegative_table28_residual"
        ),
        "double_ionization_electron_number_change": 2,
        "maximum_energy_eV": float(maximum_energy_eV),
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    deck = ElectronCollisionDeck(
        processes=(*unchanged, *ion_processes, *attachment_processes),
        payload_sha256=digest,
        source_database="NIST-evaluated SF6 product-resolved electron interactions",
        retrieved_at="2026-08-18",
        source_reference=(
            "Christophorou--Olthoff 2000 Tables 3,16,17,24,25,27,28; "
            f"doi:10.1063/1.1288407; derivation={json.dumps(derivation, sort_keys=True)}"
        ),
    )
    anchor_sum = float(sum(
        curve.cross_section_m2[0] for curve in anchors))
    total_at_100 = float(np.interp(
        100.0,
        aggregate_ionization.electron_energy_eV,
        aggregate_ionization.cross_section_m2,
    ))
    return NISTProductResolvedSF6Replay(
        aggregate_replay=aggregate,
        derived_deck=deck,
        ionization_source_anchors=anchors,
        attachment_source_curves=attachment_curves,
        ionization_threshold_eV=_SF6_POSITIVE_ION_APPEARANCE_EV,
        ionization_anchor_sum_m2=anchor_sum,
        evaluated_total_ionization_at_100eV_m2=total_at_100,
        maximum_attachment_rounding_rescale_fraction=maximum_rounding_rescale,
        maximum_sf5_source_peak_normalized_residual=sf5_peak_residual,
    )
