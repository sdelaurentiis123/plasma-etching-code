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
NIST_SF6_TABLE17_CSV_SHA256 = (
    "9566c88d59c4ab1c39833850dec3ef2535fcb8ccbb89c30612f3a0fc770d7663"
)
NIST_SF6_TABLE20_CSV_SHA256 = (
    "a4a7ed0c6d4227f8fbfaf5e69407d88ce03688c2b12c8a96427f90f5158c42f1"
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
