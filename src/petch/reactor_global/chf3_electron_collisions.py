"""Hash-locked CHF3 electron-collision source replay and swarm evidence.

Kushner and Zhang constructed a *working* cross-section set by reconciling
measurements, calculations, and swarm data.  This module preserves that
evidence class: reproducing the NIST drift curve is a source-replay check,
not independent validation of the individual chemical branches.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
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


CHF3_MASS_AMU = 70.013849489
KUSHNER_ZHANG_CHF3_CSV_SHA256 = (
    "61d6b8ba5582174fa8300c62c76f843db91d71a41ad40c5b8770d93951ab1735"
)
NIST_CHF3_TABLE6_CSV_SHA256 = (
    "d5e0de969eedd64f38e8b15b7f5bd5c90f6e0ff6ccf22e77aebefd9c3a742901"
)
NIST_CHF3_TABLE4_CSV_SHA256 = (
    "845b090424f9df9ca8d5612783fe9cbf7a275a1bc26f239940cab266531a9da1"
)
NIST_CHF3_TABLE5_CSV_SHA256 = (
    "fd6bc530ab6794f56cb733862a71b7cd40d4757b6b9d604a63e0841e498fb865"
)
_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_KUSHNER_ZHANG_CHF3_CSV = (
    _ROOT / "data" / "experimental" / "kushner_zhang_2000_chf3"
    / "cross_sections.csv"
)
DEFAULT_NIST_CHF3_TABLE6_CSV = (
    _ROOT / "data" / "experimental"
    / "christophorou_olthoff_1999_chf3"
    / "table6_drift_velocity.csv"
)
DEFAULT_NIST_CHF3_TABLE4_CSV = (
    _ROOT / "data" / "experimental"
    / "christophorou_olthoff_1999_chf3"
    / "table4_total_scattering.csv"
)
DEFAULT_NIST_CHF3_TABLE5_CSV = (
    _ROOT / "data" / "experimental"
    / "christophorou_olthoff_1999_chf3"
    / "table5_momentum_transfer.csv"
)
_EXPECTED_PROCESS_SIGNATURES = (
    ("MOM", "MOMENTUM", "CHF3", None, 0, 150),
    ("VIB14", "EXCITATION", "CHF3(v1,v4)", 0.37, 0, 100),
    ("VIB25", "EXCITATION", "CHF3(v2,v5)", 0.18, 0, 100),
    ("VIB36", "EXCITATION", "CHF3(v3,v6)", 0.13, 0, 100),
    ("NEU1", "EXCITATION", "CF3 + H", 11.0, 0, 200),
    ("NEU2", "EXCITATION", "CHF2 + F", 13.0, 0, 200),
    ("NEU3", "EXCITATION", "CF2 + H + F", 23.6, 0, 200),
    ("NEU4", "EXCITATION", "CHF + F + F", 35.0, 0, 200),
    ("NEU5", "EXCITATION", "CF + H + F + F", 19.5, 0, 200),
    ("NEU6", "EXCITATION", "CF + H + F2", 19.5, 0, 200),
    ("NEU1_ADD_ON", "EXCITATION", "CF3 + H", 11.0, 0, 200),
    ("ION1", "IONIZATION", "CF3+ + H", 15.2, 1, 200),
    ("ION2", "IONIZATION", "CHF2+ + F", 16.8, 1, 200),
    ("ION3", "IONIZATION", "CF2+ + HF", 17.6, 1, 200),
    ("ION4", "IONIZATION", "CHF+ + F + F", 19.8, 1, 200),
    ("ION5", "IONIZATION", "CF+ + HF + F", 20.9, 1, 200),
    ("ION6", "IONIZATION", "CH+ + F2 + F", 33.5, 1, 200),
    ("ION7", "IONIZATION", "F+ + CHF2", 37.0, 1, 200),
    ("ATT1", "ATTACHMENT", "CHF2 + F-", 0.0, -1, 183),
    # Ion-pair production, not net electron attachment.
    ("ATT2", "EXCITATION", "CHF2+ + F-", 11.5, 0, 183),
)


@dataclass(frozen=True)
class KushnerZhangCHF3Replay:
    """Raw author table plus one explicit numerical support closure."""

    raw_deck: ElectronCollisionDeck
    derived_deck: ElectronCollisionDeck
    process_labels: tuple[str, ...]
    maximum_energy_eV: float
    low_energy_momentum_rule: str = "constant_first_source_value_to_zero"
    high_energy_tail_rule: str = "constant_last_source_value"
    evidence_class: str = "swarm_regressed_working_set"
    supports_source_replay: bool = True
    supports_independent_branch_validation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            self.raw_deck.payload_sha256 != KUSHNER_ZHANG_CHF3_CSV_SHA256
            or len(self.process_labels) != len(self.raw_deck.processes)
            or len(self.raw_deck.processes) != len(self.derived_deck.processes)
            or not math.isfinite(self.maximum_energy_eV)
            or self.maximum_energy_eV <= 0.0
            or self.low_energy_momentum_rule
            != "constant_first_source_value_to_zero"
            or self.high_energy_tail_rule != "constant_last_source_value"
            or self.evidence_class != "swarm_regressed_working_set"
            or not self.supports_source_replay
            or self.supports_independent_branch_validation
            or self.supports_reactor_state_prediction
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid Kushner--Zhang CHF3 replay")


@dataclass(frozen=True)
class NISTCHF3DriftCurve:
    reduced_electric_field_Td: np.ndarray
    drift_velocity_m_s: np.ndarray
    gas_temperature_K: float = 298.0
    evidence_class: str = "recommended_measured_swarm_fit"
    supports_independent_grade_of_working_set: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        field = np.asarray(self.reduced_electric_field_Td, dtype=float).copy()
        velocity = np.asarray(self.drift_velocity_m_s, dtype=float).copy()
        if (
            field.ndim != 1
            or field.size != 33
            or velocity.shape != field.shape
            or np.any(~np.isfinite(field))
            or np.any(~np.isfinite(velocity))
            or np.any(field <= 0.0)
            or np.any(velocity <= 0.0)
            or np.any(np.diff(field) <= 0.0)
            or self.gas_temperature_K != 298.0
            or self.evidence_class != "recommended_measured_swarm_fit"
            or self.supports_independent_grade_of_working_set
            or self.supports_reactor_state_prediction
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid NIST CHF3 drift curve")
        field.setflags(write=False)
        velocity.setflags(write=False)
        object.__setattr__(self, "reduced_electric_field_Td", field)
        object.__setattr__(self, "drift_velocity_m_s", velocity)


@dataclass(frozen=True)
class NISTCHF3CrossSectionCurve:
    electron_energy_eV: np.ndarray
    cross_section_m2: np.ndarray
    evidence_class: str

    def __post_init__(self):
        energy = np.asarray(self.electron_energy_eV, dtype=float).copy()
        cross_section = np.asarray(self.cross_section_m2, dtype=float).copy()
        if (
            energy.ndim != 1
            or energy.size < 5
            or cross_section.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(cross_section))
            or np.any(energy <= 0.0)
            or np.any(cross_section <= 0.0)
            or np.any(np.diff(energy) <= 0.0)
            or self.evidence_class not in {
                "recommended_total_scattering",
                "suggested_calculated_momentum_transfer",
            }
        ):
            raise ValueError("invalid NIST CHF3 cross-section curve")
        energy.setflags(write=False)
        cross_section.setflags(write=False)
        object.__setattr__(self, "electron_energy_eV", energy)
        object.__setattr__(self, "cross_section_m2", cross_section)


@dataclass(frozen=True)
class NISTEvaluatedCHF3Replay:
    """Evaluated elastic backbone plus fixed working-set chemistry.

    Table 4 is total scattering, so the fixed Kushner--Zhang inelastic rows
    are subtracted below 9 eV before the solver re-adds them to its momentum
    sum.  Table 5 is already an elastic momentum-transfer calculation.  This
    deconvolution is explicit and deterministic; it is not a direct elastic
    measurement below 9 eV.
    """

    working_set: KushnerZhangCHF3Replay
    derived_deck: ElectronCollisionDeck
    total_scattering: NISTCHF3CrossSectionCurve
    momentum_transfer: NISTCHF3CrossSectionCurve
    high_energy_closure: str
    evidence_class: str = "evaluated_transport_hybrid_fixed_chemistry"
    supports_direct_transport_constraints: bool = True
    supports_independent_chemical_branch_validation: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            self.high_energy_closure not in {
                "constant_join_ratio", "linear_return_to_working_set_at_120eV",
            }
            or self.evidence_class
            != "evaluated_transport_hybrid_fixed_chemistry"
            or not self.supports_direct_transport_constraints
            or self.supports_independent_chemical_branch_validation
            or self.supports_reactor_state_prediction
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid NIST-evaluated CHF3 replay")


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _read_processes(path: Path):
    if _digest(path) != KUSHNER_ZHANG_CHF3_CSV_SHA256:
        raise RuntimeError("Kushner--Zhang CHF3 CSV checksum changed")
    groups = {}
    order = []
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            label = row["process_label"]
            if label not in groups:
                groups[label] = []
                order.append(label)
            groups[label].append(row)
    if tuple(order) != tuple(item[0] for item in _EXPECTED_PROCESS_SIGNATURES):
        raise RuntimeError("Kushner--Zhang process order changed")
    processes = []
    for signature in _EXPECTED_PROCESS_SIGNATURES:
        label, kind, product, loss, delta_e, count = signature
        rows = groups[label]
        row0 = rows[0]
        parsed_loss = (
            None if row0["energy_loss_eV"] == ""
            else float(row0["energy_loss_eV"])
        )
        loss_changed = (
            any(float(row["energy_loss_eV"]) != loss for row in rows)
            if loss is not None
            else any(row["energy_loss_eV"] != "" for row in rows)
        )
        if (
            len(rows) != count
            or any(row["kind"] != kind for row in rows)
            or any(row["target"] != "CHF3" for row in rows)
            or any(row["product"] != product for row in rows)
            or loss_changed
        ):
            raise RuntimeError(f"Kushner--Zhang signature changed for {label}")
        if parsed_loss != loss or any(
            int(row["electron_number_change"]) != delta_e for row in rows
        ):
            raise RuntimeError(f"Kushner--Zhang balance changed for {label}")
        energy = tuple(float(row["energy_eV"]) for row in rows)
        sigma = tuple(float(row["cross_section_m2"]) for row in rows)
        processes.append(ElectronCollisionProcess(
            kind=kind,
            target="CHF3",
            product=product,
            electron_energy_eV=energy,
            cross_section_m2=sigma,
            mass_ratio=(
                ELECTRON_MASS_AMU / CHF3_MASS_AMU
                if kind == "MOMENTUM" else None
            ),
            energy_loss_eV=loss,
            electron_number_change=delta_e,
            comments=(label, "Kushner--Zhang swarm-regressed working set"),
        ))
    return tuple(order), tuple(processes)


def load_kushner_zhang_2000_chf3_replay(
    path: str | Path = DEFAULT_KUSHNER_ZHANG_CHF3_CSV,
    *,
    maximum_energy_eV: float = 120.0,
) -> KushnerZhangCHF3Replay:
    """Load the exact working set and add a declared high-energy tail."""

    path = Path(path)
    maximum = float(maximum_energy_eV)
    if not math.isfinite(maximum) or not 40.0 <= maximum <= 120.0:
        raise ValueError("CHF3 replay maximum energy must be in [40, 120] eV")
    labels, processes = _read_processes(path)
    raw = ElectronCollisionDeck(
        processes=processes,
        payload_sha256=KUSHNER_ZHANG_CHF3_CSV_SHA256,
        source_database="Kushner--Zhang 2000 author tabulation",
        retrieved_at="2026-08-18",
        source_reference=(
            "https://cpseg.eecs.umich.edu/pub/data/chf3_xsec.xls; "
            "doi:10.1063/1.1289187"
        ),
    )
    derived_processes = []
    extended_labels = []
    for label, process in zip(labels, processes):
        energy = list(process.electron_energy_eV)
        sigma = list(process.cross_section_m2)
        if process.kind == "MOMENTUM" and energy[0] > 0.0:
            energy.insert(0, 0.0)
            sigma.insert(0, sigma[0])
        if energy[-1] < maximum:
            energy.append(maximum)
            sigma.append(sigma[-1])
            extended_labels.append(label)
        derived_processes.append(replace(
            process,
            electron_energy_eV=tuple(energy),
            cross_section_m2=tuple(sigma),
        ))
    derivation = {
        "schema": "petch.kushner_zhang_chf3_replay.v1",
        "raw_payload_sha256": raw.payload_sha256,
        "maximum_energy_eV": maximum,
        "low_edge_rule": "constant_first_source_value_to_zero_for_momentum",
        "tail_rule": "constant_last_source_value",
        "extended_process_labels": extended_labels,
    }
    derived_sha = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    derived = ElectronCollisionDeck(
        processes=tuple(derived_processes),
        payload_sha256=derived_sha,
        source_database="Kushner--Zhang CHF3 explicit source replay",
        retrieved_at=raw.retrieved_at,
        source_reference=(
            f"{raw.source_reference}; derivation="
            f"{json.dumps(derivation, sort_keys=True)}"
        ),
    )
    return KushnerZhangCHF3Replay(
        raw_deck=raw,
        derived_deck=derived,
        process_labels=labels,
        maximum_energy_eV=maximum,
    )


def load_nist_1999_chf3_drift_curve(
    path: str | Path = DEFAULT_NIST_CHF3_TABLE6_CSV,
) -> NISTCHF3DriftCurve:
    path = Path(path)
    if _digest(path) != NIST_CHF3_TABLE6_CSV_SHA256:
        raise RuntimeError("NIST CHF3 Table 6 CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(
        row["gas_temperature_K"] != "298"
        or row["evidence_class"] != "recommended_measured_swarm_fit"
        for row in rows
    ):
        raise RuntimeError("NIST CHF3 Table 6 metadata changed")
    return NISTCHF3DriftCurve(
        reduced_electric_field_Td=np.asarray([
            float(row["reduced_electric_field_Td"]) for row in rows
        ]),
        drift_velocity_m_s=np.asarray([
            float(row["drift_velocity_m_s"]) for row in rows
        ]),
    )


def _load_nist_cross_section(
    path: str | Path,
    *,
    expected_sha256: str,
    evidence_class: str,
) -> NISTCHF3CrossSectionCurve:
    path = Path(path)
    if _digest(path) != expected_sha256:
        raise RuntimeError("NIST CHF3 cross-section CSV checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if any(row["evidence_class"] != evidence_class for row in rows):
        raise RuntimeError("NIST CHF3 cross-section metadata changed")
    return NISTCHF3CrossSectionCurve(
        electron_energy_eV=np.asarray([
            float(row["electron_energy_eV"]) for row in rows
        ]),
        cross_section_m2=np.asarray([
            float(row["cross_section_m2"]) for row in rows
        ]),
        evidence_class=evidence_class,
    )


def load_nist_1999_chf3_total_scattering(
    path: str | Path = DEFAULT_NIST_CHF3_TABLE4_CSV,
) -> NISTCHF3CrossSectionCurve:
    return _load_nist_cross_section(
        path,
        expected_sha256=NIST_CHF3_TABLE4_CSV_SHA256,
        evidence_class="recommended_total_scattering",
    )


def load_nist_1999_chf3_momentum_transfer(
    path: str | Path = DEFAULT_NIST_CHF3_TABLE5_CSV,
) -> NISTCHF3CrossSectionCurve:
    return _load_nist_cross_section(
        path,
        expected_sha256=NIST_CHF3_TABLE5_CSV_SHA256,
        evidence_class="suggested_calculated_momentum_transfer",
    )


def derive_nist_evaluated_chf3_replay(
    working_set: KushnerZhangCHF3Replay | None = None,
    *,
    high_energy_closure: str = "constant_join_ratio",
) -> NISTEvaluatedCHF3Replay:
    """Replace only the uncertain CHF3 momentum backbone with NIST data."""

    if working_set is None:
        working_set = load_kushner_zhang_2000_chf3_replay()
    if not isinstance(working_set, KushnerZhangCHF3Replay):
        raise TypeError("a Kushner--Zhang CHF3 working set is required")
    if high_energy_closure not in {
        "constant_join_ratio", "linear_return_to_working_set_at_120eV",
    }:
        raise ValueError("unsupported CHF3 high-energy closure")
    total = load_nist_1999_chf3_total_scattering()
    momentum = load_nist_1999_chf3_momentum_transfer()
    source_processes = list(working_set.derived_deck.processes)
    source_momentum = source_processes[0]
    if source_momentum.kind != "MOMENTUM":
        raise RuntimeError("working-set momentum row moved")

    # The solver includes every inelastic cross section in sigma_m. Remove
    # those same rows from the total-scattering proxy exactly once.
    inelastic = source_processes[1:]
    inelastic_on_total_grid = np.zeros_like(total.electron_energy_eV)
    for process in inelastic:
        inelastic_on_total_grid += np.interp(
            total.electron_energy_eV,
            process.electron_energy_eV,
            process.cross_section_m2,
        )
    use_total = total.electron_energy_eV <= 9.0
    total_energy = total.electron_energy_eV[use_total]
    elastic_proxy = (
        total.cross_section_m2[use_total]
        - inelastic_on_total_grid[use_total]
    )
    if np.any(elastic_proxy <= 0.0):
        raise RuntimeError("NIST total-scattering deconvolution became nonpositive")

    source_energy = np.asarray(source_momentum.electron_energy_eV)
    source_sigma = np.asarray(source_momentum.cross_section_m2)
    join_energy = float(momentum.electron_energy_eV[-1])
    join_ratio = float(
        momentum.cross_section_m2[-1]
        / np.interp(join_energy, source_energy, source_sigma)
    )
    high_mask = (
        (source_energy > join_energy)
        & (source_energy <= working_set.maximum_energy_eV)
    )
    high_energy = np.unique(np.concatenate((
        source_energy[high_mask],
        np.array([working_set.maximum_energy_eV]),
    )))
    if high_energy_closure == "constant_join_ratio":
        high_scale = np.full(high_energy.shape, join_ratio)
    else:
        high_scale = 1.0 + (join_ratio - 1.0) * (
            (working_set.maximum_energy_eV - high_energy)
            / (working_set.maximum_energy_eV - join_energy)
        )
    energy = np.concatenate((
        np.array([0.0]),
        total_energy,
        momentum.electron_energy_eV,
        high_energy,
    ))
    sigma = np.concatenate((
        np.array([elastic_proxy[0]]),
        elastic_proxy,
        momentum.cross_section_m2,
        np.interp(high_energy, source_energy, source_sigma) * high_scale,
    ))
    evaluated_momentum = replace(
        source_momentum,
        electron_energy_eV=tuple(energy),
        cross_section_m2=tuple(sigma),
        comments=(
            "NIST_TABLE4_5_EVALUATED_MOMENTUM",
            "Table 4 total scattering deconvolved below 9 eV",
            "Table 5 suggested calculated momentum transfer at 10-30 eV",
            f"high_energy_closure={high_energy_closure}",
        ),
    )
    source_processes[0] = evaluated_momentum
    derivation = {
        "schema": "petch.nist_evaluated_chf3_replay.v1",
        "working_set_sha256": working_set.derived_deck.payload_sha256,
        "table4_sha256": NIST_CHF3_TABLE4_CSV_SHA256,
        "table5_sha256": NIST_CHF3_TABLE5_CSV_SHA256,
        "total_scattering_upper_energy_eV": 9.0,
        "momentum_transfer_interval_eV": [10.0, 30.0],
        "high_energy_closure": high_energy_closure,
        "join_ratio": join_ratio,
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    deck = ElectronCollisionDeck(
        processes=tuple(source_processes),
        payload_sha256=digest,
        source_database="NIST-evaluated CHF3 transport hybrid",
        retrieved_at="2026-08-18",
        source_reference=(
            "Christophorou--Olthoff 1999 Tables 4/5; fixed "
            "Kushner--Zhang 2000 inelastic branches; "
            f"derivation={json.dumps(derivation, sort_keys=True)}"
        ),
    )
    return NISTEvaluatedCHF3Replay(
        working_set=working_set,
        derived_deck=deck,
        total_scattering=total,
        momentum_transfer=momentum,
        high_energy_closure=high_energy_closure,
    )
