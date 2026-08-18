"""Hash-locked Lan--Jeon C4F6 electron-collision source replay.

Lan and Jeon constructed an effective C4F6 collision set by fitting electron
swarm drift measurements.  The set closes the local electron Boltzmann
problem, but several excitation/dissociation rows are effective or inherited
from c-C4F8 and the ionization row is total rather than product resolved.
Consequently this module exposes electron transport and aggregate rate
coefficients while failing closed on a C4F6 reactor state, ion composition,
wafer flux, or feature depth.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

from .argon import ELECTRON_MASS_AMU
from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


C4F6_MASS_AMU = 162.03
LAN_JEON_2014_DOI = "10.3938/jkps.64.1320"
MOMENTUM_CSV_SHA256 = (
    "ed3bd6f0281fd498b4b95a5c8a3f34025ccb221d5833897d1483251be1a80cbd"
)
INELASTIC_CSV_SHA256 = (
    "46d65c7ffcd54382b5869d873dddb48ddd16c060466d3ac02ba43478d8fd257a"
)
_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MOMENTUM_CSV = (
    _ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
    / "table1_momentum_transfer.csv"
)
DEFAULT_INELASTIC_CSV = (
    _ROOT / "data" / "experimental" / "lan_jeon_2014_c4f6"
    / "table2_inelastic.csv"
)
_EXPECTED_LABELS = (
    "Qa", "Qv1", "Qv2", "Qex1", "Qex2", "Qex3", "Qex4", "Qdiss", "Qi",
)
_EXPECTED_COUNTS = {
    "Qa": 34,
    "Qv1": 27,
    "Qv2": 35,
    "Qex1": 13,
    "Qex2": 14,
    "Qex3": 12,
    "Qex4": 12,
    "Qdiss": 9,
    "Qi": 10,
}
_PRODUCTS = {
    "Qa": "aggregate C4F6- attachment products",
    "Qv1": "C4F6 effective vibration 1",
    "Qv2": "C4F6 effective vibration 2",
    "Qex1": "C4F6 effective excitation 1",
    "Qex2": "C4F6 effective excitation 2",
    "Qex3": "C4F6 effective excitation 3",
    "Qex4": "C4F6 effective excitation 4",
    "Qdiss": "aggregate neutral C4F6 dissociation products",
    "Qi": "aggregate C4Fx+ ionization products",
}


@dataclass(frozen=True)
class LanJeon2014C4F6Replay:
    """Exact printed tables plus declared finite-support closures."""

    raw_deck: ElectronCollisionDeck
    derived_deck: ElectronCollisionDeck
    process_labels: tuple[str, ...]
    maximum_energy_eV: float
    low_energy_momentum_rule: str = "constant_first_source_value_to_zero"
    high_energy_tail_rule: str = "constant_last_source_value"
    evidence_class: str = "mixed_measured_analog_and_swarm_regressed_working_set"
    supports_source_swarm_replay: bool = True
    supports_independent_branch_validation: bool = False
    supports_resolved_primary_chemistry: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            len(self.raw_deck.processes) != 10
            or len(self.derived_deck.processes) != 10
            or self.process_labels != ("Qm",) + _EXPECTED_LABELS
            or not math.isfinite(self.maximum_energy_eV)
            or not 100.0 <= self.maximum_energy_eV <= 200.0
            or self.low_energy_momentum_rule
            != "constant_first_source_value_to_zero"
            or self.high_energy_tail_rule != "constant_last_source_value"
            or self.evidence_class
            != "mixed_measured_analog_and_swarm_regressed_working_set"
            or not self.supports_source_swarm_replay
            or self.supports_independent_branch_validation
            or self.supports_resolved_primary_chemistry
            or self.supports_reactor_state_prediction
            or self.supports_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("invalid Lan--Jeon C4F6 replay")


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_momentum(path: Path) -> ElectronCollisionProcess:
    if _digest(path) != MOMENTUM_CSV_SHA256:
        raise RuntimeError("Lan--Jeon C4F6 momentum table checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    if (
        len(rows) != 37
        or any(row["evidence_class"] != "swarm_regressed_from_measured_drift" for row in rows)
    ):
        raise RuntimeError("Lan--Jeon C4F6 momentum topology changed")
    return ElectronCollisionProcess(
        kind="MOMENTUM",
        target="C4F6",
        product="C4F6",
        electron_energy_eV=tuple(float(row["electron_energy_eV"]) for row in rows),
        cross_section_m2=tuple(float(row["cross_section_m2"]) for row in rows),
        mass_ratio=ELECTRON_MASS_AMU / C4F6_MASS_AMU,
        comments=(
            "Lan--Jeon 2014 Table 1 Qm",
            "effective momentum transfer regressed against measured drift",
        ),
    )


def _load_inelastic(path: Path) -> tuple[tuple[str, ...], tuple[ElectronCollisionProcess, ...]]:
    if _digest(path) != INELASTIC_CSV_SHA256:
        raise RuntimeError("Lan--Jeon C4F6 inelastic table checksum changed")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = tuple(csv.DictReader(stream))
    grouped: dict[str, list[dict[str, str]]] = {}
    order = []
    for row in rows:
        label = row["process_label"]
        if label not in grouped:
            grouped[label] = []
            order.append(label)
        grouped[label].append(row)
    if tuple(order) != _EXPECTED_LABELS:
        raise RuntimeError("Lan--Jeon C4F6 process order changed")
    processes = []
    for label in _EXPECTED_LABELS:
        group = grouped[label]
        if len(group) != _EXPECTED_COUNTS[label]:
            raise RuntimeError(f"Lan--Jeon C4F6 row count changed for {label}")
        kinds = {row["kind"] for row in group}
        losses = {float(row["energy_loss_eV"]) for row in group}
        evidence = {row["evidence_class"] for row in group}
        if len(kinds) != 1 or len(losses) != 1 or len(evidence) != 1:
            raise RuntimeError(f"Lan--Jeon C4F6 metadata changed for {label}")
        kind = kinds.pop()
        loss = losses.pop()
        process = ElectronCollisionProcess(
            kind=kind,
            target="C4F6",
            product=_PRODUCTS[label],
            electron_energy_eV=tuple(float(row["electron_energy_eV"]) for row in group),
            cross_section_m2=tuple(float(row["cross_section_m2"]) for row in group),
            energy_loss_eV=loss,
            comments=(
                f"Lan--Jeon 2014 Table 2 {label}",
                evidence.pop(),
                "aggregate product branching is unresolved",
            ),
        )
        processes.append(process)
    return tuple(order), tuple(processes)


def load_lan_jeon_2014_c4f6_replay(
    momentum_path: str | Path = DEFAULT_MOMENTUM_CSV,
    inelastic_path: str | Path = DEFAULT_INELASTIC_CSV,
    *,
    maximum_energy_eV: float = 200.0,
) -> LanJeon2014C4F6Replay:
    """Load the exact printed working set and make solver support explicit."""

    maximum = float(maximum_energy_eV)
    if not math.isfinite(maximum) or not 100.0 <= maximum <= 200.0:
        raise ValueError("C4F6 replay maximum energy must be in [100, 200] eV")
    momentum = _load_momentum(Path(momentum_path))
    labels, inelastic = _load_inelastic(Path(inelastic_path))
    processes = (momentum,) + inelastic
    source_digest = sha256(
        (MOMENTUM_CSV_SHA256 + INELASTIC_CSV_SHA256).encode()
    ).hexdigest()
    raw = ElectronCollisionDeck(
        processes=processes,
        payload_sha256=source_digest,
        source_database="Lan--Jeon 2014 C4F6 printed Tables 1-2",
        retrieved_at="2026-08-18",
        source_reference=f"doi:{LAN_JEON_2014_DOI}",
    )
    derived_processes = []
    extended_labels = []
    for label, process in zip(("Qm",) + labels, processes):
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
        "schema": "petch.lan-jeon-2014-c4f6-replay.v1",
        "raw_payload_sha256": raw.payload_sha256,
        "maximum_energy_eV": maximum,
        "low_edge_rule": "constant_first_source_value_to_zero_for_momentum",
        "tail_rule": "constant_last_source_value",
        "extended_process_labels": extended_labels,
        "product_branching": "unresolved_aggregate_only",
    }
    derived_digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    derived = ElectronCollisionDeck(
        processes=tuple(derived_processes),
        payload_sha256=derived_digest,
        source_database="Lan--Jeon C4F6 explicit source replay",
        retrieved_at=raw.retrieved_at,
        source_reference=(
            f"{raw.source_reference}; derivation="
            f"{json.dumps(derivation, sort_keys=True)}"
        ),
    )
    return LanJeon2014C4F6Replay(
        raw_deck=raw,
        derived_deck=derived,
        process_labels=("Qm",) + labels,
        maximum_energy_eV=maximum,
    )
