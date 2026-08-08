"""Explicit legacy-SIGLO chlorine source-replay adapter.

The raw LXCat/SIGLO bytes remain user supplied and outside the package.  This
adapter recognizes one hash-locked 2013 deck, reproduces BOLOS/BOLSIG edge
padding as a declared numerical convention, declares double-ionization
multiplicity, and maps every inelastic row into the six-species chlorine
ledger.  It is a source-replay/sensitivity provider, not current collision
evidence or a validated reactor mechanism.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path

from .chlorine import lee_lieberman_chlorine_species
from .electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from .electron_collision_deck import (
    ElectronCollisionDeck,
    load_bolsig_lxcat_file,
)


LEGACY_SIGLO_CL2_2013_SHA256 = (
    "b8b1ff807d0586795dcabf4511a2a44b8b04828626b45e0fd95a319817c01063"
)
_EXPECTED_ROWS = (
    ("ATTACHMENT", "Cl^- + Cl", 0.0),
    ("ELASTIC", None, None),
    ("EXCITATION", "Cl2(v)", 0.069),
    ("EXCITATION", "Cl2(v)", 0.139),
    ("EXCITATION", "Cl2(3PI_u)", 3.36),
    ("EXCITATION", "Cl2(1PI_u)", 4.3),
    ("EXCITATION", "Cl2(3PI_g)", 6.38),
    ("EXCITATION", "Cl2(1PI_g)", 7.01),
    ("EXCITATION", "Cl2(3SIG_u)", 7.02),
    ("EXCITATION", "Cl2(1PI_ub)", 10.54),
    ("EXCITATION", "Cl2(1SIG_ub)", 10.7),
    ("EXCITATION", "Cl^- + Cl^+", 11.0),
    ("IONIZATION", "Cl2^+", 11.49),
    ("IONIZATION", "Cl^+ + Cl", 11.49),
    ("IONIZATION", "Cl2^++", 35.5),
    ("IONIZATION", "Cl^++ + Cl", 43.5),
)


@dataclass(frozen=True)
class LegacySigloChlorineReplay:
    raw_payload_sha256: str
    derived_deck: ElectronCollisionDeck
    collision_chemistry: ElectronCollisionChemistry
    maximum_energy_eV: float
    missing_reactor_channels: tuple[str, ...] = (
        "atomic_chlorine_ionization",
        "electron_detachment_from_Clminus",
        "tracked_vibrational_and_electronic_state_kinetics",
    )
    declared_sensitivity_closures: tuple[str, ...] = (
        "constant edge-value cross-section padding to solver maximum energy",
        "elastic low-energy edge padded with first printed value",
        "Cl2++ and Cl++ channels collapse to two tracked Cl+ products",
    )
    supports_direct_swarm_grade: bool = False
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


def _verify_rows(deck: ElectronCollisionDeck) -> None:
    actual = tuple(
        (process.kind, process.product, process.energy_loss_eV)
        for process in deck.processes
    )
    if actual != _EXPECTED_ROWS:
        raise RuntimeError(
            "legacy SIGLO chlorine process topology/signatures changed"
        )


def derive_legacy_siglo_cl2_replay(
    raw_deck: ElectronCollisionDeck,
    *,
    maximum_energy_eV: float = 200.0,
) -> LegacySigloChlorineReplay:
    """Build the explicit source-replay deck from an already hash-gated deck."""
    maximum = float(maximum_energy_eV)
    if raw_deck.payload_sha256 != LEGACY_SIGLO_CL2_2013_SHA256:
        raise RuntimeError("legacy SIGLO chlorine raw hash mismatch")
    if maximum <= 100.0:
        raise ValueError("legacy SIGLO replay maximum energy must exceed 100 eV")
    _verify_rows(raw_deck)
    processes = []
    for index, process in enumerate(raw_deck.processes):
        energy = list(process.electron_energy_eV)
        cross_section = list(process.cross_section_m2)
        if process.kind in {"ELASTIC", "MOMENTUM", "EFFECTIVE"} and energy[0] > 0.0:
            energy.insert(0, 0.0)
            cross_section.insert(0, cross_section[0])
        if energy[-1] < maximum:
            energy.append(maximum)
            cross_section.append(cross_section[-1])
        processes.append(replace(
            process,
            electron_energy_eV=tuple(energy),
            cross_section_m2=tuple(cross_section),
            electron_number_change=(
                2 if index in {14, 15} else process.electron_number_change),
        ))
    derivation = {
        "schema": "petch.legacy_siglo_cl2_replay.v1",
        "raw_payload_sha256": raw_deck.payload_sha256,
        "maximum_energy_eV": maximum,
        "low_edge": "constant_first_value_for_momentum",
        "high_edge": "constant_last_value_for_all_processes",
        "double_ionization_rows": [14, 15],
    }
    derived_sha = hashlib.sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    derived = ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=derived_sha,
        source_database="legacy SIGLO Cl2 2013 explicit source replay",
        retrieved_at=raw_deck.retrieved_at,
        source_reference=(
            f"{raw_deck.source_reference}; raw_sha256={raw_deck.payload_sha256}; "
            f"derivation={json.dumps(derivation, sort_keys=True)}"
        ),
        packaged_or_redistributed=False,
    )
    common_source = (
        "legacy SIGLO Cl2 2013 source-replay mapping; no reactor fit"
    )
    mappings = []

    def mapping(index, name, reactants, products, evidence="published_compilation"):
        mappings.append(ElectronCollisionHeavyMapping(
            process_index=index,
            reaction_name=name,
            heavy_reactants=reactants,
            heavy_products=products,
            source=common_source,
            evidence_kind=evidence,
        ))

    mapping(0, "dissociative_attachment", {"Cl2": 1}, {"Cl-": 1, "Cl": 1})
    mapping(2, "vibrational_excitation_069", {"Cl2": 1}, {"Cl2": 1})
    mapping(3, "vibrational_excitation_139", {"Cl2": 1}, {"Cl2": 1})
    for index, state in zip(range(4, 9), (
        "3PI_u", "1PI_u", "3PI_g", "1PI_g", "3SIG_u",
    )):
        mapping(
            index,
            f"dissociative_excitation_{state}",
            {"Cl2": 1},
            {"Cl": 2},
        )
    mapping(9, "rydberg_excitation_1PI_ub", {"Cl2": 1}, {"Cl2": 1})
    mapping(10, "rydberg_excitation_1SIG_ub", {"Cl2": 1}, {"Cl2": 1})
    mapping(11, "ion_pair_formation", {"Cl2": 1}, {"Cl-": 1, "Cl+": 1})
    mapping(12, "molecular_ionization", {"Cl2": 1}, {"Cl2+": 1})
    mapping(13, "dissociative_ionization", {"Cl2": 1}, {"Cl+": 1, "Cl": 1})
    mapping(
        14,
        "double_molecular_ionization_collapsed_to_two_Clplus",
        {"Cl2": 1},
        {"Cl+": 2},
        evidence="sensitivity",
    )
    mapping(
        15,
        "double_atomic_ionization_collapsed_to_two_Clplus",
        {"Cl2": 1},
        {"Cl+": 2},
        evidence="sensitivity",
    )
    chemistry = ElectronCollisionChemistry(
        derived, lee_lieberman_chlorine_species(), tuple(mappings))
    return LegacySigloChlorineReplay(
        raw_payload_sha256=raw_deck.payload_sha256,
        derived_deck=derived,
        collision_chemistry=chemistry,
        maximum_energy_eV=maximum,
    )


def load_legacy_siglo_cl2_replay(
    path: str | Path,
    *,
    maximum_energy_eV: float = 200.0,
) -> LegacySigloChlorineReplay:
    raw = load_bolsig_lxcat_file(
        path,
        source_database="SIGLO LXCat legacy chlorine",
        retrieved_at="2013-06-04",
        source_reference=str(Path(path)),
        target="Cl2",
        expected_sha256=LEGACY_SIGLO_CL2_2013_SHA256,
    )
    return derive_legacy_siglo_cl2_replay(
        raw, maximum_energy_eV=maximum_energy_eV)
