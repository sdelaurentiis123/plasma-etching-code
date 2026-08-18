"""Rights-safe replay of the Song et al. (2026) O2 collision workbook.

The AIP supplemental workbook is CC BY-NC 4.0.  It is therefore never
packaged by petch: callers supply their own copy and this module verifies its
exact SHA-256 before reading it.  The parser uses only Python's OOXML/ZIP
support and preserves the source workbook byte-for-byte.

The source separates momentum transfer, vibration, five electronic states,
neutral dissociation, total ionization, and dissociative attachment.  Turning
the finite-support tables into a Boltzmann deck requires explicit tail/onset
closures; those choices are included in the derived deck hash.  This is an
electron-kinetic input, not evidence for a unique reactor state or TiO2 depth.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path, PurePosixPath
import re
from zipfile import BadZipFile, ZipFile
import xml.etree.ElementTree as ET

import numpy as np

from .argon import ELECTRON_MASS_AMU
from .electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


O2_MASS_AMU = 31.9988
SONG_2026_O2_WORKBOOK_SHA256 = (
    "6f98ac82e169d25d0a4328b1a3703f733668539adb8141d736209d199013c860"
)
SONG_2026_O2_WORKBOOK_MD5 = "71120ea0c9db8d7c9e5a64af0ea9f9ba"
SONG_2026_O2_DATASET_DOI = "10.60893/figshare.jpr.30850013.v1"
SONG_2026_O2_ARTICLE_DOI = "10.1063/5.0287254"
_M2_PER_SOURCE_UNIT = 1.0e-20  # workbook unit: 10^-16 cm^2
_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_DOC_REL_NS = (
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
)
_PACKAGE_REL_NS = (
    "http://schemas.openxmlformats.org/package/2006/relationships"
)
_CELL_REF = re.compile(r"([A-Z]+)([1-9][0-9]*)")
_EXPECTED_DIMENSIONS = {
    "total scattering cross section": "A1:C106",
    "very low total scattering cross": "A1:C58",
    "Differential cross section ": "A1:Q25",
    "Elastic scattering cross sectio": "A1:C24",
    "momentum transfer cross section": "A1:F52",
    "vibrational excitation": "A1:C1802",
    "electronic excitation": "A1:O1502",
    "Dissociation": "A1:C15",
    "ionization": "A1:I45",
    "dissociative electron attachmen": "A1:C58",
}


@dataclass(frozen=True)
class Song2026O2CrossSectionCurve:
    """One exact numeric curve read from the user-supplied workbook."""

    label: str
    electron_energy_eV: np.ndarray
    cross_section_m2: np.ndarray
    uncertainty_m2: np.ndarray | None
    evidence_class: str

    def __post_init__(self):
        energy = np.asarray(self.electron_energy_eV, dtype=float).copy()
        sigma = np.asarray(self.cross_section_m2, dtype=float).copy()
        uncertainty = (
            None
            if self.uncertainty_m2 is None
            else np.asarray(self.uncertainty_m2, dtype=float).copy()
        )
        if (
            not self.label
            or not self.evidence_class
            or energy.ndim != 1
            or energy.size < 10
            or sigma.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(sigma))
            or np.any(energy < 0.0)
            or np.any(sigma < 0.0)
            or np.any(np.diff(energy) <= 0.0)
            or (
                uncertainty is not None
                and (
                    uncertainty.shape != energy.shape
                    or np.any(~np.isfinite(uncertainty))
                    or np.any(uncertainty < 0.0)
                )
            )
        ):
            raise ValueError("invalid Song 2026 O2 cross-section curve")
        energy.setflags(write=False)
        sigma.setflags(write=False)
        if uncertainty is not None:
            uncertainty.setflags(write=False)
        object.__setattr__(self, "electron_energy_eV", energy)
        object.__setattr__(self, "cross_section_m2", sigma)
        object.__setattr__(self, "uncertainty_m2", uncertainty)


@dataclass(frozen=True)
class Song2026O2Replay:
    """Hash-locked workbook replay and its explicitly closed kinetic deck."""

    source_curves: tuple[Song2026O2CrossSectionCurve, ...]
    derived_deck: ElectronCollisionDeck
    source_workbook_sha256: str
    high_energy_tail_closure: str
    dissociation_onset_closure: str
    maximum_energy_eV: float
    source_artifact_committed: bool = False
    source_license: str = "CC BY-NC 4.0"
    supports_resolved_primary_chemistry: bool = True
    supports_direct_swarm_validation: bool = False
    supports_target_reactor_state_prediction: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        if (
            len(self.source_curves) != 10
            or not isinstance(self.derived_deck, ElectronCollisionDeck)
            or self.source_workbook_sha256 != SONG_2026_O2_WORKBOOK_SHA256
            or self.high_energy_tail_closure not in {
                "inverse_energy", "linear_to_zero_at_30eV",
            }
            or self.dissociation_onset_closure not in {
                "linear_from_physical_threshold",
                "zero_until_first_tabulated_energy",
            }
            or self.maximum_energy_eV != 120.0
            or self.source_artifact_committed
            or self.source_license != "CC BY-NC 4.0"
            or not self.supports_resolved_primary_chemistry
            or self.supports_direct_swarm_validation
            or self.supports_target_reactor_state_prediction
            or self.supports_feature_depth
        ):
            raise ValueError("invalid Song 2026 O2 replay")

    def curve(self, label: str) -> Song2026O2CrossSectionCurve:
        matches = tuple(item for item in self.source_curves if item.label == label)
        if len(matches) != 1:
            raise KeyError(label)
        return matches[0]


def _workbook_cells(path: Path) -> dict[str, dict[str, str]]:
    """Read cached OOXML values without modifying or recalculating the file."""

    try:
        archive = ZipFile(path)
    except BadZipFile as exc:
        raise RuntimeError("Song 2026 O2 source is not a valid XLSX archive") from exc
    with archive:
        names = set(archive.namelist())
        required = {
            "xl/workbook.xml",
            "xl/_rels/workbook.xml.rels",
        }
        if not required <= names:
            raise RuntimeError("Song 2026 O2 workbook structure changed")
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in names:
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{{{_MAIN_NS}}}si"):
                shared_strings.append("".join(
                    node.text or ""
                    for node in item.iter(f"{{{_MAIN_NS}}}t")
                ))
        relationships = ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        targets = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in relationships.findall(f"{{{_PACKAGE_REL_NS}}}Relationship")
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        result: dict[str, dict[str, str]] = {}
        dimensions: dict[str, str] = {}
        sheets = workbook.find(f"{{{_MAIN_NS}}}sheets")
        if sheets is None:
            raise RuntimeError("Song 2026 O2 workbook has no sheets")
        for sheet in sheets:
            name = sheet.attrib["name"]
            relationship_id = sheet.attrib[f"{{{_DOC_REL_NS}}}id"]
            target = targets.get(relationship_id)
            if target is None:
                raise RuntimeError("Song 2026 O2 workbook relationship changed")
            member = str(PurePosixPath("xl") / target.lstrip("/"))
            if target.startswith("/"):
                member = target.lstrip("/")
            if member not in names:
                raise RuntimeError("Song 2026 O2 workbook sheet path changed")
            sheet_root = ET.fromstring(archive.read(member))
            dimension = sheet_root.find(f"{{{_MAIN_NS}}}dimension")
            if dimension is None:
                raise RuntimeError("Song 2026 O2 workbook dimension missing")
            dimensions[name] = dimension.attrib["ref"]
            cells: dict[str, str] = {}
            for cell in sheet_root.findall(f".//{{{_MAIN_NS}}}c"):
                reference = cell.attrib.get("r")
                if reference is None or _CELL_REF.fullmatch(reference) is None:
                    raise RuntimeError("Song 2026 O2 cell reference changed")
                value = cell.find(f"{{{_MAIN_NS}}}v")
                if value is None or value.text is None:
                    continue
                payload = value.text
                if cell.attrib.get("t") == "s":
                    try:
                        payload = shared_strings[int(payload)]
                    except (IndexError, ValueError) as exc:
                        raise RuntimeError(
                            "Song 2026 O2 shared-string table changed"
                        ) from exc
                cells[reference] = payload
            result[name] = cells
    if dimensions != _EXPECTED_DIMENSIONS:
        raise RuntimeError("Song 2026 O2 workbook sheet topology changed")
    return result


def _numeric_column(
    cells: dict[str, str],
    *,
    column: str,
    first_row: int,
    last_row: int,
) -> np.ndarray:
    values = []
    for row in range(first_row, last_row + 1):
        reference = f"{column}{row}"
        try:
            value = float(cells[reference])
        except (KeyError, ValueError) as exc:
            raise RuntimeError(
                f"Song 2026 O2 numeric cell {reference} changed"
            ) from exc
        if not math.isfinite(value):
            raise RuntimeError(f"Song 2026 O2 cell {reference} is non-finite")
        values.append(value)
    return np.asarray(values)


def _curve(
    cells: dict[str, str],
    *,
    label: str,
    energy_column: str,
    sigma_column: str,
    first_row: int,
    last_row: int,
    uncertainty_column: str | None,
    evidence_class: str,
) -> Song2026O2CrossSectionCurve:
    uncertainty = (
        None
        if uncertainty_column is None
        else _numeric_column(
            cells,
            column=uncertainty_column,
            first_row=first_row,
            last_row=last_row,
        ) * _M2_PER_SOURCE_UNIT
    )
    return Song2026O2CrossSectionCurve(
        label=label,
        electron_energy_eV=_numeric_column(
            cells,
            column=energy_column,
            first_row=first_row,
            last_row=last_row,
        ),
        cross_section_m2=_numeric_column(
            cells,
            column=sigma_column,
            first_row=first_row,
            last_row=last_row,
        ) * _M2_PER_SOURCE_UNIT,
        uncertainty_m2=uncertainty,
        evidence_class=evidence_class,
    )


def _read_source_curves(path: Path) -> tuple[Song2026O2CrossSectionCurve, ...]:
    if not path.is_file():
        raise FileNotFoundError(path)
    if sha256(path.read_bytes()).hexdigest() != SONG_2026_O2_WORKBOOK_SHA256:
        raise RuntimeError("Song 2026 O2 workbook checksum changed")
    sheets = _workbook_cells(path)
    signatures = {
        "momentum transfer cross section": ("A1", "MTCS"),
        "vibrational excitation": ("A1", "VI"),
        "electronic excitation": ("A1", "EX"),
        "Dissociation": ("A1", "DISS"),
        "ionization": ("A1", "ION"),
        "dissociative electron attachmen": ("A1", "DEA"),
    }
    for sheet, (cell, value) in signatures.items():
        if sheets[sheet].get(cell) != value:
            raise RuntimeError(f"Song 2026 O2 {sheet} signature changed")
    curves = [
        _curve(
            sheets["momentum transfer cross section"],
            label="momentum_transfer",
            energy_column="A", sigma_column="B", uncertainty_column="C",
            first_row=3, last_row=52,
            evidence_class="recommended_momentum_transfer",
        ),
        _curve(
            sheets["vibrational excitation"],
            label="vibrational_v0_to_v1",
            energy_column="A", sigma_column="B", uncertainty_column=None,
            first_row=3, last_row=1802,
            evidence_class="recommended_vibrational_excitation",
        ),
    ]
    electronic = sheets["electronic excitation"]
    for label, energy_column, sigma_column in (
        ("electronic_a1Delta_g", "A", "B"),
        ("electronic_b1Sigma_g_plus", "D", "E"),
        ("electronic_c1Sigma_u_minus", "G", "H"),
        ("electronic_A3Delta_u", "J", "K"),
        ("electronic_A3Sigma_u_plus", "M", "N"),
    ):
        curves.append(_curve(
            electronic,
            label=label,
            energy_column=energy_column,
            sigma_column=sigma_column,
            uncertainty_column=None,
            first_row=3,
            last_row=1502,
            evidence_class="recommended_electronic_excitation",
        ))
    curves.extend((
        _curve(
            sheets["Dissociation"],
            label="neutral_dissociation_total",
            energy_column="A", sigma_column="B", uncertainty_column="C",
            first_row=3, last_row=15,
            evidence_class="recommended_neutral_dissociation",
        ),
        _curve(
            sheets["ionization"],
            label="positive_ionization_total",
            energy_column="A", sigma_column="H", uncertainty_column="I",
            first_row=3, last_row=45,
            evidence_class="recommended_total_positive_ionization",
        ),
        _curve(
            sheets["dissociative electron attachmen"],
            label="dissociative_attachment_total",
            energy_column="A", sigma_column="B", uncertainty_column="C",
            first_row=3, last_row=58,
            evidence_class="recommended_dissociative_attachment",
        ),
    ))
    return tuple(curves)


def _tail(
    energy: np.ndarray,
    sigma: np.ndarray,
    *,
    maximum_energy_eV: float,
    closure: str,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if energy[-1] >= maximum_energy_eV:
        return tuple(energy), tuple(sigma)
    if closure == "inverse_energy":
        tail_energy = np.geomspace(energy[-1], maximum_energy_eV, 48)[1:]
        tail_sigma = sigma[-1] * energy[-1] / tail_energy
    elif closure == "linear_to_zero_at_30eV":
        zero_energy = max(30.0, float(energy[-1]) + 1.0)
        if zero_energy >= maximum_energy_eV:
            tail_energy = np.asarray([maximum_energy_eV])
            tail_sigma = np.asarray([0.0])
        else:
            tail_energy = np.asarray([zero_energy, maximum_energy_eV])
            tail_sigma = np.asarray([0.0, 0.0])
    else:
        raise ValueError("unsupported O2 high-energy tail closure")
    return (
        tuple(np.concatenate((energy, tail_energy))),
        tuple(np.concatenate((sigma, tail_sigma))),
    )


def _prepend_zero(
    energy: np.ndarray,
    sigma: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if energy[0] == 0.0:
        return energy.copy(), sigma.copy()
    return np.concatenate(([0.0], energy)), np.concatenate(([0.0], sigma))


def load_song_2026_o2_replay(
    path: str | Path,
    *,
    high_energy_tail_closure: str = "inverse_energy",
    dissociation_onset_closure: str = "zero_until_first_tabulated_energy",
    maximum_energy_eV: float = 120.0,
) -> Song2026O2Replay:
    """Load the exact external workbook and derive one closed O2 deck."""

    maximum = float(maximum_energy_eV)
    if maximum != 120.0:
        raise ValueError("the audited O2 deck currently supports exactly 120 eV")
    if high_energy_tail_closure not in {
        "inverse_energy", "linear_to_zero_at_30eV",
    }:
        raise ValueError("unsupported O2 high-energy tail closure")
    if dissociation_onset_closure not in {
        "linear_from_physical_threshold",
        "zero_until_first_tabulated_energy",
    }:
        raise ValueError("unsupported O2 dissociation onset closure")
    curves = _read_source_curves(Path(path))
    source = {curve.label: curve for curve in curves}

    momentum = source["momentum_transfer"]
    momentum_energy = np.concatenate(([0.0], momentum.electron_energy_eV))
    momentum_sigma = np.concatenate((
        [momentum.cross_section_m2[0]], momentum.cross_section_m2,
    ))
    processes = [ElectronCollisionProcess(
        kind="MOMENTUM",
        target="O2",
        product="O2",
        electron_energy_eV=tuple(momentum_energy),
        cross_section_m2=tuple(momentum_sigma),
        mass_ratio=ELECTRON_MASS_AMU / O2_MASS_AMU,
        comments=(
            "Song et al. 2026 recommended momentum-transfer cross section",
            "constant first source value from 0 to 0.001 eV",
        ),
    )]

    vibration = source["vibrational_v0_to_v1"]
    vib_energy, vib_sigma = _tail(
        vibration.electron_energy_eV,
        vibration.cross_section_m2,
        maximum_energy_eV=maximum,
        closure=high_energy_tail_closure,
    )
    processes.append(ElectronCollisionProcess(
        kind="EXCITATION",
        target="O2",
        product="O2(v=1)",
        electron_energy_eV=vib_energy,
        cross_section_m2=vib_sigma,
        energy_loss_eV=.196,
        comments=(
            "Song et al. 2026 v=0 to 1 vibrational excitation",
            f"high_energy_tail={high_energy_tail_closure}",
        ),
    ))

    electronic_channels = (
        ("electronic_a1Delta_g", "O2(a1Delta_g)", .98),
        ("electronic_b1Sigma_g_plus", "O2(b1Sigma_g+)", 1.63),
        ("electronic_c1Sigma_u_minus", "O2(c1Sigma_u-)", 6.12),
        ("electronic_A3Delta_u", "O2(A'3Delta_u)", 6.27),
        ("electronic_A3Sigma_u_plus", "O2(A3Sigma_u+)", 6.47),
    )
    for label, product, loss in electronic_channels:
        curve = source[label]
        energy, sigma = _prepend_zero(
            curve.electron_energy_eV, curve.cross_section_m2
        )
        energy_tail, sigma_tail = _tail(
            energy,
            sigma,
            maximum_energy_eV=maximum,
            closure=high_energy_tail_closure,
        )
        processes.append(ElectronCollisionProcess(
            kind="EXCITATION",
            target="O2",
            product=product,
            electron_energy_eV=energy_tail,
            cross_section_m2=sigma_tail,
            energy_loss_eV=loss,
            comments=(
                "Song et al. 2026 electronic excitation",
                f"high_energy_tail={high_energy_tail_closure}",
            ),
        ))

    dissociation = source["neutral_dissociation_total"]
    if dissociation_onset_closure == "linear_from_physical_threshold":
        diss_energy = np.concatenate((
            [0.0, 5.12], dissociation.electron_energy_eV,
        ))
        diss_sigma = np.concatenate((
            [0.0, 0.0], dissociation.cross_section_m2,
        ))
    else:
        below_first = np.nextafter(dissociation.electron_energy_eV[0], 0.0)
        diss_energy = np.concatenate((
            [0.0, 5.12, below_first], dissociation.electron_energy_eV,
        ))
        diss_sigma = np.concatenate((
            [0.0, 0.0, 0.0], dissociation.cross_section_m2,
        ))
    processes.append(ElectronCollisionProcess(
        kind="EXCITATION",
        target="O2",
        product="O + O aggregate",
        electron_energy_eV=tuple(diss_energy),
        cross_section_m2=tuple(diss_sigma),
        energy_loss_eV=5.12,
        comments=(
            "Song et al. 2026 recommended neutral dissociation",
            f"onset_closure={dissociation_onset_closure}",
            "aggregate product-state branching unresolved",
        ),
    ))

    ionization = source["positive_ionization_total"]
    ion_energy = np.concatenate((
        [0.0, 12.07], ionization.electron_energy_eV,
    ))
    ion_sigma = np.concatenate(([0.0, 0.0], ionization.cross_section_m2))
    processes.append(ElectronCollisionProcess(
        kind="IONIZATION",
        target="O2",
        product="O2+/O+ aggregate",
        electron_energy_eV=tuple(ion_energy),
        cross_section_m2=tuple(ion_sigma),
        energy_loss_eV=12.07,
        comments=(
            "Song et al. 2026 recommended total positive ionization",
            "workbook minus-sign typography rejected using article Table VII",
            "positive-ion product branching unresolved",
        ),
    ))

    attachment = source["dissociative_attachment_total"]
    attach_energy = np.concatenate((
        [0.0], attachment.electron_energy_eV, [15.0, maximum],
    ))
    attach_sigma = np.concatenate((
        [0.0], attachment.cross_section_m2, [0.0, 0.0],
    ))
    processes.append(ElectronCollisionProcess(
        kind="ATTACHMENT",
        target="O2",
        product="O- + O",
        electron_energy_eV=tuple(attach_energy),
        cross_section_m2=tuple(attach_sigma),
        comments=(
            "Song et al. 2026 recommended dissociative attachment",
            "linear closure from 9.9 eV to zero at 15 eV",
        ),
    ))

    derivation = {
        "schema": "petch.song_2026_o2_replay.v1",
        "source_workbook_sha256": SONG_2026_O2_WORKBOOK_SHA256,
        "source_artifact_committed": False,
        "source_license": "CC BY-NC 4.0",
        "maximum_energy_eV": maximum,
        "high_energy_tail_closure": high_energy_tail_closure,
        "dissociation_onset_closure": dissociation_onset_closure,
        "vibrational_energy_loss_eV": .196,
        "electronic_energy_losses_eV": [.98, 1.63, 6.12, 6.27, 6.47],
        "neutral_dissociation_energy_loss_eV": 5.12,
        "positive_ionization_energy_loss_eV": 12.07,
        "attachment_tail_zero_eV": 15.0,
        "rotational_excitation_omitted": True,
    }
    digest = sha256(json.dumps(
        derivation, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    deck = ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=digest,
        source_database="Song et al. 2026 recommended O2 collision set",
        retrieved_at="2026-08-18",
        source_reference=(
            f"article doi:{SONG_2026_O2_ARTICLE_DOI}; dataset "
            f"doi:{SONG_2026_O2_DATASET_DOI}; source workbook supplied "
            f"externally under CC BY-NC 4.0; derivation="
            f"{json.dumps(derivation, sort_keys=True)}"
        ),
    )
    return Song2026O2Replay(
        source_curves=curves,
        derived_deck=deck,
        source_workbook_sha256=SONG_2026_O2_WORKBOOK_SHA256,
        high_energy_tail_closure=high_energy_tail_closure,
        dissociation_onset_closure=dissociation_onset_closure,
        maximum_energy_eV=maximum,
    )
