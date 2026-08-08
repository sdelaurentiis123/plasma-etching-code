"""Rights-safe ingestion of user-supplied BOLSIG+/LXCat collision decks.

LXCat contributors retain ownership of their data and LXCat does not
authorize third-party redistribution.  This module therefore parses and
hash-locks a deck supplied by the user; it deliberately packages no LXCat
cross-section bytes.

Parsing a structurally complete collision set is only the input gate for an
electron Boltzmann solve.  It does not validate the set against swarm data and
cannot support a reactor, wafer-flux, or feature-depth prediction by itself.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re


BOLSIG_COLLISION_KINDS = frozenset({
    "ATTACHMENT",
    "EFFECTIVE",
    "ELASTIC",
    "EXCITATION",
    "IONIZATION",
    "MOMENTUM",
    "ROTATION",
})
_MOMENTUM_KINDS = frozenset({"EFFECTIVE", "ELASTIC", "MOMENTUM"})
_INELASTIC_KINDS = frozenset({
    "ATTACHMENT", "EXCITATION", "IONIZATION", "ROTATION",
})
_SEPARATOR = re.compile(r"^-{5,}$")


@dataclass(frozen=True)
class ElectronCollisionProcess:
    """One BOLSIG-format electron collision process in SI units."""

    kind: str
    target: str
    product: str | None
    electron_energy_eV: tuple[float, ...]
    cross_section_m2: tuple[float, ...]
    mass_ratio: float | None = None
    energy_loss_eV: float | None = None
    statistical_weight_ratio: float | None = None
    automatic_superelastic: bool = False
    comments: tuple[str, ...] = ()

    def __post_init__(self):
        energies = tuple(float(value) for value in self.electron_energy_eV)
        cross_sections = tuple(float(value) for value in self.cross_section_m2)
        comments = tuple(str(value).strip() for value in self.comments)
        if (
            self.kind not in BOLSIG_COLLISION_KINDS
            or not str(self.target).strip()
            or len(energies) < 2
            or len(energies) != len(cross_sections)
            or any(not math.isfinite(value) for value in energies)
            or any(not math.isfinite(value) for value in cross_sections)
            or any(value < 0.0 for value in energies)
            or any(value < 0.0 for value in cross_sections)
            or any(right <= left for left, right in zip(energies, energies[1:]))
        ):
            raise ValueError("invalid electron collision process")
        if self.kind in _MOMENTUM_KINDS:
            if (
                self.mass_ratio is None
                or not math.isfinite(self.mass_ratio)
                or self.mass_ratio <= 0.0
                or self.energy_loss_eV is not None
                or self.statistical_weight_ratio is not None
                or self.automatic_superelastic
            ):
                raise ValueError("invalid momentum-transfer process metadata")
        elif self.kind == "ATTACHMENT":
            if (
                self.mass_ratio is not None
                or self.energy_loss_eV not in (None, 0.0)
                or self.statistical_weight_ratio is not None
                or self.automatic_superelastic
            ):
                raise ValueError("invalid attachment process metadata")
            object.__setattr__(self, "energy_loss_eV", 0.0)
        elif self.kind in _INELASTIC_KINDS:
            if (
                self.mass_ratio is not None
                or self.energy_loss_eV is None
                or not math.isfinite(self.energy_loss_eV)
                or (
                    self.kind in {"IONIZATION", "ROTATION"}
                    and self.energy_loss_eV < 0.0
                )
                or (
                    self.automatic_superelastic
                    and (
                        self.statistical_weight_ratio is None
                        or not math.isfinite(self.statistical_weight_ratio)
                        or self.statistical_weight_ratio <= 0.0
                    )
                )
                or (
                    not self.automatic_superelastic
                    and self.statistical_weight_ratio is not None
                )
            ):
                raise ValueError("invalid inelastic process metadata")
        object.__setattr__(self, "target", str(self.target).strip())
        object.__setattr__(
            self, "product",
            None if self.product is None else str(self.product).strip(),
        )
        object.__setattr__(self, "electron_energy_eV", energies)
        object.__setattr__(self, "cross_section_m2", cross_sections)
        object.__setattr__(self, "comments", comments)


@dataclass(frozen=True)
class ElectronCollisionDeck:
    """A hash-locked local collision deck with explicit use boundaries."""

    processes: tuple[ElectronCollisionProcess, ...]
    payload_sha256: str
    source_database: str
    retrieved_at: str
    source_reference: str
    packaged_or_redistributed: bool = False

    def __post_init__(self):
        processes = tuple(self.processes)
        if (
            not processes
            or any(
                not isinstance(item, ElectronCollisionProcess)
                for item in processes
            )
            or not re.fullmatch(r"[0-9a-f]{64}", self.payload_sha256)
            or not str(self.source_database).strip()
            or not str(self.retrieved_at).strip()
            or not str(self.source_reference).strip()
            or self.packaged_or_redistributed
        ):
            raise ValueError("invalid local electron collision deck")
        object.__setattr__(self, "processes", processes)

    @property
    def targets(self) -> tuple[str, ...]:
        return tuple(sorted({item.target for item in self.processes}))

    def for_target(self, target: str) -> tuple[ElectronCollisionProcess, ...]:
        target = str(target).strip()
        if not target:
            raise ValueError("target must be non-empty")
        return tuple(item for item in self.processes if item.target == target)

    def structural_kinetic_readiness(
        self,
        target: str,
        *,
        elastic_angular_closure: str | None = None,
    ) -> dict[str, object]:
        """Report kinetic-input topology, never a validation verdict.

        A BOLSIG/LXCat momentum-transfer row does not contain the higher
        Legendre moments of a differential elastic cross section.  Multi-term
        use therefore fails closed unless an angular closure is explicit.
        ``isotropic_source_reproduction`` is the declared assumption used by
        Kawaguchi et al.; it is not promoted to differential-scattering data.
        """

        if elastic_angular_closure not in {
            None, "isotropic_source_reproduction",
        }:
            raise ValueError("unsupported elastic angular closure")

        processes = self.for_target(target)
        momentum = tuple(
            item for item in processes if item.kind in _MOMENTUM_KINDS)
        inelastic = tuple(
            item for item in processes if item.kind in _INELASTIC_KINDS)
        ionization = tuple(
            item for item in processes if item.kind == "IONIZATION")
        attachment = tuple(
            item for item in processes if item.kind == "ATTACHMENT")
        issues = []
        if not processes:
            issues.append("target_missing")
        if len(momentum) != 1:
            issues.append("requires_exactly_one_momentum_or_effective_process")
        if not inelastic:
            issues.append("no_inelastic_energy_loss_process")
        if elastic_angular_closure is None:
            issues.append("multiterm_requires_explicit_elastic_angular_closure")
        return {
            "target": target,
            "process_count": len(processes),
            "momentum_process_count": len(momentum),
            "inelastic_process_count": len(inelastic),
            "ionization_process_count": len(ionization),
            "attachment_process_count": len(attachment),
            "structurally_ready_for_kinetic_input": not issues,
            "elastic_angular_closure": elastic_angular_closure,
            "angular_evidence_class": (
                "undeclared"
                if elastic_angular_closure is None
                else "source_reproduction_assumption"
            ),
            "contains_differential_elastic_cross_sections": False,
            "issues": tuple(issues),
            "supports_swarm_validation": False,
            "supports_reactor_state_prediction": False,
            "supports_wafer_flux": False,
            "supports_feature_depth": False,
        }


def _is_separator(line: str) -> bool:
    return bool(_SEPARATOR.fullmatch(line.strip()))


def _next_nonempty(lines: tuple[str, ...], index: int) -> tuple[int, str]:
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index >= len(lines):
        raise ValueError("unexpected end of BOLSIG collision deck")
    return index, lines[index].strip()


def _parameter_values(line: str) -> tuple[float, ...]:
    payload = line.split("/", 1)[0].strip()
    try:
        values = tuple(float(token) for token in payload.split())
    except ValueError as exc:
        raise ValueError("invalid BOLSIG process parameter line") from exc
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("invalid BOLSIG process parameter line")
    return values


def _signature(value: str) -> tuple[str, str | None, bool]:
    if "<->" in value:
        target, product = value.split("<->", 1)
        return target.strip(), product.strip(), True
    if "->" in value:
        target, product = value.split("->", 1)
        return target.strip(), product.strip(), False
    return value.strip(), None, False


def parse_bolsig_lxcat_bytes(
    payload: bytes,
    *,
    source_database: str,
    retrieved_at: str,
    source_reference: str,
    target: str | None = None,
    expected_sha256: str | None = None,
) -> ElectronCollisionDeck:
    """Parse a user-supplied BOLSIG/LXCat deck without redistributing it."""

    if not isinstance(payload, bytes) or not payload:
        raise ValueError("collision deck payload must be non-empty bytes")
    digest = hashlib.sha256(payload).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise RuntimeError("electron collision deck hash mismatch")
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("collision deck must be UTF-8 text") from exc
    lines = tuple(text.splitlines())
    processes: list[ElectronCollisionProcess] = []
    index = 0
    while index < len(lines):
        current = lines[index].strip()
        if current == "COMMENT":
            index += 1
            while index < len(lines) and not _is_separator(lines[index]):
                index += 1
            index += 1
            continue
        if current not in BOLSIG_COLLISION_KINDS:
            index += 1
            continue
        kind = current
        index, raw_signature = _next_nonempty(lines, index + 1)
        process_target, product, automatic_superelastic = _signature(
            raw_signature)
        index += 1
        mass_ratio = None
        energy_loss_eV = None
        statistical_weight_ratio = None
        if kind != "ATTACHMENT":
            index, parameter_line = _next_nonempty(lines, index)
            values = _parameter_values(parameter_line)
            index += 1
            if kind in _MOMENTUM_KINDS:
                if len(values) != 1:
                    raise ValueError(
                        "momentum process requires one mass-ratio value")
                mass_ratio = values[0]
            else:
                if len(values) not in (1, 2):
                    raise ValueError(
                        "inelastic process requires threshold and optional "
                        "statistical-weight ratio"
                    )
                energy_loss_eV = values[0]
                if len(values) == 2:
                    statistical_weight_ratio = values[1]
        comments = []
        while index < len(lines) and not _is_separator(lines[index]):
            value = lines[index].strip()
            if value:
                comments.append(value)
            index += 1
        if index >= len(lines):
            raise ValueError("collision process has no data-table delimiter")
        index += 1
        energies = []
        cross_sections = []
        while index < len(lines) and not _is_separator(lines[index]):
            value = lines[index].strip()
            index += 1
            if not value:
                continue
            fields = value.split()
            if len(fields) != 2:
                raise ValueError("collision data row must have two numbers")
            try:
                energy, cross_section = map(float, fields)
            except ValueError as exc:
                raise ValueError("invalid collision data row") from exc
            energies.append(energy)
            cross_sections.append(cross_section)
        if index >= len(lines):
            raise ValueError("collision process data table is unterminated")
        index += 1
        if target is None or process_target == target:
            processes.append(ElectronCollisionProcess(
                kind=kind,
                target=process_target,
                product=product,
                electron_energy_eV=tuple(energies),
                cross_section_m2=tuple(cross_sections),
                mass_ratio=mass_ratio,
                energy_loss_eV=energy_loss_eV,
                statistical_weight_ratio=statistical_weight_ratio,
                automatic_superelastic=automatic_superelastic,
                comments=tuple(comments),
            ))
    if not processes:
        suffix = "" if target is None else f" for target {target!r}"
        raise ValueError(f"collision deck contains no supported processes{suffix}")
    return ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=digest,
        source_database=source_database,
        retrieved_at=retrieved_at,
        source_reference=source_reference,
        packaged_or_redistributed=False,
    )


def load_bolsig_lxcat_file(
    path: str | Path,
    **metadata,
) -> ElectronCollisionDeck:
    """Load a local user-supplied deck; the bytes remain outside the package."""

    return parse_bolsig_lxcat_bytes(Path(path).read_bytes(), **metadata)
