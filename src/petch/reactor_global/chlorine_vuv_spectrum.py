"""Fail-closed neutral-chlorine VUV spectrum sensitivity from ADF04 data.

This module supplies deterministic atomic bookkeeping, not a fitted photon
boundary.  It combines observed NIST level separations with effective
collision strengths and radiative branching from an independent
AUTOSTRUCTURE distorted-wave calculation distributed through OPEN-ADAS.
The hybrid is necessary because the calculation's raw Cl I level energies are
not spectroscopic: its first excited term is placed above the measured
ionisation threshold.

OPEN-ADAS files are license restricted and are therefore never bundled with
petch.  Loading requires an explicit personal-research acknowledgement, a
caller-provided path, and a pinned physical-record hash.  The resulting direct
coronal spectrum is an audit sensitivity only: cascades, resonance escape,
plasma quenching, and wavelength-resolved Si photoetch yields remain open.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re
from typing import Iterable

import numpy as np


OPEN_ADAS_CL0_COLLISION_RECORDS_SHA256 = (
    "f215093e9ab5ac36a202cdb353a0e3a3b9651982158ba105531fb44cab74c4e7"
)
OPEN_ADAS_CL0_NIST_RECORDS_SHA256 = (
    "12905140158c763f9a1cc6efff4d27fb74ecbf77e2162c67726b8dc528d1f430"
)
OPEN_ADAS_PERSONAL_USE_NOTICE = (
    "OPEN-ADAS downloads are restricted to personal use and may not be "
    "redistributed with this code without written ADAS Project permission."
)

_FORTRAN_NUMBER = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[+-]\d+)?"
)
_LEVEL = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+\(\s*(\d*)\)\s*(\S*)\s*"
    r"\(\s*([0-9.]+)\)\s+([+-]?[0-9.]+)"
)
_WAVENUMBER_TO_NM = 1.0e7
_BOLTZMANN_EV_K = 8.617333262145e-5
_MAXWELLIAN_RATE_CONSTANT_CM3_K05_S = 8.629e-6


def _adas_float(token: str) -> float:
    """Parse ADAS's exponent-without-E notation (for example ``4.09-02``)."""
    match = re.fullmatch(
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))(?:([+-]\d+))?", token.strip()
    )
    if match is None:
        raise ValueError(f"invalid ADF04 number {token!r}")
    mantissa, exponent = match.groups()
    return float(f"{mantissa}e{exponent or '0'}")


def _physical_record_bytes(raw: bytes) -> bytes:
    """Return normalized ADF04 records through the first ``-1, -1`` marker."""
    lines = raw.decode("ascii").splitlines()
    for index, line in enumerate(lines):
        try:
            upper = int(line[1:6])
            lower = int(line[6:11])
        except (ValueError, IndexError):
            continue
        if upper == -1 and lower == -1:
            return ("\n".join(lines[: index + 1]) + "\n").encode("ascii")
    raise ValueError("ADF04 physical-record terminator '-1 -1' is missing")


@dataclass(frozen=True)
class Adf04Level:
    index: int
    configuration: str
    multiplicity: int
    orbital_label: str
    total_angular_momentum: float
    energy_cm_inv: float

    @property
    def statistical_weight(self) -> float:
        return 2.0 * self.total_angular_momentum + 1.0

    @property
    def state_key(self) -> tuple[str, int, str, float]:
        configuration = self.configuration
        # The NIST converter prepends the ion-core token 521; the independent
        # AUTOSTRUCTURE file does not.  The remaining compact occupancy string
        # is identical for the resolved states used here.
        if configuration.startswith("521"):
            configuration = configuration[3:]
        return (
            configuration,
            self.multiplicity,
            self.orbital_label,
            self.total_angular_momentum,
        )


@dataclass(frozen=True)
class Adf04Type3Transition:
    upper_index: int
    lower_index: int
    transition_probability_s_inv: float
    effective_collision_strengths: tuple[float, ...]
    high_energy_parameter: float


@dataclass(frozen=True)
class Adf04Type3Dataset:
    levels: tuple[Adf04Level, ...]
    temperatures_K: tuple[float, ...]
    transitions: tuple[Adf04Type3Transition, ...]
    physical_records_sha256: str


def _parse_level(line: str) -> Adf04Level | None:
    match = _LEVEL.match(line)
    if match is None:
        return None
    index, configuration, multiplicity, orbital, angular, energy = match.groups()
    if not multiplicity or not orbital:
        # Unresolved NIST blends cannot be joined to fine-structure collision
        # states and are deliberately excluded.
        return None
    return Adf04Level(
        index=int(index),
        configuration=configuration,
        multiplicity=int(multiplicity),
        orbital_label=orbital,
        total_angular_momentum=float(angular),
        energy_cm_inv=float(energy),
    )


def _first_level_terminator(lines: list[str]) -> int:
    for index, line in enumerate(lines[1:], start=1):
        try:
            if int(line[:6]) == -1:
                return index
        except ValueError:
            continue
    raise ValueError("ADF04 level terminator is missing")


def parse_adf04_type3_bytes(raw: bytes) -> Adf04Type3Dataset:
    """Parse the ADF04 type-3 subset needed for a coronal line audit.

    Transition indices occupy fixed fields.  Numeric tokenization begins only
    at column 12 because adjacent signed ADAS exponents such as
    ``4.09-02-5.33-04`` are legal and cannot be parsed with ``split()``.
    """
    physical = _physical_record_bytes(raw)
    lines = physical.decode("ascii").splitlines()
    first_terminator = _first_level_terminator(lines)
    levels = tuple(
        level for level in (_parse_level(line) for line in lines[1:first_terminator])
        if level is not None
    )
    header_index = first_terminator + 1
    while header_index < len(lines) and not lines[header_index].strip():
        header_index += 1
    header_values = [_adas_float(item) for item in _FORTRAN_NUMBER.findall(
        lines[header_index]
    )]
    if len(header_values) < 3 or int(header_values[1]) != 3:
        raise ValueError("only an ADF04 type-3 collision block is supported")
    temperatures = tuple(header_values[2:])
    transitions: list[Adf04Type3Transition] = []
    for line in lines[header_index + 1:]:
        try:
            upper = int(line[1:6])
        except (ValueError, IndexError) as error:
            raise ValueError(f"invalid fixed-width ADF04 transition: {line!r}") from error
        if upper == -1:
            break
        try:
            lower = int(line[6:11])
        except (ValueError, IndexError) as error:
            raise ValueError(f"invalid fixed-width ADF04 transition: {line!r}") from error
        values = tuple(_adas_float(item) for item in _FORTRAN_NUMBER.findall(line[11:]))
        if len(values) != len(temperatures) + 2:
            raise ValueError(
                f"transition {upper}->{lower} carries {len(values)} values; "
                f"expected {len(temperatures) + 2}"
            )
        transitions.append(Adf04Type3Transition(
            upper_index=upper,
            lower_index=lower,
            transition_probability_s_inv=values[0],
            effective_collision_strengths=values[1:-1],
            high_energy_parameter=values[-1],
        ))
    if not levels or not temperatures or not transitions:
        raise ValueError("incomplete ADF04 type-3 dataset")
    return Adf04Type3Dataset(
        levels=levels,
        temperatures_K=temperatures,
        transitions=tuple(transitions),
        physical_records_sha256=hashlib.sha256(physical).hexdigest(),
    )


def parse_adf04_levels_bytes(raw: bytes) -> tuple[Adf04Level, ...]:
    """Parse a level-only ADF04 file, excluding unresolved NIST blends."""
    physical = _physical_record_bytes(raw)
    lines = physical.decode("ascii").splitlines()
    first_terminator = _first_level_terminator(lines)
    levels = tuple(
        level for level in (_parse_level(line) for line in lines[1:first_terminator])
        if level is not None
    )
    if not levels:
        raise ValueError("ADF04 level file contains no resolved levels")
    return levels


def load_open_adas_cl0_personal_research(
    collision_path: str | Path,
    nist_level_path: str | Path,
    *,
    accept_restricted_personal_use: bool = False,
) -> tuple[Adf04Type3Dataset, tuple[Adf04Level, ...]]:
    """Load the two pinned OPEN-ADAS Cl I files without redistributing them."""
    if not accept_restricted_personal_use:
        raise PermissionError(OPEN_ADAS_PERSONAL_USE_NOTICE)
    collision_raw = Path(collision_path).read_bytes()
    nist_raw = Path(nist_level_path).read_bytes()
    collision = parse_adf04_type3_bytes(collision_raw)
    nist_physical = _physical_record_bytes(nist_raw)
    nist_hash = hashlib.sha256(nist_physical).hexdigest()
    if collision.physical_records_sha256 != OPEN_ADAS_CL0_COLLISION_RECORDS_SHA256:
        raise ValueError("OPEN-ADAS Cl I collision physical-record hash mismatch")
    if nist_hash != OPEN_ADAS_CL0_NIST_RECORDS_SHA256:
        raise ValueError("OPEN-ADAS Cl I NIST-level physical-record hash mismatch")
    return collision, parse_adf04_levels_bytes(nist_raw)


def map_observed_to_collision_levels(
    collision_levels: Iterable[Adf04Level],
    observed_levels: Iterable[Adf04Level],
) -> dict[int, int]:
    """Join repeated fine-structure states by identity and energy-order rank."""
    calculated: dict[tuple[str, int, str, float], list[Adf04Level]] = defaultdict(list)
    observed: dict[tuple[str, int, str, float], list[Adf04Level]] = defaultdict(list)
    for level in collision_levels:
        calculated[level.state_key].append(level)
    for level in observed_levels:
        observed[level.state_key].append(level)
    mapping: dict[int, int] = {}
    for key in calculated.keys() & observed.keys():
        calculated_group = sorted(
            calculated[key], key=lambda level: (level.energy_cm_inv, level.index)
        )
        observed_group = sorted(
            observed[key], key=lambda level: (level.energy_cm_inv, level.index)
        )
        for measured, model in zip(observed_group, calculated_group):
            mapping[measured.index] = model.index
    return mapping


def _interpolate_upsilon(
    transition: Adf04Type3Transition,
    temperatures_K: tuple[float, ...],
    electron_temperature_K: float,
) -> float:
    grid = np.asarray(temperatures_K, dtype=float)
    if not grid[0] <= electron_temperature_K <= grid[-1]:
        raise ValueError("electron temperature lies outside the ADF04 grid")
    return float(np.interp(
        math.log(electron_temperature_K),
        np.log(grid),
        np.asarray(transition.effective_collision_strengths, dtype=float),
    ))


@dataclass(frozen=True)
class ChlorineVuvLine:
    upper_observed_index: int
    lower_observed_index: int
    wavelength_nm: float
    transition_probability_s_inv: float
    upper_total_radiative_probability_s_inv: float
    direct_excitation_rate_coefficient_cm3_s: float
    photon_rate_coefficient_cm3_s: float


@dataclass(frozen=True)
class ChlorineDirectCoronalSpectrum:
    electron_temperature_eV: float
    nonradiative_loss_s_inv: float
    lines: tuple[ChlorineVuvLine, ...]
    matched_observed_level_count: int
    calculated_level_count: int
    observed_level_count: int
    prediction_supported: bool = False

    def band_rate_coefficient_cm3_s(
        self, minimum_wavelength_nm: float, maximum_wavelength_nm: float,
    ) -> float:
        lower = float(minimum_wavelength_nm)
        upper = float(maximum_wavelength_nm)
        if not 0.0 < lower < upper:
            raise ValueError("invalid wavelength band")
        return float(sum(
            line.photon_rate_coefficient_cm3_s for line in self.lines
            if lower <= line.wavelength_nm < upper
        ))


def chlorine_direct_coronal_spectrum(
    collision: Adf04Type3Dataset,
    observed_levels: Iterable[Adf04Level],
    *,
    electron_temperature_eV: float,
    nonradiative_loss_s_inv: float = 0.0,
) -> ChlorineDirectCoronalSpectrum:
    """Evaluate direct ground-term excitation and radiative branching.

    The neutral density multiplying the returned coefficients is the total
    ground-term Cl density.  Its two fine-structure levels are distributed by
    their observed Boltzmann weights.  Excited-state cascades and radiation
    trapping are intentionally not inferred.
    """
    temperature_eV = float(electron_temperature_eV)
    nonradiative_loss = float(nonradiative_loss_s_inv)
    if (
        not math.isfinite(temperature_eV)
        or temperature_eV <= 0.0
        or not math.isfinite(nonradiative_loss)
        or nonradiative_loss < 0.0
    ):
        raise ValueError("invalid coronal-spectrum condition")
    temperature_K = temperature_eV / _BOLTZMANN_EV_K
    observed_tuple = tuple(observed_levels)
    observed_by_index = {level.index: level for level in observed_tuple}
    calculated_by_index = {level.index: level for level in collision.levels}
    observed_to_calculated = map_observed_to_collision_levels(
        collision.levels, observed_tuple
    )
    calculated_to_observed = {
        calculated_index: observed_index
        for observed_index, calculated_index in observed_to_calculated.items()
    }
    ground_observed = sorted(
        (observed_by_index[index] for index in (1, 2)),
        key=lambda level: level.index,
    )
    if len(ground_observed) != 2 or not all(
        level.index in observed_to_calculated for level in ground_observed
    ):
        raise ValueError("both observed Cl I ground fine-structure levels are required")
    ground_weights = np.asarray([
        level.statistical_weight * math.exp(
            -level.energy_cm_inv
            * 1.2398419843320026e-4 / temperature_eV
        )
        for level in ground_observed
    ])
    ground_weights /= ground_weights.sum()

    transition_lookup = {
        (transition.upper_index, transition.lower_index): transition
        for transition in collision.transitions
    }
    radiative_out: dict[int, float] = defaultdict(float)
    for transition in collision.transitions:
        radiative_out[transition.upper_index] += max(
            0.0, transition.transition_probability_s_inv
        )

    direct_excitation: dict[int, float] = defaultdict(float)
    for population, ground in zip(ground_weights, ground_observed):
        calculated_ground = observed_to_calculated[ground.index]
        for observed_upper in observed_tuple:
            if observed_upper.energy_cm_inv <= ground.energy_cm_inv:
                continue
            calculated_upper = observed_to_calculated.get(observed_upper.index)
            transition = transition_lookup.get((calculated_upper, calculated_ground))
            if transition is None:
                continue
            upsilon = _interpolate_upsilon(
                transition, collision.temperatures_K, temperature_K
            )
            delta_energy_eV = (
                observed_upper.energy_cm_inv - ground.energy_cm_inv
            ) * 1.2398419843320026e-4
            rate = (
                _MAXWELLIAN_RATE_CONSTANT_CM3_K05_S
                * upsilon
                / (ground.statistical_weight * math.sqrt(temperature_K))
                * math.exp(-delta_energy_eV / temperature_eV)
            )
            direct_excitation[calculated_upper] += float(population * rate)

    lines: list[ChlorineVuvLine] = []
    for transition in collision.transitions:
        upper_observed_index = calculated_to_observed.get(transition.upper_index)
        lower_observed_index = calculated_to_observed.get(transition.lower_index)
        if upper_observed_index is None or lower_observed_index is None:
            continue
        upper_observed = observed_by_index[upper_observed_index]
        lower_observed = observed_by_index[lower_observed_index]
        separation = upper_observed.energy_cm_inv - lower_observed.energy_cm_inv
        excitation = direct_excitation.get(transition.upper_index, 0.0)
        total_radiative = radiative_out.get(transition.upper_index, 0.0)
        if (
            separation <= 0.0
            or excitation <= 0.0
            or transition.transition_probability_s_inv <= 0.0
            or total_radiative <= 0.0
        ):
            continue
        photon_rate = excitation * transition.transition_probability_s_inv / (
            total_radiative + nonradiative_loss
        )
        lines.append(ChlorineVuvLine(
            upper_observed_index=upper_observed_index,
            lower_observed_index=lower_observed_index,
            wavelength_nm=_WAVENUMBER_TO_NM / separation,
            transition_probability_s_inv=transition.transition_probability_s_inv,
            upper_total_radiative_probability_s_inv=total_radiative,
            direct_excitation_rate_coefficient_cm3_s=excitation,
            photon_rate_coefficient_cm3_s=photon_rate,
        ))
    return ChlorineDirectCoronalSpectrum(
        electron_temperature_eV=temperature_eV,
        nonradiative_loss_s_inv=nonradiative_loss,
        lines=tuple(sorted(lines, key=lambda line: line.wavelength_nm)),
        matched_observed_level_count=len(observed_to_calculated),
        calculated_level_count=len(collision.levels),
        observed_level_count=len(observed_tuple),
    )
