"""Measured reactive CF3+--CHF3 collision kernel.

Peko et al. measured two *destruction* mechanisms for a CF3+ projectile in
CHF3: summed dissociative charge transfer (DCT) and summed collision-induced
dissociation (CID).  The DCT curve is reconstructed from Figure 2; the paper
reports the summed CID cross section as approximately energy independent at
18 A^2 over 20--225 eV relative collision energy.

This module deliberately does not invent the missing elastic differential
cross section.  It therefore closes a measured reactive-removal floor, not a
complete molecular-ion Boltzmann operator or target-reactor IEAD.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
from hashlib import sha256
import math
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.interpolate import PchipInterpolator


ANGSTROM2_TO_M2 = 1.0e-20
PEKO_DCT_CSV_SHA256 = (
    "e10cc18dfa63e8fb79ccd386fba86bacc019b37f3714563b044b6b10c05105bb"
)
_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PEKO_DCT_CSV = (
    _ROOT / "data" / "experimental" / "peko_2002_chf3"
    / "figure2_cf3_chf3_dct_sum.csv"
)


def _scalar_or_array(value: np.ndarray) -> float | np.ndarray:
    return float(value) if value.ndim == 0 else value


@dataclass(frozen=True)
class ReactiveCrossSectionBand:
    """Central measured cross section and conservative source-error band."""

    central_m2: float | np.ndarray
    lower_m2: float | np.ndarray
    upper_m2: float | np.ndarray

    def __post_init__(self):
        central = np.asarray(self.central_m2, dtype=float)
        lower = np.asarray(self.lower_m2, dtype=float)
        upper = np.asarray(self.upper_m2, dtype=float)
        if (
            central.shape != lower.shape
            or central.shape != upper.shape
            or np.any(~np.isfinite(central))
            or np.any(~np.isfinite(lower))
            or np.any(~np.isfinite(upper))
            or np.any(lower < 0.0)
            or np.any(lower > central)
            or np.any(central > upper)
        ):
            raise ValueError("invalid reactive cross-section band")


@dataclass(frozen=True)
class ReactiveSlabSensitivity:
    """One-energy slab scale check, explicitly not a sheath solution."""

    laboratory_energy_eV: float
    relative_collision_energy_eV: float
    path_length_m: float
    chf3_number_density_m3: float
    cross_section: ReactiveCrossSectionBand
    optical_depth_central: float
    optical_depth_lower: float
    optical_depth_upper: float
    destruction_probability_central: float
    destruction_probability_lower: float
    destruction_probability_upper: float
    feed_fraction_used_as_density_proxy: bool
    supports_complete_molecular_transport: bool = False
    supports_target_iead: bool = False
    supports_absolute_depth_prediction: bool = False

    def __post_init__(self):
        values = np.asarray([
            self.laboratory_energy_eV,
            self.relative_collision_energy_eV,
            self.path_length_m,
            self.chf3_number_density_m3,
            self.optical_depth_central,
            self.optical_depth_lower,
            self.optical_depth_upper,
            self.destruction_probability_central,
            self.destruction_probability_lower,
            self.destruction_probability_upper,
        ])
        if (
            np.any(~np.isfinite(values))
            or np.any(values[:4] <= 0.0)
            or not 0.0 <= self.optical_depth_lower <= self.optical_depth_central
            or not self.optical_depth_central <= self.optical_depth_upper
            or not 0.0 <= self.destruction_probability_lower <= 1.0
            or not (
                self.destruction_probability_lower
                <= self.destruction_probability_central
                <= self.destruction_probability_upper
                <= 1.0
            )
            or self.supports_complete_molecular_transport
            or self.supports_target_iead
            or self.supports_absolute_depth_prediction
        ):
            raise ValueError("invalid reactive slab sensitivity")


@dataclass(frozen=True)
class Peko2002CF3CHF3ReactiveCollisionModel:
    """C1 measured CF3+ destruction floor for stationary CHF3 targets.

    Peko's abscissa is *relative* collision energy.  For a projectile of
    laboratory energy ``E_lab`` striking a stationary neutral,

    ``E_rel = m_target / (m_projectile + m_target) * E_lab``.

    The conversion is enforced here so a sheath energy cannot be silently
    inserted as a relative energy.  PCHIP in log-energy is deterministic,
    positive on the data support, local, and continuously differentiable.
    Extrapolation is forbidden.
    """

    relative_energy_eV: np.ndarray
    dct_sum_cross_section_A2: np.ndarray
    projectile_mass_amu: float = 69.0
    target_mass_amu: float = 70.0
    cid_sum_cross_section_A2: float = 18.0
    cid_relative_energy_support_eV: tuple[float, float] = (20.0, 225.0)
    dct_source_relative_uncertainty: float = 0.25
    cid_source_relative_uncertainty: float = 0.15
    digitization_cross_section_A2_bound: float = 0.046677327
    source: str = "peko-2002-chf3-ion-molecule"
    provenance: Mapping[str, object] = field(default_factory=lambda: {
        "dct": "Figure 2 summed DCT asterisks, checksum-pinned digitization",
        "cid": "text-reported summed Reactions 18--21, about 18 A^2",
        "collision_energy": "source relative collision energy",
        "missing": (
            "elastic/momentum-transfer cross section and angular kernel; "
            "other molecular-ion/neutral pairs"
        ),
        "coefficient_selected_from_depth_target": None,
    })

    def __post_init__(self):
        energy = np.asarray(self.relative_energy_eV, dtype=float).copy()
        dct = np.asarray(self.dct_sum_cross_section_A2, dtype=float).copy()
        scalars = np.asarray([
            self.projectile_mass_amu,
            self.target_mass_amu,
            self.cid_sum_cross_section_A2,
            *self.cid_relative_energy_support_eV,
            self.dct_source_relative_uncertainty,
            self.cid_source_relative_uncertainty,
            self.digitization_cross_section_A2_bound,
        ])
        if (
            energy.ndim != 1
            or energy.size < 8
            or dct.shape != energy.shape
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(dct))
            or np.any(energy <= 0.0)
            or np.any(dct <= 0.0)
            or np.any(np.diff(energy) <= 0.0)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars <= 0.0)
            or self.cid_relative_energy_support_eV[0]
            >= self.cid_relative_energy_support_eV[1]
            or self.dct_source_relative_uncertainty >= 1.0
            or self.cid_source_relative_uncertainty >= 1.0
            or energy[0] < self.cid_relative_energy_support_eV[0]
            or energy[-1] > self.cid_relative_energy_support_eV[1]
            or not str(self.source).strip()
        ):
            raise ValueError("invalid Peko CF3+--CHF3 collision model")
        energy.setflags(write=False)
        dct.setflags(write=False)
        object.__setattr__(self, "relative_energy_eV", energy)
        object.__setattr__(self, "dct_sum_cross_section_A2", dct)
        object.__setattr__(self, "provenance", MappingProxyType(
            dict(self.provenance)))
        object.__setattr__(self, "_dct_pchip", PchipInterpolator(
            np.log(energy), dct, extrapolate=False))

    @property
    def relative_energy_support_eV(self) -> tuple[float, float]:
        return float(self.relative_energy_eV[0]), float(
            self.relative_energy_eV[-1])

    @property
    def laboratory_energy_support_eV(self) -> tuple[float, float]:
        factor = (
            self.projectile_mass_amu + self.target_mass_amu
        ) / self.target_mass_amu
        lower, upper = self.relative_energy_support_eV
        return factor * lower, factor * upper

    @property
    def laboratory_to_relative_energy_factor(self) -> float:
        return self.target_mass_amu / (
            self.projectile_mass_amu + self.target_mass_amu)

    def relative_energy_from_laboratory_eV(
        self, laboratory_energy_eV,
    ) -> float | np.ndarray:
        energy = np.asarray(laboratory_energy_eV, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy <= 0.0):
            raise ValueError("laboratory collision energy must be positive")
        result = self.laboratory_to_relative_energy_factor * energy
        return _scalar_or_array(result)

    def _checked_relative_energy(self, relative_energy_eV) -> np.ndarray:
        energy = np.asarray(relative_energy_eV, dtype=float)
        lower, upper = self.relative_energy_support_eV
        if (
            np.any(~np.isfinite(energy))
            or np.any(energy < lower)
            or np.any(energy > upper)
        ):
            raise ValueError(
                "relative collision energy is outside measured DCT support")
        return energy

    def dct_sum_cross_section_m2(
        self, relative_energy_eV,
    ) -> float | np.ndarray:
        energy = self._checked_relative_energy(relative_energy_eV)
        result = np.asarray(self._dct_pchip(np.log(energy))) * ANGSTROM2_TO_M2
        return _scalar_or_array(result)

    def dct_sum_cross_section_derivative_m2_per_eV(
        self, relative_energy_eV,
    ) -> float | np.ndarray:
        energy = self._checked_relative_energy(relative_energy_eV)
        # The interpolant coordinate is log(E), hence the chain-rule 1/E.
        derivative = (
            np.asarray(self._dct_pchip.derivative()(np.log(energy)))
            / energy
            * ANGSTROM2_TO_M2
        )
        return _scalar_or_array(derivative)

    def reactive_destruction_cross_section(
        self, relative_energy_eV,
    ) -> ReactiveCrossSectionBand:
        energy = self._checked_relative_energy(relative_energy_eV)
        dct_A2 = np.asarray(self.dct_sum_cross_section_m2(energy)) / (
            ANGSTROM2_TO_M2)
        central_A2 = self.cid_sum_cross_section_A2 + dct_A2
        lower_A2 = (
            self.cid_sum_cross_section_A2
            * (1.0 - self.cid_source_relative_uncertainty)
            + np.maximum(
                0.0,
                dct_A2 * (1.0 - self.dct_source_relative_uncertainty)
                - self.digitization_cross_section_A2_bound,
            )
        )
        upper_A2 = (
            self.cid_sum_cross_section_A2
            * (1.0 + self.cid_source_relative_uncertainty)
            + dct_A2 * (1.0 + self.dct_source_relative_uncertainty)
            + self.digitization_cross_section_A2_bound
        )
        return ReactiveCrossSectionBand(
            central_m2=_scalar_or_array(central_A2 * ANGSTROM2_TO_M2),
            lower_m2=_scalar_or_array(lower_A2 * ANGSTROM2_TO_M2),
            upper_m2=_scalar_or_array(upper_A2 * ANGSTROM2_TO_M2),
        )

    def slab_sensitivity(
        self,
        *,
        laboratory_energy_eV: float,
        chf3_number_density_m3: float,
        path_length_m: float,
        feed_fraction_used_as_density_proxy: bool,
    ) -> ReactiveSlabSensitivity:
        energy_lab = float(laboratory_energy_eV)
        density = float(chf3_number_density_m3)
        path = float(path_length_m)
        if (
            not math.isfinite(density)
            or density <= 0.0
            or not math.isfinite(path)
            or path <= 0.0
        ):
            raise ValueError("positive CHF3 density and path length required")
        energy_rel = float(self.relative_energy_from_laboratory_eV(energy_lab))
        band = self.reactive_destruction_cross_section(energy_rel)
        central = density * path * float(band.central_m2)
        lower = density * path * float(band.lower_m2)
        upper = density * path * float(band.upper_m2)
        return ReactiveSlabSensitivity(
            laboratory_energy_eV=energy_lab,
            relative_collision_energy_eV=energy_rel,
            path_length_m=path,
            chf3_number_density_m3=density,
            cross_section=band,
            optical_depth_central=central,
            optical_depth_lower=lower,
            optical_depth_upper=upper,
            destruction_probability_central=-math.expm1(-central),
            destruction_probability_lower=-math.expm1(-lower),
            destruction_probability_upper=-math.expm1(-upper),
            feed_fraction_used_as_density_proxy=bool(
                feed_fraction_used_as_density_proxy),
        )


def load_peko_2002_cf3_chf3_reactive_collision_model(
    path: str | Path = DEFAULT_PEKO_DCT_CSV,
) -> Peko2002CF3CHF3ReactiveCollisionModel:
    """Load the checksum-pinned Figure-2 reconstruction."""

    source = Path(path)
    raw = source.read_bytes()
    if sha256(raw).hexdigest() != PEKO_DCT_CSV_SHA256:
        raise ValueError("Peko Figure-2 DCT table checksum mismatch")
    parsed = list(csv.DictReader(raw.decode("utf-8").splitlines()))
    if len(parsed) != 21:
        raise ValueError("unexpected Peko Figure-2 point count")
    digitization_bounds = {
        float(row["digitization_cross_section_A2_bound"])
        for row in parsed
    }
    uncertainties = {
        float(row["source_measurement_relative_uncertainty"])
        for row in parsed
    }
    if len(digitization_bounds) != 1 or uncertainties != {0.25}:
        raise ValueError("Peko digitization metadata changed")
    return Peko2002CF3CHF3ReactiveCollisionModel(
        relative_energy_eV=np.asarray([
            float(row["collision_energy_eV"]) for row in parsed]),
        dct_sum_cross_section_A2=np.asarray([
            float(row["dct_sum_cross_section_A2"]) for row in parsed]),
        digitization_cross_section_A2_bound=digitization_bounds.pop(),
    )
