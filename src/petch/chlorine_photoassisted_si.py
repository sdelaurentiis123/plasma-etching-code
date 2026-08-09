"""Measured shortwave-VUV photo-assisted etching response of chlorinated Si."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.interpolate import PchipInterpolator


@dataclass(frozen=True)
class Du2022ShortwavePhotoEtchYield:
    """Absolute external-VUV yield interval measured by Du et al. (2022).

    The tandem source was dominated by 104.82 and 106.67 nm photons.  The
    reported 90--244 Si atoms/photon is large because the photon initiates a
    chlorine-mediated reaction chain; it is not a one-photon/one-atom law.
    Application at 139 nm or to an RF-biased surface is not source-supported.
    """

    silicon_atoms_per_photon: float
    silicon_atom_density_m3: float = 5.0e28
    minimum_supported_wavelength_nm: float = 104.82
    maximum_supported_wavelength_nm: float = 106.67

    def __post_init__(self):
        values = (
            float(self.silicon_atoms_per_photon),
            float(self.silicon_atom_density_m3),
            float(self.minimum_supported_wavelength_nm),
            float(self.maximum_supported_wavelength_nm),
        )
        if (
            any(not math.isfinite(value) for value in values)
            or not 90.0 <= values[0] <= 244.0
            or values[1] <= 0.0
            or values[2] <= 0.0
            or values[3] < values[2]
        ):
            raise ValueError("invalid Du-2022 shortwave photo-etch card")
        for name, value in zip((
            "silicon_atoms_per_photon", "silicon_atom_density_m3",
            "minimum_supported_wavelength_nm",
            "maximum_supported_wavelength_nm",
        ), values):
            object.__setattr__(self, name, value)

    @classmethod
    def measured_lower_yield(cls) -> "Du2022ShortwavePhotoEtchYield":
        return cls(silicon_atoms_per_photon=90.0)

    @classmethod
    def measured_upper_yield(cls) -> "Du2022ShortwavePhotoEtchYield":
        return cls(silicon_atoms_per_photon=244.0)

    def etch_velocity_m_s(
        self,
        photon_flux_m2_s: float,
        *,
        photon_wavelength_nm: float,
    ) -> float:
        flux = float(photon_flux_m2_s)
        wavelength = float(photon_wavelength_nm)
        if not math.isfinite(flux) or flux < 0.0 or not math.isfinite(wavelength):
            raise ValueError("invalid photon boundary")
        if not (
            self.minimum_supported_wavelength_nm
            <= wavelength
            <= self.maximum_supported_wavelength_nm
        ):
            raise ValueError(
                "Du-2022 yield is wavelength-specific and cannot be transferred"
            )
        return float(
            flux * self.silicon_atoms_per_photon / self.silicon_atom_density_m3
        )

    def required_photon_flux_m2_s(
        self,
        etch_velocity_m_s: float,
        *,
        photon_wavelength_nm: float,
    ) -> float:
        velocity = float(etch_velocity_m_s)
        if not math.isfinite(velocity) or velocity < 0.0:
            raise ValueError("etch velocity must be finite and nonnegative")
        # Exercise the same wavelength-domain refusal in both directions.
        self.etch_velocity_m_s(0.0, photon_wavelength_nm=photon_wavelength_nm)
        return float(
            velocity * self.silicon_atom_density_m3
            / self.silicon_atoms_per_photon
        )

    @property
    def supports_prediction(self) -> bool:
        return True


@dataclass(frozen=True)
class Hirsch2020PulsedDCAntiSynergySensitivity:
    """Digitized Figure-8 PAE suppression curve for its exact pulsed-DC case.

    Hirsch et al. constructed the red curve so corrected IAE yield is nearly
    independent of Cl/ion ratio.  It is useful evidence for anti-synergy and a
    smooth sensitivity operator, but it is not a measured universal damage law
    and cannot be transferred to RF bias as a predictive closure.
    """

    duty_cycle_percent: tuple[float, ...] = (
        0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0,
        50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0,
    )
    relative_pae_yield: tuple[float, ...] = (
        1.0, 0.981081081081081, 0.9662162162162162,
        0.9108108108108108, 0.831081081081081,
        0.7405405405405405, 0.6594594594594595,
        0.581081081081081, 0.5081081081081081,
        0.43243243243243246, 0.3472972972972973,
        0.2648648648648649, 0.1918918918918919,
        0.13783783783783785, 0.0972972972972973,
        0.062162162162162166, 0.03783783783783784,
        0.021621621621621623, 0.01891891891891892,
    )

    def __post_init__(self):
        duty = np.asarray(self.duty_cycle_percent, dtype=float)
        relative = np.asarray(self.relative_pae_yield, dtype=float)
        if (
            duty.ndim != 1
            or duty.size < 3
            or duty.shape != relative.shape
            or np.any(~np.isfinite(duty))
            or np.any(~np.isfinite(relative))
            or np.any(np.diff(duty) <= 0.0)
            or duty[0] != 0.0
            or duty[-1] > 100.0
            or np.any((relative < 0.0) | (relative > 1.0))
            or np.any(np.diff(relative) > 0.0)
        ):
            raise ValueError("invalid Hirsch Figure-8 sensitivity")
        object.__setattr__(self, "duty_cycle_percent", tuple(duty.tolist()))
        object.__setattr__(self, "relative_pae_yield", tuple(relative.tolist()))

    def relative_yield(self, duty_cycle_percent):
        duty = np.asarray(duty_cycle_percent, dtype=float)
        support = np.asarray(self.duty_cycle_percent)
        if (
            np.any(~np.isfinite(duty))
            or np.any(duty < support[0])
            or np.any(duty > support[-1])
        ):
            raise ValueError("Hirsch Figure-8 duty is outside digitized support")
        interpolator = PchipInterpolator(
            support,
            np.asarray(self.relative_pae_yield),
            extrapolate=False,
        )
        result = np.asarray(interpolator(duty), dtype=float)
        if result.ndim == 0:
            return float(result)
        return result

    @property
    def supports_prediction(self) -> bool:
        return False
