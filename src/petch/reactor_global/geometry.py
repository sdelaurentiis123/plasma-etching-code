"""Cylindrical global-model geometry and sheath-edge factors."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import j1


@dataclass(frozen=True)
class ElectropositiveEdgeFactors:
    """Sheath-edge to volume-average density ratios."""

    axial: float
    radial: float

    def __post_init__(self):
        values = np.asarray([self.axial, self.radial], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0) or np.any(values > 1.0):
            raise ValueError("edge factors must lie in (0, 1]")
        object.__setattr__(self, "axial", float(self.axial))
        object.__setattr__(self, "radial", float(self.radial))


@dataclass(frozen=True)
class ElectronegativeEdgeFactors(ElectropositiveEdgeFactors):
    """Lee--Lieberman sheath-edge ratios with negative-ion correction.

    ``electronegativity`` is the volume-average ratio ``n_minus / n_e`` and
    ``electron_to_ion_temperature_ratio`` is ``T_e / T_i``.  These names keep
    the temperature ratio distinct from a surface-reaction probability, which
    is also conventionally written as gamma.
    """

    electronegativity: float
    electron_to_ion_temperature_ratio: float
    electronegative_correction: float

    def __post_init__(self):
        super().__post_init__()
        alpha = float(self.electronegativity)
        temperature_ratio = float(self.electron_to_ion_temperature_ratio)
        correction = float(self.electronegative_correction)
        if not np.isfinite(alpha) or alpha < 0.0:
            raise ValueError("electronegativity must be nonnegative and finite")
        if not np.isfinite(temperature_ratio) or temperature_ratio <= 0.0:
            raise ValueError(
                "electron-to-ion temperature ratio must be positive and finite")
        expected = (1.0 + 3.0 * alpha / temperature_ratio) / (1.0 + alpha)
        if not np.isclose(correction, expected, rtol=2e-15, atol=0.0):
            raise ValueError("electronegative correction is inconsistent")
        object.__setattr__(self, "electronegativity", alpha)
        object.__setattr__(
            self, "electron_to_ion_temperature_ratio", temperature_ratio)
        object.__setattr__(self, "electronegative_correction", correction)


@dataclass(frozen=True)
class CylindricalReactor:
    """Right circular cylinder used by a volume-averaged reactor model."""

    radius_m: float
    length_m: float

    def __post_init__(self):
        values = np.asarray([self.radius_m, self.length_m], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("reactor dimensions must be positive and finite")
        object.__setattr__(self, "radius_m", float(self.radius_m))
        object.__setattr__(self, "length_m", float(self.length_m))

    @property
    def volume_m3(self) -> float:
        return float(np.pi * self.radius_m ** 2 * self.length_m)

    @property
    def physical_area_m2(self) -> float:
        return float(
            2.0 * np.pi * self.radius_m ** 2
            + 2.0 * np.pi * self.radius_m * self.length_m)

    @property
    def diffusion_length_m(self) -> float:
        inverse_square = (
            (np.pi / self.length_m) ** 2
            + (2.405 / self.radius_m) ** 2
        )
        return float(1.0 / np.sqrt(inverse_square))

    def effective_loss_area_m2(
            self, edge_factors: ElectropositiveEdgeFactors) -> float:
        if not isinstance(edge_factors, ElectropositiveEdgeFactors):
            raise TypeError("sheath-edge factors are required")
        return float(
            edge_factors.axial * 2.0 * np.pi * self.radius_m ** 2
            + edge_factors.radial
            * 2.0 * np.pi * self.radius_m * self.length_m
        )

    def electropositive_edge_factors(
            self, *, ion_mean_free_path_m: float,
            bohm_speed_m_s: float | None = None,
            ambipolar_diffusion_m2_s: float | None = None,
            include_high_pressure_diffusion: bool = True,
            ) -> ElectropositiveEdgeFactors:
        """Evaluate Lee--Lieberman Eqs. 13--14 for ``alpha=0``.

        The high-pressure terms require both Bohm speed and ambipolar
        diffusivity. They may be omitted together only as an explicit
        low/intermediate-pressure approximation.
        """
        mean_free_path = float(ion_mean_free_path_m)
        if not np.isfinite(mean_free_path) or mean_free_path <= 0.0:
            raise ValueError("ion mean free path must be positive and finite")
        if include_high_pressure_diffusion:
            if bohm_speed_m_s is None or ambipolar_diffusion_m2_s is None:
                raise ValueError(
                    "high-pressure edge factors require Bohm speed and diffusivity")
            bohm_speed = float(bohm_speed_m_s)
            diffusivity = float(ambipolar_diffusion_m2_s)
            if (
                not np.isfinite(bohm_speed)
                or bohm_speed <= 0.0
                or not np.isfinite(diffusivity)
                or diffusivity <= 0.0
            ):
                raise ValueError("invalid Bohm speed or ambipolar diffusivity")
            axial_diffusion = (
                0.86 * self.length_m * bohm_speed
                / (np.pi * diffusivity)
            ) ** 2
            radial_diffusion = (
                0.8 * self.radius_m * bohm_speed
                / (2.405 * j1(2.405) * diffusivity)
            ) ** 2
        else:
            if bohm_speed_m_s is not None or ambipolar_diffusion_m2_s is not None:
                raise ValueError(
                    "do not supply diffusion inputs when the term is disabled")
            axial_diffusion = 0.0
            radial_diffusion = 0.0
        axial = 0.86 / np.sqrt(
            3.0
            + self.length_m / (2.0 * mean_free_path)
            + axial_diffusion
        )
        radial = 0.8 / np.sqrt(
            4.0
            + self.radius_m / mean_free_path
            + radial_diffusion
        )
        return ElectropositiveEdgeFactors(axial=axial, radial=radial)

    def electronegative_edge_factors(
            self, *, electronegativity: float,
            electron_to_ion_temperature_ratio: float,
            ion_mean_free_path_m: float,
            bohm_speed_m_s: float | None = None,
            ambipolar_diffusion_m2_s: float | None = None,
            include_high_pressure_diffusion: bool = True,
            ) -> ElectronegativeEdgeFactors:
        """Evaluate Lee--Lieberman Eqs. 13--14.

        The source assumes a spatially uniform electron density, a parabolic
        negative-ion density that vanishes at the sheath edge, and a common
        edge-to-bulk ratio for all positive-ion species.  The returned values
        must therefore be treated as a source-model closure, not a universal
        sheath law.
        """
        alpha = float(electronegativity)
        temperature_ratio = float(electron_to_ion_temperature_ratio)
        if not np.isfinite(alpha) or alpha < 0.0:
            raise ValueError("electronegativity must be nonnegative and finite")
        if not np.isfinite(temperature_ratio) or temperature_ratio <= 0.0:
            raise ValueError(
                "electron-to-ion temperature ratio must be positive and finite")
        electropositive = self.electropositive_edge_factors(
            ion_mean_free_path_m=ion_mean_free_path_m,
            bohm_speed_m_s=bohm_speed_m_s,
            ambipolar_diffusion_m2_s=ambipolar_diffusion_m2_s,
            include_high_pressure_diffusion=include_high_pressure_diffusion,
        )
        correction = (
            1.0 + 3.0 * alpha / temperature_ratio
        ) / (1.0 + alpha)
        return ElectronegativeEdgeFactors(
            axial=correction * electropositive.axial,
            radial=correction * electropositive.radial,
            electronegativity=alpha,
            electron_to_ion_temperature_ratio=temperature_ratio,
            electronegative_correction=correction,
        )
