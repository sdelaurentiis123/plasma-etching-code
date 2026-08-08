"""Neutral diffusion to a partially reactive cylindrical wall.

The volume-loss calculation solves the fundamental separable mode of the
diffusion equation with a Robin wall boundary.  It therefore retains radial
and axial geometry rather than replacing transport with a ballistic ``A/V``
loss.  Chantry's resistance-sum approximation is reported only as a diagnostic;
the returned loss frequency comes from the two transcendental eigenvalues.

The wall reaction velocity uses the Motz--Wise/Chantry correction

    h = gamma * vbar / (2 * (2 - gamma)),

so the extrapolation length is ``D/h``.  This continuum boundary is distinct
from the local incident molecular flux ``n*vbar/4``: its density is the
diffusion-field boundary value after partial reflections have been reduced to
a Robin condition.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros

from .geometry import CylindricalReactor

_FIRST_J0_ZERO = float(jn_zeros(0, 1)[0])
_BOLTZMANN_J_K = 1.380649e-23
_STANDARD_PRESSURE_PA = 101325.0
_CM2_TO_M2 = 1.0e-4

NEUFELD_1972_OMEGA_11_COEFFICIENTS = (
    1.06036,
    0.15610,
    0.19300,
    0.47635,
    1.03587,
    1.52996,
    1.76474,
    3.89411,
)
NEUFELD_1972_REDUCED_TEMPERATURE_DOMAIN = (0.3, 100.0)
NEUFELD_1972_OMEGA_11_MAXIMUM_RELATIVE_FIT_ERROR = 0.0011
NEUTRAL_DIFFUSIVITY_EVIDENCE_KINDS = frozenset({
    "measured",
    "validated_model",
    "published_model",
    "assumed",
    "sensitivity",
})
_PREDICTIVE_DIFFUSIVITY_EVIDENCE = frozenset({
    "measured",
    "validated_model",
})


def neufeld_1972_lennard_jones_omega_11(
    reduced_temperature: float,
) -> float:
    """Evaluate Neufeld--Janzen--Aziz Table-I ``Omega(1,1)*``.

    The coefficient row and four-term equation are taken from the original
    article, DOI 10.1063/1.1678363.  The authors fit the Lennard--Jones
    (12-6) collision integral over ``0.3 <= T* <= 100``; Table II reports a
    maximum deviation of 0.11% for this row.  No extrapolation is allowed.
    """
    temperature = float(reduced_temperature)
    lower, upper = NEUFELD_1972_REDUCED_TEMPERATURE_DOMAIN
    if (
        not np.isfinite(temperature)
        or not lower <= temperature <= upper
    ):
        raise ValueError(
            "reduced temperature is outside the Neufeld 1972 domain")
    a, b, c, d, e, f, g, h = NEUFELD_1972_OMEGA_11_COEFFICIENTS
    return float(
        a / temperature ** b
        + c / np.exp(d * temperature)
        + e / np.exp(f * temperature)
        + g / np.exp(h * temperature)
    )


@dataclass(frozen=True)
class ChapmanEnskogBinaryDiffusivity:
    """First-approximation binary diffusion for a Lennard--Jones pair.

    The evaluated SI state uses the conventional Chapman--Enskog expression
    printed by Malyshev et al. in cgs engineering units,

    ``D[cm2/s] = 0.002628 sqrt(T^3 (Ma+Mb)/(2 Ma Mb))``
    ``            / (p[atm] sigma_ab[A]^2 Omega(1,1)*)``.

    Arithmetic ``sigma`` and geometric ``epsilon`` mixing rules are explicit.
    A source-declared multiplicative correction may be retained, but it is
    provenance rather than an adjustable calibration knob.
    """

    species_a: str
    species_b: str
    molar_mass_a_g_mol: float
    molar_mass_b_g_mol: float
    sigma_a_angstrom: float
    sigma_b_angstrom: float
    epsilon_a_over_k_K: float
    epsilon_b_over_k_K: float
    source_correction_factor: float
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        values = np.asarray([
            self.molar_mass_a_g_mol,
            self.molar_mass_b_g_mol,
            self.sigma_a_angstrom,
            self.sigma_b_angstrom,
            self.epsilon_a_over_k_K,
            self.epsilon_b_over_k_K,
            self.source_correction_factor,
        ], dtype=float)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not str(self.species_a).strip()
            or not str(self.species_b).strip()
            or self.species_a == self.species_b
            or np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
            or self.evidence_kind not in NEUTRAL_DIFFUSIVITY_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid Chapman-Enskog binary diffusivity")
        for name, value in zip(
            (
                "molar_mass_a_g_mol",
                "molar_mass_b_g_mol",
                "sigma_a_angstrom",
                "sigma_b_angstrom",
                "epsilon_a_over_k_K",
                "epsilon_b_over_k_K",
                "source_correction_factor",
            ),
            values,
        ):
            object.__setattr__(self, name, float(value))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def mixed_sigma_angstrom(self) -> float:
        return float(0.5 * (
            self.sigma_a_angstrom + self.sigma_b_angstrom))

    @property
    def mixed_epsilon_over_k_K(self) -> float:
        return float(np.sqrt(
            self.epsilon_a_over_k_K * self.epsilon_b_over_k_K))

    @property
    def valid_temperature_K(self) -> tuple[float, float]:
        lower, upper = NEUFELD_1972_REDUCED_TEMPERATURE_DOMAIN
        epsilon = self.mixed_epsilon_over_k_K
        return float(lower * epsilon), float(upper * epsilon)

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_DIFFUSIVITY_EVIDENCE
            and self.relative_uncertainty is not None
        )

    def collision_integral(self, gas_temperature_K: float) -> float:
        temperature = float(gas_temperature_K)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("gas temperature must be positive and finite")
        return neufeld_1972_lennard_jones_omega_11(
            temperature / self.mixed_epsilon_over_k_K)

    def diffusivity_cm2_s_at_pressure(
        self,
        *,
        gas_temperature_K: float,
        pressure_Pa: float,
    ) -> float:
        temperature = float(gas_temperature_K)
        pressure = float(pressure_Pa)
        if (
            not np.isfinite(temperature)
            or temperature <= 0.0
            or not np.isfinite(pressure)
            or pressure <= 0.0
        ):
            raise ValueError("temperature and pressure must be positive")
        omega = self.collision_integral(temperature)
        mass_factor = (
            temperature ** 3
            * (self.molar_mass_a_g_mol + self.molar_mass_b_g_mol)
            / (2.0 * self.molar_mass_a_g_mol * self.molar_mass_b_g_mol)
        )
        raw = (
            0.002628 * np.sqrt(mass_factor)
            / (
                (pressure / _STANDARD_PRESSURE_PA)
                * self.mixed_sigma_angstrom ** 2
                * omega
            )
        )
        return float(self.source_correction_factor * raw)

    def evaluate(
        self,
        *,
        total_neutral_density_m3: float,
        gas_temperature_K: float,
    ) -> "NeutralDiffusivityState":
        density = float(total_neutral_density_m3)
        temperature = float(gas_temperature_K)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError(
                "total neutral density must be positive and finite")
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("gas temperature must be positive and finite")
        pressure = density * _BOLTZMANN_J_K * temperature
        diffusivity_cm2_s = self.diffusivity_cm2_s_at_pressure(
            gas_temperature_K=temperature,
            pressure_Pa=pressure,
        )
        reduced_temperature = temperature / self.mixed_epsilon_over_k_K
        omega = self.collision_integral(temperature)
        return NeutralDiffusivityState(
            diffusivity_m2_s=diffusivity_cm2_s * _CM2_TO_M2,
            total_neutral_density_m3=density,
            gas_temperature_K=temperature,
            source=self.source,
            evidence_kind=self.evidence_kind,
            relative_uncertainty=self.relative_uncertainty,
            provenance={
                **self.provenance,
                "species_pair": (self.species_a, self.species_b),
                "molar_masses_g_mol": (
                    self.molar_mass_a_g_mol,
                    self.molar_mass_b_g_mol,
                ),
                "sigma_angstrom": (
                    self.sigma_a_angstrom,
                    self.sigma_b_angstrom,
                ),
                "epsilon_over_k_K": (
                    self.epsilon_a_over_k_K,
                    self.epsilon_b_over_k_K,
                ),
                "mixed_sigma_angstrom": self.mixed_sigma_angstrom,
                "mixed_epsilon_over_k_K": self.mixed_epsilon_over_k_K,
                "reduced_temperature": reduced_temperature,
                "omega_11": omega,
                "pressure_Pa": pressure,
                "diffusivity_cm2_s": diffusivity_cm2_s,
                "source_correction_factor": self.source_correction_factor,
                "neufeld_maximum_relative_fit_error": (
                    NEUFELD_1972_OMEGA_11_MAXIMUM_RELATIVE_FIT_ERROR),
                "valid_temperature_K": self.valid_temperature_K,
            },
        )


@dataclass(frozen=True)
class ReducedNeutralDiffusivity:
    """Density-reduced binary diffusivity with a strict temperature domain."""

    reduced_diffusivity_m_inv_s: float
    reference_temperature_K: float
    valid_temperature_K: tuple[float, float]
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        values = np.asarray([
            self.reduced_diffusivity_m_inv_s,
            self.reference_temperature_K,
        ], dtype=float)
        try:
            lower, upper = (
                float(value) for value in self.valid_temperature_K)
        except (TypeError, ValueError):
            raise ValueError(
                "temperature domain must contain two numbers") from None
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not np.isfinite(lower)
            or not np.isfinite(upper)
            or lower <= 0.0
            or upper < lower
            or not lower <= self.reference_temperature_K <= upper
            or not str(self.source).strip()
            or self.evidence_kind not in NEUTRAL_DIFFUSIVITY_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid reduced neutral diffusivity")
        object.__setattr__(
            self,
            "reduced_diffusivity_m_inv_s",
            float(self.reduced_diffusivity_m_inv_s),
        )
        object.__setattr__(
            self, "reference_temperature_K",
            float(self.reference_temperature_K))
        object.__setattr__(self, "valid_temperature_K", (lower, upper))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_DIFFUSIVITY_EVIDENCE
            and self.relative_uncertainty is not None
        )

    def evaluate(
        self,
        *,
        total_neutral_density_m3: float,
        gas_temperature_K: float,
    ) -> "NeutralDiffusivityState":
        density = float(total_neutral_density_m3)
        temperature = float(gas_temperature_K)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError(
                "total neutral density must be positive and finite")
        if (
            not np.isfinite(temperature)
            or not self.valid_temperature_K[0]
            <= temperature
            <= self.valid_temperature_K[1]
        ):
            raise ValueError(
                "gas temperature is outside the diffusivity evidence domain")
        return NeutralDiffusivityState(
            diffusivity_m2_s=(
                self.reduced_diffusivity_m_inv_s / density),
            total_neutral_density_m3=density,
            gas_temperature_K=temperature,
            source=self.source,
            evidence_kind=self.evidence_kind,
            relative_uncertainty=self.relative_uncertainty,
            provenance={
                **self.provenance,
                "reduced_diffusivity_m_inv_s":
                    self.reduced_diffusivity_m_inv_s,
                "reference_temperature_K": self.reference_temperature_K,
                "valid_temperature_K": self.valid_temperature_K,
            },
        )


@dataclass(frozen=True)
class NeutralDiffusivityState:
    """Evaluated bulk diffusivity and its evidence chain."""

    diffusivity_m2_s: float
    total_neutral_density_m3: float
    gas_temperature_K: float
    source: str
    evidence_kind: str
    relative_uncertainty: float | None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        values = np.asarray([
            self.diffusivity_m2_s,
            self.total_neutral_density_m3,
            self.gas_temperature_K,
        ], dtype=float)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
            or self.evidence_kind not in NEUTRAL_DIFFUSIVITY_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid neutral diffusivity state")
        for name, value in zip(
            (
                "diffusivity_m2_s",
                "total_neutral_density_m3",
                "gas_temperature_K",
            ),
            values,
        ):
            object.__setattr__(self, name, float(value))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_DIFFUSIVITY_EVIDENCE
            and self.relative_uncertainty is not None
        )


@dataclass(frozen=True)
class CylindricalNeutralWallLoss:
    """Fundamental neutral-loss mode for one cylindrical reactor state."""

    diffusivity_m2_s: float
    mean_thermal_speed_m_s: float
    wall_reaction_probability: float
    extrapolation_length_m: float
    radial_dimensionless_root: float
    axial_dimensionless_root: float
    radial_wavenumber_m_inv: float
    axial_wavenumber_m_inv: float
    exact_loss_frequency_s_inv: float
    chantry_loss_frequency_s_inv: float
    surface_limit_frequency_s_inv: float
    absorbing_wall_frequency_s_inv: float
    radial_eigen_residual: float
    axial_eigen_residual: float

    def __post_init__(self):
        finite_nonnegative = np.asarray([
            self.diffusivity_m2_s,
            self.mean_thermal_speed_m_s,
            self.wall_reaction_probability,
            self.radial_dimensionless_root,
            self.axial_dimensionless_root,
            self.radial_wavenumber_m_inv,
            self.axial_wavenumber_m_inv,
            self.exact_loss_frequency_s_inv,
            self.chantry_loss_frequency_s_inv,
            self.surface_limit_frequency_s_inv,
            self.absorbing_wall_frequency_s_inv,
            self.radial_eigen_residual,
            self.axial_eigen_residual,
        ], dtype=float)
        if (
            np.any(~np.isfinite(finite_nonnegative))
            or np.any(finite_nonnegative < 0.0)
            or self.diffusivity_m2_s <= 0.0
            or self.mean_thermal_speed_m_s <= 0.0
            or self.wall_reaction_probability > 1.0
            or (
                not np.isfinite(self.extrapolation_length_m)
                and not (
                    np.isinf(self.extrapolation_length_m)
                    and self.wall_reaction_probability == 0.0
                )
            )
            or (
                np.isfinite(self.extrapolation_length_m)
                and self.extrapolation_length_m <= 0.0
            )
        ):
            raise ValueError("invalid cylindrical neutral wall-loss state")
        for name in (
            "diffusivity_m2_s",
            "mean_thermal_speed_m_s",
            "wall_reaction_probability",
            "extrapolation_length_m",
            "radial_dimensionless_root",
            "axial_dimensionless_root",
            "radial_wavenumber_m_inv",
            "axial_wavenumber_m_inv",
            "exact_loss_frequency_s_inv",
            "chantry_loss_frequency_s_inv",
            "surface_limit_frequency_s_inv",
            "absorbing_wall_frequency_s_inv",
            "radial_eigen_residual",
            "axial_eigen_residual",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))

    @property
    def maximum_eigen_residual(self) -> float:
        return float(max(
            self.radial_eigen_residual,
            self.axial_eigen_residual,
        ))

    @property
    def numerical_closure_passes(self) -> bool:
        return self.maximum_eigen_residual <= 1.0e-10

    @property
    def residence_time_s(self) -> float:
        if self.exact_loss_frequency_s_inv == 0.0:
            return float(np.inf)
        return float(1.0 / self.exact_loss_frequency_s_inv)


def solve_cylindrical_neutral_wall_loss(
    *,
    geometry: CylindricalReactor,
    diffusivity_m2_s: float,
    mean_thermal_speed_m_s: float,
    wall_reaction_probability: float,
) -> CylindricalNeutralWallLoss:
    """Solve the exact fundamental cylindrical Robin eigenmode.

    The cylinder extends over ``-L/2 <= z <= L/2``.  With extrapolation
    length ``lambda = D/h``, the dimensionless roots satisfy

    ``(lambda/R) x J1(x) = J0(x)`` and
    ``(2 lambda/L) y tan(y) = 1``.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("cylindrical reactor geometry is required")
    values = np.asarray([
        diffusivity_m2_s,
        mean_thermal_speed_m_s,
        wall_reaction_probability,
    ], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or diffusivity_m2_s <= 0.0
        or mean_thermal_speed_m_s <= 0.0
        or not 0.0 <= wall_reaction_probability <= 1.0
    ):
        raise ValueError("invalid neutral wall-transport inputs")

    diffusivity = float(diffusivity_m2_s)
    mean_speed = float(mean_thermal_speed_m_s)
    probability = float(wall_reaction_probability)
    absorbing_inverse_length_squared = (
        (_FIRST_J0_ZERO / geometry.radius_m) ** 2
        + (np.pi / geometry.length_m) ** 2
    )
    absorbing_frequency = float(
        diffusivity * absorbing_inverse_length_squared)

    if probability == 0.0:
        return CylindricalNeutralWallLoss(
            diffusivity_m2_s=diffusivity,
            mean_thermal_speed_m_s=mean_speed,
            wall_reaction_probability=probability,
            extrapolation_length_m=np.inf,
            radial_dimensionless_root=0.0,
            axial_dimensionless_root=0.0,
            radial_wavenumber_m_inv=0.0,
            axial_wavenumber_m_inv=0.0,
            exact_loss_frequency_s_inv=0.0,
            chantry_loss_frequency_s_inv=0.0,
            surface_limit_frequency_s_inv=0.0,
            absorbing_wall_frequency_s_inv=absorbing_frequency,
            radial_eigen_residual=0.0,
            axial_eigen_residual=0.0,
        )

    reaction_velocity = (
        probability * mean_speed / (2.0 * (2.0 - probability)))
    extrapolation_length = diffusivity / reaction_velocity
    radial_root, axial_root = _fundamental_cylindrical_roots(
        geometry.radius_m,
        geometry.length_m,
        extrapolation_length,
    )
    radial_wavenumber = radial_root / geometry.radius_m
    axial_wavenumber = 2.0 * axial_root / geometry.length_m
    exact_frequency = diffusivity * (
        radial_wavenumber ** 2 + axial_wavenumber ** 2)

    volume_to_area = geometry.volume_m3 / geometry.physical_area_m2
    absorbing_diffusion_length_squared = (
        1.0 / absorbing_inverse_length_squared)
    chantry_effective_length_squared = (
        absorbing_diffusion_length_squared
        + volume_to_area * extrapolation_length
    )
    chantry_frequency = (
        diffusivity / chantry_effective_length_squared)
    surface_frequency = (
        reaction_velocity
        * geometry.physical_area_m2
        / geometry.volume_m3
    )

    radial_ratio = extrapolation_length / geometry.radius_m
    axial_ratio = 2.0 * extrapolation_length / geometry.length_m
    radial_residual = abs(
        radial_ratio * radial_root * j1(radial_root)
        - j0(radial_root)
    )
    axial_residual = abs(
        axial_ratio * axial_root * np.tan(axial_root) - 1.0)

    return CylindricalNeutralWallLoss(
        diffusivity_m2_s=diffusivity,
        mean_thermal_speed_m_s=mean_speed,
        wall_reaction_probability=probability,
        extrapolation_length_m=extrapolation_length,
        radial_dimensionless_root=radial_root,
        axial_dimensionless_root=axial_root,
        radial_wavenumber_m_inv=radial_wavenumber,
        axial_wavenumber_m_inv=axial_wavenumber,
        exact_loss_frequency_s_inv=exact_frequency,
        chantry_loss_frequency_s_inv=chantry_frequency,
        surface_limit_frequency_s_inv=surface_frequency,
        absorbing_wall_frequency_s_inv=absorbing_frequency,
        radial_eigen_residual=radial_residual,
        axial_eigen_residual=axial_residual,
    )


@lru_cache(maxsize=4096)
def _fundamental_cylindrical_roots(
    radius_m: float,
    length_m: float,
    extrapolation_length_m: float,
) -> tuple[float, float]:
    radial_ratio = extrapolation_length_m / radius_m
    axial_ratio = 2.0 * extrapolation_length_m / length_m

    radial_root = brentq(
        lambda root: (
            radial_ratio * root * j1(root) - j0(root)
        ),
        0.0,
        _FIRST_J0_ZERO,
        xtol=1.0e-14,
        rtol=1.0e-14,
    )
    axial_root = brentq(
        lambda root: axial_ratio * root * np.tan(root) - 1.0,
        0.0,
        np.nextafter(0.5 * np.pi, 0.0),
        xtol=1.0e-14,
        rtol=1.0e-14,
    )
    return float(radial_root), float(axial_root)
