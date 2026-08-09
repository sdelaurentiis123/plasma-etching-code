"""Deterministic volume-radiation transfer from a cylinder to a wafer disk.

The kernel is the exact line-of-sight cosine law for isotropic photons emitted
uniformly in a right circular cylinder.  Fixed Gauss--Legendre quadrature
replaces photon Monte Carlo, so repeated reactor sweeps are reproducible,
smooth in every physical input, and embarrassingly parallel over conditions.

This module deliberately stops at radiative transfer.  It does not infer an
emission spectrum, resonance-line escape probability, wall reflectivity, or a
silicon photo-etch yield.  Those are separate physical closures.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import (
    c as SPEED_OF_LIGHT_M_S,
    e as ELEMENTARY_CHARGE_C,
    epsilon_0 as VACUUM_PERMITTIVITY_F_M,
    k as BOLTZMANN_CONSTANT_J_K,
    m_e as ELECTRON_MASS_KG,
    physical_constants,
)
from scipy.fft import irfft, next_fast_len, rfft
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.special import exp1, voigt_profile

from .geometry import CylindricalReactor


@dataclass(frozen=True)
class CylinderDiskRadiationResult:
    """Integrated one-pass radiation receipt for one circular wafer."""

    volume_emissivity_m3_s: float
    wafer_radius_m: float
    extinction_coefficient_m_inv: float
    wafer_photon_flux_m2_s: float
    wafer_intercept_probability: float
    geometry_flux_length_m: float
    quadrature_order: int

    def __post_init__(self):
        values = np.asarray((
            self.volume_emissivity_m3_s,
            self.wafer_radius_m,
            self.extinction_coefficient_m_inv,
            self.wafer_photon_flux_m2_s,
            self.wafer_intercept_probability,
            self.geometry_flux_length_m,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.volume_emissivity_m3_s < 0.0
            or self.wafer_radius_m <= 0.0
            or self.extinction_coefficient_m_inv < 0.0
            or self.wafer_photon_flux_m2_s < 0.0
            or not 0.0 <= self.wafer_intercept_probability <= 0.5
            or self.geometry_flux_length_m < 0.0
            or int(self.quadrature_order) != self.quadrature_order
            or self.quadrature_order < 4
        ):
            raise ValueError("invalid cylinder-to-disk radiation result")
        for name, value in zip((
            "volume_emissivity_m3_s", "wafer_radius_m",
            "extinction_coefficient_m_inv", "wafer_photon_flux_m2_s",
            "wafer_intercept_probability", "geometry_flux_length_m",
        ), values):
            object.__setattr__(self, name, float(value))
        object.__setattr__(self, "quadrature_order", int(self.quadrature_order))


@dataclass(frozen=True)
class LamellarFloorRadiationResult:
    """Direct photon receipt at the floor of an absorbing line trench."""

    wafer_plane_photon_flux_m2_s: float
    trench_floor_photon_flux_m2_s: float
    direct_floor_fraction: float
    opening_width_m: float
    optical_path_depth_m: float
    radial_quadrature_order: int
    axial_quadrature_order: int

    def __post_init__(self):
        values = np.asarray((
            self.wafer_plane_photon_flux_m2_s,
            self.trench_floor_photon_flux_m2_s,
            self.direct_floor_fraction,
            self.opening_width_m,
            self.optical_path_depth_m,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.wafer_plane_photon_flux_m2_s < 0.0
            or self.trench_floor_photon_flux_m2_s < 0.0
            or not 0.0 <= self.direct_floor_fraction <= 1.0
            or self.opening_width_m <= 0.0
            or self.optical_path_depth_m < 0.0
            or int(self.radial_quadrature_order) != self.radial_quadrature_order
            or int(self.axial_quadrature_order) != self.axial_quadrature_order
            or min(self.radial_quadrature_order, self.axial_quadrature_order) < 4
        ):
            raise ValueError("invalid lamellar-floor radiation result")
        if self.wafer_plane_photon_flux_m2_s == 0.0:
            if self.trench_floor_photon_flux_m2_s != 0.0:
                raise ValueError("a zero wafer flux requires a zero floor flux")
        elif not np.isclose(
            self.trench_floor_photon_flux_m2_s,
            self.wafer_plane_photon_flux_m2_s * self.direct_floor_fraction,
            rtol=2.0e-14,
            atol=0.0,
        ):
            raise ValueError("lamellar-floor flux ledger does not close")
        for name, value in zip((
            "wafer_plane_photon_flux_m2_s", "trench_floor_photon_flux_m2_s",
            "direct_floor_fraction", "opening_width_m", "optical_path_depth_m",
        ), values):
            object.__setattr__(self, name, float(value))
        object.__setattr__(
            self, "radial_quadrature_order", int(self.radial_quadrature_order))
        object.__setattr__(
            self, "axial_quadrature_order", int(self.axial_quadrature_order))


@dataclass(frozen=True)
class ResonanceLineData:
    """Atomic data required by a homogeneous Voigt escape calculation."""

    wavelength_nm: float
    transition_probability_s_inv: float
    lower_statistical_weight: float
    upper_statistical_weight: float
    absorber_mass_kg: float
    source: str

    def __post_init__(self):
        values = np.asarray((
            self.wavelength_nm,
            self.transition_probability_s_inv,
            self.lower_statistical_weight,
            self.upper_statistical_weight,
            self.absorber_mass_kg,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid resonance-line atomic data")
        for name, value in zip((
            "wavelength_nm", "transition_probability_s_inv",
            "lower_statistical_weight", "upper_statistical_weight",
            "absorber_mass_kg",
        ), values):
            object.__setattr__(self, name, float(value))

    @property
    def wavelength_m(self) -> float:
        return self.wavelength_nm * 1.0e-9

    @property
    def center_frequency_hz(self) -> float:
        return SPEED_OF_LIGHT_M_S / self.wavelength_m

    @property
    def absorption_oscillator_strength(self) -> float:
        # NIST Atomic Spectroscopy compendium relation between A_ki and f_ik.
        return float(
            self.upper_statistical_weight / self.lower_statistical_weight
            * ELECTRON_MASS_KG * SPEED_OF_LIGHT_M_S
            * VACUUM_PERMITTIVITY_F_M * self.wavelength_m ** 2
            * self.transition_probability_s_inv
            / (2.0 * np.pi * ELEMENTARY_CHARGE_C ** 2)
        )

    @property
    def natural_lorentz_hwhm_hz(self) -> float:
        return self.transition_probability_s_inv / (4.0 * np.pi)

    def doppler_standard_deviation_hz(self, gas_temperature_K: float) -> float:
        temperature = float(gas_temperature_K)
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("gas temperature must be finite and positive")
        return float(
            self.center_frequency_hz
            * math.sqrt(
                BOLTZMANN_CONSTANT_J_K * temperature
                / (self.absorber_mass_kg * SPEED_OF_LIGHT_M_S ** 2)
            )
        )


@dataclass(frozen=True)
class CylinderResonanceEscapeResult:
    """One-cycle escape and complete-redistribution trapping receipt."""

    escape_probability: float
    trapping_factor: float
    absorber_density_m3: float
    gas_temperature_K: float
    source_axial_scale_length_m: float
    doppler_standard_deviation_hz: float
    lorentz_hwhm_hz: float
    line_center_absorption_coefficient_m_inv: float
    source_weighted_mean_exit_path_m: float
    line_center_mean_path_optical_depth: float
    frequency_profile_normalization: float
    geometry_quadrature_order: int
    frequency_quadrature_order: int
    complete_frequency_redistribution_assumed: bool = True

    def __post_init__(self):
        values = np.asarray((
            self.escape_probability,
            self.trapping_factor,
            self.absorber_density_m3,
            self.gas_temperature_K,
            self.source_axial_scale_length_m,
            self.doppler_standard_deviation_hz,
            self.lorentz_hwhm_hz,
            self.line_center_absorption_coefficient_m_inv,
            self.source_weighted_mean_exit_path_m,
            self.line_center_mean_path_optical_depth,
            self.frequency_profile_normalization,
        ), dtype=float)
        scale = float(self.source_axial_scale_length_m)
        if (
            np.any(np.isnan(values))
            or not 0.0 < self.escape_probability <= 1.0
            or not np.isclose(
                self.trapping_factor, 1.0 / self.escape_probability,
                rtol=2.0e-14, atol=0.0)
            or self.absorber_density_m3 < 0.0
            or self.gas_temperature_K <= 0.0
            or not (math.isinf(scale) or scale > 0.0)
            or min(
                self.doppler_standard_deviation_hz,
                self.lorentz_hwhm_hz,
            ) <= 0.0
            or min(
                self.line_center_absorption_coefficient_m_inv,
                self.source_weighted_mean_exit_path_m,
                self.line_center_mean_path_optical_depth,
            ) < 0.0
            or not 0.999999 <= self.frequency_profile_normalization <= 1.000001
            or int(self.geometry_quadrature_order)
            != self.geometry_quadrature_order
            or int(self.frequency_quadrature_order)
            != self.frequency_quadrature_order
            or min(
                self.geometry_quadrature_order,
                self.frequency_quadrature_order,
            ) < 4
        ):
            raise ValueError("invalid cylinder resonance-escape result")
        for name, value in zip((
            "escape_probability", "trapping_factor", "absorber_density_m3",
            "gas_temperature_K", "source_axial_scale_length_m",
            "doppler_standard_deviation_hz", "lorentz_hwhm_hz",
            "line_center_absorption_coefficient_m_inv",
            "source_weighted_mean_exit_path_m",
            "line_center_mean_path_optical_depth",
            "frequency_profile_normalization",
        ), values):
            object.__setattr__(self, name, float(value))
        object.__setattr__(
            self, "geometry_quadrature_order",
            int(self.geometry_quadrature_order))
        object.__setattr__(
            self, "frequency_quadrature_order",
            int(self.frequency_quadrature_order))


@dataclass(frozen=True)
class CylinderPartialRedistributionResult:
    """Frequency-state trapping receipt with coherent/redistributed cycles."""

    complete_redistribution_escape_probability: float
    complete_redistribution_trapping_factor: float
    partial_redistribution_trapping_factor: float
    coherent_redistribution_probability: float
    velocity_redistribution_probability: float
    radiative_survival_probability_after_absorption: float
    velocity_changing_collision_frequency_s_inv: float
    quenching_collision_frequency_s_inv: float
    frequency_profile_normalization: float
    coherent_half_range_doppler_standard_deviations: float
    coherent_frequency_grid_points: int
    coherent_grid_points_per_lorentz_hwhm: float
    linear_solver_iterations: int
    linear_solver_relative_residual: float
    geometry_quadrature_order: int
    frequency_quadrature_order: int

    def __post_init__(self):
        values = np.asarray((
            self.complete_redistribution_escape_probability,
            self.complete_redistribution_trapping_factor,
            self.partial_redistribution_trapping_factor,
            self.coherent_redistribution_probability,
            self.velocity_redistribution_probability,
            self.radiative_survival_probability_after_absorption,
            self.velocity_changing_collision_frequency_s_inv,
            self.quenching_collision_frequency_s_inv,
            self.frequency_profile_normalization,
            self.coherent_half_range_doppler_standard_deviations,
            self.coherent_grid_points_per_lorentz_hwhm,
            self.linear_solver_relative_residual,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or not 0.0 < self.complete_redistribution_escape_probability <= 1.0
            or not np.isclose(
                self.complete_redistribution_trapping_factor,
                1.0 / self.complete_redistribution_escape_probability,
                rtol=2.0e-14,
                atol=0.0,
            )
            or self.partial_redistribution_trapping_factor < 1.0
            or not 0.0 <= self.coherent_redistribution_probability <= 1.0
            or not 0.0 <= self.velocity_redistribution_probability <= 1.0
            or not np.isclose(
                self.coherent_redistribution_probability
                + self.velocity_redistribution_probability,
                1.0,
                rtol=2.0e-14,
                atol=2.0e-14,
            )
            or not 0.0 < self.radiative_survival_probability_after_absorption <= 1.0
            or min(
                self.velocity_changing_collision_frequency_s_inv,
                self.quenching_collision_frequency_s_inv,
            ) < 0.0
            or not 0.999999 <= self.frequency_profile_normalization <= 1.000001
            or self.coherent_half_range_doppler_standard_deviations < 32.0
            or int(self.coherent_frequency_grid_points)
            != self.coherent_frequency_grid_points
            or self.coherent_frequency_grid_points < 4096
            or self.coherent_frequency_grid_points % 2
            or self.coherent_grid_points_per_lorentz_hwhm < 4.0
            or int(self.linear_solver_iterations) != self.linear_solver_iterations
            or self.linear_solver_iterations < 0
            or self.linear_solver_relative_residual > 1.0e-7
            or int(self.geometry_quadrature_order)
            != self.geometry_quadrature_order
            or int(self.frequency_quadrature_order)
            != self.frequency_quadrature_order
            or min(
                self.geometry_quadrature_order,
                self.frequency_quadrature_order,
            ) < 4
        ):
            raise ValueError("invalid partial-redistribution result")
        for name, value in zip((
            "complete_redistribution_escape_probability",
            "complete_redistribution_trapping_factor",
            "partial_redistribution_trapping_factor",
            "coherent_redistribution_probability",
            "velocity_redistribution_probability",
            "radiative_survival_probability_after_absorption",
            "velocity_changing_collision_frequency_s_inv",
            "quenching_collision_frequency_s_inv",
            "frequency_profile_normalization",
            "coherent_half_range_doppler_standard_deviations",
            "coherent_grid_points_per_lorentz_hwhm",
            "linear_solver_relative_residual",
        ), values):
            object.__setattr__(self, name, float(value))
        object.__setattr__(
            self, "coherent_frequency_grid_points",
            int(self.coherent_frequency_grid_points))
        object.__setattr__(
            self, "linear_solver_iterations", int(self.linear_solver_iterations))
        object.__setattr__(
            self, "geometry_quadrature_order",
            int(self.geometry_quadrature_order))
        object.__setattr__(
            self, "frequency_quadrature_order",
            int(self.frequency_quadrature_order))


def _mapped_legendre(order: int, upper: float) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(order)
    return 0.5 * upper * (nodes + 1.0), 0.5 * upper * weights


def _mapped_legendre_interval(
    order: int, lower: float, upper: float,
) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(order)
    half_width = 0.5 * (upper - lower)
    return half_width * nodes + 0.5 * (upper + lower), half_width * weights


def _cylinder_exit_path_distribution(
    geometry: CylindricalReactor,
    *,
    source_axial_scale_length_m: float,
    quadrature_order: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic source/direction exit paths and normalized weights."""
    nodes, weights = leggauss(quadrature_order)
    unit = 0.5 * (nodes + 1.0)
    unit_weights = 0.5 * weights
    radius = geometry.radius_m * np.sqrt(unit)
    direction_cosine = nodes
    azimuth = np.pi * (nodes + 1.0)
    scale = float(source_axial_scale_length_m)
    if math.isinf(scale):
        axial = geometry.length_m * unit
    else:
        # Inverse-CDF quadrature for an emitter proportional to
        # exp(-(L-z)/scale), avoiding resolution loss for a thin skin source.
        retained = -math.expm1(-geometry.length_m / scale)
        depth_below_top = -scale * np.log1p(-unit * retained)
        axial = geometry.length_m - depth_below_top

    r, z, mu, phi = np.meshgrid(
        radius, axial, direction_cosine, azimuth, indexing="ij")
    path_weights = np.einsum(
        "i,j,k,l->ijkl",
        unit_weights, unit_weights, unit_weights, unit_weights,
    ).ravel()
    r = r.ravel()
    z = z.ravel()
    mu = mu.ravel()
    phi = phi.ravel()
    transverse = np.sqrt(np.maximum(0.0, 1.0 - mu * mu))
    radial_path = (
        -r * np.cos(phi)
        + np.sqrt(np.maximum(
            0.0,
            geometry.radius_m ** 2 - r * r * np.sin(phi) ** 2,
        ))
    ) / transverse
    axial_path = np.where(
        mu > 0.0,
        (geometry.length_m - z) / mu,
        z / (-mu),
    )
    exit_path = np.minimum(radial_path, axial_path)
    return exit_path, path_weights


def _uniform_cylinder_projected_chords(
    geometry: CylindricalReactor,
    *,
    quadrature_order: int,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Projected surface-ray chord measure for a uniform convex cylinder.

    Re-parameterizing every oriented volume ray by its entry surface removes
    the optically thick boundary layer that converges slowly under direct
    volume quadrature.  For a chord ``C``, the path integral is
    ``(1-exp(-k*C))/k``.  End-cap and side measures retain their physical area
    ratio; common angular constants cancel in the normalized transmission.
    """
    nodes, weights = leggauss(quadrature_order)
    unit = 0.5 * (nodes + 1.0)
    unit_w = 0.5 * weights

    # Both end caps.  Point radius is uniform in disk area; mu is the inward
    # axial direction cosine and phi is relative to the local radial vector.
    radius = geometry.radius_m * np.sqrt(unit)
    mu = unit
    phi = 2.0 * np.pi * unit
    r, axial_cosine, azimuth = np.meshgrid(
        radius, mu, phi, indexing="ij")
    transverse = np.sqrt(np.maximum(0.0, 1.0 - axial_cosine ** 2))
    side_path = (
        -r * np.cos(azimuth)
        + np.sqrt(np.maximum(
            0.0,
            geometry.radius_m ** 2 - r * r * np.sin(azimuth) ** 2,
        ))
    ) / transverse
    opposite_cap_path = geometry.length_m / axial_cosine
    end_chord = np.minimum(side_path, opposite_cap_path).ravel()
    end_weight = (
        2.0 * np.pi * geometry.radius_m ** 2
        * np.einsum(
            "i,j,k->ijk", unit_w, unit_w, unit_w
        )
        * axial_cosine
    ).ravel()

    # Cylindrical side. eta is the inward radial direction cosine; psi spans
    # the tangent plane and partitions horizontal-tangent versus axial motion.
    axial = geometry.length_m * unit
    eta = unit
    psi = 2.0 * np.pi * unit
    z, radial_cosine, tangent_angle = np.meshgrid(
        axial, eta, psi, indexing="ij")
    tangent = np.sqrt(np.maximum(0.0, 1.0 - radial_cosine ** 2))
    horizontal_tangent = tangent * np.cos(tangent_angle)
    axial_direction = tangent * np.sin(tangent_angle)
    horizontal_norm_squared = (
        radial_cosine ** 2 + horizontal_tangent ** 2)
    opposite_side_path = (
        2.0 * geometry.radius_m * radial_cosine
        / horizontal_norm_squared)
    end_path = np.where(
        axial_direction > 0.0,
        (geometry.length_m - z) / axial_direction,
        z / (-axial_direction),
    )
    side_chord = np.minimum(opposite_side_path, end_path).ravel()
    side_weight = (
        2.0 * np.pi * geometry.radius_m * geometry.length_m
        * np.einsum(
            "i,j,k->ijk", unit_w, unit_w, unit_w
        )
        * radial_cosine
    ).ravel()
    return (
        np.concatenate((end_chord, side_chord)),
        np.concatenate((end_weight, side_weight)),
    )


def _uniform_cylinder_transmission(
    geometry: CylindricalReactor,
    absorption_coefficient_m_inv: np.ndarray,
    *,
    quadrature_order: int,
) -> tuple[np.ndarray, float]:
    chord, weight = _uniform_cylinder_projected_chords(
        geometry, quadrature_order=quadrature_order)
    normalization = float(np.dot(weight, chord))
    absorption = np.asarray(absorption_coefficient_m_inv, dtype=float)
    transmission = np.ones_like(absorption)
    positive = absorption > 0.0
    for start in range(0, int(np.count_nonzero(positive)), 64):
        indices = np.flatnonzero(positive)[start:start + 64]
        kappa = absorption[indices, None]
        integrated_path = -np.expm1(-kappa * chord[None, :]) / kappa
        transmission[indices] = integrated_path @ weight / normalization
    mean_exit_path = float(
        np.dot(weight, chord * chord) / (2.0 * normalization))
    return transmission, mean_exit_path


def _voigt_positive_frequency_quadrature(
    *,
    doppler_standard_deviation_hz: float,
    lorentz_hwhm_hz: float,
    quadrature_order: int,
    maximum_doppler_standard_deviations: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Integrate a symmetric Voigt density on fixed doubled intervals."""
    sigma = float(doppler_standard_deviation_hz)
    gamma = float(lorentz_hwhm_hz)
    maximum = float(maximum_doppler_standard_deviations)
    if (
        not all(math.isfinite(item) for item in (sigma, gamma, maximum))
        or min(sigma, gamma) <= 0.0
        or maximum < 32.0
    ):
        raise ValueError("invalid Voigt frequency quadrature")
    nodes, weights = leggauss(quadrature_order)
    bounds = [0.0, 1.0]
    while bounds[-1] < maximum:
        bounds.append(min(maximum, 2.0 * bounds[-1]))
    detuning = []
    integration_weights = []
    for lower, upper in zip(bounds, bounds[1:]):
        half_width = 0.5 * (upper - lower)
        detuning.append(half_width * nodes + 0.5 * (upper + lower))
        integration_weights.append(half_width * weights)
    normalized_detuning = np.concatenate(detuning)
    profile = voigt_profile(normalized_detuning * sigma, sigma, gamma)
    # Symmetry supplies the factor two; dx = sigma d(x/sigma).
    probability_weights = (
        2.0 * sigma * np.concatenate(integration_weights) * profile)
    # At the declared large detuning the Voigt tail is Lorentzian and every
    # ray is optically thin.  Retain its analytic probability instead of
    # truncating a physically important escape channel.
    positive_and_negative_tail = 2.0 * gamma / (np.pi * sigma * maximum)
    normalization = float(
        np.sum(probability_weights) + positive_and_negative_tail)
    return profile, probability_weights, normalization


def deterministic_cylinder_resonance_escape(
    geometry: CylindricalReactor,
    line: ResonanceLineData,
    *,
    absorber_density_m3: float,
    gas_temperature_K: float,
    source_axial_scale_length_m: float = math.inf,
    additional_lorentz_hwhm_hz: float = 0.0,
    geometry_quadrature_order: int = 12,
    frequency_quadrature_order: int = 40,
    maximum_doppler_standard_deviations: float = 2048.0,
) -> CylinderResonanceEscapeResult:
    r"""Evaluate complete-redistribution resonance escape without Monte Carlo.

    The absorber is uniform.  Emission is uniform radially and either uniform
    axially or exponentially localized below the top boundary.  Fixed
    Gauss--Legendre quadrature integrates source position, isotropic direction,
    and the Voigt frequency profile.  The trapping factor is the geometric
    number of radiative cycles, ``1 / P_escape``, under complete frequency
    redistribution and zero quenching.

    This is a differentiable reduced closure away from the cylinder's
    measure-zero path-branch surfaces.  It deliberately exposes source
    localization, collisional broadening, and quenching as separate closures.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("a cylindrical reactor geometry is required")
    if not isinstance(line, ResonanceLineData):
        raise TypeError("resonance-line atomic data are required")
    density = float(absorber_density_m3)
    temperature = float(gas_temperature_K)
    scale = float(source_axial_scale_length_m)
    extra_gamma = float(additional_lorentz_hwhm_hz)
    geometry_order = int(geometry_quadrature_order)
    frequency_order = int(frequency_quadrature_order)
    if (
        not math.isfinite(density)
        or density < 0.0
        or not math.isfinite(temperature)
        or temperature <= 0.0
        or not (math.isinf(scale) or (math.isfinite(scale) and scale > 0.0))
        or not math.isfinite(extra_gamma)
        or extra_gamma < 0.0
        or geometry_order != geometry_quadrature_order
        or frequency_order != frequency_quadrature_order
        or min(geometry_order, frequency_order) < 4
    ):
        raise ValueError("invalid cylinder resonance-escape condition")

    if math.isinf(scale):
        exit_path = path_weights = None
    else:
        exit_path, path_weights = _cylinder_exit_path_distribution(
            geometry,
            source_axial_scale_length_m=scale,
            quadrature_order=geometry_order,
        )
    sigma = line.doppler_standard_deviation_hz(temperature)
    gamma = line.natural_lorentz_hwhm_hz + extra_gamma
    profile, frequency_weights, normalization = (
        _voigt_positive_frequency_quadrature(
            doppler_standard_deviation_hz=sigma,
            lorentz_hwhm_hz=gamma,
            quadrature_order=frequency_order,
            maximum_doppler_standard_deviations=(
                maximum_doppler_standard_deviations),
        )
    )
    classical_electron_radius_m = physical_constants[
        "classical electron radius"][0]
    integrated_cross_section_m2_hz = (
        np.pi * classical_electron_radius_m * SPEED_OF_LIGHT_M_S
        * line.absorption_oscillator_strength
    )
    absorption = density * integrated_cross_section_m2_hz * profile
    if density == 0.0:
        escape_probability = 1.0
        if math.isinf(scale):
            _, mean_path = _uniform_cylinder_transmission(
                geometry, np.asarray([0.0]),
                quadrature_order=geometry_order)
        else:
            mean_path = float(np.dot(path_weights, exit_path))
    else:
        if math.isinf(scale):
            transmissions, mean_path = _uniform_cylinder_transmission(
                geometry, absorption, quadrature_order=geometry_order)
        else:
            mean_path = float(np.dot(path_weights, exit_path))
            transmissions = np.empty_like(absorption)
            for start in range(0, absorption.size, 64):
                stop = min(absorption.size, start + 64)
                transmissions[start:stop] = np.exp(
                    -absorption[start:stop, None] * exit_path[None, :]
                ) @ path_weights
        # The analytic far-tail probability has transmission one.
        escape_probability = float(
            np.dot(frequency_weights, transmissions)
            + (normalization - np.sum(frequency_weights))
        )
        escape_probability /= normalization
    line_center_profile = float(voigt_profile(0.0, sigma, gamma))
    line_center_absorption = (
        density * integrated_cross_section_m2_hz * line_center_profile)
    return CylinderResonanceEscapeResult(
        escape_probability=escape_probability,
        trapping_factor=1.0 / escape_probability,
        absorber_density_m3=density,
        gas_temperature_K=temperature,
        source_axial_scale_length_m=scale,
        doppler_standard_deviation_hz=sigma,
        lorentz_hwhm_hz=gamma,
        line_center_absorption_coefficient_m_inv=line_center_absorption,
        source_weighted_mean_exit_path_m=mean_path,
        line_center_mean_path_optical_depth=(
            line_center_absorption * mean_path),
        frequency_profile_normalization=normalization,
        geometry_quadrature_order=geometry_order,
        frequency_quadrature_order=frequency_order,
    )


def deterministic_cylinder_partial_redistribution(
    geometry: CylindricalReactor,
    line: ResonanceLineData,
    *,
    absorber_density_m3: float,
    gas_temperature_K: float,
    velocity_changing_collision_frequency_s_inv: float,
    quenching_collision_frequency_s_inv: float = 0.0,
    additional_lorentz_hwhm_hz: float = 0.0,
    geometry_quadrature_order: int = 12,
    frequency_quadrature_order: int = 40,
    coherent_half_range_doppler_standard_deviations: float = 32.0,
    coherent_frequency_grid_points: int | None = None,
    coherent_grid_points_per_lorentz_hwhm: float = 12.0,
    linear_solver_relative_tolerance: float = 5.0e-9,
) -> CylinderPartialRedistributionResult:
    r"""Solve Tian's partial-frequency-redistribution cycle deterministically.

    A finite frequency state represents the emitted photon's detuning.  After
    resonant absorption, radiative decay competes with quenching.  If decay
    occurs before a velocity-changing collision, the next frequency walks by
    the natural Lorentz emission kernel centred on the absorbed frequency;
    otherwise it is redrawn from the local Voigt profile.  A zero-padded FFT
    applies that coherent propagator and preconditioned GMRES solves the
    resulting convolution-plus-rank-one renewal equation.  This is the
    deterministic counterpart of Tian Eqs. 2.30--2.39 and reduces to
    ``1/P_escape`` in the complete-redistribution, zero-quench limit.

    Collision frequencies are explicit physical inputs.  This operator does
    not infer them from a trapping target and does not by itself certify a
    reactor photon boundary.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("a cylindrical reactor geometry is required")
    if not isinstance(line, ResonanceLineData):
        raise TypeError("resonance-line atomic data are required")
    density = float(absorber_density_m3)
    temperature = float(gas_temperature_K)
    velocity_frequency = float(
        velocity_changing_collision_frequency_s_inv)
    quench_frequency = float(quenching_collision_frequency_s_inv)
    extra_gamma = float(additional_lorentz_hwhm_hz)
    geometry_order = int(geometry_quadrature_order)
    frequency_order = int(frequency_quadrature_order)
    half_range = float(coherent_half_range_doppler_standard_deviations)
    requested_grid_points = coherent_frequency_grid_points
    requested_points_per_hwhm = float(
        coherent_grid_points_per_lorentz_hwhm)
    solver_tolerance = float(linear_solver_relative_tolerance)
    if (
        not math.isfinite(density)
        or density < 0.0
        or not math.isfinite(temperature)
        or temperature <= 0.0
        or not math.isfinite(velocity_frequency)
        or velocity_frequency < 0.0
        or not math.isfinite(quench_frequency)
        or quench_frequency < 0.0
        or not math.isfinite(extra_gamma)
        or extra_gamma < 0.0
        or geometry_order != geometry_quadrature_order
        or frequency_order != frequency_quadrature_order
        or min(geometry_order, frequency_order) < 4
        or not math.isfinite(half_range)
        or half_range < 32.0
        or not math.isfinite(requested_points_per_hwhm)
        or requested_points_per_hwhm < 4.0
        or not math.isfinite(solver_tolerance)
        or not 0.0 < solver_tolerance <= 1.0e-7
    ):
        raise ValueError("invalid partial-redistribution condition")

    sigma = line.doppler_standard_deviation_hz(temperature)
    gamma = line.natural_lorentz_hwhm_hz + extra_gamma
    natural_gamma_ratio = line.natural_lorentz_hwhm_hz / sigma
    required_grid_points = max(
        4096,
        int(math.ceil(
            2.0 * half_range * requested_points_per_hwhm
            / natural_gamma_ratio
        )),
    )
    if requested_grid_points is None:
        grid_points = 1 << (required_grid_points - 1).bit_length()
    else:
        if (
            int(requested_grid_points) != requested_grid_points
            or int(requested_grid_points) < 4096
            or int(requested_grid_points) % 2
        ):
            raise ValueError(
                "coherent frequency grid must be an even integer >= 4096")
        grid_points = int(requested_grid_points)
    if grid_points > 1_048_576:
        raise ValueError(
            "coherent frequency grid exceeds the certified memory bound")
    grid_spacing = 2.0 * half_range / grid_points
    actual_points_per_hwhm = natural_gamma_ratio / grid_spacing
    if actual_points_per_hwhm < 4.0:
        raise ValueError(
            "coherent frequency grid resolves the natural Lorentz HWHM "
            "with fewer than four points")

    normalized_detuning = (
        -half_range
        + (np.arange(grid_points, dtype=float) + 0.5) * grid_spacing
    )
    dimensionless_profile = (
        voigt_profile(normalized_detuning * sigma, sigma, gamma) * sigma)
    # The omitted far wing is Lorentzian at this range.  It represents photons
    # that escape on the current emission, so their excess future-cycle count
    # is exactly zero in the renewal equation.
    far_tail_probability = 2.0 * gamma / (
        np.pi * sigma * half_range)
    normalization = float(
        np.sum(dimensionless_profile) * grid_spacing
        + far_tail_probability)
    normalized_profile = dimensionless_profile / normalization
    finite_profile_probability = float(
        np.sum(normalized_profile) * grid_spacing)
    normalized_tail_probability = max(
        0.0, 1.0 - finite_profile_probability)

    classical_electron_radius_m = physical_constants[
        "classical electron radius"][0]
    integrated_cross_section_m2_hz = (
        np.pi * classical_electron_radius_m * SPEED_OF_LIGHT_M_S
        * line.absorption_oscillator_strength)

    # Transmission is smooth on the Doppler scale even when the natural
    # coherent kernel is extremely narrow.  Evaluate it on an independent
    # refinement-controlled grid and interpolate onto the propagator grid.
    minimum_escape_intervals = max(4096, frequency_order * 256)
    escape_intervals = 1 << (minimum_escape_intervals - 1).bit_length()
    escape_detuning = np.linspace(
        -half_range, half_range, escape_intervals + 1)
    escape_profile = voigt_profile(
        escape_detuning * sigma, sigma, gamma)
    escape_absorption = (
        density * integrated_cross_section_m2_hz * escape_profile)
    escape_transmission, _ = _uniform_cylinder_transmission(
        geometry,
        escape_absorption,
        quadrature_order=geometry_order,
    )
    transmission = np.interp(
        normalized_detuning, escape_detuning, escape_transmission)
    complete_escape = min(1.0, max(0.0, float(
        np.dot(normalized_profile, transmission) * grid_spacing
        + normalized_tail_probability
    )))
    complete_trapping = 1.0 / complete_escape

    redistribution_probability = (
        velocity_frequency
        / (line.transition_probability_s_inv + velocity_frequency))
    coherent_probability = 1.0 - redistribution_probability
    survival_probability = (
        line.transition_probability_s_inv
        / (line.transition_probability_s_inv + quench_frequency))
    absorbed_and_surviving = np.clip(
        (1.0 - transmission) * survival_probability, 0.0, 1.0)

    # K(x-y) is the natural Lorentz redistribution kernel in Doppler units.
    # Zero padding prevents the unphysical periodic wraparound of an FFT on a
    # finite interval.  The unknown is u=t-1; outside the interval u=0 because
    # a far-wing photon exits on its current emission.
    difference = (
        np.arange(grid_points, dtype=float) - grid_points // 2
    ) * grid_spacing
    coherent_kernel = (
        natural_gamma_ratio / np.pi
        / (difference * difference + natural_gamma_ratio ** 2)
    )
    fft_points = next_fast_len(2 * grid_points - 1)
    kernel_transform = rfft(coherent_kernel, fft_points)
    convolution_start = grid_points // 2

    def coherent_convolution(vector: np.ndarray) -> np.ndarray:
        full = irfft(
            rfft(np.asarray(vector, dtype=float), fft_points)
            * kernel_transform,
            fft_points,
        )
        return (
            full[convolution_start:convolution_start + grid_points]
            * grid_spacing
        )

    def renewal_matvec(vector: np.ndarray) -> np.ndarray:
        values = np.asarray(vector, dtype=float)
        redistributed_mean = float(
            np.dot(normalized_profile, values) * grid_spacing)
        return values - absorbed_and_surviving * (
            coherent_probability * coherent_convolution(values)
            + redistribution_probability * redistributed_mean
        )

    rhs = absorbed_and_surviving
    rhs_norm = float(np.linalg.norm(rhs))
    iteration_residuals: list[float] = []
    if rhs_norm == 0.0:
        excess_cycles = np.zeros_like(rhs)
        solver_iterations = 0
        solver_relative_residual = 0.0
    else:
        operator = LinearOperator(
            (grid_points, grid_points), matvec=renewal_matvec, dtype=float)
        approximate_diagonal = (
            1.0
            - absorbed_and_surviving * coherent_probability
            * coherent_kernel[grid_points // 2] * grid_spacing
        )
        if np.any(approximate_diagonal <= 0.0):
            raise RuntimeError(
                "partial-redistribution preconditioner is not positive")
        preconditioner = LinearOperator(
            (grid_points, grid_points),
            matvec=lambda vector: (
                np.asarray(vector, dtype=float) / approximate_diagonal),
            dtype=float,
        )
        excess_cycles, solver_info = gmres(
            operator,
            rhs,
            M=preconditioner,
            rtol=solver_tolerance,
            atol=min(1.0e-12, solver_tolerance * rhs_norm),
            restart=100,
            maxiter=500,
            callback=iteration_residuals.append,
            callback_type="pr_norm",
        )
        if solver_info != 0:
            raise RuntimeError(
                "partial-redistribution GMRES did not converge "
                f"(info={solver_info})")
        solver_iterations = len(iteration_residuals)
        solver_relative_residual = float(
            np.linalg.norm(renewal_matvec(excess_cycles) - rhs) / rhs_norm)
        if solver_relative_residual > max(1.0e-7, 10.0 * solver_tolerance):
            raise RuntimeError(
                "partial-redistribution residual exceeds tolerance")

    partial_trapping = max(1.0, float(
        1.0
        + np.dot(normalized_profile, excess_cycles) * grid_spacing
    ))
    return CylinderPartialRedistributionResult(
        complete_redistribution_escape_probability=complete_escape,
        complete_redistribution_trapping_factor=complete_trapping,
        partial_redistribution_trapping_factor=partial_trapping,
        coherent_redistribution_probability=coherent_probability,
        velocity_redistribution_probability=redistribution_probability,
        radiative_survival_probability_after_absorption=survival_probability,
        velocity_changing_collision_frequency_s_inv=velocity_frequency,
        quenching_collision_frequency_s_inv=quench_frequency,
        frequency_profile_normalization=normalization,
        coherent_half_range_doppler_standard_deviations=half_range,
        coherent_frequency_grid_points=grid_points,
        coherent_grid_points_per_lorentz_hwhm=actual_points_per_hwhm,
        linear_solver_iterations=solver_iterations,
        linear_solver_relative_residual=solver_relative_residual,
        geometry_quadrature_order=geometry_order,
        frequency_quadrature_order=frequency_order,
    )


def _disk_overlap_area_m2(
    separation_m: np.ndarray,
    first_radius_m: float,
    second_radius_m: float,
) -> np.ndarray:
    """Area shared by two offset disks, evaluated away from branch points."""
    distance = np.asarray(separation_m, dtype=float)
    large = max(first_radius_m, second_radius_m)
    small = min(first_radius_m, second_radius_m)
    overlap = np.empty_like(distance)
    contained = distance <= large - small
    overlap[contained] = np.pi * small ** 2
    partial = ~contained
    d = distance[partial]
    first_cosine = np.clip(
        (d * d + first_radius_m ** 2 - second_radius_m ** 2)
        / (2.0 * d * first_radius_m), -1.0, 1.0)
    second_cosine = np.clip(
        (d * d + second_radius_m ** 2 - first_radius_m ** 2)
        / (2.0 * d * second_radius_m), -1.0, 1.0)
    radical = np.maximum(
        0.0,
        (-d + first_radius_m + second_radius_m)
        * (d + first_radius_m - second_radius_m)
        * (d - first_radius_m + second_radius_m)
        * (d + first_radius_m + second_radius_m),
    )
    overlap[partial] = (
        first_radius_m ** 2 * np.arccos(first_cosine)
        + second_radius_m ** 2 * np.arccos(second_cosine)
        - 0.5 * np.sqrt(radical)
    )
    return overlap


def _axial_cosine_kernel_m_inv(
    horizontal_separation_m: np.ndarray,
    length_m: float,
    extinction_coefficient_m_inv: float,
) -> np.ndarray:
    r"""Return ``integral_0^L z exp(-k s) / s^3 dz`` analytically."""
    distance = np.asarray(horizontal_separation_m, dtype=float)
    upper = np.sqrt(distance * distance + length_m * length_m)
    if extinction_coefficient_m_inv == 0.0:
        return 1.0 / distance - 1.0 / upper
    kappa = extinction_coefficient_m_inv
    return (
        np.exp(-kappa * distance) / distance
        - np.exp(-kappa * upper) / upper
        - kappa * (exp1(kappa * distance) - exp1(kappa * upper))
    )


def _lamellar_direct_survival(shift_over_width: np.ndarray) -> np.ndarray:
    r"""Average ``max(0, 1-u*abs(cos(phi)))`` over source azimuth."""
    ratio = np.asarray(shift_over_width, dtype=float)
    if np.any(ratio < 0.0) or np.any(~np.isfinite(ratio)):
        raise ValueError("shift-to-opening ratio must be finite and nonnegative")
    result = np.empty_like(ratio)
    small = ratio <= 1.0
    result[small] = 1.0 - 2.0 / np.pi * ratio[small]
    large_ratio = ratio[~small]
    result[~small] = 2.0 / np.pi * (
        np.arcsin(1.0 / large_ratio)
        - large_ratio
        + np.sqrt(large_ratio * large_ratio - 1.0)
    )
    return np.clip(result, 0.0, 1.0)


def uniform_isotropic_cylinder_to_disk_transfer(
    geometry: CylindricalReactor,
    *,
    wafer_radius_m: float,
    volume_emissivity_m3_s: float,
    extinction_coefficient_m_inv: float = 0.0,
    quadrature_order: int = 24,
) -> CylinderDiskRadiationResult:
    r"""Integrate direct isotropic radiation from a cylinder to its end disk.

    For source point ``x`` and wafer point ``y``, the differential receipt is

    ``q exp(-k s) cos(theta) dV dA / (4 pi s**2)``.

    Axisymmetry removes one azimuth analytically.  The remaining source radius,
    source height, target radius, and relative azimuth are evaluated by a
    tensor-product Gauss rule.  The cylinder is convex, so every integrated ray
    remains inside the declared gas volume.  Reflections and re-emission are
    intentionally absent.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("a cylindrical reactor geometry is required")
    wafer_radius = float(wafer_radius_m)
    emissivity = float(volume_emissivity_m3_s)
    extinction = float(extinction_coefficient_m_inv)
    order = int(quadrature_order)
    if (
        not math.isfinite(wafer_radius)
        or not 0.0 < wafer_radius <= geometry.radius_m
        or not math.isfinite(emissivity)
        or emissivity < 0.0
        or not math.isfinite(extinction)
        or extinction < 0.0
        or order != quadrature_order
        or order < 4
    ):
        raise ValueError("invalid cylinder-to-disk radiation condition")

    # The horizontal double-area integral is a convolution of two disk
    # indicators.  At separation rho, its angular measure is
    # 2*pi*rho*A_overlap(rho).  Integrating z analytically leaves one smooth
    # radial integral.  Split at the containment/partial-overlap kink.
    split = geometry.radius_m - wafer_radius
    intervals = []
    if split > 0.0:
        intervals.append((0.0, split))
    intervals.append((split, geometry.radius_m + wafer_radius))
    radial_integral_m3 = 0.0
    for lower, upper in intervals:
        if upper <= lower:
            continue
        separation, separation_w = _mapped_legendre_interval(
            order, lower, upper)
        overlap = _disk_overlap_area_m2(
            separation, geometry.radius_m, wafer_radius)
        axial = _axial_cosine_kernel_m_inv(
            separation, geometry.length_m, extinction)
        radial_integral_m3 += float(np.sum(
            separation_w * separation * overlap * axial
        ))
    # The 1/2 is the horizontal convolution's 2*pi divided by the 4*pi
    # isotropic-emission denominator.
    wafer_area = np.pi * wafer_radius ** 2
    intercept_probability = radial_integral_m3 / (2.0 * geometry.volume_m3)
    geometry_flux_length = (
        geometry.volume_m3 * intercept_probability / wafer_area
    )
    photon_flux = emissivity * geometry_flux_length
    return CylinderDiskRadiationResult(
        volume_emissivity_m3_s=emissivity,
        wafer_radius_m=wafer_radius,
        extinction_coefficient_m_inv=extinction,
        wafer_photon_flux_m2_s=photon_flux,
        wafer_intercept_probability=intercept_probability,
        geometry_flux_length_m=geometry_flux_length,
        quadrature_order=order,
    )


def uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
    geometry: CylindricalReactor,
    *,
    wafer_radius_m: float,
    opening_width_m: float,
    optical_path_depth_m: float,
    volume_emissivity_m3_s: float,
    extinction_coefficient_m_inv: float = 0.0,
    radial_quadrature_order: int = 32,
    axial_quadrature_order: int = 32,
) -> LamellarFloorRadiationResult:
    r"""Integrate unreflected rays reaching an absorbing line-trench floor.

    A ray emitted at horizontal source/target separation ``rho`` and height
    ``z`` shifts by ``depth*rho/z`` between the opening and floor.  Averaging
    its projection onto the line-normal direction and its uniformly sampled
    entry point gives :func:`_lamellar_direct_survival`.  The remaining
    cylinder/disk convolution is evaluated by fixed Gauss--Legendre rules.

    This is exact geometrical transport for uniform isotropic volume emission,
    a perfectly absorbing mask/sidewall, and an infinitely long line opening.
    It is not an electromagnetic approximation at openings only a few
    wavelengths wide; callers must retain that separate diffraction gate.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("a cylindrical reactor geometry is required")
    wafer_radius = float(wafer_radius_m)
    width = float(opening_width_m)
    depth = float(optical_path_depth_m)
    emissivity = float(volume_emissivity_m3_s)
    extinction = float(extinction_coefficient_m_inv)
    radial_order = int(radial_quadrature_order)
    axial_order = int(axial_quadrature_order)
    if (
        not math.isfinite(wafer_radius)
        or not 0.0 < wafer_radius <= geometry.radius_m
        or not math.isfinite(width)
        or width <= 0.0
        or not math.isfinite(depth)
        or depth < 0.0
        or not math.isfinite(emissivity)
        or emissivity < 0.0
        or not math.isfinite(extinction)
        or extinction < 0.0
        or radial_order != radial_quadrature_order
        or axial_order != axial_quadrature_order
        or min(radial_order, axial_order) < 4
    ):
        raise ValueError("invalid cylinder-to-lamellar-floor condition")
    wafer_transfer = uniform_isotropic_cylinder_to_disk_transfer(
        geometry,
        wafer_radius_m=wafer_radius,
        volume_emissivity_m3_s=emissivity,
        extinction_coefficient_m_inv=extinction,
        quadrature_order=radial_order,
    )
    if depth == 0.0:
        fraction = 1.0
        return LamellarFloorRadiationResult(
            wafer_plane_photon_flux_m2_s=wafer_transfer.wafer_photon_flux_m2_s,
            trench_floor_photon_flux_m2_s=(
                wafer_transfer.wafer_photon_flux_m2_s * fraction),
            direct_floor_fraction=fraction,
            opening_width_m=width,
            optical_path_depth_m=depth,
            radial_quadrature_order=radial_order,
            axial_quadrature_order=axial_order,
        )

    axial_z, axial_w = _mapped_legendre(axial_order, geometry.length_m)
    split = geometry.radius_m - wafer_radius
    intervals = []
    if split > 0.0:
        intervals.append((0.0, split))
    intervals.append((split, geometry.radius_m + wafer_radius))
    shadowed_integral_m3 = 0.0
    clear_integral_m3 = 0.0
    for lower, upper in intervals:
        if upper <= lower:
            continue
        separation, separation_w = _mapped_legendre_interval(
            radial_order, lower, upper)
        overlap = _disk_overlap_area_m2(
            separation, geometry.radius_m, wafer_radius)
        rho = separation[:, None]
        z = axial_z[None, :]
        distance = np.sqrt(rho * rho + z * z)
        kernel = z / (distance ** 3)
        if extinction > 0.0:
            kernel *= np.exp(-extinction * distance)
        survival = _lamellar_direct_survival(depth * rho / (width * z))
        clear_axial = np.sum(axial_w[None, :] * kernel, axis=1)
        shadowed_axial = np.sum(
            axial_w[None, :] * kernel * survival, axis=1)
        radial_weight = separation_w * separation * overlap
        clear_integral_m3 += float(np.sum(radial_weight * clear_axial))
        shadowed_integral_m3 += float(np.sum(radial_weight * shadowed_axial))
    if clear_integral_m3 <= 0.0:
        raise RuntimeError("nonpositive direct cylinder-to-wafer integral")
    fraction = shadowed_integral_m3 / clear_integral_m3
    return LamellarFloorRadiationResult(
        wafer_plane_photon_flux_m2_s=wafer_transfer.wafer_photon_flux_m2_s,
        trench_floor_photon_flux_m2_s=(
            wafer_transfer.wafer_photon_flux_m2_s * fraction),
        direct_floor_fraction=fraction,
        opening_width_m=width,
        optical_path_depth_m=depth,
        radial_quadrature_order=radial_order,
        axial_quadrature_order=axial_order,
    )


@dataclass(frozen=True)
class KemaneciCl139nmEmissionSensitivity:
    """Kemaneci reaction-18 atomic-Cl excitation interpreted as 139-nm light.

    The Maxwellian rate is a published compilation and the wavelength/state
    assignment follows the Kemaneci source.  The survival fraction is not
    supplied by that global model: resonance trapping, collisional quenching,
    and wall loss make it an explicit sensitivity.  Therefore this card never
    grants a predictive photon boundary by itself.
    """

    radiative_survival_fraction: float = 1.0
    photon_wavelength_nm: float = 139.0

    def __post_init__(self):
        survival = float(self.radiative_survival_fraction)
        wavelength = float(self.photon_wavelength_nm)
        if (
            not math.isfinite(survival)
            or not 0.0 <= survival <= 1.0
            or not math.isfinite(wavelength)
            or wavelength <= 0.0
        ):
            raise ValueError("invalid atomic-chlorine emission sensitivity")
        object.__setattr__(self, "radiative_survival_fraction", survival)
        object.__setattr__(self, "photon_wavelength_nm", wavelength)

    @staticmethod
    def excitation_rate_coefficient_m3_s(electron_temperature_eV: float) -> float:
        from .chlorine_kemaneci import build_kemaneci_2014_forward_chlorine_network
        from .network import RateContext

        temperature = float(electron_temperature_eV)
        if not math.isfinite(temperature) or not 0.5 <= temperature <= 10.0:
            raise ValueError("Kemaneci excitation fit requires 0.5 <= Te <= 10 eV")
        network = build_kemaneci_2014_forward_chlorine_network()
        reaction = next(
            item for item in network.reactions
            if item.name == "k18_Cl_excitation_1P5_2"
        )
        return float(reaction.rate_coefficient.coefficient_si(
            RateContext(electron_temperature_eV=temperature)
        ))

    def primary_excitation_rate_m3_s(
        self,
        *,
        electron_density_m3: float,
        chlorine_atom_density_m3: float,
        electron_temperature_eV: float,
    ) -> float:
        electron_density = float(electron_density_m3)
        chlorine_density = float(chlorine_atom_density_m3)
        if (
            not math.isfinite(electron_density)
            or electron_density < 0.0
            or not math.isfinite(chlorine_density)
            or chlorine_density < 0.0
        ):
            raise ValueError("emitter densities must be finite and nonnegative")
        return float(
            self.excitation_rate_coefficient_m3_s(electron_temperature_eV)
            * electron_density * chlorine_density
        )

    def escaping_emissivity_m3_s(self, **condition: float) -> float:
        return float(
            self.radiative_survival_fraction
            * self.primary_excitation_rate_m3_s(**condition)
        )

    @property
    def supports_prediction(self) -> bool:
        return False
