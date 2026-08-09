"""Deterministic resonance transport over axisymmetric reactor moments.

The global chemistry can remain zero-dimensional while exposing a small,
fixed axisymmetric moment field for radiation: gas temperature, ground-state
absorber density, and a line-specific emitter density.  This module integrates
those fields with fixed source, direction, path, and frequency quadrature.  It
is deterministic and embarrassingly parallel over conditions and lines.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import (
    c as SPEED_OF_LIGHT_M_S,
    k as BOLTZMANN_CONSTANT_J_K,
    physical_constants,
)
from scipy.special import voigt_profile

from .geometry import CylindricalReactor
from .vuv_radiation import ResonanceLineData
from .vuv_radiation import deterministic_cylinder_resonance_escape


@dataclass(frozen=True)
class AxisymmetricRadiationMomentField:
    """Nodal axisymmetric fields consumed by deterministic ray transport."""

    radial_nodes_m: np.ndarray
    axial_nodes_m: np.ndarray
    gas_temperature_K: np.ndarray
    absorber_density_m3: np.ndarray
    emitter_density_m3: np.ndarray
    source: str

    def __post_init__(self):
        radial = np.array(self.radial_nodes_m, dtype=float, copy=True)
        axial = np.array(self.axial_nodes_m, dtype=float, copy=True)
        temperature = np.array(self.gas_temperature_K, dtype=float, copy=True)
        absorber = np.array(self.absorber_density_m3, dtype=float, copy=True)
        emitter = np.array(self.emitter_density_m3, dtype=float, copy=True)
        expected = (radial.size, axial.size)
        if (
            radial.ndim != 1
            or axial.ndim != 1
            or min(radial.size, axial.size) < 2
            or np.any(~np.isfinite(radial))
            or np.any(~np.isfinite(axial))
            or radial[0] != 0.0
            or axial[0] != 0.0
            or np.any(np.diff(radial) <= 0.0)
            or np.any(np.diff(axial) <= 0.0)
            or temperature.shape != expected
            or absorber.shape != expected
            or emitter.shape != expected
            or np.any(~np.isfinite(temperature))
            or np.any(~np.isfinite(absorber))
            or np.any(~np.isfinite(emitter))
            or np.any(temperature <= 0.0)
            or np.any(absorber < 0.0)
            or np.any(emitter < 0.0)
            or not np.any(emitter > 0.0)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid axisymmetric radiation moment field")
        for array in (radial, axial, temperature, absorber, emitter):
            array.setflags(write=False)
        object.__setattr__(self, "radial_nodes_m", radial)
        object.__setattr__(self, "axial_nodes_m", axial)
        object.__setattr__(self, "gas_temperature_K", temperature)
        object.__setattr__(self, "absorber_density_m3", absorber)
        object.__setattr__(self, "emitter_density_m3", emitter)

    @property
    def geometry(self) -> CylindricalReactor:
        return CylindricalReactor(
            radius_m=float(self.radial_nodes_m[-1]),
            length_m=float(self.axial_nodes_m[-1]),
        )


@dataclass(frozen=True)
class AxisymmetricResonanceEscapeResult:
    escape_probability: float
    trapping_factor: float
    source_weighted_temperature_K: float
    source_weighted_absorber_density_m3: float
    source_weighted_mean_exit_path_m: float
    frequency_profile_normalization_minimum: float
    frequency_profile_normalization_maximum: float
    source_quadrature_order: int
    direction_quadrature_order: int
    path_quadrature_order: int
    frequency_quadrature_order: int

    def __post_init__(self):
        values = np.asarray((
            self.escape_probability,
            self.trapping_factor,
            self.source_weighted_temperature_K,
            self.source_weighted_absorber_density_m3,
            self.source_weighted_mean_exit_path_m,
            self.frequency_profile_normalization_minimum,
            self.frequency_profile_normalization_maximum,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or not 0.0 < self.escape_probability <= 1.0
            or not np.isclose(
                self.trapping_factor, 1.0 / self.escape_probability,
                rtol=2.0e-14, atol=0.0)
            or self.source_weighted_temperature_K <= 0.0
            or self.source_weighted_absorber_density_m3 < 0.0
            or self.source_weighted_mean_exit_path_m <= 0.0
            or not 0.99999 <= self.frequency_profile_normalization_minimum
            or not self.frequency_profile_normalization_maximum <= 1.00001
            or any(
                int(order) != order or order < 4
                for order in (
                    self.source_quadrature_order,
                    self.direction_quadrature_order,
                    self.path_quadrature_order,
                    self.frequency_quadrature_order,
                )
            )
        ):
            raise ValueError("invalid axisymmetric resonance-escape result")
        for name, value in zip((
            "escape_probability", "trapping_factor",
            "source_weighted_temperature_K",
            "source_weighted_absorber_density_m3",
            "source_weighted_mean_exit_path_m",
            "frequency_profile_normalization_minimum",
            "frequency_profile_normalization_maximum",
        ), values):
            object.__setattr__(self, name, float(value))
        for name in (
            "source_quadrature_order", "direction_quadrature_order",
            "path_quadrature_order", "frequency_quadrature_order",
        ):
            object.__setattr__(self, name, int(getattr(self, name)))


@dataclass(frozen=True)
class AxisymmetricResonanceConvergenceReceipt:
    """Coarse/refined numerical receipt required before prediction use."""

    coarse: AxisymmetricResonanceEscapeResult
    refined: AxisymmetricResonanceEscapeResult
    relative_escape_change: float
    relative_trapping_change: float
    relative_tolerance: float
    converged: bool

    def __post_init__(self):
        maximum_change = max(
            float(self.relative_escape_change),
            float(self.relative_trapping_change),
        )
        if (
            not isinstance(self.coarse, AxisymmetricResonanceEscapeResult)
            or not isinstance(self.refined, AxisymmetricResonanceEscapeResult)
            or not math.isfinite(maximum_change)
            or min(
                self.relative_escape_change,
                self.relative_trapping_change,
            ) < 0.0
            or not math.isfinite(self.relative_tolerance)
            or self.relative_tolerance <= 0.0
            or bool(self.converged)
            != (maximum_change <= self.relative_tolerance)
        ):
            raise ValueError("invalid axisymmetric convergence receipt")
        object.__setattr__(
            self, "relative_escape_change", float(self.relative_escape_change))
        object.__setattr__(
            self, "relative_trapping_change",
            float(self.relative_trapping_change))
        object.__setattr__(
            self, "relative_tolerance", float(self.relative_tolerance))
        object.__setattr__(self, "converged", bool(self.converged))


def _bilinear(
    radial_nodes: np.ndarray,
    axial_nodes: np.ndarray,
    values: np.ndarray,
    radius: np.ndarray,
    axial: np.ndarray,
) -> np.ndarray:
    r = np.clip(np.asarray(radius, dtype=float), radial_nodes[0], radial_nodes[-1])
    z = np.clip(np.asarray(axial, dtype=float), axial_nodes[0], axial_nodes[-1])
    ri = np.clip(np.searchsorted(radial_nodes, r, side="right") - 1,
                 0, radial_nodes.size - 2)
    zi = np.clip(np.searchsorted(axial_nodes, z, side="right") - 1,
                 0, axial_nodes.size - 2)
    r0 = radial_nodes[ri]
    r1 = radial_nodes[ri + 1]
    z0 = axial_nodes[zi]
    z1 = axial_nodes[zi + 1]
    fr = (r - r0) / (r1 - r0)
    fz = (z - z0) / (z1 - z0)
    return (
        (1.0 - fr) * (1.0 - fz) * values[ri, zi]
        + fr * (1.0 - fz) * values[ri + 1, zi]
        + (1.0 - fr) * fz * values[ri, zi + 1]
        + fr * fz * values[ri + 1, zi + 1]
    )


def _positive_voigt_grid(
    sigma_hz: float,
    gamma_hz: float,
    order: int,
    maximum_sigma: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    nodes, weights = leggauss(order)
    bounds = [0.0, 1.0]
    while bounds[-1] < maximum_sigma:
        bounds.append(min(maximum_sigma, 2.0 * bounds[-1]))
    scaled = []
    scaled_weights = []
    for lower, upper in zip(bounds, bounds[1:]):
        half = 0.5 * (upper - lower)
        scaled.append(half * nodes + 0.5 * (upper + lower))
        scaled_weights.append(half * weights)
    y = np.concatenate(scaled)
    detuning = sigma_hz * y
    profile = voigt_profile(detuning, sigma_hz, gamma_hz)
    probability_weights = (
        2.0 * sigma_hz * np.concatenate(scaled_weights) * profile)
    tail = 2.0 * gamma_hz / (np.pi * sigma_hz * maximum_sigma)
    return detuning, probability_weights, float(
        np.sum(probability_weights) + tail)


def _axisymmetric_projected_chords(
    geometry: CylindricalReactor,
    *,
    surface_order: int,
    direction_order: int,
    endcap_radial_break_m: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Return inward chords sampled in projected surface-ray measure.

    Every oriented volume ray in a convex body can be represented uniquely by
    its exit point and inward direction.  The Jacobian is ``|n dot omega|``.
    This parameterization makes optically thick transport a surface problem
    instead of asking volume quadrature to resolve a vanishingly thin escape
    layer.  Common angular constants cancel in the normalized escape ratio;
    the retained end/side weights preserve the physical surface-area ratio.
    """
    surface_nodes, surface_weights = leggauss(surface_order)
    direction_nodes, direction_weights = leggauss(direction_order)
    surface_unit = 0.5 * (surface_nodes + 1.0)
    surface_w = 0.5 * surface_weights
    direction_unit = 0.5 * (direction_nodes + 1.0)
    direction_w = 0.5 * direction_weights

    radial_break = (
        None if endcap_radial_break_m is None
        else float(endcap_radial_break_m)
    )
    if radial_break is not None and (
        not math.isfinite(radial_break)
        or not 0.0 < radial_break < geometry.radius_m
    ):
        raise ValueError("endcap radial quadrature break lies outside disk")
    if radial_break is None:
        endcap_area_unit = surface_unit
        endcap_surface_w = surface_w
    else:
        # A partial wafer introduces a step in the terminal-boundary label at
        # r=R_wafer. Split the disk-area coordinate there exactly; otherwise a
        # single Gauss rule oscillates as nodes cross the step with order.
        area_break = (radial_break / geometry.radius_m) ** 2
        endcap_area_unit = np.concatenate((
            0.5 * area_break * (surface_nodes + 1.0),
            area_break + 0.5 * (1.0 - area_break) * (surface_nodes + 1.0),
        ))
        endcap_surface_w = np.concatenate((
            0.5 * area_break * surface_weights,
            0.5 * (1.0 - area_break) * surface_weights,
        ))

    exit_r_parts = []
    exit_z_parts = []
    direction_parts = []
    chord_parts = []
    weight_parts = []

    # The two circular end caps must remain separate because the moment field
    # need not be symmetric in z.  Disk radius is sampled uniformly in area.
    radial = geometry.radius_m * np.sqrt(endcap_area_unit)
    mu = direction_unit
    azimuth = 2.0 * np.pi * direction_unit
    radius0, axial_cosine, phi = np.meshgrid(
        radial, mu, azimuth, indexing="ij")
    transverse = np.sqrt(np.maximum(0.0, 1.0 - axial_cosine ** 2))
    dx = transverse * np.cos(phi)
    dy = transverse * np.sin(phi)
    radial_path = (
        -radius0 * dx
        + np.sqrt(np.maximum(
            0.0,
            transverse ** 2 * geometry.radius_m ** 2
            - radius0 ** 2 * dy ** 2,
        ))
    ) / (transverse ** 2)
    end_base_weight = (
        np.pi * geometry.radius_m ** 2
        * np.einsum(
            "i,j,k->ijk", endcap_surface_w, direction_w, direction_w)
        * axial_cosine
    )
    for exit_z, dz in ((0.0, 1.0), (geometry.length_m, -1.0)):
        opposite_cap_path = geometry.length_m / axial_cosine
        chord = np.minimum(radial_path, opposite_cap_path)
        exit_r_parts.append(radius0.ravel())
        exit_z_parts.append(np.full(radius0.size, exit_z))
        direction_parts.append(np.column_stack((
            dx.ravel(), dy.ravel(), np.full(radius0.size, dz)
            * axial_cosine.ravel(),
        )))
        chord_parts.append(chord.ravel())
        weight_parts.append(end_base_weight.ravel())

    # Cylindrical side.  The exit point is fixed at azimuth zero because all
    # supplied fields are axisymmetric; its integrated circumference remains
    # in the surface-area weight.
    axial = geometry.length_m * surface_unit
    eta = direction_unit
    tangent_angle = 2.0 * np.pi * direction_unit
    axial0, radial_cosine, psi = np.meshgrid(
        axial, eta, tangent_angle, indexing="ij")
    tangent = np.sqrt(np.maximum(0.0, 1.0 - radial_cosine ** 2))
    dy = tangent * np.cos(psi)
    dz = tangent * np.sin(psi)
    horizontal_norm_squared = radial_cosine ** 2 + dy ** 2
    opposite_side_path = (
        2.0 * geometry.radius_m * radial_cosine
        / horizontal_norm_squared)
    end_path = np.where(
        dz > 0.0,
        (geometry.length_m - axial0) / dz,
        axial0 / (-dz),
    )
    chord = np.minimum(opposite_side_path, end_path)
    side_weight = (
        2.0 * np.pi * geometry.radius_m * geometry.length_m
        * np.einsum("i,j,k->ijk", surface_w, direction_w, direction_w)
        * radial_cosine
    )
    exit_r_parts.append(np.full(axial0.size, geometry.radius_m))
    exit_z_parts.append(axial0.ravel())
    direction_parts.append(np.column_stack((
        -radial_cosine.ravel(), dy.ravel(), dz.ravel(),
    )))
    chord_parts.append(chord.ravel())
    weight_parts.append(side_weight.ravel())

    return (
        np.concatenate(exit_r_parts),
        np.concatenate(exit_z_parts),
        np.concatenate(direction_parts, axis=0),
        np.concatenate(chord_parts),
        np.concatenate(weight_parts),
    )


def deterministic_axisymmetric_resonance_escape(
    field: AxisymmetricRadiationMomentField,
    line: ResonanceLineData,
    *,
    additional_lorentz_hwhm_hz: float = 0.0,
    source_quadrature_order: int = 6,
    direction_quadrature_order: int = 6,
    path_quadrature_order: int = 12,
    frequency_quadrature_order: int = 24,
    maximum_doppler_standard_deviations: float = 1024.0,
) -> AxisymmetricResonanceEscapeResult:
    """Integrate Voigt absorption over axisymmetric reactor moment fields."""
    if not isinstance(field, AxisymmetricRadiationMomentField):
        raise TypeError("an axisymmetric radiation moment field is required")
    if not isinstance(line, ResonanceLineData):
        raise TypeError("resonance-line atomic data are required")
    orders = tuple(int(item) for item in (
        source_quadrature_order, direction_quadrature_order,
        path_quadrature_order, frequency_quadrature_order,
    ))
    if (
        any(value != integer or integer < 4 for value, integer in zip((
            source_quadrature_order, direction_quadrature_order,
            path_quadrature_order, frequency_quadrature_order,
        ), orders))
        or not math.isfinite(additional_lorentz_hwhm_hz)
        or additional_lorentz_hwhm_hz < 0.0
        or not math.isfinite(maximum_doppler_standard_deviations)
        or maximum_doppler_standard_deviations < 32.0
    ):
        raise ValueError("invalid axisymmetric radiation quadrature")
    source_order, direction_order, path_order, frequency_order = orders
    geometry = field.geometry

    # Exact uniform-field limit uses the projected-chord formulation, which
    # resolves the optically thick boundary layer analytically along chords.
    if (
        np.all(field.gas_temperature_K == field.gas_temperature_K.flat[0])
        and np.all(
            field.absorber_density_m3 == field.absorber_density_m3.flat[0])
        and np.all(field.emitter_density_m3 == field.emitter_density_m3.flat[0])
    ):
        homogeneous = deterministic_cylinder_resonance_escape(
            geometry,
            line,
            absorber_density_m3=float(field.absorber_density_m3.flat[0]),
            gas_temperature_K=float(field.gas_temperature_K.flat[0]),
            additional_lorentz_hwhm_hz=additional_lorentz_hwhm_hz,
            geometry_quadrature_order=max(8, direction_order),
            frequency_quadrature_order=frequency_order,
            maximum_doppler_standard_deviations=(
                maximum_doppler_standard_deviations),
        )
        return AxisymmetricResonanceEscapeResult(
            escape_probability=homogeneous.escape_probability,
            trapping_factor=homogeneous.trapping_factor,
            source_weighted_temperature_K=float(
                field.gas_temperature_K.flat[0]),
            source_weighted_absorber_density_m3=float(
                field.absorber_density_m3.flat[0]),
            source_weighted_mean_exit_path_m=(
                homogeneous.source_weighted_mean_exit_path_m),
            frequency_profile_normalization_minimum=(
                homogeneous.frequency_profile_normalization),
            frequency_profile_normalization_maximum=(
                homogeneous.frequency_profile_normalization),
            source_quadrature_order=source_order,
            direction_quadrature_order=direction_order,
            path_quadrature_order=path_order,
            frequency_quadrature_order=frequency_order,
        )

    exit_r, exit_z, inward_direction, chord, chord_weight = (
        _axisymmetric_projected_chords(
            geometry,
            surface_order=source_order,
            direction_order=direction_order,
        )
    )
    gamma = line.natural_lorentz_hwhm_hz + additional_lorentz_hwhm_hz
    classical_radius = physical_constants["classical electron radius"][0]
    cross_section_integral = (
        np.pi * classical_radius * SPEED_OF_LIGHT_M_S
        * line.absorption_oscillator_strength)

    # A common dimensional frequency grid is required because emission and
    # absorption profiles may have different local Doppler widths.  Scaling
    # it by the coldest supplied temperature resolves every local line core.
    sigma_reference = line.doppler_standard_deviation_hz(
        float(np.min(field.gas_temperature_K)))
    frequency_nodes, frequency_weights = leggauss(frequency_order)
    bounds = [0.0, 1.0]
    while bounds[-1] < maximum_doppler_standard_deviations:
        bounds.append(min(
            maximum_doppler_standard_deviations, 2.0 * bounds[-1]))
    detuning_parts = []
    frequency_weight_parts = []
    for lower, upper in zip(bounds, bounds[1:]):
        half = 0.5 * (upper - lower)
        detuning_parts.append(
            sigma_reference
            * (half * frequency_nodes + 0.5 * (upper + lower)))
        # Symmetric positive-frequency integration.
        frequency_weight_parts.append(
            2.0 * sigma_reference * half * frequency_weights)
    detuning = np.concatenate(detuning_parts)
    frequency_weight = np.concatenate(frequency_weight_parts)
    maximum_detuning = (
        sigma_reference * maximum_doppler_standard_deviations)
    # Leading asymptotic integral of both Lorentz wings.  At this cutoff its
    # opacity is negligible for the declared Tian/Lam audit envelope.
    profile_tail = 2.0 * gamma / (np.pi * maximum_detuning)

    source_measure = 0.0
    escaped_measure = 0.0
    temperature_measure = 0.0
    absorber_measure = 0.0
    exit_path_measure = 0.0
    normalization_minimum = math.inf
    normalization_maximum = -math.inf
    # Equal finite-volume segments plus exact exponential attenuation within
    # each segment resolve arbitrarily large optical depth without a spatial
    # boundary layer.  The segment count controls only field-variation error.
    segment_unit = (np.arange(path_order, dtype=float) + 0.5) / path_order
    edge_unit = np.arange(path_order + 1, dtype=float) / path_order
    for radius0, axial0, direction, exit_path, projected_weight in zip(
        exit_r, exit_z, inward_direction, chord, chord_weight
    ):
        ds = exit_path / path_order
        path_s = exit_path * segment_unit
        edge_s = exit_path * edge_unit
        x = radius0 + path_s * direction[0]
        y = path_s * direction[1]
        path_r = np.sqrt(x * x + y * y)
        path_z = axial0 + path_s * direction[2]
        edge_x = radius0 + edge_s * direction[0]
        edge_y = edge_s * direction[1]
        edge_r = np.sqrt(edge_x * edge_x + edge_y * edge_y)
        edge_z = axial0 + edge_s * direction[2]
        emitter_edge = _bilinear(
            field.radial_nodes_m, field.axial_nodes_m,
            field.emitter_density_m3, edge_r, edge_z)
        emitter = 0.5 * (emitter_edge[:-1] + emitter_edge[1:])
        temperature = _bilinear(
            field.radial_nodes_m, field.axial_nodes_m,
            field.gas_temperature_K, path_r, path_z)
        temperature_edge = _bilinear(
            field.radial_nodes_m, field.axial_nodes_m,
            field.gas_temperature_K, edge_r, edge_z)
        absorber = _bilinear(
            field.radial_nodes_m, field.axial_nodes_m,
            field.absorber_density_m3, path_r, path_z)
        sigma = line.center_frequency_hz * np.sqrt(
            BOLTZMANN_CONSTANT_J_K * temperature
            / (line.absorber_mass_kg * SPEED_OF_LIGHT_M_S ** 2))
        profile = voigt_profile(
            detuning[:, None], sigma[None, :], gamma)
        profile_normalization = (
            np.sum(frequency_weight[:, None] * profile, axis=0)
            + profile_tail)
        normalization_minimum = min(
            normalization_minimum, float(np.min(profile_normalization)))
        normalization_maximum = max(
            normalization_maximum, float(np.max(profile_normalization)))
        sigma_edge = line.center_frequency_hz * np.sqrt(
            BOLTZMANN_CONSTANT_J_K * temperature_edge
            / (line.absorber_mass_kg * SPEED_OF_LIGHT_M_S ** 2))
        profile_edge = voigt_profile(
            detuning[:, None], sigma_edge[None, :], gamma)
        profile_normalization_edge = (
            np.sum(frequency_weight[:, None] * profile_edge, axis=0)
            + profile_tail)
        normalization_minimum = min(
            normalization_minimum,
            float(np.min(profile_normalization_edge)))
        normalization_maximum = max(
            normalization_maximum,
            float(np.max(profile_normalization_edge)))
        normalized_emission_edge = (
            emitter_edge[None, :] * profile_edge
            / profile_normalization_edge[None, :])
        normalized_tail_edge = (
            emitter_edge * profile_tail / profile_normalization_edge)

        optical_depth = (
            cross_section_integral * absorber[None, :] * profile * ds)
        attenuation_before = np.exp(-np.column_stack((
            np.zeros(detuning.size),
            np.cumsum(optical_depth[:, :-1], axis=1),
        )))
        attenuation_integral = np.ones_like(optical_depth)
        positive_tau = optical_depth > 0.0
        attenuation_integral[positive_tau] = (
            -np.expm1(-optical_depth[positive_tau])
            / optical_depth[positive_tau])
        linear_moment = np.empty_like(optical_depth)
        small_tau = optical_depth < 1.0e-4
        tau_small = optical_depth[small_tau]
        linear_moment[small_tau] = (
            0.5 - tau_small / 3.0 + tau_small ** 2 / 8.0
            - tau_small ** 3 / 30.0)
        tau_large = optical_depth[~small_tau]
        linear_moment[~small_tau] = (
            1.0 - (1.0 + tau_large) * np.exp(-tau_large)
        ) / (tau_large ** 2)
        emission_near = normalized_emission_edge[:, :-1]
        emission_far = normalized_emission_edge[:, 1:]
        finite_escape_by_frequency = np.sum(
            ds * attenuation_before
            * (
                emission_near * attenuation_integral
                + (emission_far - emission_near) * linear_moment
            ),
            axis=1,
        )
        # The omitted far Lorentz wings are asymptotically optically thin.
        tail_escape = float(
            np.sum(0.5 * (
                normalized_tail_edge[:-1] + normalized_tail_edge[1:]
            )) * ds)
        chord_escaped = (
            float(np.dot(frequency_weight, finite_escape_by_frequency))
            + tail_escape)
        chord_source = float(np.sum(emitter) * ds)
        source_measure += projected_weight * chord_source
        escaped_measure += projected_weight * chord_escaped
        temperature_measure += projected_weight * float(
            np.sum(emitter * temperature) * ds)
        absorber_measure += projected_weight * float(
            np.sum(emitter * absorber) * ds)
        exit_path_measure += projected_weight * float(
            np.sum(emitter * path_s) * ds)

    if source_measure <= 0.0:
        raise ValueError("projected chords do not resolve the declared emitter")
    escape = escaped_measure / source_measure
    return AxisymmetricResonanceEscapeResult(
        escape_probability=escape,
        trapping_factor=1.0 / escape,
        source_weighted_temperature_K=temperature_measure / source_measure,
        source_weighted_absorber_density_m3=(
            absorber_measure / source_measure),
        source_weighted_mean_exit_path_m=(
            exit_path_measure / source_measure),
        frequency_profile_normalization_minimum=normalization_minimum,
        frequency_profile_normalization_maximum=normalization_maximum,
        source_quadrature_order=source_order,
        direction_quadrature_order=direction_order,
        path_quadrature_order=path_order,
        frequency_quadrature_order=frequency_order,
    )


def certify_deterministic_axisymmetric_resonance_escape(
    field: AxisymmetricRadiationMomentField,
    line: ResonanceLineData,
    *,
    additional_lorentz_hwhm_hz: float = 0.0,
    source_quadrature_order: int = 8,
    direction_quadrature_order: int = 8,
    path_quadrature_order: int = 24,
    frequency_quadrature_order: int = 24,
    refinement_factor: int = 2,
    relative_tolerance: float = 0.01,
    maximum_doppler_standard_deviations: float = 1024.0,
) -> AxisymmetricResonanceConvergenceReceipt:
    """Run independent coarse/refined transport and expose the numerical gate.

    The refined solution is the prediction candidate.  Consumers must inspect
    ``converged``; a failed receipt is returned instead of being hidden behind
    a numerically unconverged radiation boundary.
    """
    if (
        int(refinement_factor) != refinement_factor
        or refinement_factor < 2
        or not math.isfinite(relative_tolerance)
        or relative_tolerance <= 0.0
    ):
        raise ValueError("invalid axisymmetric convergence request")
    coarse_orders = tuple(int(order) for order in (
        source_quadrature_order,
        direction_quadrature_order,
        path_quadrature_order,
        frequency_quadrature_order,
    ))
    refined_orders = tuple(
        order * int(refinement_factor) for order in coarse_orders)
    common = {
        "additional_lorentz_hwhm_hz": additional_lorentz_hwhm_hz,
        "maximum_doppler_standard_deviations": (
            maximum_doppler_standard_deviations),
    }
    coarse = deterministic_axisymmetric_resonance_escape(
        field,
        line,
        source_quadrature_order=coarse_orders[0],
        direction_quadrature_order=coarse_orders[1],
        path_quadrature_order=coarse_orders[2],
        frequency_quadrature_order=coarse_orders[3],
        **common,
    )
    refined = deterministic_axisymmetric_resonance_escape(
        field,
        line,
        source_quadrature_order=refined_orders[0],
        direction_quadrature_order=refined_orders[1],
        path_quadrature_order=refined_orders[2],
        frequency_quadrature_order=refined_orders[3],
        **common,
    )
    escape_change = abs(
        refined.escape_probability - coarse.escape_probability
    ) / refined.escape_probability
    trapping_change = abs(
        refined.trapping_factor - coarse.trapping_factor
    ) / refined.trapping_factor
    return AxisymmetricResonanceConvergenceReceipt(
        coarse=coarse,
        refined=refined,
        relative_escape_change=escape_change,
        relative_trapping_change=trapping_change,
        relative_tolerance=relative_tolerance,
        converged=max(escape_change, trapping_change) <= relative_tolerance,
    )
