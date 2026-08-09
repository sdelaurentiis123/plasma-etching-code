"""Deterministic coupled space-frequency resonance-radiation transport.

Projected surface chords integrate all emission positions along each ray.  A
piecewise-constant axisymmetric zone field then produces, at every frequency,
a conservative matrix mapping an emitting zone to the zone of resonant
absorption or to escape.  Natural-Lorentz frequency walking is applied by
zero-padded FFT independently in each zone, and matrix-free GMRES solves the
coupled renewal equation.  This is a deterministic finite-volume analogue of
the absorption/re-emission cycle in Tian (2017), without photon Monte Carlo.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import (
    c as SPEED_OF_LIGHT_M_S,
    k as BOLTZMANN_CONSTANT_J_K,
    physical_constants,
)
from scipy.fft import irfft, next_fast_len, rfft
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.special import voigt_profile

from .axisymmetric_radiation import _axisymmetric_projected_chords
from .geometry import CylindricalReactor
from .vuv_radiation import ResonanceLineData


FULL_ENDCAP_ESCAPE_BOUNDARY_LABELS = (
    "lower_endcap_wafer_plane",
    "upper_endcap",
    "cylindrical_sidewall",
)
PARTIAL_WAFER_ESCAPE_BOUNDARY_LABELS = (
    "lower_endcap_wafer_plane",
    "lower_endcap_outside_wafer",
    "upper_endcap",
    "cylindrical_sidewall",
)


@dataclass(frozen=True)
class AxisymmetricRadiationZoneField:
    """Rectilinear axisymmetric cells aggregated into physical zones."""

    radial_edges_m: np.ndarray
    axial_edges_m: np.ndarray
    cell_zone_index: np.ndarray
    gas_temperature_K: np.ndarray
    absorber_density_m3: np.ndarray
    emitter_density_m3: np.ndarray
    source: str

    def __post_init__(self):
        radial = np.array(self.radial_edges_m, dtype=float, copy=True)
        axial = np.array(self.axial_edges_m, dtype=float, copy=True)
        zones = np.array(self.cell_zone_index, dtype=int, copy=True)
        temperature = np.array(self.gas_temperature_K, dtype=float, copy=True)
        absorber = np.array(self.absorber_density_m3, dtype=float, copy=True)
        emitter = np.array(self.emitter_density_m3, dtype=float, copy=True)
        zone_count = temperature.size
        if (
            radial.ndim != 1
            or axial.ndim != 1
            or min(radial.size, axial.size) < 2
            or radial[0] != 0.0
            or axial[0] != 0.0
            or np.any(~np.isfinite(radial))
            or np.any(~np.isfinite(axial))
            or np.any(np.diff(radial) <= 0.0)
            or np.any(np.diff(axial) <= 0.0)
            or zones.shape != (radial.size - 1, axial.size - 1)
            or temperature.ndim != 1
            or absorber.ndim != 1
            or emitter.ndim != 1
            or zone_count < 1
            or zone_count > 8
            or any(array.shape != (zone_count,) for array in (
                absorber, emitter))
            or np.any(~np.isfinite(temperature))
            or np.any(~np.isfinite(absorber))
            or np.any(~np.isfinite(emitter))
            or np.any(temperature <= 0.0)
            or np.any(absorber < 0.0)
            or np.any(emitter < 0.0)
            or not np.any(emitter > 0.0)
            or np.any(zones < 0)
            or np.any(zones >= zone_count)
            or set(np.unique(zones)) != set(range(zone_count))
            or not str(self.source).strip()
        ):
            raise ValueError("invalid axisymmetric radiation zone field")
        for array in (radial, axial, zones, temperature, absorber, emitter):
            array.setflags(write=False)
        object.__setattr__(self, "radial_edges_m", radial)
        object.__setattr__(self, "axial_edges_m", axial)
        object.__setattr__(self, "cell_zone_index", zones)
        object.__setattr__(self, "gas_temperature_K", temperature)
        object.__setattr__(self, "absorber_density_m3", absorber)
        object.__setattr__(self, "emitter_density_m3", emitter)

    @property
    def geometry(self) -> CylindricalReactor:
        return CylindricalReactor(
            radius_m=float(self.radial_edges_m[-1]),
            length_m=float(self.axial_edges_m[-1]),
        )

    @property
    def zone_count(self) -> int:
        return int(self.gas_temperature_K.size)

    @property
    def analytic_zone_volume_m3(self) -> np.ndarray:
        volume = np.zeros(self.zone_count)
        for radial_index in range(self.radial_edges_m.size - 1):
            annular_area = np.pi * (
                self.radial_edges_m[radial_index + 1] ** 2
                - self.radial_edges_m[radial_index] ** 2)
            for axial_index in range(self.axial_edges_m.size - 1):
                cell_volume = annular_area * (
                    self.axial_edges_m[axial_index + 1]
                    - self.axial_edges_m[axial_index])
                volume[self.cell_zone_index[radial_index, axial_index]] += (
                    cell_volume)
        return volume


@dataclass(frozen=True)
class ZonalPartialRedistributionResult:
    """Coupled spatial/frequency trapping and numerical-conservation receipt."""

    trapping_factor: float
    complete_frequency_redistribution_trapping_factor: float
    initial_emission_zone_probability: np.ndarray
    escape_boundary_labels: tuple[str, ...]
    partial_redistribution_escape_boundary_probability: np.ndarray
    complete_redistribution_escape_boundary_probability: np.ndarray
    partial_redistribution_quench_probability: float
    complete_redistribution_quench_probability: float
    terminal_probability_conservation_error_maximum: float
    zone_source_measure_relative_volume_error_maximum: float
    transition_probability_conservation_error_maximum: float
    frequency_profile_normalization_error_maximum: float
    coherent_frequency_grid_points: int
    coherent_grid_points_per_lorentz_hwhm: float
    coherent_half_range_doppler_standard_deviations: float
    linear_solver_iterations: int
    linear_solver_relative_residual: float
    surface_quadrature_order: int
    direction_quadrature_order: int
    frequency_quadrature_order: int

    def __post_init__(self):
        probability = np.array(
            self.initial_emission_zone_probability, dtype=float, copy=True)
        partial_escape = np.array(
            self.partial_redistribution_escape_boundary_probability,
            dtype=float,
            copy=True,
        )
        complete_escape = np.array(
            self.complete_redistribution_escape_boundary_probability,
            dtype=float,
            copy=True,
        )
        values = np.asarray((
            self.trapping_factor,
            self.complete_frequency_redistribution_trapping_factor,
            self.partial_redistribution_quench_probability,
            self.complete_redistribution_quench_probability,
            self.terminal_probability_conservation_error_maximum,
            self.zone_source_measure_relative_volume_error_maximum,
            self.transition_probability_conservation_error_maximum,
            self.frequency_profile_normalization_error_maximum,
            self.coherent_grid_points_per_lorentz_hwhm,
            self.coherent_half_range_doppler_standard_deviations,
            self.linear_solver_relative_residual,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.trapping_factor < 1.0
            or self.complete_frequency_redistribution_trapping_factor < 1.0
            or probability.ndim != 1
            or probability.size < 1
            or np.any(~np.isfinite(probability))
            or np.any(probability < 0.0)
            or not np.isclose(np.sum(probability), 1.0, atol=2.0e-12)
            or tuple(self.escape_boundary_labels) not in {
                FULL_ENDCAP_ESCAPE_BOUNDARY_LABELS,
                PARTIAL_WAFER_ESCAPE_BOUNDARY_LABELS,
            }
            or partial_escape.shape != (len(self.escape_boundary_labels),)
            or complete_escape.shape != (len(self.escape_boundary_labels),)
            or np.any(~np.isfinite(partial_escape))
            or np.any(~np.isfinite(complete_escape))
            or np.any(partial_escape < 0.0)
            or np.any(complete_escape < 0.0)
            or not 0.0 <= self.partial_redistribution_quench_probability <= 1.0
            or not 0.0 <= self.complete_redistribution_quench_probability <= 1.0
            or not np.isclose(
                np.sum(partial_escape)
                + self.partial_redistribution_quench_probability,
                1.0,
                atol=2.0e-6,
            )
            or not np.isclose(
                np.sum(complete_escape)
                + self.complete_redistribution_quench_probability,
                1.0,
                atol=2.0e-6,
            )
            or not 0.0 <= self.terminal_probability_conservation_error_maximum < 2.0e-6
            or not 0.0 <= self.zone_source_measure_relative_volume_error_maximum < 0.1
            or not 0.0 <= self.transition_probability_conservation_error_maximum < 1.0e-8
            or not 0.0 <= self.frequency_profile_normalization_error_maximum < 1.0e-5
            or int(self.coherent_frequency_grid_points)
            != self.coherent_frequency_grid_points
            or self.coherent_frequency_grid_points < 4096
            or self.coherent_frequency_grid_points % 2
            or self.coherent_grid_points_per_lorentz_hwhm < 4.0
            or self.coherent_half_range_doppler_standard_deviations < 32.0
            or int(self.linear_solver_iterations) != self.linear_solver_iterations
            or self.linear_solver_iterations < 0
            or not 0.0 <= self.linear_solver_relative_residual <= 1.0e-7
            or any(
                int(order) != order or order < 4
                for order in (
                    self.surface_quadrature_order,
                    self.direction_quadrature_order,
                    self.frequency_quadrature_order,
                )
            )
        ):
            raise ValueError("invalid zonal partial-redistribution result")
        probability.setflags(write=False)
        partial_escape.setflags(write=False)
        complete_escape.setflags(write=False)
        object.__setattr__(
            self, "escape_boundary_labels", tuple(self.escape_boundary_labels))
        object.__setattr__(
            self, "initial_emission_zone_probability", probability)
        object.__setattr__(
            self,
            "partial_redistribution_escape_boundary_probability",
            partial_escape,
        )
        object.__setattr__(
            self,
            "complete_redistribution_escape_boundary_probability",
            complete_escape,
        )
        for name in (
            "trapping_factor",
            "complete_frequency_redistribution_trapping_factor",
            "partial_redistribution_quench_probability",
            "complete_redistribution_quench_probability",
            "terminal_probability_conservation_error_maximum",
            "zone_source_measure_relative_volume_error_maximum",
            "transition_probability_conservation_error_maximum",
            "frequency_profile_normalization_error_maximum",
            "coherent_grid_points_per_lorentz_hwhm",
            "coherent_half_range_doppler_standard_deviations",
            "linear_solver_relative_residual",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))
        for name in (
            "coherent_frequency_grid_points", "linear_solver_iterations",
            "surface_quadrature_order", "direction_quadrature_order",
            "frequency_quadrature_order",
        ):
            object.__setattr__(self, name, int(getattr(self, name)))

    @property
    def partial_redistribution_wafer_escape_probability(self) -> float:
        index = self.escape_boundary_labels.index(
            "lower_endcap_wafer_plane")
        return float(
            self.partial_redistribution_escape_boundary_probability[index])

    @property
    def complete_redistribution_wafer_escape_probability(self) -> float:
        index = self.escape_boundary_labels.index(
            "lower_endcap_wafer_plane")
        return float(
            self.complete_redistribution_escape_boundary_probability[index])

    def partial_redistribution_wafer_flux_m2_s(
        self,
        *,
        total_line_emission_rate_s: float,
        wafer_area_m2: float,
    ) -> float:
        """Convert the conservative wafer escape fraction to flux."""
        rate = float(total_line_emission_rate_s)
        area = float(wafer_area_m2)
        if (
            not math.isfinite(rate)
            or rate < 0.0
            or not math.isfinite(area)
            or area <= 0.0
        ):
            raise ValueError("invalid line-emission rate or wafer area")
        return float(
            rate
            * self.partial_redistribution_wafer_escape_probability
            / area)


def _ray_zone_segments(
    field: AxisymmetricRadiationZoneField,
    exit_radius_m: float,
    exit_axial_m: float,
    inward_direction: np.ndarray,
    chord_length_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Split one inward chord at every radial and axial cell boundary."""
    radius0 = float(exit_radius_m)
    axial0 = float(exit_axial_m)
    dx, dy, dz = np.asarray(inward_direction, dtype=float)
    chord = float(chord_length_m)
    intersections = [0.0, chord]
    tolerance = 2.0e-12 * max(field.geometry.radius_m, field.geometry.length_m)
    if abs(dz) > 1.0e-14:
        for boundary in field.axial_edges_m[1:-1]:
            distance = (boundary - axial0) / dz
            if tolerance < distance < chord - tolerance:
                intersections.append(float(distance))
    horizontal_squared = dx * dx + dy * dy
    if horizontal_squared > 1.0e-14:
        for boundary in field.radial_edges_m[1:-1]:
            discriminant = (
                (2.0 * radius0 * dx) ** 2
                - 4.0 * horizontal_squared * (radius0 ** 2 - boundary ** 2)
            )
            if discriminant < 0.0:
                continue
            root = math.sqrt(max(0.0, discriminant))
            for distance in (
                (-2.0 * radius0 * dx - root) / (2.0 * horizontal_squared),
                (-2.0 * radius0 * dx + root) / (2.0 * horizontal_squared),
            ):
                if tolerance < distance < chord - tolerance:
                    intersections.append(float(distance))
    edges = np.asarray(sorted(intersections), dtype=float)
    keep = np.concatenate(([True], np.diff(edges) > tolerance))
    edges = edges[keep]
    midpoint = 0.5 * (edges[:-1] + edges[1:])
    x = radius0 + midpoint * dx
    y = midpoint * dy
    radius = np.sqrt(x * x + y * y)
    axial = axial0 + midpoint * dz
    radial_index = np.clip(
        np.searchsorted(field.radial_edges_m, radius, side="right") - 1,
        0, field.radial_edges_m.size - 2)
    axial_index = np.clip(
        np.searchsorted(field.axial_edges_m, axial, side="right") - 1,
        0, field.axial_edges_m.size - 2)
    zones = field.cell_zone_index[radial_index, axial_index]
    length = np.diff(edges)
    # Merge adjacent rectilinear cells that belong to the same aggregate zone.
    merged_zones = [int(zones[0])]
    merged_lengths = [float(length[0])]
    for zone, segment_length in zip(zones[1:], length[1:]):
        if int(zone) == merged_zones[-1]:
            merged_lengths[-1] += float(segment_length)
        else:
            merged_zones.append(int(zone))
            merged_lengths.append(float(segment_length))
    return np.asarray(merged_zones, dtype=int), np.asarray(merged_lengths)


def _attenuation_integral(optical_depth: np.ndarray) -> np.ndarray:
    output = np.ones_like(optical_depth)
    positive = optical_depth > 0.0
    output[positive] = -np.expm1(-optical_depth[positive]) / optical_depth[positive]
    return output


def _zonal_transition_probabilities(
    field: AxisymmetricRadiationZoneField,
    absorption_coefficient_m_inv: np.ndarray,
    *,
    wafer_radius_m: float,
    escape_boundary_labels: tuple[str, ...],
    surface_quadrature_order: int,
    direction_quadrature_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Return absorption and boundary-resolved escape transition ledgers."""
    zone_count, frequency_count = absorption_coefficient_m_inv.shape
    exit_r, exit_z, direction, chord, weight = _axisymmetric_projected_chords(
        field.geometry,
        surface_order=surface_quadrature_order,
        direction_order=direction_quadrature_order,
        endcap_radial_break_m=(
            wafer_radius_m
            if len(escape_boundary_labels) == 4 else None
        ),
    )
    absorption_measure = np.zeros((zone_count, zone_count, frequency_count))
    boundary_count = len(escape_boundary_labels)
    escape_measure = np.zeros((
        zone_count, boundary_count, frequency_count))
    source_measure = np.zeros(zone_count)
    boundary_source_measure = np.zeros((zone_count, boundary_count))
    conservation_error = 0.0
    for radius0, axial0, ray_direction, ray_length, ray_weight in zip(
        exit_r, exit_z, direction, chord, weight
    ):
        if axial0 == 0.0:
            if boundary_count == 3 or radius0 <= wafer_radius_m:
                boundary_index = 0
            else:
                boundary_index = 1
        elif axial0 == field.geometry.length_m:
            boundary_index = boundary_count - 2
        else:
            boundary_index = boundary_count - 1
        zones, lengths = _ray_zone_segments(
            field, radius0, axial0, ray_direction, ray_length)
        optical_depth = (
            absorption_coefficient_m_inv[zones]
            * lengths[:, None])
        prefix = np.vstack((
            np.zeros(frequency_count),
            np.cumsum(optical_depth, axis=0),
        ))
        for source_segment, (source_zone, source_length) in enumerate(
            zip(zones, lengths)
        ):
            source_measure[source_zone] += ray_weight * source_length
            boundary_source_measure[source_zone, boundary_index] += (
                ray_weight * source_length)
            source_integral = (
                source_length
                * _attenuation_integral(optical_depth[source_segment]))
            escaped = (
                np.exp(-prefix[source_segment]) * source_integral)
            escape_measure[source_zone, boundary_index] += (
                ray_weight * escaped)
            same_segment_absorption = source_length - source_integral
            absorption_measure[source_zone, source_zone] += (
                ray_weight * same_segment_absorption)
            for absorbing_segment in range(source_segment):
                intervening_depth = (
                    prefix[source_segment]
                    - prefix[absorbing_segment + 1])
                absorbed = (
                    source_integral * np.exp(-intervening_depth)
                    * -np.expm1(-optical_depth[absorbing_segment]))
                absorption_measure[
                    source_zone, zones[absorbing_segment]
                ] += ray_weight * absorbed
    if np.any(source_measure <= 0.0):
        raise RuntimeError("projected chords did not resolve every zone")
    absorption_probability = (
        absorption_measure / source_measure[:, None, None])
    escape_probability = escape_measure / source_measure[:, None, None]
    row_sum = (
        np.sum(absorption_probability, axis=1)
        + np.sum(escape_probability, axis=1))
    conservation_error = float(np.max(np.abs(row_sum - 1.0)))
    if conservation_error > 1.0e-8:
        raise RuntimeError(
            "zonal radiation transition ledger does not conserve probability")
    return (
        absorption_probability,
        escape_probability,
        source_measure,
        boundary_source_measure,
        conservation_error,
    )


def deterministic_zonal_partial_redistribution(
    field: AxisymmetricRadiationZoneField,
    line: ResonanceLineData,
    *,
    velocity_changing_collision_frequency_s_inv: float | np.ndarray = 0.0,
    quenching_collision_frequency_s_inv: float | np.ndarray = 0.0,
    additional_lorentz_hwhm_hz: float = 0.0,
    wafer_radius_m: float | None = None,
    surface_quadrature_order: int = 8,
    direction_quadrature_order: int = 8,
    frequency_quadrature_order: int = 24,
    coherent_half_range_doppler_standard_deviations: float = 32.0,
    coherent_frequency_grid_points: int | None = None,
    coherent_grid_points_per_lorentz_hwhm: float = 8.0,
    linear_solver_relative_tolerance: float = 5.0e-9,
) -> ZonalPartialRedistributionResult:
    """Solve coupled spatial-zone and partial-frequency redistribution."""
    if not isinstance(field, AxisymmetricRadiationZoneField):
        raise TypeError("an axisymmetric radiation zone field is required")
    if not isinstance(line, ResonanceLineData):
        raise TypeError("resonance-line atomic data are required")
    zone_count = field.zone_count
    wafer_radius = (
        field.geometry.radius_m
        if wafer_radius_m is None else float(wafer_radius_m))
    escape_boundary_labels = (
        FULL_ENDCAP_ESCAPE_BOUNDARY_LABELS
        if wafer_radius == field.geometry.radius_m
        else PARTIAL_WAFER_ESCAPE_BOUNDARY_LABELS)
    velocity = np.broadcast_to(np.asarray(
        velocity_changing_collision_frequency_s_inv, dtype=float),
        (zone_count,)).copy()
    quench = np.broadcast_to(np.asarray(
        quenching_collision_frequency_s_inv, dtype=float),
        (zone_count,)).copy()
    gamma = line.natural_lorentz_hwhm_hz + float(
        additional_lorentz_hwhm_hz)
    surface_order = int(surface_quadrature_order)
    direction_order = int(direction_quadrature_order)
    frequency_order = int(frequency_quadrature_order)
    half_range = float(coherent_half_range_doppler_standard_deviations)
    requested_points_per_hwhm = float(
        coherent_grid_points_per_lorentz_hwhm)
    solver_tolerance = float(linear_solver_relative_tolerance)
    if (
        np.any(~np.isfinite(velocity))
        or np.any(velocity < 0.0)
        or np.any(~np.isfinite(quench))
        or np.any(quench < 0.0)
        or not math.isfinite(gamma)
        or gamma <= 0.0
        or not math.isfinite(wafer_radius)
        or not 0.0 < wafer_radius <= field.geometry.radius_m
        or any(
            integer != value or integer < 4
            for integer, value in (
                (surface_order, surface_quadrature_order),
                (direction_order, direction_quadrature_order),
                (frequency_order, frequency_quadrature_order),
            )
        )
        or not math.isfinite(half_range)
        or half_range < 32.0
        or not math.isfinite(requested_points_per_hwhm)
        or requested_points_per_hwhm < 4.0
        or not math.isfinite(solver_tolerance)
        or not 0.0 < solver_tolerance <= 1.0e-7
    ):
        raise ValueError("invalid zonal partial-redistribution condition")

    sigma = np.asarray([
        line.doppler_standard_deviation_hz(temperature)
        for temperature in field.gas_temperature_K
    ])
    sigma_reference = float(np.max(sigma))
    natural_gamma_ratio = line.natural_lorentz_hwhm_hz / sigma_reference
    required_grid_points = max(4096, int(math.ceil(
        2.0 * half_range * requested_points_per_hwhm
        / natural_gamma_ratio)))
    if coherent_frequency_grid_points is None:
        grid_points = 1 << (required_grid_points - 1).bit_length()
    else:
        if (
            int(coherent_frequency_grid_points)
            != coherent_frequency_grid_points
            or int(coherent_frequency_grid_points) < 4096
            or int(coherent_frequency_grid_points) % 2
        ):
            raise ValueError("invalid coherent frequency grid")
        grid_points = int(coherent_frequency_grid_points)
    if grid_points > 1_048_576:
        raise ValueError("coherent frequency grid exceeds memory bound")
    grid_spacing = 2.0 * half_range / grid_points
    actual_points_per_hwhm = natural_gamma_ratio / grid_spacing
    if actual_points_per_hwhm < 4.0:
        raise ValueError("coherent natural linewidth is under-resolved")
    x = -half_range + (np.arange(grid_points) + 0.5) * grid_spacing

    profile = np.asarray([
        voigt_profile(x * sigma_reference, local_sigma, gamma)
        * sigma_reference
        for local_sigma in sigma
    ])
    tail = 2.0 * gamma / (np.pi * sigma_reference * half_range)
    profile_normalization = np.sum(profile, axis=1) * grid_spacing + tail
    normalized_profile = profile / profile_normalization[:, None]
    finite_profile_probability = (
        np.sum(normalized_profile, axis=1) * grid_spacing)
    normalized_tail_probability = np.clip(
        1.0 - finite_profile_probability, 0.0, 1.0)
    profile_normalization_error = float(np.max(
        np.abs(profile_normalization - 1.0)))

    minimum_escape_intervals = max(4096, frequency_order * 256)
    escape_intervals = 1 << (minimum_escape_intervals - 1).bit_length()
    escape_x = np.linspace(-half_range, half_range, escape_intervals + 1)
    classical_radius = physical_constants["classical electron radius"][0]
    integrated_cross_section = (
        np.pi * classical_radius * SPEED_OF_LIGHT_M_S
        * line.absorption_oscillator_strength)
    escape_profile = np.asarray([
        voigt_profile(
            escape_x * sigma_reference, local_sigma, gamma)
        for local_sigma in sigma
    ])
    escape_absorption = (
        field.absorber_density_m3[:, None]
        * integrated_cross_section * escape_profile)
    (
        coarse_transition,
        coarse_escape,
        source_measure,
        boundary_source_measure,
        conservation_error,
    ) = (
        _zonal_transition_probabilities(
            field,
            escape_absorption,
            wafer_radius_m=wafer_radius,
            escape_boundary_labels=escape_boundary_labels,
            surface_quadrature_order=surface_order,
            direction_quadrature_order=direction_order,
        ))
    transition = np.empty((zone_count, zone_count, grid_points))
    direct_escape = np.empty((
        zone_count, len(escape_boundary_labels), grid_points))
    for source_zone in range(zone_count):
        for absorbing_zone in range(zone_count):
            transition[source_zone, absorbing_zone] = np.interp(
                x,
                escape_x,
                coarse_transition[source_zone, absorbing_zone],
            )
        for boundary_index in range(len(escape_boundary_labels)):
            direct_escape[source_zone, boundary_index] = np.interp(
                x,
                escape_x,
                coarse_escape[source_zone, boundary_index],
            )
    interpolated_conservation_error = float(np.max(np.abs(
        np.sum(transition, axis=1)
        + np.sum(direct_escape, axis=1)
        - 1.0
    )))
    conservation_error = max(
        conservation_error, interpolated_conservation_error)
    if conservation_error > 1.0e-8:
        raise RuntimeError(
            "interpolated zonal radiation ledger does not conserve probability")

    analytic_volume = field.analytic_zone_volume_m3
    relative_measure = source_measure / analytic_volume
    relative_volume_error = float(
        np.max(np.abs(relative_measure / np.mean(relative_measure) - 1.0)))
    clear_boundary_probability = (
        boundary_source_measure / source_measure[:, None])
    # The projected-ray measure is used only to normalize each transition row.
    # Initial source fractions have an exact analytic finite-volume measure,
    # avoiding any source bias when a quadrature node lies near a zone edge.
    initial_weight = field.emitter_density_m3 * analytic_volume
    initial_probability = initial_weight / np.sum(initial_weight)

    redistribution_probability = velocity / (
        line.transition_probability_s_inv + velocity)
    coherent_probability = 1.0 - redistribution_probability
    survival_probability = line.transition_probability_s_inv / (
        line.transition_probability_s_inv + quench)

    difference = (
        np.arange(grid_points, dtype=float) - grid_points // 2
    ) * grid_spacing
    coherent_kernel = natural_gamma_ratio / np.pi / (
        difference * difference + natural_gamma_ratio ** 2)
    fft_points = next_fast_len(2 * grid_points - 1)
    kernel_transform = rfft(coherent_kernel, fft_points)
    convolution_start = grid_points // 2

    def coherent_convolution(values: np.ndarray) -> np.ndarray:
        transformed = rfft(values, fft_points, axis=1)
        full = irfft(
            transformed * kernel_transform[None, :],
            fft_points,
            axis=1,
        )
        return (
            full[:, convolution_start:convolution_start + grid_points]
            * grid_spacing)

    cycle_rhs = np.einsum(
        "ijf,j->if", transition, survival_probability, optimize=True)

    def renewal_matvec(flat_values: np.ndarray) -> np.ndarray:
        values = np.asarray(flat_values, dtype=float).reshape(
            zone_count, grid_points)
        redistributed_mean = np.einsum(
            "if,if->i", normalized_profile, values) * grid_spacing
        next_zone = survival_probability[:, None] * (
            coherent_probability[:, None] * coherent_convolution(values)
            + redistribution_probability[:, None]
            * redistributed_mean[:, None]
        )
        mapped = np.einsum(
            "ijf,jf->if", transition, next_zone, optimize=True)
        return (values - mapped).ravel()

    size = zone_count * grid_points
    operator = LinearOperator(
        (size, size), matvec=renewal_matvec, dtype=float)
    diagonal = 1.0 - np.asarray([
        transition[zone, zone]
        * survival_probability[zone]
        * coherent_probability[zone]
        * coherent_kernel[grid_points // 2] * grid_spacing
        for zone in range(zone_count)
    ])
    if np.any(diagonal <= 0.0):
        raise RuntimeError("zonal preconditioner is not positive")
    preconditioner = LinearOperator(
        (size, size),
        matvec=lambda vector: (
            np.asarray(vector).reshape(zone_count, grid_points)
            / diagonal
        ).ravel(),
        dtype=float,
    )
    solve_iterations: list[int] = []
    solve_residuals: list[float] = []

    def solve_renewal(rhs: np.ndarray, *, label: str) -> np.ndarray:
        rhs = np.asarray(rhs, dtype=float)
        if rhs.shape != (zone_count, grid_points):
            raise RuntimeError("zonal renewal right-hand side shape mismatch")
        rhs_flat = rhs.ravel()
        rhs_norm = float(np.linalg.norm(rhs_flat))
        iteration_residuals: list[float] = []
        if rhs_norm == 0.0:
            solve_iterations.append(0)
            solve_residuals.append(0.0)
            return np.zeros_like(rhs)
        solution, info = gmres(
            operator,
            rhs_flat,
            M=preconditioner,
            rtol=solver_tolerance,
            atol=min(1.0e-12, solver_tolerance * rhs_norm),
            restart=100,
            maxiter=500,
            callback=iteration_residuals.append,
            callback_type="pr_norm",
        )
        if info != 0:
            raise RuntimeError(
                "zonal partial-redistribution GMRES failed for "
                f"{label} (info={info})")
        relative_residual = float(
            np.linalg.norm(renewal_matvec(solution) - rhs_flat) / rhs_norm)
        if relative_residual > max(1.0e-7, 10.0 * solver_tolerance):
            raise RuntimeError(
                f"zonal GMRES residual exceeds tolerance for {label}")
        solve_iterations.append(len(iteration_residuals))
        solve_residuals.append(relative_residual)
        return solution.reshape(zone_count, grid_points)

    excess = solve_renewal(cycle_rhs, label="expected-cycle ledger")

    expected_by_initial_zone = (
        1.0 + np.einsum(
            "if,if->i", normalized_profile, excess) * grid_spacing)
    trapping = max(1.0, float(np.dot(
        initial_probability, expected_by_initial_zone)))

    # A coherent natural-line walk or a local Voigt redraw can leave the
    # finite frequency grid.  At >=32 Doppler widths the omitted photon is
    # optically thin.  Allocate that probability to the exact clear-geometry
    # boundary view rather than silently discarding it.
    coherent_in_grid_probability = coherent_convolution(
        np.ones((zone_count, grid_points)))
    coherent_tail_probability = np.clip(
        1.0 - coherent_in_grid_probability, 0.0, 1.0)
    local_reemission_tail_probability = (
        coherent_probability[:, None] * coherent_tail_probability
        + redistribution_probability[:, None]
        * normalized_tail_probability[:, None]
    )
    absorbed_then_far_escape = np.einsum(
        "ijf,j,jf,jb->ibf",
        transition,
        survival_probability,
        local_reemission_tail_probability,
        clear_boundary_probability,
        optimize=True,
    )
    partial_escape_by_boundary = np.empty(len(escape_boundary_labels))
    for boundary_index, boundary_label in enumerate(escape_boundary_labels):
        terminal = solve_renewal(
            direct_escape[:, boundary_index]
            + absorbed_then_far_escape[:, boundary_index],
            label=f"{boundary_label} terminal ledger",
        )
        terminal_by_initial_zone = (
            np.einsum(
                "if,if->i", normalized_profile, terminal) * grid_spacing
            + normalized_tail_probability
            * clear_boundary_probability[:, boundary_index]
        )
        partial_escape_by_boundary[boundary_index] = np.dot(
            initial_probability, terminal_by_initial_zone)
    partial_quench_rhs = np.einsum(
        "ijf,j->if",
        transition,
        1.0 - survival_probability,
        optimize=True,
    )
    partial_quench_terminal = solve_renewal(
        partial_quench_rhs, label="quench terminal ledger")
    partial_quench_by_initial_zone = np.einsum(
        "if,if->i", normalized_profile, partial_quench_terminal
    ) * grid_spacing
    partial_quench_probability = float(np.dot(
        initial_probability, partial_quench_by_initial_zone))

    # Independent complete-redistribution limit of the same spatial ledger.
    profile_weighted_transition = np.einsum(
        "if,ijf->ij", normalized_profile, transition) * grid_spacing
    complete_operator = (
        np.eye(zone_count)
        - profile_weighted_transition * survival_probability[None, :])
    complete_by_zone = np.linalg.solve(
        complete_operator, np.ones(zone_count))
    complete_trapping = max(1.0, float(np.dot(
        initial_probability, complete_by_zone)))
    profile_weighted_direct_escape = np.einsum(
        "if,ibf->ib", normalized_profile, direct_escape
    ) * grid_spacing
    complete_direct_escape = (
        profile_weighted_direct_escape
        + normalized_tail_probability[:, None] * clear_boundary_probability)
    complete_escape_by_zone = np.linalg.solve(
        complete_operator, complete_direct_escape)
    complete_escape_by_boundary = np.einsum(
        "i,ib->b", initial_probability, complete_escape_by_zone)
    complete_quench_rhs = np.einsum(
        "ij,j->i",
        profile_weighted_transition,
        1.0 - survival_probability,
    )
    complete_quench_by_zone = np.linalg.solve(
        complete_operator, complete_quench_rhs)
    complete_quench_probability = float(np.dot(
        initial_probability, complete_quench_by_zone))
    partial_terminal_error = abs(
        float(np.sum(partial_escape_by_boundary))
        + partial_quench_probability - 1.0)
    complete_terminal_error = abs(
        float(np.sum(complete_escape_by_boundary))
        + complete_quench_probability - 1.0)
    terminal_conservation_error = max(
        partial_terminal_error, complete_terminal_error)
    if terminal_conservation_error >= 2.0e-6:
        raise RuntimeError(
            "zonal terminal photon ledger does not conserve probability")

    return ZonalPartialRedistributionResult(
        trapping_factor=trapping,
        complete_frequency_redistribution_trapping_factor=(
            complete_trapping),
        initial_emission_zone_probability=initial_probability,
        escape_boundary_labels=escape_boundary_labels,
        partial_redistribution_escape_boundary_probability=(
            partial_escape_by_boundary),
        complete_redistribution_escape_boundary_probability=(
            complete_escape_by_boundary),
        partial_redistribution_quench_probability=(
            partial_quench_probability),
        complete_redistribution_quench_probability=(
            complete_quench_probability),
        terminal_probability_conservation_error_maximum=(
            terminal_conservation_error),
        zone_source_measure_relative_volume_error_maximum=(
            relative_volume_error),
        transition_probability_conservation_error_maximum=(
            conservation_error),
        frequency_profile_normalization_error_maximum=(
            profile_normalization_error),
        coherent_frequency_grid_points=grid_points,
        coherent_grid_points_per_lorentz_hwhm=(
            actual_points_per_hwhm),
        coherent_half_range_doppler_standard_deviations=half_range,
        linear_solver_iterations=max(solve_iterations),
        linear_solver_relative_residual=max(solve_residuals),
        surface_quadrature_order=surface_order,
        direction_quadrature_order=direction_order,
        frequency_quadrature_order=frequency_order,
    )
