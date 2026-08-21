"""Measured-waveform SPTS Bosch reactor with cylindrical 3-D wafer transfer."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import MappingProxyType
import math

import numpy as np

from .bosch_spts_reduced import (
    BoschSPTSReducedParameters, BoschSPTSReducedReactorSolution,
    solve_bosch_spts_reduced_reactor,
)
from .cylindrical_inventory_lift import (
    CylindricalFiniteVolumeGrid,
    DeterministicCylindricalIndependentInventoryLift,
    normalized_cylindrical_annular_skin_source,
)
from .geometry import CylindricalReactor


_SPECIES = ("F", "C4F8_film_precursor", "positive_ion")
_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(128)
_LOG_TWO = math.log(2.0)
_TRANSMISSION_COEFFICIENT_BOUND = 0.1
_TRANSMISSION_VPP_REFERENCE_V = 637.4409584828442
_TRANSMISSION_VPP_SCALE_V = 3.816305957878358
_TRANSMISSION_VPP_DOMAIN_V = (626.9533265149638, 643.534317555529)


@dataclass(frozen=True)
class BoschSPTSCylindricalParameters:
    reduced: BoschSPTSReducedParameters = field(
        default_factory=BoschSPTSReducedParameters)
    azimuthal_cell_count: int = 16
    species_source_ring_radius_m: tuple[float, float, float] | None = None
    species_source_radial_width_m: tuple[float, float, float] | None = None
    ion_edge_focus_amplitude: float = 0.0
    ion_edge_focus_onset_radius_m: float = 0.095
    ion_edge_focus_width_m: float = 0.005
    source_cosine_coefficients: tuple[tuple[float, ...], ...] = ((), (), ())
    source_sine_coefficients: tuple[tuple[float, ...], ...] = ((), (), ())

    def __post_init__(self):
        if not isinstance(self.reduced, BoschSPTSReducedParameters):
            raise TypeError("reduced parameters must be BoschSPTSReducedParameters")
        if (int(self.azimuthal_cell_count) != self.azimuthal_cell_count
                or self.azimuthal_cell_count < 8):
            raise ValueError("cylindrical Bosch transfer requires at least 8 phi cells")
        ring = self.species_source_ring_radius_m
        if ring is None:
            ring = (self.reduced.source_ring_radius_m,) * len(_SPECIES)
        ring = tuple(float(value) for value in ring)
        width = self.species_source_radial_width_m
        if width is None:
            width = (self.reduced.source_radial_width_m,) * len(_SPECIES)
        width = tuple(float(value) for value in width)
        if (len(ring) != len(_SPECIES) or len(width) != len(_SPECIES)
                or any(not math.isfinite(value)
                       or not 0.0 <= value <= self.reduced.reactor_radius_m
                       for value in ring)
                or any(not math.isfinite(value) or value <= 0.0
                       for value in width)):
            raise ValueError("invalid species-resolved radial source moments")
        if (not math.isfinite(self.ion_edge_focus_amplitude)
                or not 0.0 <= self.ion_edge_focus_amplitude <= 5.0
                or not math.isfinite(self.ion_edge_focus_onset_radius_m)
                or not 0.0 <= self.ion_edge_focus_onset_radius_m <= (
                    self.reduced.wafer_radius_m)
                or not math.isfinite(self.ion_edge_focus_width_m)
                or self.ion_edge_focus_width_m <= 0.0):
            raise ValueError("invalid ion edge-focus parameters")
        cosine = tuple(tuple(float(value) for value in row)
                       for row in self.source_cosine_coefficients)
        sine = tuple(tuple(float(value) for value in row)
                     for row in self.source_sine_coefficients)
        if (len(cosine) != len(_SPECIES) or len(sine) != len(_SPECIES)
                or any(len(c) != len(s) or len(c) > 4
                       for c, s in zip(cosine, sine))
                or any(not math.isfinite(value) or abs(value) > 5.0
                       for row in cosine + sine for value in row)):
            raise ValueError("invalid species-resolved source harmonics")
        object.__setattr__(self, "source_cosine_coefficients", cosine)
        object.__setattr__(self, "source_sine_coefficients", sine)
        object.__setattr__(self, "species_source_ring_radius_m", ring)
        object.__setattr__(self, "species_source_radial_width_m", width)

    def manifest(self):
        return {
            "schema": "petch-spts-bosch-cylindrical-parameters-v1",
            "reduced": self.reduced.manifest(),
            "azimuthal_cell_count": self.azimuthal_cell_count,
            "species_source_ring_radius_m": list(
                self.species_source_ring_radius_m),
            "species_source_radial_width_m": list(
                self.species_source_radial_width_m),
            "ion_edge_focus_amplitude": self.ion_edge_focus_amplitude,
            "ion_edge_focus_onset_radius_m": (
                self.ion_edge_focus_onset_radius_m),
            "ion_edge_focus_width_m": self.ion_edge_focus_width_m,
            "source_cosine_coefficients": [
                list(row) for row in self.source_cosine_coefficients],
            "source_sine_coefficients": [
                list(row) for row in self.source_sine_coefficients],
            "azimuthal_source_form": (
                "positive exponential Fourier modulation of each species source"),
        }


def bosch_real_zernike_modes(maximum_order):
    """Return the frozen complete real-Zernike ordering without a piston."""
    maximum_order = int(maximum_order)
    if not 0 <= maximum_order <= 10:
        raise ValueError("Bosch Zernike maximum order must be in [0, 10]")
    modes = []
    for radial_order in range(1, maximum_order + 1):
        for azimuthal_order in range(
                radial_order % 2, radial_order + 1, 2):
            if azimuthal_order == 0:
                modes.append((radial_order, azimuthal_order, "cos"))
            else:
                modes.extend((
                    (radial_order, azimuthal_order, "cos"),
                    (radial_order, azimuthal_order, "sin"),
                ))
    return tuple(modes)


def _zernike_radial(radial_order, azimuthal_order, radius):
    radius = np.asarray(radius, dtype=float)
    output = np.zeros_like(radius)
    for index in range((radial_order - azimuthal_order) // 2 + 1):
        output += (
            (-1) ** index * math.factorial(radial_order - index)
            / (
                math.factorial(index)
                * math.factorial(
                    (radial_order + azimuthal_order) // 2 - index)
                * math.factorial(
                    (radial_order - azimuthal_order) // 2 - index)
            )
            * radius ** (radial_order - 2 * index)
        )
    return output


def bosch_real_zernike_design(maximum_order, radius_fraction, phi_rad):
    """Evaluate the complete frozen real-Zernike basis on the unit disk."""
    radius = np.asarray(radius_fraction, dtype=float)
    phi = np.asarray(phi_rad, dtype=float)
    radius, phi = np.broadcast_arrays(radius, phi)
    if (
        np.any(~np.isfinite(radius))
        or np.any(~np.isfinite(phi))
        or np.any(radius < 0.0)
        or np.any(radius > 1.0 + 2.0e-14)
    ):
        raise ValueError("Bosch Zernike coordinates must lie on the unit disk")
    radius = np.minimum(radius, 1.0)
    columns = []
    for radial_order, azimuthal_order, phase in bosch_real_zernike_modes(
            maximum_order):
        radial = _zernike_radial(
            radial_order, azimuthal_order, radius)
        if phase == "cos":
            columns.append(
                radial * np.cos(azimuthal_order * phi))
        else:
            columns.append(
                radial * np.sin(azimuthal_order * phi))
    if not columns:
        return np.empty(radius.shape + (0,), dtype=float)
    return np.stack(columns, axis=-1)


@dataclass(frozen=True)
class BoschSPTSWaferIonTransmissionLaw:
    """Positive current-conserving tool fingerprint on the ion boundary."""

    static_maximum_order: int
    static_coefficients: tuple[float, ...]
    dynamic_maximum_order: int = 0
    dynamic_coefficients: tuple[float, ...] = ()
    vpp_reference_V: float = _TRANSMISSION_VPP_REFERENCE_V
    vpp_scale_V: float = _TRANSMISSION_VPP_SCALE_V
    vpp_domain_V: tuple[float, float] = _TRANSMISSION_VPP_DOMAIN_V
    coefficient_bound: float = _TRANSMISSION_COEFFICIENT_BOUND
    maximum_absolute_log_field: float = _LOG_TWO

    def __post_init__(self):
        static_order = int(self.static_maximum_order)
        dynamic_order = int(self.dynamic_maximum_order)
        static = tuple(float(value) for value in self.static_coefficients)
        dynamic = tuple(float(value) for value in self.dynamic_coefficients)
        voltage_domain = tuple(float(value) for value in self.vpp_domain_V)
        if (
            static_order != self.static_maximum_order
            or not 1 <= static_order <= 10
            or dynamic_order != self.dynamic_maximum_order
            or dynamic_order not in (0, 2)
            or len(static) != len(bosch_real_zernike_modes(static_order))
            or len(dynamic) != len(bosch_real_zernike_modes(dynamic_order))
            or len(voltage_domain) != 2
            or not voltage_domain[0] < voltage_domain[1]
            or not math.isfinite(self.vpp_reference_V)
            or not math.isfinite(self.vpp_scale_V)
            or self.vpp_scale_V <= 0.0
            or not math.isfinite(self.coefficient_bound)
            or self.coefficient_bound <= 0.0
            or not math.isfinite(self.maximum_absolute_log_field)
            or self.maximum_absolute_log_field <= 0.0
            or any(not math.isfinite(value) for value in static + dynamic)
            or any(abs(value) > self.coefficient_bound
                   for value in static + dynamic)
        ):
            raise ValueError("invalid Bosch wafer ion-transmission law")
        object.__setattr__(self, "static_maximum_order", static_order)
        object.__setattr__(self, "dynamic_maximum_order", dynamic_order)
        object.__setattr__(self, "static_coefficients", static)
        object.__setattr__(self, "dynamic_coefficients", dynamic)
        object.__setattr__(self, "vpp_domain_V", voltage_domain)
        maximum = self._certification_maximum_absolute_log_field()
        if maximum > self.maximum_absolute_log_field + 2.0e-14:
            raise ValueError(
                "Bosch ion-transmission log field exceeds its frozen bound")

    def standardized_vpp(self, vpp_rms_V):
        value = float(vpp_rms_V)
        if (
            not math.isfinite(value)
            or value < self.vpp_domain_V[0]
            or value > self.vpp_domain_V[1]
        ):
            raise ValueError("C4F8 platen Vpp RMS lies outside the frozen domain")
        return (value - self.vpp_reference_V) / self.vpp_scale_V

    def log_field(self, radius_fraction, phi_rad, standardized_vpp):
        standardized_vpp = float(standardized_vpp)
        if not math.isfinite(standardized_vpp):
            raise ValueError("standardized Vpp must be finite")
        static = bosch_real_zernike_design(
            self.static_maximum_order, radius_fraction, phi_rad)
        output = np.einsum(
            "...k,k->...", static, np.asarray(self.static_coefficients))
        if self.dynamic_coefficients:
            dynamic = bosch_real_zernike_design(
                self.dynamic_maximum_order, radius_fraction, phi_rad)
            output = output + standardized_vpp * np.einsum(
                "...k,k->...", dynamic,
                np.asarray(self.dynamic_coefficients))
        if np.any(~np.isfinite(output)):
            raise RuntimeError("Bosch ion-transmission field is nonfinite")
        return output

    def _certification_maximum_absolute_log_field(self):
        radius = np.linspace(0.0, 1.0, 129)
        phi = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
        rho, angle = np.meshgrid(radius, phi, indexing="ij")
        values = []
        for voltage in self.vpp_domain_V:
            values.append(np.max(np.abs(self.log_field(
                rho, angle, self.standardized_vpp(voltage)))))
        return float(max(values, default=0.0))

    def manifest(self):
        return {
            "schema": "petch-spts-bosch-wafer-ion-transmission-law-v1",
            "basis": "complete real Zernike basis excluding piston",
            "static_maximum_order": self.static_maximum_order,
            "static_modes": [list(item) for item in bosch_real_zernike_modes(
                self.static_maximum_order)],
            "static_coefficients": list(self.static_coefficients),
            "dynamic_maximum_order": self.dynamic_maximum_order,
            "dynamic_modes": [list(item) for item in bosch_real_zernike_modes(
                self.dynamic_maximum_order)],
            "dynamic_coefficients": list(self.dynamic_coefficients),
            "dynamic_input": "standardized C4F8 platen peak-to-peak RMS",
            "vpp_reference_V": self.vpp_reference_V,
            "vpp_scale_V": self.vpp_scale_V,
            "vpp_domain_V": list(self.vpp_domain_V),
            "coefficient_bound": self.coefficient_bound,
            "maximum_absolute_log_field": self.maximum_absolute_log_field,
            "certified_maximum_absolute_log_field": (
                self._certification_maximum_absolute_log_field()),
            "positive_form": "exp(static Zernike field + z*dynamic field)",
            "current_normalization": (
                "baseline-ion-current-weighted finite-volume wafer integral"),
            "target_depth_used": False,
        }


@dataclass(frozen=True)
class BoschSPTSCylindricalSourceResponse:
    species_names: tuple[str, ...]
    x_m: np.ndarray
    y_m: np.ndarray
    unit_point_flux_per_density_m_s: np.ndarray
    unit_wafer_average_flux_per_density_m_s: np.ndarray
    maximum_species_ledger_relative_residual: float
    maximum_linear_system_relative_residual: float
    provenance: dict

    def __post_init__(self):
        names = tuple(self.species_names)
        x = np.asarray(self.x_m, dtype=float).copy()
        y = np.asarray(self.y_m, dtype=float).copy()
        point = np.asarray(self.unit_point_flux_per_density_m_s, dtype=float).copy()
        average = np.asarray(
            self.unit_wafer_average_flux_per_density_m_s, dtype=float).copy()
        if (
            names != _SPECIES
            or x.ndim != 1 or y.shape != x.shape or x.size < 4
            or point.shape != (len(names), x.size)
            or average.shape != (len(names),)
            or any(np.any(~np.isfinite(value)) for value in (x, y, point, average))
            or np.any(point < 0.0) or np.any(average < 0.0)
            or len(set(zip(x, y))) != x.size
            or not 0.0 <= self.maximum_species_ledger_relative_residual < 1.0e-8
            or not 0.0 <= self.maximum_linear_system_relative_residual < 1.0e-10
        ):
            raise ValueError("invalid cylindrical Bosch source response")
        for value in (x, y, point, average):
            value.setflags(write=False)
        object.__setattr__(self, "species_names", names)
        object.__setattr__(self, "x_m", x)
        object.__setattr__(self, "y_m", y)
        object.__setattr__(self, "unit_point_flux_per_density_m_s", point)
        object.__setattr__(
            self, "unit_wafer_average_flux_per_density_m_s", average)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


@dataclass(frozen=True)
class BoschSPTSWaferBoundaryMapTrace:
    reactor: BoschSPTSReducedReactorSolution
    species_names: tuple[str, ...]
    x_m: np.ndarray
    y_m: np.ndarray
    point_flux_m2_s: np.ndarray
    wafer_area_average_flux_m2_s: np.ndarray
    maximum_cylindrical_species_ledger_relative_residual: float
    maximum_cylindrical_inventory_relative_residual: float
    maximum_cylindrical_linear_system_relative_residual: float
    source_jvp_supported: bool
    provenance: dict

    def __post_init__(self):
        names = tuple(self.species_names)
        x = np.asarray(self.x_m, dtype=float).copy()
        y = np.asarray(self.y_m, dtype=float).copy()
        point = np.asarray(self.point_flux_m2_s, dtype=float).copy()
        average = np.asarray(self.wafer_area_average_flux_m2_s, dtype=float).copy()
        intervals = self.reactor.interval_duration_s.size
        if (names != _SPECIES or x.ndim != 1 or y.shape != x.shape or x.size < 4
                or point.shape != (intervals, len(names), x.size)
                or average.shape != (intervals, len(names))
                or any(np.any(~np.isfinite(value)) for value in (x, y, point, average))
                or np.any(point < 0.0) or np.any(average < 0.0)
                or len(set(zip(x, y))) != x.size
                or not 0.0 <= self.maximum_cylindrical_species_ledger_relative_residual < 1e-8
                or not 0.0 <= self.maximum_cylindrical_inventory_relative_residual < 1e-10
                or not 0.0 <= self.maximum_cylindrical_linear_system_relative_residual < 1e-10
                or not bool(self.source_jvp_supported)):
            raise ValueError("invalid cylindrical Bosch wafer boundary map")
        for value in (x, y, point, average):
            value.setflags(write=False)
        object.__setattr__(self, "species_names", names)
        object.__setattr__(self, "x_m", x)
        object.__setattr__(self, "y_m", y)
        object.__setattr__(self, "point_flux_m2_s", point)
        object.__setattr__(self, "wafer_area_average_flux_m2_s", average)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def _point_interpolation_weights(grid, x_m, y_m):
    x = np.asarray(x_m, dtype=float)
    y = np.asarray(y_m, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or np.any(~np.isfinite(x + y)):
        raise ValueError("invalid wafer sample coordinates")
    radius = np.sqrt(x * x + y * y)
    if np.any(radius > grid.geometry.radius_m):
        raise ValueError("wafer sample lies outside the cylindrical reactor")
    phi = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    radial = grid.radial_centers_m
    phi_center = grid.azimuthal_centers_rad
    dphi = 2.0 * np.pi / grid.azimuthal_cell_count
    weights = np.zeros((x.size, grid.radial_cell_count, grid.azimuthal_cell_count))
    for point, (r_value, phi_value) in enumerate(zip(radius, phi)):
        if r_value <= radial[0]:
            weights[point, 0, :] = 1.0 / grid.azimuthal_cell_count
            continue
        if r_value >= radial[-1]:
            r0 = r1 = radial.size - 1
            radial_fraction = 0.0
        else:
            r1 = int(np.searchsorted(radial, r_value))
            r0 = r1 - 1
            radial_fraction = (
                (r_value - radial[r0]) / (radial[r1] - radial[r0]))
        phi_coordinate = (phi_value - phi_center[0]) / dphi
        phi_floor = math.floor(phi_coordinate)
        p0 = phi_floor % grid.azimuthal_cell_count
        p1 = (p0 + 1) % grid.azimuthal_cell_count
        phi_fraction = phi_coordinate - phi_floor
        for radial_index, radial_weight in (
                (r0, 1.0 - radial_fraction), (r1, radial_fraction)):
            weights[point, radial_index, p0] += (
                radial_weight * (1.0 - phi_fraction))
            weights[point, radial_index, p1] += radial_weight * phi_fraction
    if np.max(np.abs(np.sum(weights, axis=(1, 2)) - 1.0)) > 1.0e-13:
        raise RuntimeError("cylindrical point interpolation lost normalization")
    return weights


def _source_shapes(grid, parameters):
    base = parameters.reduced
    shapes = []
    for species in range(len(_SPECIES)):
        cosine = parameters.source_cosine_coefficients[species]
        sine = parameters.source_sine_coefficients[species]
        annular = normalized_cylindrical_annular_skin_source(
            grid, axial_skin_depth_m=base.source_axial_skin_depth_m,
            ring_radius_m=parameters.species_source_ring_radius_m[species],
            radial_width_m=parameters.species_source_radial_width_m[species],
            cosine_coefficients=cosine, sine_coefficients=sine)
        central = normalized_cylindrical_annular_skin_source(
            grid, axial_skin_depth_m=base.source_axial_skin_depth_m,
            ring_radius_m=0.0,
            radial_width_m=(base.central_source_radial_scale_m / math.sqrt(2.0)),
            cosine_coefficients=cosine, sine_coefficients=sine)
        central_fraction = base.source_central_fraction[species]
        shapes.append(
            central_fraction * central + (1.0 - central_fraction) * annular)
    return np.stack(shapes)


def _wafer_area_weights(grid, wafer_radius_m):
    clipped_outer = np.minimum(grid.radial_edges_m[1:], wafer_radius_m)
    clipped_inner = np.minimum(grid.radial_edges_m[:-1], wafer_radius_m)
    sector_area = (
        0.5 * np.maximum(clipped_outer ** 2 - clipped_inner ** 2, 0.0)
        * np.diff(grid.azimuthal_edges_rad)[0])
    area_weights = np.repeat(
        sector_area[:, None], grid.azimuthal_cell_count, axis=1)
    area_weights /= np.sum(area_weights)
    return area_weights


def _ion_edge_focus_factors(unit_lower, grid, parameters, point_radius_m):
    """Resolve a sub-grid wafer-edge ion layer at fixed total ion current."""
    unit_lower = np.asarray(unit_lower, dtype=float)
    point_radius = np.asarray(point_radius_m, dtype=float)
    amplitude = float(parameters.ion_edge_focus_amplitude)
    if amplitude == 0.0:
        return (
            np.ones_like(point_radius),
            np.ones((grid.radial_cell_count, grid.azimuthal_cell_count)),
            1.0,
            1.0,
        )
    onset = float(parameters.ion_edge_focus_onset_radius_m)
    width = float(parameters.ion_edge_focus_width_m)
    wafer_radius = float(parameters.reduced.wafer_radius_m)

    quadrature_radius = 0.5 * wafer_radius * (_GAUSS_X + 1.0)
    quadrature_sigmoid = 1.0 / (
        1.0 + np.exp(-(quadrature_radius - onset) / width))
    sigmoid_area_mean = float(np.sum(
        _GAUSS_W * quadrature_radius * quadrature_sigmoid)
        * wafer_radius / (wafer_radius ** 2))
    continuous_area_mean = 1.0 + amplitude * 2.0 * sigmoid_area_mean

    radial_sigmoid = 1.0 / (
        1.0 + np.exp(-(grid.radial_centers_m - onset) / width))
    area_normalized = (
        1.0 + amplitude * radial_sigmoid[:, None]) / continuous_area_mean
    area_weights = _wafer_area_weights(grid, wafer_radius)
    initial_current = float(np.sum(unit_lower[2] * area_weights))
    candidate_current = float(np.sum(
        unit_lower[2] * area_normalized * area_weights))
    if initial_current <= 0.0 or candidate_current <= 0.0:
        raise RuntimeError("ion edge focus encountered zero wafer current")
    current_correction = initial_current / candidate_current
    final_current = float(np.sum(
        unit_lower[2] * area_normalized * current_correction * area_weights))
    current_residual = abs(final_current - initial_current) / initial_current
    if current_residual > 2.0e-14:
        raise RuntimeError("ion edge focus failed current conservation")
    point_sigmoid = 1.0 / (
        1.0 + np.exp(-(point_radius - onset) / width))
    point_factor = (
        (1.0 + amplitude * point_sigmoid) / continuous_area_mean
        * current_correction)
    if np.any(~np.isfinite(point_factor)) or np.any(point_factor <= 0.0):
        raise RuntimeError("ion edge focus produced an invalid point factor")
    grid_factor = area_normalized * current_correction
    return (
        point_factor,
        np.broadcast_to(
            grid_factor,
            (grid.radial_cell_count, grid.azimuthal_cell_count)).copy(),
        continuous_area_mean,
        current_correction,
    )


def _ion_edge_focus_point_factor(unit_lower, grid, parameters, point_radius_m):
    point, _grid, continuous_mean, correction = _ion_edge_focus_factors(
        unit_lower, grid, parameters, point_radius_m)
    return point, continuous_mean, correction


def _wafer_ion_transmission_point_factor(
        unit_lower, grid, parameters, point_x_m, point_y_m,
        edge_focus_grid_factor, law, c4f8_platen_vpp_rms_V):
    """Evaluate a positive Zernike map with exact grid-current normalization."""
    if not isinstance(law, BoschSPTSWaferIonTransmissionLaw):
        raise TypeError("invalid Bosch wafer ion-transmission law")
    x = np.asarray(point_x_m, dtype=float)
    y = np.asarray(point_y_m, dtype=float)
    point_radius = np.hypot(x, y)
    wafer_radius = float(parameters.reduced.wafer_radius_m)
    if np.any(point_radius > wafer_radius + 2.0e-14):
        raise ValueError("ion-transmission point lies outside the wafer")
    standardized_vpp = law.standardized_vpp(c4f8_platen_vpp_rms_V)

    radial = grid.radial_centers_m[:, None]
    angle = grid.azimuthal_centers_rad[None, :]
    on_wafer = radial <= wafer_radius
    grid_log_field = np.zeros(
        (grid.radial_cell_count, grid.azimuthal_cell_count))
    if np.any(on_wafer):
        rho = np.broadcast_to(
            radial / wafer_radius, grid_log_field.shape)[on_wafer.repeat(
                grid.azimuthal_cell_count, axis=1)]
        phi = np.broadcast_to(angle, grid_log_field.shape)[on_wafer.repeat(
            grid.azimuthal_cell_count, axis=1)]
        grid_log_field[on_wafer.repeat(
            grid.azimuthal_cell_count, axis=1)] = law.log_field(
                rho, phi, standardized_vpp)
    raw_grid_factor = np.exp(grid_log_field)
    area_weights = _wafer_area_weights(grid, wafer_radius)
    baseline_grid_current = (
        np.asarray(unit_lower, dtype=float)[2]
        * np.asarray(edge_focus_grid_factor, dtype=float))
    initial_current = float(np.sum(baseline_grid_current * area_weights))
    candidate_current = float(np.sum(
        baseline_grid_current * raw_grid_factor * area_weights))
    if initial_current <= 0.0 or candidate_current <= 0.0:
        raise RuntimeError("ion-transmission map encountered zero wafer current")
    current_correction = initial_current / candidate_current
    final_current = float(np.sum(
        baseline_grid_current * raw_grid_factor * current_correction
        * area_weights))
    current_residual = abs(final_current - initial_current) / initial_current
    if current_residual > 2.0e-14:
        raise RuntimeError("ion-transmission map failed current conservation")

    point_phi = np.arctan2(y, x)
    point_log_field = law.log_field(
        np.minimum(point_radius / wafer_radius, 1.0),
        point_phi,
        standardized_vpp,
    )
    point_factor = np.exp(point_log_field) * current_correction
    if np.any(~np.isfinite(point_factor)) or np.any(point_factor <= 0.0):
        raise RuntimeError("ion-transmission map produced an invalid factor")
    return point_factor, {
        "standardized_c4f8_platen_vpp_rms": standardized_vpp,
        "exact_current_correction": current_correction,
        "relative_total_ion_current_residual": current_residual,
        "minimum_point_factor": float(np.min(point_factor)),
        "maximum_point_factor": float(np.max(point_factor)),
        "minimum_grid_factor": float(np.min(
            raw_grid_factor * current_correction)),
        "maximum_grid_factor": float(np.max(
            raw_grid_factor * current_correction)),
    }


class DeterministicBoschSPTSCylindricalReactorToWafer:
    """0-D measured waveform plus positive cylindrical 3-D inventory lift."""

    def __init__(self, parameters: BoschSPTSCylindricalParameters):
        if not isinstance(parameters, BoschSPTSCylindricalParameters):
            raise TypeError("parameters must be BoschSPTSCylindricalParameters")
        self.parameters = parameters
        base = parameters.reduced
        grid = CylindricalFiniteVolumeGrid.uniform(
            CylindricalReactor(base.reactor_radius_m, base.reactor_length_m),
            radial_cell_count=base.radial_cell_count,
            azimuthal_cell_count=parameters.azimuthal_cell_count,
            axial_cell_count=base.axial_cell_count)
        source_shapes = _source_shapes(grid, parameters)
        wall = np.column_stack((
            base.lower_wall_velocity_m_s,
            base.upper_wall_velocity_m_s,
            base.side_wall_velocity_m_s))
        # Chamber conditioning changes neutral recombination on non-wafer
        # walls.  Lower-endcap collection is deliberately unchanged so this
        # equipment closure cannot mutate the frozen wafer surface law.
        wall[:2, 1:] *= base.neutral_wall_loss_multiplier
        self._lift = DeterministicCylindricalIndependentInventoryLift(
            grid=grid, species_names=_SPECIES,
            diffusion_coefficient_m2_s=base.diffusion_coefficient_m2_s,
            wall_velocity_m_s=wall, source_shape=source_shapes,
            source=(
                "SPTS species-resolved positive cylindrical source moments; "
                "source-2 measured off"))
        self._grid = grid
        self._maximum_ledger = self._lift.maximum_unit_ledger_relative_residual
        self._maximum_linear = (
            self._lift.maximum_unit_linear_system_relative_residual)

    def source_response(
            self, *, x_m, y_m, source_cosine_coefficients=None,
            source_sine_coefficients=None, source_ring_radius_m=None,
            source_radial_width_m=None, source_central_fraction=None,
            species_source_ring_radius_m=None,
            species_source_radial_width_m=None,
            ion_edge_focus_amplitude=None,
            ion_edge_focus_onset_radius_m=None,
            ion_edge_focus_width_m=None,
            ion_transmission_law=None,
            c4f8_platen_vpp_rms_V=None):
        if ((source_cosine_coefficients is None)
                != (source_sine_coefficients is None)):
            raise ValueError("cosine and sine source coefficients must be paired")
        if ((ion_transmission_law is None)
                != (c4f8_platen_vpp_rms_V is None)):
            raise ValueError(
                "ion-transmission law and C4F8 Vpp RMS must be paired")
        if (source_ring_radius_m is not None
                and species_source_ring_radius_m is not None):
            raise ValueError("shared and species ring radii are mutually exclusive")
        if (source_radial_width_m is not None
                and species_source_radial_width_m is not None):
            raise ValueError("shared and species radial widths are mutually exclusive")
        source_geometry_changed = any(value is not None for value in (
            source_ring_radius_m, source_radial_width_m,
            source_central_fraction, species_source_ring_radius_m,
            species_source_radial_width_m))
        if source_cosine_coefficients is None and not source_geometry_changed:
            parameters = replace(
                self.parameters,
                ion_edge_focus_amplitude=(
                    self.parameters.ion_edge_focus_amplitude
                    if ion_edge_focus_amplitude is None
                    else ion_edge_focus_amplitude),
                ion_edge_focus_onset_radius_m=(
                    self.parameters.ion_edge_focus_onset_radius_m
                    if ion_edge_focus_onset_radius_m is None
                    else ion_edge_focus_onset_radius_m),
                ion_edge_focus_width_m=(
                    self.parameters.ion_edge_focus_width_m
                    if ion_edge_focus_width_m is None
                    else ion_edge_focus_width_m))
            unit_lower = self._lift.unit_lower_flux_per_density_m_s
            maximum_ledger = self._maximum_ledger
            maximum_linear = self._maximum_linear
            reused_factorization = True
        else:
            base = self.parameters.reduced
            reduced = replace(
                base,
                source_ring_radius_m=(
                    base.source_ring_radius_m if source_ring_radius_m is None
                    else source_ring_radius_m),
                source_radial_width_m=(
                    base.source_radial_width_m if source_radial_width_m is None
                    else source_radial_width_m),
                source_central_fraction=(
                    base.source_central_fraction if source_central_fraction is None
                    else source_central_fraction))
            parameters = BoschSPTSCylindricalParameters(
                reduced=reduced,
                azimuthal_cell_count=self.parameters.azimuthal_cell_count,
                species_source_ring_radius_m=(
                    species_source_ring_radius_m
                    if species_source_ring_radius_m is not None
                    else (None if source_ring_radius_m is not None
                          else self.parameters.species_source_ring_radius_m)),
                species_source_radial_width_m=(
                    species_source_radial_width_m
                    if species_source_radial_width_m is not None
                    else (None if source_radial_width_m is not None
                          else self.parameters.species_source_radial_width_m)),
                source_cosine_coefficients=(
                    self.parameters.source_cosine_coefficients
                    if source_cosine_coefficients is None
                    else source_cosine_coefficients),
                source_sine_coefficients=(
                    self.parameters.source_sine_coefficients
                    if source_sine_coefficients is None
                    else source_sine_coefficients),
                ion_edge_focus_amplitude=(
                    self.parameters.ion_edge_focus_amplitude
                    if ion_edge_focus_amplitude is None
                    else ion_edge_focus_amplitude),
                ion_edge_focus_onset_radius_m=(
                    self.parameters.ion_edge_focus_onset_radius_m
                    if ion_edge_focus_onset_radius_m is None
                    else ion_edge_focus_onset_radius_m),
                ion_edge_focus_width_m=(
                    self.parameters.ion_edge_focus_width_m
                    if ion_edge_focus_width_m is None
                    else ion_edge_focus_width_m))
            unit_lower, maximum_ledger, maximum_linear = (
                self._lift.source_shape_to_unit_lower_flux(
                    _source_shapes(self._grid, parameters)))
            reused_factorization = True
        weights = _point_interpolation_weights(self._grid, x_m, y_m)
        unit_point_flux = np.einsum("srp,qrp->sq", unit_lower, weights)
        point_radius = np.sqrt(
            np.asarray(x_m, dtype=float) ** 2 + np.asarray(y_m, dtype=float) ** 2)
        (point_factor, edge_grid_factor, continuous_area_mean,
         current_correction) = _ion_edge_focus_factors(
            unit_lower, self._grid, parameters, point_radius)
        unit_point_flux[2] *= point_factor
        if ion_transmission_law is None:
            transmission_provenance = {
                "enabled": False,
                "relative_total_ion_current_residual": 0.0,
            }
        else:
            transmission_factor, transmission_provenance = (
                _wafer_ion_transmission_point_factor(
                    unit_lower,
                    self._grid,
                    parameters,
                    np.asarray(x_m, dtype=float),
                    np.asarray(y_m, dtype=float),
                    edge_grid_factor,
                    ion_transmission_law,
                    c4f8_platen_vpp_rms_V,
                ))
            unit_point_flux[2] *= transmission_factor
            transmission_provenance = {
                "enabled": True,
                "law": ion_transmission_law.manifest(),
                "c4f8_platen_vpp_rms_V": float(c4f8_platen_vpp_rms_V),
                **transmission_provenance,
            }
        area_weights = _wafer_area_weights(
            self._grid, parameters.reduced.wafer_radius_m)
        unit_average_flux = np.einsum(
            "srp,rp->s", unit_lower, area_weights)
        return BoschSPTSCylindricalSourceResponse(
            species_names=_SPECIES, x_m=np.asarray(x_m), y_m=np.asarray(y_m),
            unit_point_flux_per_density_m_s=unit_point_flux,
            unit_wafer_average_flux_per_density_m_s=unit_average_flux,
            maximum_species_ledger_relative_residual=maximum_ledger,
            maximum_linear_system_relative_residual=maximum_linear,
            provenance={
                "model": "spts-bosch-cylindrical-source-response-v1",
                "parameters": parameters.manifest(),
                "transport_factorization_reused": reused_factorization,
                "ion_edge_focus_continuous_area_mean": continuous_area_mean,
                "ion_edge_focus_exact_current_correction": current_correction,
                "wafer_ion_transmission": transmission_provenance,
                "total_wafer_ion_current_conserved": True,
                "target_depth_used": False,
            })

    def solve(
            self, trace, *, x_m, y_m, initial_density_m3=(0.0, 0.0, 0.0),
            source_response=None):
        base = self.parameters.reduced
        reactor = solve_bosch_spts_reduced_reactor(
            trace, base, initial_density_m3=initial_density_m3)
        if source_response is None:
            source_response = self.source_response(x_m=x_m, y_m=y_m)
        if (not isinstance(source_response, BoschSPTSCylindricalSourceResponse)
                or not np.array_equal(source_response.x_m, np.asarray(x_m))
                or not np.array_equal(source_response.y_m, np.asarray(y_m))):
            raise ValueError("source response does not match requested wafer points")
        point_flux = (
            reactor.volume_average_density_m3[:, :, None]
            * source_response.unit_point_flux_per_density_m_s[None])
        average_flux = (
            reactor.volume_average_density_m3
            * source_response.unit_wafer_average_flux_per_density_m_s[None])
        return BoschSPTSWaferBoundaryMapTrace(
            reactor=reactor, species_names=_SPECIES,
            x_m=np.asarray(x_m), y_m=np.asarray(y_m),
            point_flux_m2_s=point_flux,
            wafer_area_average_flux_m2_s=average_flux,
            maximum_cylindrical_species_ledger_relative_residual=(
                source_response.maximum_species_ledger_relative_residual),
            maximum_cylindrical_inventory_relative_residual=0.0,
            maximum_cylindrical_linear_system_relative_residual=(
                source_response.maximum_linear_system_relative_residual),
            source_jvp_supported=True,
            provenance={
                "model": "spts-bosch-measured-waveform-cylindrical-wafer-v1",
                "parameters": source_response.provenance["parameters"],
                "point_interpolation": "bilinear r/periodic-phi; axis-averaged at r=0",
                "target_depth_used": False,
            })

    def density_to_point_flux_jvp(
            self, density_tangent_m3, *, x_m, y_m, source_response=None):
        tangent = np.asarray(density_tangent_m3, dtype=float)
        if (tangent.ndim != 2 or tangent.shape[1] != len(_SPECIES)
                or np.any(~np.isfinite(tangent))):
            raise ValueError("density tangent must have shape (time, 3)")
        if source_response is None:
            source_response = self.source_response(x_m=x_m, y_m=y_m)
        if (not isinstance(source_response, BoschSPTSCylindricalSourceResponse)
                or not np.array_equal(source_response.x_m, np.asarray(x_m))
                or not np.array_equal(source_response.y_m, np.asarray(y_m))):
            raise ValueError("source response does not match requested wafer points")
        return (
            tangent[:, :, None]
            * source_response.unit_point_flux_per_density_m_s[None])
