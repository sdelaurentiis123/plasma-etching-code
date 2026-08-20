"""Measured-waveform SPTS Bosch reactor with cylindrical 3-D wafer transfer."""
from __future__ import annotations

from dataclasses import dataclass, field
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


@dataclass(frozen=True)
class BoschSPTSCylindricalParameters:
    reduced: BoschSPTSReducedParameters = field(
        default_factory=BoschSPTSReducedParameters)
    azimuthal_cell_count: int = 16
    source_cosine_coefficients: tuple[tuple[float, ...], ...] = ((), (), ())
    source_sine_coefficients: tuple[tuple[float, ...], ...] = ((), (), ())

    def __post_init__(self):
        if not isinstance(self.reduced, BoschSPTSReducedParameters):
            raise TypeError("reduced parameters must be BoschSPTSReducedParameters")
        if (int(self.azimuthal_cell_count) != self.azimuthal_cell_count
                or self.azimuthal_cell_count < 8):
            raise ValueError("cylindrical Bosch transfer requires at least 8 phi cells")
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

    def manifest(self):
        return {
            "schema": "petch-spts-bosch-cylindrical-parameters-v1",
            "reduced": self.reduced.manifest(),
            "azimuthal_cell_count": self.azimuthal_cell_count,
            "source_cosine_coefficients": [
                list(row) for row in self.source_cosine_coefficients],
            "source_sine_coefficients": [
                list(row) for row in self.source_sine_coefficients],
            "azimuthal_source_form": (
                "positive exponential Fourier modulation of each species source"),
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
            ring_radius_m=base.source_ring_radius_m,
            radial_width_m=base.source_radial_width_m,
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
            source_sine_coefficients=None):
        if ((source_cosine_coefficients is None)
                != (source_sine_coefficients is None)):
            raise ValueError("cosine and sine source coefficients must be paired")
        if source_cosine_coefficients is None:
            parameters = self.parameters
            unit_lower = self._lift.unit_lower_flux_per_density_m_s
            maximum_ledger = self._maximum_ledger
            maximum_linear = self._maximum_linear
            reused_factorization = True
        else:
            parameters = BoschSPTSCylindricalParameters(
                reduced=self.parameters.reduced,
                azimuthal_cell_count=self.parameters.azimuthal_cell_count,
                source_cosine_coefficients=source_cosine_coefficients,
                source_sine_coefficients=source_sine_coefficients)
            unit_lower, maximum_ledger, maximum_linear = (
                self._lift.source_shape_to_unit_lower_flux(
                    _source_shapes(self._grid, parameters)))
            reused_factorization = True
        weights = _point_interpolation_weights(self._grid, x_m, y_m)
        unit_point_flux = np.einsum("srp,qrp->sq", unit_lower, weights)
        area_weights = _wafer_area_weights(
            self._grid, self.parameters.reduced.wafer_radius_m)
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
