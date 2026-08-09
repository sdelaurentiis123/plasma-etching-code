"""Deterministic line-resolved Cl I source-to-wafer radiation boundary.

The direct-coronal source supplies primary photons created by electron-impact
excitation.  A conservative partial-frequency-redistribution solve then tracks
each photon through ground-state absorption, same-line re-emission, branching
to other lines, and escape onto the wafer disk or the other chamber surfaces.
The operator is deterministic; the known incompleteness of the coronal source
and spatially uniform reactor state remains explicit in the returned receipt.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import atomic_mass, physical_constants

from .chlorine_vuv_spectrum import Adf04Level, ChlorineVuvLine
from .geometry import CylindricalReactor
from .vuv_radiation import ResonanceLineData
from .zonal_radiation import (
    AxisymmetricRadiationZoneField,
    ZonalPartialRedistributionResult,
    deterministic_zonal_partial_redistribution,
)


_BOLTZMANN_EV_K = physical_constants["Boltzmann constant in eV/K"][0]
_CHLORINE_ATOM_MASS_KG = 35.45 * atomic_mass


@dataclass(frozen=True)
class ChlorineLineWaferBoundaryResult:
    """One line's primary source, resonant transfer, and wafer flux receipt."""

    wavelength_nm: float
    lower_observed_index: int
    upper_observed_index: int
    lower_ground_population_fraction: float
    lower_state_absorber_density_m3: float
    primary_line_emissivity_m3_s: float
    primary_line_emission_rate_s: float
    alternate_branch_loss_frequency_s_inv: float
    wafer_radius_m: float
    wafer_photon_flux_m2_s: float
    radiation: ZonalPartialRedistributionResult
    prediction_supported: bool
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        values = np.asarray((
            self.wavelength_nm,
            self.lower_ground_population_fraction,
            self.lower_state_absorber_density_m3,
            self.primary_line_emissivity_m3_s,
            self.primary_line_emission_rate_s,
            self.alternate_branch_loss_frequency_s_inv,
            self.wafer_radius_m,
            self.wafer_photon_flux_m2_s,
        ), dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.wavelength_nm <= 0.0
            or not 0.0 < self.lower_ground_population_fraction <= 1.0
            or np.any(values[2:] < 0.0)
            or self.wafer_radius_m <= 0.0
            or not isinstance(self.radiation, ZonalPartialRedistributionResult)
            or bool(self.prediction_supported)
            or not self.known_limitations
        ):
            raise ValueError("invalid chlorine line wafer-boundary result")
        object.__setattr__(
            self, "known_limitations", tuple(self.known_limitations))


def chlorine_ground_population_fraction(
    observed_levels: tuple[Adf04Level, ...] | list[Adf04Level],
    lower_observed_index: int,
    *,
    gas_temperature_K: float,
) -> float:
    """Return the thermal population of one Cl I ground fine-structure level."""
    temperature = float(gas_temperature_K)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("gas temperature must be finite and positive")
    by_index = {level.index: level for level in observed_levels}
    if set((1, 2)) - set(by_index):
        raise ValueError("observed Cl I ground fine-structure levels are missing")
    if int(lower_observed_index) not in (1, 2):
        raise ValueError("resonance boundary requires a ground-terminating line")
    ground = (by_index[1], by_index[2])
    energy_origin = min(level.energy_cm_inv for level in ground)
    weights = np.asarray([
        level.statistical_weight * math.exp(
            -(level.energy_cm_inv - energy_origin)
            * 1.2398419843320026e-4
            / (_BOLTZMANN_EV_K * temperature)
        )
        for level in ground
    ])
    weights /= np.sum(weights)
    return float(weights[(1, 2).index(int(lower_observed_index))])


def deterministic_uniform_chlorine_line_wafer_boundary(
    geometry: CylindricalReactor,
    line: ChlorineVuvLine,
    observed_levels: tuple[Adf04Level, ...] | list[Adf04Level],
    *,
    wafer_radius_m: float,
    electron_density_m3: float,
    chlorine_atom_density_m3: float,
    gas_temperature_K: float,
    velocity_changing_collision_frequency_s_inv: float = 0.0,
    nonradiative_quenching_frequency_s_inv: float = 0.0,
    surface_quadrature_order: int = 8,
    direction_quadrature_order: int = 8,
    frequency_quadrature_order: int = 24,
    coherent_grid_points_per_lorentz_hwhm: float = 8.0,
) -> ChlorineLineWaferBoundaryResult:
    """Propagate one direct-coronal Cl I line to a finite wafer disk.

    The source rate coefficient already includes the primary radiative branch.
    Following resonant absorption, competing branches from the same upper state
    are carried as destruction of this line rather than being allowed to
    re-emit at the wrong wavelength.
    """
    if not isinstance(geometry, CylindricalReactor):
        raise TypeError("a cylindrical reactor geometry is required")
    if not isinstance(line, ChlorineVuvLine):
        raise TypeError("a chlorine VUV line is required")
    electron_density = float(electron_density_m3)
    chlorine_density = float(chlorine_atom_density_m3)
    temperature = float(gas_temperature_K)
    wafer_radius = float(wafer_radius_m)
    velocity_frequency = float(
        velocity_changing_collision_frequency_s_inv)
    nonradiative_quench = float(nonradiative_quenching_frequency_s_inv)
    if (
        not math.isfinite(electron_density)
        or electron_density < 0.0
        or not math.isfinite(chlorine_density)
        or chlorine_density < 0.0
        or not math.isfinite(temperature)
        or temperature <= 0.0
        or not math.isfinite(wafer_radius)
        or not 0.0 < wafer_radius <= geometry.radius_m
        or not math.isfinite(velocity_frequency)
        or velocity_frequency < 0.0
        or not math.isfinite(nonradiative_quench)
        or nonradiative_quench < 0.0
    ):
        raise ValueError("invalid chlorine line reactor boundary")
    levels = tuple(observed_levels)
    by_index = {level.index: level for level in levels}
    try:
        lower = by_index[line.lower_observed_index]
        upper = by_index[line.upper_observed_index]
    except KeyError as error:
        raise ValueError("line level is absent from observed Cl I levels") from error
    lower_fraction = chlorine_ground_population_fraction(
        levels,
        line.lower_observed_index,
        gas_temperature_K=temperature,
    )
    absorber_density = chlorine_density * lower_fraction
    primary_emissivity = (
        electron_density
        * chlorine_density
        * line.photon_rate_coefficient_cm3_s
        * 1.0e-6
    )
    # Emitter density sets the initial spatial source weights in the radiation
    # solver.  A uniform nonzero proxy is sufficient when the primary source is
    # zero because amplitude is applied separately below.
    emitter_weight = max(primary_emissivity, 1.0)
    field = AxisymmetricRadiationZoneField(
        radial_edges_m=np.asarray([0.0, geometry.radius_m]),
        axial_edges_m=np.asarray([0.0, geometry.length_m]),
        cell_zone_index=np.asarray([[0]]),
        gas_temperature_K=np.asarray([temperature]),
        absorber_density_m3=np.asarray([absorber_density]),
        emitter_density_m3=np.asarray([emitter_weight]),
        source=(
            "uniform direct-coronal Cl I primary source and ground-state "
            "absorber moment"),
    )
    resonance_line = ResonanceLineData(
        wavelength_nm=line.wavelength_nm,
        transition_probability_s_inv=line.transition_probability_s_inv,
        lower_statistical_weight=lower.statistical_weight,
        upper_statistical_weight=upper.statistical_weight,
        absorber_mass_kg=_CHLORINE_ATOM_MASS_KG,
        source=(
            "OPEN-ADAS transition probability plus observed NIST Cl I "
            "fine-structure levels"),
    )
    alternate_branch = max(
        0.0,
        line.upper_total_radiative_probability_s_inv
        - line.transition_probability_s_inv,
    )
    radiation = deterministic_zonal_partial_redistribution(
        field,
        resonance_line,
        wafer_radius_m=wafer_radius,
        velocity_changing_collision_frequency_s_inv=velocity_frequency,
        quenching_collision_frequency_s_inv=(
            alternate_branch + nonradiative_quench),
        surface_quadrature_order=surface_quadrature_order,
        direction_quadrature_order=direction_quadrature_order,
        frequency_quadrature_order=frequency_quadrature_order,
        coherent_grid_points_per_lorentz_hwhm=(
            coherent_grid_points_per_lorentz_hwhm),
    )
    primary_rate = primary_emissivity * geometry.volume_m3
    wafer_area = math.pi * wafer_radius ** 2
    wafer_flux = radiation.partial_redistribution_wafer_flux_m2_s(
        total_line_emission_rate_s=primary_rate,
        wafer_area_m2=wafer_area,
    )
    return ChlorineLineWaferBoundaryResult(
        wavelength_nm=line.wavelength_nm,
        lower_observed_index=line.lower_observed_index,
        upper_observed_index=line.upper_observed_index,
        lower_ground_population_fraction=lower_fraction,
        lower_state_absorber_density_m3=absorber_density,
        primary_line_emissivity_m3_s=primary_emissivity,
        primary_line_emission_rate_s=primary_rate,
        alternate_branch_loss_frequency_s_inv=alternate_branch,
        wafer_radius_m=wafer_radius,
        wafer_photon_flux_m2_s=wafer_flux,
        radiation=radiation,
        prediction_supported=False,
        known_limitations=(
            "direct excitation only; excited-state cascades are absent",
            "uniform global electron and chlorine densities replace spatial fields",
            "velocity-changing and nonradiative collision rates require external inputs",
            "OPEN-ADAS collision strengths are a distorted-wave calculation, not a target-tool spectrum measurement",
        ),
    )
