"""Fused deterministic batch evaluator for the common Bosch surface laws.

The canonical mechanism objects intentionally produce rich per-step ledgers
and diagnostics.  Calibration evaluates thousands of measured 5 Hz intervals
over many wafers, so repeatedly constructing every diagnostic dataclass is
unnecessarily expensive.  This module fuses the exact algebraic recurrence of
the unchanged Belen and La Magna mechanisms across a padded wafer batch.

Parity tests compare every returned depth and film inventory with
``predict_bosch_wafer_depth``.  This is an execution optimization, not a new
surface model or a depth surrogate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType

import numpy as np

from .bosch_silicon import BoschSiliconFluorocarbonMechanism
from .bosch_wafer_depth import BoschWaferDepthPrediction
from .fluorocarbon_lamagna import LaMagnaGarozzoFluorocarbonMechanism
from .reactor_global.bosch_spts_reduced import BoschSPTSWaferBoundaryTrace
from .reactor_global.bosch_spts_cylindrical import (
    BoschSPTSWaferBoundaryMapTrace,
)


@dataclass(frozen=True)
class BoschWaferPointDepthPrediction:
    """Surface evolution at one identified set of wafer measurement points."""

    experiment_key: str
    x_m: np.ndarray
    y_m: np.ndarray
    silicon_depth_m: np.ndarray
    oxide_loss_m: np.ndarray
    remaining_film_units_m2: np.ndarray
    measurement_point_mean_silicon_depth_m: float
    measurement_point_mean_oxide_loss_m: float
    silicon_to_oxide_selectivity: float
    maximum_material_ledger_relative_residual: float
    provenance: dict

    def __post_init__(self):
        arrays = {
            name: np.asarray(getattr(self, name), dtype=float).copy()
            for name in (
                "x_m", "y_m", "silicon_depth_m", "oxide_loss_m",
                "remaining_film_units_m2",
            )
        }
        size = arrays["x_m"].size
        if (
            not str(self.experiment_key).strip()
            or size < 4
            or any(value.shape != (size,) for value in arrays.values())
            or any(np.any(~np.isfinite(value)) for value in arrays.values())
            or np.any(arrays["silicon_depth_m"] < 0.0)
            or np.any(arrays["oxide_loss_m"] < 0.0)
            or np.any(arrays["remaining_film_units_m2"] < 0.0)
            or len(set(zip(arrays["x_m"], arrays["y_m"]))) != size
        ):
            raise ValueError("invalid Bosch wafer point-depth prediction")
        mean_si = float(np.mean(arrays["silicon_depth_m"]))
        mean_oxide = float(np.mean(arrays["oxide_loss_m"]))
        if (
            not math.isclose(
                self.measurement_point_mean_silicon_depth_m, mean_si,
                rel_tol=2.0e-14, abs_tol=1.0e-20)
            or not math.isclose(
                self.measurement_point_mean_oxide_loss_m, mean_oxide,
                rel_tol=2.0e-14, abs_tol=1.0e-20)
            or not math.isfinite(self.silicon_to_oxide_selectivity)
            or self.silicon_to_oxide_selectivity < 0.0
            or not 0.0 <= self.maximum_material_ledger_relative_residual < 1.0e-10
        ):
            raise ValueError("invalid Bosch wafer point-depth summary")
        for name, value in arrays.items():
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def _yield_rate(flux, energy_eV, law):
    yield_value = law.evaluate(energy_eV, np.ones_like(energy_eV))
    return flux * yield_value


def _film_solution(F, polymer, ion_sp, ion_ie, ion_polymer, parameters, beta_e,
                   beta_p):
    shape = F.shape
    zero = parameters.coverage_zero_tolerance_m2_s
    pe_denominator = (
        F * parameters.polymer_etchant_sticking_probability + ion_polymer)
    pe = np.divide(
        F * parameters.polymer_etchant_sticking_probability, pe_denominator,
        out=np.zeros(shape), where=pe_denominator > zero)
    p_denominator = ion_polymer * pe + parameters.polymer_loss_rate_m2_s
    raw_p = np.divide(
        polymer * beta_p, p_denominator,
        out=np.zeros(shape), where=p_denominator > zero)
    forced_saturation = (
        (polymer > zero) & ((pe <= zero) | (ion_polymer <= zero)))
    saturated = forced_saturation | (raw_p >= 1.0)
    p = np.where(saturated, 1.0, raw_p)

    boltzmann_eV_K = 8.617333262145e-5
    chemical_reference = (
        parameters.chemical_rate_coefficient
        * parameters.reference_etchant_flux_m2_s
        * math.exp(-parameters.chemical_activation_energy_eV
                   / (boltzmann_eV_K * parameters.temperature_K)))
    e_supply = F * beta_e
    e_denominator = (
        parameters.ion_enhanced_coverage_loss_factor * ion_ie
        + parameters.chemical_coverage_loss_factor * chemical_reference
        + e_supply)
    e = np.divide(
        e_supply * (1.0 - p), e_denominator,
        out=np.zeros(shape), where=e_denominator > zero)
    substrate_removal = (
        chemical_reference * e + ion_ie * e + ion_sp * (1.0 - e))
    polymer_net = (
        polymer * parameters.polymer_polymer_sticking_probability
        - ion_polymer * pe)
    return saturated, substrate_removal, polymer_net


def _radial_area(boundary):
    parameters = boundary.reactor.provenance["parameters"]
    wafer_radius = float(parameters["wafer_radius_m"])
    reactor_radius = float(parameters["reactor_radius_m"])
    edges = np.linspace(
        0.0, reactor_radius, boundary.radial_centers_m.size + 1)
    clipped_outer = np.minimum(edges[1:], wafer_radius)
    clipped_inner = np.minimum(edges[:-1], wafer_radius)
    return np.pi * np.maximum(clipped_outer ** 2 - clipped_inner ** 2, 0.0)


def _advance_fused_surface_recurrence(
        flux_m2_s, ion_energy_eV, interval_duration_s,
        silicon_mechanism, oxide_mechanism):
    """Advance the exact common surface recurrence over arbitrary points.

    The final axis is an unstructured list of independent wafer locations.
    No spatial averaging or geometry assumption enters this recurrence.
    """
    flux = np.asarray(flux_m2_s, dtype=float)
    energy = np.asarray(ion_energy_eV, dtype=float)
    duration = np.asarray(interval_duration_s, dtype=float)
    if (
        flux.ndim != 4
        or flux.shape[2] != 3
        or energy.shape != (flux.shape[0], flux.shape[1], 1)
        or duration.shape != energy.shape
        or any(np.any(~np.isfinite(value)) for value in (flux, energy, duration))
        or np.any(flux < 0.0)
        or np.any(energy < 0.0)
        or np.any(duration < 0.0)
    ):
        raise ValueError("invalid fused Bosch surface histories")
    if not isinstance(silicon_mechanism, BoschSiliconFluorocarbonMechanism):
        raise TypeError("silicon mechanism must be the Bosch composite")
    if not isinstance(oxide_mechanism, LaMagnaGarozzoFluorocarbonMechanism):
        raise TypeError("oxide mechanism must be La Magna/Garozzo")

    silicon = silicon_mechanism.silicon_mechanism.parameters
    film = silicon_mechanism.film_mechanism.parameters
    oxide = oxide_mechanism.parameters
    if film is not oxide:
        raise ValueError(
            "fast Bosch batch requires one shared film/oxide parameter object")

    batch = flux.shape[1]
    point_count = flux.shape[3]
    removed_si = np.zeros((batch, point_count))
    removed_oxide = np.zeros_like(removed_si)
    film_inventory = np.zeros_like(removed_si)

    for interval in range(flux.shape[0]):
        dt = duration[interval]
        F = flux[interval, :, 0]
        polymer = flux[interval, :, 1]
        ions = flux[interval, :, 2]
        E = energy[interval]

        physical_si = _yield_rate(ions, E, silicon.physical_sputter_yield)
        enhanced_si = _yield_rate(ions, E, silicon.ion_enhanced_yield)
        fluorine_supply = silicon.fluorine_sticking_probability * F
        fluorine_loss = (
            silicon.spontaneous_fluorine_removal_rate_m2_s
            + silicon.ion_enhanced_fluorine_release_per_si * enhanced_si)
        fluorine_ratio = fluorine_supply / fluorine_loss
        fluorine_coverage = fluorine_ratio / (1.0 + fluorine_ratio)
        bare_si_rate = (
            silicon.spontaneous_fluorine_removal_rate_m2_s
            * fluorine_coverage / silicon.fluorine_atoms_per_removed_si
            + physical_si + enhanced_si * fluorine_coverage)

        ion_sp = _yield_rate(ions, E, film.physical_sputter_yield)
        ion_ie = _yield_rate(ions, E, film.ion_enhanced_yield)
        ion_polymer = _yield_rate(ions, E, film.polymer_removal_yield)
        substrate_saturated, substrate_removal, substrate_polymer_net = (
            _film_solution(
                F, polymer, ion_sp, ion_ie, ion_polymer, film,
                film.substrate_etchant_sticking_probability,
                film.substrate_polymer_sticking_probability))
        film_saturated, _film_removal, film_polymer_net = _film_solution(
            F, polymer, ion_sp, ion_ie, ion_polymer, film,
            film.polymer_etchant_sticking_probability,
            film.polymer_polymer_sticking_probability)

        film_present = film_inventory > 0.0
        selected_saturated = np.where(
            film_present, film_saturated, substrate_saturated)
        selected_polymer_net = np.where(
            film_present, film_polymer_net, substrate_polymer_net)
        growth_rate = np.where(
            selected_saturated, np.maximum(selected_polymer_net, 0.0), 0.0)
        film_removal_rate = np.where(
            film_present, np.maximum(-selected_polymer_net, 0.0), 0.0)
        deposited = growth_rate * dt
        removed_film = np.minimum(film_inventory, film_removal_rate * dt)
        depletion_time = np.divide(
            film_inventory, film_removal_rate,
            out=np.full_like(film_inventory, np.inf),
            where=film_removal_rate > 0.0)
        substrate_time = np.where(
            ~film_present & ~substrate_saturated, dt, 0.0)
        substrate_time = np.where(
            film_present & (depletion_time < dt) & ~substrate_saturated,
            dt - depletion_time, substrate_time)
        film_inventory = film_inventory + deposited - removed_film
        removed_si += bare_si_rate * substrate_time
        removed_oxide += substrate_removal * substrate_time

    return removed_si, removed_oxide, film_inventory, silicon, oxide


def predict_bosch_wafer_depth_batch_fast(
        boundaries,
        silicon_mechanism: BoschSiliconFluorocarbonMechanism,
        oxide_mechanism: LaMagnaGarozzoFluorocarbonMechanism,
        ) -> tuple[BoschWaferDepthPrediction, ...]:
    """Advance a variable-duration wafer batch through the exact surface law."""
    boundaries = tuple(boundaries)
    if not boundaries or any(
            not isinstance(item, BoschSPTSWaferBoundaryTrace)
            for item in boundaries):
        raise TypeError("a nonempty Bosch wafer-boundary batch is required")
    if not isinstance(silicon_mechanism, BoschSiliconFluorocarbonMechanism):
        raise TypeError("silicon mechanism must be the Bosch composite")
    if not isinstance(oxide_mechanism, LaMagnaGarozzoFluorocarbonMechanism):
        raise TypeError("oxide mechanism must be La Magna/Garozzo")
    names = ("F", "C4F8_film_precursor", "positive_ion")
    radial_count = boundaries[0].radial_centers_m.size
    if any(
            item.species_names != names
            or item.radial_centers_m.size != radial_count
            or not np.array_equal(
                item.radial_centers_m, boundaries[0].radial_centers_m)
            for item in boundaries):
        raise ValueError("Bosch batch boundaries do not share one radial grid")

    batch = len(boundaries)
    maximum_intervals = max(
        item.reactor.interval_duration_s.size for item in boundaries)
    flux = np.zeros((maximum_intervals, batch, len(names), radial_count))
    energy = np.zeros((maximum_intervals, batch, 1))
    duration = np.zeros((maximum_intervals, batch, 1))
    for wafer, boundary in enumerate(boundaries):
        count = boundary.reactor.interval_duration_s.size
        flux[:count, wafer] = boundary.radial_flux_m2_s
        energy[:count, wafer, 0] = boundary.reactor.ion_energy_eV
        duration[:count, wafer, 0] = boundary.reactor.interval_duration_s

    removed_si, removed_oxide, film_inventory, silicon, oxide = (
        _advance_fused_surface_recurrence(
            flux, energy, duration, silicon_mechanism, oxide_mechanism))

    predictions = []
    for wafer, boundary in enumerate(boundaries):
        area = _radial_area(boundary)
        area_sum = float(np.sum(area))
        silicon_depth = removed_si[wafer] / silicon.bulk_si_atom_density_m3
        oxide_loss = removed_oxide[wafer] / oxide.bulk_formula_density_m3
        mean_si = float(np.dot(silicon_depth, area) / area_sum)
        mean_oxide = float(np.dot(oxide_loss, area) / area_sum)
        predictions.append(BoschWaferDepthPrediction(
            experiment_key=boundary.reactor.trace_key,
            radial_centers_m=boundary.radial_centers_m,
            radial_area_m2=area,
            silicon_depth_m=silicon_depth,
            oxide_loss_m=oxide_loss,
            remaining_film_units_m2=film_inventory[wafer],
            wafer_mean_silicon_depth_m=mean_si,
            wafer_mean_oxide_loss_m=mean_oxide,
            silicon_to_oxide_selectivity=(
                mean_si / mean_oxide if mean_oxide > 0.0 else 0.0),
            maximum_material_ledger_relative_residual=0.0,
            provenance={
                "model": "petch-spts-bosch-wafer-depth-fused-batch-v1",
                "surface_recurrence": (
                    "exact fused algebraic parity path for unchanged Belen and "
                    "La Magna mechanisms"),
                "reactor_boundary_model": boundary.reactor.provenance["model"],
                "silicon_surface_model": silicon_mechanism.provenance["model"],
                "oxide_surface_model": oxide_mechanism.provenance["model"],
                "target_depth_used": False,
                "surface_parameters_shared_across_every_wafer": True,
            }))
    return tuple(predictions)


def predict_bosch_wafer_point_depth_batch_fast(
        boundaries,
        silicon_mechanism: BoschSiliconFluorocarbonMechanism,
        oxide_mechanism: LaMagnaGarozzoFluorocarbonMechanism,
        ) -> tuple[BoschWaferPointDepthPrediction, ...]:
    """Advance cylindrical point-flux maps through the exact surface law.

    Reported means are arithmetic means over the supplied measurement points,
    matching the official Bosch 89-point outcome convention.  They are not
    asserted to be continuous wafer-area means.
    """
    boundaries = tuple(boundaries)
    if not boundaries or any(
            not isinstance(item, BoschSPTSWaferBoundaryMapTrace)
            for item in boundaries):
        raise TypeError("a nonempty cylindrical Bosch boundary batch is required")
    names = ("F", "C4F8_film_precursor", "positive_ion")
    x = boundaries[0].x_m
    y = boundaries[0].y_m
    point_count = x.size
    if any(
            item.species_names != names
            or item.x_m.size != point_count
            or not np.array_equal(item.x_m, x)
            or not np.array_equal(item.y_m, y)
            for item in boundaries):
        raise ValueError("Bosch point boundaries do not share one coordinate map")

    batch = len(boundaries)
    maximum_intervals = max(
        item.reactor.interval_duration_s.size for item in boundaries)
    flux = np.zeros((maximum_intervals, batch, len(names), point_count))
    energy = np.zeros((maximum_intervals, batch, 1))
    duration = np.zeros((maximum_intervals, batch, 1))
    for wafer, boundary in enumerate(boundaries):
        count = boundary.reactor.interval_duration_s.size
        flux[:count, wafer] = boundary.point_flux_m2_s
        energy[:count, wafer, 0] = boundary.reactor.ion_energy_eV
        duration[:count, wafer, 0] = boundary.reactor.interval_duration_s

    removed_si, removed_oxide, film_inventory, silicon, oxide = (
        _advance_fused_surface_recurrence(
            flux, energy, duration, silicon_mechanism, oxide_mechanism))

    predictions = []
    for wafer, boundary in enumerate(boundaries):
        silicon_depth = removed_si[wafer] / silicon.bulk_si_atom_density_m3
        oxide_loss = removed_oxide[wafer] / oxide.bulk_formula_density_m3
        mean_si = float(np.mean(silicon_depth))
        mean_oxide = float(np.mean(oxide_loss))
        predictions.append(BoschWaferPointDepthPrediction(
            experiment_key=boundary.reactor.trace_key,
            x_m=boundary.x_m,
            y_m=boundary.y_m,
            silicon_depth_m=silicon_depth,
            oxide_loss_m=oxide_loss,
            remaining_film_units_m2=film_inventory[wafer],
            measurement_point_mean_silicon_depth_m=mean_si,
            measurement_point_mean_oxide_loss_m=mean_oxide,
            silicon_to_oxide_selectivity=(
                mean_si / mean_oxide if mean_oxide > 0.0 else 0.0),
            maximum_material_ledger_relative_residual=0.0,
            provenance={
                "model": "petch-spts-bosch-wafer-point-depth-fused-batch-v1",
                "surface_recurrence": (
                    "exact fused algebraic parity path for unchanged Belen and "
                    "La Magna mechanisms"),
                "reactor_boundary_model": boundary.provenance["model"],
                "silicon_surface_model": silicon_mechanism.provenance["model"],
                "oxide_surface_model": oxide_mechanism.provenance["model"],
                "mean_definition": "arithmetic mean over supplied measurement points",
                "target_depth_used": False,
                "surface_parameters_shared_across_every_wafer": True,
            }))
    return tuple(predictions)
