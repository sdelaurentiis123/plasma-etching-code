"""Common-engine SPTS Bosch wafer-depth and selectivity prediction.

This layer consumes the measured-waveform reactor/wafer boundary and advances
the same surface laws used by the feature engine.  Every radial annulus is
advanced in one vectorized deterministic calculation; no depth regression or
wafer-specific multiplier appears here.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
import math

import numpy as np

from .bosch_silicon import BoschSiliconFluorocarbonMechanism
from .fluorocarbon_lamagna import (
    LaMagnaFluorocarbonParameters, LaMagnaGarozzoFluorocarbonMechanism,
)
from .reactor_global.bosch_spts_reduced import BoschSPTSWaferBoundaryTrace
from .silicon_sf6o2 import belen_2005_reference_silicon_mechanism
from .surface_kinetics import EnergeticFlux, SurfaceFluxes


@dataclass(frozen=True)
class BoschWaferDepthPrediction:
    experiment_key: str
    radial_centers_m: np.ndarray
    radial_area_m2: np.ndarray
    silicon_depth_m: np.ndarray
    oxide_loss_m: np.ndarray
    remaining_film_units_m2: np.ndarray
    wafer_mean_silicon_depth_m: float
    wafer_mean_oxide_loss_m: float
    silicon_to_oxide_selectivity: float
    maximum_material_ledger_relative_residual: float
    provenance: dict

    def __post_init__(self):
        radial = np.asarray(self.radial_centers_m, dtype=float).copy()
        area = np.asarray(self.radial_area_m2, dtype=float).copy()
        silicon = np.asarray(self.silicon_depth_m, dtype=float).copy()
        oxide = np.asarray(self.oxide_loss_m, dtype=float).copy()
        film = np.asarray(self.remaining_film_units_m2, dtype=float).copy()
        if (not self.experiment_key or radial.ndim != 1 or radial.size < 4
                or any(value.shape != radial.shape for value in (
                    area, silicon, oxide, film))
                or any(np.any(~np.isfinite(value)) for value in (
                    radial, area, silicon, oxide, film))
                or np.any(radial <= 0.0) or np.any(area < 0.0)
                or np.any(silicon < 0.0) or np.any(oxide < 0.0)
                or np.any(film < 0.0) or np.sum(area) <= 0.0
                or not math.isfinite(self.wafer_mean_silicon_depth_m)
                or self.wafer_mean_silicon_depth_m < 0.0
                or not math.isfinite(self.wafer_mean_oxide_loss_m)
                or self.wafer_mean_oxide_loss_m < 0.0
                or not math.isfinite(self.silicon_to_oxide_selectivity)
                or self.silicon_to_oxide_selectivity < 0.0
                or not math.isfinite(self.maximum_material_ledger_relative_residual)
                or not 0.0 <= self.maximum_material_ledger_relative_residual < 1.0e-10):
            raise ValueError("invalid Bosch wafer-depth prediction")
        for value in (radial, area, silicon, oxide, film):
            value.setflags(write=False)
        object.__setattr__(self, "radial_centers_m", radial)
        object.__setattr__(self, "radial_area_m2", area)
        object.__setattr__(self, "silicon_depth_m", silicon)
        object.__setattr__(self, "oxide_loss_m", oxide)
        object.__setattr__(self, "remaining_film_units_m2", film)
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def build_bosch_reference_surface_mechanisms(
        *, reference_etchant_flux_m2_s=1.0e21,
        fluorine_sticking_probability=0.5):
    """Return unchanged reference Si and SiO2 laws for the Bosch boundary.

    ``reference_etchant_flux_m2_s`` is the La Magna chemical-rate reference,
    not a fitted incident flux.  It is fixed across every wafer and both the
    Si film clock and oxide mask use the same parameter object.
    """
    if (not np.isfinite(reference_etchant_flux_m2_s)
            or reference_etchant_flux_m2_s <= 0.0):
        raise ValueError("reference etchant flux must be positive")
    film_parameters = LaMagnaFluorocarbonParameters.viennaps_4_6_1_reference(
        reference_etchant_flux_m2_s=reference_etchant_flux_m2_s,
        etchant_species=("F",),
        polymer_species=("C4F8_film_precursor",),
        projectile_species=("positive_ion",),
        neutral_transport_mode="species_specific")
    oxide = LaMagnaGarozzoFluorocarbonMechanism(film_parameters)
    silicon = BoschSiliconFluorocarbonMechanism(
        belen_2005_reference_silicon_mechanism(
            fluorine_sticking_probability=fluorine_sticking_probability,
            projectile_species=("positive_ion",)),
        LaMagnaGarozzoFluorocarbonMechanism(film_parameters))
    return silicon, oxide


def predict_bosch_wafer_depth(
        boundary: BoschSPTSWaferBoundaryTrace,
        silicon_mechanism: BoschSiliconFluorocarbonMechanism,
        oxide_mechanism: LaMagnaGarozzoFluorocarbonMechanism):
    if not isinstance(boundary, BoschSPTSWaferBoundaryTrace):
        raise TypeError("boundary must be BoschSPTSWaferBoundaryTrace")
    if not isinstance(silicon_mechanism, BoschSiliconFluorocarbonMechanism):
        raise TypeError("silicon mechanism must be the Bosch composite")
    if not isinstance(oxide_mechanism, LaMagnaGarozzoFluorocarbonMechanism):
        raise TypeError("oxide mechanism must be La Magna/Garozzo")

    species_index = {name: index for index, name in enumerate(boundary.species_names)}
    if set(species_index) != {"F", "C4F8_film_precursor", "positive_ion"}:
        raise ValueError("unexpected Bosch wafer-boundary species")
    radial_count = boundary.radial_centers_m.size
    silicon_state = silicon_mechanism.initial_state((radial_count,))
    oxide_state = oxide_mechanism.initial_state((radial_count,))
    maximum_ledger = 0.0

    for interval, duration in enumerate(boundary.reactor.interval_duration_s):
        radial = boundary.radial_flux_m2_s[interval]
        fluxes = SurfaceFluxes(
            {
                "F": radial[species_index["F"]],
                "C4F8_film_precursor": radial[
                    species_index["C4F8_film_precursor"]],
            },
            (EnergeticFlux(
                "positive_ion", radial[species_index["positive_ion"]],
                [boundary.reactor.ion_energy_eV[interval]], [1.0], [1.0]),))
        silicon_step = silicon_mechanism.advance(
            silicon_state, fluxes, float(duration), strict=True)
        oxide_step = oxide_mechanism.advance(
            oxide_state, fluxes, float(duration), strict=True)
        silicon_state = silicon_step.state
        oxide_state = oxide_step.state
        for exchange in (
                silicon_step.material_exchange, oxide_step.material_exchange):
            for inventory, removed in exchange.removed_units_m2.items():
                residual = np.abs(exchange.residual_units_m2(inventory))
                scale = np.maximum(np.asarray(removed, dtype=float), 1.0)
                maximum_ledger = max(
                    maximum_ledger, float(np.max(residual / scale)))

    silicon_density = silicon_mechanism.silicon_mechanism.parameters.bulk_si_atom_density_m3
    oxide_density = oxide_mechanism.parameters.bulk_formula_density_m3
    silicon_depth = silicon_state.removed_si_atoms_m2 / silicon_density
    oxide_loss = oxide_state.removed_formula_units_m2 / oxide_density

    parameters = boundary.reactor.provenance["parameters"]
    wafer_radius = float(parameters["wafer_radius_m"])
    reactor_radius = float(parameters["reactor_radius_m"])
    radial_count = boundary.radial_centers_m.size
    edges = np.linspace(0.0, reactor_radius, radial_count + 1)
    clipped_outer = np.minimum(edges[1:], wafer_radius)
    clipped_inner = np.minimum(edges[:-1], wafer_radius)
    area = np.pi * np.maximum(clipped_outer ** 2 - clipped_inner ** 2, 0.0)
    wafer_area = np.sum(area)
    mean_si = float(np.dot(silicon_depth, area) / wafer_area)
    mean_oxide = float(np.dot(oxide_loss, area) / wafer_area)
    selectivity = mean_si / mean_oxide if mean_oxide > 0.0 else 0.0
    return BoschWaferDepthPrediction(
        experiment_key=boundary.reactor.trace_key,
        radial_centers_m=boundary.radial_centers_m,
        radial_area_m2=area,
        silicon_depth_m=silicon_depth,
        oxide_loss_m=oxide_loss,
        remaining_film_units_m2=silicon_state.polymer_film_units_m2,
        wafer_mean_silicon_depth_m=mean_si,
        wafer_mean_oxide_loss_m=mean_oxide,
        silicon_to_oxide_selectivity=selectivity,
        maximum_material_ledger_relative_residual=maximum_ledger,
        provenance={
            "model": "petch-spts-bosch-wafer-depth-v1",
            "reactor_boundary_model": boundary.reactor.provenance["model"],
            "silicon_surface_model": silicon_mechanism.provenance["model"],
            "oxide_surface_model": oxide_mechanism.provenance["model"],
            "target_depth_used": False,
            "surface_parameters_shared_across_every_wafer": True,
        })
