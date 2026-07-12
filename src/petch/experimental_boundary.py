"""Adapters from experimental diagnostics to the common plasma-boundary contract.

Experimental measurements and missing closures remain separate. In particular, a self-bias voltage
is not silently converted into an IEDF, and an integrated radical/ion ratio is not silently promoted
to a species-resolved flux vector.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import (
    MaxwellianFluxVelocityDensity, PlasmaBoundaryState, SpeciesBoundaryState,
)
from .experimental_data import (
    Jeon2022ElectronBiasControl, Jeon2022PlasmaControl,
    jeon_2022_bohm_ion_flux_m2_s,
)


@dataclass(frozen=True)
class Jeon2022BoundaryClosure:
    ion_name: str
    ion_mass_amu: float
    ion_normal_energy_eV: np.ndarray
    ion_normal_energy_weight: np.ndarray
    ion_tangential_temperature_eV: float
    neutral_flux_fraction: Mapping[str, float]
    neutral_mass_amu: Mapping[str, float]
    neutral_temperature_K: float
    provenance: Mapping[str, object]
    supports_prediction_within_declared_domain: bool = False

    def __post_init__(self):
        energy = np.asarray(self.ion_normal_energy_eV, dtype=float).copy()
        weight = np.asarray(self.ion_normal_energy_weight, dtype=float).copy()
        fraction = dict(self.neutral_flux_fraction)
        mass = dict(self.neutral_mass_amu)
        if (not self.ion_name or self.ion_mass_amu <= 0.0 or energy.ndim != 1
                or weight.shape != energy.shape or energy.size == 0
                or np.any(~np.isfinite(energy)) or np.any(energy < 0.0)
                or np.any(~np.isfinite(weight)) or np.any(weight < 0.0) or weight.sum() <= 0.0
                or self.ion_tangential_temperature_eV <= 0.0
                or self.neutral_temperature_K <= 0.0 or not fraction
                or set(fraction) != set(mass) or any(not name for name in fraction)
                or any(not np.isfinite(value) or value < 0.0 for value in fraction.values())
                or any(not np.isfinite(value) or value <= 0.0 for value in mass.values())
                or not np.isclose(sum(fraction.values()), 1.0, rtol=0.0, atol=2e-13)
                or not isinstance(self.supports_prediction_within_declared_domain, bool)):
            raise ValueError("invalid Jeon boundary closure")
        weight /= weight.sum()
        energy.setflags(write=False); weight.setflags(write=False)
        object.__setattr__(self, "ion_normal_energy_eV", energy)
        object.__setattr__(self, "ion_normal_energy_weight", weight)
        object.__setattr__(self, "neutral_flux_fraction", MappingProxyType(fraction))
        object.__setattr__(self, "neutral_mass_amu", MappingProxyType(mass))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def _same_jeon_condition(left, right):
    return (left.condition_family == right.condition_family
            and left.c4f8_fraction == right.c4f8_fraction
            and left.pulse_off_ms == right.pulse_off_ms)


def build_jeon_2022_boundary_state(
        plasma_control: Jeon2022PlasmaControl,
        electron_bias_control: Jeon2022ElectronBiasControl,
        closure: Jeon2022BoundaryClosure, *, reference_plane_m, n_transverse_ion=3,
        n_transverse_neutral=5, n_normal_neutral=8):
    """Build one dimensional boundary state from measurements plus an explicit closure."""
    if (not isinstance(plasma_control, Jeon2022PlasmaControl)
            or not isinstance(electron_bias_control, Jeon2022ElectronBiasControl)
            or not isinstance(closure, Jeon2022BoundaryClosure)):
        raise TypeError("Jeon boundary construction requires typed evidence and closure")
    if not _same_jeon_condition(plasma_control, electron_bias_control):
        raise ValueError("Jeon plasma controls refer to different experimental conditions")
    if not np.isfinite(reference_plane_m):
        raise ValueError("reference_plane_m must be finite")
    if int(n_transverse_ion) <= 0 or int(n_transverse_neutral) <= 0 or int(n_normal_neutral) <= 0:
        raise ValueError("boundary quadrature orders must be positive")

    ion_flux = jeon_2022_bohm_ion_flux_m2_s(electron_bias_control)
    nodes, gh_weight = np.polynomial.hermite.hermgauss(int(n_transverse_ion))
    transverse = np.sqrt(closure.ion_tangential_temperature_eV) * nodes
    transverse_weight = gh_weight / np.sqrt(np.pi)
    ix, iy, iz = np.meshgrid(
        np.arange(nodes.size), np.arange(nodes.size),
        np.arange(closure.ion_normal_energy_eV.size), indexing="ij")
    ion_velocity = np.column_stack((
        transverse[ix.ravel()], transverse[iy.ravel()],
        np.sqrt(closure.ion_normal_energy_eV[iz.ravel()])))
    ion_weight = (transverse_weight[ix.ravel()] * transverse_weight[iy.ravel()]
                  * closure.ion_normal_energy_weight[iz.ravel()])
    ion = SpeciesBoundaryState(
        closure.ion_name, 1, closure.ion_mass_amu, ion_flux,
        ion_velocity, ion_weight,
        provenance={
            "role": "explicit_iedf_closure",
            "closure": dict(closure.provenance),
            "supports_prediction": closure.supports_prediction_within_declared_domain,
        })

    temperature_eV = closure.neutral_temperature_K * 8.617333262145e-5
    hermite_node, hermite_weight = np.polynomial.hermite.hermgauss(int(n_transverse_neutral))
    laguerre_node, laguerre_weight = np.polynomial.laguerre.laggauss(int(n_normal_neutral))
    nx, ny, nz = np.meshgrid(
        np.arange(hermite_node.size), np.arange(hermite_node.size),
        np.arange(laguerre_node.size), indexing="ij")
    neutral_velocity = np.column_stack((
        np.sqrt(temperature_eV) * hermite_node[nx.ravel()],
        np.sqrt(temperature_eV) * hermite_node[ny.ravel()],
        np.sqrt(temperature_eV * laguerre_node[nz.ravel()])))
    neutral_weight = (hermite_weight[nx.ravel()] * hermite_weight[ny.ravel()]
                      * laguerre_weight[nz.ravel()] / np.pi)
    total_neutral_flux = plasma_control.neutral_to_ion_flux_ratio * ion_flux
    neutral_species = tuple(SpeciesBoundaryState(
        name, 0, closure.neutral_mass_amu[name], total_neutral_flux * fraction,
        neutral_velocity, neutral_weight,
        density_model=MaxwellianFluxVelocityDensity(temperature_eV),
        provenance={
            "role": "species_composition_closure",
            "integrated_flux_source": plasma_control.source_location,
            "closure": dict(closure.provenance),
            "supports_prediction": closure.supports_prediction_within_declared_domain,
        }) for name, fraction in closure.neutral_flux_fraction.items())
    return PlasmaBoundaryState(
        (ion,) + neutral_species, float(reference_plane_m),
        provenance={
            "source": "Jeong_2022_diagnostics_plus_explicit_closure",
            "electron_density": electron_bias_control.source_location,
            "self_bias_energy_scale_v": electron_bias_control.self_bias_magnitude_v,
            "self_bias_is_not_iedf": True,
            "integrated_neutral_to_ion_ratio": plasma_control.neutral_to_ion_flux_ratio,
            "closure_supports_prediction": closure.supports_prediction_within_declared_domain,
        })
