"""Deterministic species-resolved reactor-to-feature boundary construction.

This adapter preserves absolute species flux, ion charge/mass, joint ion
energy-angle structure, and thermal-neutral identities.  It performs no
particle Monte Carlo and no chemistry aggregation.  A reactor calculation can
therefore feed a microscopic surface deck without first collapsing all ions
into an anonymous projectile or all radicals into one scalar channel.
"""
from __future__ import annotations

import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import (
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
)
from .iadf_two_component import TwoComponentIADF, build_two_component_boundary


BOLTZMANN_EV_K = 8.617333262145e-5


def _finite_mapping(values, *, positive: bool, label: str):
    converted = {str(name): float(value) for name, value in dict(values).items()}
    if (
        not converted
        or any(not name for name in converted)
        or any(
            not math.isfinite(value)
            or (value <= 0.0 if positive else value < 0.0)
            for value in converted.values()
        )
    ):
        raise ValueError(f"invalid {label}")
    return MappingProxyType(converted)


def _charge_mapping(values, names):
    converted = {str(name): int(value) for name, value in dict(values).items()}
    if (
        set(converted) != set(names)
        or any(value <= 0 for value in converted.values())
    ):
        raise ValueError("invalid ion charge-number mapping")
    return MappingProxyType(converted)


def _energy_measures(energy_eV, energy_weight, names):
    energy = {}
    weight = {}
    supplied_weight = {} if energy_weight is None else dict(energy_weight)
    if supplied_weight and set(supplied_weight) != set(names):
        raise ValueError("ion energy weights do not cover the ion species")
    if set(energy_eV) != set(names):
        raise ValueError("ion energy measures do not cover the ion species")
    for name in names:
        values = np.atleast_1d(np.asarray(energy_eV[name], dtype=float)).copy()
        if values.ndim != 1 or values.size == 0 or np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(f"invalid ion energy measure for {name}")
        if supplied_weight:
            weights = np.asarray(supplied_weight[name], dtype=float).copy()
        else:
            weights = np.full(values.shape, 1.0 / values.size)
        if (
            weights.shape != values.shape
            or np.any(~np.isfinite(weights))
            or np.any(weights < 0.0)
            or weights.sum() <= 0.0
        ):
            raise ValueError(f"invalid ion energy weights for {name}")
        weights /= weights.sum()
        values.setflags(write=False)
        weights.setflags(write=False)
        energy[name] = values
        weight[name] = weights
    return MappingProxyType(energy), MappingProxyType(weight)


def _thermal_neutral_quadrature(temperature_eV, n_transverse, n_normal):
    hermite_node, hermite_weight = np.polynomial.hermite.hermgauss(
        int(n_transverse)
    )
    laguerre_node, laguerre_weight = np.polynomial.laguerre.laggauss(
        int(n_normal)
    )
    ix, iy, iz = np.meshgrid(
        np.arange(hermite_node.size),
        np.arange(hermite_node.size),
        np.arange(laguerre_node.size),
        indexing="ij",
    )
    velocity = np.column_stack((
        np.sqrt(temperature_eV) * hermite_node[ix.ravel()],
        np.sqrt(temperature_eV) * hermite_node[iy.ravel()],
        np.sqrt(temperature_eV * laguerre_node[iz.ravel()]),
    ))
    weight = (
        hermite_weight[ix.ravel()]
        * hermite_weight[iy.ravel()]
        * laguerre_weight[iz.ravel()]
        / np.pi
    )
    return velocity, weight


def build_species_resolved_feature_boundary(
    *,
    ion_flux_m2_s: Mapping[str, float],
    ion_mass_amu: Mapping[str, float],
    ion_charge_number: Mapping[str, int],
    ion_energy_eV: Mapping[str, float | np.ndarray],
    ion_iadf: TwoComponentIADF,
    neutral_flux_m2_s: Mapping[str, float],
    neutral_mass_amu: Mapping[str, float],
    neutral_temperature_K: float,
    reference_plane_m: float,
    ion_energy_weight: Mapping[str, np.ndarray] | None = None,
    ion_polar_order: int = 12,
    ion_azimuthal_order: int = 8,
    neutral_transverse_order: int = 3,
    neutral_normal_order: int = 4,
    provenance: Mapping[str, object] | None = None,
) -> PlasmaBoundaryState:
    """Build one absolute, deterministic reactor-to-feature boundary.

    Ion energy measures may be scalar, discrete RF-phase distributions, or
    deterministic collisional-sheath quadratures.  The supplied IADF supplies
    the conditional angular law at every energy node.  Neutral fluxes are
    one-way thermal wall fluxes with analytic half-Maxwellian densities.
    """
    ions = _finite_mapping(ion_flux_m2_s, positive=False, label="ion flux")
    ion_mass = _finite_mapping(ion_mass_amu, positive=True, label="ion mass")
    charges = _charge_mapping(ion_charge_number, ions)
    if set(ion_mass) != set(ions):
        raise ValueError("ion masses do not cover the ion species")
    energies, energy_weights = _energy_measures(
        dict(ion_energy_eV), ion_energy_weight, ions
    )
    neutrals = _finite_mapping(
        neutral_flux_m2_s, positive=False, label="neutral flux"
    )
    neutral_mass = _finite_mapping(
        neutral_mass_amu, positive=True, label="neutral mass"
    )
    if set(neutral_mass) != set(neutrals) or set(ions) & set(neutrals):
        raise ValueError("neutral identities are incomplete or overlap ions")
    if not isinstance(ion_iadf, TwoComponentIADF):
        raise TypeError("ion_iadf must be a TwoComponentIADF")
    values = np.asarray([
        neutral_temperature_K,
        reference_plane_m,
        ion_polar_order,
        ion_azimuthal_order,
        neutral_transverse_order,
        neutral_normal_order,
    ], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or neutral_temperature_K <= 0.0
        or int(ion_polar_order) < 2
        or int(ion_azimuthal_order) < 1
        or int(neutral_transverse_order) < 1
        or int(neutral_normal_order) < 1
    ):
        raise ValueError("invalid feature-boundary quadrature closure")

    species = []
    for name in sorted(ions):
        temporary = build_two_component_boundary(
            ion_iadf,
            ions[name],
            energies[name],
            energy_weight=energy_weights[name],
            ion_mass_amu=ion_mass[name],
            name=name,
            n_polar=int(ion_polar_order),
            azimuthal_order=int(ion_azimuthal_order),
            reference_plane_m=float(reference_plane_m),
            extra_provenance={
                "role": "reactor_species_resolved_positive_ion",
                "charge_number": charges[name],
                "energy_measure_node_count": int(energies[name].size),
            },
        ).species[0]
        species.append(SpeciesBoundaryState(
            name=temporary.name,
            charge_number=charges[name],
            mass_amu=temporary.mass_amu,
            flux_m2_s=temporary.flux_m2_s,
            velocity_sqrt_eV=temporary.velocity_sqrt_eV,
            weight=temporary.weight,
            phase_rad=temporary.phase_rad,
            position_m=temporary.position_m,
            density_model=temporary.density_model,
            density_model_2d=temporary.density_model_2d,
            provenance=temporary.provenance,
        ))

    neutral_temperature_eV = BOLTZMANN_EV_K * float(neutral_temperature_K)
    neutral_velocity, neutral_weight = _thermal_neutral_quadrature(
        neutral_temperature_eV,
        int(neutral_transverse_order),
        int(neutral_normal_order),
    )
    for name in sorted(neutrals):
        species.append(SpeciesBoundaryState(
            name=name,
            charge_number=0,
            mass_amu=neutral_mass[name],
            flux_m2_s=neutrals[name],
            velocity_sqrt_eV=neutral_velocity,
            weight=neutral_weight,
            density_model=MaxwellianFluxVelocityDensity(
                neutral_temperature_eV
            ),
            provenance={
                "role": "reactor_species_resolved_thermal_neutral",
                "temperature_K": float(neutral_temperature_K),
                "one_way_thermal_flux_supplied_by_reactor": True,
            },
        ))
    metadata = {
        "provider": "species_resolved_feature_boundary",
        "deterministic_quadrature": True,
        "monte_carlo": False,
        "ion_species_count": len(ions),
        "neutral_species_count": len(neutrals),
        "positive_ion_total_flux_m2_s": float(sum(ions.values())),
        "neutral_total_flux_m2_s": float(sum(neutrals.values())),
        "neutral_temperature_K": float(neutral_temperature_K),
    }
    if provenance:
        metadata.update(dict(provenance))
    return PlasmaBoundaryState(
        species=tuple(species),
        reference_plane_m=float(reference_plane_m),
        provenance=metadata,
    )


__all__ = ["build_species_resolved_feature_boundary"]
