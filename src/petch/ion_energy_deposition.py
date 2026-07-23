"""Analytic ion energy deposition for derived (knob-free) yield laws.

Nuclear stopping follows the ZBL universal reduced-stopping fit (Ziegler, Biersack,
Littmark 1985); electronic stopping follows Lindhard-Scharff velocity-proportional
stopping.  Compound targets combine by Bragg additivity.  Everything downstream is
integrals of these two curves: CSDA path length, projected range, and -- the quantity
the yield law actually needs -- the NUCLEAR energy deposited within the reactive
surface layer along a slant path.  Inputs are atomic numbers, masses, and densities;
there are no adjustable parameters in this module.

Accuracy expectation at sub-10-keV energies is the ZBL analytic level (~10-20 percent
against full SRIM transport); the module carries validation anchors as tests rather
than claiming transport-code fidelity.  Iteration 1 uses the straight-slant-path
approximation for the depth mapping; cascade straggling and detour corrections are
declared refinements, not silent assumptions.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos, log, radians, sqrt

import numpy as np

# eV cm^2 (per atom) conversion prefactor of the ZBL reduced nuclear stopping.
_ZBL_PREFACTOR_EV_CM2 = 8.462e-15


@dataclass(frozen=True)
class TargetComponent:
    atomic_number: int
    mass_amu: float
    stoichiometry: float


@dataclass(frozen=True)
class Target:
    """Amorphous compound target with Bragg-additive stopping."""

    components: tuple[TargetComponent, ...]
    atom_density_m3: float

    def __post_init__(self):
        if (not self.components
                or not np.isfinite(self.atom_density_m3)
                or self.atom_density_m3 <= 0.0):
            raise ValueError("invalid stopping target")


SIO2 = Target(
    components=(TargetComponent(14, 28.086, 1.0), TargetComponent(8, 15.999, 2.0)),
    atom_density_m3=6.6e28)
AMORPHOUS_CARBON = Target(
    components=(TargetComponent(6, 12.011, 1.0),),
    atom_density_m3=1.0e29)
FLUOROCARBON_FILM = Target(
    components=(TargetComponent(6, 12.011, 1.0), TargetComponent(9, 18.998, 1.5)),
    atom_density_m3=7.5e28)


def _reduced_energy(energy_eV, z1, m1, z2, m2):
    return (32.53 * m2 * (energy_eV / 1000.0)
            / (z1 * z2 * (m1 + m2) * (z1 ** 0.23 + z2 ** 0.23)))


def _zbl_reduced_nuclear(reduced):
    reduced = np.asarray(reduced, dtype=float)
    small = reduced <= 30.0
    result = np.empty_like(reduced)
    with np.errstate(divide="ignore", invalid="ignore"):
        result[small] = (
            np.log1p(1.1383 * reduced[small])
            / (2.0 * (reduced[small]
                      + 0.01321 * reduced[small] ** 0.21226
                      + 0.19593 * np.sqrt(reduced[small]))))
        result[~small] = np.log(reduced[~small]) / (2.0 * reduced[~small])
    result[reduced <= 0.0] = 0.0
    return result


def _lindhard_k(z1, m1, z2, m2):
    return (0.0793 * z1 ** (2.0 / 3.0) * sqrt(z2) * (m1 + m2) ** 1.5
            / ((z1 ** (2.0 / 3.0) + z2 ** (2.0 / 3.0)) ** 0.75
               * m1 ** 1.5 * sqrt(m2)))


def stopping_cross_sections_eV_cm2(energy_eV, z1, m1, target):
    """Per-atom nuclear and electronic stopping cross sections, Bragg-additive."""
    energy = np.asarray(energy_eV, dtype=float)
    nuclear = np.zeros_like(energy)
    electronic = np.zeros_like(energy)
    total_stoich = sum(c.stoichiometry for c in target.components)
    for component in target.components:
        z2, m2 = component.atomic_number, component.mass_amu
        reduced = _reduced_energy(energy, z1, m1, z2, m2)
        conversion = (_ZBL_PREFACTOR_EV_CM2 * z1 * z2 * m1
                      / ((m1 + m2) * (z1 ** 0.23 + z2 ** 0.23)))
        weight = component.stoichiometry / total_stoich
        nuclear += weight * conversion * _zbl_reduced_nuclear(reduced)
        k = _lindhard_k(z1, m1, z2, m2)
        # d(reduced)/d(path) conversion reused for the electronic channel.
        electronic += weight * conversion * k * np.sqrt(np.maximum(reduced, 0.0))
    return nuclear, electronic


def csda_path_nm(energy_eV, z1, m1, target, *, minimum_energy_eV=10.0, steps=400):
    """Continuous-slowing-down path length by trapezoid integration of 1/(N S)."""
    energy = float(energy_eV)
    if energy <= minimum_energy_eV:
        return 0.0
    grid = np.geomspace(minimum_energy_eV, energy, int(steps))
    nuclear, electronic = stopping_cross_sections_eV_cm2(grid, z1, m1, target)
    total = (nuclear + electronic) * (target.atom_density_m3 * 1e-6)  # eV/cm
    path_cm = np.trapz(1.0 / total, grid)
    return float(path_cm * 1e7)


def projected_range_nm(energy_eV, z1, m1, target):
    """Projected range from CSDA path with the standard mass-ratio detour factor."""
    mean_m2 = (sum(c.stoichiometry * c.mass_amu for c in target.components)
               / sum(c.stoichiometry for c in target.components))
    detour = 1.0 + mean_m2 / (3.0 * m1)
    return csda_path_nm(energy_eV, z1, m1, target) / detour


def nuclear_energy_in_layer_eV(energy_eV, cosine_incidence, layer_depth_nm, z1, m1,
                               target, *, minimum_energy_eV=10.0, steps=300):
    """Nuclear energy (eV) deposited within the top layer along a straight slant path.

    The ion enters at incidence cosine ``mu`` and slows along a straight path; the
    depth coordinate advances as path * mu.  Nuclear losses within depth <=
    layer_depth_nm are the chemistry-driving deposit (Sigmund); electronic losses and
    deeper nuclear losses are counted as buried.  Straight-path is the declared
    iteration-1 approximation.
    """
    energy = float(energy_eV)
    mu = float(np.clip(cosine_incidence, 1e-3, 1.0))
    if energy <= minimum_energy_eV or layer_depth_nm <= 0.0:
        return 0.0
    path_limit_cm = (layer_depth_nm * 1e-7) / mu
    density_cm3 = target.atom_density_m3 * 1e-6
    grid = np.geomspace(minimum_energy_eV, energy, int(steps))[::-1]
    deposited_nuclear = 0.0
    travelled = 0.0
    for index in range(len(grid) - 1):
        e_high, e_low = grid[index], grid[index + 1]
        e_mid = 0.5 * (e_high + e_low)
        nuclear, electronic = stopping_cross_sections_eV_cm2(
            np.array([e_mid]), z1, m1, target)
        total = float(nuclear[0] + electronic[0]) * density_cm3
        segment = (e_high - e_low) / total
        if travelled + segment >= path_limit_cm:
            fraction = max(0.0, (path_limit_cm - travelled) / segment)
            deposited_nuclear += (e_high - e_low) * fraction * (
                float(nuclear[0]) / float(nuclear[0] + electronic[0]))
            return float(deposited_nuclear)
        travelled += segment
        deposited_nuclear += (e_high - e_low) * (
            float(nuclear[0]) / float(nuclear[0] + electronic[0]))
    return float(deposited_nuclear)


def derived_yield_energy_factor(energy_eV, cosine_incidence, *, layer_depth_nm,
                                reference_energy_eV, reference_cosine=1.0,
                                threshold_energy_eV=0.0, z1=18, m1=39.948,
                                target=SIO2):
    """Energy-angle factor proportional to nuclear energy deposited in the layer,
    normalized to unity at the reference energy and normal incidence.

    Anchoring at the reference keeps the mechanism's absolute yield table intact:
    this factor replaces only the EXTRAPOLATION shape (the fitted knee and part of
    the empirical angular law), not the calibrated magnitude.
    """
    if threshold_energy_eV and energy_eV <= threshold_energy_eV:
        return 0.0
    reference = nuclear_energy_in_layer_eV(
        reference_energy_eV, reference_cosine, layer_depth_nm, z1, m1, target)
    if reference <= 0.0:
        raise ValueError("reference energy deposits nothing in the layer")
    value = nuclear_energy_in_layer_eV(
        energy_eV, cosine_incidence, layer_depth_nm, z1, m1, target)
    return float(value / reference)
