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

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# eV cm^2 (per atom) conversion prefactor of the ZBL reduced nuclear stopping.
_ZBL_PREFACTOR_EV_CM2 = 8.462e-15
_RESIDUAL_PATH_TABLE_CACHE = {}


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
    path_cm = _trapezoid(1.0 / total, grid)
    return float(path_cm * 1e7)


def _csda_path_table(z1, m1, target, minimum_energy_eV, maximum_energy_eV,
                     points):
    """Return a dense monotone ``energy -> CSDA path`` table.

    The table is the cumulative form of the same quadrature used by
    :func:`csda_path_nm`.  It exists so residual-energy inversion is a
    deterministic interpolation problem rather than a separate fitted
    attenuation law.
    """
    maximum = max(float(maximum_energy_eV), float(minimum_energy_eV) * 1.001)
    # Round the cache ceiling upward by decades so repeated calls with nearby
    # event energies share one table without clipping their support.
    decade_ceiling = 10.0 ** np.ceil(np.log10(maximum))
    key = (int(z1), float(m1), target, float(minimum_energy_eV),
           float(decade_ceiling), int(points))
    cached = _RESIDUAL_PATH_TABLE_CACHE.get(key)
    if cached is not None:
        return cached

    energy = np.geomspace(minimum_energy_eV, decade_ceiling, int(points))
    nuclear, electronic = stopping_cross_sections_eV_cm2(
        energy, z1, m1, target)
    stopping_eV_cm = (
        (nuclear + electronic) * (target.atom_density_m3 * 1e-6))
    inverse_stopping = 1.0 / stopping_eV_cm
    increments_cm = (
        0.5 * (inverse_stopping[1:] + inverse_stopping[:-1])
        * np.diff(energy))
    path_nm = np.concatenate((
        np.zeros(1, dtype=float),
        np.cumsum(increments_cm, dtype=float) * 1e7,
    ))
    cached = (energy, path_nm)
    _RESIDUAL_PATH_TABLE_CACHE[key] = cached
    return cached


def residual_energy_after_layer_eV(
        energy_eV, cosine_incidence, layer_depth_nm, z1, m1, target,
        *, minimum_energy_eV=10.0, table_points=4096):
    """Ion energy remaining after crossing a layer by CSDA path inversion.

    ``layer_depth_nm`` is the layer thickness along its surface normal.
    The traversed material path is therefore ``depth / cosine_incidence``.
    Starting from the CSDA path-to-rest at the incident energy, this function
    subtracts that material path and inverts the same ZBL/Lindhard stopping
    integral.  It has no attenuation length or fitted transmission constant.

    The analytic stopping model is not extended below ``minimum_energy_eV``:
    when the residual path enters that unresolved end-of-range interval, the
    ion is declared stopped and the returned energy is zero.  A zero-thickness
    layer is an exact identity, including below the cutoff.
    """
    energy = np.asarray(energy_eV, dtype=float)
    cosine = np.asarray(cosine_incidence, dtype=float)
    depth = np.asarray(layer_depth_nm, dtype=float)
    energy_b, cosine_b, depth_b = np.broadcast_arrays(energy, cosine, depth)
    if (np.any(~np.isfinite(energy_b))
            or np.any(~np.isfinite(cosine_b))
            or np.any(~np.isfinite(depth_b))):
        raise ValueError("energy, incidence cosine, and layer depth must be finite")
    if np.any(energy_b < 0.0):
        raise ValueError("ion energy must be nonnegative")
    if np.any((cosine_b < 0.0) | (cosine_b > 1.0)):
        raise ValueError("incidence cosine must lie in [0, 1]")
    if np.any(depth_b < 0.0):
        raise ValueError("layer depth must be nonnegative")
    if minimum_energy_eV <= 0.0:
        raise ValueError("minimum energy must be positive")
    if int(table_points) < 128:
        raise ValueError("residual path table requires at least 128 points")

    output = np.zeros_like(energy_b, dtype=float)
    identity = depth_b == 0.0
    output[identity] = energy_b[identity]
    active = (
        (~identity)
        & (energy_b > minimum_energy_eV)
        & (cosine_b > 0.0))
    if not np.any(active):
        return float(output) if output.ndim == 0 else output

    incident_max = float(np.max(energy_b[active]))
    energy_grid, path_grid_nm = _csda_path_table(
        z1, m1, target, minimum_energy_eV, incident_max, table_points)
    incident_path_nm = np.interp(
        energy_b[active], energy_grid, path_grid_nm)
    traversed_path_nm = depth_b[active] / cosine_b[active]
    residual_path_nm = incident_path_nm - traversed_path_nm
    survives = residual_path_nm > 0.0
    remaining = np.zeros_like(residual_path_nm)
    remaining[survives] = np.interp(
        residual_path_nm[survives], path_grid_nm, energy_grid)
    # Interpolation roundoff must never create energy.
    output[active] = np.minimum(remaining, energy_b[active])
    return float(output) if output.ndim == 0 else output


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


_TARGETS = {
    "sio2": SIO2,
    "amorphous_carbon": AMORPHOUS_CARBON,
    "fluorocarbon_film": FLUOROCARBON_FILM,
}
_TABLE_CACHE = {}
_TABLE_ENERGY_POINTS = 72
_TABLE_COSINE_POINTS = 40
_TABLE_MAXIMUM_ENERGY_EV = 2.0e4


def _factor_table(layer_depth_nm, reference_energy_eV, target_name, z1, m1):
    """Dense (log-energy x cosine) table of the deposition factor, built once.

    The deposited-energy factor is smooth in both arguments, so bilinear
    interpolation on this grid is far below the quadrature receipts; building the
    table costs a few thousand scalar integrals exactly once per configuration.
    """
    key = (float(layer_depth_nm), float(reference_energy_eV), str(target_name),
           int(z1), float(m1))
    cached = _TABLE_CACHE.get(key)
    if cached is not None:
        return cached
    target = _TARGETS[str(target_name)]
    reference = nuclear_energy_in_layer_eV(
        reference_energy_eV, 1.0, layer_depth_nm, z1, m1, target)
    if reference <= 0.0:
        raise ValueError("reference energy deposits nothing in the layer")
    energies = np.geomspace(12.0, _TABLE_MAXIMUM_ENERGY_EV, _TABLE_ENERGY_POINTS)
    cosines = np.linspace(1e-3, 1.0, _TABLE_COSINE_POINTS)
    table = np.empty((_TABLE_ENERGY_POINTS, _TABLE_COSINE_POINTS))
    for i, energy in enumerate(energies):
        for j, mu in enumerate(cosines):
            table[i, j] = nuclear_energy_in_layer_eV(
                energy, mu, layer_depth_nm, z1, m1, target) / reference
    cached = (np.log(energies), cosines, table)
    _TABLE_CACHE[key] = cached
    return cached


def cached_layer_factor(energy_eV, cosine_incidence, layer_depth_nm,
                        reference_energy_eV, threshold_energy_eV, target_name,
                        z1, m1):
    """Vectorized deposition factor via bilinear interpolation of the dense table."""
    log_e_grid, mu_grid, table = _factor_table(
        layer_depth_nm, reference_energy_eV, target_name, z1, m1)
    energy = np.asarray(energy_eV, dtype=float)
    cosine = np.asarray(cosine_incidence, dtype=float)
    energy_b, cosine_b = np.broadcast_arrays(energy, cosine)
    log_e = np.log(np.clip(energy_b, 12.0, _TABLE_MAXIMUM_ENERGY_EV))
    mu = np.clip(cosine_b, mu_grid[0], 1.0)
    ei = np.clip(np.searchsorted(log_e_grid, log_e) - 1, 0, len(log_e_grid) - 2)
    mj = np.clip(np.searchsorted(mu_grid, mu) - 1, 0, len(mu_grid) - 2)
    te = (log_e - log_e_grid[ei]) / (log_e_grid[ei + 1] - log_e_grid[ei])
    tm = (mu - mu_grid[mj]) / (mu_grid[mj + 1] - mu_grid[mj])
    te = np.clip(te, 0.0, 1.0)
    tm = np.clip(tm, 0.0, 1.0)
    value = ((1 - te) * (1 - tm) * table[ei, mj]
             + te * (1 - tm) * table[ei + 1, mj]
             + (1 - te) * tm * table[ei, mj + 1]
             + te * tm * table[ei + 1, mj + 1])
    if threshold_energy_eV:
        value = np.where(energy_b <= threshold_energy_eV, 0.0, value)
    return value.reshape(energy_b.shape)
