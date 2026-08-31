"""Deterministic direct transport through an extruded mask footprint.

The operator in this module is deliberately narrower than the common 3-D
surface engine.  It evaluates the collision-free, no-wall-return component of
an incident angular distribution through a mask whose lateral footprint is
constant over its height.  That component is both useful and identifiable
before a material-specific wall-reaction law is available:

* energetic ions often remain close to the surface normal, so the direct map
  can dominate their entrance transmission;
* thermal radicals have a cosine incident distribution, so the direct map is
  a rigorous no-return contribution while diffuse/reactive wall return remains
  a separate declared operator.

All angular integration is deterministic Gaussian quadrature.  No random or
Monte-Carlo sampling enters the result.  The footprint is one unique periodic
cell (no duplicate endpoint); a nonperiodic calculation may instead select
solid exterior boundaries.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.ndimage import map_coordinates


@dataclass(frozen=True)
class AngularOrdinate:
    """One angular characteristic, expressed as transverse slopes."""

    tangent_x: float
    tangent_y: float
    weight: float

    def __post_init__(self):
        values = np.asarray(
            [self.tangent_x, self.tangent_y, self.weight], dtype=float)
        if np.any(~np.isfinite(values)) or self.weight < 0.0:
            raise ValueError("angular ordinate must be finite with nonnegative weight")
        object.__setattr__(self, "tangent_x", float(self.tangent_x))
        object.__setattr__(self, "tangent_y", float(self.tangent_y))
        object.__setattr__(self, "weight", float(self.weight))


@dataclass(frozen=True)
class ExtrudedMaskTransmission:
    """Direct-transmission field and the numerical contract that produced it."""

    transmission: np.ndarray
    mask_height: float
    grid_spacing: float
    ordinates: tuple[AngularOrdinate, ...]
    periodic_lateral: bool
    subdivisions_per_crossed_cell: float

    def __post_init__(self):
        field = np.asarray(self.transmission, dtype=float).copy()
        ordinates = tuple(self.ordinates)
        if (
            field.ndim != 2
            or field.size == 0
            or np.any(~np.isfinite(field))
            or np.any((field < -1e-12) | (field > 1.0 + 1e-12))
            or not np.isfinite(self.mask_height)
            or self.mask_height <= 0.0
            or not np.isfinite(self.grid_spacing)
            or self.grid_spacing <= 0.0
            or not ordinates
            or any(not isinstance(item, AngularOrdinate) for item in ordinates)
            or not np.isfinite(self.subdivisions_per_crossed_cell)
            or self.subdivisions_per_crossed_cell < 1.0
        ):
            raise ValueError("invalid extruded-mask transmission result")
        total_weight = sum(item.weight for item in ordinates)
        if not np.isclose(total_weight, 1.0, rtol=0.0, atol=2e-14):
            raise ValueError("angular quadrature weights must sum to one")
        field = np.clip(field, 0.0, 1.0)
        field.setflags(write=False)
        object.__setattr__(self, "transmission", field)
        object.__setattr__(self, "mask_height", float(self.mask_height))
        object.__setattr__(self, "grid_spacing", float(self.grid_spacing))
        object.__setattr__(self, "ordinates", ordinates)
        object.__setattr__(
            self, "subdivisions_per_crossed_cell",
            float(self.subdivisions_per_crossed_cell))


def gaussian_transverse_angle_ordinates(
        component_sigma_rad: float, *, order_per_component: int = 7):
    """Tensor Gauss-Hermite rule for two independent Gaussian angles.

    ``component_sigma_rad`` is the standard deviation of either signed
    transverse angular component.  Tangents, rather than small-angle values,
    are returned so propagation stays geometric away from zero angle.
    """
    if (
        not np.isfinite(component_sigma_rad)
        or component_sigma_rad < 0.0
        or int(order_per_component) != order_per_component
        or order_per_component <= 0
    ):
        raise ValueError("invalid Gaussian angular quadrature")
    nodes, weights = hermgauss(int(order_per_component))
    weights = weights / np.sqrt(np.pi)
    angles = np.sqrt(2.0) * float(component_sigma_rad) * nodes
    ordinates = tuple(
        AngularOrdinate(
            np.tan(theta_x), np.tan(theta_y), weights[ix] * weights[iy])
        for ix, theta_x in enumerate(angles)
        for iy, theta_y in enumerate(angles)
    )
    return ordinates


def cosine_flux_hemisphere_ordinates(
        *, polar_cosine_order: int = 8, azimuth_count: int = 16):
    """Product rule for an incident cosine-flux hemisphere.

    The flux density in polar cosine mu is ``2*mu`` on [0, 1].  Uniform
    midpoint azimuths avoid placing a preferred direction on rectilinear mask
    edges while preserving opposite-direction pairs.
    """
    if (
        int(polar_cosine_order) != polar_cosine_order
        or polar_cosine_order <= 0
        or int(azimuth_count) != azimuth_count
        or azimuth_count < 4
    ):
        raise ValueError("invalid cosine-flux angular quadrature")
    raw_mu, raw_weight = leggauss(int(polar_cosine_order))
    cosine = 0.5 * (raw_mu + 1.0)
    cosine_weight = 0.5 * raw_weight * 2.0 * cosine
    return tuple(
        AngularOrdinate(
            np.sqrt(1.0 - mu * mu) / mu * np.cos(phi),
            np.sqrt(1.0 - mu * mu) / mu * np.sin(phi),
            weight / int(azimuth_count),
        )
        for mu, weight in zip(cosine, cosine_weight)
        for phi in (
            2.0 * np.pi * (np.arange(int(azimuth_count)) + 0.5)
            / int(azimuth_count)
        )
    )


def direct_extruded_mask_transmission(
        opening, *, mask_height, grid_spacing, ordinates,
        periodic_lateral=True, subdivisions_per_crossed_cell=2.0):
    """Integrate direct characteristics through a constant mask footprint.

    Parameters use one common length unit.  ``opening`` is true in gas and
    false in mask over the unique lateral cell.  For each floor node and
    ordinate, the complete line from mask bottom to mask top is tested.  The
    path subdivision count is set by the largest number of grid cells crossed
    in x, y, or z, preventing one-cell tracks from being skipped.
    """
    gas = np.asarray(opening, dtype=bool)
    ordinates = tuple(ordinates)
    if (
        gas.ndim != 2
        or gas.size == 0
        or not np.any(gas)
        or not np.isfinite(mask_height)
        or mask_height <= 0.0
        or not np.isfinite(grid_spacing)
        or grid_spacing <= 0.0
        or not ordinates
        or any(not isinstance(item, AngularOrdinate) for item in ordinates)
        or not np.isfinite(subdivisions_per_crossed_cell)
        or subdivisions_per_crossed_cell < 1.0
    ):
        raise ValueError("invalid extruded-mask direct-transport inputs")
    total_weight = sum(item.weight for item in ordinates)
    if not np.isclose(total_weight, 1.0, rtol=0.0, atol=2e-14):
        raise ValueError("angular quadrature weights must sum to one")

    nx, ny = gas.shape
    x, y = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    occupancy = gas.astype(np.uint8)
    transmission = np.zeros(gas.shape, dtype=float)
    mode = "wrap" if periodic_lateral else "constant"
    for ordinate in ordinates:
        shift_x = ordinate.tangent_x * float(mask_height) / float(grid_spacing)
        shift_y = ordinate.tangent_y * float(mask_height) / float(grid_spacing)
        subdivisions = int(np.ceil(float(subdivisions_per_crossed_cell) * max(
            float(mask_height) / float(grid_spacing),
            abs(shift_x), abs(shift_y), 1.0)))
        alive = gas.copy()
        for index in range(1, subdivisions + 1):
            fraction = index / subdivisions
            alive &= map_coordinates(
                occupancy,
                [x + fraction * shift_x, y + fraction * shift_y],
                order=0,
                mode=mode,
                cval=0.0,
                prefilter=False,
            ).astype(bool)
            if not np.any(alive):
                break
        transmission += ordinate.weight * alive
    transmission[~gas] = 0.0
    return ExtrudedMaskTransmission(
        transmission,
        float(mask_height),
        float(grid_spacing),
        ordinates,
        bool(periodic_lateral),
        float(subdivisions_per_crossed_cell),
    )


__all__ = [
    "AngularOrdinate",
    "ExtrudedMaskTransmission",
    "cosine_flux_hemisphere_ordinates",
    "direct_extruded_mask_transmission",
    "gaussian_transverse_angle_ordinates",
]
