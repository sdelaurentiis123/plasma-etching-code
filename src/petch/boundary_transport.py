"""Transport adapters consuming the unified plasma boundary state."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from .charging_nodal import trace_nodal


@dataclass(frozen=True)
class BoundaryLaunches2D:
    x: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    normalized_weight: np.ndarray
    flux_weight_m2_s: np.ndarray


def boundary_launches_2d(species: SpeciesBoundaryState, x_min, x_max, n_position):
    """Tensor a weighted species phase-space measure with uniform lateral midpoint quadrature.

    The adapter contains no source physics. It preserves correlations already present in the boundary
    state. The ignored 3-D `vy` coordinate is returned for total-energy/surface-reaction accounting.
    """
    n_position = int(n_position)
    if n_position <= 0 or not x_max > x_min:
        raise ValueError("require n_position>0 and x_max>x_min")
    position = float(x_min) + (np.arange(n_position) + 0.5) * (float(x_max) - float(x_min)) / n_position
    sample_count = species.weight.size
    x = np.tile(position, sample_count)
    velocity = np.repeat(species.velocity_sqrt_eV, n_position, axis=0)
    normalized_weight = np.repeat(species.weight / n_position, n_position)
    flux_weight = normalized_weight * species.flux_m2_s
    for array in (x, velocity, normalized_weight, flux_weight):
        array.setflags(write=False)
    return BoundaryLaunches2D(
        x=x, vx=velocity[:, 0], vy=velocity[:, 1], vz=velocity[:, 2],
        normalized_weight=normalized_weight, flux_weight_m2_s=flux_weight,
    )


def trace_boundary_state_floor_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, target_mask, *,
        x_min=0.0, x_max=None, n_position=256, max_steps=None, dt_cap=0.15, dt_field=0.10):
    """Trace one boundary-state species and return normalized/absolute target flux.

    `target_mask` is a cell mask, making this adapter independent of named regions. Normalized flux is
    reported per unit target width relative to incident flux per unit source width.
    """
    solid = np.asarray(solid, dtype=bool); target_mask = np.asarray(target_mask, dtype=bool)
    if target_mask.shape != solid.shape or np.any(target_mask & ~solid):
        raise ValueError("target_mask must select solid cells")
    nx, nz = solid.shape
    if x_max is None:
        x_max = float(nx)
    species = boundary.get(species_name)
    launches = boundary_launches_2d(species, x_min, x_max, n_position)
    z = np.full(launches.x.shape, 1e-3)
    if max_steps is None:
        max_steps = 200 * nz
    hit_x, hit_z, *_ = trace_nodal(
        nodal_potential, solid, launches.x, z, launches.vx, launches.vz,
        float(species.charge_number), nx, nz, int(max_steps), dt_cap, dt_field)
    valid = hit_x >= 0
    target = np.zeros(hit_x.shape, dtype=bool)
    target[valid] = target_mask[hit_x[valid], hit_z[valid]]
    source_width = float(x_max) - float(x_min)
    target_width = float(np.count_nonzero(target_mask[:, np.max(np.where(target_mask)[1])]))
    if target_width <= 0.0:
        raise ValueError("target mask has zero horizontal measure")
    hit_probability = float(np.sum(launches.normalized_weight[target]))
    normalized_flux = hit_probability * source_width / target_width
    return dict(
        normalized_flux=normalized_flux,
        absolute_flux_m2_s=normalized_flux * species.flux_m2_s,
        hit_probability=hit_probability,
        incident_flux_m2_s=species.flux_m2_s,
    )
