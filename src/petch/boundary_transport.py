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


def adjoint_boundary_state_floor_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, floor_cells, *,
        proposal_species=None, n_face_position=8, max_steps=None, dt_cap=0.15, dt_field=0.10):
    """Generic Liouville adjoint floor gather using the boundary state's joint flux density.

    This function contains no species source law. The supplied species quadrature is used as the surface
    proposal, and its density model scores both surface and time-reversed plasma-exit states.
    """
    solid = np.asarray(solid, dtype=bool); nx, nz = solid.shape
    species = boundary.get(species_name)
    if species.density_model is None:
        raise ValueError("adjoint transport requires a boundary density model")
    proposal = species if proposal_species is None else proposal_species
    if proposal.density_model is None:
        raise ValueError("adjoint surface proposal requires a density model")
    if max_steps is None:
        max_steps = 200 * nz
    n_face_position = int(n_face_position)
    if n_face_position <= 0 or not floor_cells:
        raise ValueError("positive face quadrature and nonempty floor_cells are required")
    per_cell = np.zeros(len(floor_cells))
    base_velocity = proposal.velocity_sqrt_eV
    log_surface_density = proposal.log_flux_density(
        base_velocity, proposal.phase_rad, proposal.position_m)
    if not np.all(np.isfinite(log_surface_density)):
        raise ValueError("surface proposal samples must lie inside their density support")
    for ci, (cx, cz) in enumerate(floor_cells):
        face_u = (np.arange(n_face_position) + 0.5) / n_face_position
        x0 = np.tile(cx + face_u, base_velocity.shape[0])
        z0 = np.full(x0.shape, cz - 1e-3)
        # Time reverse the incident quadrature to launch outward from the floor.
        vx0 = np.repeat(-base_velocity[:, 0], n_face_position)
        vz0 = np.repeat(-base_velocity[:, 2], n_face_position)
        hit_x, _, _, _, survivor, exit_vx, exit_vz = trace_nodal(
            nodal_potential, solid, x0, z0, vx0, vz0, float(species.charge_number),
            nx, nz, int(max_steps), dt_cap, dt_field)
        escaped = (hit_x < 0) & (survivor < 0.5) & (exit_vz < 0.0)
        sample_index = np.repeat(np.arange(base_velocity.shape[0]), n_face_position)
        exit_forward = np.column_stack((
            -exit_vx,
            base_velocity[sample_index, 1],
            -exit_vz,
        ))
        log_exit_density = species.log_flux_density(exit_forward)
        surface_normal = base_velocity[sample_index, 2]
        exit_normal = np.maximum(-exit_vz, 1e-300)
        log_ratio = log_exit_density - log_surface_density[sample_index]
        with np.errstate(over="ignore", invalid="ignore"):
            density_ratio = np.exp(log_ratio)
        contribution = np.where(
            escaped & np.isfinite(log_ratio),
            surface_normal / exit_normal * density_ratio, 0.0)
        quadrature_weight = proposal.weight[sample_index] / n_face_position
        per_cell[ci] = float(np.sum(quadrature_weight * contribution))
    normalized_flux = float(per_cell.mean())
    return dict(normalized_flux=normalized_flux,
                absolute_flux_m2_s=normalized_flux * species.flux_m2_s,
                per_cell=per_cell)
