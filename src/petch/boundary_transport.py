"""Transport adapters consuming the unified plasma boundary state."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from .boundary_state import qmc_boundary_proposal
from .charging_nodal import trace_nodal
from .adaptive_quadrature import adaptive_surface_quadrature


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


def adjoint_boundary_state_face_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, cells, normals, *,
        proposal_species=None, n_face_position=8, max_steps=None, dt_cap=0.15, dt_field=0.10,
        want_energy=False, face_quadrature_offset=0.5, face_position_samples=None):
    """Generic Liouville adjoint gather on arbitrary axis-aligned material faces.

    This function contains no species source law. The supplied species quadrature is used as the surface
    proposal, and its density model scores both surface and time-reversed plasma-exit states. Velocities
    use the boundary convention (positive ``vz`` points into the feature). For outward material normal
    ``n``, a forward surface state is incident only when ``-v.n > 0``. Its time reverse is traced to the
    plasma plane and weighted by the Liouville flux Jacobian ``(-v.n)/vz_exit``.

    The returned values are fluxes per unit face length divided by incident flux per unit horizontal
    source length. No orientation correction is hidden in the caller.
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
    face_quadrature_offset = float(face_quadrature_offset)
    cells = [tuple(map(int, cell)) for cell in cells]
    normals = np.asarray(normals, dtype=float)
    if n_face_position <= 0 or not cells:
        raise ValueError("positive face quadrature and nonempty cells are required")
    if not 0.0 <= face_quadrature_offset < 1.0:
        raise ValueError("face_quadrature_offset must lie in [0,1)")
    if normals.shape != (len(cells), 2):
        raise ValueError("one 2-D outward normal is required per face")
    if (np.any(~np.isfinite(normals)) or
            np.any(np.abs(np.linalg.norm(normals, axis=1) - 1.0) > 1e-12)):
        raise ValueError("face normals must be finite unit vectors")
    per_face = np.zeros(len(cells)); energy_numerator = np.zeros(len(cells))
    base_velocity = proposal.velocity_sqrt_eV
    if face_position_samples is not None:
        face_position_samples = np.asarray(face_position_samples, dtype=float)
        if (face_position_samples.shape != (base_velocity.shape[0],)
                or np.any((face_position_samples < 0.0) | (face_position_samples >= 1.0))):
            raise ValueError("face_position_samples must provide one coordinate in [0,1) per velocity")
    log_surface_density = proposal.log_flux_density(
        base_velocity, proposal.phase_rad, proposal.position_m)
    if not np.all(np.isfinite(log_surface_density)):
        raise ValueError("surface proposal samples must lie inside their density support")
    for ci, ((cx, cz), (normal_x, normal_z)) in enumerate(zip(cells, normals)):
        if face_position_samples is None:
            face_u = (np.arange(n_face_position) + face_quadrature_offset) / n_face_position
            sample_index = np.repeat(np.arange(base_velocity.shape[0]), n_face_position)
            quadrature_weight = proposal.weight[sample_index] / n_face_position
        else:
            face_u = face_position_samples
            sample_index = np.arange(base_velocity.shape[0])
            quadrature_weight = proposal.weight
        # Face centre plus tangent coordinate, displaced one epsilon into the adjacent gas.
        face_s = face_u - 0.5
        x_center = cx + 0.5 + (0.5 + 1e-3) * normal_x
        z_center = cz + 0.5 + (0.5 + 1e-3) * normal_z
        if face_position_samples is None:
            x0 = np.tile(x_center - normal_z * face_s, base_velocity.shape[0])
            z0 = np.tile(z_center + normal_x * face_s, base_velocity.shape[0])
        else:
            x0 = x_center - normal_z * face_s
            z0 = z_center + normal_x * face_s
        # Time reverse the forward incident quadrature to launch outward from the material face.
        vx0 = -base_velocity[sample_index, 0]
        vz0 = -base_velocity[sample_index, 2]
        hit_x, _, _, _, survivor, exit_vx, exit_vz = trace_nodal(
            nodal_potential, solid, x0, z0, vx0, vz0, float(species.charge_number),
            nx, nz, int(max_steps), dt_cap, dt_field)
        escaped = (hit_x < 0) & (survivor < 0.5) & (exit_vz < 0.0)
        exit_forward = np.column_stack((
            -exit_vx,
            base_velocity[sample_index, 1],
            -exit_vz,
        ))
        # RF phase is a trajectory label for the electrostatic feature solve and is preserved under
        # time reversal. This retains phase-energy-angle correlations supplied by a sheath model.
        exit_phase = (None if proposal.phase_rad is None
                      else proposal.phase_rad[sample_index])
        log_exit_density = species.log_flux_density(exit_forward, exit_phase)
        surface_normal = -(base_velocity[sample_index, 0] * normal_x
                           + base_velocity[sample_index, 2] * normal_z)
        exit_normal = np.maximum(-exit_vz, 1e-300)
        log_ratio = log_exit_density - log_surface_density[sample_index]
        with np.errstate(over="ignore", invalid="ignore"):
            density_ratio = np.exp(log_ratio)
        contribution = np.where(
            escaped & (surface_normal > 0.0) & np.isfinite(log_ratio),
            surface_normal / exit_normal * density_ratio, 0.0)
        per_face[ci] = float(np.sum(quadrature_weight * contribution))
        if want_energy:
            # The adjoint launches at the material face, so its initial kinetic energy is exactly the
            # time-reversed forward impact energy. The tracer's ``impact_energy`` is zero for the desired
            # escaping adjoint paths because those paths do not hit material a second time.
            impact_energy_3d = np.sum(base_velocity[sample_index] ** 2, axis=1)
            energy_numerator[ci] = float(np.sum(
                quadrature_weight * contribution * impact_energy_3d))
    normalized_flux = float(per_face.mean())
    result = dict(normalized_flux=normalized_flux,
                  absolute_flux_m2_s=normalized_flux * species.flux_m2_s,
                  per_face=per_face, per_cell=per_face)
    if want_energy:
        result["mean_impact_energy_eV_per_face"] = np.divide(
            energy_numerator, per_face, out=np.zeros_like(per_face), where=per_face > 0.0)
        result["mean_impact_energy_eV"] = (
            float(energy_numerator.sum() / per_face.sum()) if per_face.sum() > 0.0 else 0.0)
    return result


def adaptive_adjoint_boundary_state_face_flux(
        boundary, species_name, nodal_potential, solid, cells, normals, *, proposal_species=None,
        n_face_position=4, base_log2=6, max_log2=12, n_replicates=4, seed=0,
        absolute_tolerance=1e-3, relative_tolerance=5e-3,
        element_absolute_tolerance=None, element_relative_tolerance=0.0, refine_fraction=0.5,
        max_steps=None, dt_cap=0.15, dt_field=0.10):
    """Universally adapt randomized phase-space quadrature on arbitrary material faces."""
    physical = boundary.get(species_name)
    template = physical if proposal_species is None else proposal_species
    cells = [tuple(cell) for cell in cells]; normals = np.asarray(normals, dtype=float)

    def evaluator(indices, log2_samples, replicate_seed):
        proposal = qmc_boundary_proposal(
            template, log2_samples, replicate_seed, name=f"{species_name}-adaptive-proposal")
        # Jointly refine surface position with velocity. Randomly permuted stratification avoids a
        # fixed velocity-position correlation and works for stratified mixtures of any component count.
        rng = np.random.default_rng(replicate_seed + 7919)
        count = proposal.weight.size
        face_position = (rng.permutation(count) + rng.random(count)) / count
        result = adjoint_boundary_state_face_flux(
            boundary, species_name, nodal_potential, solid,
            [cells[index] for index in indices], normals[indices], proposal_species=proposal,
            n_face_position=n_face_position, max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field,
            face_position_samples=face_position)
        return result["per_face"]

    return adaptive_surface_quadrature(
        evaluator, len(cells), weights=np.full(len(cells), 1.0 / len(cells)),
        base_log2=base_log2, max_log2=max_log2, n_replicates=n_replicates, seed=seed,
        absolute_tolerance=absolute_tolerance, relative_tolerance=relative_tolerance,
        element_absolute_tolerance=element_absolute_tolerance,
        element_relative_tolerance=element_relative_tolerance, refine_fraction=refine_fraction)


def adjoint_boundary_state_floor_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, floor_cells, *,
        proposal_species=None, n_face_position=8, max_steps=None, dt_cap=0.15, dt_field=0.10,
        want_energy=False):
    """Backward-compatible horizontal-floor specialization of the arbitrary-face gather."""
    return adjoint_boundary_state_face_flux(
        boundary, species_name, nodal_potential, solid, floor_cells,
        [(0.0, -1.0)] * len(floor_cells), proposal_species=proposal_species,
        n_face_position=n_face_position, max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field,
        want_energy=want_energy)
