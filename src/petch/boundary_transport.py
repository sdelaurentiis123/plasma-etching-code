"""Transport adapters consuming the unified plasma boundary state."""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from .boundary_state import qmc_boundary_proposal, qmc_boundary_proposal_with_auxiliary
from .charging_nodal import trace_nodal
from .adaptive_quadrature import AdaptiveQuadratureResult, adaptive_surface_quadrature


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
        x_min=0.0, x_max=None, n_position=256, max_steps=None, dt_cap=0.15, dt_field=0.10,
        fixed_dt=0.0):
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
        float(species.charge_number), nx, nz, int(max_steps), dt_cap, dt_field, fixed_dt)
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


def forward_boundary_state_cell_flux_qmc(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, cells, *,
        normals=None, proposal_species=None, log2_samples=12, seed=0, x_min=0.0, x_max=None,
        max_steps=None, dt_cap=0.15, dt_field=0.10, fixed_dt=0.0):
    """Forward QMC current to arbitrary material cells from the same boundary-density contract.

    The score is total current to each unit-depth cell divided by incident current per unit horizontal
    source length. It is therefore directly comparable to the sum of adjoint unit-face scores belonging
    to that cell. No target orientation or feature label enters the estimator.
    """
    solid = np.asarray(solid, dtype=bool); nx, nz = solid.shape
    if x_max is None:
        x_max = float(nx)
    if not float(x_max) > float(x_min):
        raise ValueError("x_max must exceed x_min")
    species = boundary.get(species_name)
    template = species if proposal_species is None else proposal_species
    proposal = qmc_boundary_proposal(
        template, int(log2_samples), int(seed), name=f"{species_name}-forward-proposal")
    velocity = proposal.velocity_sqrt_eV; count = velocity.shape[0]
    rng = np.random.default_rng(int(seed) + 15485863)
    lateral_u = (rng.permutation(count) + rng.random(count)) / count
    x0 = float(x_min) + (float(x_max) - float(x_min)) * lateral_u
    z0 = np.full(count, 1e-3)
    traced = trace_nodal(
        nodal_potential, solid, x0, z0, velocity[:, 0], velocity[:, 2],
        float(species.charge_number), nx, nz,
        int(200 * nz if max_steps is None else max_steps), dt_cap, dt_field, fixed_dt)
    hit_x, hit_z, hit_nx, hit_nz = traced[0], traced[1], traced[7], traced[8]
    hit_x_position, hit_z_position = traced[9], traced[10]
    log_physical = species.log_flux_density(
        velocity, proposal.phase_rad, proposal.position_m)
    log_proposal = proposal.log_flux_density(
        velocity, proposal.phase_rad, proposal.position_m)
    with np.errstate(over="ignore", invalid="ignore"):
        density_ratio = np.exp(log_physical - log_proposal)
    score = proposal.weight * density_ratio * (float(x_max) - float(x_min))
    cells = [tuple(map(int, cell)) for cell in cells]
    if normals is not None:
        normals = np.asarray(normals, dtype=int)
        if normals.shape != (len(cells), 2):
            raise ValueError("normals must match cells")
    result = np.zeros(len(cells)); hit_count = np.zeros(len(cells), dtype=np.int64); lookup = {}
    endpoint_result = None if normals is None else np.zeros((len(cells), 2))
    for index, cell in enumerate(cells):
        key = cell if normals is None else (cell, tuple(normals[index]))
        lookup.setdefault(key, []).append(index)
    valid = hit_x >= 0
    for hx, hz, hnx, hnz, xpos, zpos, value in zip(
            hit_x[valid], hit_z[valid], hit_nx[valid], hit_nz[valid],
            hit_x_position[valid], hit_z_position[valid], score[valid]):
        cell = (int(hx), int(hz))
        key = cell if normals is None else (cell, (int(hnx), int(hnz)))
        indices = lookup.get(key)
        if indices is not None:
            # A cell current is stored once even when the cell owns multiple exposed faces.
            result[indices[0]] += value
            hit_count[indices[0]] += 1
            if endpoint_result is not None:
                face_index = indices[0]
                normal = tuple(normals[face_index])
                u = ((float(zpos) - cells[face_index][1]) if normal[0] != 0
                     else (float(xpos) - cells[face_index][0]))
                u = min(max(u, 0.0), 1.0)
                endpoint_result[face_index, 0] += value * (1.0 - u)
                endpoint_result[face_index, 1] += value * u
    return dict(
        per_cell=result, normalized_total=float(result.sum()),
        absolute_per_cell_m2_s=result * species.flux_m2_s,
        per_face_endpoint=endpoint_result,
        hit_count_per_cell=hit_count, samples=count,
        equal_weight_sampling=bool(np.allclose(
            proposal.weight, np.full(count, 1.0 / count), rtol=1e-12, atol=1e-15)),
        seed=int(seed))


def adjoint_boundary_state_face_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, cells, normals, *,
        proposal_species=None, n_face_position=8, max_steps=None, dt_cap=0.15, dt_field=0.10,
        want_energy=False, face_quadrature_offset=0.5, face_position_samples=None,
        fixed_dt=0.0, face_offset=1e-3):
    """Generic Liouville adjoint gather on arbitrary axis-aligned material faces.

    This function contains no species source law. The supplied species quadrature is used as the surface
    proposal is interpreted in each face's local tangent/outward-normal frame: proposal ``vx`` is the
    in-plane tangent and positive proposal ``vz`` is the inward normal speed. Rotating that local rule
    per face is essential for support completeness: on a vertical wall, global ``vz`` is tangential and
    must span both signs. The time-reversed global state is traced to the plasma plane and weighted by
    the Liouville flux Jacobian ``v_normal_in/vz_exit``. The proposal density is numerical and remains
    evaluated in its local coordinates; the rotation has unit Jacobian.

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
    face_offset = float(face_offset)
    cells = [tuple(map(int, cell)) for cell in cells]
    normals = np.asarray(normals, dtype=float)
    if n_face_position <= 0 or not cells:
        raise ValueError("positive face quadrature and nonempty cells are required")
    if not 0.0 <= face_quadrature_offset < 1.0:
        raise ValueError("face_quadrature_offset must lie in [0,1)")
    if not np.isfinite(face_offset) or not 0.0 < face_offset < 0.5:
        raise ValueError("face_offset must lie strictly between zero and half a cell")
    if normals.shape != (len(cells), 2):
        raise ValueError("one 2-D outward normal is required per face")
    if (np.any(~np.isfinite(normals)) or
            np.any(np.abs(np.linalg.norm(normals, axis=1) - 1.0) > 1e-12)):
        raise ValueError("face normals must be finite unit vectors")
    per_face = np.zeros(len(cells)); energy_numerator = np.zeros(len(cells))
    effective_sample_size = np.zeros(len(cells)); max_sample_fraction = np.zeros(len(cells))
    dominant_sample_index = np.full(len(cells), -1, dtype=int)
    dominant_surface_velocity = np.zeros((len(cells), 3))
    dominant_exit_velocity = np.zeros((len(cells), 3))
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
    if face_position_samples is None:
        face_u = (np.arange(n_face_position) + face_quadrature_offset) / n_face_position
        sample_index = np.repeat(np.arange(base_velocity.shape[0]), n_face_position)
        quadrature_weight = proposal.weight[sample_index] / n_face_position
        sample_face_u = np.tile(face_u, base_velocity.shape[0])
        face_s = sample_face_u - 0.5
    else:
        face_u = face_position_samples
        sample_index = np.arange(base_velocity.shape[0])
        quadrature_weight = proposal.weight
        sample_face_u = face_u
        face_s = sample_face_u - 0.5
    cell_array = np.asarray(cells, dtype=float)
    x_center = cell_array[:, 0] + 0.5 + (0.5 + face_offset) * normals[:, 0]
    z_center = cell_array[:, 1] + 0.5 + (0.5 + face_offset) * normals[:, 1]
    x0 = (x_center[:, None] - normals[:, 1, None] * face_s[None, :]).ravel()
    z0 = (z_center[:, None] + normals[:, 0, None] * face_s[None, :]).ravel()
    samples_per_face = sample_index.size; face_count = len(cells)
    local_tangent = base_velocity[sample_index, 0]
    local_normal_in = base_velocity[sample_index, 2]
    # In-plane unit tangent t=(-nz,nx). A forward incident surface velocity is
    # v = local_tangent*t - local_normal_in*n; trace its exact time reverse.
    surface_vx = (-normals[:, 1, None] * local_tangent[None, :]
                  - normals[:, 0, None] * local_normal_in[None, :])
    surface_vz = (normals[:, 0, None] * local_tangent[None, :]
                  - normals[:, 1, None] * local_normal_in[None, :])
    vx0 = -surface_vx.ravel()
    vz0 = -surface_vz.ravel()
    traced = trace_nodal(
        nodal_potential, solid, x0, z0, vx0, vz0, float(species.charge_number),
        nx, nz, int(max_steps), dt_cap, dt_field, fixed_dt)
    hit_x, survivor, exit_vx, exit_vz = traced[0], traced[4], traced[5], traced[6]
    escaped = (hit_x < 0) & (survivor < 0.5) & (exit_vz < 0.0)
    tiled_sample_index = np.tile(sample_index, face_count)
    exit_forward = np.column_stack((
        -exit_vx, base_velocity[tiled_sample_index, 1], -exit_vz))
    exit_phase = (None if proposal.phase_rad is None
                  else proposal.phase_rad[tiled_sample_index])
    log_exit_density = species.log_flux_density(exit_forward, exit_phase)
    surface_normal = np.tile(local_normal_in, face_count)
    exit_normal = np.maximum(-exit_vz, 1e-300)
    log_ratio = log_exit_density - np.tile(log_surface_density[sample_index], face_count)
    with np.errstate(over="ignore", invalid="ignore"):
        density_ratio = np.exp(log_ratio)
    contribution = np.where(
        escaped & (surface_normal > 0.0) & np.isfinite(log_ratio),
        surface_normal / exit_normal * density_ratio, 0.0).reshape(face_count, samples_per_face)
    weighted_contribution = contribution * quadrature_weight[None, :]
    per_face[:] = weighted_contribution.sum(axis=1)
    per_face_endpoint = np.stack((
        np.sum(weighted_contribution * (1.0 - sample_face_u)[None, :], axis=1),
        np.sum(weighted_contribution * sample_face_u[None, :], axis=1)), axis=1)
    # Canonical endpoint order is top-to-bottom on vertical faces and left-to-right on horizontal
    # faces, matching ``material_face_nodes`` and exact forward hit deposition. The local tangent
    # parameter runs oppositely on left- and bottom-facing material edges.
    reverse_endpoint = (normals[:, 0] < 0.0) | (normals[:, 1] > 0.0)
    per_face_endpoint[reverse_endpoint] = per_face_endpoint[reverse_endpoint, ::-1]
    square_sum = np.sum(weighted_contribution * weighted_contribution, axis=1)
    positive = square_sum > 0.0
    effective_sample_size[positive] = per_face[positive] ** 2 / square_sum[positive]
    dominant = np.argmax(weighted_contribution, axis=1)
    dominant_sample_index[positive] = sample_index[dominant[positive]]
    max_sample_fraction[positive] = (
        weighted_contribution[np.where(positive)[0], dominant[positive]] / per_face[positive])
    positive_face = np.where(positive)[0]
    positive_sample = dominant[positive]
    dominant_surface_velocity[positive, 0] = surface_vx[positive_face, positive_sample]
    dominant_surface_velocity[positive, 1] = base_velocity[dominant_sample_index[positive], 1]
    dominant_surface_velocity[positive, 2] = surface_vz[positive_face, positive_sample]
    exit_matrix = exit_forward.reshape(face_count, samples_per_face, 3)
    dominant_exit_velocity[positive] = exit_matrix[np.where(positive)[0], dominant[positive]]
    if want_energy:
        impact_energy = np.sum(base_velocity[sample_index] ** 2, axis=1)
        energy_numerator[:] = np.sum(weighted_contribution * impact_energy[None, :], axis=1)
    normalized_flux = float(per_face.mean())
    result = dict(normalized_flux=normalized_flux,
                  absolute_flux_m2_s=normalized_flux * species.flux_m2_s,
                  per_face=per_face, per_cell=per_face,
                  per_face_endpoint=per_face_endpoint,
                  effective_sample_size=effective_sample_size,
                  max_sample_fraction=max_sample_fraction,
                  dominant_sample_index=dominant_sample_index,
                  dominant_surface_velocity=dominant_surface_velocity,
                  dominant_exit_velocity=dominant_exit_velocity)
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
        max_steps=None, dt_cap=0.15, dt_field=0.10, fixed_dt=0.0,
        initial_log2_samples=None, face_offset=1e-3):
    """Universally adapt randomized phase-space quadrature on arbitrary material faces."""
    physical = boundary.get(species_name)
    template = physical if proposal_species is None else proposal_species
    cells = [tuple(cell) for cell in cells]; normals = np.asarray(normals, dtype=float)
    endpoint_cache = {}

    def evaluate_full(indices, log2_samples, replicate_seed):
        proposal, auxiliary = qmc_boundary_proposal_with_auxiliary(
            template, log2_samples, 1, replicate_seed,
            name=f"{species_name}-adaptive-proposal")
        face_position = auxiliary[:, 0]
        result = adjoint_boundary_state_face_flux(
            boundary, species_name, nodal_potential, solid,
            [cells[index] for index in indices], normals[indices], proposal_species=proposal,
            n_face_position=n_face_position, max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field,
            face_position_samples=face_position, fixed_dt=fixed_dt, face_offset=face_offset)
        for local_index, element_index in enumerate(indices):
            endpoint_cache[(int(replicate_seed), int(log2_samples), int(element_index))] = (
                result["per_face_endpoint"][local_index].copy())
        return result

    def evaluator(indices, log2_samples, replicate_seed):
        return evaluate_full(indices, log2_samples, replicate_seed)["per_face"]

    adaptive = adaptive_surface_quadrature(
        evaluator, len(cells), weights=np.full(len(cells), 1.0 / len(cells)),
        base_log2=base_log2, max_log2=max_log2, n_replicates=n_replicates, seed=seed,
        absolute_tolerance=absolute_tolerance, relative_tolerance=relative_tolerance,
        element_absolute_tolerance=element_absolute_tolerance,
        element_relative_tolerance=element_relative_tolerance, refine_fraction=refine_fraction,
        initial_log2_samples=initial_log2_samples)
    endpoint_replicates = np.empty((n_replicates, len(cells), 2))
    for element_index, level in enumerate(adaptive.log2_samples):
        for replicate in range(n_replicates):
            endpoint_replicates[replicate, element_index] = endpoint_cache[
                (int(seed + replicate), int(level), int(element_index))]
    return replace(
        adaptive, auxiliary_mean=endpoint_replicates.mean(axis=0),
        auxiliary_replicates=endpoint_replicates)


def adaptive_forward_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, cells, *, proposal_species=None,
        normals=None,
        base_log2=8, max_log2=16, n_replicates=4, seed=0,
        absolute_tolerance=1e-3, relative_tolerance=5e-3,
        element_absolute_tolerance=None, element_relative_tolerance=0.0,
        refine_fraction=0.5, initial_log2_samples=None,
        x_min=0.0, x_max=None, max_steps=None, dt_cap=0.15, dt_field=0.10,
        fixed_dt=0.0):
    """Adapt the complementary forward estimator on arbitrary physical material cells."""
    cells = [tuple(map(int, cell)) for cell in cells]
    if proposal_species is not None:
        raise ValueError("adaptive forward zero-hit bounds currently require direct physical sampling")
    if n_replicates < 2:
        raise ValueError("at least two forward replicates are required")
    nx = np.asarray(solid).shape[0]
    source_width = float(nx if x_max is None else x_max) - float(x_min)
    evaluations = 0; converged = False; rounds = 0
    estimates = np.empty((n_replicates, len(cells)))
    endpoint_estimates = (None if normals is None
                          else np.empty((n_replicates, len(cells), 2)))
    hit_counts = np.empty((n_replicates, len(cells)), dtype=np.int64)
    replicate_samples = np.empty(n_replicates, dtype=np.int64)
    equal_weight_sampling = True
    start_level = int(base_log2)
    if initial_log2_samples is not None:
        initial = np.asarray(initial_log2_samples, dtype=int)
        if initial.shape != (len(cells),):
            raise ValueError("initial forward levels have wrong shape")
        start_level = int(initial.max())
    for level in range(start_level, int(max_log2) + 1):
        rounds += 1
        for replicate in range(n_replicates):
            result = forward_boundary_state_cell_flux_qmc(
                boundary, species_name, nodal_potential, solid, cells,
                normals=normals,
                log2_samples=level, seed=int(seed + replicate), x_min=x_min, x_max=x_max,
                max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field, fixed_dt=fixed_dt)
            estimates[replicate] = result["per_cell"]
            if endpoint_estimates is not None:
                endpoint_estimates[replicate] = result["per_face_endpoint"]
            hit_counts[replicate] = result["hit_count_per_cell"]
            replicate_samples[replicate] = result["samples"]
            equal_weight_sampling &= result["equal_weight_sampling"]
        evaluations += n_replicates * (2 ** level)
        element_mean = estimates.mean(axis=0)
        element_stderr = estimates.std(axis=0, ddof=1) / np.sqrt(n_replicates)
        if equal_weight_sampling:
            # A handful of randomized-QMC replicates can accidentally agree on the same rare hit
            # count and dramatically understate uncertainty. Direct physical sampling has equal
            # weights, so the ordinary Bernoulli hit-count standard error is a conservative floor.
            # This does not change the mean and does not assume a region or species.
            sample_count = int(replicate_samples.sum())
            hit_probability = hit_counts.sum(axis=0) / sample_count
            count_stderr = source_width * np.sqrt(
                hit_probability * (1.0 - hit_probability) / sample_count)
            element_stderr = np.maximum(element_stderr, count_stderr)
        # Zero observed hits do not imply zero flux. For direct physical sampling, the rule-of-three
        # 95% upper probability bound is 3/N across the pooled independent histories.
        zero = np.all(estimates == 0.0, axis=0)
        zero_sample_count = (int(replicate_samples.sum()) if equal_weight_sampling
                             else n_replicates * (2 ** level))
        zero_upper = 3.0 * source_width / zero_sample_count
        element_stderr[zero] = np.maximum(element_stderr[zero], zero_upper)
        totals = estimates.mean(axis=1)
        total_mean = float(totals.mean())
        total_stderr = float(totals.std(ddof=1) / np.sqrt(n_replicates))
        total_ok = total_stderr <= absolute_tolerance + relative_tolerance * abs(total_mean)
        element_ok = (element_absolute_tolerance is None or np.all(
            element_stderr <= element_absolute_tolerance
            + element_relative_tolerance * np.abs(element_mean)))
        if total_ok and element_ok:
            converged = True
            break
    levels = np.full(len(cells), level, dtype=np.int64)
    return AdaptiveQuadratureResult(
        element_mean=element_mean.copy(), element_stderr=element_stderr.copy(),
        element_replicates=estimates.copy(), log2_samples=levels,
        total_mean=total_mean, total_stderr=total_stderr, converged=converged,
        rounds=rounds, evaluations=evaluations,
        auxiliary_mean=(None if endpoint_estimates is None
                        else endpoint_estimates.mean(axis=0)),
        auxiliary_replicates=endpoint_estimates)


def bidirectional_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, cells, normals, *,
        proposal_species=None, adjoint_options=None, forward_options=None,
        element_absolute_tolerance=1e-3, element_relative_tolerance=0.05,
        method_hint=None, switch_factor=2.0, consistency_sigma=5.0,
        support_sigma=2.0, support_ratio=0.5):
    """Select forward or adjoint current per physical cell solely by measured uncertainty.

    This is not a named-region switch. Both unbiased estimators use the same physical boundary density;
    the lower-uncertainty estimate wins independently for each cell. The result refuses cells for which
    neither direction meets the requested mixed tolerance. When both estimators claim certification,
    they must also agree within ``consistency_sigma`` combined standard errors; otherwise a precise but
    biased transport path could silently control the charging fixed point.
    """
    cells = [tuple(map(int, cell)) for cell in cells]
    normals = np.asarray(normals, dtype=float)
    adjoint_kwargs = {} if adjoint_options is None else dict(adjoint_options)
    forward_kwargs = {} if forward_options is None else dict(forward_options)
    adjoint_fixed_dt = float(adjoint_kwargs.get("fixed_dt", 0.0))
    forward_fixed_dt = float(forward_kwargs.get("fixed_dt", 0.0))
    if np.ptp(np.asarray(nodal_potential, dtype=float)) > 0.0:
        if adjoint_fixed_dt <= 0.0 or forward_fixed_dt <= 0.0:
            raise ValueError(
                "bidirectional transport in a nonuniform field requires a positive fixed_dt; "
                "state-dependent timesteps do not define a reversible adjoint map")
        if not np.isclose(adjoint_fixed_dt, forward_fixed_dt, rtol=0.0, atol=0.0):
            raise ValueError("forward and adjoint estimators require the same fixed_dt")
    adjoint_kwargs.setdefault("element_absolute_tolerance", element_absolute_tolerance)
    adjoint_kwargs.setdefault("element_relative_tolerance", element_relative_tolerance)
    forward_kwargs.setdefault("element_absolute_tolerance", element_absolute_tolerance)
    forward_kwargs.setdefault("element_relative_tolerance", element_relative_tolerance)
    adjoint = adaptive_adjoint_boundary_state_face_flux(
        boundary, species_name, nodal_potential, solid, cells, normals,
        proposal_species=proposal_species, **adjoint_kwargs)

    unique_cells = list(dict.fromkeys(cells))
    if method_hint is not None:
        method_hint = np.asarray(method_hint)
        if method_hint.shape != (len(unique_cells),):
            raise ValueError("method_hint must match the number of unique cells")
    if switch_factor < 1.0:
        raise ValueError("switch_factor must be at least one")
    if consistency_sigma <= 0.0:
        raise ValueError("consistency_sigma must be positive")
    if support_sigma <= 0.0:
        raise ValueError("support_sigma must be positive")
    if not 0.0 < support_ratio < 1.0:
        raise ValueError("support_ratio must lie in (0, 1)")
    forward = adaptive_forward_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, cells, normals=normals,
        proposal_species=None, **forward_kwargs)
    if (adjoint.element_replicates.shape[0] < 4
            or forward.element_replicates.shape[0] < 4):
        raise ValueError(
            "bidirectional certification requires at least four independent replicates")
    if adjoint.element_replicates.shape[0] != forward.element_replicates.shape[0]:
        raise ValueError("bidirectional estimators require the same replicate count")
    first_face = {cell: cells.index(cell) for cell in unique_cells}

    # The forward controller resolves oriented faces, but a dielectric cell owns the sum of all its
    # exposed-face currents. Individually certified faces can still have an uncertified pooled cell
    # current. Raise the common forward sample level until every physical-cell sum is resolved or the
    # declared maximum is reached.
    forward_pool_refinement_rounds = 0
    max_forward_level = int(forward_kwargs.get("max_log2", 16))
    while True:
        pooled_forward_ok = True
        for cell in unique_cells:
            face_indices = np.array([i for i, item in enumerate(cells) if item == cell], dtype=int)
            replicates = forward.element_replicates[:, face_indices].sum(axis=1)
            mean = float(replicates.mean())
            stderr = max(
                float(replicates.std(ddof=1) / np.sqrt(replicates.size)),
                float(np.sqrt(np.sum(forward.element_stderr[face_indices] ** 2))))
            pooled_forward_ok &= stderr <= (
                element_absolute_tolerance + element_relative_tolerance * abs(mean))
        current_forward_level = int(forward.log2_samples.max())
        if pooled_forward_ok or current_forward_level >= max_forward_level:
            break
        next_level = current_forward_level + 1
        refined_forward_options = dict(forward_kwargs)
        refined_forward_options["initial_log2_samples"] = np.full(
            len(cells), next_level, dtype=int)
        forward = adaptive_forward_boundary_state_cell_flux(
            boundary, species_name, nodal_potential, solid, cells, normals=normals,
            proposal_species=None, **refined_forward_options)
        forward_pool_refinement_rounds += 1

    # A replicate-only importance-sampling controller cannot detect a mode missed by every scramble:
    # all replicates then agree on a spuriously precise near-zero value. Direct forward transport is
    # the independent support audit. When two nominally precise estimates disagree, refine both the
    # selected adjoint faces and the common forward history budget without using species, orientation,
    # region, aspect ratio, or geometry names.
    cross_refinement_rounds = 0
    forward_cross_refinement_rounds = 0
    max_adjoint_level = int(adjoint_kwargs.get("max_log2", 12))
    while True:
        inconsistent_cells = []
        for cell_index, cell in enumerate(unique_cells):
            face_indices = np.array([i for i, item in enumerate(cells) if item == cell], dtype=int)
            replicates = adjoint.element_replicates[:, face_indices].sum(axis=1)
            adjoint_mean = float(replicates.mean())
            adjoint_stderr = float(replicates.std(ddof=1) / np.sqrt(replicates.size))
            forward_replicates = forward.element_replicates[:, face_indices].sum(axis=1)
            forward_mean = float(forward_replicates.mean())
            forward_stderr = max(
                float(forward_replicates.std(ddof=1) / np.sqrt(forward_replicates.size)),
                float(np.sqrt(np.sum(forward.element_stderr[face_indices] ** 2))))
            forward_ok = forward_stderr <= (
                element_absolute_tolerance + element_relative_tolerance * abs(forward_mean))
            adjoint_ok = adjoint_stderr <= (
                element_absolute_tolerance + element_relative_tolerance * abs(adjoint_mean))
            combined = np.hypot(forward_stderr, adjoint_stderr)
            discrepancy = (abs(forward_mean - adjoint_mean) / combined
                           if combined > 0.0
                           else (0.0 if forward_mean == adjoint_mean else np.inf))
            support_separated = (
                adjoint_mean + support_sigma * adjoint_stderr
                < support_ratio * max(
                    forward_mean - support_sigma * forward_stderr, 0.0))
            if (forward_ok and adjoint_ok
                    and (discrepancy > consistency_sigma or support_separated)):
                inconsistent_cells.append(cell)
        if not inconsistent_cells:
            break
        inconsistent_set = set(inconsistent_cells)
        levels = adjoint.log2_samples.copy()
        selected_faces = np.array([
            index for index, cell in enumerate(cells) if cell in inconsistent_set
            and levels[index] < max_adjoint_level], dtype=int)
        current_forward_level = int(forward.log2_samples.max())
        refine_forward = current_forward_level < max_forward_level
        if selected_faces.size == 0 and not refine_forward:
            break
        if selected_faces.size:
            levels[selected_faces] = np.minimum(
                levels[selected_faces] + 2, max_adjoint_level)
            refined_options = dict(adjoint_kwargs)
            refined_options["initial_log2_samples"] = levels
            adjoint = adaptive_adjoint_boundary_state_face_flux(
                boundary, species_name, nodal_potential, solid, cells, normals,
                proposal_species=proposal_species, **refined_options)
        if refine_forward:
            refined_forward_options = dict(forward_kwargs)
            refined_forward_options["initial_log2_samples"] = np.full(
                len(cells), current_forward_level + 1, dtype=int)
            refined_forward = adaptive_forward_boundary_state_cell_flux(
                boundary, species_name, nodal_potential, solid, cells, normals=normals,
                proposal_species=None, **refined_forward_options)
            if int(refined_forward.log2_samples.max()) <= current_forward_level:
                max_forward_level = current_forward_level
            else:
                forward = refined_forward
                forward_cross_refinement_rounds += 1
        cross_refinement_rounds += 1

    # Refinement replaces the replicate ensembles. Endpoint moments must be read from the final
    # ensembles, never from the pre-refinement cache, or nodal deposition becomes history-dependent
    # even while face totals are certified.
    adjoint_endpoint_replicates = (
        adjoint.auxiliary_replicates if adjoint.auxiliary_replicates is not None
        else np.repeat(0.5 * adjoint.element_replicates[:, :, None], 2, axis=2))
    forward_endpoint_replicates = (
        forward.auxiliary_replicates if forward.auxiliary_replicates is not None
        else np.repeat(0.5 * forward.element_replicates[:, :, None], 2, axis=2))
    per_face = np.zeros(len(cells)); per_face_stderr = np.zeros(len(cells))
    selected_face_mean = np.zeros(len(cells)); selected_face_stderr = np.zeros(len(cells))
    selected_face_replicates = np.zeros_like(adjoint.element_replicates)
    selected_endpoint_replicates = np.zeros_like(adjoint_endpoint_replicates)
    method = np.empty(len(unique_cells), dtype="U7")
    forward_cell_mean = np.zeros(len(unique_cells)); forward_cell_stderr = np.zeros(len(unique_cells))
    adjoint_cell_mean = np.zeros(len(unique_cells)); adjoint_cell_stderr = np.zeros(len(unique_cells))
    estimator_discrepancy_sigma = np.zeros(len(unique_cells))
    method_within_tolerance = np.zeros(len(unique_cells), dtype=bool)
    estimator_consistent = np.ones(len(unique_cells), dtype=bool)
    cell_converged = np.zeros(len(unique_cells), dtype=bool)
    adjoint_zero_unresolved = np.zeros(len(unique_cells), dtype=bool)
    adjoint_support_unresolved = np.zeros(len(unique_cells), dtype=bool)
    converged = True
    for cell_index, cell in enumerate(unique_cells):
        face_indices = np.array([i for i, item in enumerate(cells) if item == cell], dtype=int)
        adjoint_replicates = adjoint.element_replicates[:, face_indices].sum(axis=1)
        adjoint_mean = float(adjoint_replicates.mean())
        adjoint_stderr = float(adjoint_replicates.std(ddof=1) / np.sqrt(adjoint_replicates.size))
        forward_replicates = forward.element_replicates[:, face_indices].sum(axis=1)
        forward_mean = float(forward_replicates.mean())
        forward_stderr = max(
            float(forward_replicates.std(ddof=1) / np.sqrt(forward_replicates.size)),
            float(np.sqrt(np.sum(forward.element_stderr[face_indices] ** 2))))
        forward_cell_mean[cell_index] = forward_mean
        forward_cell_stderr[cell_index] = forward_stderr
        adjoint_cell_mean[cell_index] = adjoint_mean
        adjoint_cell_stderr[cell_index] = adjoint_stderr
        combined_stderr = np.hypot(forward_stderr, adjoint_stderr)
        estimator_discrepancy_sigma[cell_index] = (
            abs(forward_mean - adjoint_mean) / combined_stderr
            if combined_stderr > 0.0 else (0.0 if forward_mean == adjoint_mean else np.inf))
        forward_metric = forward_stderr / max(abs(forward_mean), element_absolute_tolerance)
        adjoint_metric = adjoint_stderr / max(abs(adjoint_mean), element_absolute_tolerance)
        forward_allowed = element_absolute_tolerance + element_relative_tolerance * abs(forward_mean)
        adjoint_allowed = element_absolute_tolerance + element_relative_tolerance * abs(adjoint_mean)
        forward_ok = forward_stderr <= forward_allowed
        # Replicate variance cannot distinguish an exact zero from a rare event missed by every
        # surface-launched adjoint history. The independent forward estimator supplies the missing
        # support check. Unlike direct physical sampling, weighted adjoint histories do not have a
        # finite universal rule-of-three upper bound.
        adjoint_zero_unresolved[cell_index] = (
            adjoint_mean == 0.0 and adjoint_stderr == 0.0
            and forward_mean + support_sigma * forward_stderr > 0.0)
        adjoint_support_unresolved[cell_index] = (
            adjoint_mean + support_sigma * adjoint_stderr
            < support_ratio * max(
                forward_mean - support_sigma * forward_stderr, 0.0))
        adjoint_precision_ok = adjoint_stderr <= adjoint_allowed
        adjoint_ok = (adjoint_precision_ok
                      and not adjoint_zero_unresolved[cell_index]
                      and not adjoint_support_unresolved[cell_index])
        preferred = "forward" if forward_metric < adjoint_metric else "adjoint"
        # Hysteresis is only allowed between two estimators with the same certification status. An
        # old method hint must never retain an uncertified estimate when the complementary estimator
        # meets the requested tolerance.
        if forward_ok != adjoint_ok:
            chosen = "forward" if forward_ok else "adjoint"
        elif method_hint is not None and method_hint[cell_index] in ("forward", "adjoint"):
            hinted = str(method_hint[cell_index])
            hinted_metric = forward_metric if hinted == "forward" else adjoint_metric
            other_metric = adjoint_metric if hinted == "forward" else forward_metric
            chosen = preferred if hinted_metric > switch_factor * other_metric else hinted
        else:
            chosen = preferred
        if chosen == "forward":
            mean, stderr = forward_mean, forward_stderr
            selected_face_mean[face_indices] = forward.element_mean[face_indices]
            selected_face_stderr[face_indices] = forward.element_stderr[face_indices]
            selected_face_replicates[:, face_indices] = forward.element_replicates[:, face_indices]
            selected_endpoint_replicates[:, face_indices] = (
                forward_endpoint_replicates[:, face_indices])
        else:
            mean, stderr = adjoint_mean, adjoint_stderr
            selected_face_mean[face_indices] = adjoint.element_mean[face_indices]
            selected_face_stderr[face_indices] = adjoint.element_stderr[face_indices]
            selected_face_replicates[:, face_indices] = adjoint.element_replicates[:, face_indices]
            selected_endpoint_replicates[:, face_indices] = (
                adjoint_endpoint_replicates[:, face_indices])
        allowed = element_absolute_tolerance + element_relative_tolerance * abs(mean)
        method_within_tolerance[cell_index] = stderr <= allowed
        # A support audit that invalidates the adjoint estimator does not invalidate an independently
        # certified direct-physical forward estimate. Require agreement only when both estimators are
        # actually admissible; otherwise use the admissible direction and retain the diagnostic flag.
        estimator_consistent[cell_index] = (
            not (forward_ok and adjoint_ok)
            or estimator_discrepancy_sigma[cell_index] <= consistency_sigma)
        cell_converged[cell_index] = (
            method_within_tolerance[cell_index] and estimator_consistent[cell_index])
        converged &= cell_converged[cell_index]
        index = first_face[cell]
        per_face[index] = mean; per_face_stderr[index] = stderr; method[cell_index] = chosen
    return dict(
        per_face=per_face, per_face_stderr=per_face_stderr,
        selected_face_mean=selected_face_mean, selected_face_stderr=selected_face_stderr,
        selected_face_replicates=selected_face_replicates,
        selected_endpoint_mean=selected_endpoint_replicates.mean(axis=0),
        selected_endpoint_stderr=np.maximum(
            selected_endpoint_replicates.std(axis=0, ddof=1)
            / np.sqrt(selected_endpoint_replicates.shape[0]),
            selected_face_stderr[:, None]),
        selected_endpoint_replicates=selected_endpoint_replicates,
        unique_cells=np.asarray(unique_cells, dtype=int), method=method,
        forward_cell_mean=forward_cell_mean, forward_cell_stderr=forward_cell_stderr,
        adjoint_cell_mean=adjoint_cell_mean, adjoint_cell_stderr=adjoint_cell_stderr,
        estimator_discrepancy_sigma=estimator_discrepancy_sigma,
        method_within_tolerance=method_within_tolerance,
        estimator_consistent=estimator_consistent, cell_converged=cell_converged,
        adjoint_zero_unresolved=adjoint_zero_unresolved,
        adjoint_support_unresolved=adjoint_support_unresolved,
        consistency_sigma=float(consistency_sigma),
        support_sigma=float(support_sigma),
        support_ratio=float(support_ratio),
        cross_refinement_rounds=cross_refinement_rounds,
        forward_cross_refinement_rounds=forward_cross_refinement_rounds,
        forward_pool_refinement_rounds=forward_pool_refinement_rounds,
        converged=bool(converged), adjoint=adjoint, forward=forward)


def adjoint_boundary_state_floor_flux(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, floor_cells, *,
        proposal_species=None, n_face_position=8, max_steps=None, dt_cap=0.15, dt_field=0.10,
        want_energy=False, fixed_dt=0.0):
    """Backward-compatible horizontal-floor specialization of the arbitrary-face gather."""
    return adjoint_boundary_state_face_flux(
        boundary, species_name, nodal_potential, solid, floor_cells,
        [(0.0, -1.0)] * len(floor_cells), proposal_species=proposal_species,
        n_face_position=n_face_position, max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field,
        want_energy=want_energy, fixed_dt=fixed_dt)
