"""Transport adapters consuming the unified plasma boundary state."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from .boundary_state import qmc_boundary_proposal
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


def forward_boundary_state_cell_flux_qmc(
        boundary: PlasmaBoundaryState, species_name, nodal_potential, solid, cells, *,
        proposal_species=None, log2_samples=12, seed=0, x_min=0.0, x_max=None,
        max_steps=None, dt_cap=0.15, dt_field=0.10):
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
    hit_x, hit_z, *_ = trace_nodal(
        nodal_potential, solid, x0, z0, velocity[:, 0], velocity[:, 2],
        float(species.charge_number), nx, nz,
        int(200 * nz if max_steps is None else max_steps), dt_cap, dt_field)
    log_physical = species.log_flux_density(
        velocity, proposal.phase_rad, proposal.position_m)
    log_proposal = proposal.log_flux_density(
        velocity, proposal.phase_rad, proposal.position_m)
    with np.errstate(over="ignore", invalid="ignore"):
        density_ratio = np.exp(log_physical - log_proposal)
    score = proposal.weight * density_ratio * (float(x_max) - float(x_min))
    cells = [tuple(map(int, cell)) for cell in cells]
    result = np.zeros(len(cells)); lookup = {}
    for index, cell in enumerate(cells):
        lookup.setdefault(cell, []).append(index)
    valid = hit_x >= 0
    for hx, hz, value in zip(hit_x[valid], hit_z[valid], score[valid]):
        indices = lookup.get((int(hx), int(hz)))
        if indices is not None:
            # A cell current is stored once even when the cell owns multiple exposed faces.
            result[indices[0]] += value
    return dict(
        per_cell=result, normalized_total=float(result.sum()),
        absolute_per_cell_m2_s=result * species.flux_m2_s,
        samples=count, seed=int(seed))


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
        face_s = np.tile(face_u - 0.5, base_velocity.shape[0])
    else:
        face_u = face_position_samples
        sample_index = np.arange(base_velocity.shape[0])
        quadrature_weight = proposal.weight
        face_s = face_u - 0.5
    cell_array = np.asarray(cells, dtype=float)
    x_center = cell_array[:, 0] + 0.5 + (0.5 + 1e-3) * normals[:, 0]
    z_center = cell_array[:, 1] + 0.5 + (0.5 + 1e-3) * normals[:, 1]
    x0 = (x_center[:, None] - normals[:, 1, None] * face_s[None, :]).ravel()
    z0 = (z_center[:, None] + normals[:, 0, None] * face_s[None, :]).ravel()
    samples_per_face = sample_index.size; face_count = len(cells)
    vx0 = np.tile(-base_velocity[sample_index, 0], face_count)
    vz0 = np.tile(-base_velocity[sample_index, 2], face_count)
    hit_x, _, _, _, survivor, exit_vx, exit_vz = trace_nodal(
        nodal_potential, solid, x0, z0, vx0, vz0, float(species.charge_number),
        nx, nz, int(max_steps), dt_cap, dt_field)
    escaped = (hit_x < 0) & (survivor < 0.5) & (exit_vz < 0.0)
    tiled_sample_index = np.tile(sample_index, face_count)
    exit_forward = np.column_stack((
        -exit_vx, base_velocity[tiled_sample_index, 1], -exit_vz))
    exit_phase = (None if proposal.phase_rad is None
                  else proposal.phase_rad[tiled_sample_index])
    log_exit_density = species.log_flux_density(exit_forward, exit_phase)
    surface_normal = -(
        normals[:, 0, None] * base_velocity[sample_index, 0][None, :]
        + normals[:, 1, None] * base_velocity[sample_index, 2][None, :]).ravel()
    exit_normal = np.maximum(-exit_vz, 1e-300)
    log_ratio = log_exit_density - np.tile(log_surface_density[sample_index], face_count)
    with np.errstate(over="ignore", invalid="ignore"):
        density_ratio = np.exp(log_ratio)
    contribution = np.where(
        escaped & (surface_normal > 0.0) & np.isfinite(log_ratio),
        surface_normal / exit_normal * density_ratio, 0.0).reshape(face_count, samples_per_face)
    weighted_contribution = contribution * quadrature_weight[None, :]
    per_face[:] = weighted_contribution.sum(axis=1)
    square_sum = np.sum(weighted_contribution * weighted_contribution, axis=1)
    positive = square_sum > 0.0
    effective_sample_size[positive] = per_face[positive] ** 2 / square_sum[positive]
    dominant = np.argmax(weighted_contribution, axis=1)
    dominant_sample_index[positive] = sample_index[dominant[positive]]
    max_sample_fraction[positive] = (
        weighted_contribution[np.where(positive)[0], dominant[positive]] / per_face[positive])
    dominant_surface_velocity[positive] = base_velocity[dominant_sample_index[positive]]
    exit_matrix = exit_forward.reshape(face_count, samples_per_face, 3)
    dominant_exit_velocity[positive] = exit_matrix[np.where(positive)[0], dominant[positive]]
    if want_energy:
        impact_energy = np.sum(base_velocity[sample_index] ** 2, axis=1)
        energy_numerator[:] = np.sum(weighted_contribution * impact_energy[None, :], axis=1)
    normalized_flux = float(per_face.mean())
    result = dict(normalized_flux=normalized_flux,
                  absolute_flux_m2_s=normalized_flux * species.flux_m2_s,
                  per_face=per_face, per_cell=per_face,
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
        max_steps=None, dt_cap=0.15, dt_field=0.10, initial_log2_samples=None):
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
        element_relative_tolerance=element_relative_tolerance, refine_fraction=refine_fraction,
        initial_log2_samples=initial_log2_samples)


def adaptive_forward_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, cells, *, proposal_species=None,
        base_log2=8, max_log2=16, n_replicates=4, seed=0,
        absolute_tolerance=1e-3, relative_tolerance=5e-3,
        element_absolute_tolerance=None, element_relative_tolerance=0.0,
        refine_fraction=0.5, initial_log2_samples=None,
        x_min=0.0, x_max=None, max_steps=None, dt_cap=0.15, dt_field=0.10):
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
                log2_samples=level, seed=int(seed + replicate), x_min=x_min, x_max=x_max,
                max_steps=max_steps, dt_cap=dt_cap, dt_field=dt_field)
            estimates[replicate] = result["per_cell"]
        evaluations += n_replicates * (2 ** level)
        element_mean = estimates.mean(axis=0)
        element_stderr = estimates.std(axis=0, ddof=1) / np.sqrt(n_replicates)
        # Zero observed hits do not imply zero flux. For direct physical sampling, the rule-of-three
        # 95% upper probability bound is 3/N across the pooled independent histories.
        zero = np.all(estimates == 0.0, axis=0)
        zero_upper = 3.0 * source_width / (n_replicates * (2 ** level))
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
        rounds=rounds, evaluations=evaluations)


def bidirectional_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, cells, normals, *,
        proposal_species=None, adjoint_options=None, forward_options=None,
        element_absolute_tolerance=1e-3, element_relative_tolerance=0.05,
        method_hint=None, switch_factor=2.0):
    """Select forward or adjoint current per physical cell solely by measured uncertainty.

    This is not a named-region switch. Both unbiased estimators use the same physical boundary density;
    the lower-uncertainty estimate wins independently for each cell. The result refuses cells for which
    neither direction meets the requested mixed tolerance.
    """
    cells = [tuple(map(int, cell)) for cell in cells]
    normals = np.asarray(normals, dtype=float)
    adjoint_kwargs = {} if adjoint_options is None else dict(adjoint_options)
    forward_kwargs = {} if forward_options is None else dict(forward_options)
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
    forward = adaptive_forward_boundary_state_cell_flux(
        boundary, species_name, nodal_potential, solid, unique_cells,
        proposal_species=None, **forward_kwargs)
    first_face = {cell: cells.index(cell) for cell in unique_cells}
    per_face = np.zeros(len(cells)); per_face_stderr = np.zeros(len(cells))
    method = np.empty(len(unique_cells), dtype="U7")
    converged = True
    for cell_index, cell in enumerate(unique_cells):
        face_indices = np.array([i for i, item in enumerate(cells) if item == cell], dtype=int)
        adjoint_replicates = adjoint.element_replicates[:, face_indices].sum(axis=1)
        adjoint_mean = float(adjoint_replicates.mean())
        adjoint_stderr = float(adjoint_replicates.std(ddof=1) / np.sqrt(adjoint_replicates.size))
        forward_mean = float(forward.element_mean[cell_index])
        forward_stderr = float(forward.element_stderr[cell_index])
        forward_metric = forward_stderr / max(abs(forward_mean), element_absolute_tolerance)
        adjoint_metric = adjoint_stderr / max(abs(adjoint_mean), element_absolute_tolerance)
        preferred = "forward" if forward_metric < adjoint_metric else "adjoint"
        if method_hint is not None and method_hint[cell_index] in ("forward", "adjoint"):
            hinted = str(method_hint[cell_index])
            hinted_metric = forward_metric if hinted == "forward" else adjoint_metric
            other_metric = adjoint_metric if hinted == "forward" else forward_metric
            chosen = preferred if hinted_metric > switch_factor * other_metric else hinted
        else:
            chosen = preferred
        if chosen == "forward":
            mean, stderr = forward_mean, forward_stderr
        else:
            mean, stderr = adjoint_mean, adjoint_stderr
        allowed = element_absolute_tolerance + element_relative_tolerance * abs(mean)
        converged &= stderr <= allowed
        index = first_face[cell]
        per_face[index] = mean; per_face_stderr[index] = stderr; method[cell_index] = chosen
    return dict(
        per_face=per_face, per_face_stderr=per_face_stderr,
        unique_cells=np.asarray(unique_cells, dtype=int), method=method,
        converged=bool(converged), adjoint=adjoint, forward=forward)


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
