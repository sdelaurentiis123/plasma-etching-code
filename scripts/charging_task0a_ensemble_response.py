"""Paired scrambled-QMC response audit for the real AR4 charging operator.

This driver intentionally evaluates currents without taking a charging step.  Estimator selection is
loaded from a separately certified checkpoint and never reselected.  Every +/- pair uses identical
Sobol seeds, levels, proposals, and forward launch coordinates (common random numbers).

The default direction set is deliberately small but diagnostic: the two worst residual coordinates,
five fixed random unit vectors, and the first three right-singular vectors of the archived Jacobian.
The reported condition number is therefore the condition of this declared 10-dimensional response
subspace, not a mislabeled full 47 x 47 Jacobian condition number.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path

# Four independent processes are the supported local parallel path. Numba's workqueue layer aborts
# under concurrent Python threads, so ThreadPoolExecutor must not be used for this campaign.
# Thread-count environment variables must be fixed before Numba is imported.  Do that only for
# direct campaign execution: importing this module for its audit helpers must not mutate the host
# test process after another module has already initialized Numba's worker pool.
if __name__ == "__main__":
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[variable] = "1"

import numpy as np

from petch.boundary_state import (
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    folded_normal_tangential_proposal,
    maxwellian_electron_boundary_state,
    mixture_boundary_proposal,
    qmc_boundary_proposal,
    qmc_boundary_proposal_with_auxiliary,
)
from petch.boundary_transport import (
    adjoint_boundary_state_face_flux,
    forward_boundary_state_cell_flux_qmc,
)
from petch.charging_backward import _gas_faces
from petch.charging_nodal import material_face_nodes
from petch.charging_poisson import NodalPoissonSystem
from petch.sheath import CollisionlessRFSheath


RADII = (1.0, 0.5, 0.25, 0.1, 0.05, 0.025)
SEEDS = (401, 409, 419, 421, 431, 433, 439, 443)
_WORKER_MODEL = None


@dataclass(frozen=True)
class AuditModel:
    solid: np.ndarray
    boundary: PlasmaBoundaryState
    proposals: dict
    poisson: NodalPoissonSystem
    base_charge: np.ndarray
    cells: np.ndarray
    normals: np.ndarray
    face_nodes: tuple
    dielectric_nodes: np.ndarray
    capacitance: np.ndarray
    method: dict
    current_scale: float


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _physical_boundary_and_proposals():
    sheath = CollisionlessRFSheath(
        40.0, 10.0, 2e6, 4.0, 40.0, thickness_m=5e-4)
    ion = collisionless_sheath_boundary_state(
        sheath, 1e19, n_phase=16, tangential_temperature_eV=0.2,
        n_transverse=3, normal_energy_bins=64).get("ion")
    electron = maxwellian_electron_boundary_state(
        4.0, 1e19, n_transverse=3, n_normal=6).get("electron")
    boundary = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    ion_tail = SpeciesBoundaryState(
        "ion_tail", 1, 40.0, 1.0, [[0.0, 0.0, np.sqrt(40.0)]], [1.0],
        density_model=MaxwellianFluxVelocityDensity(40.0))
    electron_tail = SpeciesBoundaryState(
        "electron_tail", -1, electron.mass_amu, 1.0,
        [[0.0, 0.0, np.sqrt(16.0)]], [1.0],
        density_model=MaxwellianFluxVelocityDensity(16.0))
    ion_proposal = mixture_boundary_proposal((
        ion, folded_normal_tangential_proposal(ion, +1),
        folded_normal_tangential_proposal(ion, -1), ion_tail,
    ), (0.4, 0.25, 0.25, 0.1), name="ion_grazing_proposal")
    electron_proposal = mixture_boundary_proposal((
        electron, folded_normal_tangential_proposal(electron, +1),
        folded_normal_tangential_proposal(electron, -1), electron_tail,
    ), (0.4, 0.25, 0.25, 0.1), name="electron_grazing_proposal")
    return boundary, {"ion": ion_proposal, "electron": electron_proposal}


def load_model(checkpoint_path):
    with np.load(checkpoint_path) as saved:
        required = {
            "solid", "surface_charge_node_c_per_m", "dielectric_nodes",
            "method_hint_ion", "method_hint_electron",
        }
        missing = required - set(saved.files)
        if missing:
            raise ValueError(f"checkpoint lacks Task 0A fields: {sorted(missing)}")
        solid = saved["solid"].astype(bool)
        base_charge = saved["surface_charge_node_c_per_m"].astype(float)
        dielectric_nodes = saved["dielectric_nodes"].astype(int)
        method = {
            "ion": saved["method_hint_ion"].astype("U7"),
            "electron": saved["method_hint_electron"].astype("U7"),
        }
    cells, normals = _gas_faces(solid, solid)
    cells = np.asarray(cells, dtype=int); normals = np.asarray(normals, dtype=int)
    unique_cells = list(dict.fromkeys(map(tuple, cells)))
    if any(values.shape != (len(unique_cells),) for values in method.values()):
        raise ValueError("frozen estimator map does not match checkpoint surface cells")
    face_nodes = tuple(
        material_face_nodes(tuple(cell), tuple(normal))
        for cell, normal in zip(cells, normals))
    expected_nodes = np.asarray(sorted({
        node for endpoints in face_nodes for node in endpoints if node[1] != 0
    }), dtype=int)
    if not np.array_equal(dielectric_nodes, expected_nodes):
        raise ValueError("checkpoint dielectric DOF ordering differs from the engine ordering")
    epsilon_r = np.ones_like(solid, dtype=float); epsilon_r[solid] = 3.9
    grounded = np.zeros((solid.shape[0] + 1, solid.shape[1] + 1), dtype=bool)
    grounded[:, -1] = True
    poisson = NodalPoissonSystem(
        epsilon_r, grounded, np.zeros_like(grounded, dtype=float))
    capacitance = poisson.diagonal_surface_capacitance(dielectric_nodes)
    boundary, proposals = _physical_boundary_and_proposals()
    current_scale = max(
        species.flux_m2_s * abs(species.charge_number)
        for species in boundary.species)
    return AuditModel(
        solid, boundary, proposals, poisson, base_charge, cells, normals, face_nodes,
        dielectric_nodes, capacitance, method, current_scale)


def direction_set(jacobian_path, count_random=5, count_dominant=3, residual_override=None):
    with np.load(jacobian_path) as saved:
        jacobian = np.asarray(saved["jacobian"], dtype=float)
        archived_residual = np.asarray(saved["residual"], dtype=float)
    residual = (archived_residual if residual_override is None
                else np.asarray(residual_override, dtype=float))
    if jacobian.shape != (residual.size, residual.size):
        raise ValueError("archived Jacobian and residual dimensions differ")
    directions = []
    for rank, index in enumerate(np.argsort(-np.abs(residual))[:2]):
        vector = np.zeros(residual.size); vector[index] = 1.0
        directions.append((f"worst_coordinate_{rank}_dof_{index:02d}", "worst_coordinate", vector))
    rng = np.random.default_rng(20260713)
    for index in range(int(count_random)):
        vector = rng.normal(size=residual.size); vector /= np.linalg.norm(vector)
        directions.append((f"random_{index}", "random", vector))
    _, _, vh = np.linalg.svd(jacobian, full_matrices=False)
    for index in range(int(count_dominant)):
        directions.append((f"dominant_{index}", "dominant", vh[index].copy()))
    return residual, directions


def _checkpoint_residual(checkpoint_path):
    with np.load(checkpoint_path) as saved:
        if "mean_log_ratio_history" in saved and len(saved["mean_log_ratio_history"]):
            return np.asarray(saved["mean_log_ratio_history"][-1], dtype=float)
        ion = np.asarray(saved["ion_current"], dtype=float)
        electron = np.asarray(saved["electron_current"], dtype=float)
        return np.log((ion + 1e-12) / (electron + 1e-12))


def _method_per_face(model, species_name):
    unique_cells = list(dict.fromkeys(map(tuple, model.cells)))
    by_cell = dict(zip(unique_cells, model.method[species_name]))
    return np.asarray([by_cell[tuple(cell)] for cell in model.cells])


def _species_evaluation(
        model, species_name, potential, log2_samples, seed, max_steps=None,
        fixed_dt=0.01):
    method_face = _method_per_face(model, species_name)
    endpoint = np.zeros((len(model.cells), 2))
    detail = {}
    if np.any(method_face == "forward"):
        forward = forward_boundary_state_cell_flux_qmc(
            model.boundary, species_name, potential, model.solid, model.cells,
            normals=model.normals, log2_samples=log2_samples, seed=seed,
            fixed_dt=fixed_dt, max_steps=max_steps, source_offset=1e-6,
            return_trajectory_outcomes=True, return_trajectory_contributions=True)
        selected = method_face == "forward"
        endpoint[selected] = forward["per_face_endpoint"][selected]
        detail["forward"] = forward
    adjoint_faces = np.where(method_face == "adjoint")[0]
    if adjoint_faces.size:
        proposal, auxiliary = qmc_boundary_proposal_with_auxiliary(
            model.proposals[species_name], int(log2_samples), 1, int(seed),
            name=f"{species_name}-task0a")
        adjoint = adjoint_boundary_state_face_flux(
            model.boundary, species_name, potential, model.solid,
            model.cells[adjoint_faces], model.normals[adjoint_faces],
            proposal_species=proposal, face_position_samples=auxiliary[:, 0],
            fixed_dt=fixed_dt, max_steps=max_steps, face_offset=1e-6,
            return_trajectory_outcomes=True, return_trajectory_contributions=True)
        endpoint[adjoint_faces] = adjoint["per_face_endpoint"]
        detail["adjoint"] = adjoint
        detail["adjoint_faces"] = adjoint_faces
    return endpoint, detail


def _assemble_nodes(model, species_name, endpoint):
    species = model.boundary.get(species_name)
    scale = species.flux_m2_s * abs(species.charge_number) / model.current_scale
    node_index = {tuple(node): index for index, node in enumerate(model.dielectric_nodes)}
    current = np.zeros(len(model.dielectric_nodes))
    for face_value, endpoints in zip(endpoint, model.face_nodes):
        for endpoint_index, node in enumerate(endpoints):
            if node[1] != 0:
                current[node_index[tuple(node)]] += scale * face_value[endpoint_index]
    return current


def _evaluate(model, direction, radius, sign, log2_samples, seed):
    charge = model.base_charge.copy()
    for index, node in enumerate(model.dielectric_nodes):
        charge[tuple(node)] += sign * radius * direction[index] * model.capacitance[index]
    potential, _ = model.poisson.solve(charge)
    endpoints = {}; details = {}; currents = {}
    for species_name in ("ion", "electron"):
        endpoints[species_name], details[species_name] = _species_evaluation(
            model, species_name, potential, log2_samples, seed)
        currents[species_name] = _assemble_nodes(
            model, species_name, endpoints[species_name])
    residual = np.log(
        (currents["ion"] + 1e-12) / (currents["electron"] + 1e-12))
    return dict(
        potential=potential, endpoints=endpoints, details=details,
        ion_current=currents["ion"], electron_current=currents["electron"],
        residual=residual)


def _deposit_forward_subset(endpoint, result, selected_histories, allowed_faces):
    face = result["trajectory_face_index"]
    score = result["trajectory_score"]
    u = result["trajectory_endpoint_u"]
    selected = np.where(
        selected_histories & (face >= 0)
        & np.where(face >= 0, allowed_faces[np.maximum(face, 0)], False))[0]
    for history in selected:
        index = face[history]
        endpoint[index, 0] += score[history] * (1.0 - u[history])
        endpoint[index, 1] += score[history] * u[history]


def _switch_endpoints(model, minus, plus, species_name):
    minus_endpoint = np.zeros((len(model.cells), 2))
    plus_endpoint = np.zeros((len(model.cells), 2))
    switched = 0; histories = 0
    minus_detail = minus["details"][species_name]
    plus_detail = plus["details"][species_name]
    if "forward" in minus_detail:
        minus_forward = minus_detail["forward"]; plus_forward = plus_detail["forward"]
        changed = np.any(
            minus_forward["trajectory_outcomes"] != plus_forward["trajectory_outcomes"], axis=1)
        switched += int(np.count_nonzero(changed)); histories += int(changed.size)
        forward_faces = _method_per_face(model, species_name) == "forward"
        _deposit_forward_subset(minus_endpoint, minus_forward, changed, forward_faces)
        _deposit_forward_subset(plus_endpoint, plus_forward, changed, forward_faces)
    if "adjoint" in minus_detail:
        minus_adjoint = minus_detail["adjoint"]; plus_adjoint = plus_detail["adjoint"]
        changed = np.any(
            minus_adjoint["trajectory_outcomes"] != plus_adjoint["trajectory_outcomes"], axis=2)
        switched += int(np.count_nonzero(changed)); histories += int(changed.size)
        faces = minus_detail["adjoint_faces"]
        face_u = minus_adjoint["trajectory_face_u"]
        for local, face in enumerate(faces):
            mask = changed[local]
            for target, source in (
                    (minus_endpoint, minus_adjoint), (plus_endpoint, plus_adjoint)):
                contribution = source["trajectory_contribution"][local, mask]
                u = face_u[mask]
                target[face, 0] += np.sum(contribution * (1.0 - u))
                target[face, 1] += np.sum(contribution * u)
    return minus_endpoint, plus_endpoint, switched, histories


def evaluate_pair(model, direction, radius, log2_samples, seed):
    start = time.perf_counter()
    minus = _evaluate(model, direction, radius, -1.0, log2_samples, seed)
    plus = _evaluate(model, direction, radius, +1.0, log2_samples, seed)
    denominator = 2.0 * radius
    response = (plus["residual"] - minus["residual"]) / denominator
    current_response = (
        (plus["ion_current"] - plus["electron_current"])
        - (minus["ion_current"] - minus["electron_current"])) / denominator
    switch_currents = {}; switched = 0; histories = 0
    for species_name in ("ion", "electron"):
        minus_endpoint, plus_endpoint, count, total = _switch_endpoints(
            model, minus, plus, species_name)
        switched += count; histories += total
        switch_currents[species_name] = (
            _assemble_nodes(model, species_name, plus_endpoint)
            - _assemble_nodes(model, species_name, minus_endpoint)) / denominator
    midpoint_ion = 0.5 * (plus["ion_current"] + minus["ion_current"])
    midpoint_electron = 0.5 * (plus["electron_current"] + minus["electron_current"])
    switch_response = _switch_log_response(
        switch_currents["ion"], switch_currents["electron"],
        midpoint_ion, midpoint_electron, 1e-4)
    return dict(
        response=response, switch_response=switch_response,
        switch_ion_current_response=switch_currents["ion"],
        switch_electron_current_response=switch_currents["electron"],
        midpoint_ion_current=midpoint_ion, midpoint_electron_current=midpoint_electron,
        current_response=current_response,
        global_current_response=float(np.sum(current_response)),
        face_switch_fraction=float(switched / histories) if histories else 0.0,
        switched_histories=switched, total_histories=histories,
        minus_residual=minus["residual"], plus_residual=plus["residual"],
        elapsed_seconds=time.perf_counter() - start)


def _switch_log_response(ion_response, electron_response, midpoint_ion, midpoint_electron,
                         active_threshold):
    result = np.zeros_like(np.asarray(ion_response, dtype=float))
    active = np.asarray(midpoint_ion) + np.asarray(midpoint_electron) >= float(active_threshold)
    result[active] = (
        np.asarray(ion_response)[active] / np.maximum(np.asarray(midpoint_ion)[active], 1e-12)
        - np.asarray(electron_response)[active]
        / np.maximum(np.asarray(midpoint_electron)[active], 1e-12))
    return result


def _initialize_worker(checkpoint):
    global _WORKER_MODEL
    _WORKER_MODEL = load_model(checkpoint)


def _evaluate_job(job):
    level, seed, radius, name, kind, vector = job
    return job, evaluate_pair(_WORKER_MODEL, vector, radius, level, seed)


def _job_path(output_dir, level, seed, radius, direction_name):
    radius_text = str(radius).replace(".", "p")
    return output_dir / "raw" / (
        f"l{level:02d}_s{seed}_r{radius_text}_{direction_name}.npz")


def run(args):
    output_dir = Path(args.output_dir); (output_dir / "raw").mkdir(parents=True, exist_ok=True)
    # Validate once in the parent before launching expensive work. Each process owns its sparse
    # Poisson factorization; sharing SuperLU state across threads is not assumed safe.
    load_model(args.checkpoint)
    checkpoint_residual = _checkpoint_residual(args.checkpoint)
    residual, directions = direction_set(
        args.jacobian, residual_override=checkpoint_residual)
    if args.direction_kind:
        directions = [item for item in directions if item[1] in set(args.direction_kind)]
    if args.direction:
        directions = [item for item in directions if item[0] in set(args.direction)]
        if not directions:
            raise ValueError("--direction did not match the declared audit directions")
    config = dict(
        schema="petch.charging.task0a.v1", checkpoint=Path(args.checkpoint).name,
        checkpoint_sha256=_sha256(args.checkpoint),
        jacobian=Path(args.jacobian).name, jacobian_sha256=_sha256(args.jacobian),
        levels=list(args.levels), seeds=list(args.seeds), radii=list(args.radii),
        directions=[dict(name=name, kind=kind, vector=vector.tolist())
                    for name, kind, vector in directions],
        common_random_numbers=True, estimator_method_map="checkpoint_frozen",
        hard_visibility=True, fixed_dt=0.01, face_offset=1e-6,
        nested_sobol_levels=True, active_threshold_sweep=[1e-5, 1e-4, 1e-3],
        process_workers=args.workers, process_thread_count=1)
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config["config_hash"] = hashlib.sha256(canonical).hexdigest()
    (output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    jobs = [(level, seed, radius, direction)
            for level in args.levels for seed in args.seeds for radius in args.radii
            for direction in directions]
    pending = []
    for level, seed, radius, (name, kind, vector) in jobs:
        path = _job_path(output_dir, level, seed, radius, name)
        if path.exists() and not args.force:
            continue
        pending.append((level, seed, radius, name, kind, vector))
    if args.max_jobs:
        pending = pending[:args.max_jobs]
    completed = 0
    executor = ProcessPoolExecutor(
        max_workers=args.workers, initializer=_initialize_worker,
        initargs=(args.checkpoint,))
    try:
        evaluations = executor.map(_evaluate_job, pending)
        for index, (job, result) in enumerate(evaluations):
            level, seed, radius, name, kind, vector = job
            path = _job_path(output_dir, level, seed, radius, name)
            np.savez_compressed(
                path, schema="petch.charging.task0a.pair.v1", config_hash=config["config_hash"],
                direction_name=name, direction_kind=kind, direction=vector,
                level=level, seed=seed, radius_v=radius, **result)
            completed += 1
            print(
                f"[{index + 1}/{len(pending)}] {path.name} "
                f"L={np.linalg.norm(result['response']):.5g} "
                f"switch={result['face_switch_fraction']:.3g} "
                f"seconds={result['elapsed_seconds']:.1f}", flush=True)
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
    print(f"wrote {completed} paired evaluations; checkpoint residual rms="
          f"{np.sqrt(np.mean(residual ** 2)):.6g}")


def _load_rows(output_dir):
    rows = []
    for path in sorted((output_dir / "raw").glob("*.npz")):
        with np.load(path) as saved:
            rows.append({key: saved[key] for key in saved.files})
    return rows


def _condition(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return float(singular[0] / singular[-1]), singular


def _switch_response_for(row, active_threshold):
    required = {
        "switch_ion_current_response", "switch_electron_current_response",
        "midpoint_ion_current", "midpoint_electron_current",
    }
    if required <= set(row):
        return _switch_log_response(
            row["switch_ion_current_response"], row["switch_electron_current_response"],
            row["midpoint_ion_current"], row["midpoint_electron_current"], active_threshold)
    if np.isclose(active_threshold, 1e-4):
        return np.asarray(row["switch_response"], dtype=float)
    raise ValueError("raw pair predates active-threshold sensitivity fields")


def response_decomposition_metrics(total_matrix, switch_matrix):
    """Signed norm decomposition; no clamping can turn an adverse switch effect into zero."""
    total = np.asarray(total_matrix, dtype=float)
    switch = np.asarray(switch_matrix, dtype=float)
    nonswitch = total - switch
    total_norm2 = float(np.sum(total * total))
    component_norm2 = float(np.sum(switch * switch) + np.sum(nonswitch * nonswitch))
    condition, singular = _condition(total)
    no_switch_condition, no_switch_singular = _condition(nonswitch)
    return dict(
        condition=condition, no_switch_condition=no_switch_condition,
        switch_component_energy_fraction=float(np.sum(switch * switch))
        / max(component_norm2, 1e-300),
        switch_signed_projection_fraction=float(np.sum(switch * total))
        / max(total_norm2, 1e-300),
        switch_nonswitch_interference_fraction=float(2.0 * np.sum(switch * nonswitch))
        / max(total_norm2, 1e-300),
        signed_log_condition_change=float(np.log(condition / no_switch_condition)),
        sigma_max=float(singular[0]), sigma_min=float(singular[-1]),
        no_switch_sigma_min=float(no_switch_singular[-1]))


def _bootstrap_metrics(response_replicates, switch_replicates, samples=400):
    count = response_replicates.shape[0]
    rng = np.random.default_rng(20260713)
    values = []
    for _ in range(int(samples)):
        selected = rng.integers(0, count, size=count)
        values.append(response_decomposition_metrics(
            response_replicates[selected].mean(axis=0),
            switch_replicates[selected].mean(axis=0)))
    keys = ("condition", "no_switch_condition", "signed_log_condition_change",
            "switch_component_energy_fraction", "switch_signed_projection_fraction")
    return {key: np.asarray([item[key] for item in values]) for key in keys}


def analyze(args):
    output_dir = Path(args.output_dir)
    config = json.loads((output_dir / "config.json").read_text())
    rows = _load_rows(output_dir)
    expected = len(config["levels"]) * len(config["seeds"]) * len(config["radii"]) * len(config["directions"])
    if len(rows) != expected and not args.allow_incomplete:
        raise RuntimeError(f"Task 0A incomplete: found {len(rows)} of {expected} paired jobs")
    grouped = {}
    for row in rows:
        key = (int(row["level"]), float(row["radius_v"]), str(row["direction_name"]))
        grouped.setdefault(key, []).append(row)
    table = []
    for (level, radius, name), items in sorted(grouped.items()):
        response = np.stack([item["response"] for item in items])
        switch_response = np.stack([_switch_response_for(item, 1e-4) for item in items])
        mean = response.mean(axis=0)
        stderr = (response.std(axis=0, ddof=1) / np.sqrt(len(items))
                  if len(items) > 1 else np.full(mean.shape, np.nan))
        table.append(dict(
            level=level, radius_v=radius, direction=name,
            direction_kind=str(items[0]["direction_kind"]), replicates=len(items),
            frozen_lipschitz=float(np.linalg.norm(response[0])),
            ensemble_lipschitz=float(np.linalg.norm(mean)),
            between_scramble_error=float(np.linalg.norm(stderr)),
            signal_to_sampling_error=float(np.linalg.norm(mean) / np.linalg.norm(stderr))
            if np.all(np.isfinite(stderr)) and np.linalg.norm(stderr) > 0 else np.inf,
            local_current_change_l2=float(np.linalg.norm(np.mean(
                np.stack([item["current_response"] for item in items]), axis=0))),
            global_current_change=float(np.mean([
                float(item["global_current_response"]) for item in items])),
            face_switch_fraction=float(np.mean([
                float(item["face_switch_fraction"]) for item in items])),
            switch_response_fraction=float(
                np.linalg.norm(switch_response.mean(axis=0)) / max(np.linalg.norm(mean), 1e-300)),
            switch_signed_projection_fraction=float(np.dot(
                switch_response.mean(axis=0), mean) / max(np.dot(mean, mean), 1e-300))))
    table_path = output_dir / "ensemble_response.csv"
    if table:
        with open(table_path, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(table)
    condition_rows = []
    direction_names = [item["name"] for item in config["directions"]]
    thresholds = config.get("active_threshold_sweep", [1e-4])
    for level in config["levels"]:
        for radius in config["radii"]:
            by_direction = []; available = []
            for name in direction_names:
                items = sorted(
                    grouped.get((int(level), float(radius), name), []),
                    key=lambda item: int(item["seed"]))
                if items:
                    by_direction.append(items)
                    available.append(name)
            if len(by_direction) < 2:
                continue
            replicate_count = min(map(len, by_direction))
            response_replicates = np.stack([
                np.column_stack([items[index]["response"] for items in by_direction])
                for index in range(replicate_count)])
            frozen_condition, _ = _condition(response_replicates[0])
            for threshold in thresholds:
                try:
                    switch_replicates = np.stack([
                        np.column_stack([
                            _switch_response_for(items[index], threshold)
                            for items in by_direction])
                        for index in range(replicate_count)])
                except ValueError:
                    continue
                metrics = response_decomposition_metrics(
                    response_replicates.mean(axis=0), switch_replicates.mean(axis=0))
                bootstrap = _bootstrap_metrics(response_replicates, switch_replicates)
                condition_rows.append(dict(
                    level=level, radius_v=radius, active_threshold=threshold,
                    directions=len(by_direction), replicates=replicate_count,
                    complete_direction_set=len(by_direction) == len(direction_names),
                    frozen_directional_condition=frozen_condition,
                    ensemble_directional_condition=metrics["condition"],
                    ensemble_condition_stderr=float(bootstrap["condition"].std(ddof=1)),
                    ensemble_condition_ci_low=float(np.quantile(bootstrap["condition"], 0.025)),
                    ensemble_condition_ci_high=float(np.quantile(bootstrap["condition"], 0.975)),
                    no_switch_directional_condition=metrics["no_switch_condition"],
                    no_switch_condition_stderr=float(
                        bootstrap["no_switch_condition"].std(ddof=1)),
                    signed_log_condition_change=metrics["signed_log_condition_change"],
                    signed_log_condition_change_stderr=float(
                        bootstrap["signed_log_condition_change"].std(ddof=1)),
                    switch_component_energy_fraction=metrics[
                        "switch_component_energy_fraction"],
                    switch_component_energy_stderr=float(
                        bootstrap["switch_component_energy_fraction"].std(ddof=1)),
                    switch_signed_projection_fraction=metrics[
                        "switch_signed_projection_fraction"],
                    switch_signed_projection_stderr=float(
                        bootstrap["switch_signed_projection_fraction"].std(ddof=1)),
                    sigma_max=metrics["sigma_max"], sigma_min=metrics["sigma_min"],
                    direction_names=";".join(available)))
    condition_path = output_dir / "ensemble_condition.csv"
    if condition_rows:
        with open(condition_path, "w", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(condition_rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(condition_rows)
    scaling_rows = _switch_scaling_rows(table)
    _write_rows(output_dir / "switch_scaling.csv", scaling_rows)
    stability_rows = _level_stability_rows(grouped, config)
    _write_rows(output_dir / "level_stability.csv", stability_rows)
    highest_level = max((row["level"] for row in table), default=-1)
    highest = [row for row in table if row["level"] == highest_level]
    min_signal = min((row["signal_to_sampling_error"] for row in highest), default=0.0)
    stability_sigma = max((row["difference_sigma"] for row in stability_rows), default=np.inf)
    precision_pass = bool(min_signal >= 3.0 and stability_sigma <= 2.0)
    summary = dict(
        schema="petch.charging.task0a.summary.v1", config_hash=config["config_hash"],
        complete=len(rows) == expected, paired_jobs=len(rows), expected_paired_jobs=expected,
        exact_operator="hard visibility; no smoothing",
        condition_scope=f"declared {len(direction_names)}-direction response subspace",
        frozen_ensemble_comparison="same checkpoint, directions, radii, levels, and first scramble",
        checkpoint_residual_rms=float(np.sqrt(np.mean(_checkpoint_residual(
            output_dir.parent / "charging_task0_inputs" / config["checkpoint"]) ** 2)))
            if (output_dir.parent / "charging_task0_inputs" / config["checkpoint"]).exists()
            else None,
        highest_level=highest_level, minimum_high_level_signal_to_error=min_signal,
        maximum_level_difference_sigma=stability_sigma,
        precision_gate_pass=precision_pass,
        decision_gate=("evaluate smoother/conditioning branch" if precision_pass else
                       "increase paired-replicate precision; no solver promotion"),
        switch_attribution=(
            "signed projection and component energy; signed log(cond_total/cond_no_switch); "
            "no clamping"))
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _plot(output_dir, table, condition_rows, scaling_rows, summary["complete"])
    print(json.dumps(summary, indent=2))


def _write_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _switch_scaling_rows(table):
    result = []
    keys = sorted({(row["level"], row["direction"]) for row in table})
    for level, direction in keys:
        rows = sorted(
            [row for row in table if row["level"] == level and row["direction"] == direction],
            key=lambda row: row["radius_v"])
        positive = [row for row in rows if row["face_switch_fraction"] > 0.0]
        if len(positive) < 2:
            continue
        slope, intercept = np.polyfit(
            np.log([row["radius_v"] for row in positive]),
            np.log([row["face_switch_fraction"] for row in positive]), 1)
        result.append(dict(
            level=level, direction=direction, direction_kind=positive[0]["direction_kind"],
            radius_count=len(positive), switch_fraction_loglog_slope=float(slope),
            log_intercept=float(intercept)))
    return result


def _level_stability_rows(grouped, config):
    levels = sorted(config["levels"])
    if len(levels) < 2:
        return []
    low, high = levels[-2:]
    result = []
    for radius in config["radii"]:
        for item in config["directions"]:
            low_items = grouped.get((int(low), float(radius), item["name"]), [])
            high_items = grouped.get((int(high), float(radius), item["name"]), [])
            if len(low_items) < 2 or len(high_items) < 2:
                continue
            low_response = np.stack([row["response"] for row in low_items])
            high_response = np.stack([row["response"] for row in high_items])
            difference = high_response.mean(axis=0) - low_response.mean(axis=0)
            error = np.sqrt(
                low_response.var(axis=0, ddof=1) / len(low_items)
                + high_response.var(axis=0, ddof=1) / len(high_items))
            result.append(dict(
                low_level=low, high_level=high, radius_v=radius, direction=item["name"],
                difference_l2=float(np.linalg.norm(difference)),
                combined_error_l2=float(np.linalg.norm(error)),
                difference_sigma=float(np.linalg.norm(difference) / max(np.linalg.norm(error), 1e-300))))
    return result


def _plot(output_dir, table, condition_rows, scaling_rows, complete):
    if not table:
        return
    import matplotlib.pyplot as plt
    levels = sorted({row["level"] for row in table})
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for level in levels:
        radii = sorted({row["radius_v"] for row in table if row["level"] == level})
        aggregate = []
        for radius in radii:
            subset = [row for row in table if row["level"] == level and row["radius_v"] == radius]
            aggregate.append((
                radius, np.median([row["frozen_lipschitz"] for row in subset]),
                np.median([row["ensemble_lipschitz"] for row in subset]),
                np.median([row["signal_to_sampling_error"] for row in subset]),
                np.mean([row["face_switch_fraction"] for row in subset])))
        values = np.asarray(aggregate)
        axes[0, 0].plot(values[:, 0], values[:, 1], "o--", label=f"frozen 2^{level}")
        axes[0, 0].plot(values[:, 0], values[:, 2], "o-", label=f"ensemble 2^{level}")
        axes[0, 1].plot(values[:, 0], values[:, 3], "o-", label=f"2^{level}")
        axes[1, 0].plot(values[:, 0], values[:, 4], "o-", label=f"2^{level}")
        condition = [row for row in condition_rows
                     if row["level"] == level and np.isclose(row["active_threshold"], 1e-4)]
        if condition:
            radius = np.asarray([row["radius_v"] for row in condition])
            ensemble = np.asarray([row["ensemble_directional_condition"] for row in condition])
            axes[1, 1].plot(radius, ensemble, "o-", label=f"ensemble 2^{level}")
            axes[1, 1].plot(
                radius, [row["frozen_directional_condition"] for row in condition],
                "x--", label=f"frozen 2^{level}")
            axes[1, 1].fill_between(
                radius, [row["ensemble_condition_ci_low"] for row in condition],
                [row["ensemble_condition_ci_high"] for row in condition], alpha=0.15)
    for axis in axes.ravel():
        axis.set_xscale("log"); axis.invert_xaxis(); axis.grid(alpha=0.25); axis.legend()
        axis.set_xlabel("perturbation radius (V)")
    axes[0, 0].set_yscale("log"); axes[0, 0].set_ylabel("apparent response norm (V$^{-1}$)")
    axes[0, 1].set_yscale("log"); axes[0, 1].set_ylabel("signal / between-scramble error")
    axes[1, 0].set_yscale("log"); axes[1, 0].set_ylabel("face-outcome switching fraction")
    axes[1, 1].set_yscale("log"); axes[1, 1].set_ylabel("directional condition number")
    figure.suptitle("Task 0A paired ensemble response" + ("" if complete else " — INCOMPLETE PILOT"))
    figure.savefig(output_dir / "ensemble_response.png", dpi=180)
    plt.close(figure)
    if condition_rows:
        figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
        for level in levels:
            rows = [row for row in condition_rows
                    if row["level"] == level and np.isclose(row["active_threshold"], 1e-4)]
            if not rows:
                continue
            radius = [row["radius_v"] for row in rows]
            axes[0].plot(radius, [row["ensemble_directional_condition"] for row in rows],
                         "o-", label=f"full 2^{level}")
            axes[0].plot(radius, [row["no_switch_directional_condition"] for row in rows],
                         "x--", label=f"without switch component 2^{level}")
            axes[1].errorbar(
                radius, [row["signed_log_condition_change"] for row in rows],
                yerr=[row["signed_log_condition_change_stderr"] for row in rows],
                fmt="o-", label=f"2^{level}")
        for axis in axes:
            axis.set_xscale("log"); axis.grid(alpha=0.25); axis.legend()
            axis.set_xlabel("perturbation radius (V)")
        axes[0].set_yscale("log"); axes[0].set_ylabel("directional condition number")
        axes[0].set_title("Signed switch-component removal")
        axes[1].axhline(0.0, color="black", lw=0.8)
        axes[1].set_ylabel("log(cond full / cond without switches)")
        axes[1].set_title("positive = switches worsen conditioning")
        figure.savefig(output_dir / "condition_decomposition.png", dpi=180)
        plt.close(figure)
    if scaling_rows:
        figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
        for level in levels:
            rows = [row for row in table if row["level"] == level]
            radii = sorted({row["radius_v"] for row in rows})
            fraction = [np.median([
                row["face_switch_fraction"] for row in rows if row["radius_v"] == radius])
                for radius in radii]
            axis.plot(radii, fraction, "o-", label=f"median 2^{level}")
        anchor_radius = max(radii); anchor = fraction[-1]
        reference_radius = np.asarray(sorted(radii))
        axis.plot(reference_radius, anchor * reference_radius / anchor_radius,
                  ":", color="black", label="linear slope = 1")
        axis.set_xscale("log"); axis.set_yscale("log"); axis.grid(alpha=0.25)
        axis.set_xlabel("perturbation radius (V)")
        axis.set_ylabel("paired face-outcome switching fraction")
        median_slope = np.median([row["switch_fraction_loglog_slope"] for row in scaling_rows])
        axis.set_title(f"Switch-set scaling; median fitted slope = {median_slope:.2f}")
        axis.legend()
        figure.savefig(output_dir / "switch_scaling.png", dpi=180)
        plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--checkpoint", required=True)
    run_parser.add_argument("--jacobian", required=True)
    run_parser.add_argument("--output-dir", required=True)
    run_parser.add_argument("--levels", type=int, nargs="+", default=[9, 11])
    run_parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    run_parser.add_argument("--radii", type=float, nargs="+", default=list(RADII))
    run_parser.add_argument("--direction", action="append")
    run_parser.add_argument(
        "--direction-kind", action="append",
        choices=("worst_coordinate", "random", "dominant"))
    run_parser.add_argument("--max-jobs", type=int, default=0)
    run_parser.add_argument("--workers", type=int, default=1)
    run_parser.add_argument("--force", action="store_true")
    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("--output-dir", required=True)
    analyze_parser.add_argument("--allow-incomplete", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    if arguments.command == "run":
        run(arguments)
    else:
        analyze(arguments)
