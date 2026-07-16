#!/usr/bin/env python3
"""Score one frozen Nozawa charge state with independent exact current ensembles.

This is a diagnostic, not a continuation.  The authoritative face charge is
held fixed while nested scrambled-Sobol ensembles estimate the kinetic current
at increasing ``n_position``.  The output distinguishes a coherent physical
drift from a constant-step stochastic noise floor without re-etching or
re-integrating the preceding charge trajectory.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from petch.charging_coevolution_3d import (
    _face_material_ids,
    _patch_balances,
    _q1_patch_balance_diagnostics,
    _surface_gas_normals,
    integrate_surface_charging_to_saturation_3d,
)
from petch.charging_coupled_3d import current_balance_metrics_3d
from petch.charging_poisson_3d import (
    CompatibleQ1SurfaceChargeProjector3D,
    lump_mixed_surface_density_3d,
)
from petch.feature_step_3d import extract_mesh_3d
from petch.nozawa_replay_3d import make_nozawa_1995_replay_setup


SCHEMA = "petch-nozawa-1995-fixed-charge-ensemble-audit-v1"


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json_value(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if np.isnan(value):
            return "nan"
        if np.isposinf(value):
            return "+inf"
        if np.isneginf(value):
            return "-inf"
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__}")


def _atomic_json(path: Path, payload) -> None:
    encoded = (json.dumps(_json_value(payload), indent=2, sort_keys=True) + "\n").encode()
    with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(encoded)
    os.replace(temporary, path)


def _atomic_npz(path: Path, **arrays) -> None:
    with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", suffix=".npz", delete=False) as stream:
        temporary = Path(stream.name)
        np.savez_compressed(stream, **arrays)
    os.replace(temporary, path)


def _bootstrap_metrics(
        positive_face, negative_face, potential_rate, physical_area,
        patch_groups, patch_scales, *, count=2000, seed=20260716):
    rng = np.random.default_rng(seed)
    b1 = np.empty(count)
    b2_max = np.empty((count, len(patch_scales)))
    b2_rms = np.empty_like(b2_max)
    replicates = len(positive_face)
    for sample in range(count):
        chosen = rng.integers(0, replicates, size=replicates)
        positive = np.mean(positive_face[chosen], axis=0)
        negative = np.mean(negative_face[chosen], axis=0)
        b1[sample] = np.max(np.abs(np.mean(potential_rate[chosen], axis=0)))
        balances = _patch_balances(
            positive, negative, physical_area, patch_groups, patch_scales)
        b2_max[sample] = [
            item.b2_maximum_ion_normalized_imbalance for item in balances]
        b2_rms[sample] = [
            item.b2_rms_ion_normalized_imbalance for item in balances]
    return {
        "b1_potential_rate_v_s_95ci": np.quantile(
            b1, (0.025, 0.975), method="nearest"),
        # Nearest-order bootstrap bounds remain defined when sparse ion sampling
        # produces legitimate +inf B2 values in part or all of the resamples.
        "b2_patch_max_95ci": np.quantile(
            b2_max, (0.025, 0.975), axis=0, method="nearest").T,
        "b2_patch_rms_95ci": np.quantile(
            b2_rms, (0.025, 0.975), axis=0, method="nearest").T,
    }


def _evaluate_level(
        *, condition, sigma, n_position, seed, epochs, timestep_s,
        trajectory_emergency_max_steps, transport_device, output):
    setup = make_nozawa_1995_replay_setup(
        condition, mode="charge_audit", n_position=n_position, seed=seed,
        charging_timestep_s=timestep_s, maximum_charging_steps=1,
        trajectory_emergency_max_steps=trajectory_emergency_max_steps,
        transport_device=transport_device)
    process = setup.process
    geometry = process.geometry
    vertices, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    if sigma.shape != (len(faces),):
        raise ValueError(
            f"checkpoint has {len(sigma)} faces but the replay geometry has {len(faces)}")
    normals = _surface_gas_normals(vertices, faces, centroids, geometry)
    material = _face_material_ids(centroids, geometry)
    poisson = process.charging_system_builder(geometry)
    face_conductor_id = poisson.classify_surface_floating_conductors(
        centroids, normals, grid_origin=process.potential_origin,
        grid_spacing=process.potential_spacing)
    raw_charge = lump_mixed_surface_density_3d(
        poisson, vertices, faces, sigma, face_conductor_id,
        grid_origin=process.potential_origin,
        grid_spacing=process.potential_spacing,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    raw_potential, raw_poisson = poisson.solve(raw_charge)
    projector = CompatibleQ1SurfaceChargeProjector3D.from_mixed_poisson_system(
        poisson, vertices, faces, face_conductor_id,
        grid_origin=process.potential_origin,
        grid_spacing=process.potential_spacing,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    original_face_charge = sigma * physical_area
    compatible_face_charge = projector.project_face_charge(original_face_charge)
    compatible_sigma = compatible_face_charge / physical_area
    solver = dict(process.solver_options)
    options = dict(process.charging_options)
    options.pop("terminal_window_s", None)
    # The state is projected once above. A zero-step current audit cannot create a new charge
    # state, so rebuilding the identical dense mixed projector inside every independent replicate
    # is pure overhead. The kinetic current and Poisson field are exactly the same.
    options.update(
        maximum_steps=0, stop_on_saturation=False,
        compatible_q1_charge_state=False)

    positive_face = []
    negative_face = []
    positive_node = []
    negative_node = []
    potential_rate = []
    sampling_seeds = []
    trajectory_horizons = []
    result_reference = None
    for replicate, epoch in enumerate(epochs):
        result = integrate_surface_charging_to_saturation_3d(
            poisson, compatible_sigma, process.boundary, vertices, faces, areas,
            face_centroids=centroids, face_gas_normals=normals,
            face_material_id=material, source_bounds=process.source_bounds,
            source_z=process.source_z, potential_origin=process.potential_origin,
            potential_spacing=process.potential_spacing,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            mesh_origin_m=geometry.mesh_origin_m, n_position=n_position,
            seed=seed, trajectory_fixed_dt=solver["trajectory_fixed_dt"],
            trajectory_max_steps=solver["trajectory_max_steps"],
            trajectory_adaptive_horizon=solver["trajectory_adaptive_horizon"],
            trajectory_emergency_max_steps=solver[
                "trajectory_emergency_max_steps"],
            periodic_lateral=solver["periodic_lateral"],
            transport_device=solver["transport_device"],
            charged_surface_response=solver.get("charged_surface_response"),
            initial_sampling_epoch=int(epoch), **options)
        step = result.final_step
        record = result.history[-1]
        if result_reference is None:
            result_reference = result
        elif not np.array_equal(result.potential_v, result_reference.potential_v):
            raise RuntimeError("fixed-state audits produced different electrostatic potentials")
        positive_face.append(step.positive_face_current_density_a_m2)
        negative_face.append(step.negative_face_current_density_a_m2)
        positive_node.append(step.positive_current_node_a)
        negative_node.append(step.negative_current_node_a)
        potential_rate.append(
            (step.potential_after_v - step.potential_before_v) / timestep_s)
        sampling_seeds.append(record["sampling_seed"])
        trajectory_horizons.append(record["transport_trajectory_final_max_steps"])
        replicate_directory = output / f"replicate_{replicate:02d}"
        replicate_directory.mkdir(parents=True, exist_ok=False)
        _atomic_npz(
            replicate_directory / "current_audit.npz",
            positive_face_current_density_a_m2=step.positive_face_current_density_a_m2,
            negative_face_current_density_a_m2=step.negative_face_current_density_a_m2,
            positive_current_node_a=step.positive_current_node_a,
            negative_current_node_a=step.negative_current_node_a,
            potential_before_v=step.potential_before_v,
            potential_rate_v_s=potential_rate[-1],
            physical_face_area_m2=areas * geometry.mesh_length_unit_m ** 2,
            sampling_epoch=np.asarray(epoch),
            sampling_seed=np.asarray(record["sampling_seed"]))

    positive_face = np.stack(positive_face)
    negative_face = np.stack(negative_face)
    positive_node = np.stack(positive_node)
    negative_node = np.stack(negative_node)
    potential_rate = np.stack(potential_rate)
    patch_scales = np.asarray(
        [item.patch_scale_m for item in result_reference.patch_balance], dtype=float)
    patch_groups = tuple(item.group for item in result_reference.patch_balance)
    mean_positive_face = np.mean(positive_face, axis=0)
    mean_negative_face = np.mean(negative_face, axis=0)
    balances = _patch_balances(
        mean_positive_face, mean_negative_face, physical_area,
        patch_groups, patch_scales)
    raw_net_face_current_a = (
        mean_positive_face - mean_negative_face) * physical_area
    resolved_net_face_current_a = projector.project_face_charge(
        raw_net_face_current_a)
    patch_null_sensitivity = tuple(
        max(projector.unresolved_linear_functional_fraction(
            (group == label).astype(float)) for label in np.unique(group))
        for group in patch_groups)
    q1_patch = _q1_patch_balance_diagnostics(
        mean_positive_face, mean_negative_face, physical_area,
        patch_groups, patch_scales, resolved_net_face_current_a,
        patch_null_sensitivity)
    reduced_net = np.stack([
        poisson.reduce_charge(positive - negative)
        for positive, negative in zip(positive_node, negative_node)])
    mean_reduced_net = np.mean(reduced_net, axis=0)
    standard_error_reduced_net = np.std(reduced_net, axis=0, ddof=1) / np.sqrt(
        len(reduced_net))
    signal_l2 = float(np.linalg.norm(mean_reduced_net))
    uncertainty_l2 = float(np.linalg.norm(standard_error_reduced_net))
    mean_positive_reduced = np.mean(
        [poisson.reduce_charge(item) for item in positive_node], axis=0)
    mean_negative_reduced = np.mean(
        [poisson.reduce_charge(item) for item in negative_node], axis=0)
    node_metrics = current_balance_metrics_3d(
        mean_positive_reduced, mean_negative_reduced)
    mean_potential_rate = np.mean(potential_rate, axis=0)
    bootstrap = _bootstrap_metrics(
        positive_face, negative_face, potential_rate, physical_area,
        patch_groups, patch_scales)
    conductor_balance = []
    for conductor_id in poisson.floating_conductor_ids:
        selected = face_conductor_id == conductor_id
        positive_total = float(np.sum(
            mean_positive_face[selected] * physical_area[selected]))
        negative_total = float(np.sum(
            mean_negative_face[selected] * physical_area[selected]))
        conductor_balance.append({
            "conductor_id": int(conductor_id),
            "positive_current_a": positive_total,
            "negative_current_a": negative_total,
            "ion_normalized_imbalance": (
                abs(positive_total - negative_total) / positive_total
                if positive_total > 0.0 else float("inf")),
        })
    audit_path = output / "ensemble_current_audit.npz"
    compatible_path = output / "compatible_state.npz"
    _atomic_npz(
        audit_path,
        ensemble_mean_positive_face_current_density_a_m2=mean_positive_face,
        ensemble_mean_negative_face_current_density_a_m2=mean_negative_face,
        ensemble_mean_positive_current_node_a=np.mean(positive_node, axis=0),
        ensemble_mean_negative_current_node_a=np.mean(negative_node, axis=0),
        ensemble_net_reduced_node_current_standard_error_a=(
            standard_error_reduced_net),
        ensemble_mean_potential_rate_v_s=mean_potential_rate,
        replicate_potential_rate_v_s=potential_rate,
        potential_before_v=result_reference.potential_v,
        physical_face_area_m2=physical_area,
        patch_scales_m=patch_scales,
        patch_group_by_scale=np.stack(patch_groups),
        sampling_epochs=np.asarray(epochs),
        sampling_seeds=np.asarray(sampling_seeds))
    _atomic_npz(
        compatible_path,
        sigma_c_per_m2=compatible_sigma,
        charge_node_c=result_reference.charge_node_c,
        potential_v=result_reference.potential_v,
        source_sampling_epoch=np.asarray(epochs[0]),
        q1_face_coupling_rank=np.asarray(
            result_reference.diagnostics["q1_face_coupling_rank"]),
        q1_face_coupling_nullity=np.asarray(
            result_reference.diagnostics["q1_face_coupling_nullity"]))
    summary = {
        "n_position": n_position,
        "replicate_count": len(epochs),
        "sampling_epochs": list(epochs),
        "sampling_seeds": sampling_seeds,
        "potential_min_v": float(np.min(result_reference.potential_v)),
        "potential_max_v": float(np.max(result_reference.potential_v)),
        "mixed_compatible_projection": {
            "initial_unresolved_face_charge_fraction": projector.unresolved_fraction(
                original_face_charge),
            "initial_unresolved_face_charge_l1_c": float(np.sum(np.abs(
                original_face_charge - compatible_face_charge))),
            "q1_face_coupling_rank": projector.rank,
            "q1_face_coupling_nullity": projector.nullity,
            "maximum_raw_to_compatible_potential_change_v": float(np.max(np.abs(
                result_reference.potential_v - raw_potential))),
            "raw_floating_conductor_charge_c": raw_poisson.floating_conductor_charge_c,
            "compatible_floating_conductor_charge_c": (
                result_reference.final_step.poisson_before.floating_conductor_charge_c),
            "compatible_state_artifact": {
                "name": compatible_path.name,
                "sha256": _hash(compatible_path),
            },
        },
        "ensemble_b1_potential_rate_max_v_s": float(
            np.max(np.abs(mean_potential_rate))),
        "ensemble_b1_potential_rate_v_s_95ci": bootstrap[
            "b1_potential_rate_v_s_95ci"],
        "ensemble_patch_b2_max": [
            item.b2_maximum_ion_normalized_imbalance for item in balances],
        "ensemble_patch_b2_max_95ci": bootstrap["b2_patch_max_95ci"],
        "ensemble_patch_b2_rms": [
            item.b2_rms_ion_normalized_imbalance for item in balances],
        "ensemble_patch_b2_rms_95ci": bootstrap["b2_patch_rms_95ci"],
        "ensemble_q1_resolved_patch_b2_max": [
            item["q1_resolved_maximum_ion_normalized_imbalance"]
            for item in q1_patch],
        "ensemble_q1_resolved_patch_b2_rms": [
            item["q1_resolved_rms_ion_normalized_imbalance"]
            for item in q1_patch],
        "ensemble_q1_unresolved_patch_b2_max": [
            item["q1_unresolved_maximum_ion_normalized_imbalance"]
            for item in q1_patch],
        "ensemble_q1_unresolved_face_current_fraction": (
            projector.unresolved_fraction(raw_net_face_current_a)),
        "floating_conductor_current_balance": conductor_balance,
        "ensemble_node_rms_relative_current_imbalance": (
            node_metrics.rms_relative_imbalance),
        "ensemble_node_worst_relative_current_imbalance": (
            node_metrics.maximum_relative_imbalance),
        "net_reduced_node_current_signal_l2_a": signal_l2,
        "net_reduced_node_current_standard_error_l2_a": uncertainty_l2,
        "net_current_signal_to_standard_error_l2": (
            signal_l2 / uncertainty_l2 if uncertainty_l2 > 0.0 else float("inf")),
        "mean_signed_net_charge_rate_a": float(np.sum(mean_reduced_net)),
        "maximum_trajectory_horizon_steps": int(max(trajectory_horizons)),
        "exact_hard_visibility": True,
        "artifact": {
            "name": audit_path.name,
            "sha256": _hash(audit_path),
        },
    }
    _atomic_json(output / "summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--condition", default="fig10_l06s06_05")
    parser.add_argument("--n-position", type=int, nargs="+", default=(4, 8, 16))
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--first-unused-epoch", type=int, required=True)
    parser.add_argument("--seed", type=int, default=1701)
    parser.add_argument("--timestep-s", type=float, default=5e-8)
    parser.add_argument("--trajectory-emergency-max-steps", type=int, default=65536)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.replicates < 2 or args.first_unused_epoch < 0
            or any(item < 4 for item in args.n_position)
            or len(set(args.n_position)) != len(args.n_position)):
        parser.error("need >=2 replicates, unique n_position >=4, and a nonnegative epoch")
    checkpoint = args.checkpoint.resolve()
    with np.load(checkpoint, allow_pickle=False) as data:
        required = {
            "sigma_c_per_m2", "accepted_steps", "physical_time_s",
            "resume_sampling_epoch"}
        if not required.issubset(data.files):
            parser.error("checkpoint is missing required restart fields")
        sigma = np.asarray(data["sigma_c_per_m2"], dtype=float).copy()
        checkpoint_metadata = {
            "accepted_steps": int(data["accepted_steps"]),
            "physical_time_s": float(data["physical_time_s"]),
            "recorded_resume_sampling_epoch": int(data["resume_sampling_epoch"]),
        }
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    epochs = tuple(
        args.first_unused_epoch + replicate for replicate in range(args.replicates))
    levels = []
    for n_position in args.n_position:
        level_output = output / f"n_position_{n_position}"
        level_output.mkdir()
        levels.append(_evaluate_level(
            condition=args.condition, sigma=sigma, n_position=n_position,
            seed=args.seed, epochs=epochs, timestep_s=args.timestep_s,
            trajectory_emergency_max_steps=args.trajectory_emergency_max_steps,
            transport_device=args.transport_device, output=level_output))
    summary = {
        "schema": SCHEMA,
        "condition_id": args.condition,
        "checkpoint": {
            "name": checkpoint.name,
            "sha256": _hash(checkpoint),
            **checkpoint_metadata,
        },
        "first_unused_epoch_used": args.first_unused_epoch,
        "nested_scramble_epochs_across_sample_levels": True,
        "timestep_s_used_only_to_report_instantaneous_dvdt": args.timestep_s,
        "levels": levels,
        "decision_rule": (
            "A coherent drift requires sample-level stability and a net-current "
            "signal-to-standard-error ratio comfortably above one; otherwise the "
            "state is estimator-noise limited and must not be advanced by another "
            "low-sample constant-step march."),
    }
    _atomic_json(output / "summary.json", summary)
    print(json.dumps(_json_value(summary), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
