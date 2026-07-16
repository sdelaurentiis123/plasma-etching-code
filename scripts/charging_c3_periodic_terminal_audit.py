#!/usr/bin/env python3
"""Independently audit one bounded periodic C3 terminal-window trajectory.

This script separates two questions that must never be conflated:

* ``integrity_pass`` verifies topology, provenance, conservative time integration, Poisson
  replay, fresh-scramble sequencing, and artifact closure.
* ``contract_converged`` reports the unchanged CCA-R2 B1/B2 gates.

A trajectory can therefore be a completely valid physical calculation without yet being a
converged charging state.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import _geometry_and_poisson  # noqa: E402
from petch.charging_coevolution_3d import (  # noqa: E402
    _patch_balances,
    _q1_patch_balance_diagnostics,
)
from petch.charging_coupled_3d import current_balance_metrics_3d  # noqa: E402
from petch.charging_poisson_3d import (  # noqa: E402
    CompatibleQ1SurfaceChargeProjector3D,
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json_value(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__}")


def _relative_l1(error, reference) -> float:
    denominator = max(float(np.sum(np.abs(reference))), np.finfo(float).tiny)
    return float(np.sum(np.abs(error)) / denominator)


def _block_rows(history, block_size):
    output = []
    for start in range(0, len(history), block_size):
        block = history[start:start + block_size]
        output.append(dict(
            first_evaluation=int(block[0]["evaluation"]),
            last_evaluation=int(block[-1]["evaluation"]),
            physical_time_us=float(block[-1]["physical_time_s"] * 1e6),
            node_rms_mean=float(np.mean([
                item["rms_relative_current_imbalance_node"] for item in block])),
            node_rms_last=float(block[-1]["rms_relative_current_imbalance_node"]),
            node_worst_mean=float(np.mean([
                item["max_relative_current_imbalance_node"] for item in block])),
            potential_rate_median_v_s=float(np.median([
                item["potential_rate_max_v_s"] for item in block])),
            raw_patch_b2_median=float(np.median([
                item["maximum_patch_relative_imbalance"] for item in block])),
            q1_resolved_patch_b2_median=float(np.median([
                max(item["patch_q1_resolved_max_ion_normalized_imbalance"])
                for item in block])),
            unresolved_current_fraction_median=float(np.median([
                item["unresolved_face_current_fraction"] for item in block]))))
    return output


def _plot(history, result, output):
    time_us = np.asarray([item["physical_time_s"] for item in history]) * 1e6
    node_rms = np.asarray([
        item["rms_relative_current_imbalance_node"] for item in history])
    node_worst = np.asarray([
        item["max_relative_current_imbalance_node"] for item in history])
    rate = np.asarray([item["potential_rate_max_v_s"] for item in history])
    raw_b2 = np.asarray([item["maximum_patch_relative_imbalance"] for item in history])
    resolved_b2 = np.asarray([
        max(item["patch_q1_resolved_max_ion_normalized_imbalance"])
        for item in history])
    unresolved = np.asarray([
        item["unresolved_face_current_fraction"] for item in history])

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), constrained_layout=True)
    axis = axes[0, 0]
    axis.plot(time_us, node_rms, alpha=0.55, label="node RMS")
    axis.plot(time_us, node_worst, alpha=0.55, label="worst node")
    axis.axhline(0.08, color="black", linestyle="--", linewidth=1, label="0.08")
    axis.set_yscale("log")
    axis.set_title("Independent Q1-node current diagnostic")
    axis.set_xlabel("corrected physical time (µs)")
    axis.set_ylabel("relative imbalance")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[0, 1]
    axis.plot(time_us, rate, alpha=0.5, label="instantaneous |dV/dt|")
    axis.axhline(1e3, color="black", linestyle="--", linewidth=1, label="B1 gate")
    axis.scatter(
        [time_us[-1]], [result["final_potential_rate_max_v_s"]], marker="D", s=50,
        label="50 µs integrated B1")
    axis.set_yscale("log")
    axis.set_title("Voltage motion: noisy samples vs integrated drift")
    axis.set_xlabel("corrected physical time (µs)")
    axis.set_ylabel("V/s")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 0]
    axis.plot(time_us, raw_b2, alpha=0.5, label="raw face-patch B2")
    axis.plot(time_us, resolved_b2, alpha=0.65, label="Q1-resolved diagnostic")
    axis.axhline(0.08, color="black", linestyle="--", linewidth=1, label="B2 gate")
    axis.scatter(
        [time_us[-1]], [max(result["patch_b2_max_ion_normalized"])], marker="D", s=50,
        label="50 µs raw B2")
    axis.scatter(
        [time_us[-1]], [max(result["q1_resolved_patch_b2_max"])], marker="s", s=45,
        label="50 µs Q1-resolved")
    axis.set_yscale("log")
    axis.set_title("Patch balance separates raw and field-visible modes")
    axis.set_xlabel("corrected physical time (µs)")
    axis.set_ylabel("ion-normalized imbalance")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)

    axis = axes[1, 1]
    axis.plot(time_us, 100.0 * unresolved, alpha=0.65)
    axis.set_ylim(0.0, 100.0)
    axis.set_title("Instantaneous current in the Q1-null face subspace")
    axis.set_xlabel("corrected physical time (µs)")
    axis.set_ylabel("unresolved norm (%)")
    axis.grid(alpha=0.25)

    figure.suptitle("Corrected periodic C3: valid trajectory, saturation still open")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial-checkpoint", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--grid-dx-um", type=float, default=0.25)
    parser.add_argument("--block-size", type=int, default=50)
    args = parser.parse_args()
    initial_path = args.initial_checkpoint.resolve()
    run = args.run_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.block_size <= 0:
        parser.error("--block-size must be positive")

    required = ("summary.json", "face_checkpoint.npz", "current_audit.npz")
    if any(not (run / name).exists() for name in required):
        parser.error("run directory is missing a required terminal artifact")
    summary = json.loads((run / "summary.json").read_text())
    result = summary["result"]
    history = summary["history"]
    config = summary["config"]
    config_hash = sha256(json.dumps(
        config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    artifact_hashes = {
        name: _hash(run / record["name"]) == record["sha256"]
        for name, record in summary["artifacts"].items()}
    with np.load(initial_path) as initial, np.load(
            run / "face_checkpoint.npz") as final, np.load(
            run / "current_audit.npz") as current:
        initial_face_charge = np.asarray(initial["face_charge_c"], dtype=float)
        final_face_charge = np.asarray(final["face_charge_c"], dtype=float)
        final_sigma = np.asarray(final["sigma_c_per_m2"], dtype=float)
        final_node = np.asarray(final["charge_node_c"], dtype=float)
        final_potential = np.asarray(final["potential_v"], dtype=float)
        vertices = np.asarray(final["vertices"])
        faces = np.asarray(final["faces"], dtype=int)
        areas = np.asarray(final["areas"], dtype=float)
        terminal_positive_density = np.asarray(
            current["terminal_window_positive_face_current_density_a_m2"], dtype=float)
        terminal_negative_density = np.asarray(
            current["terminal_window_negative_face_current_density_a_m2"], dtype=float)
        positive_node = np.asarray(current["positive_current_node_a"], dtype=float)
        negative_node = np.asarray(current["negative_current_node_a"], dtype=float)
        patch_groups = np.asarray(current["patch_group_by_scale"], dtype=int)
        patch_scales = np.asarray(current["patch_scales_m"], dtype=float)
        terminal_window_s = float(np.asarray(current["terminal_window_s"]).item())

    geometry, poisson = _geometry_and_poisson(args.grid_dx_um)
    projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        poisson, vertices, faces, grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    terminal_positive = terminal_positive_density * physical_area
    terminal_negative = terminal_negative_density * physical_area
    terminal_net = terminal_positive - terminal_negative
    resolved_terminal_net = projector.project_face_charge(terminal_net)
    predicted_delta = resolved_terminal_net * terminal_window_s
    measured_delta = final_face_charge - initial_face_charge
    transported_scale = max(
        float(np.sum(np.abs(terminal_positive) + np.abs(terminal_negative)))
        * terminal_window_s, np.finfo(float).tiny)

    initial_reduced = projector.node_charge_from_face_charge(initial_face_charge)
    final_reduced = projector.node_charge_from_face_charge(final_face_charge)
    predicted_node_delta = projector.node_charge_from_face_charge(
        terminal_net) * terminal_window_s
    recomputed_node = poisson.canonicalize_reduced_charge(final_reduced)
    recomputed_potential, poisson_diagnostics = poisson.solve(recomputed_node)
    manual_b1 = float(np.max(np.abs(
        final_potential - np.asarray(np.load(initial_path)["potential_v"], dtype=float)))
        / terminal_window_s)

    final_node_metrics = current_balance_metrics_3d(
        poisson.reduce_charge(positive_node), poisson.reduce_charge(negative_node))
    patch = _patch_balances(
        terminal_positive_density, terminal_negative_density,
        physical_area, tuple(patch_groups), tuple(patch_scales))
    sensitivities = tuple(
        max(projector.unresolved_linear_functional_fraction(
            (groups == group).astype(float)) for group in np.unique(groups))
        for groups in patch_groups)
    q1_patch = _q1_patch_balance_diagnostics(
        terminal_positive_density, terminal_negative_density,
        physical_area, tuple(patch_groups), tuple(patch_scales),
        resolved_terminal_net, sensitivities)

    timestep = float(config["timestep_s"])
    epochs = np.asarray([item["sampling_epoch"] for item in history], dtype=int)
    expected_epochs = int(config["initial_sampling_epoch"]) + np.arange(len(history))
    physical_times = np.asarray([item["physical_time_s"] for item in history])
    expected_times = np.arange(len(history)) * timestep
    manual_raw_b2 = [item.b2_maximum_ion_normalized_imbalance for item in patch]
    manual_q1_b2 = [
        item["q1_resolved_maximum_ion_normalized_imbalance"] for item in q1_patch]
    integrity_checks = dict(
        config_hash=config_hash == summary["config_hash"],
        artifact_hashes=all(artifact_hashes.values()),
        corrected_topology=(
            config["particle_periodic_lateral"]
            and tuple(config["poisson_periodic_axes"]) == (0, 1)
            and tuple(poisson.periodic_axes) == (0, 1)),
        fresh_unique_sequential_epochs=(
            config["scramble_mode"] == "fresh"
            and np.array_equal(epochs, expected_epochs)
            and len(np.unique(epochs)) == len(epochs)),
        exact_step_times=np.allclose(
            physical_times, expected_times, rtol=2e-13, atol=2e-18),
        accepted_without_rejection=(
            result["accepted_steps"] == config["maximum_steps"]
            and result["rejected_steps"] == 0
            and len(history) == result["accepted_steps"] + 1),
        periodic_potential_seam=(
            np.array_equal(final_potential[0], final_potential[-1])
            and np.array_equal(final_potential[:, 0], final_potential[:, -1])),
        periodic_charge_seam=(
            np.array_equal(final_node[0], final_node[-1])
            and np.array_equal(final_node[:, 0], final_node[:, -1])),
        compatible_face_state=projector.unresolved_fraction(final_face_charge) < 1e-12,
        face_inventory=np.allclose(
            final_face_charge, final_sigma * physical_area,
            rtol=2e-13, atol=2e-31),
        terminal_face_update=_relative_l1(
            measured_delta - predicted_delta, measured_delta) < 5e-10,
        terminal_node_update=_relative_l1(
            (final_reduced - initial_reduced) - predicted_node_delta,
            final_reduced - initial_reduced) < 5e-10,
        terminal_global_charge=(
            abs(float(np.sum(measured_delta) - np.sum(terminal_net) * terminal_window_s))
            / transported_scale < 5e-12),
        checkpoint_node_replay=_relative_l1(
            recomputed_node - final_node, recomputed_node) < 5e-13,
        checkpoint_potential_replay=(
            np.max(np.abs(recomputed_potential - final_potential)) < 1e-10),
        poisson_free_residual=poisson_diagnostics.max_abs_free_charge_residual_c < 1e-27,
        poisson_global_balance=abs(poisson_diagnostics.charge_balance_c) < 1e-26,
        summary_b1=np.isclose(
            manual_b1, result["final_potential_rate_max_v_s"], rtol=2e-12),
        summary_raw_b2=np.allclose(
            manual_raw_b2, result["patch_b2_max_ion_normalized"],
            rtol=2e-12, atol=2e-14),
        summary_q1_b2=np.allclose(
            manual_q1_b2, result["q1_resolved_patch_b2_max"],
            rtol=2e-12, atol=2e-14),
        summary_node_metrics=(
            np.isclose(
                final_node_metrics.rms_relative_imbalance,
                result["retained_node_rms_relative_current_imbalance"], rtol=2e-12)
            and np.isclose(
                final_node_metrics.maximum_relative_imbalance,
                result["retained_node_max_relative_current_imbalance"], rtol=2e-12)),
        per_step_charge_ledger=(
            max(item["charge_conservation_relative_error"] for item in history) < 5e-13),
        per_step_surface_transfer_ledger=(
            max(item["surface_transfer_relative_charge_balance_error"]
                for item in history) < 5e-13),
        hard_visibility_tail_bound=(
            max(item["response_tail_closure_l1_current_error_bound_relative"]
                for item in history) <= 2.1e-10))
    integrity_pass = all(bool(value) for value in integrity_checks.values())
    contract_converged = bool(
        result["final_potential_rate_max_v_s"] <= config["potential_rate_tolerance_v_s"]
        and all(value <= 0.08 for value in result["patch_b2_max_ion_normalized"]))

    block_rows = _block_rows(history, args.block_size)
    csv_path = output / "periodic_terminal_blocks.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(block_rows[0]))
        writer.writeheader()
        writer.writerows(block_rows)
    plot_path = output / "periodic_terminal_audit.png"
    _plot(history, result, plot_path)

    artifact = dict(
        schema="petch.charging.c3.periodic-terminal-audit.v1",
        status=(
            "integrity verified; saturation gates remain open"
            if integrity_pass and not contract_converged else
            "integrity verified; saturation gates pass"
            if integrity_pass else "integrity audit failed"),
        integrity_pass=integrity_pass,
        contract_converged=contract_converged,
        integrity_checks=integrity_checks,
        measured=dict(
            corrected_physical_time_s=result["physical_time_s"],
            wall_clock_s=summary["run_manifest"]["wall_clock_s"],
            b1_terminal_potential_rate_v_s=result["final_potential_rate_max_v_s"],
            b1_gate_v_s=config["potential_rate_tolerance_v_s"],
            raw_patch_b2_max=result["patch_b2_max_ion_normalized"],
            q1_resolved_patch_b2_max=result["q1_resolved_patch_b2_max"],
            b2_gate=0.08,
            retained_node_rms=result["retained_node_rms_relative_current_imbalance"],
            retained_node_worst=result["retained_node_max_relative_current_imbalance"],
            maximum_charge_conservation_relative_error=result[
                "maximum_charge_conservation_relative_error"],
            maximum_surface_transfer_relative_error=result[
                "maximum_surface_transfer_relative_charge_balance_error"],
            maximum_replay_fraction=result["maximum_transport_lineage_replay_fraction"],
            maximum_bounce_extensions=result[
                "maximum_response_bounce_budget_extension_count"],
            maximum_horizon_extensions=result[
                "maximum_transport_trajectory_horizon_extension_count"],
            face_terminal_update_relative_l1=_relative_l1(
                measured_delta - predicted_delta, measured_delta),
            node_terminal_update_relative_l1=_relative_l1(
                (final_reduced - initial_reduced) - predicted_node_delta,
                final_reduced - initial_reduced),
            global_terminal_update_error_c=float(
                np.sum(measured_delta) - np.sum(terminal_net) * terminal_window_s),
            endpoint_poisson_potential_max_error_v=float(
                np.max(np.abs(recomputed_potential - final_potential))),
            endpoint_face_q1_null_fraction=projector.unresolved_fraction(
                final_face_charge),
            endpoint_x_seam_max_v=float(np.max(np.abs(
                final_potential[0] - final_potential[-1]))),
            endpoint_y_seam_max_v=float(np.max(np.abs(
                final_potential[:, 0] - final_potential[:, -1])))),
        interpretation=dict(
            convergence=(
                "The corrected physical trajectory is valid but is not stationary under the "
                "unchanged CCA-R2 contract."),
            discretization=(
                "Raw patch B2 sees a large exact P0-face/Q1-field null component. The separately "
                "reported Q1-resolved diagnostic is smaller but also remains above 0.08. R3 is "
                "still change control, not an active replacement gate."),
            next_step=(
                "Do not credit pre-repair time. If more dynamics are run, continue only from the "
                "corrected checkpoint in another explicit bounded window and require timestep, "
                "sample, and grid refinement before any convergence claim.")),
        provenance=dict(
            initial_checkpoint=dict(name=initial_path.name, sha256=_hash(initial_path)),
            run_summary=dict(name="summary.json", sha256=_hash(run / "summary.json")),
            artifact_hash_matches=artifact_hashes,
            run_source_sha256=summary["run_manifest"]["source_sha256"],
            audit_script_sha256=_hash(Path(__file__).resolve())),
        artifacts=dict(
            blocks_csv=dict(name=csv_path.name, sha256=_hash(csv_path)),
            plot=dict(name=plot_path.name, sha256=_hash(plot_path))))
    audit_path = output / "periodic_terminal_audit.json"
    audit_path.write_text(json.dumps(_json_value(artifact), indent=2) + "\n")
    print(json.dumps(_json_value({
        "status": artifact["status"],
        "integrity_pass": integrity_pass,
        "contract_converged": contract_converged,
        "measured": artifact["measured"],
        "artifacts": artifact["artifacts"]}), indent=2))
    return 0 if integrity_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
