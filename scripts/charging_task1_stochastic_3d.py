#!/usr/bin/env python3
"""Fresh-scramble physical-time charging cross-check with a separately frozen method map."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(variable, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import (  # noqa: E402
    _boundary, _depth_charge, _electron_proposal, _geometry_and_poisson, _ion_proposal,
)
from petch.charging_coupled_3d import advance_dielectric_charging_3d  # noqa: E402
from petch.feature_step_3d import _surface_gas_normals  # noqa: E402
from petch.threed import extract_mesh_3d  # noqa: E402


def _metrics(positive, negative):
    total = positive + negative
    active = total > max(1e-15 * float(np.max(total)), 1e-300)
    relative = np.abs(positive[active] - negative[active]) / total[active]
    return float(np.sqrt(np.mean(relative ** 2))), float(np.max(relative))


def _rolling_mean(value, window):
    value = np.asarray(value, dtype=float)
    if value.size < window:
        return np.full_like(value, np.nan)
    result = np.full_like(value, np.nan)
    result[window - 1:] = np.convolve(value, np.ones(window) / window, mode="valid")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/charging_task1_3d_stochastic")
    parser.add_argument("--timestep-s", type=float, default=1.25e-7)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--base-seed", type=int, default=1009)
    parser.add_argument("--forward-level", type=int, default=10)
    parser.add_argument("--ion-proposal-level", type=int, default=8)
    parser.add_argument("--electron-proposal-level", type=int, default=9)
    parser.add_argument("--trajectory-dt", type=float, default=0.005)
    parser.add_argument("--trajectory-max-steps", type=int, default=50000)
    args = parser.parse_args()
    if args.steps <= 0 or args.timestep_s <= 0.0:
        parser.error("positive steps and timestep are required")

    archive_bytes = args.initial_state.read_bytes()
    archive_sha256 = hashlib.sha256(archive_bytes).hexdigest()
    with np.load(args.initial_state) as archived:
        charge = np.asarray(archived["refined_charge_node_c"], dtype=float).copy()
        method_hint = np.asarray(archived["refined_method_hint_Ar+"]).astype("U7")

    geometry, poisson = _geometry_and_poisson()
    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(verts, faces, centroids, geometry)
    source_z = 2.0
    boundary = _boundary(source_z * geometry.mesh_length_unit_m)
    common = dict(
        poisson_system=poisson, boundary=boundary, verts=verts, faces=faces, areas=areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m, n_position=256,
        trajectory_fixed_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        phase_space_log2_samples=args.forward_level, periodic_lateral=True,
        transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-4,
        adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
        transport_device="cpu")
    config = dict(
        model="exact_hard_visibility_3d_fresh_scramble_physical_time",
        initial_state_sha256=archive_sha256, timestep_s=args.timestep_s,
        steps=args.steps, base_seed=args.base_seed, seed_stride=104729,
        forward_level=args.forward_level, ion_proposal_level=args.ion_proposal_level,
        electron_proposal_level=args.electron_proposal_level,
        method_map_source="separate certified pilot archive",
        trajectory_dt=args.trajectory_dt, trajectory_max_steps=args.trajectory_max_steps)
    config_hash = hashlib.sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    history = []
    charge_history = [charge.copy()]
    maximum_conservation = 0.0
    absolute_throughput = 0.0
    for step in range(args.steps):
        seed = args.base_seed + 104729 * step
        proposals = {
            "Ar+": _ion_proposal(boundary, args.ion_proposal_level, seed + 17),
            "electron": _electron_proposal(
                boundary, args.electron_proposal_level, seed + 31),
        }
        result = advance_dielectric_charging_3d(
            charge_node_c=charge, duration_s=args.timestep_s, seed=seed,
            adjoint_proposals=proposals,
            bidirectional_options=dict(
                forward_log2_samples=args.forward_level,
                adjoint_log2_samples=args.ion_proposal_level,
                n_replicates=4, method_hint={"Ar+": method_hint},
                require_certification=False,
                element_absolute_tolerance=0.02, element_relative_tolerance=0.1,
                face_quadrature_points=3),
            **common)
        node_rms, node_max = _metrics(
            result.positive_current_node_a, result.negative_current_node_a)
        face_rms, face_max = _metrics(
            result.positive_face_current_density_a_m2,
            result.negative_face_current_density_a_m2)
        history.append(dict(
            step=step, physical_time_s=step * args.timestep_s, seed=seed,
            rms_relative_current_imbalance_node=node_rms,
            max_relative_current_imbalance_node=node_max,
            rms_relative_current_imbalance_face=face_rms,
            max_relative_current_imbalance_face=face_max,
            minimum_potential_v=float(np.min(result.potential_before_v)),
            maximum_potential_v=float(np.max(result.potential_before_v)),
            total_stored_charge_c=float(np.sum(charge))))
        maximum_conservation = max(
            maximum_conservation,
            abs(result.diagnostics["charge_conservation_residual_c"]))
        absolute_throughput += abs(result.diagnostics["incident_charge_c"])
        charge = result.charge_node_c.copy()
        charge_history.append(charge.copy())
        if step % 25 == 0 or step + 1 == args.steps:
            print(
                f"step={step:4d} rms={node_rms:.6f} max={node_max:.6f} "
                f"V=[{result.potential_before_v.min():.3f},{result.potential_before_v.max():.3f}]",
                flush=True)

    node_rms = np.array([item["rms_relative_current_imbalance_node"] for item in history])
    node_max = np.array([item["max_relative_current_imbalance_node"] for item in history])
    tail = min(args.steps, max(20, args.steps // 5))
    tail_rms_mean = float(np.mean(node_rms[-tail:]))
    tail_rms_stderr = float(np.std(node_rms[-tail:], ddof=1) / np.sqrt(tail))
    tail_max_mean = float(np.mean(node_max[-tail:]))
    tail_max_stderr = float(np.std(node_max[-tail:], ddof=1) / np.sqrt(tail))
    region_charge = _depth_charge(charge, geometry.dx)
    summary = dict(
        config_hash=config_hash, config=config,
        exact_operator_statement="hard visibility; no smoothing or analytic current replacement",
        frozen_method_counts={name: int(np.sum(method_hint == name))
                              for name in ("forward", "adjoint")},
        final=dict(history[-1]), tail_steps=tail,
        tail_statistics=dict(
            rms_mean=tail_rms_mean, rms_stderr=tail_rms_stderr,
            worst_node_mean=tail_max_mean, worst_node_stderr=tail_max_stderr),
        charge_regions_c=region_charge,
        maximum_step_charge_conservation_residual_c=maximum_conservation,
        total_absolute_charge_throughput_c=absolute_throughput,
        conservation_relative_to_throughput=maximum_conservation / max(absolute_throughput, 1e-300),
        statistically_inside_contract=bool(
            tail_max_mean + 2.0 * tail_max_stderr <= 0.08),
        conclusion="stochastic trajectory requires an independent high-sample final-state audit")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "transient_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=history[0].keys())
        writer.writeheader(); writer.writerows(history)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savez_compressed(
        args.output_dir / "trajectory.npz", charge_history_node_c=np.stack(charge_history),
        final_charge_node_c=charge, method_hint_Ar=method_hint,
        vertices=verts, faces=faces, centroids=centroids, areas=areas)

    time_us = np.arange(args.steps) * args.timestep_s * 1e6
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].plot(time_us, node_rms, alpha=0.35, label="fresh-scramble RMS")
    axes[0].plot(time_us, _rolling_mean(node_rms, tail), label=f"{tail}-step mean")
    axes[0].axhline(0.08, color="0.4", linestyle="--", label="node contract")
    axes[0].set(xlabel="physical time (µs)", ylabel="relative imbalance",
                title="Stochastic nodal RMS")
    axes[0].legend(fontsize=8); axes[0].grid(alpha=0.25)
    axes[1].plot(time_us, node_max, alpha=0.35, label="fresh-scramble worst node")
    axes[1].plot(time_us, _rolling_mean(node_max, tail), label=f"{tail}-step mean")
    axes[1].axhline(0.08, color="0.4", linestyle="--", label="node contract")
    axes[1].set(xlabel="physical time (µs)", ylabel="relative imbalance",
                title="Stochastic worst node")
    axes[1].legend(fontsize=8); axes[1].grid(alpha=0.25)
    fig.savefig(args.output_dir / "stochastic_transient.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
