#!/usr/bin/env python3
"""Safeguarded pseudo-transient continuation using only the conservative current direction."""
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


def _metrics(step):
    total = step.positive_current_node_a + step.negative_current_node_a
    active = total > max(1e-15 * float(np.max(total)), 1e-300)
    relative = np.abs(
        step.positive_current_node_a[active] - step.negative_current_node_a[active]) / total[active]
    return float(np.sqrt(np.mean(relative ** 2))), float(np.max(relative))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/charging_task2a_ptc_3d")
    parser.add_argument("--initial-timestep-s", type=float, default=1.25e-7)
    parser.add_argument("--maximum-timestep-s", type=float, default=5e-7)
    parser.add_argument("--minimum-timestep-s", type=float, default=1e-11)
    parser.add_argument("--accepted-steps", type=int, default=80)
    parser.add_argument("--growth", type=float, default=1.25)
    parser.add_argument("--allowed-growth", type=float, default=0.005)
    args = parser.parse_args()

    payload = args.initial_state.read_bytes()
    initial_sha256 = hashlib.sha256(payload).hexdigest()
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
        mesh_length_unit_m=geometry.mesh_length_unit_m, n_position=256, seed=79,
        trajectory_fixed_dt=0.005, trajectory_max_steps=50000,
        phase_space_log2_samples=10, periodic_lateral=True,
        transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-4,
        adjoint_proposals={
            "Ar+": _ion_proposal(boundary, 8, 79),
            "electron": _electron_proposal(boundary, 10, 83)},
        adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
        bidirectional_options=dict(
            forward_log2_samples=10, adjoint_log2_samples=8, n_replicates=4,
            method_hint={"Ar+": method_hint}, require_certification=False,
            element_absolute_tolerance=0.02, element_relative_tolerance=0.1,
            face_quadrature_points=3), transport_device="cpu")
    dt = float(args.initial_timestep_s)
    accepted = 0
    rejected = 0
    pseudo_time = 0.0
    history = []
    charge_history = [charge.copy()]
    best_charge = charge.copy()
    best_maximum = float("inf")
    best_rms = float("inf")

    current = advance_dielectric_charging_3d(
        charge_node_c=charge, duration_s=dt, **common)
    current_rms, current_max = _metrics(current)
    while accepted < args.accepted_steps and dt >= args.minimum_timestep_s:
        candidate_charge = current.charge_node_c
        trial = advance_dielectric_charging_3d(
            charge_node_c=candidate_charge, duration_s=dt, **common)
        trial_rms, trial_max = _metrics(trial)
        allowed = 1.0 + float(args.allowed_growth)
        accept = trial_rms <= allowed * current_rms and trial_max <= allowed * current_max
        if accept:
            charge = candidate_charge.copy()
            pseudo_time += dt
            accepted += 1
            history.append(dict(
                accepted_step=accepted, rejected_steps=rejected,
                pseudo_time_s=pseudo_time, timestep_s=dt,
                rms_relative_current_imbalance_node=trial_rms,
                max_relative_current_imbalance_node=trial_max,
                minimum_potential_v=float(np.min(trial.potential_before_v)),
                maximum_potential_v=float(np.max(trial.potential_before_v)),
                charge_conservation_residual_c=float(
                    current.diagnostics["charge_conservation_residual_c"])))
            charge_history.append(charge.copy())
            if (trial_max < best_maximum
                    or (np.isclose(trial_max, best_maximum) and trial_rms < best_rms)):
                best_charge = charge.copy(); best_maximum = trial_max; best_rms = trial_rms
            if trial_max < 0.9 * current_max:
                dt = min(float(args.maximum_timestep_s), dt * float(args.growth))
            current = advance_dielectric_charging_3d(
                charge_node_c=charge, duration_s=dt, **common)
            current_rms, current_max = _metrics(current)
            if current_max <= 0.08:
                break
        else:
            rejected += 1
            dt *= 0.5
            current = advance_dielectric_charging_3d(
                charge_node_c=charge, duration_s=dt, **common)
            current_rms, current_max = _metrics(current)
        if (accepted + rejected) % 10 == 0:
            print(
                f"accepted={accepted} rejected={rejected} dt={dt:.3e} "
                f"rms={current_rms:.6f} max={current_max:.6f}", flush=True)

    config = dict(
        model="hard_visibility_current_direction_ptc", initial_state_sha256=initial_sha256,
        initial_timestep_s=args.initial_timestep_s,
        maximum_timestep_s=args.maximum_timestep_s,
        minimum_timestep_s=args.minimum_timestep_s, accepted_step_limit=args.accepted_steps,
        growth=args.growth, allowed_growth=args.allowed_growth,
        method_map_source="separate certified pilot archive")
    config_hash = hashlib.sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    summary = dict(
        config_hash=config_hash, config=config, accepted_steps=accepted, rejected_steps=rejected,
        final=(None if not history else history[-1]),
        best=dict(rms=best_rms, worst_node=best_maximum),
        converged=bool(best_maximum <= 0.08),
        exhausted_minimum_timestep=bool(dt < args.minimum_timestep_s),
        exact_operator_statement="hard visibility; unchanged kinetic residual and current direction",
        best_charge_regions_c=_depth_charge(best_charge, geometry.dx),
        conclusion="PTC candidate requires an independent high-sample exact-map audit")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if history:
        with (args.output_dir / "history.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=history[0].keys())
            writer.writeheader(); writer.writerows(history)
    np.savez_compressed(
        args.output_dir / "states.npz", best_charge_node_c=best_charge,
        final_charge_node_c=charge, charge_history_node_c=np.stack(charge_history),
        method_hint_Ar=method_hint)
    if history:
        step = np.array([item["accepted_step"] for item in history])
        fig, axis = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
        axis.plot(step, [item["rms_relative_current_imbalance_node"] for item in history],
                  label="RMS")
        axis.plot(step, [item["max_relative_current_imbalance_node"] for item in history],
                  label="worst node")
        axis.axhline(0.08, color="0.4", linestyle="--", label="node contract")
        axis.set(xlabel="accepted pseudo-step", ylabel="relative imbalance",
                 title="Safeguarded current-direction PTC")
        axis.grid(alpha=0.25); axis.legend()
        fig.savefig(args.output_dir / "ptc_history.png", dpi=180)
        plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
