"""Bounded real-trench sensitivity to sourced Ar+-on-SiO2 kinetic electron emission."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from petch.charged_surface_response_3d import Sobolewski2021ArKineticSEE3D
from petch.charging_coupled_3d import advance_dielectric_charging_3d, current_balance_metrics_3d

import charging_surface_response_preflight_3d as preflight
import charging_task1_physical_time_3d as task1


DEFAULT_STATE = (
    ROOT / "results/charging_task1_3d_refined_15us_level13_audit/final_states.npz")


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _evaluate(charge, method, geometry, poisson, verts, faces, centroids, areas, normals, run):
    source_z = 2.0
    boundary = task1._boundary(source_z * geometry.mesh_length_unit_m)
    proposals = {
        "Ar+": task1._ion_proposal(boundary, 9, 79),
        "electron": task1._electron_proposal(boundary, 9, 83),
    }
    bidirectional = dict(
        forward_log2_samples=11, adjoint_log2_samples=9, n_replicates=4,
        max_forward_log2_samples=11, max_adjoint_log2_samples=11,
        max_face_quadrature_points=7, element_absolute_tolerance=0.02,
        element_relative_tolerance=0.1, face_quadrature_points=3,
        method_hint={"Ar+": method}, require_certification=False)
    response = None
    if run["emission_energy_eV"] is not None:
        response = Sobolewski2021ArKineticSEE3D(
            "plasma_exposed_SiO2", emission_energy_eV=run["emission_energy_eV"],
            emission_energy_evidence=(
                "bounded 1/3/5 eV sensitivity: Huang and Kushner 2026 state ion-induced "
                "secondaries have average energy of a few eV; exact SiO2 spectrum unresolved"),
            angular_log2_samples=run["angular_log2_samples"], angular_seed=113)
    result = advance_dielectric_charging_3d(
        poisson, charge, boundary, verts, faces, areas,
        source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        duration_s=1.0e-15, mesh_length_unit_m=geometry.mesh_length_unit_m,
        n_position=256, seed=79, trajectory_fixed_dt=run["trajectory_fixed_dt"],
        trajectory_max_steps=100000, phase_space_log2_samples=11,
        periodic_lateral=True,
        transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
        face_centroids=centroids, face_gas_normals=normals,
        adjoint_face_quadrature_points=3, adjoint_ray_offset=1.0e-4,
        adjoint_proposals=proposals,
        adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
        bidirectional_options=bidirectional, transport_device="cpu",
        charged_surface_response=response,
        face_material_id=(
            None if response is None
            else np.full(len(faces), "plasma_exposed_SiO2")),
        response_launch_offset=run["response_launch_offset"],
        response_fixed_dt=run["response_fixed_dt"], response_max_bounces=2)
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    regions = preflight._surface_regions(centroids)
    row = dict(run=run["name"], emission_energy_eV=run["emission_energy_eV"],
               angular_log2_samples=run["angular_log2_samples"],
               trajectory_fixed_dt=run["trajectory_fixed_dt"],
               response_fixed_dt=run["response_fixed_dt"],
               response_launch_offset=run["response_launch_offset"])
    for region in ("top", "upper_wall", "lower_wall", "floor", "all"):
        selected = (np.ones(len(faces), dtype=bool) if region == "all" else regions == region)
        positive = float(np.dot(
            result.positive_face_current_density_a_m2[selected], physical_area[selected]))
        negative = float(np.dot(
            result.negative_face_current_density_a_m2[selected], physical_area[selected]))
        row[f"{region}_signed_imbalance"] = (positive - negative) / (positive + negative)
    node = current_balance_metrics_3d(
        result.positive_current_node_a, result.negative_current_node_a)
    row["node_rms_imbalance"] = node.rms_relative_imbalance
    row["node_worst_imbalance"] = node.maximum_relative_imbalance
    row["emitted_rate_s"] = 0.0
    row["landed_rate_s"] = 0.0
    row["escaped_rate_s"] = 0.0
    row["cascade_relative_charge_error"] = 0.0
    if response is not None:
        cascade = result.surface_transfer
        first_transfer = cascade.transfers[0]
        row["emitted_rate_s"] = float(sum(
            np.sum(item.event_rate_s) for item in first_transfer.outgoing))
        row["landed_rate_s"] = float(sum(
            item.landed_rate_s for item in cascade.flights_by_bounce[0]))
        row["escaped_rate_s"] = float(sum(
            item.escaped_rate_s for item in cascade.flights_by_bounce[0]))
        row["cascade_relative_charge_error"] = cascade.relative_charge_balance_error
    return row


def _plot(path, rows):
    import matplotlib.pyplot as plt

    display = [row for row in rows if row["run"] in ("absorber", "energy_1eV", "energy_3eV", "energy_5eV")]
    labels = ["absorber", "1 eV", "3 eV", "5 eV"]
    regions = ("top", "upper_wall", "lower_wall", "floor")
    figure, axes = plt.subplots(1, 3, figsize=(14.4, 4.4), constrained_layout=True)
    x = np.arange(len(regions)); width = 0.19
    for index, (row, label) in enumerate(zip(display, labels)):
        axes[0].bar(
            x + (index - 1.5) * width,
            [row[f"{region}_signed_imbalance"] for region in regions], width, label=label)
    axes[0].axhline(0.08, color="0.35", ls="--", lw=1)
    axes[0].axhline(-0.08, color="0.35", ls="--", lw=1)
    axes[0].set_xticks(x, regions, rotation=25)
    axes[0].set(title="Regional signed current", ylabel="(positive − negative) / total")
    axes[0].legend(fontsize=8)

    axes[1].plot(labels, [row["node_rms_imbalance"] for row in display], "o-", label="RMS")
    axes[1].plot(labels, [row["node_worst_imbalance"] for row in display], "s-", label="worst")
    axes[1].axhline(0.08, color="0.35", ls="--", lw=1, label="0.08 contract")
    axes[1].set(title="Unchanged nodal contract", ylabel="relative imbalance")
    axes[1].legend(fontsize=8)

    emitted = np.array([row["emitted_rate_s"] for row in display[1:]])
    landed = np.array([row["landed_rate_s"] for row in display[1:]])
    escaped = np.array([row["escaped_rate_s"] for row in display[1:]])
    axes[2].bar(labels[1:], landed / emitted, label="landed")
    axes[2].bar(labels[1:], escaped / emitted, bottom=landed / emitted, label="escaped")
    axes[2].set_ylim(0.0, 1.0)
    axes[2].set(title="First-flight SEE routing", ylabel="fraction of emitted rate")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(axis="y", alpha=0.22)
    figure.suptitle("Ar+ kinetic-SEE sensitivity — sourced yield, declared few-eV energy")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results/charging_ion_see_sensitivity_3d")
    args = parser.parse_args()
    geometry, poisson = task1._geometry_and_poisson(0.125)
    verts, faces, centroids, areas = task1.extract_mesh_3d(geometry.phi, geometry.dx)
    normals = task1._surface_gas_normals(verts, faces, centroids, geometry)
    with np.load(args.state) as archive:
        charge = np.asarray(archive["refined_charge_node_c"], dtype=float).copy()
        method = np.asarray(archive["refined_method_hint_Ar+"]).astype("U7")
    runs = [
        dict(name="absorber", emission_energy_eV=None, angular_log2_samples=0,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.005, response_launch_offset=1e-4),
        dict(name="energy_1eV", emission_energy_eV=1.0, angular_log2_samples=4,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.005, response_launch_offset=1e-4),
        dict(name="energy_3eV", emission_energy_eV=3.0, angular_log2_samples=4,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.005, response_launch_offset=1e-4),
        dict(name="energy_5eV", emission_energy_eV=5.0, angular_log2_samples=4,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.005, response_launch_offset=1e-4),
        dict(name="angular_level3", emission_energy_eV=3.0, angular_log2_samples=3,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.005, response_launch_offset=1e-4),
        dict(name="refined_flight", emission_energy_eV=3.0, angular_log2_samples=4,
             trajectory_fixed_dt=0.005, response_fixed_dt=0.0025, response_launch_offset=5e-5),
    ]
    rows = [
        _evaluate(charge, method, geometry, poisson, verts, faces, centroids, areas, normals, run)
        for run in runs]
    config = dict(
        schema="petch.charging.ion_see_sensitivity.v1", state=args.state.name,
        state_sha256=_sha256(args.state), primary_forward_log2_samples=11,
        material_assignment="all feature faces declared plasma_exposed_SiO2 for bounded sensitivity",
        response_provenance=dict(Sobolewski2021ArKineticSEE3D(
            "plasma_exposed_SiO2", 3.0, "bounded literature sensitivity").provenance),
        runs=runs,
        exact_operator=(
            "hard visibility; full nodal-field charged re-impact; no electron-impact SEE, "
            "surface conduction, bulk leakage, photon/metastable/potential emission"))
    encoded = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config_hash = hashlib.sha256(encoded).hexdigest()
    summary = dict(config_hash=config_hash, config=config, rows=rows,
                   conclusion="bounded sensitivity only; convergence gate unchanged")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output_dir / "ion_see_sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    _plot(args.output_dir / "ion_see_sensitivity.png", rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
