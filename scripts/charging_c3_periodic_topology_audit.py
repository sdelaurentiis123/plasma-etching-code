#!/usr/bin/env python3
"""Audit and repair the C3 lateral particle/Poisson topology mismatch.

The input checkpoint may have been evolved with periodically wrapped trajectories but a
nonperiodic Q1 field.  This script never credits that history as corrected physical time.  It
constructs the intended periodic Poisson operator, projects the saved face inventory into that
operator's exact compatible charge space, and emits a warm *proposal* for bounded cold/warm
verification under the corrected operator.
"""
from __future__ import annotations

import argparse
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
from petch.charging_poisson_3d import (  # noqa: E402
    CompatibleQ1SurfaceChargeProjector3D,
    NodalPoissonSystem3D,
    lump_triangle_sheet_charge_3d,
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


def _seam_metrics(potential):
    x = np.asarray(potential[0, :, :] - potential[-1, :, :], dtype=float)
    y = np.asarray(potential[:, 0, :] - potential[:, -1, :], dtype=float)
    return dict(
        x_max_abs_v=float(np.max(np.abs(x))),
        x_rms_v=float(np.sqrt(np.mean(x ** 2))),
        y_max_abs_v=float(np.max(np.abs(y))),
        y_rms_v=float(np.sqrt(np.mean(y ** 2))),
        x_by_z_max_abs_v=np.max(np.abs(x), axis=0),
        y_by_z_max_abs_v=np.max(np.abs(y), axis=0))


def _plot(archived, corrected, output):
    old = _seam_metrics(archived)
    new = _seam_metrics(corrected)
    z = np.arange(archived.shape[2])
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), constrained_layout=True)
    axes[0].plot(z, old["x_by_z_max_abs_v"], marker="o", label="x seam")
    axes[0].plot(z, old["y_by_z_max_abs_v"], marker="s", label="y seam")
    axes[0].set_title("Archived mismatched field")
    axes[0].set_xlabel("z node")
    axes[0].set_ylabel("endpoint |delta V| (V)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)

    axes[1].plot(z, new["x_by_z_max_abs_v"], marker="o", label="x seam")
    axes[1].plot(z, new["y_by_z_max_abs_v"], marker="s", label="y seam")
    axes[1].set_title("Correct periodic Q1 field")
    axes[1].set_xlabel("z node")
    axes[1].set_ylabel("endpoint |delta V| (V)")
    axes[1].set_ylim(-0.05, 0.5)
    axes[1].legend()
    axes[1].grid(alpha=0.25)

    floor = 0
    axes[2].plot(archived[:, :, floor].ravel(), marker="o", label="archived")
    axes[2].plot(corrected[:, :, floor].ravel(), marker="s", label="periodicized")
    axes[2].set_title("Floor-node potential proposal")
    axes[2].set_xlabel("flattened full-grid node")
    axes[2].set_ylabel("V")
    axes[2].legend()
    axes[2].grid(alpha=0.25)
    figure.suptitle("C3 topology repair: endpoints become one field unknown")
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--grid-dx-um", type=float, default=0.25)
    args = parser.parse_args()
    checkpoint_path = args.checkpoint.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with np.load(checkpoint_path) as checkpoint:
        sigma = np.asarray(checkpoint["sigma_c_per_m2"], dtype=float)
        face_charge = np.asarray(checkpoint["face_charge_c"], dtype=float)
        archived_node_charge = np.asarray(checkpoint["charge_node_c"], dtype=float)
        archived_potential = np.asarray(checkpoint["potential_v"], dtype=float)
        vertices = np.asarray(checkpoint["vertices"])
        faces = np.asarray(checkpoint["faces"], dtype=int)
        centroids = np.asarray(checkpoint["centroids"], dtype=float)
        areas = np.asarray(checkpoint["areas"], dtype=float)
        face_material_id = np.asarray(checkpoint["face_material_id"], dtype=int)
        method_hint = np.asarray(checkpoint["method_hint_Ar"])
        resume_sampling_epoch = int(np.asarray(
            checkpoint["resume_sampling_epoch"]).item())
        scramble_mode = str(np.asarray(checkpoint["scramble_mode"]).item())
        scramble_base_seed = int(np.asarray(checkpoint["scramble_base_seed"]).item())
        sampling_seed_stride = int(np.asarray(
            checkpoint["sampling_seed_stride"]).item())

    geometry, periodic = _geometry_and_poisson(args.grid_dx_um)
    if archived_node_charge.shape != periodic.shape:
        raise ValueError("checkpoint and declared Poisson grid shapes disagree")
    legacy = NodalPoissonSystem3D(
        periodic.epsilon_r, periodic.spacing_m, periodic.dirichlet_mask,
        periodic.dirichlet_voltage, periodic_axes=())
    physical_area = areas * geometry.mesh_length_unit_m ** 2
    if not np.allclose(face_charge, sigma * physical_area, rtol=2e-13, atol=2e-31):
        raise ValueError("checkpoint face inventory is internally inconsistent")

    direct_node_charge = lump_triangle_sheet_charge_3d(
        periodic.shape, vertices, faces, sigma,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    legacy_potential, _ = legacy.solve(direct_node_charge)
    legacy_projector = CompatibleQ1SurfaceChargeProjector3D.from_triangles(
        legacy.shape, vertices, faces, grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    periodic_projector = CompatibleQ1SurfaceChargeProjector3D.from_poisson_system(
        periodic, vertices, faces, grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m)
    periodic_face_charge = periodic_projector.project_face_charge(face_charge)
    periodic_sigma = periodic_face_charge / physical_area
    periodic_node_charge = periodic.canonicalize_charge(lump_triangle_sheet_charge_3d(
        periodic.shape, vertices, faces, periodic_sigma,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=args.grid_dx_um,
        coordinate_length_unit_m=geometry.mesh_length_unit_m))
    periodic_potential, poisson_diagnostics = periodic.solve(periodic_node_charge)

    reduced_before = periodic_projector.node_charge_from_face_charge(face_charge)
    reduced_after = periodic_projector.node_charge_from_face_charge(periodic_face_charge)
    reduced_scale = max(float(np.sum(np.abs(reduced_before))), np.finfo(float).tiny)
    potential_scale = max(float(np.linalg.norm(archived_potential)), np.finfo(float).tiny)
    charge_scale = max(float(np.sum(np.abs(face_charge))), np.finfo(float).tiny)
    direct_scale = max(float(np.sum(np.abs(direct_node_charge))), np.finfo(float).tiny)

    warm_path = output / "periodicized_warm_proposal.npz"
    np.savez_compressed(
        warm_path,
        sigma_c_per_m2=periodic_sigma, face_charge_c=periodic_face_charge,
        charge_node_c=periodic_node_charge, potential_v=periodic_potential,
        vertices=vertices, faces=faces, centroids=centroids, areas=areas,
        face_material_id=face_material_id, method_hint_Ar=method_hint,
        resume_sampling_epoch=np.asarray(resume_sampling_epoch),
        scramble_mode=np.asarray(scramble_mode),
        scramble_base_seed=np.asarray(scramble_base_seed),
        sampling_seed_stride=np.asarray(sampling_seed_stride),
        compatible_q1_charge_state=np.asarray(True),
        poisson_periodic_axes=np.asarray(periodic.periodic_axes, dtype=int),
        source_checkpoint_sha256=np.asarray(_hash(checkpoint_path)),
        status=np.asarray(
            "warm proposal only; prior physical time invalid under corrected topology"))

    plot_path = output / "periodic_topology_repair.png"
    _plot(archived_potential, periodic_potential, plot_path)
    artifact = dict(
        schema="petch.charging.c3.periodic-topology-audit.v1",
        status=(
            "topology mismatch confirmed and repaired in the discrete operator; warm state is a "
            "proposal only, not a convergence claim"),
        mismatch=dict(
            particle_topology="periodic lateral x/y",
            archived_poisson_topology="independent lateral endpoint nodes (natural Neumann)",
            archived_checkpoint_seam=_seam_metrics(archived_potential),
            recomputed_legacy_seam=_seam_metrics(legacy_potential),
            archived_vs_recomputed_legacy_potential_relative_l2=float(
                np.linalg.norm(archived_potential - legacy_potential) / potential_scale),
            archived_vs_direct_node_charge_relative_l1=float(
                np.sum(np.abs(archived_node_charge - direct_node_charge)) / direct_scale)),
        repair=dict(
            poisson_periodic_axes=list(periodic.periodic_axes),
            full_node_shape=list(periodic.shape),
            independent_node_shape=list(periodic.reduced_shape),
            corrected_seam=_seam_metrics(periodic_potential),
            legacy_face_coupling=dict(
                shape=list(legacy_projector.coupling.shape), rank=legacy_projector.rank,
                nullity=legacy_projector.nullity,
                condition_number=legacy_projector.condition_number),
            periodic_face_coupling=dict(
                shape=list(periodic_projector.coupling.shape), rank=periodic_projector.rank,
                nullity=periodic_projector.nullity,
                condition_number=periodic_projector.condition_number),
            effective_nodal_load_relative_l1_error=float(
                np.sum(np.abs(reduced_after - reduced_before)) / reduced_scale),
            global_face_charge_error_c=float(
                periodic_face_charge.sum() - face_charge.sum()),
            global_face_charge_relative_error=float(
                abs(periodic_face_charge.sum() - face_charge.sum()) / charge_scale),
            input_periodic_null_fraction=periodic_projector.unresolved_fraction(face_charge),
            output_periodic_null_fraction=periodic_projector.unresolved_fraction(
                periodic_face_charge),
            maximum_potential_change_v=float(np.max(np.abs(
                periodic_potential - archived_potential))),
            potential_relative_l2_change=float(
                np.linalg.norm(periodic_potential - archived_potential) / potential_scale),
            poisson_free_charge_residual_max_c=(
                poisson_diagnostics.max_abs_free_charge_residual_c),
            poisson_global_charge_balance_c=poisson_diagnostics.charge_balance_c),
        next_test=dict(
            action=(
                "run short timestep-refined fresh-scramble branches from zero charge and this "
                "periodicized warm proposal; require convergence toward the same corrected branch"),
            no_replay_statement=(
                "the archived 2.94875 ms is not replayed or credited; only its instantaneous "
                "inventory is used as a proposal")),
        provenance=dict(
            checkpoint_name=checkpoint_path.name,
            checkpoint_sha256=_hash(checkpoint_path),
            script_sha256=_hash(Path(__file__).resolve()),
            charging_poisson_sha256=_hash(ROOT / "src/petch/charging_poisson_3d.py")),
        artifacts=dict(
            warm_proposal=dict(name=warm_path.name, sha256=_hash(warm_path)),
            plot=dict(name=plot_path.name, sha256=_hash(plot_path))))
    audit_path = output / "audit.json"
    audit_path.write_text(json.dumps(_json_value(artifact), indent=2) + "\n")
    print(json.dumps(_json_value({
        "status": artifact["status"], "mismatch": artifact["mismatch"],
        "repair": artifact["repair"], "artifacts": artifact["artifacts"]}), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
