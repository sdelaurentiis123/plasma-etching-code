#!/usr/bin/env python3
"""Bounded fixed-geometry C3 trench trajectory with face-authoritative checkpoints."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
from time import perf_counter

for variable in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import (  # noqa: E402
    _boundary, _electron_proposal, _geometry_and_poisson, _ion_proposal,
)
from petch.charged_surface_response_3d import (  # noqa: E402
    GrazingSpecularIonReflection3D,
)
from petch.charging_coevolution_3d import (  # noqa: E402
    SurfaceChargingSaturationError,
    integrate_surface_charging_to_saturation_3d,
)
from petch.charging_poisson_3d import lump_triangle_sheet_charge_3d  # noqa: E402
from petch.feature_step_3d import (  # noqa: E402
    _face_material_ids, _surface_gas_normals,
)
from petch.threed import extract_mesh_3d  # noqa: E402


def json_value(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict) or hasattr(value, "items"):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__}")


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "results/charging_coevolution_c3_trench_pilot")
    parser.add_argument("--grid-dx", type=float, default=0.25)
    parser.add_argument("--timestep-s", type=float, default=1.25e-7)
    parser.add_argument("--maximum-steps", type=int, default=2)
    parser.add_argument("--timestep-policy", choices=("fixed", "ser"), default="fixed")
    parser.add_argument("--maximum-timestep-s", type=float, default=5e-7)
    parser.add_argument("--minimum-timestep-s", type=float, default=1e-11)
    parser.add_argument("--ser-activation-rms", type=float, default=0.5)
    parser.add_argument("--ser-maximum-growth", type=float, default=2.0)
    parser.add_argument("--ser-allowed-residual-growth", type=float, default=0.005)
    parser.add_argument("--potential-rate-tolerance-v-s", type=float, default=1e3)
    parser.add_argument("--patch-scales-um", type=float, nargs="+", default=(0.25, 0.5))
    parser.add_argument("--initial-face-state", type=Path)
    parser.add_argument(
        "--method-map", type=Path,
        default=ROOT / "results/charging_task1_3d_bidir_cont6/final_states.npz")
    parser.add_argument("--method-key", default="refined_method_hint_Ar+")
    parser.add_argument("--forward-level", type=int, default=10)
    parser.add_argument("--adjoint-level", type=int, default=8)
    parser.add_argument("--n-position", type=int, default=256)
    parser.add_argument("--seed", type=int, default=79)
    parser.add_argument("--scramble-mode", choices=("frozen", "fresh"), default="frozen")
    parser.add_argument("--sampling-seed-stride", type=int, default=1000003)
    parser.add_argument(
        "--initial-sampling-epoch", type=int,
        help=("fresh-scramble epoch attached to an input checkpoint; when checkpoint metadata "
              "exists it is loaded automatically and an explicit value must agree"))
    parser.add_argument("--trajectory-dt", type=float, default=0.0003125)
    parser.add_argument("--trajectory-max-steps", type=int, default=128000)
    parser.add_argument(
        "--transport-device", choices=("cpu", "cuda", "cuda:0"), default="cpu")
    parser.add_argument("--response-tail-tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    if args.maximum_steps < 0:
        parser.error("--maximum-steps must be nonnegative")
    if len(args.patch_scales_um) < 2:
        parser.error("at least two --patch-scales-um values are required")

    geometry, poisson = _geometry_and_poisson(args.grid_dx)
    vertices, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(vertices, faces, centroids, geometry)
    material = _face_material_ids(centroids, geometry)
    source_z = 2.0
    boundary = _boundary(source_z * geometry.mesh_length_unit_m)
    with np.load(args.method_map) as archived:
        method_hint = np.asarray(archived[args.method_key]).astype("U7")
        archived_faces = np.asarray(archived["faces"], dtype=int)
        archived_vertices = np.asarray(archived["vertices"], dtype=float)
    if (method_hint.shape != (len(faces),) or not np.array_equal(faces, archived_faces)
            or not np.array_equal(vertices, archived_vertices)):
        parser.error("the separately selected estimator map does not match this exact mesh")

    initial_sigma = np.zeros(len(faces))
    initial_state_sha256 = None
    checkpoint_sampling = {}
    if args.initial_face_state is not None:
        initial_state_sha256 = file_hash(args.initial_face_state)
        with np.load(args.initial_face_state) as archived:
            initial_sigma = np.asarray(archived["sigma_c_per_m2"], dtype=float).copy()
            state_faces = np.asarray(archived["faces"], dtype=int)
            state_vertices = np.asarray(archived["vertices"], dtype=float)
            for name in (
                    "resume_sampling_epoch", "scramble_mode", "scramble_base_seed",
                    "sampling_seed_stride"):
                if name in archived:
                    checkpoint_sampling[name] = np.asarray(archived[name]).item()
        if (initial_sigma.shape != (len(faces),) or not np.array_equal(faces, state_faces)
                or not np.array_equal(vertices, state_vertices)):
            parser.error("initial face-charge checkpoint does not match this exact mesh")

    if args.scramble_mode == "fresh":
        checkpoint_epoch = checkpoint_sampling.get("resume_sampling_epoch")
        if checkpoint_epoch is not None:
            checkpoint_epoch = int(checkpoint_epoch)
            if (args.initial_sampling_epoch is not None
                    and args.initial_sampling_epoch != checkpoint_epoch):
                parser.error("--initial-sampling-epoch conflicts with checkpoint metadata")
            args.initial_sampling_epoch = checkpoint_epoch
            if str(checkpoint_sampling.get("scramble_mode", "fresh")) != "fresh":
                parser.error("fresh continuation requires a fresh-scramble checkpoint")
            if ("scramble_base_seed" in checkpoint_sampling
                    and int(checkpoint_sampling["scramble_base_seed"]) != args.seed):
                parser.error("--seed conflicts with fresh checkpoint metadata")
            if ("sampling_seed_stride" in checkpoint_sampling
                    and int(checkpoint_sampling["sampling_seed_stride"])
                    != args.sampling_seed_stride):
                parser.error("--sampling-seed-stride conflicts with checkpoint metadata")
        elif args.initial_face_state is not None and args.initial_sampling_epoch is None:
            parser.error(
                "fresh continuation from a legacy checkpoint requires --initial-sampling-epoch")
        elif args.initial_sampling_epoch is None:
            args.initial_sampling_epoch = 0
    else:
        if args.initial_sampling_epoch not in (None, 0):
            parser.error("--initial-sampling-epoch is only valid in fresh mode")
        args.initial_sampling_epoch = 0

    reflection = GrazingSpecularIonReflection3D.literature_bounded_sensitivity(1, "Ar+")
    config = dict(
        status="bounded coarse real-trench pilot; not a convergence claim",
        geometry=dict(
            cell_width_um=1.0, cell_length_um=0.5, domain_height_um=2.0,
            dx_um=args.grid_dx, opening_width_um=0.5, mask_thickness_um=0.25,
            substrate_top_um=1.25, etched_depth_um=0.75),
        timestep_s=args.timestep_s, maximum_steps=args.maximum_steps,
        timestep_policy=args.timestep_policy,
        maximum_timestep_s=args.maximum_timestep_s,
        minimum_timestep_s=args.minimum_timestep_s,
        ser_activation_rms=args.ser_activation_rms,
        ser_maximum_growth=args.ser_maximum_growth,
        ser_allowed_residual_growth=args.ser_allowed_residual_growth,
        potential_rate_tolerance_v_s=args.potential_rate_tolerance_v_s,
        patch_scales_um=list(args.patch_scales_um), n_position=args.n_position,
        seed=args.seed,
        scramble_mode=args.scramble_mode,
        sampling_seed_stride=args.sampling_seed_stride,
        initial_sampling_epoch=args.initial_sampling_epoch,
        adjoint_proposal_seeds={"Ar+": args.seed, "electron": args.seed + 4},
        trajectory_dt=args.trajectory_dt,
        trajectory_max_steps=args.trajectory_max_steps,
        transport_device=args.transport_device,
        forward_level=args.forward_level, adjoint_level=args.adjoint_level,
        initial_face_state_sha256=initial_state_sha256,
        response_tail_tolerance=args.response_tail_tolerance,
        method_map_sha256=file_hash(args.method_map), method_key=args.method_key,
        estimator_map_source="separate pre-C3 pilot; estimator choice only, no nodal charge",
        exact_operator="hard visibility with bounded grazing-ion reflection")
    encoded = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config_hash = sha256(encoded).hexdigest()
    source_paths = (
        ROOT / "src/petch/charging_coevolution_3d.py",
        ROOT / "src/petch/charging_coupled_3d.py",
        ROOT / "src/petch/boundary_transport_3d.py",
        ROOT / "src/petch/charged_surface_cascade_3d.py",
        ROOT / "src/petch/charged_surface_response_3d.py",
        Path(__file__).resolve(),
    )

    def run_manifest(wall_clock_s):
        return dict(
            engine_git_revision=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            source_sha256={
                path.relative_to(ROOT).as_posix(): file_hash(path) for path in source_paths},
            hardware=platform.platform(), python=platform.python_version(),
            wall_clock_s=wall_clock_s,
            reflection_provenance=json_value(reflection.provenance))

    started = perf_counter()
    def proposal_factory(evaluation_seed):
        return {
            "Ar+": _ion_proposal(boundary, args.adjoint_level, evaluation_seed),
            "electron": _electron_proposal(
                boundary, args.adjoint_level, evaluation_seed + 4),
        }

    try:
        result = integrate_surface_charging_to_saturation_3d(
            poisson, initial_sigma, boundary, vertices, faces, areas,
            face_centroids=centroids, face_gas_normals=normals,
            face_material_id=material,
            source_bounds=(0.0, 1.0, 0.0, 0.5), source_z=source_z,
            potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
            patch_scales_m=tuple(value * 1e-6 for value in args.patch_scales_um),
            potential_rate_tolerance_v_s=args.potential_rate_tolerance_v_s,
            timestep_s=args.timestep_s, maximum_steps=args.maximum_steps,
            current_balance_tolerance=0.08, timestep_policy=args.timestep_policy,
            maximum_timestep_s=args.maximum_timestep_s,
            minimum_timestep_s=args.minimum_timestep_s,
            ser_activation_rms=args.ser_activation_rms,
            ser_maximum_growth=args.ser_maximum_growth,
            ser_allowed_residual_growth=args.ser_allowed_residual_growth,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            mesh_origin_m=geometry.mesh_origin_m,
            n_position=args.n_position, seed=args.seed,
            scramble_mode=args.scramble_mode,
            sampling_seed_stride=args.sampling_seed_stride,
            initial_sampling_epoch=args.initial_sampling_epoch,
            trajectory_fixed_dt=args.trajectory_dt,
            trajectory_max_steps=args.trajectory_max_steps,
            phase_space_log2_samples=args.forward_level,
            periodic_lateral=True,
            transport_estimator={"Ar+": "bidirectional", "electron": "adjoint"},
            adjoint_face_quadrature_points=3, adjoint_ray_offset=1e-4,
            adjoint_proposals=proposal_factory(args.seed),
            fresh_adjoint_proposal_factory=(
                proposal_factory if args.scramble_mode == "fresh" else None),
            adjoint_proposal_frames={"Ar+": "source_aligned", "electron": "surface_local"},
            bidirectional_options=dict(
                forward_log2_samples=args.forward_level,
                adjoint_log2_samples=args.adjoint_level, n_replicates=4,
                method_hint={"Ar+": method_hint}, require_certification=False,
                element_absolute_tolerance=0.02, element_relative_tolerance=0.1,
                face_quadrature_points=3),
            transport_device=args.transport_device, charged_surface_response=reflection,
            response_launch_offset=1e-5, response_max_bounces=16,
            response_relative_tail_tolerance=args.response_tail_tolerance)
    except SurfaceChargingSaturationError as error:
        wall_clock_s = perf_counter() - started
        history = tuple(dict(item) for item in error.history)
        last = history[-1] if history else {}
        failure_sigma = np.asarray(error.sigma_c_per_m2)
        failure_charge = lump_triangle_sheet_charge_3d(
            poisson.shape, vertices, faces, failure_sigma,
            grid_origin=(0.0, 0.0, 0.0), grid_spacing=geometry.dx,
            coordinate_length_unit_m=geometry.mesh_length_unit_m)
        failure_potential, _ = poisson.solve(failure_charge)
        failure_summary = dict(
            schema="petch.charging.coevolution.c3.trench-pilot.v1",
            config_hash=config_hash, config=config, run_manifest=run_manifest(wall_clock_s),
            result=dict(
                converged=False, failed=True, error_type=type(error).__name__,
                error_message=str(error), accepted_steps=error.accepted_steps,
                rejected_steps=error.rejected_steps,
                physical_time_s=error.physical_time_s,
                pseudo_time_s=error.pseudo_time_s,
                resume_sampling_epoch=(
                    args.initial_sampling_epoch + error.accepted_steps
                    if args.scramble_mode == "fresh" else 0),
                retained_node_rms_relative_current_imbalance=last.get(
                    "rms_relative_current_imbalance_node"),
                retained_node_max_relative_current_imbalance=last.get(
                    "max_relative_current_imbalance_node"),
                final_potential_rate_max_v_s=last.get("potential_rate_max_v_s"),
                patch_scales_m=last.get("patch_scales_m", ()),
                patch_b2_max_ion_normalized=last.get("patch_max_relative_imbalance", ()),
                maximum_charge_conservation_relative_error=max(
                    (item["charge_conservation_relative_error"] for item in history),
                    default=None),
                maximum_surface_transfer_relative_charge_balance_error=max(
                    (item["surface_transfer_relative_charge_balance_error"] for item in history),
                    default=None),
                maximum_response_tail_closure_relative_absolute_charge_rate=max(
                    (item["response_tail_closure_relative_absolute_charge_rate"]
                     for item in history), default=None),
                maximum_response_tail_closure_l1_current_error_bound_relative=max(
                    (item["response_tail_closure_l1_current_error_bound_relative"]
                     for item in history), default=None),
                maximum_transport_lineage_replay_count=max(
                    (item.get("transport_lineage_replay_count", 0) for item in history),
                    default=None),
                maximum_transport_lineage_replay_fraction=max(
                    (item.get("transport_lineage_replay_fraction", 0.0) for item in history),
                    default=None),
                minimum_potential_v=float(np.min(failure_potential)),
                maximum_potential_v=float(np.max(failure_potential))),
            history=json_value(history),
            conclusion="solver refused; replayable face state and failure diagnostics retained")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
        (args.output_dir / "summary.json").write_text(
            json.dumps(failure_summary, indent=2) + "\n")
        np.savez_compressed(
            args.output_dir / "face_checkpoint.npz",
            sigma_c_per_m2=failure_sigma,
            face_charge_c=failure_sigma * areas * geometry.mesh_length_unit_m ** 2,
            charge_node_c=failure_charge, potential_v=failure_potential,
            vertices=vertices, faces=faces, centroids=centroids, areas=areas,
            face_material_id=material, method_hint_Ar=method_hint,
            resume_sampling_epoch=np.asarray(
                args.initial_sampling_epoch + error.accepted_steps
                if args.scramble_mode == "fresh" else 0),
            scramble_mode=np.asarray(args.scramble_mode),
            scramble_base_seed=np.asarray(args.seed),
            sampling_seed_stride=np.asarray(args.sampling_seed_stride))
        print(json.dumps(failure_summary, indent=2), flush=True)
        raise SystemExit(2) from error
    wall_clock_s = perf_counter() - started
    summary = dict(
        schema="petch.charging.coevolution.c3.trench-pilot.v1",
        config_hash=config_hash, config=config,
        run_manifest=run_manifest(wall_clock_s),
        result=dict(
            converged=result.converged, accepted_steps=result.accepted_steps,
            rejected_steps=result.rejected_steps, physical_time_s=result.physical_time_s,
            pseudo_time_s=result.pseudo_time_s,
            resume_sampling_epoch=result.diagnostics["resume_sampling_epoch"],
            retained_node_rms_relative_current_imbalance=(
                result.diagnostics["retained_node_rms_relative_current_imbalance"]),
            retained_node_max_relative_current_imbalance=(
                result.diagnostics["retained_node_max_relative_current_imbalance"]),
            final_potential_rate_max_v_s=result.diagnostics[
                "final_potential_rate_max_v_s"],
            patch_scales_m=list(result.diagnostics["patch_scales_m"]),
            patch_b2_max_ion_normalized=[
                item.b2_maximum_ion_normalized_imbalance for item in result.patch_balance],
            patch_symmetric_max=[
                item.maximum_relative_imbalance for item in result.patch_balance],
            maximum_charge_conservation_relative_error=max(
                item["charge_conservation_relative_error"] for item in result.history),
            maximum_surface_transfer_relative_charge_balance_error=max(
                item["surface_transfer_relative_charge_balance_error"]
                for item in result.history),
            maximum_response_tail_closure_relative_absolute_charge_rate=max(
                item["response_tail_closure_relative_absolute_charge_rate"]
                for item in result.history),
            maximum_response_tail_closure_l1_current_error_bound_relative=max(
                item["response_tail_closure_l1_current_error_bound_relative"]
                for item in result.history),
            maximum_transport_lineage_replay_count=max(
                item.get("transport_lineage_replay_count", 0) for item in result.history),
            maximum_transport_lineage_replay_fraction=max(
                item.get("transport_lineage_replay_fraction", 0.0)
                for item in result.history),
            minimum_potential_v=float(np.min(result.potential_v)),
            maximum_potential_v=float(np.max(result.potential_v))),
        history=json_value(result.history),
        conclusion=("C3 saturation gates pass" if result.converged else
                    "bounded pilot only; continue from face checkpoint with refinement evidence"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    checkpoint_path = args.output_dir / "face_checkpoint.npz"
    current_audit_path = args.output_dir / "current_audit.npz"
    np.savez_compressed(
        checkpoint_path,
        sigma_c_per_m2=result.sigma_c_per_m2,
        face_charge_c=result.face_charge_c,
        charge_node_c=result.charge_node_c,
        potential_v=result.potential_v,
        vertices=vertices, faces=faces, centroids=centroids, areas=areas,
        face_material_id=material, method_hint_Ar=method_hint,
        resume_sampling_epoch=np.asarray(result.diagnostics["resume_sampling_epoch"]),
        scramble_mode=np.asarray(args.scramble_mode),
        scramble_base_seed=np.asarray(args.seed),
        sampling_seed_stride=np.asarray(args.sampling_seed_stride))
    np.savez_compressed(
        current_audit_path,
        positive_face_current_density_a_m2=(
            result.final_step.positive_face_current_density_a_m2),
        negative_face_current_density_a_m2=(
            result.final_step.negative_face_current_density_a_m2),
        net_face_current_density_a_m2=result.final_step.face_current_density_a_m2,
        positive_current_node_a=result.final_step.positive_current_node_a,
        negative_current_node_a=result.final_step.negative_current_node_a,
        potential_before_v=result.final_step.potential_before_v,
        potential_after_v=result.final_step.potential_after_v,
        potential_rate_v_s=(
            result.final_step.potential_after_v - result.final_step.potential_before_v
        ) / args.timestep_s,
        patch_scales_m=np.asarray(
            [item.patch_scale_m for item in result.patch_balance], dtype=float),
        patch_group_by_scale=np.stack(
            [item.group for item in result.patch_balance], axis=0),
        physical_face_area_m2=areas * geometry.mesh_length_unit_m ** 2)
    summary["artifacts"] = {
        "face_checkpoint": {
            "name": checkpoint_path.name, "sha256": file_hash(checkpoint_path)},
        "current_audit": {
            "name": current_audit_path.name, "sha256": file_hash(current_audit_path)},
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
