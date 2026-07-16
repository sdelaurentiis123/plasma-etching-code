#!/usr/bin/env python3
"""Run the source-scale Hwang--Giapis Figure 13 notch-profile replay.

This is a quantitative development replay, not a held-out validation claim:
the experimental contour was inspected during implementation and its
measurement uncertainty was not reported.  No profile datum is fitted.

The expensive 20 nm global charging/transport boundary is cached separately
from the 5 nm local profile refinement.  Reruns therefore reuse the exact
charged ion lineage unless the declared global configuration changes.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from petch.charging2d import solve_edge_array_charging
from petch.experimental_boundary import (
    HWANG_GIAPIS_1997_EEDF_SHA256,
    HWANG_GIAPIS_1997_IEDF_SHA256,
    build_hwang_giapis_1997_boundary_state,
)
from petch.hwang_giapis_notch_profile_2d import (
    HwangGiapisLocalNotchCheckpoint2D,
    evolve_hwang_giapis_local_notch_event_driven_2d,
    hwang_giapis_local_boundary_from_edge_array_result,
)
from petch.hwang_giapis_notch_validation_2d import (
    HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256,
    load_hwang_giapis_1997_fig13_profile,
    score_hwang_giapis_1997_fig13_profile,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "hwang_giapis_1997_fig13_validation"
EVIDENCE = (
    ROOT / "data" / "experimental" / "hwang_giapis_1997"
    / "fig13_notch_profile.csv")
IEDF_EVIDENCE = (
    ROOT / "data" / "experimental" / "hwang_giapis_1997"
    / "fig4a_ion_energy_distribution.csv")
EEDF_EVIDENCE = (
    ROOT / "data" / "experimental" / "hwang_giapis_1997"
    / "fig4b_electron_energy_distribution.csv")
GLOBAL_SCHEMA = "petch-hwang-giapis-fig13-global-boundary-v6"
LEGACY_GLOBAL_SCHEMAS = {
    "petch-hwang-giapis-fig13-global-boundary-v5",
    GLOBAL_SCHEMA,
}
AUDIT_SCHEMA = "petch-hwang-giapis-fig13-profile-audit-v1"


def _git_revision():
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=False,
        capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def _stable_hash(value):
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256(encoded).hexdigest()


def _file_sha256(path):
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _relative_or_name(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return path.name


def _continuation_artifacts(source):
    source = Path(source).resolve()
    if source.is_file():
        source = source.parent
    metadata_path = source / "global_boundary.json"
    arrays_path = source / "global_boundary_arrays.npz"
    if not metadata_path.is_file() or not arrays_path.is_file():
        raise ValueError(
            "global warm start must contain global_boundary.json and "
            "global_boundary_arrays.npz")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") not in LEGACY_GLOBAL_SCHEMAS:
        raise ValueError("global warm start uses an unsupported schema")
    return source, metadata_path, arrays_path, metadata


def _continuation_provenance(source):
    source, metadata_path, arrays_path, metadata = (
        _continuation_artifacts(source))
    return {
        "directory": _relative_or_name(source),
        "metadata_path": _relative_or_name(metadata_path),
        "metadata_sha256": _file_sha256(metadata_path),
        "arrays_path": _relative_or_name(arrays_path),
        "arrays_sha256": _file_sha256(arrays_path),
        "source_schema": metadata["schema"],
        "source_config_sha256": metadata.get("config_sha256"),
    }


def _load_continuation_state(config):
    provenance = config.get("continuation_source")
    if provenance is None:
        return {}
    source = ROOT / provenance["directory"]
    _, metadata_path, arrays_path, metadata = _continuation_artifacts(source)
    if (_file_sha256(metadata_path) != provenance["metadata_sha256"]
            or _file_sha256(arrays_path) != provenance["arrays_sha256"]):
        raise ValueError("global warm-start artifact checksum changed")
    with np.load(arrays_path) as arrays:
        if "continuation_surface_potential_v" in arrays:
            surface = np.asarray(
                arrays["continuation_surface_potential_v"], dtype=float)
        else:
            # Version 5 saved the tail-averaged Laplace field.  On
            # nonconducting surface cells this is exactly the declared
            # potential state, so it is a valid warm proposal.  Version 6
            # additionally saves the raw stochastic endpoint.
            surface = np.asarray(arrays["potential_v"], dtype=float)
    continuation = metadata.get("continuation_state", {})
    observables = metadata.get("global_observables", {})
    gain = observables.get("stochastic_gain", {})
    source_config = metadata.get("config", {})
    age = continuation.get(
        "stochastic_gain_age",
        gain.get(
            "ending_age",
            int(gain.get("starting_age", 0))
            + int(source_config.get("n_iter", 0))))
    return {
        "initial_surface_potential_v": surface,
        "initial_edge_potential_v": float(continuation.get(
            "edge_potential_v",
            observables["edge_poly_potential_v"])),
        "initial_neighbor_potential_v": float(continuation.get(
            "neighbor_potential_v",
            observables["neighbor_poly_potential_v"])),
        "stochastic_gain_age": int(age),
    }


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _global_config(args):
    config = {
        "AR": 2.6,
        "W": 25,
        "mouth": 185,
        "Te": 4.0,
        "V_dc": 37.0,
        "V_rf": 30.0,
        "iadf_hwhm_deg": 4.3,
        "n_per_iter": int(args.global_samples),
        "n_iter": int(args.global_iterations),
        "seed": int(args.seed),
        "poly_um": 0.3,
        "feature_w_um": 0.5,
        "rf_bursts": True,
        "sheath_um": 89.0,
        "boundary_um": 3.7,
        "ion_angle_energy_corr": "anticorrelated",
        "source_model": "primary_digitized_PlasmaBoundaryState",
        "open_width_um": 2.0,
        "right_buffer_um": 0.5,
        "domain_model": "hwang_mirror_cell",
        "edge_open_model": "none",
        "electron_model": "mc",
        "relax": float(args.global_relax),
        "stochastic_gain_exponent": float(args.global_gain_exponent),
        "stochastic_gain_offset": float(args.global_gain_offset),
        "return_final_ion_lineage": True,
        "final_audit_seed": int(
            args.global_final_seed
            if args.global_final_seed is not None
            else args.seed + 1_000_003),
        "plasma_boundary_iedf_path": str(IEDF_EVIDENCE.relative_to(ROOT)),
        "plasma_boundary_iedf_sha256": HWANG_GIAPIS_1997_IEDF_SHA256,
        "plasma_boundary_eedf_path": str(EEDF_EVIDENCE.relative_to(ROOT)),
        "plasma_boundary_eedf_sha256": HWANG_GIAPIS_1997_EEDF_SHA256,
        "plasma_boundary_eadf_cosine_power": 0.6,
        "plasma_boundary_dimensional_projection": (
            "source_faithful_2d_energy_angle_laws"),
        "plasma_boundary_reference_plane_m": 3.7e-6,
    }
    if args.global_final_samples is not None:
        config["final_audit_samples"] = int(args.global_final_samples)
    if args.ion_tangential_temperature_eV is not None:
        config["plasma_boundary_ion_tangential_temperature_eV"] = float(
            args.ion_tangential_temperature_eV)
    if args.global_warm_start is not None:
        config["continuation_source"] = _continuation_provenance(
            args.global_warm_start)
    return config


def _global_solver_kwargs(config, *, progress_callback=None):
    kwargs = {
        key: value for key, value in config.items()
        if (not key.startswith("plasma_boundary_")
            and not key.startswith("continuation_"))
    }
    kwargs.update(_load_continuation_state(config))
    kwargs["progress_callback"] = progress_callback
    kwargs["plasma_boundary"] = build_hwang_giapis_1997_boundary_state(
        IEDF_EVIDENCE, EEDF_EVIDENCE,
        reference_plane_m=float(config["plasma_boundary_reference_plane_m"]),
        ion_tangential_temperature_eV=float(config.get(
            "plasma_boundary_ion_tangential_temperature_eV", 0.5)))
    return kwargs


def _save_global_boundary(output, config, result, elapsed_s):
    lineage = result["final_ion_lineage"]
    floor_line = np.asarray(result["Vfloor"], dtype=float)
    edge_exclusion = max(1, floor_line.size // 5)
    floor_interior = floor_line[
        edge_exclusion:floor_line.size - edge_exclusion]
    np.savez_compressed(
        output / "global_boundary_arrays.npz",
        potential_v=result["V"],
        continuation_surface_potential_v=(
            result["continuation_state"]["surface_potential_v"]),
        hit_type=lineage["hit_type"],
        hit_ix=lineage["hit_ix"],
        hit_iz=lineage["hit_iz"],
        impact_energy_eV=lineage["impact_energy_eV"],
        hit_x_grid=lineage["hit_x_grid"],
        hit_z_grid=lineage["hit_z_grid"],
        hit_vx_sqrt_eV=lineage["hit_vx_sqrt_eV"],
        hit_vz_sqrt_eV=lineage["hit_vz_sqrt_eV"])
    metadata = {
        "schema": GLOBAL_SCHEMA,
        "config": config,
        "config_sha256": _stable_hash(config),
        "elapsed_s": float(elapsed_s),
        "geometry": result["geom"],
        "poly_potential_v": float(result["V_poly_edge"]),
        "continuation_state": {
            key: value
            for key, value in result["continuation_state"].items()
            if not isinstance(value, np.ndarray)
        },
        "lineage": {
            key: value for key, value in lineage.items()
            if not isinstance(value, np.ndarray)
        },
        "global_observables": {
            "floor_flux": result["floor_flux"],
            "floor_flux_tail": result["floor_flux_tail"],
            "floor_potential_center_v": result["V_floor_center"],
            "floor_potential_peak_v": result["V_foot_peak"],
            "floor_potential_mean_v": float(np.mean(floor_line)),
            "floor_potential_median_v": float(np.median(floor_line)),
            "floor_potential_interior_mean_v": float(
                np.mean(floor_interior)),
            "floor_potential_interior_median_v": float(
                np.median(floor_interior)),
            "edge_poly_potential_v": result["V_poly_edge"],
            "neighbor_poly_potential_v": result["V_poly_neighbor"],
            "foot_ion_flux": result["foot_ion_flux"],
            "foot_ion_mean_energy_eV": result["foot_ion_Emean"],
            "tail_residual": result["diag"]["residual"],
            "snapshot_residual": result["diag"]["residual_snapshot"],
            "region_current_balance": (
                result["diag"]["region_current_balance"]),
            "potential_stationarity": (
                result["diag"]["potential_stationarity"]),
            "stochastic_gain": result["diag"]["stochastic_gain"],
            "potential_history": result["diag"]["potential_history"],
            "trace": result["diag"]["trace"],
        },
    }
    (output / "global_boundary.json").write_text(
        json.dumps(_json_value(metadata), indent=2) + "\n", encoding="utf-8")


def _load_global_boundary(output, config):
    metadata_path = output / "global_boundary.json"
    arrays_path = output / "global_boundary_arrays.npz"
    if not metadata_path.exists() or not arrays_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (metadata.get("schema") != GLOBAL_SCHEMA
            or metadata.get("config_sha256") != _stable_hash(config)):
        return None
    arrays = np.load(arrays_path)
    lineage = {
        name: arrays[name] for name in (
            "hit_type", "hit_ix", "hit_iz", "impact_energy_eV",
            "hit_x_grid", "hit_z_grid", "hit_vx_sqrt_eV",
            "hit_vz_sqrt_eV")
    }
    lineage.update(metadata["lineage"])
    return {
        "V": arrays["potential_v"],
        "V_poly_edge": metadata["poly_potential_v"],
        "geom": metadata["geometry"],
        "final_ion_lineage": lineage,
        "_cached_metadata": metadata,
    }


def _score_dict(score):
    return {
        "rmse_um": score.rmse_um,
        "mean_absolute_error_um": score.mean_absolute_error_um,
        "maximum_absolute_error_um": score.maximum_absolute_error_um,
        "maximum_depth_error_um": score.maximum_depth_error_um,
        "digitization_bound_coverage_fraction": (
            score.digitization_bound_coverage_fraction),
        "strict_validation_pass": score.strict_validation_pass,
        "claim_status": score.claim_status,
    }


def _completed_local_run(output, label, config_sha256):
    metadata_path = output / f"{label}.json"
    arrays_path = output / f"{label}.npz"
    if not metadata_path.exists() or not arrays_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("config_sha256") != config_sha256:
        return None
    metadata.setdefault("landed_photoresist_count", 0.0)
    arrays = np.load(arrays_path)
    metadata["_depth_um"] = arrays["notch_depth_by_height_um"]
    return metadata


def _save_local_run(output, label, config, result, score, elapsed_s):
    count_closure = (
        result.landed_poly_count + result.landed_oxide_count
        + result.landed_photoresist_count + result.escaped_count
        - result.launched_count)
    remaining_reactive = float(np.sum(result.reactive_collision_inventory))
    reactive_closure = (
        result.threshold_removed_reactive_collisions
        + result.detached_reactive_collisions + remaining_reactive
        - result.direct_reactive_collisions
        - result.scattered_reactive_collisions)
    cell_area_um2 = float(config["cell_size_um"]) ** 2
    metadata = {
        "label": label,
        "config": config,
        "config_sha256": _stable_hash(config),
        "elapsed_s": float(elapsed_s),
        "maximum_notch_depth_um": result.maximum_notch_depth_um,
        "removed_cell_count": result.removed_cell_count,
        "threshold_removed_cell_count": result.threshold_removed_cell_count,
        "detached_cell_count": result.detached_cell_count,
        "threshold_removed_area_um2": (
            result.threshold_removed_cell_count * cell_area_um2),
        "detached_area_um2": result.detached_cell_count * cell_area_um2,
        "threshold_removed_reactive_collisions": (
            result.threshold_removed_reactive_collisions),
        "detached_reactive_collisions": (
            result.detached_reactive_collisions),
        "remaining_reactive_collision_inventory": remaining_reactive,
        "reactive_collision_ledger_error": reactive_closure,
        "direct_reactive_collisions": result.direct_reactive_collisions,
        "scattered_reactive_collisions": result.scattered_reactive_collisions,
        "launched_count": result.launched_count,
        "landed_poly_count": result.landed_poly_count,
        "landed_oxide_count": result.landed_oxide_count,
        "landed_photoresist_count": result.landed_photoresist_count,
        "escaped_count": result.escaped_count,
        "particle_count_ledger_error": count_closure,
        "score": _score_dict(score),
        "provenance": dict(result.provenance),
    }
    np.savez_compressed(
        output / f"{label}.npz",
        notch_depth_by_height_um=result.notch_depth_by_height_um,
        poly_cell=result.poly_cell,
        reactive_collision_inventory=result.reactive_collision_inventory,
        cumulative_oxide_ion_count=result.cumulative_oxide_ion_count)
    (output / f"{label}.json").write_text(
        json.dumps(_json_value(metadata), indent=2) + "\n",
        encoding="utf-8")
    return metadata


def _load_local_checkpoint(output, label, config_sha256):
    metadata_path = output / f"{label}_checkpoint.json"
    arrays_path = output / f"{label}_checkpoint.npz"
    if not metadata_path.exists() or not arrays_path.exists():
        return None
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("config_sha256") != config_sha256:
        return None
    arrays = np.load(arrays_path)
    return HwangGiapisLocalNotchCheckpoint2D(
        poly_cell=arrays["poly_cell"],
        reactive_collision_inventory=arrays[
            "reactive_collision_inventory"],
        cumulative_oxide_ion_count=arrays["cumulative_oxide_ion_count"],
        remaining_campaign_fraction=metadata[
            "remaining_campaign_fraction"],
        direct_reactive_collisions=metadata[
            "direct_reactive_collisions"],
        scattered_reactive_collisions=metadata[
            "scattered_reactive_collisions"],
        landed_poly_count=metadata["landed_poly_count"],
        landed_oxide_count=metadata["landed_oxide_count"],
        escaped_count=metadata["escaped_count"],
        threshold_removed_cell_count=metadata[
            "threshold_removed_cell_count"],
        detached_cell_count=metadata["detached_cell_count"],
        threshold_removed_reactive_collisions=metadata[
            "threshold_removed_reactive_collisions"],
        detached_reactive_collisions=metadata[
            "detached_reactive_collisions"],
        front_events=metadata["front_events"],
        landed_photoresist_count=metadata.get(
            "landed_photoresist_count", 0.0))


def _save_local_checkpoint(output, label, config_sha256, checkpoint):
    if not isinstance(checkpoint, HwangGiapisLocalNotchCheckpoint2D):
        raise TypeError("local checkpoint writer received an invalid state")
    arrays_path = output / f"{label}_checkpoint.npz"
    arrays_temporary = output / f"{label}_checkpoint.tmp.npz"
    np.savez_compressed(
        arrays_temporary,
        poly_cell=checkpoint.poly_cell,
        reactive_collision_inventory=(
            checkpoint.reactive_collision_inventory),
        cumulative_oxide_ion_count=(
            checkpoint.cumulative_oxide_ion_count))
    arrays_temporary.replace(arrays_path)
    metadata = {
        "config_sha256": config_sha256,
        "front_events": checkpoint.front_events,
        "remaining_campaign_fraction": (
            checkpoint.remaining_campaign_fraction),
        "direct_reactive_collisions": (
            checkpoint.direct_reactive_collisions),
        "scattered_reactive_collisions": (
            checkpoint.scattered_reactive_collisions),
        "landed_poly_count": checkpoint.landed_poly_count,
        "landed_oxide_count": checkpoint.landed_oxide_count,
        "landed_photoresist_count": (
            checkpoint.landed_photoresist_count),
        "escaped_count": checkpoint.escaped_count,
        "threshold_removed_cell_count": (
            checkpoint.threshold_removed_cell_count),
        "detached_cell_count": checkpoint.detached_cell_count,
        "threshold_removed_reactive_collisions": (
            checkpoint.threshold_removed_reactive_collisions),
        "detached_reactive_collisions": (
            checkpoint.detached_reactive_collisions),
        "timestamp_unix_s": time.time(),
    }
    metadata_path = output / f"{label}_checkpoint.json"
    metadata_temporary = output / f"{label}_checkpoint.tmp.json"
    metadata_temporary.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    metadata_temporary.replace(metadata_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--global-samples", type=int, default=3000)
    parser.add_argument("--global-iterations", type=int, default=2000)
    parser.add_argument(
        "--global-relax", type=float, default=8.0,
        help=(
            "base voltage-space stochastic step; the decreasing gain changes "
            "the path, not the current-balance fixed point"))
    parser.add_argument(
        "--global-gain-exponent", type=float, default=0.75,
        help=(
            "Robbins--Monro gain exponent in (0.5, 1]; unlike the legacy "
            "finite gain, this removes the terminal random walk"))
    parser.add_argument(
        "--global-gain-offset", type=float, default=25.0,
        help="positive early-iteration offset for the decreasing gain")
    parser.add_argument(
        "--global-final-samples", type=int,
        help=(
            "independent final ion/electron audit population; defaults to "
            "four times --global-samples"))
    parser.add_argument(
        "--global-final-seed", type=int,
        help=(
            "dedicated nested final-audit seed; defaults to run seed + "
            "1,000,003"))
    parser.add_argument(
        "--global-warm-start", type=Path,
        help=(
            "continue from a prior global-boundary directory; the source "
            "artifacts are checksum-pinned and the stochastic gain age is "
            "continued rather than reset"))
    parser.add_argument(
        "--ion-tangential-temperature-eV", type=float,
        help=(
            "optional lower-boundary effective transverse scale; the "
            "source-paper 0.5 eV value remains the default"))
    parser.add_argument(
        "--trajectory-fixed-dt", type=float, default=2.5e-4)
    parser.add_argument(
        "--trajectory-max-steps", type=int, default=8192,
        help=(
            "charged-particle trajectory horizon in fixed steps; double this "
            "when halving --trajectory-fixed-dt to preserve physical time"))
    parser.add_argument(
        "--maximum-front-events", type=int, default=2000)
    parser.add_argument(
        "--include-scatter-off", action="store_true",
        help="also run the source mechanism with SiO2 forward scattering disabled")
    parser.add_argument(
        "--include-local-charging-off", action="store_true",
        help=(
            "also run the open-area-only control with exposed-SiO2 charging "
            "and forward scattering disabled"))
    parser.add_argument(
        "--force-global", action="store_true",
        help="discard the matching cached global boundary and recompute it")
    parser.add_argument(
        "--global-only", action="store_true",
        help="build or verify the cached global boundary, then stop")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    observations = load_hwang_giapis_1997_fig13_profile(EVIDENCE)
    config = _global_config(args)
    global_result = None if args.force_global else _load_global_boundary(
        output, config)
    if global_result is None:
        print("global: solving 20 nm edge/open-area charging boundary", flush=True)

        def global_progress(record):
            heartbeat = {
                "phase": "global_charging",
                "config_sha256": _stable_hash(config),
                "timestamp_unix_s": time.time(),
                **record,
            }
            heartbeat_path = output / "global_heartbeat.json"
            heartbeat_temporary = output / "global_heartbeat.tmp.json"
            heartbeat_temporary.write_text(
                json.dumps(heartbeat, indent=2) + "\n",
                encoding="utf-8")
            heartbeat_temporary.replace(heartbeat_path)

        start = time.perf_counter()
        global_result = solve_edge_array_charging(
            **_global_solver_kwargs(
                config, progress_callback=global_progress))
        global_elapsed = time.perf_counter() - start
        _save_global_boundary(output, config, global_result, global_elapsed)
        global_metadata = json.loads(
            (output / "global_boundary.json").read_text(encoding="utf-8"))
        print(f"global: complete in {global_elapsed:.1f} s", flush=True)
    else:
        global_metadata = global_result.pop("_cached_metadata")
        print("global: reused checksum-matched cached boundary", flush=True)
    boundary = hwang_giapis_local_boundary_from_edge_array_result(global_result)
    if args.global_only:
        print(
            "global: local boundary contains "
            f"{boundary.target_event_count} weighted target-wall entries "
            f"from {boundary.source_particle_count} final ion samples",
            flush=True)
        return
    results = {}
    profile_arrays = {}
    profile_modes = [("event_driven_scatter_on", True, True)]
    if args.include_scatter_off or args.include_local_charging_off:
        profile_modes.append(("event_driven_scatter_off", False, True))
    if args.include_local_charging_off:
        profile_modes.append((
            "event_driven_local_charging_off", False, False))
    for label, scatter, exposed_oxide_charging in profile_modes:
        local_config = {
            "global_config_sha256": global_metadata["config_sha256"],
            "integrator": "event_driven_next_50_collision_threshold",
            "integrator_revision": (
                "event-driven-v4-gas-side-entry-connected-material-ledger"),
            "trajectory_fixed_dt": float(args.trajectory_fixed_dt),
            "trajectory_max_steps": int(args.trajectory_max_steps),
            "maximum_front_events": int(args.maximum_front_events),
            "include_forward_scatter": bool(scatter),
            "reactive_collisions_per_cell": 50.0,
            "cell_size_um": 0.005,
            "published_fluence_ions_per_0p5um": 18.7e6,
            "detached_fragment_policy": (
                "remove_zero_face_connected_components"),
        }
        if not exposed_oxide_charging:
            local_config["include_exposed_oxide_charging"] = False
        config_sha256 = _stable_hash(local_config)
        completed = _completed_local_run(output, label, config_sha256)
        if completed is None:
            print(f"local: {label}", flush=True)
            initial_checkpoint = _load_local_checkpoint(
                output, label, config_sha256)
            initial_event = (
                -1 if initial_checkpoint is None
                else initial_checkpoint.front_events)
            if initial_checkpoint is not None:
                print(
                    f"local: resuming event={initial_event} "
                    f"fluence={100.0 * (1.0 - initial_checkpoint.remaining_campaign_fraction):.2f}%",
                    flush=True)
            last_reported = {"event": initial_event - 5}
            last_saved = {"event": initial_event}

            def progress(front_event, campaign_fraction, maximum_depth_um):
                if (front_event == 0
                        or front_event >= last_reported["event"] + 5
                        or campaign_fraction >= 1.0 - 1e-12):
                    last_reported["event"] = front_event
                    heartbeat = {
                        "label": label,
                        "front_event": front_event,
                        "campaign_fraction": campaign_fraction,
                        "maximum_depth_um": maximum_depth_um,
                        "timestamp_unix_s": time.time(),
                    }
                    (output / "heartbeat.json").write_text(
                        json.dumps(heartbeat, indent=2) + "\n",
                        encoding="utf-8")
                    print(
                        f"  event={front_event} "
                        f"fluence={100.0 * campaign_fraction:.2f}% "
                        f"max_depth={maximum_depth_um:.4f} um",
                        flush=True)

            def checkpoint(state):
                if (state.front_events >= last_saved["event"] + 25
                        or state.remaining_campaign_fraction <= 2e-14):
                    _save_local_checkpoint(
                        output, label, config_sha256, state)
                    last_saved["event"] = state.front_events

            start = time.perf_counter()
            result = evolve_hwang_giapis_local_notch_event_driven_2d(
                boundary.entries, boundary.sidewall_potential_v,
                poly_potential_v=boundary.poly_potential_v,
                trajectory_fixed_dt=float(args.trajectory_fixed_dt),
                trajectory_max_steps=int(args.trajectory_max_steps),
                include_forward_scatter=scatter,
                include_exposed_oxide_charging=exposed_oxide_charging,
                maximum_front_events=int(args.maximum_front_events),
                initial_checkpoint=initial_checkpoint,
                progress_callback=progress,
                checkpoint_callback=checkpoint)
            elapsed = time.perf_counter() - start
            score = score_hwang_giapis_1997_fig13_profile(
                observations, result.notch_depth_by_height_um)
            completed = _save_local_run(
                output, label, local_config, result, score, elapsed)
            completed["_depth_um"] = result.notch_depth_by_height_um
            print(
                f"local: {label} max={result.maximum_notch_depth_um:.4f} um "
                f"RMSE={score.rmse_um:.4f} um in {elapsed:.1f} s",
                flush=True)
        else:
            print(f"local: reused {label}", flush=True)
        profile_arrays[f"{label}_depth_um"] = completed.pop("_depth_um")
        results[label] = completed
    experiment_height = np.asarray([
        item.height_above_oxide_um for item in observations])
    experiment_depth = np.asarray([
        item.notch_depth_um for item in observations])
    np.savez_compressed(
        output / "profile_results.npz",
        experimental_height_um=experiment_height,
        experimental_depth_um=experiment_depth,
        local_cell_center_height_um=(
            (np.arange(60) + 0.5) * 0.005),
        **profile_arrays)
    causal_mechanism_audit = {
        "status": "not_run",
        "source_reference": {
            "local_oxide_charging_reactive_collision_gain_fraction": 0.55,
            "scattering_reactive_collision_gain_fraction": 0.086,
            "citation": (
                "Hwang & Giapis 1997, Sec. IV D, comparison of "
                "Figure 14 profiles (a), (b), and (c)"),
        },
    }
    required_controls = {
        "event_driven_scatter_on",
        "event_driven_scatter_off",
        "event_driven_local_charging_off",
    }
    if required_controls.issubset(results):
        full = results["event_driven_scatter_on"]
        no_scatter = results["event_driven_scatter_off"]
        open_area_only = results["event_driven_local_charging_off"]
        full_reactive = (
            full["direct_reactive_collisions"]
            + full["scattered_reactive_collisions"])
        no_scatter_reactive = (
            no_scatter["direct_reactive_collisions"]
            + no_scatter["scattered_reactive_collisions"])
        open_area_reactive = (
            open_area_only["direct_reactive_collisions"]
            + open_area_only["scattered_reactive_collisions"])
        local_gain = (
            np.nan if open_area_reactive == 0.0
            else no_scatter_reactive / open_area_reactive - 1.0)
        scatter_gain = (
            np.nan if no_scatter_reactive == 0.0
            else full_reactive / no_scatter_reactive - 1.0)
        causal_mechanism_audit = {
            "status": "complete",
            "full_reactive_collisions": full_reactive,
            "scatter_off_reactive_collisions": no_scatter_reactive,
            "local_charging_off_reactive_collisions": open_area_reactive,
            "local_oxide_charging_reactive_collision_gain_fraction": (
                local_gain),
            "scattering_reactive_collision_gain_fraction": scatter_gain,
            "maximum_depth_ordering_um": {
                "full": full["maximum_notch_depth_um"],
                "scatter_off": no_scatter["maximum_notch_depth_um"],
                "local_charging_off": (
                    open_area_only["maximum_notch_depth_um"]),
            },
            "expected_depth_ordering_pass": bool(
                full["maximum_notch_depth_um"]
                >= no_scatter["maximum_notch_depth_um"]
                >= open_area_only["maximum_notch_depth_um"]),
            "source_reference": {
                "local_oxide_charging_reactive_collision_gain_fraction": 0.55,
                "scattering_reactive_collision_gain_fraction": 0.086,
                "citation": (
                    "Hwang & Giapis 1997, Sec. IV D, comparison of "
                    "Figure 14 profiles (a), (b), and (c)"),
            },
        }
    audit = {
        "schema": AUDIT_SCHEMA,
        "campaign": "Hwang--Giapis 1997 Figure 13 source-scale notch replay",
        "claim_status": (
            "QUANTITATIVE_DEVELOPMENT_REPLAY_ONLY: no target fitting; "
            "measurement uncertainty unreported; contour not held out"),
        "git_revision": _git_revision(),
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "experimental_evidence": {
            "path": str(EVIDENCE.relative_to(ROOT)),
            "sha256": HWANG_GIAPIS_1997_FIG13_PROFILE_SHA256,
            "point_count": len(observations),
            "digitization_uncertainty_um": 0.005,
            "measurement_uncertainty": "not reported",
        },
        "global_boundary": global_metadata,
        "global_to_local_reduction": {
            "target_event_count": boundary.target_event_count,
            "source_particle_count": boundary.source_particle_count,
            "weighted_launched_count": boundary.entries.launched_count,
            "poly_potential_v": boundary.poly_potential_v,
            "sidewall_potential_min_v": float(
                np.min(boundary.sidewall_potential_v)),
            "sidewall_potential_max_v": float(
                np.max(boundary.sidewall_potential_v)),
            "provenance": dict(boundary.provenance),
        },
        "profile_runs": results,
        "causal_mechanism_audit": causal_mechanism_audit,
        "discarded_equal_batch_diagnostic": {
            "batches_64_maximum_notch_depth_um": 0.280,
            "batches_64_rmse_um": 0.0384,
            "batches_128_maximum_notch_depth_um": 0.385,
            "batches_128_rmse_um": 0.0514,
            "conclusion": (
                "equal-batch front integration is not refinement-convergent "
                "and is excluded from validation scoring"),
        },
        "no_profile_parameters_fitted": True,
    }
    (output / "audit.json").write_text(
        json.dumps(_json_value(audit), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output / 'audit.json'}", flush=True)


if __name__ == "__main__":
    main()
