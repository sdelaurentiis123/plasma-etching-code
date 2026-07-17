#!/usr/bin/env python3
"""Bounded current-operator CUDA profile on the Krueger base geometry.

The already-complete 10/5 nm 0.5 s refinement pair is not repeated.  This
harness warms the exact initial 5 nm operator once, then profiles one to four
positive feature steps through the same ``advance_feature_step_3d`` path used
by the campaign.  It reads no held-out observations, performs no calibration,
and refuses more than 0.1 s of physical evolution.

Run the process under an external wall timeout as well; the internal wall gate
is checked between indivisible engine calls.
"""
from __future__ import annotations

import argparse
import cProfile
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import pstats
import subprocess
import sys
from time import perf_counter

import numpy as np
import warp as wp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import krueger_2024_multiresolution_audit as multires  # noqa: E402
import krueger_2024_trench_pilot as pilot  # noqa: E402
from petch.feature_step_3d import make_rectangular_trench_geometry_3d  # noqa: E402
from petch import threed  # noqa: E402


def _sha256(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _git_revision():
    return subprocess.run(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def _git_clean():
    return not subprocess.run(
        ("git", "status", "--porcelain"), cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout


def _display_path(filename):
    path = Path(filename)
    try:
        return str(path.resolve().relative_to(ROOT))
    except (OSError, ValueError):
        return path.name


def top_profile_rows(stats, limit=80):
    """Return deterministic cumulative-time rows from one ``pstats.Stats``."""
    rows = []
    ordered = sorted(
        stats.stats.items(),
        key=lambda item: (-float(item[1][3]), str(item[0])),
    )
    for (filename, line, function), (primitive, calls, self_s, cumulative_s, _) in (
            ordered[:int(limit)]):
        rows.append({
            "path": _display_path(filename),
            "line": int(line),
            "function": str(function),
            "primitive_calls": int(primitive),
            "total_calls": int(calls),
            "self_time_s": float(self_s),
            "cumulative_time_s": float(cumulative_s),
        })
    return rows


def _write_json(path, payload):
    target = Path(path)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(target)


def _maximum_radiosity_residual(result):
    return max((
        float(value["relative_balance_error"])
        for value in result.diagnostics["neutral_radiosity"].values()
    ), default=0.0)


def run(args):
    if not _git_clean():
        raise RuntimeError("CUDA profiler requires a clean checksum-bound worktree")
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    requested_device = wp.get_device(args.device)
    level_set_device = wp.get_device(threed.DEVICE)
    if requested_device.is_cuda != level_set_device.is_cuda:
        raise RuntimeError(
            "mixed engine devices are not a valid CUDA profile: "
            f"transport={requested_device}, level_set={level_set_device}; "
            "launch with PETCH_DEVICE set to the declared --device")
    dx = float(args.dx_nm) / 1000.0
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.13,
        cell_length=0.02,
        domain_height=2.8,
        dx=dx,
        opening_width=0.09,
        mask_thickness=0.85,
        substrate_top=1.8,
        etched_depth=0.0,
    )
    start_metrics = pilot.measure_krueger_metrics(geometry, substrate_top_um=1.8)
    start_metrics, initial_completion = multires._complete_analytic_initial_metrics(
        start_metrics)

    overall_started = perf_counter()
    initialization_started = perf_counter()
    frozen = multires._advance(
        geometry, None, None,
        duration_s=0.0,
        seed=int(args.seed),
        device=args.device,
        topology_policy="refuse",
        remap_backend=args.surface_state_remap_backend,
    )
    wp.synchronize_device(args.device)
    wp.synchronize_device(threed.DEVICE)
    initialization_wall = perf_counter() - initialization_started
    geometry = frozen.geometry
    state = frozen.next_surface_state
    fingerprint = frozen.next_surface_state_mesh_fingerprint

    positive_warmup_wall = []
    for index in range(int(args.positive_warmup_steps)):
        warmup_started = perf_counter()
        warmup = multires._advance(
            geometry,
            state,
            fingerprint,
            duration_s=float(args.step_duration_s),
            seed=int(args.seed) + index,
            device=args.device,
            topology_policy="refuse",
            remap_backend=args.surface_state_remap_backend,
        )
        wp.synchronize_device(args.device)
        wp.synchronize_device(threed.DEVICE)
        positive_warmup_wall.append(perf_counter() - warmup_started)
        geometry = warmup.geometry
        state = warmup.next_surface_state
        fingerprint = warmup.next_surface_state_mesh_fingerprint

    profiler = cProfile.Profile()
    step_receipts = []
    profiler.enable()
    profiled_started = perf_counter()
    for index in range(int(args.profile_steps)):
        if perf_counter() - overall_started >= float(args.max_wall_s):
            raise TimeoutError("CUDA profile reached its internal wall budget")
        step_started = perf_counter()
        result = multires._advance(
            geometry,
            state,
            fingerprint,
            duration_s=float(args.step_duration_s),
            seed=int(args.seed) + int(args.positive_warmup_steps) + index,
            device=args.device,
            topology_policy="refuse",
            remap_backend=args.surface_state_remap_backend,
        )
        wp.synchronize_device(args.device)
        wp.synchronize_device(threed.DEVICE)
        step_wall = perf_counter() - step_started
        step_receipts.append({
            "step": index + 1,
            "wall_time_s": step_wall,
            "active_face_count": int(result.active_face_area.size),
            "maximum_displacement_cells": float(
                result.diagnostics["max_displacement_mesh_units"] / geometry.dx),
            "maximum_material_ledger_residual_units_m2": float(
                pilot._maximum_ledger_residual(result.surface.material_exchange)),
            "maximum_radiosity_relative_balance_error": (
                _maximum_radiosity_residual(result)),
            "topology_event": result.diagnostics["topology_event"],
        })
        geometry = result.geometry
        state = result.next_surface_state
        fingerprint = result.next_surface_state_mesh_fingerprint
    wp.synchronize_device(args.device)
    wp.synchronize_device(threed.DEVICE)
    profiled_wall = perf_counter() - profiled_started
    profiler.disable()

    profile_path = output / "profile.pstats"
    profiler.dump_stats(profile_path)
    stats = pstats.Stats(profiler)
    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats(
        "cumulative").print_stats(int(args.text_rows))
    (output / "profile.txt").write_text(stream.getvalue(), encoding="utf-8")
    device = wp.get_device(args.device)
    end_metrics = pilot.measure_krueger_metrics(geometry, substrate_top_um=1.8)
    payload = {
        "schema": "petch.krueger_2024_cuda_profile.v1",
        "scope": "bounded_base_geometry_performance_profile_no_held_out_data",
        "status": "complete",
        "git_revision": _git_revision(),
        "worktree_clean": True,
        "source_sha256": {
            "profiler": _sha256(Path(__file__)),
            "multiresolution_worker": _sha256(Path(multires.__file__)),
        },
        "held_out_profile_data_read": False,
        "calibration_performed": False,
        "configuration": {
            "dx_nm": float(args.dx_nm),
            "profile_steps": int(args.profile_steps),
            "positive_warmup_steps": int(args.positive_warmup_steps),
            "step_duration_s": float(args.step_duration_s),
            "total_profiled_physical_time_s": (
                int(args.profile_steps) * float(args.step_duration_s)),
            "total_executed_physical_time_s": (
                (int(args.profile_steps) + int(args.positive_warmup_steps))
                * float(args.step_duration_s)),
            "seed": int(args.seed),
            "device": str(args.device),
            "operator": multires.OPERATOR,
            "surface_state_remap_backend": str(args.surface_state_remap_backend),
            "calibration_parameters": multires.CALIBRATION,
        },
        "hardware": {
            "warp_version": str(wp.config.version),
            "device": str(device),
            "device_name": str(getattr(device, "name", device)),
            "device_architecture": str(getattr(device, "arch", "unknown")),
            "transport_device": str(requested_device),
            "level_set_device": str(level_set_device),
            "unified_device_selection": True,
        },
        "runtime": {
            "operator_initialization_wall_time_s": initialization_wall,
            "positive_warmup_step_wall_time_s": positive_warmup_wall,
            "profiled_wall_time_s": profiled_wall,
            "mean_profiled_step_wall_time_s": profiled_wall / int(args.profile_steps),
            "overall_wall_time_s": perf_counter() - overall_started,
            "internal_wall_limit_s": float(args.max_wall_s),
        },
        "start_metrics": start_metrics,
        "initial_metric_completion": initial_completion,
        "end_metrics": end_metrics,
        "step_receipts": step_receipts,
        "top_cumulative_functions": top_profile_rows(stats),
        "artifacts": {
            "pstats": profile_path.name,
            "text": "profile.txt",
        },
    }
    _write_json(output / "audit.json", payload)
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results" / "krueger_2024_cuda_profile")
    parser.add_argument("--dx-nm", type=float, default=5.0)
    parser.add_argument("--profile-steps", type=int, default=2)
    parser.add_argument("--positive-warmup-steps", type=int, default=1)
    parser.add_argument("--step-duration-s", type=float, default=0.025)
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument(
        "--surface-state-remap-backend",
        choices=multires.REMAP_BACKENDS,
        default="legacy_knn")
    parser.add_argument("--max-wall-s", type=float, default=300.0)
    parser.add_argument("--text-rows", type=int, default=80)
    args = parser.parse_args(argv)
    total_physical = (
        (int(args.profile_steps) + int(args.positive_warmup_steps))
        * float(args.step_duration_s))
    if (not np.isfinite(args.dx_nm) or args.dx_nm <= 0.0
            or int(args.profile_steps) < 1 or int(args.profile_steps) > 4
            or int(args.positive_warmup_steps) < 1
            or int(args.positive_warmup_steps) > 2
            or not np.isfinite(args.step_duration_s) or args.step_duration_s <= 0.0
            or total_physical > 0.1
            or not np.isfinite(args.max_wall_s) or args.max_wall_s <= 0.0
            or int(args.text_rows) < 1 or int(args.text_rows) > 200):
        parser.error(
            "profiler requires 1--2 positive warm-up steps, 1--4 timed steps, "
            "<=0.1 s total physical time, and positive bounded controls")
    return args


def main(argv=None):
    payload = run(parse_args(argv))
    print(json.dumps({
        "status": payload["status"],
        "mean_profiled_step_wall_time_s": payload["runtime"][
            "mean_profiled_step_wall_time_s"],
        "output": payload["artifacts"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
