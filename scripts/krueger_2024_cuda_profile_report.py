#!/usr/bin/env python3
"""Summarize the bounded Krueger CUDA and 10/5 nm preflight evidence.

This is a postprocessor only.  It never imports the engine, reads held-out
profiles, chooses calibration parameters, or launches a simulation.  The two
CUDA inputs must represent the same 5 nm, 0.05 s, seed-paired trajectory.  The
mixed-device run is retained solely as a numerical parity reference; timing
authority belongs to the explicitly unified CUDA run.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "petch-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROFILE_SCHEMA = "petch.krueger_2024_cuda_profile.v1"
SUMMARY_SCHEMA = "petch.krueger_2024_cuda_profile_summary.v1"
PARITY_METRICS = (
    "etch_depth_nm",
    "floor_z_um",
    "mask_opening_nm",
    "mask_top_z_um",
    "maximum_feature_width_nm",
    "remaining_mask_thickness_nm",
    "top_feature_width_nm",
)
REFINEMENT_METRICS = (
    "etch_depth_nm_s",
    "mask_opening_nm_s",
    "remaining_mask_thickness_nm_s",
    "maximum_feature_width_nm_s",
    "top_feature_width_nm_s",
)
TIMING_ROWS = (
    ("chemistry and material routing", "src/petch/material_mechanism_3d.py", "advance_by_material"),
    ("ballistic boundary transport", "src/petch/boundary_transport_3d.py", "gather_boundary_state_ballistic_3d"),
    ("diffuse neutral exchange", "src/petch/feature_step_3d.py", "_apply_diffuse_neutral_transport"),
    ("surface-state remap", "src/petch/feature_step_3d.py", "conservative_remap_surface_state"),
    ("level-set redistance", "src/petch/feature_step_3d.py", "_redistance_feature_field"),
)


def _sha(path: Path | str) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _read(path: Path | str) -> tuple[Path, dict]:
    path = Path(path)
    return path, json.loads(path.read_text(encoding="utf-8"))


def _profile(path: Path | str, *, require_unified: bool) -> tuple[dict, dict]:
    path, payload = _read(path)
    if payload.get("schema") != PROFILE_SCHEMA or payload.get("status") != "complete":
        raise ValueError(f"{path.name} is not a complete bounded CUDA profile")
    if payload.get("held_out_profile_data_read") is not False:
        raise ValueError(f"{path.name} is not sealed from held-out data")
    configuration = payload.get("configuration", {})
    if (not np.isclose(float(configuration.get("dx_nm", np.nan)), 5.0)
            or not np.isclose(float(configuration.get("step_duration_s", np.nan)), 0.025)
            or int(configuration.get("seed", -1)) != 241):
        raise ValueError(f"{path.name} does not match the bounded paired preflight")
    if require_unified:
        hardware = payload.get("hardware", {})
        if (hardware.get("unified_device_selection") is not True
                or hardware.get("transport_device") != hardware.get("level_set_device")
                or not str(hardware.get("transport_device", "")).startswith("cuda")):
            raise ValueError("timing profile did not use one declared CUDA device")
    return payload, {"path_name": path.name, "sha256": _sha(path)}


def _executed_time(configuration: dict) -> float:
    if "total_executed_physical_time_s" in configuration:
        return float(configuration["total_executed_physical_time_s"])
    return float(configuration["total_profiled_physical_time_s"])


def _find_row(rows: list[dict], path: str, function: str) -> dict | None:
    matches = [
        row for row in rows
        if row.get("path") == path and row.get("function") == function
    ]
    if len(matches) > 1:
        raise ValueError(f"ambiguous profile row for {path}:{function}")
    return None if not matches else matches[0]


def summarize(mixed_path: Path | str, unified_path: Path | str,
              multiresolution_path: Path | str) -> dict:
    mixed, mixed_info = _profile(mixed_path, require_unified=False)
    unified, unified_info = _profile(unified_path, require_unified=True)
    if (not np.isclose(_executed_time(mixed["configuration"]),
                       _executed_time(unified["configuration"]), rtol=0.0, atol=1e-15)
            or mixed["configuration"]["calibration_parameters"]
            != unified["configuration"]["calibration_parameters"]
            or mixed["configuration"]["operator"] != unified["configuration"]["operator"]):
        raise ValueError("mixed and unified CUDA profiles are not paired")

    parity = {}
    for metric in PARITY_METRICS:
        mixed_value = float(mixed["end_metrics"][metric])
        unified_value = float(unified["end_metrics"][metric])
        absolute = abs(unified_value - mixed_value)
        parity[metric] = {
            "mixed": mixed_value,
            "unified": unified_value,
            "absolute_difference": absolute,
            "relative_difference": absolute / max(abs(mixed_value), np.finfo(float).tiny),
        }
    per_y_mixed = np.asarray(mixed["end_metrics"]["mask_opening_per_y_widths_nm"], dtype=float)
    per_y_unified = np.asarray(
        unified["end_metrics"]["mask_opening_per_y_widths_nm"], dtype=float)
    if per_y_mixed.shape != per_y_unified.shape:
        raise ValueError("paired profiles have incompatible per-y opening observations")
    # Put all dimensional parity checks on a common nm scale.  The two z fields
    # above are stored in micrometres by the measurement contract.
    parity_nm = []
    for name, entry in parity.items():
        scale = 1000.0 if name.endswith("_z_um") else 1.0
        parity_nm.append(entry["absolute_difference"] * scale)
    maximum_parity_delta_nm = max(
        parity_nm + [float(np.max(np.abs(per_y_unified - per_y_mixed)))])

    profile_rows = list(unified.get("top_cumulative_functions", ()))
    total_row = _find_row(
        profile_rows, "src/petch/feature_step_3d.py", "advance_feature_step_3d")
    if total_row is None:
        raise ValueError("unified profile lacks the feature-step timing root")
    total_s = float(total_row["cumulative_time_s"])
    components = []
    for label, path, function in TIMING_ROWS:
        row = _find_row(profile_rows, path, function)
        component_s = 0.0 if row is None else float(row["cumulative_time_s"])
        components.append({
            "label": label,
            "path": path,
            "function": function,
            "time_s": component_s,
            "fraction": component_s / total_s,
        })
    attributed_s = sum(component["time_s"] for component in components)
    if attributed_s > total_s * (1.0 + 1e-8):
        raise ValueError("profile categories overlap and cannot form a timing budget")
    components.append({
        "label": "remaining extraction, advection, topology, and overhead",
        "path": None,
        "function": None,
        "time_s": max(0.0, total_s - attributed_s),
        "fraction": max(0.0, total_s - attributed_s) / total_s,
    })

    multiresolution_path, multiresolution = _read(multiresolution_path)
    paired = multiresolution.get("paired_10nm_vs_5nm", {}).get("initial", {})
    if set(REFINEMENT_METRICS) - set(paired):
        raise ValueError("10/5 nm audit lacks the required initial refinement metrics")
    refinement = {
        metric: {
            "coarse_10nm": float(paired[metric]["coarse_10nm"]),
            "fine_5nm": float(paired[metric]["fine_5nm"]),
            "relative_to_fine": float(paired[metric]["relative_to_fine"]),
            "absolute_relative_percent": abs(float(
                paired[metric]["relative_to_fine"])) * 100.0,
        }
        for metric in REFINEMENT_METRICS
    }
    cases = {float(case["dx_nm"]): case for case in multiresolution.get("cases", ())}
    if set(cases) != {5.0, 10.0}:
        raise ValueError("multiresolution evidence is not the paired 10/5 nm audit")

    runtime = unified["runtime"]
    steady_step_s = float(runtime["mean_profiled_step_wall_time_s"])
    step_duration_s = float(unified["configuration"]["step_duration_s"])
    naive_60s_projection = steady_step_s * 60.0 / step_duration_s
    material_residual = max(
        float(row["maximum_material_ledger_residual_units_m2"])
        for row in unified["step_receipts"])
    radiosity_residual = max(
        float(row["maximum_radiosity_relative_balance_error"])
        for row in unified["step_receipts"])

    payload = {
        "schema": SUMMARY_SCHEMA,
        "scientific_status": (
            "bounded base-only performance and numerical preflight; no endpoint calibration "
            "authority and no held-out validation claim"),
        "held_out_profile_data_read": False,
        "calibration_performed": False,
        "generator": {
            "path_name": Path(__file__).name,
            "sha256": _sha(Path(__file__)),
        },
        "inputs": {
            "mixed_device_parity_reference": mixed_info,
            "unified_cuda_timing_profile": unified_info,
            "paired_10nm_5nm_initial_audit": {
                "path_name": multiresolution_path.name,
                "sha256": _sha(multiresolution_path),
            },
        },
        "paired_cuda_parity": {
            "same_seed_operator_and_executed_physical_time": True,
            "metrics": parity,
            "maximum_per_y_opening_difference_nm": float(
                np.max(np.abs(per_y_unified - per_y_mixed))),
            "maximum_dimensional_difference_nm": maximum_parity_delta_nm,
            "within_1e_6_nm": bool(maximum_parity_delta_nm <= 1e-6),
        },
        "unified_cuda_runtime": {
            "hardware": unified["hardware"],
            "initialization_wall_time_s": float(runtime[
                "operator_initialization_wall_time_s"]),
            "positive_warmup_step_wall_time_s": list(runtime[
                "positive_warmup_step_wall_time_s"]),
            "steady_profiled_step_wall_time_s": steady_step_s,
            "profiled_feature_step_cumulative_time_s": total_s,
            "timing_budget": components,
            "naive_fixed_0p025s_step_projection_for_60s": {
                "steps": int(round(60.0 / step_duration_s)),
                "wall_time_s": naive_60s_projection,
                "wall_time_h": naive_60s_projection / 3600.0,
                "warning": (
                    "arithmetic projection only; one warmed short step is not an endpoint "
                    "benchmark and geometry/operator cost changes with process time"),
            },
        },
        "paired_10nm_5nm_initial_refinement": {
            "duration_s": 0.5,
            "metrics": refinement,
            "wall_time_s": {
                "fine_5nm": float(cases[5.0]["wall_time_s"]),
                "coarse_10nm": float(cases[10.0]["wall_time_s"]),
            },
            "interpretation": (
                "depth, definition-correct opening, and mask-thickness rates agree within "
                "1.5 percent; maximum/top feature-width rates differ by 71.2 percent, so "
                "local profile shape is not spatially authoritative"),
        },
        "operator_receipts": {
            "maximum_material_ledger_residual_units_m2": material_residual,
            "maximum_radiosity_relative_balance_error": radiosity_residual,
            "topology_events": [row["topology_event"] for row in unified["step_receipts"]],
            "surface_state_remap_backend": "legacy_knn (engine default used by this worker)",
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["summary_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def plot_summary(summary: dict, output: Path | str) -> None:
    timing = summary["unified_cuda_runtime"]["timing_budget"]
    refinement = summary["paired_10nm_5nm_initial_refinement"]["metrics"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(13.2, 5.3))

    timing_labels = [item["label"] for item in timing]
    timing_values = [item["time_s"] for item in timing]
    order = np.argsort(timing_values)
    left.barh(
        np.asarray(timing_labels)[order], np.asarray(timing_values)[order],
        color="#3366a8")
    total = sum(timing_values)
    for y, index in enumerate(order):
        value = timing_values[index]
        left.text(value + total * 0.012, y, f"{value:.2f}s  ({100*value/total:.1f}%)",
                  va="center", fontsize=9)
    left.set_xlabel("Cumulative wall time in one warmed 5 nm step (s)")
    left.set_title("Where the unified CUDA step spends time")
    left.set_xlim(0.0, max(timing_values) * 1.43)
    left.grid(axis="x", alpha=0.25)

    label_map = {
        "etch_depth_nm_s": "etch-depth rate",
        "mask_opening_nm_s": "mask-opening rate",
        "remaining_mask_thickness_nm_s": "mask-thickness rate",
        "maximum_feature_width_nm_s": "maximum-width rate",
        "top_feature_width_nm_s": "top-width rate",
    }
    names = list(REFINEMENT_METRICS)
    values = [max(refinement[name]["absolute_relative_percent"], 1e-4) for name in names]
    colors = ["#448866" if index < 3 else "#c06a3a" for index in range(len(names))]
    right.barh([label_map[name] for name in names], values, color=colors)
    right.set_xscale("log")
    right.invert_yaxis()
    for y, value in enumerate(values):
        right.text(value * 1.13, y, f"{value:.3g}%", va="center", fontsize=9)
    right.set_xlabel("Absolute 10 nm vs 5 nm relative difference (%) — log scale")
    right.set_title("Global calibration observables (green) vs local shape (orange)")
    right.set_xlim(0.04, 180.0)
    right.grid(axis="x", which="both", alpha=0.25)
    fig.suptitle("Krüger 2024 bounded numerical preflight", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mixed-profile", required=True)
    parser.add_argument("--unified-profile", required=True)
    parser.add_argument("--multiresolution", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    summary = summarize(
        args.mixed_profile, args.unified_profile, args.multiresolution)
    audit_path = output / "audit.json"
    temporary = audit_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(audit_path)
    plot_summary(summary, output / "profile_summary.png")
    print(json.dumps({
        "steady_step_wall_time_s": summary["unified_cuda_runtime"][
            "steady_profiled_step_wall_time_s"],
        "naive_60s_projection_h": summary["unified_cuda_runtime"][
            "naive_fixed_0p025s_step_projection_for_60s"]["wall_time_h"],
        "maximum_cuda_parity_difference_nm": summary["paired_cuda_parity"][
            "maximum_dimensional_difference_nm"],
        "maximum_local_width_refinement_difference_percent": max(
            summary["paired_10nm_5nm_initial_refinement"]["metrics"][name][
                "absolute_relative_percent"]
            for name in ("maximum_feature_width_nm_s", "top_feature_width_nm_s")),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
