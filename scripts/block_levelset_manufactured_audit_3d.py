#!/usr/bin/env python3
"""Bounded storage/work audit for fixed-dx block-sparse Krueger geometry.

The audit uses analytic base geometry only.  It reads no experimental outcome,
does not evolve a profile, and cannot authorize AMR or calibration.  Its job is
to decide whether fixed-spacing sparse-volume kernels can plausibly meet the
predeclared 3x geometry-memory and 2x level-set-work gates before those kernels
are built.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from petch.block_levelset_3d import build_block_sparse_levelset_3d
from petch.feature_step_3d import make_rectangular_trench_geometry_3d


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "block_levelset_manufactured_audit_3d"
DX_NM = (10.0, 5.0)
DEPTH_UM = (0.0, 0.45, 0.90)
BAND_WIDTH_CELLS = (4, 8, 12)
REQUIRED_BAND_WIDTH_CELLS = 8
MEMORY_GATE = 3.0
WORK_GATE = 2.0
MAXIMUM_CASES = 72


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(
        payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def candidate_block_shapes(shape):
    interval = tuple(int(value) - 1 for value in shape)
    x, y, z = interval
    candidates = (
        (min(8, x), y, min(8, z)),
        (min(8, x), y, min(16, z)),
        (x, y, min(16, z)),
        (x, y, min(32, z)),
    )
    return tuple(dict.fromkeys(candidates))


def _geometry(dx_nm, depth_um):
    return make_rectangular_trench_geometry_3d(
        cell_width=0.13, cell_length=0.02, domain_height=2.8,
        dx=float(dx_nm) / 1000.0, opening_width=0.09,
        mask_thickness=0.85, substrate_top=1.8,
        etched_depth=float(depth_um))


def decide(records):
    required = [
        record for record in records
        if record["dx_nm"] == 5.0
        and record["etched_depth_um"] == 0.9
        and record["band_width_cells"] == REQUIRED_BAND_WIDTH_CELLS
        and record["integrity_pass"]]
    if not required:
        return {
            "status": "blocked_missing_required_deep_5nm_cases",
            "pass": False,
        }
    best_core = max(required, key=lambda item: item["core_memory_reduction"])
    best_indexed = max(
        required, key=lambda item: item["indexed_halo_only_memory_reduction"])
    best_work = max(required, key=lambda item: item["unique_node_work_upper_bound"])
    passed = bool(
        best_core["core_memory_reduction"] >= MEMORY_GATE
        and best_indexed["indexed_halo_only_memory_reduction"] >= MEMORY_GATE
        and best_work["unique_node_work_upper_bound"] >= WORK_GATE)
    return {
        "status": (
            "fixed_dx_sparse_candidate_pass"
            if passed else "fixed_dx_sparse_no_go_for_krueger"),
        "pass": passed,
        "gates": {
            "minimum_core_memory_reduction": MEMORY_GATE,
            "minimum_indexed_packed_memory_reduction": MEMORY_GATE,
            "minimum_unique_node_work_upper_bound": WORK_GATE,
        },
        "best_required_band_deep_5nm": {
            "core": best_core,
            "indexed": best_indexed,
            "work_upper_bound": best_work,
        },
        "scientific_action": (
            "implement sparse fixed-dx evolution kernels"
            if passed else
            "do not wire fixed-dx sparse evolution for Krueger; retain the verified block "
            "authority as AMR infrastructure and target surface/transport symmetry or adaptive "
            "surface work, which dominates measured runtime"),
    }


def _best_series(records, name):
    output = {}
    for dx_nm in DX_NM:
        for depth_um in DEPTH_UM:
            selected = [
                item for item in records
                if item["dx_nm"] == dx_nm
                and item["etched_depth_um"] == depth_um
                and item["band_width_cells"] == REQUIRED_BAND_WIDTH_CELLS]
            output[(dx_nm, depth_um)] = max(float(item[name]) for item in selected)
    return output


def render_plot(records, output):
    core = _best_series(records, "core_memory_reduction")
    indexed = _best_series(records, "indexed_halo_only_memory_reduction")
    work = _best_series(records, "unique_node_work_upper_bound")
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    figure.suptitle("Fixed-dx sparse-volume gate on analytic Krueger geometry")
    colors = {10.0: "#2878B5", 5.0: "#D95F02"}
    for dx_nm in DX_NM:
        depth = np.asarray(DEPTH_UM)
        axes[0].plot(
            depth, [core[(dx_nm, value)] for value in depth], marker="o",
            color=colors[dx_nm], label=f"{dx_nm:g} nm core")
        axes[0].plot(
            depth, [indexed[(dx_nm, value)] for value in depth], marker="s",
            linestyle="--", color=colors[dx_nm],
            label=f"{dx_nm:g} nm indexed halo")
        axes[1].plot(
            depth, [work[(dx_nm, value)] for value in depth], marker="o",
            color=colors[dx_nm], label=f"{dx_nm:g} nm")
    axes[0].axhline(MEMORY_GATE, color="#B23A48", linestyle=":",
                    label="3× gate")
    axes[0].set_title("Best memory reduction, safe 8-cell band")
    axes[0].set_xlabel("Analytic etched depth (µm)")
    axes[0].set_ylabel("Dense bytes / sparse bytes")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].axhline(WORK_GATE, color="#B23A48", linestyle=":",
                    label="2× gate")
    axes[1].set_title("Optimistic unique-node work ceiling")
    axes[1].set_xlabel("Analytic etched depth (µm)")
    axes[1].set_ylabel("Dense nodes / unique packed nodes")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    figure.text(
        0.5, 0.01,
        "As the trench deepens, its surface band occupies most of this narrow 3-D cell; "
        "volume sparsity vanishes before transport/chemistry cost is touched.",
        ha="center", fontsize=9)
    figure.tight_layout(rect=(0.02, 0.08, 0.98, 0.93))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def run(args):
    started = perf_counter()
    records = []
    status = "complete"
    for dx_nm in DX_NM:
        for depth_um in DEPTH_UM:
            geometry = _geometry(dx_nm, depth_um)
            shapes = candidate_block_shapes(geometry.phi.shape)
            for band_width in BAND_WIDTH_CELLS:
                for block_shape in shapes:
                    if perf_counter() - started >= float(args.max_wall_s):
                        status = "wall_budget_checkpoint"
                        break
                    case_started = perf_counter()
                    sparse = build_block_sparse_levelset_3d(
                        geometry, block_cell_shape=block_shape,
                        band_width_cells=band_width, periodic_axes=(0, 1))
                    receipt = sparse.storage_receipt(halo_cells=1)
                    combined, owner, material = sparse.to_truncated_dense()
                    near = np.abs(geometry.phi) <= 2.0 * geometry.dx
                    maximum_band_error = float(np.max(
                        np.abs(combined[near] - geometry.phi[near]), initial=0.0))
                    material_error = max(float(np.max(np.abs(
                        material[material_id][near]
                        - geometry.material_levelsets[material_id][near]), initial=0.0))
                        for material_id in sparse.material_ids)
                    tolerance = max(
                        sparse.periodic_endpoint_max_abs_difference,
                        1.0e-12 * geometry.dx)
                    integrity_pass = bool(
                        maximum_band_error <= tolerance
                        and material_error <= tolerance
                        and np.array_equal(owner[near], geometry.material_id[near]))
                    records.append({
                        "dx_nm": float(dx_nm),
                        "etched_depth_um": float(depth_um),
                        "band_width_cells": int(band_width),
                        "block_cell_shape": list(block_shape),
                        "wall_time_s": perf_counter() - case_started,
                        "integrity_pass": integrity_pass,
                        "maximum_exact_band_error_mesh_units": maximum_band_error,
                        "maximum_material_band_error_mesh_units": material_error,
                        **receipt,
                        "unique_node_work_upper_bound": float(
                            receipt["dense_node_count"]
                            / receipt["indexed_unique_node_count"]),
                    })
                if status != "complete":
                    break
            if status != "complete":
                break
        if status != "complete":
            break
    decision = decide(records) if status == "complete" else {
        "status": "incomplete_wall_budget", "pass": False}
    revision = subprocess.run(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, check=True,
        capture_output=True, text=True).stdout.strip()
    payload = {
        "schema": "petch.block_levelset_manufactured_audit_3d.v1",
        "status": decision["status"],
        "scientific_scope": (
            "analytic geometry and numerical storage only; no evolution, calibration, "
            "experimental outcome, or held-out data"),
        "configuration": {
            "dx_nm": list(DX_NM),
            "etched_depth_um": list(DEPTH_UM),
            "band_width_cells": list(BAND_WIDTH_CELLS),
            "required_band_width_cells": REQUIRED_BAND_WIDTH_CELLS,
            "halo_cells": 1,
            "maximum_cases": MAXIMUM_CASES,
            "max_wall_s": float(args.max_wall_s),
        },
        "generator": {
            "path_name": Path(__file__).name,
            "sha256": _sha(__file__),
            "git_revision": revision,
        },
        "case_count": len(records),
        "wall_time_s": perf_counter() - started,
        "cases": records,
        "decision": decision,
    }
    output = Path(args.output_dir)
    _write_json(output / "audit.json", payload)
    if status == "complete":
        render_plot(records, output / "memory_gate.png")
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-wall-s", type=float, default=120.0)
    args = parser.parse_args(argv)
    case_count = len(DX_NM) * len(DEPTH_UM) * len(BAND_WIDTH_CELLS) * 4
    if case_count > MAXIMUM_CASES:
        parser.error("declared sparse audit exceeds its hard case budget")
    if args.max_wall_s <= 0.0 or args.max_wall_s > 180.0:
        parser.error("max-wall-s must be in (0, 180]")
    return args


if __name__ == "__main__":
    run(parse_args())
