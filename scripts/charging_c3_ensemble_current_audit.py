#!/usr/bin/env python3
"""Combine independent exact current audits at one fixed C3 charge state."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def _atomic_npz(path: Path, **arrays) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-audits", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.current_audits) < 2:
        parser.error("at least two independent current audits are required")
    paths = [path.resolve() for path in args.current_audits]
    records = []
    for path in paths:
        with np.load(path) as audit:
            required = {
                "positive_face_current_density_a_m2",
                "negative_face_current_density_a_m2", "positive_current_node_a",
                "negative_current_node_a", "potential_before_v", "physical_face_area_m2",
                "patch_scales_m", "patch_group_by_scale"}
            if not required.issubset(audit.files):
                parser.error(f"incomplete current audit: {path}")
            records.append({name: np.asarray(audit[name]).copy() for name in required})
    reference = records[0]
    exact_names = ("physical_face_area_m2", "patch_scales_m", "patch_group_by_scale")
    for index, record in enumerate(records[1:], start=1):
        for name in exact_names:
            if not np.array_equal(record[name], reference[name]):
                parser.error(f"audit {index} disagrees on {name}")
        if record["potential_before_v"].shape != reference["potential_before_v"].shape:
            parser.error(f"audit {index} has a different state shape")
    potential = np.asarray([record["potential_before_v"] for record in records])
    potential_spread = float(np.max(np.ptp(potential, axis=0)))
    potential_scale = max(float(np.max(np.abs(potential[0]))), np.finfo(float).tiny)
    if potential_spread / potential_scale > 5e-12:
        parser.error("current audits do not score the same fixed potential state")

    positive_face = np.asarray([
        record["positive_face_current_density_a_m2"] for record in records])
    negative_face = np.asarray([
        record["negative_face_current_density_a_m2"] for record in records])
    positive_node = np.asarray([record["positive_current_node_a"] for record in records])
    negative_node = np.asarray([record["negative_current_node_a"] for record in records])
    count = len(records)
    net_face = positive_face - negative_face
    net_node = positive_node - negative_node
    standard_error_face = np.std(net_face, axis=0, ddof=1) / np.sqrt(count)
    standard_error_node = np.std(net_node, axis=0, ddof=1) / np.sqrt(count)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    audit_path = output / "ensemble_current_audit.npz"
    _atomic_npz(
        audit_path,
        ensemble_mean_positive_face_current_density_a_m2=np.mean(positive_face, axis=0),
        ensemble_mean_negative_face_current_density_a_m2=np.mean(negative_face, axis=0),
        ensemble_net_face_current_density_standard_error_a_m2=standard_error_face,
        ensemble_mean_positive_current_node_a=np.mean(positive_node, axis=0),
        ensemble_mean_negative_current_node_a=np.mean(negative_node, axis=0),
        ensemble_net_node_current_standard_error_a=standard_error_node,
        replicate_count=np.asarray(count),
        potential_before_v=potential[0],
        maximum_potential_spread_v=np.asarray(potential_spread),
        physical_face_area_m2=reference["physical_face_area_m2"],
        patch_scales_m=reference["patch_scales_m"],
        patch_group_by_scale=reference["patch_group_by_scale"])
    net_mean = np.mean(net_node, axis=0)
    signal_l2 = float(np.linalg.norm(net_mean))
    uncertainty_l2 = float(np.linalg.norm(standard_error_node))
    manifest = dict(
        schema="petch.charging.c3.fixed-state-ensemble-current-audit.v1",
        replicate_count=count,
        maximum_potential_spread_v=potential_spread,
        maximum_relative_potential_spread=potential_spread / potential_scale,
        net_node_current_signal_l2_a=signal_l2,
        net_node_current_standard_error_l2_a=uncertainty_l2,
        signal_to_standard_error_l2=(
            signal_l2 / uncertainty_l2 if uncertainty_l2 > 0.0 else float("inf")),
        provenance=[
            dict(source_directory=path.parent.name, name=path.name, sha256=_hash(path))
            for path in paths
        ],
        artifact=dict(name=audit_path.name, sha256=_hash(audit_path)))
    manifest_path = output / "summary.json"
    _atomic_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
