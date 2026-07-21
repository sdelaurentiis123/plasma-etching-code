#!/usr/bin/env python3
"""Emit the R5 authority launch manifest binding a candidate run to its source epoch.

Reads the candidate's completed audit, the shipped source archive, and the current git
state; refuses a dirty tree.  The manifest carries the executable checksums the R5 freeze
verifies against current source, the declared deterministic-exchange operator inputs from
the audit configuration, and the standard sealed-data attestations.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent.parent
EXECUTABLES = (
    "scripts/krueger_2024_trench_pilot.py",
    "src/petch/boundary_transport_3d.py",
    "src/petch/feature_step_3d.py",
    "src/petch/deterministic_exchange_2d.py",
    "src/petch/extruded_exchange_3d.py",
    "src/petch/surface_common_refinement_3d.py",
)


def _sha_file(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", required=True)
    parser.add_argument("--archive", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    audit = json.loads(Path(args.audit).read_text(encoding="utf-8"))
    if audit.get("status") != "complete":
        raise ValueError("launch manifest requires a complete candidate audit")
    config = audit["configuration"]
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
        check=True).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True,
        check=True).stdout.strip())
    if dirty:
        raise ValueError("refusing to bind a launch manifest to a dirty tree")

    payload = {
        "schema": "petch.krueger-2024.r5-base-authority.v1",
        "protocol_id": "K24-PETCH-R5",
        "authority_candidate": True,
        "held_out_profile_data_read": False,
        "calibration_performed_by_this_launch": False,
        "source_epoch": {
            "git_revision": revision,
            "git_dirty": False,
            "archive_sha256": _sha_file(args.archive),
            "archive_name": Path(args.archive).name,
        },
        "fixed_parameters": {
            name: float(config[name])
            for name in ("effective_mask_crosslinked_growth_fraction",
                         "oxide_etch_yield_scale")},
        "numerical_operator": {
            "duration_s": float(config["duration_s"]),
            "dx_um": float(config["dx_um"]),
            "radiosity_backend": str(config["radiosity_backend"]),
            "deterministic_exchange": {
                "exchange_method": str(
                    config["deterministic_exchange"]["exchange_method"]),
                "exchange_relative_tolerance": float(
                    config["deterministic_exchange"]["exchange_relative_tolerance"]),
                "exchange_geometry_tolerance": float(
                    config["deterministic_exchange"]["exchange_geometry_tolerance"]),
                "maximum_refinement_level": int(
                    config["deterministic_exchange"]["maximum_refinement_level"]),
                "extrusion_projection_guard_cells": float(
                    config["deterministic_exchange"][
                        "extrusion_projection_guard_cells"]),
            },
            "surface_state_remap_backend": str(config["surface_state_remap_backend"]),
        },
        "candidate_audit": {
            "path_name": Path(args.audit).name,
            "sha256": _sha_file(args.audit),
            "config_hash": str(audit["config_hash"]),
        },
        "executable_source_sha256": {
            name: _sha_file(ROOT / name) for name in EXECUTABLES},
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["launch_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(output)
    print(payload["launch_sha256"])


if __name__ == "__main__":
    main()
