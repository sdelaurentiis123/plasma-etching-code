#!/usr/bin/env python3
"""Deterministic Jacobian calibration proposer for the analytic-occlusion base runs.

The deterministic engine returns bit-reproducible endpoints, so base calibration no longer
needs noise-aware regression: two probe runs (one per knob, base-relative) measure the
2x2 endpoint Jacobian exactly, and a damped Newton step proposes the next candidate pair.
Consumes completed base endpoint audits only (no held-out access); every proposal embeds
the sha256 of every input audit.  Fail-closed: refuses mixed epochs, incomplete runs, or a
singular/ill-conditioned Jacobian.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

BOUNDS = {"fraction": (0.84, 0.93), "yield_scale": (0.50, 0.60)}
TARGETS = {"mask_opening_nm": 45.0, "etch_depth_nm": 825.0}
MAXIMUM_CONDITION = 1.0e4
MAXIMUM_STEP = {"fraction": 0.02, "yield_scale": 0.02}


def endpoint(audit_path):
    audit = json.loads(Path(audit_path).read_text())
    if audit.get("status") != "complete":
        raise ValueError(f"not a complete endpoint: {audit_path}")
    config = audit["configuration"]
    final = audit["history"][-1]["metrics"]
    return {
        "fraction": float(config["effective_mask_crosslinked_growth_fraction"]),
        "yield_scale": float(config["oxide_etch_yield_scale"]),
        "mask_opening_nm": float(final["mask_opening_nm"]),
        "etch_depth_nm": float(final["etch_depth_nm"]),
        "git_commit": str(audit.get("git", {}).get("commit", "")),
        "config_hash": str(audit.get("config_hash", "")),
        "sha256": sha256(Path(audit_path).read_bytes()).hexdigest(),
        "path_name": Path(audit_path).name,
    }


def propose(base, probe_fraction, probe_yield, damping):
    for name, probe in (("fraction", probe_fraction), ("yield_scale", probe_yield)):
        if probe["git_commit"] != base["git_commit"]:
            raise ValueError(f"{name} probe and base come from different engine epochs")
    df = probe_fraction["fraction"] - base["fraction"]
    dy = probe_yield["yield_scale"] - base["yield_scale"]
    if df == 0.0 or dy == 0.0:
        raise ValueError("probe runs must perturb their declared knob")
    if (probe_fraction["yield_scale"] != base["yield_scale"]
            or probe_yield["fraction"] != base["fraction"]):
        raise ValueError("each probe must hold the other knob at the base value")
    observables = ("mask_opening_nm", "etch_depth_nm")
    jacobian = np.empty((2, 2))
    for row, name in enumerate(observables):
        jacobian[row, 0] = (probe_fraction[name] - base[name]) / df
        jacobian[row, 1] = (probe_yield[name] - base[name]) / dy
    condition = float(np.linalg.cond(jacobian))
    if not np.isfinite(condition) or condition > MAXIMUM_CONDITION:
        raise ValueError(
            f"endpoint Jacobian is ill-conditioned (cond={condition:.3g}); "
            "calibration requires larger or better-separated probes")
    residual = np.array([
        TARGETS["mask_opening_nm"] - base["mask_opening_nm"],
        TARGETS["etch_depth_nm"] - base["etch_depth_nm"]])
    step = float(damping) * np.linalg.solve(jacobian, residual)
    step[0] = float(np.clip(step[0], -MAXIMUM_STEP["fraction"], MAXIMUM_STEP["fraction"]))
    step[1] = float(np.clip(step[1], -MAXIMUM_STEP["yield_scale"],
                            MAXIMUM_STEP["yield_scale"]))
    candidate = {
        "fraction": float(np.clip(
            base["fraction"] + step[0], *BOUNDS["fraction"])),
        "yield_scale": float(np.clip(
            base["yield_scale"] + step[1], *BOUNDS["yield_scale"])),
    }
    return {
        "candidate": candidate,
        "jacobian": jacobian.tolist(),
        "jacobian_condition": condition,
        "residual_nm": residual.tolist(),
        "newton_step": step.tolist(),
        "damping": float(damping),
        "targets": TARGETS,
        "engine_commit": base["git_commit"],
        "inputs": [
            {key: item[key] for key in (
                "path_name", "sha256", "config_hash", "fraction", "yield_scale",
                "mask_opening_nm", "etch_depth_nm")}
            for item in (base, probe_fraction, probe_yield)],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-audit", required=True)
    parser.add_argument("--probe-fraction-audit", required=True)
    parser.add_argument("--probe-yield-audit", required=True)
    parser.add_argument("--damping", type=float, default=1.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    proposal = propose(
        endpoint(args.base_audit), endpoint(args.probe_fraction_audit),
        endpoint(args.probe_yield_audit), args.damping)
    payload = json.dumps(proposal, indent=2, sort_keys=True)
    Path(args.output).write_text(payload + "\n", encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()
