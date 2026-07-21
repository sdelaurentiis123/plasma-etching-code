#!/usr/bin/env python3
"""Deterministic Gaussian-process feasibility proposer for the R3 base calibration.

Consumes completed 10 nm base endpoint audits (base runs only; no held-out access), fits an
independent GP per declared observable over (mask fraction, oxide yield scale), and proposes
the batch of candidate pairs maximizing the joint probability that BOTH smoothed observables
land within the declared +/-5 nm base tolerance.  Fixed seed; observation noise set from the
measured late-time neck jitter.  Receipted: the proposal embeds every input endpoint hash.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern

BOUNDS = {"fraction": (0.84, 0.93), "yield_scale": (0.52, 0.58)}
TARGETS = {"mask_opening_nm": 45.0, "etch_depth_nm": 825.0}
TOLERANCE_NM = 5.0
NOISE_NM = {"mask_opening_nm": 1.5, "etch_depth_nm": 3.0}
SMOOTH_LAST_STEPS = 5


def smoothed_endpoint(audit_path):
    audit = json.loads(Path(audit_path).read_text())
    if audit.get("status") != "complete":
        raise ValueError(f"not a complete endpoint: {audit_path}")
    config = audit["configuration"]
    history = audit["history"][-SMOOTH_LAST_STEPS:]
    return {
        "fraction": float(config["effective_mask_crosslinked_growth_fraction"]),
        "yield_scale": float(config["oxide_etch_yield_scale"]),
        "mask_opening_nm": float(np.median(
            [h["metrics"]["mask_opening_nm"] for h in history])),
        "etch_depth_nm": float(np.median(
            [h["metrics"]["etch_depth_nm"] for h in history])),
        "sha256": sha256(Path(audit_path).read_bytes()).hexdigest(),
        "path_name": Path(audit_path).name,
    }


def propose(endpoints, batch):
    x = np.array([[e["fraction"], e["yield_scale"]] for e in endpoints])
    span = np.array([BOUNDS["fraction"], BOUNDS["yield_scale"]])
    xn = (x - span[:, 0]) / (span[:, 1] - span[:, 0])
    models = {}
    for name in TARGETS:
        y = np.array([e[name] for e in endpoints])
        kernel = ConstantKernel(np.var(y) + 1.0) * Matern(
            length_scale=[0.5, 0.5], length_scale_bounds=(0.05, 5.0), nu=2.5)
        gp = GaussianProcessRegressor(
            kernel=kernel, alpha=NOISE_NM[name] ** 2, normalize_y=True,
            random_state=241, n_restarts_optimizer=8)
        gp.fit(xn, y)
        models[name] = gp

    grid_axis = np.linspace(0.0, 1.0, 61)
    gx, gy = np.meshgrid(grid_axis, grid_axis, indexing="ij")
    candidates = np.column_stack([gx.ravel(), gy.ravel()])
    feasibility = np.ones(len(candidates))
    for name, gp in models.items():
        mean, std = gp.predict(candidates, return_std=True)
        std = np.maximum(std, 1e-6)
        target = TARGETS[name]
        feasibility *= (norm.cdf((target + TOLERANCE_NM - mean) / std)
                        - norm.cdf((target - TOLERANCE_NM - mean) / std))

    order = np.argsort(-feasibility)
    chosen = []
    for index in order:
        point = candidates[index]
        if any(np.linalg.norm(point - c) < 0.12 for c in chosen):
            continue
        if any(np.linalg.norm(point - p) < 0.02 for p in xn):
            continue
        chosen.append(point)
        if len(chosen) == batch:
            break
    proposals = []
    for point in chosen:
        raw = span[:, 0] + point * (span[:, 1] - span[:, 0])
        proposals.append({
            "effective_mask_crosslinked_growth_fraction": float(raw[0]),
            "oxide_etch_yield_scale": float(raw[1]),
            "joint_feasibility": float(
                feasibility[np.argmin(np.linalg.norm(candidates - point, axis=1))]),
        })
    return proposals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint", action="append", required=True)
    parser.add_argument("--batch", type=int, default=2)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    endpoints = [smoothed_endpoint(path) for path in args.endpoint]
    proposals = propose(endpoints, args.batch)
    payload = {
        "schema": "petch.krueger-2024.r3-gp-proposal.v1",
        "protocol_id": "K24-PETCH-R3",
        "held_out_profile_data_read": False,
        "bounds": BOUNDS, "targets": TARGETS, "tolerance_nm": TOLERANCE_NM,
        "noise_nm": NOISE_NM, "smoothing_last_steps": SMOOTH_LAST_STEPS,
        "seed": 241,
        "inputs": [{k: e[k] for k in ("path_name", "sha256", "fraction",
                                      "yield_scale", "mask_opening_nm",
                                      "etch_depth_nm")} for e in endpoints],
        "proposals": proposals,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode()).hexdigest()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["proposals"], indent=2))


if __name__ == "__main__":
    main()
