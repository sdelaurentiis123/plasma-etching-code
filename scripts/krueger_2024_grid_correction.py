#!/usr/bin/env python3
"""Derive the single permitted Krüger R1.9 fine-grid base correction.

The correction is earned only when the first 5 nm base endpoint misses either preregistered base
target by more than one fine cell.  It reuses the checksum-bound two-parameter response matrix from
the final 10 nm axisymmetric derivation and consumes no transfer/held-out observation.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md"
TARGET = np.asarray([45.0, 825.0], dtype=float)
PARAMETERS = (
    "effective_mask_crosslinked_growth_fraction",
    "oxide_etch_yield_scale",
)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _verify_embedded_sha256(payload, field):
    canonical = dict(payload)
    claimed = canonical.pop(field, None)
    actual = sha256(json.dumps(
        canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if claimed != actual:
        raise ValueError(f"{field} does not match its artifact content")


def _endpoint(path, expected_dx):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("configuration", {})
    if (payload.get("status") != "complete"
            or config.get("boundary_case") != "base"
            or not np.isclose(float(config.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(config.get("dx_um", np.nan)), expected_dx)
            or config.get("ion_azimuthal_closure") != "axisymmetric_uniform"
            or int(config.get("ion_azimuthal_order", -1)) != 16):
        raise ValueError(f"not a complete R1.9 base endpoint at dx={expected_dx}: {path}")
    pair = np.asarray([float(config[name]) for name in PARAMETERS], dtype=float)
    response = np.asarray([
        float(payload["final_metrics"]["mask_opening_nm"]),
        float(payload["final_metrics"]["etch_depth_nm"]),
    ], dtype=float)
    return {
        "path_name": path.name,
        "sha256": _sha(path),
        "config_hash": payload["config_hash"],
        "dx_um": float(expected_dx),
        "parameters": pair,
        "response": response,
    }


def derive(derivation_path, base_10nm_path, base_5nm_path):
    derivation_path = Path(derivation_path)
    derivation = json.loads(derivation_path.read_text(encoding="utf-8"))
    if (derivation.get("schema")
            != "petch.krueger-2024.base-axisymmetric-secant.v1"
            or derivation.get("held_out_profile_data_read") is not False
            or derivation.get("protocol_sha256") != _sha(PROTOCOL)):
        raise ValueError("grid correction requires the sealed final 10 nm base derivation")
    _verify_embedded_sha256(derivation, "proposal_sha256")
    base_10 = _endpoint(base_10nm_path, 0.01)
    base_5 = _endpoint(base_5nm_path, 0.005)
    proposed_10 = np.asarray([
        float(derivation["proposed_configuration"][name]) for name in PARAMETERS],
        dtype=float)
    if (not np.allclose(base_10["parameters"], proposed_10, rtol=0.0, atol=1e-14)
            or not np.allclose(base_5["parameters"], proposed_10, rtol=0.0, atol=1e-14)):
        raise ValueError("10/5 nm endpoints do not use the sealed calibration pair")

    error = base_5["response"] - TARGET
    if np.all(np.abs(error) <= 5.0):
        raise ValueError("5 nm endpoint already passes; no grid correction is earned")
    jacobian = np.asarray(derivation["derivation"]["jacobian"], dtype=float)
    if jacobian.shape != (2, 2) or not np.all(np.isfinite(jacobian)):
        raise ValueError("sealed calibration response matrix is invalid")
    condition_number = float(np.linalg.cond(jacobian))
    if condition_number > 1e6:
        raise ValueError("sealed calibration response matrix is ill-conditioned")
    update = np.linalg.solve(jacobian, -error)
    proposed = proposed_10 + update
    if (not np.all(np.isfinite(proposed)) or proposed[0] < 0.0 or proposed[0] > 1.0
            or proposed[1] <= 0.0):
        raise ValueError("fine-grid correction lies outside physical parameter bounds")

    payload = {
        "schema": "petch.krueger-2024.base-grid-correction.v1",
        "protocol_id": "K24-PETCH-R1.9",
        "protocol_sha256": _sha(PROTOCOL),
        "scientific_status": "single permitted 5 nm base-only correction; not validation",
        "calibration_derivation": {
            "path_name": derivation_path.name,
            "sha256": _sha(derivation_path),
            "proposal_sha256": derivation["proposal_sha256"],
        },
        "base_10nm_endpoint": {
            key: value.tolist() if isinstance(value, np.ndarray) else value
            for key, value in base_10.items()
        },
        "base_5nm_endpoint": {
            key: value.tolist() if isinstance(value, np.ndarray) else value
            for key, value in base_5.items()
        },
        "derivation": {
            "response_order": ["mask_opening_nm", "etch_depth_nm"],
            "parameter_order": list(PARAMETERS),
            "frozen_jacobian": jacobian.tolist(),
            "condition_number": condition_number,
            "fine_grid_target_error": error.tolist(),
            "parameter_update": update.tolist(),
        },
        "proposed_configuration": {
            PARAMETERS[0]: float(proposed[0]),
            PARAMETERS[1]: float(proposed[1]),
            "ion_azimuthal_closure": "axisymmetric_uniform",
            "ion_azimuthal_order": 16,
        },
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--calibration-derivation", required=True)
    parser.add_argument("--base-10nm", required=True)
    parser.add_argument("--base-5nm", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = derive(
        args.calibration_derivation, args.base_10nm, args.base_5nm)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps(payload["proposed_configuration"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
