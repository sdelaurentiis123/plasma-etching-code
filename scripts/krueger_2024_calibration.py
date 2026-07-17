#!/usr/bin/env python3
"""Produce a checksum-bound Krüger R1 calibration proposal from the two base endpoints.

This tool consumes only the preregistered base-case mask opening.  It never loads the transfer
observation table or any held-out profile.  The output is a proposal, not a reveal: the interpolated
fraction must still be confirmed at 10 nm and subjected to the R1.1 5 nm grid procedure.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "krueger_2024"
PROTOCOL = ROOT / "KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md"


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load_endpoint(path, expected_fraction):
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("configuration", {})
    geometry = config.get("geometry", {})
    fraction = float(config.get("effective_mask_crosslinked_growth_fraction", 0.0))
    if (not isinstance(geometry, dict)
            or config.get("boundary_case") != "base"
            or config.get("charging")
            != "disabled_for_Krueger_2024_calibration_and_transfer"
            or not np.isclose(float(config.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(config.get("dx_um", np.nan)), 0.01)
            or not np.isclose(float(geometry.get("substrate_top_um", np.nan)), 1.8)
            or not np.isclose(float(geometry.get("domain_height_um", np.nan)), 2.8)
            or not np.isclose(fraction, expected_fraction)):
        raise ValueError("calibration endpoint does not match protocol K24-PETCH-R1.1")
    status = payload.get("status")
    terminal = payload.get("terminal_event")
    if status == "complete":
        if not np.isclose(
                float(payload["history"][-1]["physical_time_s"]), 60.0,
                rtol=0.0, atol=1e-10):
            raise ValueError("completed endpoint did not reach 60 s")
        opening = float(payload["final_metrics"]["mask_opening_nm"])
        if not np.isfinite(opening) or opening < 0.0:
            raise ValueError("completed endpoint has invalid mask opening")
        endpoint_kind = "resolved_opening"
    elif (status == "terminal_feature_clogged" and isinstance(terminal, dict)
          and terminal.get("kind") == "feature_clogged"):
        # Once the gas path through the mask closes, its physical opening is exactly zero.  The
        # pre-event checkpoint remains untouched; the topology event supplies a censored endpoint.
        opening = 0.0
        endpoint_kind = "resolved_closure_event"
    else:
        raise ValueError("calibration endpoint is incomplete")
    return {
        "audit_path_name": path.name,
        "audit_sha256": _sha(path),
        "fraction": fraction,
        "mask_opening_nm": opening,
        "endpoint_kind": endpoint_kind,
        "status": status,
    }


def _target_opening():
    rows = (DATA / "base_case_metrics.csv").read_text(encoding="utf-8").splitlines()
    header = rows[0].split(",")
    records = [dict(zip(header, row.split(","))) for row in rows[1:] if row.strip()]
    selected = [row for row in records if row["metric"] == "mask_opening"]
    if (len(selected) != 1 or selected[0]["split"] != "calibration"
            or selected[0]["evidence_type"] != "experiment"):
        raise ValueError("missing sole Krüger mask-opening calibration target")
    return float(selected[0]["value"]), _sha(DATA / "base_case_metrics.csv")


def propose(endpoint_zero, endpoint_one):
    zero = _load_endpoint(endpoint_zero, 0.0)
    one = _load_endpoint(endpoint_one, 1.0)
    target, target_sha = _target_opening()
    y0 = zero["mask_opening_nm"]
    y1 = one["mask_opening_nm"]
    if (y0 - target) * (y1 - target) > 0.0 or y0 == y1:
        raise ValueError("crosslink endpoints do not bracket the calibration opening")
    fraction = float(np.clip((target - y0) / (y1 - y0), 0.0, 1.0))
    output = {
        "protocol_id": "K24-PETCH-R1.1",
        "protocol_sha256": _sha(PROTOCOL),
        "calibration_target": {
            "observable": "mask_opening",
            "value_nm": target,
            "source_table_sha256": target_sha,
        },
        "endpoints": [zero, one],
        "interpolation": {
            "method": "preregistered_monotone_secant",
            "proposed_effective_mask_crosslinked_growth_fraction": fraction,
            "opening_slope_nm_per_fraction": y1 - y0,
        },
        "operator_scope": "reduced_chemistry_with_charging_disabled",
        "status": (
            "uncharged_development_proposal_requires_charging_causality_audit_"
            "and_10nm_confirmation_and_5nm_grid_procedure"),
        "held_out_profile_data_read": False,
    }
    canonical = json.dumps(output, sort_keys=True, separators=(",", ":"))
    output["proposal_sha256"] = sha256(canonical.encode("utf-8")).hexdigest()
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-zero", required=True)
    parser.add_argument("--endpoint-one", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = propose(args.endpoint_zero, args.endpoint_one)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result["interpolation"], sort_keys=True))


if __name__ == "__main__":
    main()
