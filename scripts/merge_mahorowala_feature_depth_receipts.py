#!/usr/bin/env python3
"""Strictly merge independently sharded Mahorowala feature receipts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


MATCH_KEYS = (
    "schema",
    "reactor_receipt",
    "dx_um",
    "product_wall_limit",
    "chlorine_specular_reflection",
    "observed_depth_used_for_conditioning",
    "geometry",
    "formal_feature_depth_pass",
    "evidence_blockers",
)


def merge(paths: list[Path]):
    if not paths:
        raise ValueError("at least one feature receipt is required")
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    reference = payloads[0]
    for payload in payloads[1:]:
        for key in MATCH_KEYS:
            if payload[key] != reference[key]:
                raise ValueError(f"feature receipt mismatch at {key}")
    rows = [row for payload in payloads for row in payload["rows"]]
    run_ids = [int(row["run"]) for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("duplicate feature run across shards")
    rows.sort(key=lambda row: int(row["run"]))
    usable = [
        float(row["signed_error_percent"])
        for row in rows if row["signed_error_percent"] is not None
    ]
    merged = {key: reference[key] for key in MATCH_KEYS}
    merged["schema"] = "petch.mahorowala_1998_deterministic_feature_depth.v1"
    merged["rows"] = rows
    merged["mape_percent"] = (
        None if not usable else float(np.mean(np.abs(usable))))
    merged["shard_receipts"] = [str(path) for path in paths]
    merged["shard_merge_recomputed_mape"] = True
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("receipts", nargs="+", type=Path)
    arguments = parser.parse_args()
    result = merge(arguments.receipts)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "run_count": len(result["rows"]),
        "mape_percent": result["mape_percent"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
