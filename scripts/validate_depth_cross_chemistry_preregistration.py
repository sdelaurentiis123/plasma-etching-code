#!/usr/bin/env python3
"""Validate the value-blind, cross-chemistry absolute-depth commitment."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

from petch.validation_contract import (
    PreRegisteredValidationProtocol,
    ValidationParameter,
    ValidationTargetCommitment,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATH = (
    ROOT / "data" / "experimental" / "depth_cross_chemistry_v1"
    / "preregistration.json"
)


def _expand_target_set(specification):
    common = dict(specification["common"])
    kind = specification["kind"]
    if kind == "rows":
        for row in specification["rows"]:
            yield ValidationTargetCommitment(**(common | dict(row)))
        return
    if kind != "cartesian":
        raise ValueError(f"unknown target-set kind: {kind}")
    axes = specification["axes"]
    names = tuple(axes)
    for values in product(*(axes[name] for name in names)):
        variables = dict(zip(names, values))
        fields = {
            name: value.format(**variables) if isinstance(value, str) else value
            for name, value in common.items()
        }
        yield ValidationTargetCommitment(**fields)


def load_preregistration(path=DEFAULT_PATH):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    expected_keys = {
        "schema", "protocol", "parameters", "target_sets",
        "expected_commit_sha256",
    }
    if set(payload) != expected_keys or payload["schema"] != (
        "petch.value-blind-validation-preregistration.v1"
    ):
        raise ValueError("unexpected preregistration schema")
    targets = tuple(
        target
        for target_set in payload["target_sets"]
        for target in _expand_target_set(target_set)
    )
    protocol = PreRegisteredValidationProtocol(
        targets=targets,
        parameters=tuple(
            ValidationParameter(**item) for item in payload["parameters"]),
        **payload["protocol"],
    )
    return protocol, payload["expected_commit_sha256"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=Path, default=DEFAULT_PATH)
    parser.add_argument(
        "--show-commit", action="store_true",
        help="print the computed commitment without requiring it to match",
    )
    args = parser.parse_args()
    protocol, expected = load_preregistration(args.path)
    if not args.show_commit and protocol.commit_sha256 != expected:
        raise SystemExit(
            "preregistration commitment mismatch: "
            f"expected {expected}, computed {protocol.commit_sha256}"
        )
    families = sorted({
        item.chemistry_family
        for item in protocol.targets
        if item.split == "held_out_transfer"
    })
    print(json.dumps({
        "protocol_id": protocol.protocol_id,
        "commit_sha256": protocol.commit_sha256,
        "target_count": len(protocol.targets),
        "calibration_count": sum(
            item.split == "calibration" for item in protocol.targets),
        "boundary_input_count": sum(
            item.split == "boundary_input" for item in protocol.targets),
        "held_out_count": sum(
            item.split == "held_out_transfer" for item in protocol.targets),
        "held_out_chemistry_families": families,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
