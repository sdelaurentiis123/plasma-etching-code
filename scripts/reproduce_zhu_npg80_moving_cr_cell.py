#!/usr/bin/env python3
"""Reproduce one frozen Oxford moving-Cr blind-board trajectory.

The production audit intentionally executes a content-addressed board in
parallel.  This narrow driver makes a failed cell independently reproducible
without changing the job specification, consuming other missing cells, or
writing a board cache.  Its JSON output is diagnostic only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_zhu_npg80_moving_cr_profiles import (
    ANALOG_BOARD,
    PREREGISTRATION,
    REACTOR_DOSE,
    _load,
    _run_trajectory,
    _scenario_inputs,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width-nm", type=float, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--selectivity", type=float, required=True)
    parser.add_argument("--duration-s", type=float, default=1200.0)
    parser.add_argument("--dx-nm", type=float, default=10.0)
    parser.add_argument("--transport-device", default="cpu")
    args = parser.parse_args()

    preregistration = _load(PREREGISTRATION)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    matching = tuple(
        scenario for scenario in scenarios
        if str(scenario["name"]) == str(args.scenario)
    )
    if len(matching) != 1:
        available = tuple(str(scenario["name"]) for scenario in scenarios)
        raise ValueError(
            f"scenario {args.scenario!r} is not unique; available={available}"
        )
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    )
    profiles = _run_trajectory(
        width_nm=args.width_nm,
        scenario=matching[0],
        rates_nm_min=rates,
        selectivity=args.selectivity,
        duration_s=args.duration_s,
        dx_nm=args.dx_nm,
        preregistration=preregistration,
        transport_device=args.transport_device,
    )
    print(json.dumps(profiles, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
