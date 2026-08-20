#!/usr/bin/env python3
"""Run one moving-Cr trajectory with an externally triggerable stack trace.

This driver exists for deterministic production cells which stop making
progress inside the feature step.  It deliberately calls the same `_execute`
path as the board, so a successful run creates the normal content-addressed
cache.  It does not change any physics, retry budget, or certification rule.

While the process is running, ``kill -USR1 <pid>`` writes every Python thread
stack to stderr.  A periodic watchdog does the same automatically.
"""
from __future__ import annotations

import argparse
import faulthandler
import os
import signal
import sys
from time import monotonic

from scripts import audit_zhu_npg80_moving_cr_profiles as board


def _production_job(*, width_nm: float, scenario_name: str, selectivity: float):
    preregistration = board._load(board.PREREGISTRATION)
    analog = board._load(board.ANALOG_BOARD)
    reactor = board._load(board.REACTOR_DOSE)
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"
        ]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"
        ]),
    )
    scenarios = board._scenario_inputs(preregistration, reactor)
    matching = [item for item in scenarios if item["name"] == scenario_name]
    if len(matching) != 1:
        names = ", ".join(item["name"] for item in scenarios)
        raise ValueError(f"scenario {scenario_name!r} not in: {names}")
    allowed_widths = tuple(float(value) for value in
                           preregistration["inferred_geometry_board"][
                               "width_nm"
                           ])
    allowed_selectivities = tuple(float(value) for value in
                                  preregistration["surface_response_axes"][
                                      "tio2_to_cr_selectivity"
                                  ])
    if float(width_nm) not in allowed_widths:
        raise ValueError(f"width {width_nm} not in board {allowed_widths}")
    if float(selectivity) not in allowed_selectivities:
        raise ValueError(
            f"selectivity {selectivity} not in board {allowed_selectivities}"
        )
    return (
        float(width_nm), matching[0], rates, float(selectivity), 1200.0,
        board.PRODUCTION_MESH_SPACING_NM,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width-nm", type=float, default=320.0)
    parser.add_argument("--scenario", default="ion_low_tail_0p0")
    parser.add_argument("--selectivity", type=float, default=14.0)
    parser.add_argument("--transport-device", default="cuda:0")
    parser.add_argument(
        "--stack-interval-s", type=float, default=180.0,
        help="seconds between automatic all-thread Python stack dumps",
    )
    args = parser.parse_args()
    if args.stack_interval_s <= 0.0:
        parser.error("--stack-interval-s must be positive")

    faulthandler.enable(file=sys.stderr, all_threads=True)
    if hasattr(signal, "SIGUSR1"):
        faulthandler.register(
            signal.SIGUSR1, file=sys.stderr, all_threads=True, chain=False
        )
    faulthandler.dump_traceback_later(
        args.stack_interval_s, repeat=True, file=sys.stderr, exit=False
    )

    job = _production_job(
        width_nm=args.width_nm,
        scenario_name=args.scenario,
        selectivity=args.selectivity,
    )
    spec = board._job_spec(job)
    cache = board._cache_path(spec)
    print(f"pid={os.getpid()}", flush=True)
    print(f"job_spec={board._render(spec).strip()}", flush=True)
    print(f"cache={cache}", flush=True)
    started = monotonic()
    try:
        profiles, path, device = board._execute((job, args.transport_device))
    finally:
        faulthandler.cancel_dump_traceback_later()
    print(
        f"COMPLETE wall_s={monotonic() - started:.3f} device={device} "
        f"profiles={len(profiles)} cache={path}",
        flush=True,
    )


if __name__ == "__main__":
    main()
