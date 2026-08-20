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
import json

import numpy as np

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
    parser.add_argument(
        "--trace-symmetry-shapes", action="store_true",
        help="print sparse-event and expanded-event counts before strip projection",
    )
    parser.add_argument(
        "--trace-advection", action="store_true",
        help="print the CFL substep count and extended-speed extrema per feature step",
    )
    parser.add_argument(
        "--explore-post-mask", action="store_true",
        help=(
            "diagnostic only: bypass the campaign's sub-cell Cr stop and do "
            "not write a production cache"
        ),
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

    if args.trace_symmetry_shapes:
        # Patch only the diagnostic process.  The wrapper reports the exact
        # cardinality presented to the unchanged production implementation.
        # `symmetrize_transport_across_strips` resolves this private helper
        # through its defining module's globals, so this also observes the
        # function imported by feature_step_3d.
        from petch import boundary_transport_3d as boundary_transport

        real_symmetrize = boundary_transport._symmetrize_face_resolved

        def traced_symmetrize(population, areas, groups):
            event_face = np.asarray(population.event_face, dtype=int)
            group_size_by_face = np.ones(population.face_count, dtype=np.int64)
            for members in groups.values():
                group_size_by_face[np.asarray(members, dtype=int)] = len(members)
            expanded = int(np.sum(
                group_size_by_face[event_face], dtype=np.int64
            ))
            largest = max((len(members) for members in groups.values()), default=0)
            print(
                "SYMMETRY_SHAPE "
                f"population={population.name} faces={population.face_count} "
                f"events={event_face.size} groups={len(groups)} "
                f"largest_group={largest} expanded_events={expanded}",
                flush=True,
            )
            return real_symmetrize(population, areas, groups)

        boundary_transport._symmetrize_face_resolved = traced_symmetrize

    if args.trace_advection:
        from petch import feature_step_3d as feature_step

        real_advect_materials = feature_step._advect_exposed_material_levelsets

        def traced_advect_materials(
                material_levelsets, etchable_material_ids, extended_velocity,
                dx, duration_s, substeps, *, periodic_lateral=False):
            velocity = np.asarray(extended_velocity, dtype=float)
            print(
                "ADVECTION_SHAPE "
                f"grid={velocity.shape} duration_s={float(duration_s):.17g} "
                f"dx={float(dx):.17g} substeps={int(substeps)} "
                f"speed_abs_max={float(np.max(np.abs(velocity))):.17g} "
                f"speed_min={float(np.min(velocity)):.17g} "
                f"speed_max={float(np.max(velocity)):.17g}",
                flush=True,
            )
            return real_advect_materials(
                material_levelsets, etchable_material_ids, extended_velocity,
                dx, duration_s, substeps,
                periodic_lateral=periodic_lateral,
            )

        feature_step._advect_exposed_material_levelsets = traced_advect_materials

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
        if args.explore_post_mask:
            # The v3 campaign guard was added after a grid-aligned sliver made
            # mesh extraction fail.  Mesh extraction and strip projection have
            # since both been repaired independently.  Bypass only these two
            # driver-level stop predicates to discover what the unchanged core
            # now does after local Cr perforation.  Never write this exploratory
            # path under a production content-addressed cache key.
            real_minimum = board._mask_interior_minimum_thickness_nm
            real_metrics = board._mask_metrics

            def bypassed_minimum(*_args, **_kwargs):
                return float("inf")

            def diagnostic_metrics(*metric_args, **metric_kwargs):
                metrics = dict(real_metrics(*metric_args, **metric_kwargs))
                metrics["mask_below_vertical_resolution_at_center"] = False
                metrics["mask_exhausted_at_center"] = False
                return metrics

            board._mask_interior_minimum_thickness_nm = bypassed_minimum
            board._mask_metrics = diagnostic_metrics
            preregistration = board._load(board.PREREGISTRATION)
            width, scenario, rates, selectivity, duration, dx = job
            profiles = board._run_trajectory(
                width_nm=width,
                scenario=scenario,
                rates_nm_min=rates,
                selectivity=selectivity,
                duration_s=duration,
                dx_nm=dx,
                preregistration=preregistration,
                transport_device=args.transport_device,
            )
            path = None
            device = args.transport_device
        else:
            profiles, path, device = board._execute(
                (job, args.transport_device)
            )
    finally:
        faulthandler.cancel_dump_traceback_later()
    if args.explore_post_mask:
        print(
            f"COMPLETE_EXPLORATION wall_s={monotonic() - started:.3f} "
            f"device={device} profiles={len(profiles)} cache_written=false",
            flush=True,
        )
        print(json.dumps(profiles, indent=2, sort_keys=True), flush=True)
    else:
        print(
            f"COMPLETE wall_s={monotonic() - started:.3f} device={device} "
            f"profiles={len(profiles)} cache={path}",
            flush=True,
        )


if __name__ == "__main__":
    main()
