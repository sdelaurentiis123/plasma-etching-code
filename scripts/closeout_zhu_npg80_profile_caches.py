#!/usr/bin/env python3
"""Inventory completed Oxford profile trajectories without computing missing jobs.

The production campaigns are expensive and can be stopped intentionally.  This
script turns the cache state into a deterministic, checksum-bound receipt.  It
never calls either board's ``build`` function, so a closeout cannot silently
launch a missing trajectory or promote an incomplete grid to a complete board.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_zhu_npg80_conditional_profiles import (
    _load,
    _scenario_inputs,
)
from scripts import audit_zhu_npg80_gds_square_profiles as exact
from scripts import audit_zhu_npg80_moving_cr_profiles as legacy


LEGACY_OUTPUT = legacy.OUTPUT.with_name("closeout_partial.json")
EXACT_OUTPUT = exact.OUTPUT.with_name("closeout_partial.json")


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _legacy_jobs():
    preregistration = _load(legacy.PREREGISTRATION)
    analog = _load(legacy.ANALOG_BOARD)
    reactor = _load(legacy.REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    )
    selectivities = tuple(float(value) for value in
                          preregistration["surface_response_axes"][
                              "tio2_to_cr_selectivity"])
    widths = tuple(float(value) for value in
                   preregistration["inferred_geometry_board"]["width_nm"])
    return [
        (width, scenario, rates, selectivity, 1200.0,
         legacy.PRODUCTION_MESH_SPACING_NM)
        for width in widths
        for scenario in scenarios
        for selectivity in selectivities
    ]


def _range(profiles, section, field):
    values = [float(profile[section][field]) for profile in profiles]
    return [float(min(values)), float(max(values))]


def summarize_profiles(width, profiles):
    """Return an envelope only when every expected trajectory is present."""
    return {
        "width_nm": float(width),
        "profile_count": len(profiles),
        "etch_depth_nm": _range(profiles, "profile", "etched_depth_nm"),
        "top_cd_nm": _range(profiles, "profile", "top_cd_nm"),
        "middle_cd_nm": _range(profiles, "profile", "middle_cd_nm"),
        "bottom_cd_nm": _range(profiles, "profile", "bottom_cd_nm"),
        "sidewall_angle_from_wafer_deg": _range(
            profiles, "profile", "sidewall_angle_from_wafer_deg"),
        "bow_nm": _range(profiles, "profile", "bow_nm"),
        "cr_center_remaining_thickness_nm": _range(
            profiles, "cr_mask", "center_remaining_thickness_nm"),
        "cr_center_exhausted_fraction": float(np.mean([
            bool(profile["cr_mask"]["mask_exhausted_at_center"])
            for profile in profiles
        ])),
        "terminal_reasons": sorted(set(
            str(profile["terminal_reason"]) for profile in profiles)),
        "maximum_transport_relative_particle_balance_error": max(
            float(profile[
                "maximum_transport_relative_particle_balance_error"])
            for profile in profiles),
        "maximum_state_remap_relative_conservation_residual": max(
            float(profile[
                "maximum_state_remap_relative_conservation_residual"])
            for profile in profiles),
    }


def inventory(*, name, jobs, job_spec, cache_path, preregistration,
              gds_path=None):
    valid = []
    missing = []
    corrupt = []
    grouped_profiles = {}
    trajectory_counts = {}
    for job in jobs:
        spec = job_spec(job)
        path = cache_path(spec)
        width = float(spec["width_nm"])
        trajectory_counts.setdefault(width, 0)
        if not path.exists():
            missing.append(spec)
            continue
        try:
            payload = _load(path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            corrupt.append({"path": str(path.relative_to(ROOT)),
                            "error": str(error)})
            continue
        if payload.get("job_spec") != spec:
            corrupt.append({"path": str(path.relative_to(ROOT)),
                            "error": "job_spec mismatch"})
            continue
        profiles = payload.get("profiles")
        if not isinstance(profiles, list) or not profiles:
            corrupt.append({"path": str(path.relative_to(ROOT)),
                            "error": "profiles must be a nonempty list"})
            continue
        receipt = {
            "path": str(path.relative_to(ROOT)),
            "sha256": _hash(path),
            "profile_count": len(profiles),
            "transport_device": payload.get("execution", {}).get(
                "transport_device", "unrecorded"),
            "job_spec": spec,
        }
        valid.append(receipt)
        trajectory_counts[width] += 1
        grouped_profiles.setdefault(width, []).extend(profiles)

    expected_per_width = len(jobs) // len(trajectory_counts)
    complete_widths = sorted(
        width for width, count in trajectory_counts.items()
        if count == expected_per_width)
    envelopes = [
        summarize_profiles(width, grouped_profiles[width])
        for width in complete_widths
    ]
    output = {
        "schema": "petch.zhu-npg80-profile-cache-closeout.v1",
        "campaign": name,
        "result_status": (
            "complete" if len(valid) == len(jobs) and not corrupt
            else "partial_grid_do_not_interpolate"
        ),
        "board_complete": len(valid) == len(jobs) and not corrupt,
        "target_sem_used": False,
        "target_depth_used": False,
        "coefficient_selected_from_target": None,
        "absolute_oxford_profile_prediction_certified": False,
        "preregistration": {
            "path": str(Path(preregistration).relative_to(ROOT)),
            "sha256": _hash(preregistration),
        },
        "physics_grid": {
            "expected_trajectory_count": len(jobs),
            "valid_trajectory_count": len(valid),
            "missing_trajectory_count": len(missing),
            "corrupt_trajectory_count": len(corrupt),
            "expected_trajectories_per_width": expected_per_width,
            "completed_width_nm": complete_widths,
            "trajectory_count_by_width_nm": {
                str(width): count
                for width, count in sorted(trajectory_counts.items())
            },
        },
        "claim_boundary": {
            "envelopes_include_only_complete_width_slices": True,
            "missing_grid_cells_are_not_inferred": True,
            "surface_law_is_cross_machine_rate_normalized": True,
            "mask_polarity_remains_operator_unconfirmed": name == "exact_gds",
            "target_sem_is_required_for_blind_validation_not_as_input": True,
        },
        "complete_width_envelopes": envelopes,
        "valid_trajectory_receipts": valid,
        "missing_job_specs": missing,
        "corrupt_cache_receipts": corrupt,
    }
    if gds_path is not None:
        output["gds"] = {
            "path": str(Path(gds_path).relative_to(ROOT)),
            "sha256": _hash(gds_path),
        }
    return output


def build():
    exact_jobs = exact._jobs(smoke=False)[0]
    exact_document = _load(exact.PREREGISTRATION)
    gds_path = exact.DATA / exact_document["geometry_source"]["file"]
    return {
        LEGACY_OUTPUT: inventory(
            name="legacy_square_width_prior",
            jobs=_legacy_jobs(),
            job_spec=legacy._job_spec,
            cache_path=legacy._cache_path,
            preregistration=legacy.PREREGISTRATION,
        ),
        EXACT_OUTPUT: inventory(
            name="exact_gds",
            jobs=exact_jobs,
            job_spec=exact._job_spec,
            cache_path=exact._cache_path,
            preregistration=exact.PREREGISTRATION,
            gds_path=gds_path,
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    outputs = build()
    for path, payload in outputs.items():
        rendered = _render(payload)
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(rendered, encoding="utf-8")
            print(path.relative_to(ROOT))
        elif not path.exists() or path.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"profile cache closeout is stale: {path}")
        else:
            print(f"PASS {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
