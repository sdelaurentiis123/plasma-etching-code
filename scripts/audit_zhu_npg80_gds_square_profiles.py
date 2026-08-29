#!/usr/bin/env python3
"""Run the moving-Cr Oxford profile model on Freddie's exact GDS squares.

The GDS fixes pitch and square critical dimensions.  It does not fix mask
polarity, absorbed RF power, the ion-energy distribution, or the TiO2/Cr
surface law.  This board therefore preserves the previously frozen physical
sensitivity axes while replacing the public-geometry prior with the supplied
layout.  No target SEM or target depth is consumed.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from hashlib import sha256
import multiprocessing
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_zhu_npg80_conditional_profiles import (
    _hash,
    _load,
    _scenario_inputs,
)
from scripts.audit_zhu_npg80_moving_cr_profiles import (
    _process_pool_options,
    _render,
    _run_trajectory,
)


DATA = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
PREREGISTRATION = DATA / "gds_square_profile_preregistration.json"
ANALOG_BOARD = DATA / "janissen_tio2_analog_board.json"
REACTOR_DOSE = (
    ROOT / "results" / "curated" / "zhu_npg80_daughter_wafer_dose_v1"
    / "audit.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_gds_square_profiles_v1"
    / "audit.json"
)
CACHE_DIR = OUTPUT.parent / "trajectories"
MODEL_REVISION = "two-material-moving-tio2-cr-exact-gds-squares-v1"
PRODUCTION_MESH_SPACING_NM = 10.0


def _physics_preregistration(document):
    """Adapt the exact-layout schema to the frozen moving-mask kernel API."""
    exact = document["exact_layout_geometry"]
    adapted = deepcopy(document)
    adapted["inferred_geometry_board"] = {
        "pitch_nm": float(exact["pitch_nm"]),
        "width_nm": [float(value) for value in exact["square_width_nm"]],
        "evidence_class": "operator_supplied_exact_GDSII",
        "target_layout_confirmed": True,
    }
    return adapted


def _job_spec(job):
    width, scenario, rates, selectivity, duration, dx = job
    return {
        "model_revision": MODEL_REVISION,
        "preregistration_sha256": _hash(PREREGISTRATION),
        "gds_sha256": _load(PREREGISTRATION)["geometry_source"]["sha256"],
        "width_nm": float(width),
        "scenario": dict(scenario),
        "rates_nm_min": [float(value) for value in rates],
        "selectivity": float(selectivity),
        "duration_s": float(duration),
        "mesh_spacing_nm": float(dx),
    }


def _cache_path(spec):
    digest = sha256(_render(spec).encode("utf-8")).hexdigest()[:16]
    return CACHE_DIR / (
        f"w{int(round(spec['width_nm'])):03d}_s{spec['selectivity']:.3f}_"
        f"{spec['scenario']['name']}_{digest}.json"
    )


def _execute(payload):
    job, transport_device = payload
    document = _load(PREREGISTRATION)
    preregistration = _physics_preregistration(document)
    width, scenario, rates, selectivity, duration, dx = job
    spec = _job_spec(job)
    try:
        profiles = _run_trajectory(
            width_nm=width,
            scenario=scenario,
            rates_nm_min=rates,
            selectivity=selectivity,
            duration_s=duration,
            dx_nm=dx,
            preregistration=preregistration,
            transport_device=transport_device,
        )
    except Exception as error:
        raise RuntimeError(
            f"exact-GDS square trajectory failed for {_render(spec).strip()}"
        ) from error
    path = _cache_path(spec)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render({
        "job_spec": spec,
        "execution": {"transport_device": str(transport_device)},
        "profiles": profiles,
    }), encoding="utf-8")
    return profiles, path, str(transport_device)


def _jobs(*, smoke=False, shard_count=1, shard_index=0):
    shard_count = int(shard_count)
    shard_index = int(shard_index)
    if shard_count <= 0 or shard_index < 0 or shard_index >= shard_count:
        raise ValueError("require 0 <= shard_index < shard_count")
    document = _load(PREREGISTRATION)
    preregistration = _physics_preregistration(document)
    analog = _load(ANALOG_BOARD)
    reactor = _load(REACTOR_DOSE)
    scenarios = _scenario_inputs(preregistration, reactor)
    rates = (
        float(analog["source_feature_depth_board"][
            "minimum_implied_rate_nm_min"]),
        float(analog["source_feature_depth_board"][
            "maximum_implied_rate_nm_min"]),
    )
    selectivities = tuple(
        float(value) for value in preregistration["surface_response_axes"][
            "tio2_to_cr_selectivity"]
    )
    if smoke:
        return [(
            float(preregistration["inferred_geometry_board"]["width_nm"][0]),
            scenarios[0],
            (rates[0],),
            selectivities[0],
            12.0,
            20.0,
        )], rates, selectivities, scenarios
    widths = [
        width for index, width in enumerate(
            preregistration["inferred_geometry_board"]["width_nm"])
        if index % shard_count == shard_index
    ]
    jobs = [
        (width, scenario, rates, selectivity, 1200.0,
         PRODUCTION_MESH_SPACING_NM)
        for width in widths
        for scenario in scenarios
        for selectivity in selectivities
    ]
    return jobs, rates, selectivities, scenarios


def build(*, smoke=False, transport_device="cpu", workers=None,
          shard_count=1, shard_index=0):
    document = _load(PREREGISTRATION)
    jobs, rates, selectivities, scenarios = _jobs(
        smoke=smoke, shard_count=shard_count, shard_index=shard_index)
    groups = []
    receipts = []
    missing = []
    for job in jobs:
        spec = _job_spec(job)
        path = _cache_path(spec)
        if path.exists():
            cached = _load(path)
            if cached.get("job_spec") != spec:
                raise RuntimeError(f"exact-GDS cache mismatch: {path}")
            groups.append(cached["profiles"])
            receipts.append({
                "path": str(path.relative_to(ROOT)),
                "sha256": _hash(path),
                "transport_device": cached.get("execution", {}).get(
                    "transport_device", "unrecorded"),
            })
        else:
            groups.append(None)
            receipts.append(None)
            missing.append((len(groups) - 1, job))
    if missing:
        if workers is None:
            workers = 1 if str(transport_device).startswith("cuda") else 4
        workers = int(workers)
        if workers <= 0:
            raise ValueError("workers must be a positive integer")
        payloads = [(job, str(transport_device)) for _, job in missing]
        if workers == 1:
            computed = [_execute(payload) for payload in payloads]
        else:
            with ProcessPoolExecutor(
                max_workers=min(workers, len(missing)),
                **_process_pool_options(transport_device),
            ) as pool:
                computed = list(pool.map(_execute, payloads))
        for (index, _), (profiles, path, device) in zip(missing, computed):
            groups[index] = profiles
            receipts[index] = {
                "path": str(path.relative_to(ROOT)),
                "sha256": _hash(path),
                "transport_device": str(device),
            }
    profiles = [profile for group in groups for profile in group]
    return {
        "schema": "petch.zhu-npg80-exact-gds-square-profile-board.v1",
        "condition_id": document["condition_id"],
        "model_revision": MODEL_REVISION,
        "smoke_only": bool(smoke),
        "target_sem_used": False,
        "target_depth_used": False,
        "coefficient_selected_from_target": None,
        "geometry": {
            "gds_sha256": document["geometry_source"]["sha256"],
            "pitch_nm": document["exact_layout_geometry"]["pitch_nm"],
            "square_width_nm": document["exact_layout_geometry"][
                "square_width_nm"],
            "mask_polarity_assumption": document["exact_layout_geometry"][
                "mask_polarity_assumption"],
            "mask_polarity_confirmed_by_operator": False,
        },
        "mesh_spacing_nm": 20.0 if smoke else PRODUCTION_MESH_SPACING_NM,
        "moving_materials": ["ALD TiO2", "Cr hard mask"],
        "conditional_axes": {
            "tio2_rate_nm_min": list(rates),
            "tio2_to_cr_selectivity": list(selectivities),
            "transport_scenarios": scenarios,
        },
        "claim_boundary": document["claim_boundary"],
        "execution": {
            "trajectory_transport_devices": sorted(set(
                receipt["transport_device"] for receipt in receipts
            )),
            "execution_device_not_part_of_physics_spec": True,
            "shard_count": int(shard_count),
            "shard_index": int(shard_index),
            "complete_layout_board": int(shard_count) == 1,
        },
        "trajectory_receipts": receipts,
        "profiles": profiles,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument(
        "--transport-device",
        default=os.environ.get("PETCH_TRANSPORT_DEVICE", "cpu"),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=(
            int(os.environ["PETCH_PROFILE_WORKERS"])
            if "PETCH_PROFILE_WORKERS" in os.environ else None
        ),
    )
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    args = parser.parse_args()
    if args.smoke:
        print(_render(build(
            smoke=True,
            transport_device=args.transport_device,
            workers=args.workers,
            shard_count=args.shard_count,
            shard_index=args.shard_index,
        )))
        return
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build(
        smoke=False,
        transport_device=args.transport_device,
        workers=args.workers,
        shard_count=args.shard_count,
        shard_index=args.shard_index,
    ))
    if args.write:
        output = OUTPUT
        if args.shard_count != 1:
            output = OUTPUT.with_name(
                f"audit_shard_{args.shard_index}_of_{args.shard_count}.json")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
        print(output.relative_to(ROOT))
        return
    if args.shard_count != 1 or args.shard_index != 0:
        raise SystemExit("--check requires the complete unsharded board")
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("exact-GDS square profile audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
