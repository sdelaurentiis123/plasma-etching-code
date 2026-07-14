#!/usr/bin/env python3
"""Run durable C3 charging segments until saturation or a declared campaign bound.

The physics runner remains the authority for each segment.  This supervisor only provides
process-level durability: every completed or safely refused segment owns a face checkpoint,
and ``campaign_status.json`` is atomically replaced after inspecting its provenance-bearing
summary.  A charged-cascade work-horizon refusal is the sole auto-recovery currently classified;
it doubles the initial bounce budget up to the declared emergency ceiling and resumes the exact
saved state/epoch.  Conservation, geometry-certification, corrupt-state, and unknown failures
remain hard stops.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


RECOVERABLE_BOUNCE_TEXT = (
    "charged surface-response cascade reached its bounce cap with explicit unresolved charge")


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--initial-face-state", type=Path, required=True)
    parser.add_argument("--method-map", type=Path, required=True)
    parser.add_argument("--base-physical-time-s", type=float, required=True)
    parser.add_argument("--steps-per-segment", type=int, default=50)
    parser.add_argument("--maximum-segments", type=int, default=100)
    parser.add_argument("--initial-response-max-bounces", type=int, default=512)
    parser.add_argument("--emergency-response-max-bounces", type=int, default=1024)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--transport-device", default="cuda:0")
    args = parser.parse_args()
    if (args.steps_per_segment <= 0 or args.maximum_segments <= 0
            or args.initial_response_max_bounces <= 0
            or args.emergency_response_max_bounces < args.initial_response_max_bounces
            or args.base_physical_time_s < 0.0):
        parser.error("invalid campaign bounds")

    root = Path(__file__).resolve().parents[1]
    runner = root / "scripts/charging_coevolution_c3_trench.py"
    campaign = args.campaign_dir.resolve()
    campaign.mkdir(parents=True, exist_ok=True)
    status_path = campaign / "campaign_status.json"
    checkpoint = args.initial_face_state.resolve()
    method_map = args.method_map.resolve()
    bounce_budget = args.initial_response_max_bounces
    cumulative_time = float(args.base_physical_time_s)
    records = []
    status = {
        "schema": "petch.charging.c3.unattended-supervisor.v1",
        "status": "starting",
        "started_utc": _utc_now(),
        "updated_utc": _utc_now(),
        "engine_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "supervisor_sha256": _hash(Path(__file__).resolve()),
        "initial_face_state": checkpoint.name,
        "initial_face_state_sha256": _hash(checkpoint),
        "method_map_sha256": _hash(method_map),
        "base_physical_time_s": cumulative_time,
        "steps_per_segment": args.steps_per_segment,
        "maximum_segments": args.maximum_segments,
        "initial_response_max_bounces": bounce_budget,
        "emergency_response_max_bounces": args.emergency_response_max_bounces,
        "seed": args.seed,
        "transport_device": args.transport_device,
        "python": platform.python_version(),
        "records": records,
    }
    _atomic_json(status_path, status)

    for segment in range(args.maximum_segments):
        output = campaign / f"segment_{segment:04d}"
        command = [
            sys.executable, str(runner), "--output-dir", str(output),
            "--initial-face-state", str(checkpoint),
            "--method-map", str(method_map), "--method-key", "refined_method_hint_Ar+",
            "--timestep-s", "1.25e-7", "--maximum-steps", str(args.steps_per_segment),
            "--timestep-policy", "fixed", "--forward-level", "11",
            "--adjoint-level", "9", "--electron-estimator", "forward",
            "--n-position", "256", "--seed", str(args.seed),
            "--scramble-mode", "fresh", "--sampling-seed-stride", "1000003",
            "--trajectory-dt", "0.000078125", "--trajectory-max-steps", "4096000",
            "--transport-device", args.transport_device,
            "--response-max-bounces", str(bounce_budget),
            "--response-tail-tolerance", "1e-10",
            "--response-launch-offset", "5e-6",
        ]
        output.mkdir(parents=True, exist_ok=True)
        status.update(status="running", updated_utc=_utc_now(), active_segment=segment,
                      cumulative_physical_time_s=cumulative_time,
                      active_response_max_bounces=bounce_budget)
        _atomic_json(status_path, status)
        with (output / "process.log").open("w") as stream:
            completed = subprocess.run(
                command, cwd=root, stdout=stream, stderr=subprocess.STDOUT, check=False)
        summary_path = output / "summary.json"
        if not summary_path.exists():
            status.update(
                status="hard_failure", updated_utc=_utc_now(),
                failure="physics runner exited without a replayable summary",
                returncode=completed.returncode)
            _atomic_json(status_path, status)
            return 2
        summary = json.loads(summary_path.read_text())
        result = summary["result"]
        segment_time = float(result.get("physical_time_s", 0.0))
        cumulative_time += segment_time
        next_checkpoint = output / "face_checkpoint.npz"
        record = {
            "segment": segment,
            "returncode": completed.returncode,
            "summary": summary_path.relative_to(campaign).as_posix(),
            "summary_sha256": _hash(summary_path),
            "checkpoint": next_checkpoint.relative_to(campaign).as_posix(),
            "checkpoint_sha256": _hash(next_checkpoint),
            "segment_physical_time_s": segment_time,
            "cumulative_physical_time_s": cumulative_time,
            "accepted_steps": result.get("accepted_steps"),
            "resume_sampling_epoch": result.get("resume_sampling_epoch"),
            "converged": bool(result.get("converged", False)),
            "failed": bool(result.get("failed", False)),
            "response_max_bounces": bounce_budget,
            "node_rms": result.get("retained_node_rms_relative_current_imbalance"),
            "node_worst": result.get("retained_node_max_relative_current_imbalance"),
            "potential_rate_max_v_s": result.get("final_potential_rate_max_v_s"),
            "patch_b2_max": result.get("patch_b2_max_ion_normalized"),
        }
        records.append(record)
        checkpoint = next_checkpoint
        if record["converged"] and completed.returncode == 0:
            status.update(
                status="converged", updated_utc=_utc_now(), active_segment=None,
                cumulative_physical_time_s=cumulative_time,
                final_checkpoint=record["checkpoint"])
            _atomic_json(status_path, status)
            return 0
        if completed.returncode == 0:
            status.update(updated_utc=_utc_now(), cumulative_physical_time_s=cumulative_time)
            _atomic_json(status_path, status)
            continue

        message = str(result.get("error_message", ""))
        if (RECOVERABLE_BOUNCE_TEXT in message
                and bounce_budget < args.emergency_response_max_bounces):
            bounce_budget = min(2 * bounce_budget, args.emergency_response_max_bounces)
            record["recovery"] = "resume exact checkpoint with doubled bounce budget"
            record["next_response_max_bounces"] = bounce_budget
            status.update(
                status="recovering_bounce_budget", updated_utc=_utc_now(),
                cumulative_physical_time_s=cumulative_time,
                active_response_max_bounces=bounce_budget)
            _atomic_json(status_path, status)
            continue
        status.update(
            status="hard_failure", updated_utc=_utc_now(),
            cumulative_physical_time_s=cumulative_time,
            failure=message or f"unclassified runner exit {completed.returncode}")
        _atomic_json(status_path, status)
        return 2

    status.update(
        status="campaign_bound_reached", updated_utc=_utc_now(), active_segment=None,
        cumulative_physical_time_s=cumulative_time,
        final_checkpoint=records[-1]["checkpoint"] if records else checkpoint.name)
    _atomic_json(status_path, status)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
