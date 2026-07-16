#!/usr/bin/env python3
"""Stop an unattended C3 tmux campaign after a declared completed-segment budget."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import time


TERMINAL_STATUSES = {
    "converged", "campaign_bound_reached", "hard_failure",
    "stopped_at_declared_compute_bound", "stopped_after_discretization_diagnosis",
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--tmux-session", required=True)
    parser.add_argument("--completed-segments", type=int, required=True)
    parser.add_argument("--poll-seconds", type=float, default=5.0)
    args = parser.parse_args()
    if (args.completed_segments <= 0 or args.poll_seconds <= 0.0
            or not args.tmux_session.strip()):
        parser.error("invalid campaign-bound watcher controls")
    campaign = args.campaign_dir.resolve()
    status_path = campaign / "campaign_status.json"
    while True:
        if not status_path.exists():
            time.sleep(args.poll_seconds)
            continue
        status = json.loads(status_path.read_text())
        records = status.get("records", [])
        if status.get("status") in TERMINAL_STATUSES:
            return 0
        if len(records) >= args.completed_segments:
            break
        time.sleep(args.poll_seconds)

    stopped_utc = _utc_now()
    subprocess.run(
        ["tmux", "send-keys", "-t", args.tmux_session, "C-c"], check=False)
    for _ in range(60):
        alive = subprocess.run(
            ["tmux", "has-session", "-t", args.tmux_session],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        if alive.returncode != 0:
            break
        time.sleep(0.5)
    status = json.loads(status_path.read_text())
    records = status.get("records", [])
    if len(records) < args.completed_segments:
        raise RuntimeError("campaign record count regressed while applying its compute bound")
    final = records[args.completed_segments - 1]
    final_checkpoint = campaign / final["checkpoint"]
    stop_record = dict(
        schema="petch.charging.c3.declared-compute-bound-stop.v1",
        stopped_utc=stopped_utc,
        completed_segment_budget=args.completed_segments,
        completed_segments=len(records),
        tmux_session=args.tmux_session,
        final_checkpoint=final["checkpoint"],
        final_checkpoint_sha256=_hash(final_checkpoint),
        cumulative_physical_time_s=final["cumulative_physical_time_s"])
    _atomic_json(campaign / "bounded_stop.json", stop_record)
    status.update(
        status="stopped_at_declared_compute_bound",
        updated_utc=_utc_now(), active_segment=None,
        completed_segment_budget=args.completed_segments,
        cumulative_physical_time_s=final["cumulative_physical_time_s"],
        final_checkpoint=final["checkpoint"],
        bounded_stop="bounded_stop.json")
    _atomic_json(status_path, status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
