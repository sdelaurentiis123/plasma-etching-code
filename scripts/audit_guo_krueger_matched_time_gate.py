#!/usr/bin/env python3
"""Grade the nonlinear 4 s Guo/Krueger medium/fine time-step pair."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


DEFAULT_MEDIUM = Path("/private/tmp/krueger_guo_transient_dt125_dx10")
DEFAULT_FINE = Path("/private/tmp/krueger_guo_transient_dt0625_dx10_to4s")
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "curated"
    / "guo_krueger_finite_fluence_4s_time_gate"
    / "audit.json"
)
MATCH_START_S = 0.5
MATCH_STOP_S = 4.0


def interpolate_metrics(history, physical_time_s):
    """Linearly sample geometry metrics from accepted-step history."""
    time = float(physical_time_s)
    if not math.isfinite(time):
        raise ValueError("matched time must be finite")
    rows = sorted(history, key=lambda row: float(row["physical_time_s"]))
    if not rows or time < float(rows[0]["physical_time_s"]) or time > float(
        rows[-1]["physical_time_s"]
    ):
        raise ValueError("matched time leaves the accepted trajectory")
    for row in rows:
        if math.isclose(
            time, float(row["physical_time_s"]), rel_tol=0.0, abs_tol=1e-13
        ):
            return {
                "etch_depth_nm": float(row["metrics"]["etch_depth_nm"]),
                "mask_opening_nm": float(row["metrics"]["mask_opening_nm"]),
            }
    for left, right in zip(rows[:-1], rows[1:]):
        left_time = float(left["physical_time_s"])
        right_time = float(right["physical_time_s"])
        if left_time < time < right_time:
            fraction = (time - left_time) / (right_time - left_time)
            return {
                name: (
                    float(left["metrics"][name])
                    + fraction
                    * (
                        float(right["metrics"][name])
                        - float(left["metrics"][name])
                    )
                )
                for name in ("etch_depth_nm", "mask_opening_nm")
            }
    raise RuntimeError("accepted trajectory is not a continuous time bracket")


def _load(path):
    audit_path = Path(path) / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if float(audit["history"][-1]["physical_time_s"]) < MATCH_STOP_S:
        raise ValueError(f"{audit_path} has not reached the 4 s gate")
    return audit


def _maximum(history, name):
    return max(float(row.get(name, 0.0)) for row in history)


def _history_through_gate(audit):
    return [
        row
        for row in audit["history"]
        if float(row["physical_time_s"]) <= MATCH_STOP_S + 1.0e-13
    ]


def build_report(medium_path=DEFAULT_MEDIUM, fine_path=DEFAULT_FINE):
    medium = _load(medium_path)
    fine = _load(fine_path)
    if not math.isclose(
        float(medium["configuration"]["dx_um"]),
        float(fine["configuration"]["dx_um"]),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise ValueError("matched-time pair changes spatial resolution")
    medium_dt = (
        float(medium["configuration"]["duration_s"])
        / int(medium["configuration"]["n_steps"])
    )
    fine_dt = (
        float(fine["configuration"]["duration_s"])
        / int(fine["configuration"]["n_steps"])
    )
    if not math.isclose(medium_dt / fine_dt, 2.0, rel_tol=0.0, abs_tol=1e-14):
        raise ValueError("matched-time pair is not a twofold time refinement")

    matched_times = sorted(
        {
            float(row["physical_time_s"])
            for row in fine["history"]
            if MATCH_START_S <= float(row["physical_time_s"]) <= MATCH_STOP_S
        }
        | {MATCH_START_S, MATCH_STOP_S}
    )
    trajectory = []
    for time in matched_times:
        coarse = interpolate_metrics(medium["history"], time)
        refined = interpolate_metrics(fine["history"], time)
        trajectory.append(
            {
                "physical_time_s": time,
                "medium": coarse,
                "fine": refined,
                "signed_depth_relative": (
                    coarse["etch_depth_nm"] / refined["etch_depth_nm"] - 1.0
                ),
                "signed_mask_opening_relative": (
                    coarse["mask_opening_nm"]
                    / refined["mask_opening_nm"]
                    - 1.0
                ),
                "mask_opening_difference_nm": (
                    coarse["mask_opening_nm"] - refined["mask_opening_nm"]
                ),
            }
        )

    terminal = trajectory[-1]
    maximum_depth_relative = max(
        abs(row["signed_depth_relative"]) for row in trajectory
    )
    rms_depth_relative = math.sqrt(
        sum(row["signed_depth_relative"] ** 2 for row in trajectory)
        / len(trajectory)
    )
    maximum_mask_difference_nm = max(
        abs(row["mask_opening_difference_nm"]) for row in trajectory
    )
    gate_histories = tuple(
        _history_through_gate(audit) for audit in (medium, fine)
    )
    maximum_ledger = max(
        _maximum(history, "maximum_material_ledger_residual_units_m2")
        for history in gate_histories
    )
    maximum_radiosity = max(
        _maximum(
            history,
            "maximum_neutral_radiosity_relative_balance_error",
        )
        for history in gate_histories
    )
    topology_event_count = sum(
        sum(
            float(event.get("physical_time_lower_s", math.inf))
            < MATCH_STOP_S
            for event in audit["topology_events"]
        )
        for audit in (medium, fine)
    )
    gates = {
        "terminal_depth_abs_relative_le_5pct": {
            "limit": 0.05,
            "value": abs(terminal["signed_depth_relative"]),
        },
        "terminal_mask_abs_relative_le_5pct": {
            "limit": 0.05,
            "value": abs(terminal["signed_mask_opening_relative"]),
        },
        "maximum_matched_depth_abs_relative_le_5pct": {
            "limit": 0.05,
            "value": maximum_depth_relative,
        },
        "maximum_matched_mask_difference_le_5nm": {
            "limit": 5.0,
            "value": maximum_mask_difference_nm,
        },
        "material_ledger_exact": {
            "limit": 0.0,
            "value": maximum_ledger,
        },
        "neutral_radiosity_balance_le_1e_9": {
            "limit": 1.0e-9,
            "value": maximum_radiosity,
        },
        "topology_event_count_zero": {
            "limit": 0,
            "value": topology_event_count,
        },
    }
    for gate in gates.values():
        gate["passed"] = gate["value"] <= gate["limit"]

    return {
        "schema": "petch-guo-krueger-4s-matched-time-gate-v1",
        "claim": (
            "Twofold time-step comparison through the nonlinear 2.9 s mouth "
            "transition; numerical authorization only, not a physical-boundary "
            "validation."
        ),
        "medium_path": str(medium_path),
        "fine_path": str(fine_path),
        "dx_nm": 1000.0 * float(medium["configuration"]["dx_um"]),
        "medium_nominal_step_s": medium_dt,
        "fine_nominal_step_s": fine_dt,
        "match_interval_s": [MATCH_START_S, MATCH_STOP_S],
        "matched_point_count": len(trajectory),
        "terminal": terminal,
        "maximum_matched_depth_abs_relative": maximum_depth_relative,
        "rms_matched_depth_relative": rms_depth_relative,
        "maximum_matched_mask_difference_nm": maximum_mask_difference_nm,
        "gates": gates,
        "all_gates_passed": all(gate["passed"] for gate in gates.values()),
        "trajectory": trajectory,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--medium", type=Path, default=DEFAULT_MEDIUM)
    parser.add_argument("--fine", type=Path, default=DEFAULT_FINE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = build_report(args.medium, args.fine)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"stale audit: {args.output}")
        print(f"verified {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
