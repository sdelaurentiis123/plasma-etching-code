#!/usr/bin/env python3
"""Grade preregistered deterministic-extruded Guo/Krueger prefix pairs.

This is a numerical gate only.  It neither uses nor grades the experimental
depth endpoint, and passing it cannot promote the published Krueger boundary
to a physically identified reactor boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TIME_COARSE = Path(
    "/private/tmp/krueger_guo_deterministic_extruded_dt0625_dx10_prefix"
)
DEFAULT_TIME_FINE = Path(
    "/private/tmp/krueger_guo_deterministic_extruded_dt03125_dx10_prefix"
)
DEFAULT_SPACE_FINE = Path(
    "/private/tmp/krueger_guo_deterministic_extruded_dt03125_dx5_prefix"
)
DEFAULT_OUTPUT = (
    ROOT
    / "results"
    / "curated"
    / "guo_krueger_deterministic_prefix"
    / "time_gate.json"
)
MATCH_START_S = 0.0625
MATCH_STOP_S = 0.5


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(directory: Path) -> tuple[dict, Path]:
    path = Path(directory) / "audit.json"
    audit = json.loads(path.read_text(encoding="utf-8"))
    if audit.get("status") != "complete":
        raise ValueError(f"{path} is not a completed trajectory")
    if float(audit["history"][-1]["physical_time_s"]) < MATCH_STOP_S:
        raise ValueError(f"{path} has not reached {MATCH_STOP_S} s")
    if (
        audit["configuration"].get("radiosity_backend")
        != "deterministic_extruded_2d"
    ):
        raise ValueError(f"{path} does not use deterministic_extruded_2d")
    return audit, path


def _interpolate_metrics(history: list[dict], time: float) -> dict[str, float]:
    rows = sorted(history, key=lambda row: float(row["physical_time_s"]))
    for row in rows:
        if math.isclose(
            time, float(row["physical_time_s"]), rel_tol=0.0, abs_tol=1.0e-13
        ):
            return {
                name: float(row["metrics"][name])
                for name in ("etch_depth_nm", "mask_opening_nm")
            }
    for left, right in zip(rows[:-1], rows[1:]):
        left_time = float(left["physical_time_s"])
        right_time = float(right["physical_time_s"])
        if left_time < time < right_time:
            fraction = (time - left_time) / (right_time - left_time)
            return {
                name: float(left["metrics"][name])
                + fraction
                * (
                    float(right["metrics"][name])
                    - float(left["metrics"][name])
                )
                for name in ("etch_depth_nm", "mask_opening_nm")
            }
    raise ValueError(f"time {time} s leaves the accepted trajectory")


def _configuration_difference(
    reference: dict, refined: dict, allowed: set[str]
) -> dict[str, dict[str, object]]:
    difference: dict[str, dict[str, object]] = {}
    for key in sorted(set(reference) | set(refined)):
        if reference.get(key) != refined.get(key):
            difference[key] = {
                "reference": reference.get(key),
                "refined": refined.get(key),
                "allowed": key in allowed,
            }
    return difference


def _maximum(history: list[dict], name: str) -> float:
    return max(float(row.get(name, 0.0)) for row in history)


def _health(audits: tuple[dict, dict]) -> dict[str, float]:
    histories = tuple(
        [
            row
            for row in audit["history"]
            if 0.0 < float(row["physical_time_s"]) <= MATCH_STOP_S + 1.0e-13
        ]
        for audit in audits
    )
    ratios = []
    for history in histories:
        for row in history:
            resolved = abs(float(row["max_velocity_m_s"]))
            raw = abs(float(row["raw_maximum_face_velocity_m_s"]))
            ratios.append(raw / resolved if resolved > 0.0 else math.inf)
    return {
        "maximum_material_ledger_residual_units_m2": max(
            _maximum(history, "maximum_material_ledger_residual_units_m2")
            for history in histories
        ),
        "maximum_neutral_radiosity_relative_balance_error": max(
            _maximum(
                history,
                "maximum_neutral_radiosity_relative_balance_error",
            )
            for history in histories
        ),
        "maximum_extrusion_projection_deviation_mesh_units": max(
            float(audit["extrusion_projection_max_deviation_mesh_units"])
            for audit in audits
        ),
        "maximum_raw_to_resolved_speed_ratio": max(ratios),
        "rejected_trial_count": float(
            sum(
                len(row.get("rejected_trials", ()))
                for history in histories
                for row in history
            )
        ),
        "accepted_topology_event_count": float(
            sum(len(audit.get("topology_events", ())) for audit in audits)
        ),
        "maximum_asymmetric_cell_count": max(
            float(row["metrics"].get("asymmetry_cell_count", 0.0))
            for history in histories
            for row in history
        ),
        "maximum_mirrored_node_sign_mismatch_pair_count": max(
            float(row["metrics"].get(
                "mirrored_node_sign_mismatch_pair_count", 0.0))
            for history in histories
            for row in history
        ),
        "maximum_mirrored_material_label_mismatch_pair_count": max(
            float(row["metrics"].get(
                "mirrored_material_label_mismatch_pair_count", 0.0))
            for history in histories
            for row in history
        ),
        "maximum_subcell_interface_asymmetry_cells": max(
            float(row["metrics"].get(
                "maximum_subcell_interface_asymmetry_cells", 0.0))
            for history in histories
            for row in history
        ),
    }


def build_report(
    reference_directory: Path,
    refined_directory: Path,
    comparison: str,
) -> dict[str, object]:
    reference, reference_path = _load(reference_directory)
    refined, refined_path = _load(refined_directory)
    if comparison not in {"time", "space"}:
        raise ValueError("comparison must be 'time' or 'space'")

    allowed = {"maximum_accepted_steps"}
    if comparison == "time":
        allowed.add("n_steps")
        reference_dt = (
            float(reference["configuration"]["duration_s"])
            / int(reference["configuration"]["n_steps"])
        )
        refined_dt = (
            float(refined["configuration"]["duration_s"])
            / int(refined["configuration"]["n_steps"])
        )
        if not math.isclose(
            reference_dt / refined_dt, 2.0, rel_tol=0.0, abs_tol=1.0e-14
        ):
            raise ValueError("time pair is not an exact twofold refinement")
        if reference["configuration"]["dx_um"] != refined["configuration"]["dx_um"]:
            raise ValueError("time pair changes spatial resolution")
        limits = {
            "terminal_depth_abs_relative": 0.02,
            "terminal_mask_opening_abs_relative": 0.02,
            "maximum_matched_depth_abs_relative": 0.05,
            "maximum_matched_mask_opening_difference_nm": 5.0,
        }
    else:
        allowed.add("dx_um")
        if reference["configuration"]["n_steps"] != refined["configuration"]["n_steps"]:
            raise ValueError("space pair changes the nominal time step")
        if not math.isclose(
            float(reference["configuration"]["dx_um"])
            / float(refined["configuration"]["dx_um"]),
            2.0,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ):
            raise ValueError("space pair is not an exact twofold refinement")
        limits = {
            "terminal_depth_abs_relative": 0.05,
            "terminal_mask_opening_difference_nm": 5.0,
            "maximum_matched_depth_abs_relative": 0.075,
            "maximum_matched_mask_opening_difference_nm": 7.5,
        }

    differences = _configuration_difference(
        reference["configuration"], refined["configuration"], allowed
    )
    unauthorized = {
        key: value for key, value in differences.items() if not value["allowed"]
    }
    if unauthorized:
        raise ValueError(f"pair changes unauthorized configuration: {unauthorized}")

    matched_times = sorted(
        {
            float(row["physical_time_s"])
            for row in refined["history"]
            if MATCH_START_S
            <= float(row["physical_time_s"])
            <= MATCH_STOP_S + 1.0e-13
        }
        | {MATCH_START_S, MATCH_STOP_S}
    )
    trajectory = []
    for time in matched_times:
        reference_metrics = _interpolate_metrics(reference["history"], time)
        refined_metrics = _interpolate_metrics(refined["history"], time)
        depth_relative = (
            reference_metrics["etch_depth_nm"] / refined_metrics["etch_depth_nm"]
            - 1.0
        )
        mouth_difference = (
            reference_metrics["mask_opening_nm"]
            - refined_metrics["mask_opening_nm"]
        )
        trajectory.append(
            {
                "physical_time_s": time,
                "reference": reference_metrics,
                "refined": refined_metrics,
                "signed_depth_relative": depth_relative,
                "mask_opening_difference_nm": mouth_difference,
            }
        )

    terminal = trajectory[-1]
    metrics = {
        "terminal_depth_abs_relative": abs(terminal["signed_depth_relative"]),
        "terminal_mask_opening_abs_relative": abs(
            terminal["reference"]["mask_opening_nm"]
            / terminal["refined"]["mask_opening_nm"]
            - 1.0
        ),
        "terminal_mask_opening_difference_nm": abs(
            terminal["mask_opening_difference_nm"]
        ),
        "maximum_matched_depth_abs_relative": max(
            abs(row["signed_depth_relative"]) for row in trajectory
        ),
        "maximum_matched_mask_opening_difference_nm": max(
            abs(row["mask_opening_difference_nm"]) for row in trajectory
        ),
    }
    health = _health((reference, refined))
    health_limits = {
        "maximum_material_ledger_residual_units_m2": 0.0,
        "maximum_neutral_radiosity_relative_balance_error": 1.0e-9,
        "maximum_extrusion_projection_deviation_mesh_units": 1.0e-9,
        "maximum_raw_to_resolved_speed_ratio": 2.0,
        "rejected_trial_count": 0.0,
        "accepted_topology_event_count": 0.0,
        "maximum_asymmetric_cell_count": 0.0,
        "maximum_mirrored_node_sign_mismatch_pair_count": 0.0,
        "maximum_mirrored_material_label_mismatch_pair_count": 0.0,
    }
    gates: dict[str, dict[str, object]] = {}
    for name, limit in limits.items():
        value = metrics[name]
        gates[name] = {
            "limit": limit,
            "value": value,
            "passed": value <= limit,
        }
    for name, limit in health_limits.items():
        value = health[name]
        gates[name] = {
            "limit": limit,
            "value": value,
            "passed": value <= limit,
        }

    return {
        "schema": "petch.guo-krueger.deterministic-prefix-gate.v2",
        "comparison": comparison,
        "claim": (
            "Preregistered deterministic-extruded numerical refinement gate; "
            "the experimental endpoint is not an input and physical prediction "
            "authority is not granted."
        ),
        "reference": {
            "directory": str(reference_directory),
            "audit_sha256": _sha256(reference_path),
            "config_hash": reference["config_hash"],
            "dx_nm": 1000.0 * float(reference["configuration"]["dx_um"]),
            "nominal_step_s": (
                float(reference["configuration"]["duration_s"])
                / int(reference["configuration"]["n_steps"])
            ),
        },
        "refined": {
            "directory": str(refined_directory),
            "audit_sha256": _sha256(refined_path),
            "config_hash": refined["config_hash"],
            "dx_nm": 1000.0 * float(refined["configuration"]["dx_um"]),
            "nominal_step_s": (
                float(refined["configuration"]["duration_s"])
                / int(refined["configuration"]["n_steps"])
            ),
        },
        "configuration_differences": differences,
        "match_interval_s": [MATCH_START_S, MATCH_STOP_S],
        "matched_point_count": len(trajectory),
        "metrics": metrics,
        "health": health,
        "gates": gates,
        "all_gates_passed": all(gate["passed"] for gate in gates.values()),
        "trajectory": trajectory,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison", choices=("time", "space"), default="time")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--refined", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.reference is None:
        args.reference = DEFAULT_TIME_COARSE if args.comparison == "time" else DEFAULT_TIME_FINE
    if args.refined is None:
        args.refined = DEFAULT_TIME_FINE if args.comparison == "time" else DEFAULT_SPACE_FINE
    report = build_report(args.reference, args.refined, args.comparison)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"stale audit: {args.output}")
        print(f"verified {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")
    print(
        f"{args.comparison} gate: "
        f"{'PASS' if report['all_gates_passed'] else 'FAIL'} -> {args.output}"
    )


if __name__ == "__main__":
    main()
