#!/usr/bin/env python3
"""Decompose corrected C3 voltage motion and score a stochastic-tail confirmation.

The dominant direction is learned only from the four consecutive fixed-physical-time
reference windows.  The pseudo-time tail and its independent fixed-step confirmation are
then projected onto that archived direction.  This keeps the diagnostic separate from the
integrator and prevents the confirmation from defining the direction used to score itself.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


REFERENCE_NAMES = (
    "periodicized_warm_proposal.npz",
    "c3_periodic_warm_terminal50us_l11/face_checkpoint.npz",
    "c3_periodic_warm_50to100us_l11/face_checkpoint.npz",
    "c3_periodic_warm_100to150us_l11/face_checkpoint.npz",
    "c3_periodic_warm_150to200us_l11/face_checkpoint.npz",
)
TAIL_NAME = "c3_periodic_decreasing_gain_tail_400_l11/face_checkpoint.npz"
CONFIRMATION_NAME = (
    "c3_periodic_decreasing_gain_fixed_confirmation_50us_l11/face_checkpoint.npz"
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _potential(path: Path) -> np.ndarray:
    with np.load(path) as data:
        return np.asarray(data["potential_v"], dtype=float).ravel()


def _decomposition(delta: np.ndarray, mode: np.ndarray) -> dict[str, float]:
    projection = float(delta @ mode)
    orthogonal = delta - projection * mode
    norm = float(np.linalg.norm(delta))
    return {
        "maximum_absolute_voltage_change_v": float(np.max(np.abs(delta))),
        "l2_voltage_change_v": norm,
        "dominant_mode_projection_v": projection,
        "orthogonal_l2_voltage_change_v": float(np.linalg.norm(orthogonal)),
        "dominant_mode_cosine": projection / norm if norm else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    paths = tuple(root / name for name in REFERENCE_NAMES)
    tail_path = root / TAIL_NAME
    confirmation_path = root / CONFIRMATION_NAME
    required = paths + (tail_path, confirmation_path)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        parser.error(f"missing required checkpoints: {missing}")

    reference = tuple(_potential(path) for path in paths)
    if len({item.shape for item in reference}) != 1:
        parser.error("reference checkpoint potential arrays do not share one shape")
    displacements = np.stack([
        reference[index + 1] - reference[index]
        for index in range(len(reference) - 1)
    ])
    _, singular_values, right_vectors = np.linalg.svd(displacements, full_matrices=False)
    dominant_mode = right_vectors[0].copy()
    if float(displacements[0] @ dominant_mode) > 0.0:
        dominant_mode *= -1.0

    tail = _potential(tail_path)
    confirmation = _potential(confirmation_path)
    if tail.shape != reference[-1].shape or confirmation.shape != tail.shape:
        parser.error("tail and confirmation potentials do not match the reference shape")

    physical_rows = []
    for index, delta in enumerate(displacements):
        physical_rows.append({
            "window_us": [50 * index, 50 * (index + 1)],
            **_decomposition(delta, dominant_mode),
        })
    tail_delta = tail - reference[-1]
    confirmation_delta = confirmation - tail
    combined_delta = confirmation - reference[-1]
    tail_result = _decomposition(tail_delta, dominant_mode)
    confirmation_result = _decomposition(confirmation_delta, dominant_mode)
    combined_result = _decomposition(combined_delta, dominant_mode)
    same_direction = (
        tail_result["dominant_mode_projection_v"]
        * confirmation_result["dominant_mode_projection_v"] > 0.0
    )

    output = {
        "schema": "petch.charging.c3.closure-mode-audit.v1",
        "status": "tail did not overshoot; corrected physical drift remains resolved"
        if same_direction else
        "confirmation reversed the tail direction; tail overshoot requires review",
        "method": {
            "reference": "SVD of four consecutive 50 microsecond fixed physical windows",
            "orientation": "first physical-window projection is negative",
            "confirmation_excluded_from_mode_fit": True,
        },
        "dominant_reference_energy_fraction": float(
            singular_values[0] ** 2 / np.sum(singular_values ** 2)),
        "reference_singular_values_v": singular_values.tolist(),
        "physical_windows": physical_rows,
        "decreasing_gain_tail_from_200us": tail_result,
        "fixed_confirmation_from_tail": confirmation_result,
        "combined_from_200us": combined_result,
        "decision": {
            "confirmation_matches_tail_dominant_direction": bool(same_direction),
            "tail_overshoot_detected": bool(not same_direction),
            "canonical_warm_checkpoint": CONFIRMATION_NAME
            if same_direction else REFERENCE_NAMES[-1],
            "formal_saturation_claim": False,
            "reason": (
                "The confirmation preserves the charging direction rather than restoring the "
                "pre-tail state. It validates the branch but also demonstrates continued mean "
                "motion, so it cannot certify stationarity."
            ),
        },
        "checkpoint_sha256": {
            str(path.relative_to(root)): _hash(path) for path in required
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
