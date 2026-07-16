#!/usr/bin/env python3
"""Pair exact C3 projective/PTC candidate scores by unused sampling epoch."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t


METRICS = (
    "residual_current_norm_a",
    "rms_relative_current_imbalance_node",
    "max_relative_current_imbalance_node",
    "potential_rate_max_v_s",
    "maximum_patch_relative_imbalance",
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def _parse_pattern(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("candidate patterns must be LABEL=GLOB")
    label, pattern = value.split("=", 1)
    if not label or not pattern:
        raise argparse.ArgumentTypeError("candidate patterns must be LABEL=GLOB")
    return label, pattern


def _load(pattern: str) -> tuple[dict[int, dict[str, float]], list[dict[str, str]]]:
    paths = sorted(Path().glob(pattern))
    if len(paths) < 2:
        raise ValueError(f"pattern needs at least two summaries: {pattern}")
    records: dict[int, dict[str, float]] = {}
    provenance = []
    for path in paths:
        payload = json.loads(path.read_text())
        history = payload.get("history", [])
        if len(history) != 1 or int(payload["result"]["accepted_steps"]) != 0:
            raise ValueError(f"score is not a zero-step exact audit: {path}")
        record = history[0]
        epoch = int(record["sampling_epoch"])
        if epoch in records:
            raise ValueError(f"duplicate sampling epoch {epoch} for {pattern}")
        records[epoch] = {name: float(record[name]) for name in METRICS}
        provenance.append(dict(
            source_directory=path.parent.name,
            name=path.name,
            sha256=_hash(path)))
    return records, provenance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-pattern", required=True)
    parser.add_argument(
        "--candidate-pattern", action="append", type=_parse_pattern, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    baseline, baseline_provenance = _load(args.baseline_pattern)
    baseline_epochs = set(baseline)
    candidates = {}
    provenance = {"baseline": baseline_provenance}
    for label, pattern in args.candidate_pattern:
        if label in candidates:
            parser.error(f"duplicate candidate label: {label}")
        records, candidate_provenance = _load(pattern)
        if set(records) != baseline_epochs:
            parser.error(f"candidate {label} does not have the baseline sampling epochs")
        candidates[label] = records
        provenance[label] = candidate_provenance

    epochs = sorted(baseline_epochs)
    count = len(epochs)
    critical = float(student_t.ppf(0.975, count - 1))
    rows = []
    results = {}
    for label, records in candidates.items():
        results[label] = {}
        for metric in METRICS:
            base = np.asarray([baseline[epoch][metric] for epoch in epochs])
            values = np.asarray([records[epoch][metric] for epoch in epochs])
            difference = values - base
            difference_standard_error = float(np.std(difference, ddof=1) / math.sqrt(count))
            row = dict(
                candidate=label,
                metric=metric,
                replicate_count=count,
                candidate_mean=float(np.mean(values)),
                candidate_standard_error=float(np.std(values, ddof=1) / math.sqrt(count)),
                baseline_mean=float(np.mean(base)),
                paired_difference_mean=float(np.mean(difference)),
                paired_difference_95ci_half_width=critical * difference_standard_error,
            )
            rows.append(row)
            results[label][metric] = {
                key: value for key, value in row.items()
                if key not in ("candidate", "metric")
            }

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "paired_metrics.csv"
    temporary = csv_path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, csv_path)
    summary = dict(
        schema="petch.charging.c3.projective-ptc-paired-score-audit.v1",
        scoring_operator="exact hard visibility at zero accepted steps",
        common_random_numbers=True,
        sampling_epochs=epochs,
        replicate_count=count,
        metrics=results,
        provenance=provenance,
        artifact=dict(name=csv_path.name, sha256=_hash(csv_path)))
    summary_path = output / "summary.json"
    _atomic_json(summary_path, summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
