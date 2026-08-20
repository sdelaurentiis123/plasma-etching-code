#!/usr/bin/env python3
"""Prove whether an axisymmetric Bosch map can meet the frozen shape gate.

This audit opens only the extracted calibration-only 89-point asset.  It
computes oracle lower bounds that are more favorable than any physical
axisymmetric reactor: one shared free value per measured radius, and a separate
free radial curve for every calibration wafer.  If those oracles fail the
frozen gate, the axisymmetric model form is falsified independently of any
equipment-parameter optimizer.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.bosch_process_data import load_bosch_wafer_measurements_89pt


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
CALIBRATION = DATA / "calibration_Si_Oxide_etch_89_points.csv"
CALIBRATION_MANIFEST = (
    DATA / "calibration_Si_Oxide_etch_89_points_manifest.json")
V1_PREREGISTRATION = (
    ROOT / "results" / "curated" / "zenodo_bosch_reactor_depth_holdout_v1"
    / "preregistration.json")
OUTPUT = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_cylindrical_depth_extension_v2"
    / "axisymmetry_model_form_audit.json")


def _sha256(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _rmse_percent(actual, expected):
    return 100.0 * float(np.sqrt(np.mean((actual - expected) ** 2)))


def build_audit():
    preregistration = json.loads(V1_PREREGISTRATION.read_text())
    manifest = json.loads(CALIBRATION_MANIFEST.read_text())
    allowed = tuple(preregistration["split_rule"]["calibration_experiment_keys"])
    heldout = frozenset(preregistration["split_rule"]["heldout_experiment_keys"])

    with CALIBRATION.open(newline="", encoding="utf-8") as stream:
        present = frozenset(row["experiment_key"] for row in csv.DictReader(stream))
    if present & heldout:
        raise RuntimeError("calibration model-form audit encountered a heldout key")
    if manifest["heldout_rows_copied_to_fit_asset"]:
        raise RuntimeError("calibration target firewall is not closed")
    measurements = load_bosch_wafer_measurements_89pt(
        CALIBRATION, allowed_experiment_keys=allowed)
    if len(measurements) != manifest["output_experiment_key_count"]:
        raise RuntimeError("calibration manifest count mismatch")

    x = measurements[0].x_um
    y = measurements[0].y_um
    if any(
            not np.array_equal(item.x_um, x)
            or not np.array_equal(item.y_um, y)
            for item in measurements):
        raise RuntimeError("calibration wafers do not share one coordinate map")
    silicon = np.stack([item.silicon_depth_um for item in measurements])
    oxide = np.stack([item.oxide_loss_um for item in measurements])
    normalized = silicon / np.mean(silicon, axis=1, keepdims=True)
    radius_squared = x * x + y * y
    unique_radius_squared = np.unique(radius_squared)

    shared_axisymmetric_oracle = np.empty(x.size)
    per_wafer_axisymmetric_oracle = np.empty_like(normalized)
    for radius in unique_radius_squared:
        point = radius_squared == radius
        shared_axisymmetric_oracle[point] = float(np.mean(normalized[:, point]))
        per_wafer_axisymmetric_oracle[:, point] = np.mean(
            normalized[:, point], axis=1, keepdims=True)
    per_wafer_rmse = np.asarray([
        _rmse_percent(actual, expected)
        for actual, expected in zip(normalized, per_wafer_axisymmetric_oracle)
    ])
    shared_coordinate_oracle = np.mean(normalized, axis=0)
    gate = float(
        preregistration["heldout_score"]["absolute_acceptance"]
        ["normalized_radial_shape_rmse_percent_max"])

    return {
        "schema": "petch-bosch-axisymmetry-model-form-audit-v1",
        "source_record": "https://zenodo.org/records/17122442",
        "calibration_asset": str(CALIBRATION.relative_to(ROOT)),
        "calibration_asset_sha256": _sha256(CALIBRATION),
        "calibration_manifest_sha256": _sha256(CALIBRATION_MANIFEST),
        "v1_preregistration_sha256": _sha256(V1_PREREGISTRATION),
        "heldout_outcomes_read": False,
        "calibration_wafer_count": len(measurements),
        "points_per_wafer": int(x.size),
        "unique_measured_radius_count": int(unique_radius_squared.size),
        "normalization": "each silicon map divided by its 89-point arithmetic mean",
        "weighting": "equal weight per measured point and calibration wafer",
        "calibration_measurement_summary": {
            "silicon_mean_depth_um_min": float(np.min(np.mean(silicon, axis=1))),
            "silicon_mean_depth_um_median": float(np.median(np.mean(silicon, axis=1))),
            "silicon_mean_depth_um_max": float(np.max(np.mean(silicon, axis=1))),
            "oxide_mean_loss_um_min": float(np.min(np.mean(oxide, axis=1))),
            "oxide_mean_loss_um_median": float(np.median(np.mean(oxide, axis=1))),
            "oxide_mean_loss_um_max": float(np.max(np.mean(oxide, axis=1))),
        },
        "oracle_lower_bounds": {
            "shared_axisymmetric_normalized_rmse_percent": _rmse_percent(
                normalized, shared_axisymmetric_oracle[None, :]),
            "separate_axisymmetric_per_wafer_rmse_percent_min": float(
                np.min(per_wafer_rmse)),
            "separate_axisymmetric_per_wafer_rmse_percent_median": float(
                np.median(per_wafer_rmse)),
            "separate_axisymmetric_per_wafer_rmse_percent_max": float(
                np.max(per_wafer_rmse)),
            "shared_unconstrained_coordinate_map_normalized_rmse_percent": (
                _rmse_percent(normalized, shared_coordinate_oracle[None, :])),
        },
        "frozen_shape_gate_percent": gate,
        "shared_axisymmetric_model_can_pass_frozen_gate": bool(
            _rmse_percent(normalized, shared_axisymmetric_oracle[None, :]) <= gate),
        "interpretation": (
            "Even a nonphysical oracle shared radial curve exceeds the frozen "
            "shape gate, while a shared unconstrained x-y map is below it. "
            "Axisymmetry is therefore a falsified model form for this gate."),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--print", action="store_true", dest="print_payload")
    args = parser.parse_args(argv)
    payload = (json.dumps(build_audit(), indent=2, sort_keys=True) + "\n")
    if args.print_payload:
        print(payload, end="")
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != payload:
            raise SystemExit("committed Bosch axisymmetry audit is stale")
        return 0
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
