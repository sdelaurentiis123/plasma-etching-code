#!/usr/bin/env python3
"""Freeze the chronological Bosch reactor-depth split without reading outcomes."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
DEFAULT_SUMMARY = DATA / "process_wafer_summary.csv"
DEFAULT_MANIFEST = (
    ROOT / "results" / "curated" / "zenodo_bosch_reactor_depth_holdout_v1"
    / "preregistration.json"
)

CALIBRATION_DATES = (
    "2024-07-02", "2024-07-05", "2024-07-09", "2024-07-11",
    "2024-07-19", "2024-08-01", "2024-08-05", "2024-08-07",
)
HELDOUT_DATES = ("2024-08-21", "2024-08-22")


def build_preregistration(summary_path: Path):
    payload = summary_path.read_bytes()
    with summary_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    required = {"experiment_key", "process_date", "wafer_number"}
    if len(rows) != 96 or not rows or required - set(rows[0]):
        raise ValueError("unexpected Bosch process summary")
    if len({row["experiment_key"] for row in rows}) != len(rows):
        raise ValueError("duplicate Bosch process identifier")

    calibration = sorted(
        row["experiment_key"] for row in rows
        if row["process_date"] in CALIBRATION_DATES)
    heldout = sorted(
        row["experiment_key"] for row in rows
        if row["process_date"] in HELDOUT_DATES)
    unused = sorted(
        row["experiment_key"] for row in rows
        if row["process_date"] not in set(CALIBRATION_DATES + HELDOUT_DATES))
    if len(calibration) != 76 or len(heldout) != 20 or unused:
        raise ValueError("Bosch chronological split no longer has its frozen membership")

    return {
        "schema": "petch-zenodo-bosch-reactor-depth-preregistration-v1",
        "source_record": "https://zenodo.org/records/17122442",
        "process_summary_sha256": sha256(payload).hexdigest(),
        "split_rule": {
            "calibration_dates": list(CALIBRATION_DATES),
            "heldout_dates": list(HELDOUT_DATES),
            "reason": (
                "chronological transfer across the final two experiment days and the "
                "intervening chamber-time gap; membership uses process date only"
            ),
            "calibration_process_record_count": len(calibration),
            "heldout_process_record_count": len(heldout),
            "calibration_experiment_keys": calibration,
            "heldout_experiment_keys": heldout,
        },
        "target_firewall": {
            "outcome_files_opened_during_preregistration": False,
            "fit_may_load_only_calibration_experiment_keys": True,
            "heldout_outcomes_may_be_opened_only_after_model_receipt_is_hash_sealed": True,
            "preexposure_blind": False,
            "status": (
                "execution-held-out, not pre-exposure blind: measurement CSVs already "
                "existed in the repository before this split was frozen"
            ),
        },
        "frozen_physics_path": [
            "5 Hz measured gas/pressure/source/platen time series",
            "deterministic phase-resolved zero-dimensional reactor state",
            "measured-waveform-conditioned sheath and ion-energy boundary",
            "axisymmetric wafer transfer",
            "unchanged Belen silicon surface law plus a conserved C4F8 film-memory state",
        ],
        "forbidden_shortcuts": [
            "direct regression from experiment key, lot number, wafer number, date, or target depth",
            "opening heldout silicon depth, oxide loss, selectivity, or radial maps during fitting",
            "per-heldout-wafer flux multipliers or radial-profile coefficients",
            "changing common transport, sheath, surface, or state-evolution operators after reveal",
        ],
        "calibration_policy": {
            "selection": "leave-one-lot-out cross-validation within calibration dates",
            "allowed": [
                "tool-specific absorbed-power coupling",
                "tool-specific residence/wall-loss closure",
                "one initial chamber-memory state per declared conditioning class",
                "axisymmetric equipment transfer coefficients shared by every wafer",
            ],
            "not_allowed": [
                "heldout target conditioning",
                "per-wafer free parameters",
                "surface-yield changes selected from heldout depth",
            ],
        },
        "heldout_score": {
            "bootstrap_unit": "wafer, not individual spatial point",
            "required_observables": [
                "wafer-mean silicon etch depth",
                "89-point silicon depth map",
                "wafer-mean oxide mask loss",
                "silicon-to-oxide selectivity",
                "within-lot depth drift",
            ],
            "absolute_acceptance": {
                "wafer_mean_si_depth_mae_um_max": 1.0,
                "wafer_mean_si_depth_mape_percent_max": 3.0,
                "pointwise_si_depth_rmse_um_max": 1.5,
                "normalized_radial_shape_rmse_percent_max": 2.0,
                "wafer_mean_oxide_loss_mae_um_max": 0.08,
                "selectivity_mape_percent_max": 12.0,
            },
            "baseline_requirement": (
                "the physics path must beat a calibration-global-mean depth predictor and "
                "a calibration-mean radial-map predictor on the same heldout wafers"
            ),
        },
        "claim_boundary": (
            "This gate validates reactor/wafer absolute depth, selectivity, radial transfer, "
            "and drift for one Bosch tool. It does not by itself validate feature charging, "
            "sidewall angle, scallop geometry, or ARDE."
        ),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    rendered = json.dumps(
        build_preregistration(args.summary), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("committed Bosch preregistration is stale")
        print("Bosch reactor-depth preregistration is current")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
