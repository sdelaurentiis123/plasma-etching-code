#!/usr/bin/env python3
"""Audit what Yoshie's blanket rates can and cannot identify.

This is deliberately a diagnostic, not a calibration routine.  It compares
the seven checksum-verified Figure-4 poly-Si blanket rates with all 49
preregistered Figure-5/6 bulk-Si feature rates under the matching nominal
timings.  No model parameter is changed and no feature value is fit.

The comparison tests the tempting scale-only inference

    feature rate = blanket rate * geometry transport factor.

That inference silently requires transferable material kinetics and a
stationary cycle history.  Yoshie's experiment changes both the substrate
(370 nm poly-Si film versus a Si substrate) and the exposure history
(75 cycles versus 450/675 cycles), so the numerical audit reports rather
than assumes transferability.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

import numpy as np

from petch.experimental_data import (
    load_yoshie_2023_blanket_rates,
    load_yoshie_2023_feature_depths,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "yoshie_2023"
RESULT = ROOT / "results" / "yoshie_2023_blanket_transfer" / "audit.json"


def build_report(data_dir=DATA):
    data_dir = Path(data_dir)
    blankets = load_yoshie_2023_blanket_rates(
        data_dir / "figure4_blanket_poly_si_rates.csv")
    features = load_yoshie_2023_feature_depths(
        data_dir / "figures5_6_feature_depths.csv")
    blanket_by_condition = {
        (int(row.cycle_duration_s), row.timing): row
        for row in blankets
    }
    feature_by_condition = {}
    for row in features:
        feature_by_condition.setdefault(
            (int(row.cycle_duration_s), row.timing), []).append(row)

    conditions = []
    direct_predictions = []
    direct_observations = []
    for condition in sorted(feature_by_condition):
        cycle, timing = condition
        blanket = blanket_by_condition[condition]
        rows = sorted(
            feature_by_condition[condition],
            key=lambda item: item.initial_width_nm)
        rates = np.asarray(
            [item.etch_rate_nm_per_bias_min for item in rows], dtype=float)
        ratios = rates / blanket.blanket_etch_rate_nm_per_bias_min
        ratio_lower = np.asarray([
            max(
                item.etch_rate_nm_per_bias_min
                - item.digitization_rate_uncertainty_nm_per_bias_min,
                0.0,
            ) / (
                blanket.blanket_etch_rate_nm_per_bias_min
                + blanket.digitization_rate_uncertainty_nm_per_bias_min
            )
            for item in rows
        ])
        ratio_upper = np.asarray([
            (
                item.etch_rate_nm_per_bias_min
                + item.digitization_rate_uncertainty_nm_per_bias_min
            ) / (
                blanket.blanket_etch_rate_nm_per_bias_min
                - blanket.digitization_rate_uncertainty_nm_per_bias_min
            )
            for item in rows
        ])
        direct_predictions.extend(
            [blanket.blanket_etch_rate_nm_per_bias_min] * len(rows))
        direct_observations.extend(rates.tolist())
        conditions.append({
            "cycle_duration_s": cycle,
            "timing": timing,
            "blanket_rate_nm_per_bias_min":
                blanket.blanket_etch_rate_nm_per_bias_min,
            "feature_rate_min_nm_per_bias_min": float(np.min(rates)),
            "feature_rate_median_nm_per_bias_min": float(median(rates)),
            "feature_rate_max_nm_per_bias_min": float(np.max(rates)),
            "feature_to_blanket_ratio_min": float(np.min(ratios)),
            "feature_to_blanket_ratio_median": float(median(ratios)),
            "feature_to_blanket_ratio_max": float(np.max(ratios)),
            "ratio_lower_bound_with_digitization_allowance":
                float(np.min(ratio_lower)),
            "ratio_upper_bound_with_digitization_allowance":
                float(np.max(ratio_upper)),
            "feature_points_above_blanket": int(np.sum(ratios > 1.0)),
            "feature_points": len(rows),
            "wide_to_narrow_feature_rate_ratio":
                float(rates[-1] / rates[0]),
        })

    predicted = np.asarray(direct_predictions)
    observed = np.asarray(direct_observations)
    direct_relative_error = np.abs(predicted - observed) / observed
    cycle_ranks = {}
    for cycle in (4, 8):
        selected = [
            item for item in conditions if item["cycle_duration_s"] == cycle
        ]
        cycle_ranks[str(cycle)] = {
            "blanket_descending": [
                item["timing"] for item in sorted(
                    selected,
                    key=lambda item: item[
                        "blanket_rate_nm_per_bias_min"],
                    reverse=True,
                )
            ],
            "feature_median_descending": [
                item["timing"] for item in sorted(
                    selected,
                    key=lambda item: item[
                        "feature_rate_median_nm_per_bias_min"],
                    reverse=True,
                )
            ],
        }

    return {
        "audit_id": "YOSHIE-2023-BLANKET-TRANSFER-R1",
        "operation": (
            "read-only comparison; no surface, transport, reactor, or "
            "geometry parameter was fitted"
        ),
        "source": {
            "citation": (
                "T. Yoshie et al., Applied Surface Science 638 (2023) "
                "157981, DOI 10.1016/j.apsusc.2023.157981"
            ),
            "blanket_material": "370 nm poly-Si film on a Si substrate",
            "feature_material": (
                "bulk Si substrate under an approximately 500 nm SiO2 mask"
            ),
            "blanket_cycle_count": 75,
            "feature_cycle_counts": {"4_s_cycle": 675, "8_s_cycle": 450},
            "wafer_boundary_limit": (
                "species-resolved neutral/ion fluxes and the ion "
                "energy-angle distribution were not reported"
            ),
        },
        "conditions": conditions,
        "rank_order_by_cycle": cycle_ranks,
        "scale_only_null": {
            "definition": (
                "assign each held-out feature the matching measured blanket "
                "rate, before any feature-transport correction"
            ),
            "mean_absolute_percentage_error": float(
                np.mean(direct_relative_error)),
            "rmse_nm_per_bias_min": float(
                np.sqrt(np.mean((predicted - observed) ** 2))),
            "feature_points": len(features),
        },
        "identifiability_verdict": {
            "blanket_boundary_evidence_tier":
                "B_facility_conditioned",
            "direct_scale_transfer_certified": False,
            "why": [
                (
                    "the blanket and feature substrates are not the same "
                    "material state (poly-Si film versus bulk Si)"
                ),
                (
                    "the chamber and surface histories differ by hundreds "
                    "of gas-modulation cycles"
                ),
                (
                    "the 8 s timing-I feature/blanket ratio remains above "
                    "2.5 after the committed digitization allowances"
                ),
                (
                    "the 8 s timing rank changes from II>III>IV>I on the "
                    "blanket to II>I>III>IV in the feature medians"
                ),
                (
                    "one scalar blanket rate cannot identify time-resolved "
                    "Ar+, SFx+, F, and CFx fluxes, their IEADs, and the "
                    "fluorination/polymer state simultaneously"
                ),
            ],
            "permitted_use": (
                "condition a reactor boundary only inside an independently "
                "validated material- and cycle-history-resolved mechanism; "
                "do not multiply feature depth by a blanket-derived scale"
            ),
        },
    }


def canonical_payload(report):
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.check and args.write:
        parser.error("--check and --write are mutually exclusive")
    payload = canonical_payload(build_report())
    if args.check:
        if not RESULT.exists() or RESULT.read_text() != payload:
            raise SystemExit(f"stale Yoshie blanket-transfer audit: {RESULT}")
        print(f"verified {RESULT.relative_to(ROOT)}")
    elif args.write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(payload)
        print(f"wrote {RESULT.relative_to(ROOT)}")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
