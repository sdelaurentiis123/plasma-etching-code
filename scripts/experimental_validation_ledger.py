#!/usr/bin/env python3
"""Build an evidence-strict status ledger for every bundled experimental benchmark."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from petch.experimental_data import (
    load_bosch_wafer_measurements,
    load_bosch_wafer_measurements_89pt,
    load_deboer_2002_figure9_depths,
    load_jeong_2023_etch_depths,
    load_jeong_2023_radical_densities,
    load_jeon_2022_trench_depths,
    load_krueger_2024_evidence,
)
from petch.notching_validation_3d import (
    NOZAWA_1995_NOTCH_CURVES_SHA256,
    load_nozawa_1995_notch_observations,
)
from petch.reactor_boundary import load_krueger_2024_reactor_flux_deck


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental"
RESULTS = ROOT / "results"
READINESS_COLUMNS = (
    "experimental evidence",
    "reactor boundary",
    "surface mechanism",
    "charging/topology",
    "current-operator replay",
    "independent held-out test",
)
READINESS_LABEL = {
    0: "missing",
    1: "partial",
    2: "development-ready",
    3: "claim-ready",
}


def _sha(path):
    path = Path(path)
    return sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _json(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.is_file() else None


def _artifact(path):
    path = Path(path)
    return {
        "path": str(path.relative_to(ROOT)),
        "exists": path.is_file(),
        "sha256": _sha(path),
    }


def _git_revision():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def _readiness(values):
    if len(values) != len(READINESS_COLUMNS):
        raise ValueError("readiness vector has the wrong length")
    return {
        name: {"level": int(value), "label": READINESS_LABEL[int(value)]}
        for name, value in zip(READINESS_COLUMNS, values)
    }


def build_ledger():
    head = _git_revision()
    nozawa_rows = load_nozawa_1995_notch_observations(
        DATA / "nozawa_1995" / "digitized_notch_curves.csv",
        expected_sha256=NOZAWA_1995_NOTCH_CURVES_SHA256,
    )
    jeon_rows = load_jeon_2022_trench_depths(
        DATA / "jeon_2022" / "digitized_trench_depths.csv")
    jeong_depth = load_jeong_2023_etch_depths(
        DATA / "jeong_2023" / "digitized_figure7_depths.csv")
    jeong_radical = load_jeong_2023_radical_densities(
        DATA / "jeong_2023" / "digitized_figure6_radicals.csv")
    deboer_rows = load_deboer_2002_figure9_depths(
        DATA / "deboer_2002" / "digitized_figure9.csv")
    krueger = load_krueger_2024_evidence(DATA / "krueger_2024")
    krueger_deck = load_krueger_2024_reactor_flux_deck(DATA / "krueger_2024")
    bosch_9 = load_bosch_wafer_measurements(
        DATA / "zenodo_17122442" / "Si_Oxide_etch_9_points.csv")
    bosch_89 = load_bosch_wafer_measurements_89pt(
        DATA / "zenodo_17122442" / "Si_Oxide_etch_89_points.csv")

    jeong_audit_path = RESULTS / "jeong_2023_reactor_closure_audit" / "audit.json"
    jeong_audit = _json(jeong_audit_path)
    deboer_path = RESULTS / "deboer_2002_belen_reflection_pilot" / "audit.json"
    deboer = _json(deboer_path)
    hwang_profile_path = (
        RESULTS / "hwang_giapis_1997_fig13_validation_v4_connected_l16" / "audit.json")
    hwang_profile = _json(hwang_profile_path)
    hwang_global_path = (
        RESULTS / "hwang_giapis_1997_fig13_source2d_v10_unclipped_cont1000_from8000"
        / "global_boundary.json")
    hwang_global = _json(hwang_global_path)
    nozawa_charge_path = (
        RESULTS / "nozawa_1995_mixed_compatible_high_precision_audit_50us_20260716"
        / "summary.json")
    nozawa_charge = _json(nozawa_charge_path)

    hwang_profile_score = (
        hwang_profile["profile_runs"]["event_driven_scatter_on"]["score"]
        if hwang_profile else {})
    hwang_global_obs = hwang_global.get("global_observables", {}) if hwang_global else {}
    hwang_stationarity = hwang_global_obs.get("potential_stationarity", {})
    nozawa_finest = nozawa_charge["levels"][-1] if nozawa_charge else {}

    entries = [
        {
            "benchmark": "Hwang--Giapis 1997 Figure 13 notch contour",
            "target": "source-faithful charged-notch profile",
            "status": "current_boundary_open_profile_score_stale",
            "earned_claim": (
                "The historical profile is a strong quantitative development replay only. "
                "Its RMSE cannot be promoted because the contour was not held out, measurement "
                "uncertainty is unreported, and the profile predates the current boundary operator."),
            "readiness": _readiness((3, 3, 2, 2, 1, 0)),
            "metrics": {
                "historical_profile_rmse_um": hwang_profile_score.get("rmse_um"),
                "historical_strict_validation_pass": hwang_profile_score.get(
                    "strict_validation_pass"),
                "historical_profile_git_revision": (
                    hwang_profile.get("git_revision") if hwang_profile else None),
                "historical_profile_is_current_head": (
                    hwang_profile is not None
                    and hwang_profile.get("git_revision") == head),
                "current_global_git_revision": (
                    hwang_global.get("git_revision") if hwang_global else None),
                "current_global_is_current_head": (
                    hwang_global is not None
                    and hwang_global.get("git_revision") == head),
                "current_surface_rms_drift_v": hwang_stationarity.get(
                    "surface_rms_drift_v"),
                "current_surface_max_abs_drift_v": hwang_stationarity.get(
                    "surface_max_abs_drift_v"),
                "current_neighbor_conductor_relative_imbalance": (
                    hwang_global_obs.get("region_current_balance", {})
                    .get("poly_neighbor", {}).get("symmetric_relative_imbalance")),
            },
            "blocking_evidence": [
                "current global charging state has not passed its stationarity refinement gate",
                "current operator has not yet produced a new profile replay",
                "published contour has no reported measurement uncertainty and is development data",
            ],
            "next_bounded_action": (
                "Audit/continue only the current global boundary to stationarity; rerun one profile "
                "afterward, then reserve a different notch experiment for formal validation."),
            "artifacts": [_artifact(hwang_profile_path), _artifact(hwang_global_path)],
        },
        {
            "benchmark": "Nozawa 1995 topology and notch-depth curves",
            "target": "open-area, shared-pad, and individual-pad notch transfer",
            "status": "engine_channels_present_terminal_evidence_and_stationarity_open",
            "earned_claim": (
                "The common 3-D hard-visibility smoke and fixed-geometry charge audits are "
                "operational evidence, not experimental notch validation."),
            "readiness": _readiness((3, 3, 2, 1, 1, 0)),
            "metrics": {
                "digitized_condition_count": len(nozawa_rows),
                "calibration_condition_count": sum(
                    item.split == "calibration" for item in nozawa_rows),
                "held_out_condition_count": sum(
                    item.split == "held_out_transfer" for item in nozawa_rows),
                "latest_q1_resolved_patch_b2_rms": nozawa_finest.get(
                    "ensemble_q1_resolved_patch_b2_rms"),
                "latest_q1_unresolved_face_current_fraction": nozawa_finest.get(
                    "ensemble_q1_unresolved_face_current_fraction"),
                "latest_b1_potential_rate_max_v_s": nozawa_finest.get(
                    "ensemble_b1_potential_rate_max_v_s"),
                "remote_pad_terminal_engine_available": True,
            },
            "blocking_evidence": [
                "shared/individual pad perimeters and conductor groupings must be source-resolved",
                "one legal terminal-current coefficient calibration must be preregistered",
                "terminal-aware conductor balance needs an explicit contract; face B2 is retained",
                "current 3-D charge state has not passed signed stationarity gates",
            ],
            "next_bounded_action": (
                "Digitize/verify pad perimeter and connected-line count from the primary figures; "
                "construct shared-vs-individual conductor IDs and run a tiny manufactured terminal "
                "preflight before any long charging/profile replay."),
            "artifacts": [_artifact(nozawa_charge_path)],
        },
        {
            "benchmark": "Jeong 2023 SiO2 energy/flux/width transfer",
            "target": "one-anchor prediction of 17 held-out depth points",
            "status": "strict_validation_blocked_by_reactor_boundary",
            "earned_claim": (
                "The feature response and failure mechanism are diagnosed. The historical moving "
                "profiles are stale and no held-out predictive transfer is currently earned."),
            "readiness": _readiness((3, 0, 2, 1, 1, 0)),
            "metrics": {
                "depth_point_count": len(jeong_depth),
                "radical_input_count": len(jeong_radical),
                "held_out_count": sum(
                    item.split == "held_out_transfer" for item in jeong_depth),
                "required_endpoint_gain_reduction_fraction": (
                    jeong_audit.get("inverse_missing_response_diagnostic", {})
                    .get("rows") is not None and 0.6162467122993285),
                "collisionless_maximum_flattening_fraction": (
                    jeong_audit.get("bounded_source_backed_tests", {})
                    .get("collisionless_virtual_sheath", {})
                    .get("maximum_density_induced_yield_flattening_fraction")
                    if jeong_audit else None),
            },
            "blocking_evidence": (
                jeong_audit.get("decision", {}).get("blocking_boundary_evidence", [])
                if jeong_audit else ["reactor closure audit missing"]),
            "next_bounded_action": (
                "Supply measured or independently validated species-resolved ion fluxes/IEADs, "
                "complete radical wall fluxes, and hot-neutral distributions through the new "
                "reactor deck. Do not refit the exposed flux sweep."),
            "artifacts": [_artifact(jeong_audit_path)],
        },
        {
            "benchmark": "Jeon 2022 pulsed C4F8/Ar SiO2 ARDE",
            "target": "width and pulse-off transfer across radical regimes",
            "status": "development_only_missing_reactor_boundary_and_absolute_time_context",
            "earned_claim": (
                "The checksummed 54-point evidence and dimensionless targets are ready; existing "
                "runs are explicitly nonpredictive baselines."),
            "readiness": _readiness((3, 0, 2, 0, 1, 0)),
            "metrics": {
                "digitized_depth_count": len(jeon_rows),
                "calibration_count": sum(
                    item.split == "calibration" for item in jeon_rows),
                "held_out_count": sum(
                    item.split == "held_out_transfer" for item in jeon_rows),
            },
            "blocking_evidence": [
                "etch duration is not reported for an absolute-rate replay",
                "measured IEAD/IADF and species-resolved surface fluxes are absent",
                "existing self-bias monoenergy and aggregate-neutral adapters are nonpredictive",
            ],
            "next_bounded_action": (
                "Score only within-condition dimensionless width ratios until a reactor flux deck "
                "or new process metadata supplies absolute boundary conditions and duration."),
            "artifacts": [
                _artifact(RESULTS / "jeon_2022_predictive_validation"
                          / "calibration_baseline_dx001.json"),
            ],
        },
        {
            "benchmark": "de Boer 2002 SF6/O2 silicon ARDE",
            "target": "Figure-9 depth transfer",
            "status": "development_data_exhausted_requires_new_held_out_experiment",
            "earned_claim": (
                "A source-correct Belen mechanism plus certified reflection improves the exposed "
                "development RMSE, but the twelve transfer points have been consumed and cannot "
                "become a fresh held-out validation."),
            "readiness": _readiness((3, 1, 2, 0, 2, 0)),
            "metrics": {
                "digitized_point_count": len(deboer_rows),
                "development_transfer_rmse_um": (
                    deboer.get("development_transfer_score", {}).get("rmse_um")
                    if deboer else None),
                "within_digitization_count": (
                    deboer.get("development_transfer_score", {})
                    .get("within_digitization_count") if deboer else None),
            },
            "blocking_evidence": [
                "measurement uncertainty is unreported",
                "all Figure-9 transfer points were exposed during mechanism development",
                "source omits a complete reactor/IEAD/temperature boundary for a strict replay",
            ],
            "next_bounded_action": (
                "Freeze the current Belen/reflection mechanism and preregister a different "
                "SF6/O2 experiment with independent flux/temperature/profile evidence."),
            "artifacts": [_artifact(deboer_path)],
        },
        {
            "benchmark": "Krüger 2024 C4F6/Ar/O2 SiO2",
            "target": "base profile calibration and O2/power held-out trends",
            "status": "best_reactor_fed_candidate_but_ion_IEAD_and_mechanism_incomplete",
            "earned_claim": (
                "The published HPEM neutral flux vector now enters the common boundary with exact "
                "provenance. It is not a prediction because the ion row is unresolved and the "
                "complete MCFPM reaction network is not implemented."),
            "readiness": _readiness((3, 2, 1, 0, 0, 1)),
            "metrics": {
                "experimental_calibration_metric_count": len(
                    krueger.calibration_metrics),
                "published_reactor_flux_count": len(krueger.boundary_fluxes),
                "held_out_transfer_observation_count": len(
                    krueger.transfer_observations),
                "resolved_neutral_flux_count": sum(
                    item.role == "neutral" for item in krueger_deck.species_fluxes),
                "unresolved_reactor_species": krueger_deck.unresolved_species,
                "reactor_deck_supports_prediction": (
                    krueger_deck.supports_predictive_boundary),
            },
            "blocking_evidence": [
                "Table I reports one aggregate ion flux rather than species-resolved ion flux/IEAD",
                "HPEM fluxes are source-model outputs, not measurements",
                "petch lacks the paper's full complexes/redeposition/crosslinking/mask network",
            ],
            "next_bounded_action": (
                "Obtain the source HPEM species/IEAD deck or reproduce it with a validated reactor "
                "provider; then port only the reaction channels required by the held-out trends."),
            "artifacts": [
                _artifact(DATA / "krueger_2024" / "base_case_boundary_fluxes.csv"),
            ],
        },
        {
            "benchmark": "Zenodo 17122442 Bosch wafer maps",
            "target": "reactor/wafer-scale etch and uniformity transfer",
            "status": "reactor_surrogate_validation_target_not_feature_profile_validation",
            "earned_claim": (
                "Large checksummed wafer datasets are ingestion-ready. They validate a future "
                "reactor/uniformity layer, not microscopic notch/bow/twist physics by themselves."),
            "readiness": _readiness((3, 0, 0, 0, 0, 1)),
            "metrics": {
                "nine_point_measurement_count": len(bosch_9),
                "eighty_nine_point_measurement_count": len(bosch_89),
            },
            "blocking_evidence": [
                "equipment settings have not yet been mapped to species-resolved wafer flux decks",
                "dataset observables are wafer etch/uniformity rather than feature contours",
            ],
            "next_bounded_action": (
                "Use these data to validate a reactor-to-wafer surrogate that emits local flux "
                "decks; couple those decks to petch only after cross-validation on unseen wafers."),
            "artifacts": [
                _artifact(DATA / "zenodo_17122442" / "Si_Oxide_etch_9_points.csv"),
                _artifact(DATA / "zenodo_17122442" / "Si_Oxide_etch_89_points.csv"),
            ],
        },
    ]
    return {
        "schema": "petch-experimental-validation-ledger-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": head,
        "claim_summary": {
            "formal_predictive_validations_passed": 0,
            "quantitative_development_replays": 2,
            "benchmarks_with_checksum_verified_experimental_evidence": len(entries),
            "interpretation": (
                "The engine has substantial verified operators and several quantitative "
                "development matches, but no current, independent, uncertainty-closed held-out "
                "experimental prediction is yet earned."),
        },
        "readiness_scale": READINESS_LABEL,
        "readiness_columns": READINESS_COLUMNS,
        "benchmarks": entries,
        "cross_benchmark_next_order": [
            "No long run until its boundary, topology, mechanism, and score contract are complete.",
            "Finish current Hwang boundary stationarity, then one current-operator profile replay.",
            "Resolve Nozawa pad topology/terminal evidence before shared-vs-individual held-outs.",
            "Use the Krueger flux-deck interface to acquire/ingest species-resolved ion IEAD output.",
            "Reserve new experiments for formal validation; exposed Jeong/de Boer data stay development.",
            "Validate a reactor/wafer surrogate independently on Bosch wafers before coupling it.",
        ],
    }


def _write_markdown(payload, path):
    lines = [
        "# Experimental validation ledger",
        "",
        f"Generated: {payload['created_utc']}",
        f"Git revision: `{payload['git_revision']}`",
        "",
        "## Bottom line",
        "",
        payload["claim_summary"]["interpretation"],
        "",
        "| Benchmark | Current status | Earned claim | Next bounded action |",
        "| --- | --- | --- | --- |",
    ]
    for item in payload["benchmarks"]:
        lines.append(
            f"| {item['benchmark']} | `{item['status']}` | "
            f"{item['earned_claim']} | {item['next_bounded_action']} |")
    lines += ["", "## Blocking evidence", ""]
    for item in payload["benchmarks"]:
        lines += [f"### {item['benchmark']}", ""]
        lines += [f"- {value}" for value in item["blocking_evidence"]]
        lines.append("")
    lines += ["## Cross-benchmark execution order", ""]
    lines += [
        f"{index}. {value}"
        for index, value in enumerate(payload["cross_benchmark_next_order"], start=1)
    ]
    path.write_text("\n".join(lines) + "\n")


def _plot(payload, path):
    names = [item["benchmark"].split(" 20")[0] for item in payload["benchmarks"]]
    matrix = np.asarray([
        [item["readiness"][column]["level"] for column in READINESS_COLUMNS]
        for item in payload["benchmarks"]
    ])
    figure, axis = plt.subplots(figsize=(12.5, 6.8), constrained_layout=True)
    image = axis.imshow(matrix, vmin=0, vmax=3, cmap="RdYlGn", aspect="auto")
    axis.set_xticks(np.arange(len(READINESS_COLUMNS)), READINESS_COLUMNS, rotation=28, ha="right")
    axis.set_yticks(np.arange(len(names)), names)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = int(matrix[row, column])
            axis.text(column, row, str(value), ha="center", va="center",
                      color="black" if value in {1, 2} else "white", fontweight="bold")
    axis.set_title(
        "Experimental validation readiness — an operator existing is not a validation pass",
        fontweight="bold")
    colorbar = figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, ticks=[0, 1, 2, 3])
    colorbar.ax.set_yticklabels(
        ["missing", "partial", "development", "claim-ready"])
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=RESULTS / "experimental_validation_ledger")
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    payload = build_ledger()
    json_path = args.output / "ledger.json"
    markdown_path = args.output / "README.md"
    plot_path = args.output / "readiness.png"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    _write_markdown(payload, markdown_path)
    _plot(payload, plot_path)
    print(json.dumps({
        "ledger": str(json_path),
        "report": str(markdown_path),
        "figure": str(plot_path),
        "formal_predictive_validations_passed": 0,
    }, indent=2))


if __name__ == "__main__":
    main()
