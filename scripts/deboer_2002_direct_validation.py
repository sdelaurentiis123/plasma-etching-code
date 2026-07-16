#!/usr/bin/env python3
"""Historical first direct score against de Boer et al. (2002), Figure 9.

This runner deliberately does *not* score the historical 1/.43/.29/.20 Clausing-model curve as
experimental evidence.  It evaluates the common feature engine against directly digitized trench
depths.  One open-feature marker per time series supplies an absolute-rate boundary condition and
one narrow 12.5-minute marker selects ``s_F`` from a frozen candidate set.  All other markers are
width/time transfers that were held out for the first score.

That first score has now occurred and exposed a material model miss.  Consequently every Figure-9
point used here is development data forever; this script is retained for exact replay and must not
be presented as a new held-out validation.  A different experiment is required for validation.

The engine supplies the normalized instantaneous floor-rate function ``f(depth / opening)``.  Each
profile is evolved without an extra fit by integrating

    d depth / dt = v_open * f(depth / opening).

Because most Si-F chemistry values remain calibration closures and the source reports no
measurement uncertainty, the maximum earned claim is a calibrated transfer demonstration—not a
predictive-chemistry validation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from deboer_feature3d import build_deboer_si_mechanism, floor_rate  # noqa: E402
from petch.experimental_data import load_deboer_2002_figure9_depths  # noqa: E402


PROTOCOL_SHA256 = "4df3eb33b83c82f4e54d6f479d508cc7f7ce2bd8eb27cd191f6edf5ec197ca06"
RESOLUTION_ADDENDUM_SHA256 = (
    "71853f2b607ae9288292217abf601dfba8f956e784f55b5df08b297368122caf")
OPERATOR_ADDENDUM_SHA256 = (
    "a28229a85fd7ed4c51522418acb8ef9be1b3c9a649c900971295ecf16a18757d")


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True,
            text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _key(
        *, s_f: float, dx_um: float, seed: int, aspect_ratio: float,
        floor_average: str = "face") -> str:
    average = "" if floor_average == "face" else f"|average={floor_average}"
    return (
        f"sF={s_f:.12g}|dx={dx_um:.12g}|seed={seed}|AR={aspect_ratio:.12g}{average}")


def _rate_curve(
        *, s_f: float, dx_um: float, seed: int, aspect_ratios: list[float], cache_path: Path,
        cache: dict, floor_average: str = "face") -> dict:
    mechanism = build_deboer_si_mechanism(s_F=s_f)
    raw = []
    for aspect_ratio in aspect_ratios:
        key = _key(
            s_f=s_f, dx_um=dx_um, seed=seed, aspect_ratio=aspect_ratio,
            floor_average=floor_average)
        if key not in cache["rates_m_s"]:
            started = time.monotonic()
            cache["rates_m_s"][key] = floor_rate(
                aspect_ratio, mechanism, dx_um=dx_um, seed=seed,
                floor_average=floor_average)
            cache["evaluations"].append({
                "key": key,
                "elapsed_s": time.monotonic() - started,
                "completed_utc": datetime.now(timezone.utc).isoformat(),
            })
            _atomic_json(cache_path, cache)
            print(
                f"computed s_F={s_f:.3f} dx={dx_um:g} seed={seed} "
                f"AR={aspect_ratio:g} average={floor_average}: "
                f"{cache['rates_m_s'][key]:.6e} m/s",
                flush=True)
        raw.append(float(cache["rates_m_s"][key]))
    raw_array = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(raw_array)) or np.any(raw_array <= 0.0):
        raise RuntimeError("common engine returned a non-positive or non-finite floor rate")
    normalized = raw_array / raw_array[0]
    return {
        "aspect_ratio": [float(value) for value in aspect_ratios],
        "raw_rate_m_s": raw_array.tolist(),
        "normalized_rate": normalized.tolist(),
        "floor_average": floor_average,
    }


def _integrated_curve(curve: dict, *, dense_count: int = 20001) -> tuple[np.ndarray, ...]:
    aspect_ratio = np.asarray(curve["aspect_ratio"], dtype=float)
    normalized = np.asarray(curve["normalized_rate"], dtype=float)
    if (aspect_ratio[0] != 0.0 or np.any(np.diff(aspect_ratio) <= 0.0)
            or np.any(normalized <= 0.0)):
        raise ValueError("rate curve must be positive and strictly ordered from aspect ratio zero")
    dense_ar = np.linspace(0.0, float(aspect_ratio[-1]), dense_count)
    # Interpolating log(rate) guarantees a positive physical rate between evaluated nodes.
    log_rate = PchipInterpolator(aspect_ratio, np.log(normalized))(dense_ar)
    dense_rate = np.exp(log_rate)
    transit_integral = cumulative_trapezoid(1.0 / dense_rate, dense_ar, initial=0.0)
    if np.any(np.diff(transit_integral) <= 0.0):
        raise RuntimeError("profile transit integral is not strictly increasing")
    return dense_ar, dense_rate, transit_integral


def predict_depths(rows, curve: dict) -> list[dict]:
    """Predict every source marker using only its series' declared boundary anchor."""
    dense_ar, _, transit_integral = _integrated_curve(curve)
    anchors = {row.series_time_min: row for row in rows if row.split == "boundary_input"}
    open_rates = {}
    for series_time, anchor in anchors.items():
        anchor_ar = anchor.etch_depth_um / anchor.mask_opening_um
        if anchor_ar > dense_ar[-1]:
            raise RuntimeError("boundary anchor lies beyond evaluated aspect-ratio support")
        anchor_integral = float(np.interp(anchor_ar, dense_ar, transit_integral))
        open_rates[series_time] = anchor.mask_opening_um / series_time * anchor_integral

    predictions = []
    for row in rows:
        target_integral = row.series_time_min * open_rates[row.series_time_min] / row.mask_opening_um
        if target_integral > transit_integral[-1]:
            raise RuntimeError(
                f"prediction for {row.series_time_min:g} min/{row.mask_opening_um:g} um "
                "exceeds the preregistered aspect-ratio support")
        predicted_ar = float(np.interp(target_integral, transit_integral, dense_ar))
        predicted_depth = predicted_ar * row.mask_opening_um
        predictions.append({
            "series_time_min": row.series_time_min,
            "mask_opening_um": row.mask_opening_um,
            "observed_depth_um": row.etch_depth_um,
            "predicted_depth_um": predicted_depth,
            "residual_um": predicted_depth - row.etch_depth_um,
            "predicted_aspect_ratio": predicted_ar,
            "digitization_uncertainty_y_um": row.digitization_uncertainty_y_um,
            "measurement_uncertainty_um": row.measurement_uncertainty_um,
            "split": row.split,
            "role": row.role,
            "open_rate_um_min": open_rates[row.series_time_min],
        })
    return predictions


def _score(predictions: list[dict], *, split: str) -> dict:
    selected = [row for row in predictions if row["split"] == split]
    residual = np.asarray([row["residual_um"] for row in selected], dtype=float)
    observed = np.asarray([row["observed_depth_um"] for row in selected], dtype=float)
    covered = [
        abs(row["residual_um"]) <= row["digitization_uncertainty_y_um"] for row in selected]
    return {
        "count": len(selected),
        "rmse_um": float(np.sqrt(np.mean(residual ** 2))),
        "mae_um": float(np.mean(np.abs(residual))),
        "mean_bias_um": float(np.mean(residual)),
        "relative_rmse": float(np.sqrt(np.mean(residual ** 2)) / np.mean(observed)),
        "digitization_only_coverage_fraction": float(np.mean(covered)),
        "digitization_only_covered_count": int(sum(covered)),
    }


def _write_predictions(path: Path, predictions: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(predictions[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(predictions)


def _validate_protocol(protocol: dict, rows) -> None:
    def keys(records):
        return {(float(item["series_time_min"]), float(item["mask_opening_um"])) for item in records}

    actual = {
        split: {(row.series_time_min, row.mask_opening_um) for row in rows if row.split == split}
        for split in ("boundary_input", "calibration", "held_out_transfer")}
    if keys(protocol["boundary_inputs"]) != actual["boundary_input"]:
        raise ValueError("protocol boundary inputs do not match the checksummed evidence split")
    if keys([protocol["calibration"]["row"]]) != actual["calibration"]:
        raise ValueError("protocol calibration point does not match the checksummed evidence split")
    if keys(protocol["held_out_rows"]) != actual["held_out_transfer"]:
        raise ValueError("protocol held-out rows do not match the checksummed evidence split")


def _run_resolution_audit(*, protocol: dict, rows, output: Path, cache_path: Path, cache: dict) -> None:
    addendum_path = (
        ROOT / "data" / "experimental" / "deboer_2002"
        / "direct_validation_resolution_addendum.json")
    if _sha256(addendum_path) != RESOLUTION_ADDENDUM_SHA256:
        raise ValueError("resolution-audit addendum changed after preregistration")
    addendum = json.loads(addendum_path.read_text(encoding="utf-8"))
    parent_path = output / "audit.json"
    if not parent_path.exists():
        raise ValueError("the frozen parent validation must complete before its resolution audit")
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    if (parent.get("protocol_sha256") != PROTOCOL_SHA256
            or parent.get("selected_s_F") != addendum["constraints"]["selected_s_F"]):
        raise ValueError("resolution addendum does not match the completed parent selection")

    config = addendum["fine_grid"]
    fine_curve = _rate_curve(
        s_f=float(addendum["constraints"]["selected_s_F"]),
        dx_um=float(config["dx_um"]), seed=int(config["seed"]),
        aspect_ratios=[float(value) for value in config["aspect_ratio_nodes"]],
        cache_path=cache_path, cache=cache)
    fine_predictions = predict_depths(rows, fine_curve)

    parent_rows = list(csv.DictReader((output / "predictions.csv").open(encoding="utf-8")))
    fine_by_key = {
        (row["series_time_min"], row["mask_opening_um"]): row for row in fine_predictions}
    combined = []
    for row in parent_rows:
        key = (float(row["series_time_min"]), float(row["mask_opening_um"]))
        fine = fine_by_key[key]
        typed = {
            field: (float(value) if field not in {"split", "role", "measurement_uncertainty_um"}
                    and value != "" else value)
            for field, value in row.items()
        }
        typed["fine_grid_predicted_depth_um"] = fine["predicted_depth_um"]
        typed["fine_grid_delta_from_primary_um"] = (
            fine["predicted_depth_um"] - float(row["predicted_depth_um"]))
        typed["three_level_grid_span_um"] = max(
            float(row["predicted_depth_um"]), float(row["grid_check_predicted_depth_um"]),
            fine["predicted_depth_um"]) - min(
                float(row["predicted_depth_um"]), float(row["grid_check_predicted_depth_um"]),
                fine["predicted_depth_um"])
        combined.append(typed)

    held = [row for row in combined if row["split"] == "held_out_transfer"]
    primary_coarse_delta = np.asarray([row["grid_check_delta_um"] for row in held], dtype=float)
    fine_primary_delta = np.asarray([
        row["fine_grid_delta_from_primary_um"] for row in held], dtype=float)
    coarse_rms = float(np.sqrt(np.mean(primary_coarse_delta ** 2)))
    fine_rms = float(np.sqrt(np.mean(fine_primary_delta ** 2)))
    grid_closed = fine_rms <= coarse_rms / 2.0
    fine_score = _score(fine_predictions, split="held_out_transfer")
    result = {
        "campaign": "deboer_2002_figure9_resolution_audit",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "parent_protocol_sha256": PROTOCOL_SHA256,
        "resolution_addendum_sha256": RESOLUTION_ADDENDUM_SHA256,
        "selected_s_F_unchanged": float(addendum["constraints"]["selected_s_F"]),
        "fine_curve": fine_curve,
        "primary_curve": parent["primary_curve"],
        "coarse_curve": parent["refinement_curves"]["grid_check"],
        "held_out_fine_score": fine_score,
        "primary_vs_coarse_prediction_rms_delta_um": coarse_rms,
        "fine_vs_primary_prediction_rms_delta_um": fine_rms,
        "fine_to_coarse_delta_ratio": fine_rms / coarse_rms,
        "grid_closed_by_preregistered_rule": grid_closed,
        "numerical_conclusion": (
            "PASS: fine-level change contracted by at least two"
            if grid_closed else
            "FAIL: selected-model predictions are not yet spatially grid-closed"),
        "experimental_validation_gate": "NOT_EVALUABLE_WITHOUT_MEASUREMENT_UNCERTAINTY",
    }
    _write_predictions(output / "predictions_with_fine_grid.csv", combined)
    _atomic_json(output / "resolution_audit.json", result)
    print(json.dumps(result, indent=2), flush=True)


def _run_operator_audit(*, protocol: dict, rows, output: Path, cache_path: Path, cache: dict) -> None:
    """Audit the corrected, area-authoritative observable without reopening held-out data."""
    addendum_path = (
        ROOT / "data" / "experimental" / "deboer_2002"
        / "direct_validation_operator_addendum.json")
    if _sha256(addendum_path) != OPERATOR_ADDENDUM_SHA256:
        raise ValueError("corrected-operator addendum changed after preregistration")
    addendum = json.loads(addendum_path.read_text(encoding="utf-8"))
    if addendum["parent_protocol_sha256"] != PROTOCOL_SHA256:
        raise ValueError("corrected-operator addendum points to the wrong parent protocol")

    aspect_ratios = [
        float(value) for value in addendum["operator_corrections"]["aspect_ratio_nodes"]]
    selection = addendum["selection_robustness"]
    selection_config = selection["grid"]
    candidate_results = []
    candidate_curves = {}
    for s_f in (float(value) for value in selection["candidate_values"]):
        curve = _rate_curve(
            s_f=s_f, dx_um=float(selection_config["dx_um"]),
            seed=int(selection_config["seed"]), aspect_ratios=aspect_ratios,
            cache_path=cache_path, cache=cache, floor_average="area")
        predictions = predict_depths(rows, curve)
        score = _score(predictions, split="calibration")
        candidate_results.append({"s_F": s_f, "calibration_score": score, "curve": curve})
        candidate_curves[s_f] = curve
        print(
            f"corrected candidate s_F={s_f:.3f}: calibration RMSE "
            f"{score['rmse_um']:.4f} um", flush=True)
    corrected_best = min(
        candidate_results,
        key=lambda item: (item["calibration_score"]["rmse_um"], item["s_F"]))["s_F"]
    declared_s_f = float(addendum["selected_model_audit"]["s_F"])
    selection_robust = corrected_best == declared_s_f
    if not selection_robust:
        result = {
            "campaign": "deboer_2002_figure9_corrected_operator_audit",
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "operator_addendum_sha256": OPERATOR_ADDENDUM_SHA256,
            "candidate_results": candidate_results,
            "selection_robust": False,
            "corrected_best_s_F": corrected_best,
            "numerical_conclusion": (
                "STOP: parent calibration selection is not robust to the corrected observable"),
        }
        _atomic_json(output / "operator_audit.json", result)
        print(json.dumps(result, indent=2), flush=True)
        return

    curves = {"coarse": candidate_curves[declared_s_f]}
    selected = addendum["selected_model_audit"]
    for label in ("medium", "fine", "independent_scramble"):
        config = selected[label]
        curves[label] = _rate_curve(
            s_f=declared_s_f, dx_um=float(config["dx_um"]), seed=int(config["seed"]),
            aspect_ratios=aspect_ratios, cache_path=cache_path, cache=cache,
            floor_average="area")
    predictions = {label: predict_depths(rows, curve) for label, curve in curves.items()}
    by_level = {
        label: {(row["series_time_min"], row["mask_opening_um"]): row for row in values}
        for label, values in predictions.items()}
    combined = []
    for fine in predictions["fine"]:
        key = (fine["series_time_min"], fine["mask_opening_um"])
        row = dict(fine)
        row["predicted_depth_um"] = by_level["fine"][key]["predicted_depth_um"]
        for label in ("coarse", "medium", "independent_scramble"):
            row[f"{label}_predicted_depth_um"] = by_level[label][key]["predicted_depth_um"]
        row["medium_minus_coarse_um"] = (
            row["medium_predicted_depth_um"] - row["coarse_predicted_depth_um"])
        row["fine_minus_medium_um"] = (
            row["predicted_depth_um"] - row["medium_predicted_depth_um"])
        row["independent_minus_medium_um"] = (
            row["independent_scramble_predicted_depth_um"]
            - row["medium_predicted_depth_um"])
        level_predictions = [
            row["predicted_depth_um"], row["coarse_predicted_depth_um"],
            row["medium_predicted_depth_um"],
            row["independent_scramble_predicted_depth_um"]]
        row["numerical_envelope_um"] = max(level_predictions) - min(level_predictions)
        combined.append(row)

    held = [row for row in combined if row["split"] == "held_out_transfer"]
    def rms(field):
        values = np.asarray([row[field] for row in held], dtype=float)
        return float(np.sqrt(np.mean(values ** 2)))

    medium_coarse_rms = rms("medium_minus_coarse_um")
    fine_medium_rms = rms("fine_minus_medium_um")
    scramble_rms = rms("independent_minus_medium_um")
    fine_medium_max = float(max(abs(row["fine_minus_medium_um"]) for row in held))
    contraction = fine_medium_rms < medium_coarse_rms
    invariant_to_digitization = fine_medium_max <= float(
        addendum["decision_rule"]["observable_invariance_required_for_experimental_claim_um"])
    result = {
        "campaign": "deboer_2002_figure9_corrected_operator_audit",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "parent_protocol_sha256": PROTOCOL_SHA256,
        "operator_addendum_sha256": OPERATOR_ADDENDUM_SHA256,
        "operator": {
            "floor_average": "area",
            "fixed_physical_geometry": True,
            "aspect_ratio_nodes": aspect_ratios,
        },
        "candidate_results": candidate_results,
        "selection_robust": True,
        "selected_s_F": declared_s_f,
        "curves": curves,
        "scores": {
            label: {
                "calibration": _score(values, split="calibration"),
                "held_out": _score(values, split="held_out_transfer"),
            } for label, values in predictions.items()
        },
        "numerical_refinement": {
            "medium_minus_coarse_held_out_rms_um": medium_coarse_rms,
            "fine_minus_medium_held_out_rms_um": fine_medium_rms,
            "fine_minus_medium_held_out_max_abs_um": fine_medium_max,
            "independent_scramble_minus_medium_held_out_rms_um": scramble_rms,
            "grid_change_contracts": contraction,
            "observable_invariant_within_digitization_bound": invariant_to_digitization,
        },
        "numerical_conclusion": (
            "PASS: corrected prediction contracts and closes inside the digitization bound"
            if contraction and invariant_to_digitization else
            "OPEN: corrected prediction has not closed inside the digitization bound"),
        "experimental_validation_gate": "NOT_EVALUABLE_WITHOUT_MEASUREMENT_UNCERTAINTY",
        "earned_claim": (
            "calibrated transfer demonstration only; predictive chemistry, absolute rate, "
            "charging, and formal experimental validation are not claimed"),
    }
    _write_predictions(output / "operator_audit_predictions.csv", combined)
    _atomic_json(output / "operator_audit.json", result)
    print(json.dumps({
        "selection_robust": result["selection_robust"],
        "selected_s_F": result["selected_s_F"],
        "fine_held_out_score": result["scores"]["fine"]["held_out"],
        "numerical_refinement": result["numerical_refinement"],
        "numerical_conclusion": result["numerical_conclusion"],
        "experimental_validation_gate": result["experimental_validation_gate"],
    }, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results" / "deboer_2002_direct_validation")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--resolution-audit", action="store_true")
    parser.add_argument("--operator-audit", action="store_true")
    args = parser.parse_args()

    protocol_path = ROOT / "data" / "experimental" / "deboer_2002" / "direct_validation_protocol.json"
    evidence_path = ROOT / "data" / "experimental" / "deboer_2002" / "digitized_figure9.csv"
    if _sha256(protocol_path) != PROTOCOL_SHA256:
        raise ValueError("direct-validation protocol changed after preregistration")
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    rows = load_deboer_2002_figure9_depths(evidence_path)
    _validate_protocol(protocol, rows)
    if args.validate_only:
        print("protocol, source checksum, pixel replay, and frozen split: PASS")
        return

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache_path = output / "rate_cache.json"
    if cache_path.exists():
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
        if cache.get("protocol_sha256") != PROTOCOL_SHA256:
            raise ValueError("rate cache belongs to a different protocol")
    else:
        cache = {"protocol_sha256": PROTOCOL_SHA256, "rates_m_s": {}, "evaluations": []}

    if args.resolution_audit:
        _run_resolution_audit(
            protocol=protocol, rows=rows, output=output, cache_path=cache_path, cache=cache)
        return
    if args.operator_audit:
        _run_operator_audit(
            protocol=protocol, rows=rows, output=output, cache_path=cache_path, cache=cache)
        return

    engine = protocol["engine"]
    aspect_ratios = [float(value) for value in engine["aspect_ratio_nodes"]]
    candidate_values = [float(value) for value in protocol["calibration"]["candidate_values"]]
    primary_config = engine["primary"]
    candidate_results = []
    for s_f in candidate_values:
        curve = _rate_curve(
            s_f=s_f, dx_um=float(primary_config["dx_um"]), seed=int(primary_config["seed"]),
            aspect_ratios=aspect_ratios, cache_path=cache_path, cache=cache)
        predictions = predict_depths(rows, curve)
        calibration_score = _score(predictions, split="calibration")
        candidate_results.append({
            "s_F": s_f,
            "curve": curve,
            "calibration_score": calibration_score,
            "calibration_predictions": [
                row for row in predictions if row["split"] == "calibration"],
        })
        print(
            f"candidate s_F={s_f:.3f}: calibration RMSE "
            f"{calibration_score['rmse_um']:.4f} um", flush=True)

    best = min(candidate_results, key=lambda item: (item["calibration_score"]["rmse_um"], item["s_F"]))
    best_s_f = float(best["s_F"])
    primary_predictions = predict_depths(rows, best["curve"])
    print(f"selected s_F={best_s_f:.3f} using calibration only", flush=True)

    refinements = {}
    refinement_predictions = {}
    for label in ("grid_check", "sampling_check"):
        config = engine[label]
        curve = _rate_curve(
            s_f=best_s_f, dx_um=float(config["dx_um"]), seed=int(config["seed"]),
            aspect_ratios=aspect_ratios, cache_path=cache_path, cache=cache)
        predictions = predict_depths(rows, curve)
        refinements[label] = curve
        refinement_predictions[label] = predictions

    primary_by_key = {
        (row["series_time_min"], row["mask_opening_um"]): row for row in primary_predictions}
    for label, predictions in refinement_predictions.items():
        for row in predictions:
            key = (row["series_time_min"], row["mask_opening_um"])
            primary_by_key[key][f"{label}_predicted_depth_um"] = row["predicted_depth_um"]
            primary_by_key[key][f"{label}_delta_um"] = (
                row["predicted_depth_um"] - primary_by_key[key]["predicted_depth_um"])
    for row in primary_predictions:
        row["numerical_envelope_um"] = max(
            abs(row["grid_check_delta_um"]), abs(row["sampling_check_delta_um"]))

    held_out_score = _score(primary_predictions, split="held_out_transfer")
    held_out_rows = [row for row in primary_predictions if row["split"] == "held_out_transfer"]
    numerical = {
        label: {
            "held_out_prediction_rmse_delta_um": float(np.sqrt(np.mean([
                row[f"{label}_delta_um"] ** 2 for row in held_out_rows]))),
            "held_out_prediction_max_abs_delta_um": float(max(
                abs(row[f"{label}_delta_um"]) for row in held_out_rows)),
        } for label in ("grid_check", "sampling_check")
    }
    manifest = {
        "campaign": "deboer_2002_figure9_direct_validation",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_id": protocol["protocol_id"],
        "protocol_sha256": PROTOCOL_SHA256,
        "evidence_csv_sha256": protocol["source"]["csv_sha256"],
        "git_revision": _git_revision(),
        "worktree_dirty": bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=ROOT, capture_output=True,
            text=True).stdout.strip()),
        "python": sys.version,
        "platform": platform.platform(),
        "claim": protocol["claim_boundary"],
        "selected_s_F": best_s_f,
        "candidate_results": candidate_results,
        "primary_curve": best["curve"],
        "refinement_curves": refinements,
        "calibration_score": _score(primary_predictions, split="calibration"),
        "held_out_score": held_out_score,
        "numerical_refinement": numerical,
        "measurement_uncertainty_reported": False,
        "experimental_validation_gate": "NOT_EVALUABLE_WITHOUT_MEASUREMENT_UNCERTAINTY",
        "interpretation": (
            "The score tests calibrated ARDE-shape transfer through the common engine. "
            "It does not validate predictive Si-F chemistry, absolute etch rate, or charging."),
    }
    _write_predictions(output / "predictions.csv", primary_predictions)
    _atomic_json(output / "audit.json", manifest)
    print(json.dumps({
        "selected_s_F": best_s_f,
        "held_out_score": held_out_score,
        "numerical_refinement": numerical,
        "experimental_validation_gate": manifest["experimental_validation_gate"],
        "output": str(output),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
