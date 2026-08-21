#!/usr/bin/env python3
"""Write the Bosch v8 chronological prediction without opening heldout outcomes.

This process-only writer is deliberately separate from the future scorer.  It
loads the frozen process traces, the calibration-only measurement asset for the
official coordinate map and empirical baselines, and the already certified v7
wall/v8 ion-boundary parameters.  It never opens the mixed 89-point outcome
CSV that contains the heldout answers.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.bosch_process_data import (  # noqa: E402
    load_bosch_process_traces,
    load_bosch_wafer_measurements_89pt,
)
from scripts.audit_bosch_dynamic_wall_calibration import (  # noqa: E402
    CALIBRATION,
    CYLINDRICAL_IMPLEMENTATION,
    DICTIONARY,
    PROCESS,
    SURFACE_IMPLEMENTATION,
)
from scripts.audit_bosch_recipe_path_memory_calibration import (  # noqa: E402
    _recipe_path_steps,
)
from scripts.audit_bosch_recipe_path_spatial_residual import (  # noqa: E402
    SUMMARY,
    _load_features,
)
from scripts.audit_bosch_wafer_boundary_map_calibration import (  # noqa: E402
    CAPACITY_OUTPUT,
    EXACT_REPLAY_OUTPUT,
    INTERPOLATION_OUTPUT,
    OUTPUT as V8_FIT,
    PREREGISTRATION as V8_PREREGISTRATION,
    RESPONSE_TABLE,
    V7_FIT,
    _exact_map_predictions,
    _law_from_coefficients,
    _v7_coefficients,
)


V1_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_reactor_depth_holdout_v1" / "preregistration.json"
)
V2_PREREGISTRATION = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_cylindrical_depth_extension_v2" / "preregistration.json"
)
DIRECTORY = (
    ROOT / "results" / "curated"
    / "zenodo_bosch_wafer_boundary_map_depth_extension_v8"
)
OUTPUT = DIRECTORY / "heldout_prediction.json"
SEAL = DIRECTORY / "heldout_prediction_seal.json"
SCRIPT = Path(__file__).resolve()


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


@dataclass(frozen=True)
class _HeldoutQuery:
    experiment_key: str
    x_um: np.ndarray
    y_um: np.ndarray


def _certified_inputs():
    v1 = _load_json(V1_PREREGISTRATION)
    v2 = _load_json(V2_PREREGISTRATION)
    v8 = _load_json(V8_PREREGISTRATION)
    capacity = _load_json(CAPACITY_OUTPUT)
    fit = _load_json(V8_FIT)
    interpolation = _load_json(INTERPOLATION_OUTPUT)
    replay = _load_json(EXACT_REPLAY_OUTPUT)
    expected_replay_hashes = {
        "preregistration": _hash(V8_PREREGISTRATION),
        "response_table": _hash(RESPONSE_TABLE),
        "calibration_fit": _hash(V8_FIT),
        "interpolation_validation": _hash(INTERPOLATION_OUTPUT),
        "v7_fit": _hash(V7_FIT),
    }
    if (
        v8["target_firewall"]["heldout_outcomes_read_at_freeze"] is not False
        or capacity["heldout_outcomes_read"] is not False
        or fit["heldout_outcomes_read"] is not False
        or fit["heldout_prediction_written"] is not False
        or interpolation["heldout_outcomes_read"] is not False
        or replay["heldout_outcomes_read"] is not False
        or replay["heldout_prediction_written"] is not False
        or replay["all_replay_and_refinement_gates_pass"] is not True
        or replay["input_hashes"] != expected_replay_hashes
    ):
        raise RuntimeError("Bosch v8 did not earn a heldout prediction seal")
    heldout_keys = tuple(v1["split_rule"]["heldout_experiment_keys"])
    calibration_keys = tuple(v1["split_rule"]["calibration_experiment_keys"])
    lot_type_by_date = v2["conditioning_state"]["declared_lot_types"]
    if (
        len(heldout_keys) != 20
        or len(set(heldout_keys)) != 20
        or set(heldout_keys) & set(calibration_keys)
        or not set(key[:10] for key in heldout_keys).issubset(lot_type_by_date)
    ):
        raise RuntimeError("Bosch chronological split is stale")
    return (
        v1, v2, capacity, fit, interpolation, replay,
        heldout_keys, calibration_keys, lot_type_by_date,
    )


def build_prediction(*, workers=8):
    (
        v1, _v2, capacity, fit, interpolation, replay,
        heldout_keys, calibration_keys, lot_type_by_date,
    ) = _certified_inputs()

    # CALIBRATION is a brokered calibration-only asset.  The mixed source CSV
    # is intentionally absent from every input and hash below.
    calibration = load_bosch_wafer_measurements_89pt(
        CALIBRATION, allowed_experiment_keys=calibration_keys)
    first = calibration[0]
    if any(
        not np.array_equal(item.x_um, first.x_um)
        or not np.array_equal(item.y_um, first.y_um)
        for item in calibration
    ):
        raise RuntimeError("Bosch calibration wafers do not share one map")
    calibration_silicon = np.stack([
        item.silicon_depth_um for item in calibration])
    calibration_oxide = np.stack([
        item.oxide_loss_um for item in calibration])
    calibration_mean_map = np.mean(calibration_silicon, axis=0)

    all_traces = load_bosch_process_traces(PROCESS, DICTIONARY)
    trace_by_key = {trace.experiment_key: trace for trace in all_traces}
    if not set(heldout_keys).issubset(trace_by_key):
        raise RuntimeError("Bosch heldout process traces are incomplete")
    heldout_traces = tuple(sorted(
        (trace_by_key[key] for key in heldout_keys),
        key=lambda item: (item.process_date, item.wafer_number),
    ))
    ordered_keys = tuple(trace.experiment_key for trace in heldout_traces)
    if set(ordered_keys) != set(heldout_keys):
        raise RuntimeError("Bosch heldout process ordering lost a key")

    v7_payload = _load_json(V7_FIT)
    steps, static_law, path_law = _recipe_path_steps(
        _v7_coefficients(v7_payload), heldout_traces, lot_type_by_date)
    wall_multipliers = np.asarray([
        steps[key].combined_wall_loss_multiplier for key in ordered_keys
    ])
    selected = fit["selected_candidate"]
    ion_law = _law_from_coefficients(
        int(selected["static_maximum_order"]),
        int(selected["dynamic_maximum_order"]),
        selected["coefficients"],
    )
    voltages = _load_features(ordered_keys)[
        "c4f8_platen_peak_to_peak_rms"]
    queries = tuple(
        _HeldoutQuery(
            experiment_key=key,
            x_um=np.asarray(first.x_um, dtype=float),
            y_um=np.asarray(first.y_um, dtype=float),
        )
        for key in ordered_keys
    )
    predictions = _exact_map_predictions(
        queries,
        heldout_traces,
        wall_multipliers,
        voltages,
        (ion_law,) * len(queries),
        workers=workers,
    )

    rows = []
    for index, (key, prediction) in enumerate(zip(ordered_keys, predictions)):
        silicon_um = np.asarray(prediction.silicon_depth_m) * 1.0e6
        oxide_um = np.asarray(prediction.oxide_loss_m) * 1.0e6
        silicon_mean = float(np.mean(silicon_um))
        oxide_mean = float(np.mean(oxide_um))
        rows.append({
            "experiment_key": key,
            "process_date": heldout_traces[index].process_date,
            "wafer_number": int(heldout_traces[index].wafer_number),
            "c4f8_platen_peak_to_peak_rms_V": float(voltages[index]),
            "combined_wall_loss_multiplier": float(wall_multipliers[index]),
            "predicted_silicon_depth_um": silicon_um.tolist(),
            "predicted_silicon_mean_depth_um": silicon_mean,
            "predicted_oxide_loss_um": oxide_um.tolist(),
            "predicted_oxide_mean_loss_um": oxide_mean,
            "predicted_selectivity": silicon_mean / oxide_mean,
        })

    return {
        "schema": "petch-zenodo-bosch-v8-heldout-prediction-v1",
        "status": "hash-sealable process-only chronological prediction",
        "heldout_process_record_count": len(rows),
        "heldout_experiment_keys": list(ordered_keys),
        "official_map": {
            "point_count": 89,
            "x_um": first.x_um.tolist(),
            "y_um": first.y_um.tolist(),
        },
        "calibration_only_baselines": {
            "calibration_wafer_count": len(calibration),
            "global_mean_silicon_depth_um": float(np.mean(
                np.mean(calibration_silicon, axis=1))),
            "mean_silicon_depth_map_um": calibration_mean_map.tolist(),
            "mean_silicon_depth_map_mean_um": float(np.mean(
                calibration_mean_map)),
            "global_mean_oxide_loss_um": float(np.mean(
                np.mean(calibration_oxide, axis=1))),
        },
        "frozen_model": {
            "v7_wall_parameters": v7_payload["parameters"],
            "v7_static_law": static_law.manifest(),
            "v7_recipe_path_law": path_law.manifest(),
            "v8_selected_candidate": selected,
            "v8_exact_replay_hash": _hash(EXACT_REPLAY_OUTPUT),
            "v8_replay_all_gates_pass": replay[
                "all_replay_and_refinement_gates_pass"],
            "v8_interpolation_pass": interpolation["passes"],
            "v8_candidate_count": len(capacity["candidates"]),
        },
        "predictions": rows,
        "target_firewall": {
            "mixed_outcome_csv_opened": False,
            "heldout_numeric_outcomes_read": False,
            "heldout_outcomes_read": False,
            "heldout_prediction_written": True,
            "target_depth_used": False,
            "target_radial_map_used": False,
            "eligible_for_separate_outcome_score_after_hash_commit": True,
            "preexposure_blind": v1["target_firewall"]["preexposure_blind"],
        },
        "input_hashes": {
            "v1_preregistration": _hash(V1_PREREGISTRATION),
            "v2_preregistration": _hash(V2_PREREGISTRATION),
            "v8_preregistration": _hash(V8_PREREGISTRATION),
            "v7_fit": _hash(V7_FIT),
            "v8_full_capacity": _hash(CAPACITY_OUTPUT),
            "v8_calibration_fit": _hash(V8_FIT),
            "v8_interpolation": _hash(INTERPOLATION_OUTPUT),
            "v8_exact_replay": _hash(EXACT_REPLAY_OUTPUT),
            "v8_response_table": _hash(RESPONSE_TABLE),
            "calibration_only_measurements": _hash(CALIBRATION),
            "process_data": _hash(PROCESS),
            "process_dictionary": _hash(DICTIONARY),
            "process_summary": _hash(SUMMARY),
        },
    }


def build_seal(prediction_bytes):
    prediction_hash = sha256(prediction_bytes).hexdigest()
    return {
        "schema": "petch-zenodo-bosch-v8-heldout-prediction-seal-v1",
        "prediction_file": str(OUTPUT.relative_to(ROOT)),
        "prediction_sha256": prediction_hash,
        "heldout_outcomes_read": False,
        "heldout_prediction_written": True,
        "eligible_for_separate_outcome_score_after_commit": True,
        "code_hashes": {
            "prediction_writer": _hash(SCRIPT),
            "v8_calibration_audit": _hash(
                ROOT / "scripts"
                / "audit_bosch_wafer_boundary_map_calibration.py"),
            "cylindrical_reactor": _hash(CYLINDRICAL_IMPLEMENTATION),
            "surface_recurrence": _hash(SURFACE_IMPLEMENTATION),
        },
        "artifact_hashes": {
            "v8_preregistration": _hash(V8_PREREGISTRATION),
            "v8_full_capacity": _hash(CAPACITY_OUTPUT),
            "v8_calibration_fit": _hash(V8_FIT),
            "v8_interpolation": _hash(INTERPOLATION_OUTPUT),
            "v8_exact_replay": _hash(EXACT_REPLAY_OUTPUT),
            "v8_response_table": _hash(RESPONSE_TABLE),
        },
        "scoring_authority": {
            "mixed_outcome_csv_may_be_opened_only_after_this_seal_and_prediction_are_committed_and_pushed": True,
            "model_or_prediction_changes_after_outcome_reveal_forbidden": True,
        },
    }


def build(*, workers=8):
    prediction = build_prediction(workers=workers)
    prediction_bytes = _render(prediction).encode("utf-8")
    seal = build_seal(prediction_bytes)
    return prediction_bytes, _render(seal).encode("utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(argv)
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    prediction, seal = build(workers=args.workers)
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != prediction:
            raise SystemExit("Bosch v8 heldout prediction is stale")
        if not SEAL.exists() or SEAL.read_bytes() != seal:
            raise SystemExit("Bosch v8 heldout prediction seal is stale")
        print("Bosch v8 heldout prediction and seal are current")
        return 0
    OUTPUT.write_bytes(prediction)
    SEAL.write_bytes(seal)
    print(OUTPUT.relative_to(ROOT))
    print(SEAL.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
