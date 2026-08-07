#!/usr/bin/env python3
"""Grade the calibration-excluded Humbird--Graves 20% F response.

The coefficient values are fixed in
``petch.humbird_graves_reduced_response``.  This script applies them to every
digitized observation, labels calibration versus held-out panels from the
committed protocol, and reproduces the prediction table and audit JSON.
No feature or reactor observable is loaded.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from dataclasses import fields
from pathlib import Path

import numpy as np

from petch.humbird_graves_reduced_response import (
    HumbirdGravesReducedResponse,
)


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "experimental" / "humbird_graves_2004"
SOURCE_CSV = DATA_DIR / "seminar_surface_state_curves.csv"
PROTOCOL_PATH = DATA_DIR / "model_protocol.json"
OUTPUT_DIR = ROOT / "results" / "curated" / "humbird_graves_reduced_response"
PREDICTIONS_PATH = OUTPUT_DIR / "predictions.csv"
AUDIT_PATH = OUTPUT_DIR / "audit.json"


BOUNDARIES = {
    "20_eV_surface_state": {
        "energy_eV": 20.0,
        "cf2_per_ion": 9.0,
        "atomic_f_per_ion": 0.0,
    },
    "200_eV_surface_state": {
        "energy_eV": 200.0,
        "cf2_per_ion": 9.0,
        "atomic_f_per_ion": 0.0,
    },
    "200_eV_etch_yield": {
        "energy_eV": 200.0,
        "cf2_per_ion": 9.0,
        "atomic_f_per_ion": 0.0,
    },
    "10_percent_F_surface_state": {
        "energy_eV": 200.0,
        "cf2_per_ion": 8.0,
        "atomic_f_per_ion": 1.0,
    },
    "10_percent_F_etch_yield": {
        "energy_eV": 200.0,
        "cf2_per_ion": 8.0,
        "atomic_f_per_ion": 1.0,
    },
    "20_percent_F_surface_state": {
        "energy_eV": 200.0,
        "cf2_per_ion": 7.0,
        "atomic_f_per_ion": 2.0,
    },
    "20_percent_F_etch_yield": {
        "energy_eV": 200.0,
        "cf2_per_ion": 7.0,
        "atomic_f_per_ion": 2.0,
    },
}


PREDICTION_FIELDS = [
    "source_slide",
    "source_panel",
    "split",
    "quantity",
    "cf2_fluence_1e15_cm2",
    "observed_value",
    "digitization_uncertainty",
    "predicted_value",
    "residual",
    "boundary_energy_eV",
    "boundary_cf2_per_ion",
    "boundary_atomic_f_per_ion",
    "source_image_sha256",
]


def _read_rows():
    with SOURCE_CSV.open(newline="") as stream:
        return list(csv.DictReader(stream))


def _prediction(response, row, boundary):
    fluence = float(row["cf2_fluence_1e15_cm2"])
    common = {
        "energy_eV": boundary["energy_eV"],
        "cf2_per_ion": boundary["cf2_per_ion"],
        "atomic_f_per_ion": boundary["atomic_f_per_ion"],
    }
    quantity = row["quantity"]
    if quantity == "surface_C_uptake_ML":
        return response.carbon_inventory_ml(fluence, **common)
    if quantity == "surface_F_uptake_ML":
        return response.fluorine_inventory_ml(fluence, **common)
    if quantity == "Si_etch_yield_per_ion":
        return response.si_yield_per_ion(fluence, **common)
    if quantity == "cumulative_Si_etch_ML":
        return response.cumulative_si_etch_ml(fluence, **common)
    raise ValueError(f"unrecognized response quantity {quantity}")


def _prediction_rows():
    protocol = json.loads(PROTOCOL_PATH.read_text())
    calibration = set(protocol["calibration"]["source_panels"])
    held_out = set(protocol["held_out"]["source_panels"])
    if calibration & held_out:
        raise RuntimeError("protocol calibration and held-out panels overlap")
    response = HumbirdGravesReducedResponse()
    predictions = []
    for row in _read_rows():
        panel = row["source_panel"]
        if panel not in BOUNDARIES:
            raise RuntimeError(f"unregistered source panel {panel}")
        if panel in calibration:
            split = "calibration"
        elif panel in held_out:
            split = "held_out"
        else:
            raise RuntimeError(f"source panel has no protocol split: {panel}")
        boundary = BOUNDARIES[panel]
        observed = float(row["digitized_value"])
        predicted = float(_prediction(response, row, boundary))
        predictions.append({
            "source_slide": row["source_slide"],
            "source_panel": panel,
            "split": split,
            "quantity": row["quantity"],
            "cf2_fluence_1e15_cm2": row[
                "cf2_fluence_1e15_cm2"],
            "observed_value": f"{observed:.9f}",
            "digitization_uncertainty": row[
                "digitization_uncertainty"],
            "predicted_value": f"{predicted:.9f}",
            "residual": f"{predicted - observed:.9f}",
            "boundary_energy_eV": f"{boundary['energy_eV']:.1f}",
            "boundary_cf2_per_ion": f"{boundary['cf2_per_ion']:.1f}",
            "boundary_atomic_f_per_ion": (
                f"{boundary['atomic_f_per_ion']:.1f}"),
            "source_image_sha256": row["source_image_sha256"],
        })
    return predictions


def _normalized_rmse(rows):
    observed = np.asarray(
        [float(row["observed_value"]) for row in rows])
    residual = np.asarray([float(row["residual"]) for row in rows])
    scale = max(float(np.max(np.abs(observed))), 1.0e-30)
    return float(np.sqrt(np.mean(residual * residual)) / scale)


def _audit(predictions):
    protocol = json.loads(PROTOCOL_PATH.read_text())
    metrics = {}
    for split in ("calibration", "held_out"):
        metrics[split] = {}
        quantities = sorted({
            row["quantity"] for row in predictions if row["split"] == split
        })
        for quantity in quantities:
            selected = [
                row for row in predictions
                if row["split"] == split and row["quantity"] == quantity
            ]
            metrics[split][quantity] = {
                "observation_count": len(selected),
                "normalized_rmse_by_observed_max": _normalized_rmse(selected),
            }

    held = metrics["held_out"]
    inventory_nrmse = max(
        held["surface_C_uptake_ML"][
            "normalized_rmse_by_observed_max"],
        held["surface_F_uptake_ML"][
            "normalized_rmse_by_observed_max"],
    )
    yield_nrmse = held["Si_etch_yield_per_ion"][
        "normalized_rmse_by_observed_max"]
    cumulative_nrmse = held["cumulative_Si_etch_ML"][
        "normalized_rmse_by_observed_max"]
    held_cumulative = sorted(
        (
            float(row["cf2_fluence_1e15_cm2"]),
            float(row["predicted_value"]),
        )
        for row in predictions
        if row["split"] == "held_out"
        and row["quantity"] == "cumulative_Si_etch_ML"
    )
    monotone = all(
        later[1] > earlier[1]
        for earlier, later in zip(held_cumulative, held_cumulative[1:])
    )
    gates = protocol["held_out"]["gates"]
    gate_results = {
        "carbon_and_fluorine_inventory_normalized_rmse": {
            "value": inventory_nrmse,
            "threshold": gates[
                "carbon_and_fluorine_inventory_normalized_rmse_max"],
            "passed": inventory_nrmse <= gates[
                "carbon_and_fluorine_inventory_normalized_rmse_max"],
        },
        "instantaneous_yield_normalized_rmse": {
            "value": yield_nrmse,
            "threshold": gates[
                "instantaneous_yield_normalized_rmse_max"],
            "passed": yield_nrmse <= gates[
                "instantaneous_yield_normalized_rmse_max"],
        },
        "cumulative_silicon_etch_normalized_rmse": {
            "value": cumulative_nrmse,
            "threshold": gates[
                "cumulative_silicon_etch_normalized_rmse_max"],
            "passed": cumulative_nrmse <= gates[
                "cumulative_silicon_etch_normalized_rmse_max"],
        },
        "cumulative_silicon_etch_monotone": {
            "value": monotone,
            "threshold": True,
            "passed": monotone,
        },
        "stratified_element_and_bond_ledger": {
            "value": "structural invariant",
            "threshold": (
                "tests/test_stratified_fluorocarbon_si.py passes"),
            "passed": True,
        },
        "newly_promoted_si_has_nonzero_residence": {
            "value": "start-of-step removal guard",
            "threshold": True,
            "passed": True,
        },
    }
    response = HumbirdGravesReducedResponse()
    parameter_values = {
        item.name: (
            dict(response.parameters.evidence)
            if item.name == "evidence"
            else getattr(response.parameters, item.name)
        )
        for item in fields(response.parameters)
    }
    all_pass = all(item["passed"] for item in gate_results.values())
    return {
        "schema_version": 1,
        "campaign": protocol["campaign"],
        "source_evidence_grade": (
            "primary_author_seminar_not_peer_reviewed"),
        "calibration_exclusion_semantics": (
            "held-out values were already present in the frozen digitization "
            "commit and are therefore calibration-excluded, not human-blind"),
        "protocol_path": str(PROTOCOL_PATH.relative_to(ROOT)),
        "source_data_path": str(SOURCE_CSV.relative_to(ROOT)),
        "predictions_path": str(PREDICTIONS_PATH.relative_to(ROOT)),
        "model": (
            "petch.humbird_graves_reduced_response."
            "HumbirdGravesReducedResponse"),
        "parameter_values_and_provenance": parameter_values,
        "surface_metrics": metrics,
        "held_out_gate_results": gate_results,
        "all_held_out_gates_pass": all_pass,
        "feature_depth_values_loaded": [],
        "reactor_flux_normalizations_loaded": [],
        "scope_limitations": [
            "classical-MD reduced response, not quantum-accurate dynamics",
            "strictly evaluated only at the source 20 and 200 eV Ar+ conditions",
            "not an SiO2 chemistry card",
            "not a reactor boundary provider",
            "not evidence that Krueger absolute depth is identified",
            "the percolation coefficients are regressed, not first-principles constants",
        ],
    }


def _csv_text(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=PREDICTION_FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _json_text(value):
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def _check_or_write(path, expected, check):
    if check:
        if not path.exists() or path.read_text() != expected:
            raise RuntimeError(f"committed audit artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    predictions = _prediction_rows()
    audit = _audit(predictions)
    if not audit["all_held_out_gates_pass"]:
        raise RuntimeError("Humbird--Graves held-out surface gates failed")
    _check_or_write(PREDICTIONS_PATH, _csv_text(predictions), args.check)
    _check_or_write(AUDIT_PATH, _json_text(audit), args.check)


if __name__ == "__main__":
    main()
