#!/usr/bin/env python3
"""Audit DFT-trained NNP/MD outputs on the Karahashi beam board.

An et al. publish atomistic SiO2 bombardment results generated with a neural
potential trained on DFT energy/force/stress labels plus explicit ZBL
short-range repulsion. Their released Figure 3 data overlap nine of petch's
independently digitized Karahashi F-free CFx+ beam points exactly in species,
energy, target, and normal incidence.

This script transcribes no curve and fits no coefficient. It checksum-binds
the released-output transcription, selects only exact overlap points, and
scores both the NNP/MD result and the frozen Guo closure on the same support.
The comparison is no-yield-fit, but not blind: the An paper was developed and
published with experimental yield comparisons.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from statistics import fmean

from petch.experimental_data import (
    KARAHASHI_2007_FIGURE4_SHA256,
    load_karahashi_2007_reactive_ion_yields,
)


ROOT = Path(__file__).resolve().parents[1]
MODEL_DIRECTORY = ROOT / "data" / "surface_interactions" / "an_2026_nnp"
EXPERIMENT_DIRECTORY = ROOT / "data" / "experimental" / "karahashi_2007"
GUO_AUDIT = (
    ROOT / "results" / "curated" / "guo_karahashi_transfer" / "audit.json")
OUTPUT = (
    ROOT / "results" / "curated" / "an_nnp_karahashi_transfer"
    / "audit.json")

MODEL_DATA_SHA256 = (
    "c897427ffb8055ffe71f042530efc8acda30c05bcdfe52df0fb4542af1c295d4")
AUTHOR_FIGURE_DATA_SHA256 = (
    "7fdb3b72bc55eb47e8e2fbf6b218ddcacdb869136b2978c4a90bebe3a3dafb9a")
AUTHOR_REPOSITORY_COMMIT = (
    "4bcd035090b9f652cda10150c4da4b662143b34e")


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _metrics(rows):
    errors = [
        row["predicted_yield_sio2_per_ion"]
        - row["observed_yield_sio2_per_ion"]
        for row in rows
    ]
    absolute_relative_errors = [
        abs(error) / row["observed_yield_sio2_per_ion"]
        for row, error in zip(rows, errors)
    ]
    return {
        "point_count": len(rows),
        "mean_absolute_error_sio2_per_ion": fmean(map(abs, errors)),
        "root_mean_square_error_sio2_per_ion": math.sqrt(
            fmean(error * error for error in errors)),
        "mean_absolute_relative_error": fmean(absolute_relative_errors),
        "maximum_absolute_relative_error": max(absolute_relative_errors),
        "within_plotted_source_interval_count": sum(
            row["prediction_within_plotted_source_interval"] for row in rows),
        "correct_net_etch_sign_count": sum(
            row["predicted_yield_sio2_per_ion"] > 0.0 for row in rows),
    }


def _load_model_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = []
        for row in csv.DictReader(handle):
            rows.append({
                "substrate": row["substrate"],
                "species": row["species"],
                "energy_eV": float(row["energy_eV"]),
                "yield_si_per_ion": float(row["yield_si_per_ion"]),
                "reported_regime": row["reported_regime"],
                "figure3_plot_status": row["figure3_plot_status"],
                "source_record": row["source_record"],
            })
    return rows


def _comparison_row(observation, prediction, metadata=None):
    error = prediction - observation.yield_sio2_per_ion
    row = {
        "species": observation.species,
        "energy_eV": observation.energy_eV,
        "observed_yield_sio2_per_ion": observation.yield_sio2_per_ion,
        "plotted_lower_yield": observation.plotted_lower_yield,
        "plotted_upper_yield": observation.plotted_upper_yield,
        "digitization_yield_uncertainty": (
            observation.digitization_yield_uncertainty),
        "predicted_yield_sio2_per_ion": prediction,
        "signed_error_sio2_per_ion": error,
        "signed_relative_error": error / observation.yield_sio2_per_ion,
        "prediction_within_plotted_source_interval": bool(
            observation.plotted_lower_yield
            <= prediction
            <= observation.plotted_upper_yield
        ),
        "correct_net_etch_sign": bool(prediction > 0.0),
    }
    if metadata:
        row.update(metadata)
    return row


def build_audit(root: Path | None = None):
    root = root or ROOT
    model_directory = (
        root / "data" / "surface_interactions" / "an_2026_nnp")
    experiment_directory = (
        root / "data" / "experimental" / "karahashi_2007")
    model_path = model_directory / "figure3_reported_etch_yields.csv"
    model_manifest_path = model_directory / "source_manifest.json"
    experiment_path = experiment_directory / "figure4_reactive_ion_yields.csv"
    experiment_manifest_path = (
        experiment_directory / "digitization_manifest.json")
    guo_audit_path = (
        root / "results" / "curated" / "guo_karahashi_transfer"
        / "audit.json")

    model_manifest = json.loads(model_manifest_path.read_text(encoding="utf-8"))
    assert _sha256(model_path) == MODEL_DATA_SHA256
    assert (
        model_manifest["released_figure_data"]["transcription_sha256"]
        == MODEL_DATA_SHA256
    )
    assert (
        model_manifest["released_figure_data"]["source_sha256"]
        == AUTHOR_FIGURE_DATA_SHA256
    )
    assert (
        model_manifest["source"]["author_repository_commit"]
        == AUTHOR_REPOSITORY_COMMIT
    )

    observations = {
        (row.species, row.energy_eV): row
        for row in load_karahashi_2007_reactive_ion_yields(experiment_path)
    }
    model_rows = _load_model_rows(model_path)
    plotted_sio2_rows = [
        row for row in model_rows
        if row["substrate"] == "SiO2"
        and row["figure3_plot_status"] == "plotted"
    ]

    point_rows = []
    for model_row in plotted_sio2_rows:
        key = (model_row["species"], model_row["energy_eV"])
        observation = observations.get(key)
        if observation is None:
            continue
        point_rows.append(_comparison_row(
            observation,
            model_row["yield_si_per_ion"],
            {
                "reported_model_regime": model_row["reported_regime"],
                "model_source_record": model_row["source_record"],
            },
        ))
    point_rows.sort(key=lambda row: (row["species"], row["energy_eV"]))

    guo_audit = json.loads(guo_audit_path.read_text(encoding="utf-8"))
    guo_by_key = {
        (row["species"], row["energy_eV"]): row
        for row in guo_audit["point_board"]
    }
    guo_same_support_rows = []
    for row in point_rows:
        key = (row["species"], row["energy_eV"])
        guo_row = guo_by_key[key]
        observation = observations[key]
        guo_same_support_rows.append(_comparison_row(
            observation,
            guo_row["predicted_yield_sio2_per_ion"],
        ))

    all_metrics = _metrics(point_rows)
    guo_same_support_metrics = _metrics(guo_same_support_rows)
    high_energy_rows = [
        row for row in point_rows if row["energy_eV"] >= 750.0]
    thousand_eV_rows = [
        row for row in point_rows if row["energy_eV"] == 1000.0]
    low_energy_rows = [
        row for row in point_rows if row["energy_eV"] < 750.0]
    per_species_metrics = {
        species: _metrics([
            row for row in point_rows if row["species"] == species])
        for species in ("CF+", "CF2+", "CF3+")
    }

    return {
        "audit_id": "AN-NNP-KARAHASHI-TRANSFER-R1",
        "status": "completed_no_yield_fit_not_blind",
        "question": (
            "Do DFT-trained NNP/ZBL molecular-dynamics outputs transfer "
            "quantitatively to exact-overlap mass-selected CFx+ SiO2 beam "
            "measurements without fitting an etch yield?"
        ),
        "calibration_firewall": {
            "karahashi_yields_used_in_nnp_training_loss": False,
            "petch_fitted_parameters": [],
            "selection_operation": (
                "exact match on substrate, ion species, setpoint energy, and "
                "normal incidence"),
            "interpolation_or_extrapolation": "none",
            "comparison_is_blind": False,
            "why_not_blind": (
                "The authors compare against experimental yields in the "
                "published model-development paper; post-hoc model selection "
                "informed by those comparisons cannot be excluded."
            ),
        },
        "source_receipts": {
            "author_repository_commit": AUTHOR_REPOSITORY_COMMIT,
            "author_figure3_yaml_sha256": AUTHOR_FIGURE_DATA_SHA256,
            "transcribed_model_data_path": str(model_path.relative_to(root)),
            "transcribed_model_data_sha256": _sha256(model_path),
            "model_source_manifest_path": str(
                model_manifest_path.relative_to(root)),
            "model_source_manifest_sha256": _sha256(model_manifest_path),
            "karahashi_data_path": str(experiment_path.relative_to(root)),
            "karahashi_data_sha256": _sha256(experiment_path),
            "loader_pinned_karahashi_data_sha256": (
                KARAHASHI_2007_FIGURE4_SHA256),
            "karahashi_digitization_manifest_path": str(
                experiment_manifest_path.relative_to(root)),
            "karahashi_digitization_manifest_sha256": _sha256(
                experiment_manifest_path),
            "same_support_guo_audit_path": str(
                guo_audit_path.relative_to(root)),
            "same_support_guo_audit_sha256": _sha256(guo_audit_path),
        },
        "experiment_and_exact_overlap": {
            "target": "SiO2",
            "species": ["CF+", "CF2+", "CF3+"],
            "mass_selected": True,
            "neutral_radical_flux": "none",
            "gas_phase_reactions": "none",
            "incidence": "normal",
            "point_count": len(point_rows),
            "energy_support_eV": [250.0, 1000.0],
            "model_points_interpolated": 0,
            "model_points_extrapolated": 0,
        },
        "model_evidence_boundary": {
            "surface_evidence_level": "DFT_trained_NNP_plus_ZBL_MD",
            "experimental_yield_regressed": False,
            "potential_is_exact_first_principles": False,
            "ion_neutralized_before_impact": True,
            "radicals_or_neutral_coflux_included": False,
            "angular_support": "normal_incidence_only",
            "released_sio2_energy_maximum_eV": 1000.0,
            "krueger_iead_maximum_eV": 4821.0,
            "executable_code_or_potential_imported": False,
            "license_file_found_at_pinned_author_commit": False,
        },
        "all_exact_overlap_metrics": all_metrics,
        "per_species_metrics": per_species_metrics,
        "diagnostic_subsets_not_preregistered": {
            "energy_at_least_750_eV": _metrics(high_energy_rows),
            "energy_exactly_1000_eV": _metrics(thousand_eV_rows),
            "energy_below_750_eV": _metrics(low_energy_rows),
        },
        "same_support_baseline": {
            "baseline": "frozen_Guo_Kwon_closure",
            "nn_p_zbl_metrics": all_metrics,
            "guo_metrics": guo_same_support_metrics,
            "mape_reduction_fraction": (
                1.0
                - all_metrics["mean_absolute_relative_error"]
                / guo_same_support_metrics[
                    "mean_absolute_relative_error"]
            ),
            "support_is_identical": True,
        },
        "point_board": point_rows,
        "same_support_guo_point_board": guo_same_support_rows,
        "physics_diagnosis": {
            "atomistic_transfer_validated_over_full_overlap": False,
            "high_energy_prompt_event_kernel_supported": (
                "promising_diagnostic_not_preregistered_validation"),
            "reason": (
                "The DFT-trained NNP/ZBL calculation greatly improves the "
                "same-support error and preserves etch sign, but the full "
                "board retains 67-72% errors at two low-energy points and "
                "was not blind. Atomic spatial resolution alone does not "
                "close the product-escape timescale or evolving film."
            ),
            "dominant_missing_closures": [
                (
                    "finite, depth-dependent diffusion and escape probability "
                    "for volatile reaction products"),
                (
                    "dose-resolved mixed-layer and fluorocarbon-film state "
                    "carried between impacts"),
                (
                    "incident neutral/radical coflux and its competition with "
                    "ion-driven removal"),
                (
                    "species-resolved extension from 1000 eV through the "
                    "measured 2000 eV beam range before any Krueger IEAD use"),
                (
                    "independent angular response at a pinned ion energy"),
            ],
            "krueger_depth_implication": (
                "This model cannot identify Krueger's unpublished ion mixture "
                "or neutral boundary, and its released outputs stop at 1000 "
                "eV while Krueger's IEAD reaches 4821 eV. It is evidence for "
                "the prompt event kernel, not permission to tune or claim the "
                "825 nm depth."
            ),
            "architecture_constraint": (
                "Couple a DFT/NNP-derived prompt collision/chemistry kernel to "
                "petch's atom-balanced finite mixed layer, film evolution, "
                "and a separately evidenced product-escape process. Keep "
                "reactor flux/species inference outside the surface model."
            ),
        },
    }


def canonical_payload(report):
    return json.dumps(report, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = canonical_payload(build_audit())
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(payload, encoding="utf-8")
    elif arguments.check or OUTPUT.exists():
        if not OUTPUT.exists() or OUTPUT.read_text(
                encoding="utf-8") != payload:
            raise RuntimeError("committed An-NNP/Karahashi audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
