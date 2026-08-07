#!/usr/bin/env python3
"""Build the claim-scoped cross-chemistry absolute-depth evidence ledger."""

from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TINACBA = (
    ROOT / "results" / "curated" / "tinacba_2021_sf5_depth" / "audit.json")
VELLA = ROOT / "results" / "curated" / "vella_hao_ale_depth" / "audit.json"
KRUEGER = (
    ROOT / "results" / "curated" / "depth_identifiability" / "audit.json")
YOSHIE_VISION = (
    ROOT / "results" / "curated" / "yoshie_2023_digitization_visual_audit"
    / "audit.json")
YOSHIE_FEATURES = (
    ROOT / "data" / "experimental" / "yoshie_2023"
    / "figures5_6_feature_depths.csv")
OUTPUT = (
    ROOT / "results" / "curated" / "cross_chemistry_depth_evidence"
    / "audit.json")


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build_audit() -> dict:
    tinacba = _load(TINACBA)
    vella = _load(VELLA)
    krueger = _load(KRUEGER)
    yoshie_vision = _load(YOSHIE_VISION)
    with YOSHIE_FEATURES.open(newline="", encoding="utf-8") as stream:
        yoshie_rows = list(csv.DictReader(stream))
    if not yoshie_rows or any(
        row["split"] != "held_out_transfer" for row in yoshie_rows
    ):
        raise ValueError("Yoshie feature board is not entirely held out")
    if yoshie_vision["verdict"] != "passed":
        raise ValueError("Yoshie visual digitization audit has not passed")

    tinacba_comparison = tinacba["comparison"]
    vella_board = vella["absolute_depth_board"]
    krueger_inputs = krueger["inputs"]
    krueger_verdict = krueger["verdict"]
    krueger_error = (
        krueger_inputs["simulated_depth_nm"]
        / krueger_inputs["target_depth_nm"]
        - 1.0
    )
    return {
        "schema": "petch.cross-chemistry-depth-evidence.v1",
        "audit_id": "CROSS-CHEMISTRY-ABSOLUTE-DEPTH-2026-08-06-R1",
        "claim_rule": (
            "surface depth-per-dose, development replay, evaluated feature "
            "depth, and value-blind held-out feature prediction are distinct "
            "evidence classes and may not be relabeled as one another"),
        "source_receipts": {
            "tinacba_2021_sf5_depth_sha256": _sha256(TINACBA),
            "vella_hao_ale_depth_sha256": _sha256(VELLA),
            "krueger_depth_identifiability_sha256": _sha256(KRUEGER),
            "yoshie_visual_audit_sha256": _sha256(YOSHIE_VISION),
            "yoshie_feature_table_sha256": _sha256(YOSHIE_FEATURES),
        },
        "boards": [
            {
                "board": "Tinacba 2021 SF5+ beam",
                "chemistry_and_material": "mass-selected SF5+ on Si and SiO2",
                "evidence_class": "retrospective_surface_depth_per_dose",
                "absolute_depth_or_depth_per_dose": True,
                "point_count": tinacba_comparison["point_count"],
                "mean_absolute_relative_error": tinacba_comparison[
                    "mean_absolute_relative_depth_error"],
                "maximum_absolute_relative_error": tinacba_comparison[
                    "maximum_absolute_relative_depth_error"],
                "fit_to_compared_depth_or_yield": tinacba["provider"][
                    "beam_depth_or_yield_fit_used"],
                "boundary_strength": (
                    "mass-selected ion; measured energy and sample-position "
                    "dose; normal incidence"),
                "uses_common_material_router": True,
                "feature_profile_test": False,
                "value_blind_held_out": False,
                "formal_pass_label": False,
                "scope": tinacba["claim"],
            },
            {
                "board": "Vella-Hao / Kounis-Melas Si ALE",
                "chemistry_and_material": "Cl2 adsorption plus Ar+ on Si",
                "evidence_class": "retrospective_no_depth_fit_cross_source",
                "absolute_depth_or_depth_per_dose": True,
                "point_count": len(vella_board["comparison_energy_eV"]),
                "mean_absolute_error_nm": vella_board["mae_nm"],
                "maximum_absolute_relative_error": vella_board[
                    "maximum_absolute_relative_error"],
                "fit_to_compared_depth_or_yield": vella_board[
                    "fit_to_depth_used"],
                "nominal_predeclared_gate": {
                    "limit": vella_board[
                        "nominal_gate_maximum_relative_error"],
                    "passed": vella_board["nominal_gate_passed"],
                },
                "boundary_strength": (
                    "measured positive-ion fluence; inferred mean energy; "
                    "no measured species-resolved IEAD"),
                "uses_atom_conservative_event_ledgers": True,
                "feature_profile_test": False,
                "value_blind_held_out": False,
                "formal_pass_label": True,
                "scope": vella["claim"],
            },
            {
                "board": "Krueger 2024 Ar/C4F6/O2 trench",
                "chemistry_and_material": "Ar/C4F6/O2 on masked SiO2",
                "evidence_class": "evaluated_feature_depth_with_missing_boundary",
                "absolute_depth_or_depth_per_dose": True,
                "point_count": 1,
                "target_depth_nm": krueger_inputs["target_depth_nm"],
                "predicted_depth_nm": krueger_inputs["simulated_depth_nm"],
                "signed_relative_error": krueger_error,
                "fit_to_compared_depth_or_yield": False,
                "published_boundary_identifies_depth": krueger_verdict[
                    "published_inputs_identify_absolute_depth"],
                "feature_profile_test": True,
                "value_blind_held_out": False,
                "formal_pass_label": False,
                "scope": (
                    "published aggregate boundary; base feature was the "
                    "source mechanism's optimization target; absolute depth "
                    "misses and remains underidentified"),
            },
            {
                "board": "Yoshie 2023 cyclic SF6/C4F8 trench",
                "chemistry_and_material": (
                    "alternating C4F8/SF6 with Ar bias on Si"),
                "evidence_class": "value_blind_held_out_dataset_ready",
                "absolute_depth_or_depth_per_dose": True,
                "point_count": len(yoshie_rows),
                "digitization_visual_gate_passed": True,
                "feature_profile_test": True,
                "value_blind_held_out": True,
                "model_predictions_completed": 0,
                "formal_pass_label": False,
                "boundary_strength": (
                    "same-reactor blanket rates, phase-resolved bulk electron "
                    "density and OES; no species-resolved wafer flux or IEAD"),
                "required_model_state": (
                    "phase-resolved reactor boundary plus persistent "
                    "C/F/S-containing surface memory"),
                "scope": (
                    "49 held-out feature depths are ready but none has been "
                    "predicted or graded"),
            },
        ],
        "summary": {
            "independent_non_krueger_chemistry_families_with_absolute_surface_depth_evidence": 2,
            "formal_held_out_feature_or_profile_predictions_completed": 0,
            "formal_held_out_feature_or_profile_passes": 0,
            "value_blind_held_out_feature_depth_targets_ready": len(yoshie_rows),
            "krueger_absolute_feature_depth_passed": False,
            "current_strongest_defensible_statement": (
                "absolute surface depth-per-dose transfers without a target "
                "fit in two non-Krueger chemistry families; no chemistry yet "
                "has a completed value-blind held-out feature-depth pass"),
        },
    }


def audit_text() -> str:
    return json.dumps(build_audit(), indent=2) + "\n"


def main() -> None:
    payload = audit_text()
    if OUTPUT.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed cross-chemistry depth audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
