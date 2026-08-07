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
LEVINSON = (
    ROOT / "results" / "curated"
    / "levinson_1997_feature_identifiability" / "audit.json")
WOO = (
    ROOT / "results" / "curated"
    / "woo_2024_c4f6_board" / "audit.json")
MAHOROWALA = (
    ROOT / "data" / "experimental" / "mahorowala_1998_cl2"
    / "audit_manifest.json")
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
    levinson = _load(LEVINSON)
    woo = _load(WOO)
    mahorowala = _load(MAHOROWALA)
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
        "schema": "petch.cross-chemistry-depth-evidence.v2",
        "audit_id": "CROSS-CHEMISTRY-ABSOLUTE-DEPTH-2026-08-06-R2",
        "claim_rule": (
            "surface depth-per-dose, development replay, evaluated feature "
            "depth, a source fixed-time rate/profile board, and value-blind "
            "held-out feature prediction are distinct evidence classes and "
            "may not be relabeled as one another"),
        "source_receipts": {
            "tinacba_2021_sf5_depth_sha256": _sha256(TINACBA),
            "vella_hao_ale_depth_sha256": _sha256(VELLA),
            "krueger_depth_identifiability_sha256": _sha256(KRUEGER),
            "yoshie_visual_audit_sha256": _sha256(YOSHIE_VISION),
            "yoshie_feature_table_sha256": _sha256(YOSHIE_FEATURES),
            "levinson_1997_feature_identifiability_sha256":
                _sha256(LEVINSON),
            "woo_2024_c4f6_board_sha256": _sha256(WOO),
            "mahorowala_1998_cl2_fixed_time_sha256": _sha256(MAHOROWALA),
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
                "board": "Levinson 1997 Ar+/Cl2/Si MIBE features",
                "chemistry_and_material": (
                    "monoenergetic Ar+ plus molecular Cl2 on masked Si"),
                "evidence_class": (
                    "controlled_beam_feature_profiles_missing_time_fluence"),
                "absolute_depth_or_depth_per_dose": False,
                "point_count": len(levinson["figure11_cases"]),
                "reported_absolute_feature_depths": True,
                "surface_model_independent_of_feature_profiles": (
                    levinson["surface_closure"][
                        "independent_of_feature_profiles"]),
                "fit_to_compared_depth_or_yield": False,
                "boundary_strength": (
                    "monoenergetic normal Ar+ and isotropic Cl2 with measured "
                    "sample-position current, but Figure 11 omits the "
                    "case-specific current/fluence and exposure time"),
                "feature_profile_test": False,
                "original_pixels_archived": levinson["source"][
                    "source_pixels_archived"],
                "value_blind_held_out": False,
                "formal_pass_label": False,
                "scope": levinson["verdict"],
            },
            {
                "board": "Woo 2024 CF4/C4F6/He patterned SiO2",
                "chemistry_and_material": (
                    "CF4/C4F6/He ICP on ACL-masked SiO2 lines"),
                "evidence_class": (
                    "absolute_patterned_rate_board_missing_kinetic_boundary"),
                "absolute_depth_or_depth_per_dose": True,
                "point_count": woo["quantitative_rate_points"],
                "reported_absolute_feature_depths": True,
                "fit_to_compared_depth_or_yield": False,
                "boundary_strength": (
                    "same-reactor Te, aggregate ion current, self-bias, "
                    "relative OES, XPS, and patterned rates; no "
                    "species-resolved ion flux, IEAD, or absolute neutrals"),
                "feature_profile_test": False,
                "original_pixels_archived": True,
                "value_blind_held_out": False,
                "model_predictions_completed": 0,
                "formal_pass_label": False,
                "source_internal_consistency_passed": False,
                "scope": woo["verdict"],
            },
            {
                "board": "Mahorowala 1998 Cl2/oxide-mask Si features",
                "chemistry_and_material": (
                    "pure Cl2 on oxide-masked poly-Si"),
                "evidence_class": (
                    "fixed_time_absolute_rate_and_profile_board_missing_"
                    "species_iead_boundary"),
                "absolute_depth_or_depth_per_dose": True,
                "source_run_count": mahorowala["transcription"][
                    "source_run_count"],
                "point_count": mahorowala["transcription"][
                    "usable_quantitative_run_count"],
                "reported_absolute_feature_depths": True,
                "derived_poly_si_depth_range_nm": mahorowala[
                    "transcription"]["derived_poly_si_depth_range_nm"],
                "fit_to_compared_depth_or_yield": False,
                "boundary_strength": (
                    "fixed 75 s exposure and absolute rates for every usable "
                    "run; reactor knobs and corresponding SEM panels are "
                    "published, but species-resolved wafer flux, measured "
                    "IEAD, and measured IAD are not"),
                "source_feature_profiles_available": True,
                "source_profile_panel_count": mahorowala[
                    "source_locations"]["profile_montage"]["panel_count"],
                "feature_profile_test": False,
                "original_pixels_archived": False,
                "original_resolution_visual_audit_passed": (
                    mahorowala["transcription"]["visual_audit_status"]
                    == "passed_original_resolution"),
                "value_blind_held_out": False,
                "model_predictions_completed": 0,
                "formal_pass_label": False,
                "scope": (
                    "absolute fixed-time chlorine rate/profile targets are "
                    "ready to grade an independently identified reactor "
                    "boundary; the scored rates may not be used to select "
                    "that boundary and then be counted as predictions"),
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
            "controlled_feature_boards_blocked_by_missing_time_or_fluence": 1,
            "fixed_time_chlorine_absolute_rate_targets_ready": mahorowala[
                "transcription"]["usable_quantitative_run_count"],
            "fixed_time_chlorine_sem_profile_targets_ready": mahorowala[
                "source_locations"]["profile_montage"]["panel_count"],
            "c4f6_same_reactor_absolute_patterned_rate_targets_ready": (
                woo["quantitative_rate_points"]),
            "value_blind_held_out_feature_depth_targets_ready": len(yoshie_rows),
            "krueger_absolute_feature_depth_passed": False,
            "current_strongest_defensible_statement": (
                "absolute surface depth-per-dose transfers without a target "
                "fit exist in two non-Krueger chemistry families and an "
                "11-point fixed-time chlorine rate/profile board is ready; "
                "no chemistry yet has a completed value-blind held-out "
                "feature-depth pass"),
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
