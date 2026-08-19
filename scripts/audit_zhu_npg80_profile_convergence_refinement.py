#!/usr/bin/env python3
"""Resolve the Oxford sentinel's 10 nm bottom-CD convergence failure."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_zhu_npg80_conditional_profiles import (
    ANALOG_BOARD,
    FILM_THICKNESS_UM,
    PREREGISTRATION as PROFILE_PREREGISTRATION,
    REACTOR_DOSE,
    _run_profile_dose_trajectory,
    _scenario_inputs,
)
from scripts.audit_zhu_npg80_profile_convergence import (
    _comparison,
    _hash,
    _load,
    _maximum_xy_asymmetry,
    _render,
)


PREREGISTRATION = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "profile_convergence_refinement_preregistration.json"
)
OUTPUT = (
    ROOT / "results" / "curated"
    / "zhu_npg80_profile_convergence_refinement_v1" / "audit.json"
)
CASE_DIR = OUTPUT.parent / "cases"
JI_PHYSICAL_SHAPE_WITNESS = (
    ROOT / "research_sources" / "library" / "ji-2024-tio2-hierarchical.md"
)


def _case_path(name):
    return CASE_DIR / f"{name}.json"


def run_case(name, *, device):
    preregistration = _load(PREREGISTRATION)
    profile_preregistration = _load(PROFILE_PREREGISTRATION)
    reactor = _load(REACTOR_DOSE)
    sentinel = preregistration["sentinel"]
    case = preregistration["cases"][name]
    required_device = case.get("required_execution_device")
    if required_device is not None and str(device) != str(required_device):
        raise ValueError(
            f"case {name} requires execution device {required_device!r}"
        )
    scenarios = {
        row["name"]: row
        for row in _scenario_inputs(profile_preregistration, reactor)
    }
    profile = _run_profile_dose_trajectory(
        width_nm=float(sentinel["width_nm"]),
        scenario=scenarios[sentinel["transport_scenario"]],
        target_rates_nm_min=(float(sentinel["blanket_rate_nm_min"]),),
        duration_s=float(sentinel["duration_s"]),
        dx_nm=float(case["mesh_spacing_nm"]),
        preregistration=profile_preregistration,
        maximum_step_s=float(case["maximum_step_s"]),
        transport_device=device,
    )[0]
    if profile["tio2_clearance_detected"]:
        raise RuntimeError("refinement sentinel unexpectedly cleared the film")
    payload = {
        "schema": "petch.zhu-npg80-profile-convergence-refinement-case.v1",
        "case_name": name,
        "target_sem_used": False,
        "target_depth_used": False,
        "execution_device": str(device),
        "case_specification": dict(case),
        "inputs": {
            "refinement_preregistration_sha256": _hash(PREREGISTRATION),
            "profile_preregistration_sha256": _hash(PROFILE_PREREGISTRATION),
            "reactor_dose_sha256": _hash(REACTOR_DOSE),
            "analog_board_sha256": _hash(ANALOG_BOARD),
        },
        "profile": profile,
    }
    path = _case_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render(payload), encoding="utf-8")
    print(path.relative_to(ROOT))


def build_audit():
    preregistration = _load(PREREGISTRATION)
    reference_path = ROOT / preregistration["reference_case"]["source"]
    reference = _load(reference_path)
    cases = {
        name: _load(_case_path(name))
        for name in preregistration["cases"]
    }
    device_replicate = _comparison(
        reference, cases["fine_dt4_cuda_replicate"]
    )
    timestep = _comparison(
        cases["fine_dt4_cuda_replicate"], cases["fine_dt2"]
    )
    grid = _comparison(cases["fine_dt2"], cases["ultrafine_dt2"])
    ultrafine_asymmetry = _maximum_xy_asymmetry(cases["ultrafine_dt2"])
    gates = preregistration["gates"]
    all_cases = (reference, *cases.values())
    maximum_balance = max(
        float(case["profile"][
            "maximum_transport_relative_particle_balance_error"])
        for case in all_cases
    )
    maximum_remap = max(
        float(case["profile"][
            "maximum_state_remap_relative_conservation_residual"])
        for case in all_cases
    )
    gate_results = {
        "cpu_cuda_depth": device_replicate["depth_relative_change"]
        <= gates["maximum_cpu_cuda_depth_relative_change"],
        "cpu_cuda_cd": device_replicate["maximum_cd_absolute_change_nm"]
        <= gates["maximum_cpu_cuda_cd_absolute_change_nm"],
        "fine_timestep_depth": timestep["depth_relative_change"]
        <= gates["maximum_fine_timestep_depth_relative_change"],
        "fine_timestep_cd": timestep["maximum_cd_absolute_change_nm"]
        <= gates["maximum_fine_timestep_cd_absolute_change_nm"],
        "ultrafine_grid_depth": grid["depth_relative_change"]
        <= gates["maximum_ultrafine_grid_depth_relative_change"],
        "ultrafine_grid_cd": grid["maximum_cd_absolute_change_nm"]
        <= gates["maximum_ultrafine_grid_cd_absolute_change_nm"],
        "ultrafine_xy_symmetry": ultrafine_asymmetry
        <= gates["maximum_ultrafine_xy_asymmetry_nm"],
        "particle_balance": maximum_balance
        <= gates["maximum_transport_relative_particle_balance_error"],
        "state_remap_conservation": maximum_remap
        <= gates["maximum_state_remap_relative_conservation_residual"],
    }
    bottom_cd_certified = all(
        gate_results[name]
        for name in (
            "cpu_cuda_cd", "fine_timestep_cd", "ultrafine_grid_cd",
            "ultrafine_xy_symmetry",
        )
    )
    ultrafine_metrics = cases["ultrafine_dt2"]["profile"]["profile"]
    fine_metrics = cases["fine_dt2"]["profile"]["profile"]
    fine_section = fine_metrics["cross_section"]
    ultrafine_section = ultrafine_metrics["cross_section"]
    fine_fraction = [
        float(row["relief_fraction_from_top"]) for row in fine_section
    ]
    ultrafine_fraction = [
        float(row["relief_fraction_from_top"]) for row in ultrafine_section
    ]
    if fine_fraction != ultrafine_fraction:
        raise RuntimeError("refinement cross-section sampling fractions changed")
    section_change = [
        abs(float(candidate["mean_width_nm"]) - float(reference["mean_width_nm"]))
        for reference, candidate in zip(fine_section, ultrafine_section)
    ]
    body_change = [
        change for fraction, change in zip(fine_fraction, section_change)
        if fraction <= 0.75
    ]
    junction_change = [
        change for fraction, change in zip(fine_fraction, section_change)
        if fraction >= 0.80
    ]
    first_over_bottom_gate = next((
        fraction for fraction, change in zip(fine_fraction, section_change)
        if change > gates["maximum_ultrafine_grid_cd_absolute_change_nm"]
    ), None)
    remaining_film_nm = (
        1.0e3 * float(FILM_THICKNESS_UM)
        - float(ultrafine_metrics["etched_depth_nm"])
    )
    bottom_minus_top_nm = (
        float(ultrafine_metrics["bottom_cd_nm"])
        - float(ultrafine_metrics["top_cd_nm"])
    )
    qualitative_footing_present = bottom_minus_top_nm > 0.0
    if bottom_cd_certified:
        footing_classification = (
            "numerically certified geometric-transport result within the declared "
            "conditional model"
            if qualitative_footing_present else
            "no positive bottom-minus-top footing in the numerically certified "
            "conditional result"
        )
    else:
        footing_classification = (
            "bottom-CD magnitude remains numerically unresolved; do not classify "
            "the conditional footing quantitatively"
        )
    return {
        "schema": "petch.zhu-npg80-profile-convergence-refinement.v1",
        "condition_id": preregistration["condition_id"],
        "target_sem_used": False,
        "target_depth_used": False,
        "surface_law_status": (
            "cross-machine conditional analog; numerical refinement does not "
            "validate the Oxford surface law"
        ),
        "inputs": {
            "refinement_preregistration": {
                "path": str(PREREGISTRATION.relative_to(ROOT)),
                "sha256": _hash(PREREGISTRATION),
            },
            "reference_case": {
                "path": str(reference_path.relative_to(ROOT)),
                "sha256": _hash(reference_path),
            },
            "external_physical_shape_witness": {
                "path": str(JI_PHYSICAL_SHAPE_WITNESS.relative_to(ROOT)),
                "sha256": _hash(JI_PHYSICAL_SHAPE_WITNESS),
            },
            "cases": {
                name: {
                    "path": str(_case_path(name).relative_to(ROOT)),
                    "sha256": _hash(_case_path(name)),
                }
                for name in cases
            },
        },
        "sentinel": preregistration["sentinel"],
        "cpu_cuda_replicate_comparison": device_replicate,
        "fine_timestep_comparison": timestep,
        "ultrafine_grid_comparison": grid,
        "ultrafine_maximum_xy_asymmetry_nm": float(ultrafine_asymmetry),
        "maximum_transport_relative_particle_balance_error": maximum_balance,
        "maximum_state_remap_relative_conservation_residual": maximum_remap,
        "conditional_shape_diagnostic": {
            "ultrafine_top_cd_nm": float(ultrafine_metrics["top_cd_nm"]),
            "ultrafine_bottom_cd_nm": float(ultrafine_metrics["bottom_cd_nm"]),
            "ultrafine_bottom_minus_top_cd_nm": bottom_minus_top_nm,
            "qualitative_footing_present": qualitative_footing_present,
            "classification": footing_classification,
            "surface_mechanism_growth_enabled": False,
            "surface_product_redeposition_enabled": False,
            "material_motion_invariant": {
                "evolving_material": "ALD TiO2 external gas-solid boundary",
                "pinned_materials": ["Cr mask", "fused-silica substrate"],
                "allowed_normal_velocity": "removal-only and nonnegative",
                "negative_removal_velocity_rejected_by_mechanism": True,
                "wider_lower_section_requires_slow_lower-wall_recession_not_growth": True,
            },
            "candidate_mechanism_if_present": (
                "differential geometric ion dose from deterministic angular "
                "transport and feature self-shadowing"
            ),
            "external_physical_shape_witness": {
                "bibkey": "ji-2024-tio2-hierarchical",
                "observed_shape": (
                    "upper triangle over a wider rectangular lower TiO2 section"
                ),
                "reported_mechanisms": (
                    "Cr-mask lateral shrink plus lower-feature passivation"
                ),
                "same_feed_constituents": True,
                "same_reactor_or_condition": False,
                "coefficient_transfer_allowed": False,
            },
            "interpretation_limit": (
                "numerical certification can establish that the shape is not a "
                "discretization artifact of this conditional model; it cannot "
                "validate the omitted Oxford TiO2/Cr surface response or prove "
                "which physical mechanism made Freddie's withheld profile"
            ),
        },
        "post_result_localization_diagnostic": {
            "changes_frozen_gate_result": False,
            "exploratory_not_preregistered": True,
            "sentinel_cleared_tio2": bool(
                cases["ultrafine_dt2"]["profile"]["tio2_clearance_detected"]
            ),
            "remaining_tio2_below_etched_floor_nm": remaining_film_nm,
            "independent_pillar_bottom_exists": False,
            "reason": (
                "the pre-clear sidewall joins a continuous unetched TiO2 film; "
                "the deepest width samples measure that junction rather than a "
                "freestanding pillar base"
            ),
            "body_fraction_limit_from_top": 0.75,
            "body_maximum_10nm_to_5nm_width_change_nm": max(body_change),
            "near_floor_fraction_start_from_top": 0.80,
            "near_floor_maximum_10nm_to_5nm_width_change_nm": max(junction_change),
            "first_fraction_exceeding_frozen_5nm_cd_gate": first_over_bottom_gate,
            "body_profile_within_frozen_timestep_cd_tolerance": (
                max(body_change)
                <= gates["maximum_fine_timestep_cd_absolute_change_nm"]
            ),
            "physical_classification": (
                "the relief-to-unetched-film junction is physical in a partial "
                "etch, but its quantitative flare magnitude is not grid-certified"
            ),
            "nonphysical_growth_or_interface_drift_detected": False,
            "target_pillar_bottom_certified": False,
        },
        "frozen_gates": gates,
        "gate_results": gate_results,
        "all_numerical_gates_pass": all(gate_results.values()),
        "bottom_cd_numerically_certified_for_sentinel": bottom_cd_certified,
        "supports_absolute_target_profile_prediction": False,
        "supports_atomic_accuracy": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=("fine_dt4_cuda_replicate", "fine_dt2", "ultrafine_dt2"),
    )
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    selected = sum(bool(value) for value in (args.case, args.assemble, args.check))
    if selected != 1:
        parser.error("select exactly one of --case, --assemble, or --check")
    if args.case:
        run_case(args.case, device=args.device)
        return
    payload = build_audit()
    rendered = _render(payload)
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("profile convergence refinement audit is stale")
        print(f"PASS {OUTPUT.relative_to(ROOT)}")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
