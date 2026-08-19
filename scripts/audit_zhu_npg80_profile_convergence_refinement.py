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
        "frozen_gates": gates,
        "gate_results": gate_results,
        "all_numerical_gates_pass": all(gate_results.values()),
        "bottom_cd_numerically_certified_for_sentinel": all(
            gate_results[name]
            for name in (
                "cpu_cuda_cd", "fine_timestep_cd", "ultrafine_grid_cd",
                "ultrafine_xy_symmetry",
            )
        ),
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
