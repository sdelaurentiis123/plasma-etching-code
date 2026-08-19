#!/usr/bin/env python3
"""Run and certify one target-free Oxford evolving-profile convergence sentinel."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


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


PREREGISTRATION = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "profile_convergence_preregistration.json"
)
OUTPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_profile_convergence_v1"
    / "audit.json"
)
CASE_DIR = OUTPUT.parent / "cases"


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _render(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _case_path(name):
    return CASE_DIR / f"{name}.json"


def _case_inputs(name):
    preregistration = _load(PREREGISTRATION)
    profile_preregistration = _load(PROFILE_PREREGISTRATION)
    reactor = _load(REACTOR_DOSE)
    sentinel = preregistration["sentinel"]
    case = preregistration["cases"][name]
    scenarios = {
        row["name"]: row
        for row in _scenario_inputs(profile_preregistration, reactor)
    }
    scenario = scenarios[sentinel["transport_scenario"]]
    return preregistration, profile_preregistration, sentinel, case, scenario


def run_case(name, *, device):
    (preregistration, profile_preregistration, sentinel, case,
     scenario) = _case_inputs(name)
    profile = _run_profile_dose_trajectory(
        width_nm=float(sentinel["width_nm"]),
        scenario=scenario,
        target_rates_nm_min=(float(sentinel["blanket_rate_nm_min"]),),
        duration_s=float(sentinel["duration_s"]),
        dx_nm=float(case["mesh_spacing_nm"]),
        preregistration=profile_preregistration,
        maximum_step_s=float(case["maximum_step_s"]),
        transport_device=device,
    )[0]
    if profile["tio2_clearance_detected"]:
        raise RuntimeError("convergence sentinel unexpectedly cleared the film")
    payload = {
        "schema": "petch.zhu-npg80-profile-convergence-case.v1",
        "case_name": name,
        "target_sem_used": False,
        "target_depth_used": False,
        "execution_device": str(device),
        "case_specification": dict(case),
        "inputs": {
            "convergence_preregistration_sha256": _hash(PREREGISTRATION),
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


def _cd_vector(case):
    profile = case["profile"]["profile"]
    return np.asarray([
        profile["top_cd_nm"],
        profile["middle_cd_nm"],
        profile["bottom_cd_nm"],
    ], dtype=float)


def _maximum_xy_asymmetry(case):
    return max(
        abs(float(row["width_x_nm"]) - float(row["width_y_nm"]))
        for row in case["profile"]["profile"]["cross_section"]
    )


def _comparison(reference, candidate):
    reference_depth = float(reference["profile"]["profile"]["etched_depth_nm"])
    candidate_depth = float(candidate["profile"]["profile"]["etched_depth_nm"])
    return {
        "reference_case": reference["case_name"],
        "candidate_case": candidate["case_name"],
        "depth_relative_change": abs(candidate_depth - reference_depth)
        / max(abs(candidate_depth), np.finfo(float).tiny),
        "cd_absolute_change_nm": np.abs(
            _cd_vector(candidate) - _cd_vector(reference)
        ).tolist(),
        "maximum_cd_absolute_change_nm": float(np.max(np.abs(
            _cd_vector(candidate) - _cd_vector(reference)
        ))),
    }


def build_audit():
    preregistration = _load(PREREGISTRATION)
    cases = {
        name: _load(_case_path(name))
        for name in preregistration["cases"]
    }
    timestep = _comparison(cases["coarse_dt8"], cases["coarse_dt4"])
    grid = _comparison(cases["coarse_dt4"], cases["fine_dt4"])
    fine_asymmetry = _maximum_xy_asymmetry(cases["fine_dt4"])
    gates = preregistration["gates"]
    maximum_balance = max(
        float(case["profile"][
            "maximum_transport_relative_particle_balance_error"])
        for case in cases.values()
    )
    maximum_remap = max(
        float(case["profile"][
            "maximum_state_remap_relative_conservation_residual"])
        for case in cases.values()
    )
    gate_results = {
        "timestep_depth": timestep["depth_relative_change"]
        <= gates["maximum_timestep_depth_relative_change"],
        "timestep_cd": timestep["maximum_cd_absolute_change_nm"]
        <= gates["maximum_timestep_cd_absolute_change_nm"],
        "grid_depth": grid["depth_relative_change"]
        <= gates["maximum_grid_depth_relative_change"],
        "grid_cd": grid["maximum_cd_absolute_change_nm"]
        <= gates["maximum_grid_cd_absolute_change_nm"],
        "fine_grid_xy_symmetry": fine_asymmetry
        <= gates["maximum_fine_grid_xy_asymmetry_nm"],
        "particle_balance": maximum_balance
        <= gates["maximum_transport_relative_particle_balance_error"],
        "state_remap_conservation": maximum_remap
        <= gates["maximum_state_remap_relative_conservation_residual"],
    }
    return {
        "schema": "petch.zhu-npg80-profile-convergence.v1",
        "condition_id": preregistration["condition_id"],
        "target_sem_used": False,
        "target_depth_used": False,
        "surface_law_status": (
            "cross-machine conditional analog; numerical convergence does not "
            "validate the Oxford surface law"
        ),
        "inputs": {
            "convergence_preregistration": {
                "path": str(PREREGISTRATION.relative_to(ROOT)),
                "sha256": _hash(PREREGISTRATION),
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
        "timestep_comparison": timestep,
        "grid_comparison": grid,
        "fine_grid_maximum_xy_asymmetry_nm": float(fine_asymmetry),
        "maximum_transport_relative_particle_balance_error": maximum_balance,
        "maximum_state_remap_relative_conservation_residual": maximum_remap,
        "frozen_gates": gates,
        "gate_results": gate_results,
        "all_numerical_gates_pass": all(gate_results.values()),
        "supports_absolute_target_profile_prediction": False,
        "supports_atomic_accuracy": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case", choices=("coarse_dt8", "coarse_dt4", "fine_dt4")
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
            raise SystemExit("profile convergence audit is stale")
        print(f"PASS {OUTPUT.relative_to(ROOT)}")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
