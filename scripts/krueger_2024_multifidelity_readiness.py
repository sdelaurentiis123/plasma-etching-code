#!/usr/bin/env python3
"""Audit whether Krueger evidence can legally issue a new calibration proposal.

The answer is derived from base-only receipts.  The tool cannot read the
oxygen/power transfer table, cannot launch a simulation, and cannot emit a new
parameter pair.  Its purpose is to distinguish a calibrated multi-fidelity
model from a collection of coarse and short-time numbers that do not yet span
the required current-epoch endpoint response.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "data" / "experimental" / "krueger_2024" / "base_case_metrics.csv"
CURRENT_EPOCH_SOURCES = (
    "src/petch/feature_step_3d.py",
    "src/petch/surface_mesh_3d.py",
    "src/petch/material_mechanism_3d.py",
    "src/petch/surface_kinetics.py",
    "src/petch/boundary_transport_3d.py",
    "src/petch/amorphous_carbon_mask.py",
)
PARAMETERS = (
    "effective_mask_crosslinked_growth_fraction",
    "oxide_etch_yield_scale",
)


def _sha(path: Path | str) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _canonical_sha(payload: dict, field: str) -> str:
    content = dict(payload)
    content.pop(field, None)
    return sha256(json.dumps(
        content, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _read(path: Path | str) -> tuple[Path, dict]:
    path = Path(path)
    return path, json.loads(path.read_text(encoding="utf-8"))


def _targets() -> tuple[dict, dict]:
    with TARGETS.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        row["metric"]: float(row["value"])
        for row in rows
        if row["split"] == "calibration" and row["evidence_type"] == "experiment"
        and row["metric"] in {"mask_opening", "etch_depth"}
    }
    if set(selected) != {"mask_opening", "etch_depth"}:
        raise ValueError("base target table does not contain exactly two calibration observations")
    return selected, {
        "path_name": TARGETS.name,
        "sha256": _sha(TARGETS),
        "selection": "split=calibration and evidence_type=experiment only",
    }


def current_source_manifest() -> dict:
    manifest = {}
    for relative in CURRENT_EPOCH_SOURCES:
        path = ROOT / relative
        if not path.is_file():
            raise ValueError(f"current operator source is missing: {relative}")
        manifest[relative] = _sha(path)
    return manifest


def _git_revision_and_clean() -> tuple[str, bool]:
    revision = subprocess.run(
        ("git", "rev-parse", "HEAD"), cwd=ROOT, check=True,
        capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(
        ("git", "status", "--porcelain"), cwd=ROOT, check=True,
        capture_output=True, text=True).stdout
    return revision, not bool(dirty)


def _verified_remap_selection(remap_audit_path):
    """Return a selected remapper only from a paired, complete 5 nm base-only receipt."""
    remap_audit_path, audit = _read(remap_audit_path)
    configuration = audit.get("configuration", {})
    if (audit.get("schema") != "petch.krueger_2024_remap_backend_audit.v1"
            or audit.get("held_out_profile_data_read") is not False
            or audit.get("scientific_scope")
            != "bounded base-boundary operator selection; no experimental outcomes read"
            or not np.isclose(float(configuration.get("dx_nm", np.nan)), 5.0)
            or int(configuration.get("steps", -1)) != 2
            or not np.isclose(
                float(configuration.get("total_physical_time_s", np.nan)), 0.05)):
        raise ValueError("remap receipt is not the bounded 5 nm base-only comparison")
    required = ("indexed_knn", "common_refinement")
    if not set(required).issubset(configuration.get("backends", ())):
        raise ValueError("remap receipt does not contain both required paired backends")
    cases = audit.get("cases", {})
    initial_geometry = set()
    initial_state = set()
    first_geometry = set()
    maximum_residual = {}
    maximum_ledger = {}
    for backend in required:
        case = cases.get(backend, {})
        steps = case.get("steps", ())
        if case.get("status") != "complete" or len(steps) != 2:
            raise ValueError(f"remap backend did not complete its paired receipt: {backend}")
        initial_geometry.add(case.get("initial_geometry_sha256"))
        initial_state.add(case.get("initial_state", {}).get("sha256"))
        first_geometry.add(steps[0].get("geometry_sha256"))
        if any(step.get("topology_event") is not None for step in steps):
            raise ValueError(f"remap receipt crossed a topology event: {backend}")
        if any(step.get("remap", {}).get("surface_state_remap_backend") != backend
               for step in steps):
            raise ValueError(f"remap receipt mislabeled its operator: {backend}")
        maximum_residual[backend] = max(float(
            step["maximum_remap_relative_conservation_residual"])
            for step in steps)
        maximum_ledger[backend] = max(abs(float(
            step["operator"]["maximum_material_ledger_residual_units_m2"]))
            for step in steps)
        if maximum_residual[backend] > 1e-12 or maximum_ledger[backend] > 1e-20:
            raise ValueError(f"remap receipt does not close conservation: {backend}")
    if (len(initial_geometry) != 1 or None in initial_geometry
            or len(initial_state) != 1 or None in initial_state
            or len(first_geometry) != 1 or None in first_geometry):
        raise ValueError("remap comparison is not paired through the first profile step")
    return remap_audit_path, {
        "selected_backend": "common_refinement",
        "selection_basis": (
            "paired 5 nm completion against indexed nearest-surface transfer; common "
            "refinement retained as the conservative geometric authority"),
        "maximum_remap_relative_conservation_residual": maximum_residual,
        "maximum_material_ledger_residual_units_m2": maximum_ledger,
    }


def derive(launch_path, evaluation_path, run_audit_path, multiresolution_path,
           cuda_summary_path, *, current_revision, current_sources,
           remap_audit_path=None):
    launch_path, launch = _read(launch_path)
    evaluation_path, evaluation = _read(evaluation_path)
    run_audit_path, run_audit = _read(run_audit_path)
    multiresolution_path, multiresolution = _read(multiresolution_path)
    cuda_summary_path, cuda_summary = _read(cuda_summary_path)

    if (launch.get("schema") != "petch.krueger-2024.r19-response-check-launch.v1"
            or launch.get("protocol_id") != "K24-PETCH-R1.9"
            or launch.get("authority") is not False
            or launch.get("held_out_profile_data_read") is not False
            or launch.get("ten_nm_sequence_closes_after_this_evaluation") is not True):
        raise ValueError("launch is not the sealed final 10 nm response check")
    if (evaluation.get("schema")
            != "petch.krueger-2024.r19-response-check-evaluation.v1"
            or evaluation.get("held_out_profile_data_read") is not False
            or evaluation.get("ten_nm_sequence_closed") is not True
            or evaluation.get("evaluation_sha256")
            != _canonical_sha(evaluation, "evaluation_sha256")):
        raise ValueError("response evaluation is not a valid closed receipt")
    if evaluation["inputs"]["launch_manifest"]["sha256"] != _sha(launch_path):
        raise ValueError("response evaluation does not bind the supplied launch")
    if evaluation["inputs"]["run_audit"]["sha256"] != _sha(run_audit_path):
        raise ValueError("response evaluation does not bind the supplied run audit")
    if (evaluation.get("decision") != "reject_response_model"
            or not evaluation.get("trajectory_contract", {}).get("pass")):
        raise ValueError("readiness branch expects the completed rejected R1.9 response model")

    target, target_info = _targets()
    if launch.get("base_only_input_sha256", {}).get(TARGETS.name) != target_info["sha256"]:
        raise ValueError("R1.9 launch is not bound to the current base calibration table")
    if (not np.isclose(float(launch["target_nm"]["mask_opening"]), target["mask_opening"])
            or not np.isclose(float(launch["target_nm"]["etch_depth"]), target["etch_depth"])):
        raise ValueError("R1.9 target differs from the declared base calibration observations")
    configuration = run_audit.get("configuration", {})
    if (run_audit.get("status") != "complete"
            or configuration.get("boundary_case") != "base"
            or not np.isclose(float(configuration.get("duration_s", np.nan)), 60.0)
            or not np.isclose(float(configuration.get("dx_um", np.nan)), 0.01)
            or any(not np.isclose(float(configuration.get(name, np.nan)),
                                  float(launch["candidate"][name]), rtol=0.0, atol=2e-14)
                   for name in PARAMETERS)):
        raise ValueError("R1.9 audit is not the sealed complete base response")

    if (cuda_summary.get("schema") != "petch.krueger_2024_cuda_profile_summary.v1"
            or cuda_summary.get("held_out_profile_data_read") is not False
            or cuda_summary.get("calibration_performed") is not False
            or cuda_summary.get("summary_sha256")
            != _canonical_sha(cuda_summary, "summary_sha256")):
        raise ValueError("CUDA/refinement summary is invalid or not base-only")
    multires_reference = cuda_summary["inputs"]["paired_10nm_5nm_initial_audit"]
    if multires_reference["sha256"] != _sha(multiresolution_path):
        raise ValueError("CUDA/refinement summary does not bind the supplied 10/5 audit")
    paired = multiresolution.get("paired_10nm_vs_5nm", {}).get("initial", {})
    if not {"etch_depth_nm_s", "mask_opening_nm_s"}.issubset(paired):
        raise ValueError("10/5 audit lacks paired initial calibration rates")

    old_sources = dict(launch.get("executable_source_sha256", {}))
    shared = sorted(set(old_sources) & set(current_sources))
    changed = {
        name: {"r19": old_sources[name], "current": current_sources[name]}
        for name in shared if old_sources[name] != current_sources[name]
    }
    added = sorted(set(current_sources) - set(old_sources))
    operator_epoch_compatible = bool(shared and not changed and not added)
    r19_epoch_canonical = json.dumps(
        old_sources, sort_keys=True, separators=(",", ":")).encode("utf-8")
    current_epoch_canonical = json.dumps(
        current_sources, sort_keys=True, separators=(",", ":")).encode("utf-8")

    remap_backend = cuda_summary["operator_receipts"]["surface_state_remap_backend"]
    remap_selection = None
    remap_audit = None
    if remap_audit_path is not None:
        remap_audit, remap_selection = _verified_remap_selection(remap_audit_path)
        remap_backend = remap_selection["selected_backend"]
    blockers = []
    if "legacy_knn" in remap_backend:
        blockers.append({
            "code": "remap_backend_not_selected",
            "evidence": remap_backend,
            "required_resolution": (
                "run one bounded, same-state legacy/indexed/common-refinement comparison and "
                "freeze the chosen backend in every operator manifest"),
        })
    if not operator_epoch_compatible:
        blockers.append({
            "code": "r19_response_belongs_to_prior_operator_epoch",
            "evidence": {"changed_shared_sources": changed, "new_current_sources": added},
            "required_resolution": "do not reuse the rejected R1.9 endpoint Jacobian as current",
        })
    blockers.extend((
        {
            "code": "no_current_epoch_high_fidelity_endpoint_anchor",
            "evidence": {"count": 0},
            "required_resolution": (
                "after remap selection, evaluate one clean t=0 5 nm or certified-AMR base "
                "anchor at the fixed R1.9 parameter pair"),
        },
        {
            "code": "no_paired_endpoint_discrepancy",
            "evidence": {
                "paired_endpoint_count": 0,
                "available_pair": "0.5 s initial rates only",
            },
            "required_resolution": (
                "a short initial rate pair may size numerical uncertainty but may not replace a "
                "60 s low/high endpoint discrepancy in this path-nonlinear process"),
        },
        {
            "code": "first_order_current_epoch_response_not_identified",
            "evidence": {
                "current_epoch_low_endpoint_count": 0,
                "default_affine_low_requirement_for_two_parameters": 3,
                "default_paired_discrepancy_requirement": (
                    "center plus two independent directions"),
            },
            "required_resolution": (
                "if the clean fine anchor misses, identify a current-epoch physical/empirical "
                "direction under an amended bounded protocol before proposing a correction"),
        },
    ))

    payload = {
        "schema": "petch.krueger_2024_multifidelity_readiness.v1",
        "status": "blocked_before_parameter_proposal",
        "scientific_status": (
            "base-only model-management readiness audit; no simulation, parameter proposal, "
            "calibration authority, or held-out reveal"),
        "authority": False,
        "held_out_profile_data_read": False,
        "calibration_performed": False,
        "generator": {"path_name": Path(__file__).name, "sha256": _sha(Path(__file__))},
        "base_targets": {**target_info, "values": target},
        "current_center": {
            "parameters": {name: float(launch["candidate"][name]) for name in PARAMETERS},
            "r19_coarse_response_nm": dict(evaluation["actual_nm"]),
            "r19_target_error_nm": dict(evaluation["target_error_nm"]),
            "role": (
                "fixed parameter location for the first clean high-fidelity anchor; not a "
                "newly optimized or authoritative pair"),
        },
        "operator_epochs": {
            "r19": {
                "sha256": sha256(r19_epoch_canonical).hexdigest(),
                "sources": old_sources,
            },
            "current": {
                "git_revision": str(current_revision),
                "sha256": sha256(current_epoch_canonical).hexdigest(),
                "sources": dict(current_sources),
            },
            "compatible": operator_epoch_compatible,
            "changed_shared_sources": changed,
            "new_current_sources": added,
        },
        "short_refinement_evidence": {
            "physical_interval_s": 0.5,
            "depth_rate_absolute_relative_percent": abs(float(
                paired["etch_depth_nm_s"]["relative_to_fine"])) * 100.0,
            "opening_rate_absolute_relative_percent": abs(float(
                paired["mask_opening_nm_s"]["relative_to_fine"])) * 100.0,
            "maximum_width_rate_absolute_relative_percent": abs(float(
                paired["maximum_feature_width_nm_s"]["relative_to_fine"])) * 100.0,
            "use": "initial numerical scale and cost evidence only",
            "not_earned": "late endpoint discrepancy or local profile-shape authority",
        },
        "model_management_contract": {
            "high_fidelity_remains_in_loop": True,
            "center_value_consistency_required": True,
            "first_order_discrepancy_required_by_default": True,
            "trial_acceptance": "actual/predicted high-fidelity merit reduction ratio",
            "stale_operator_response_reuse": "refused",
            "short_rate_as_endpoint_discrepancy": "refused",
        },
        "remap_operator_selection": (
            remap_selection if remap_selection is not None else {
                "selected_backend": None,
                "status": "awaiting paired bounded 5 nm receipt",
            }),
        "blockers": blockers,
        "next_bounded_sequence": ([
            "explicitly select the remap backend with a same-state short comparison",
        ] if remap_selection is None else []) + [
            "run one clean current-epoch high-fidelity base anchor at the fixed R1.9 pair",
            "if the anchor meets base tolerances, freeze without another parameter step",
            "if it misses, update the discrepancy receipt and earn a current-epoch direction "
            "before one safeguarded high-fidelity correction",
        ],
        "inputs": {
            "launch": {"path_name": launch_path.name, "sha256": _sha(launch_path)},
            "evaluation": {
                "path_name": evaluation_path.name, "sha256": _sha(evaluation_path)},
            "r19_run": {"path_name": run_audit_path.name, "sha256": _sha(run_audit_path)},
            "multiresolution": {
                "path_name": multiresolution_path.name, "sha256": _sha(multiresolution_path)},
            "cuda_summary": {
                "path_name": cuda_summary_path.name, "sha256": _sha(cuda_summary_path)},
            **({"remap_audit": {
                "path_name": remap_audit.name, "sha256": _sha(remap_audit)}}
               if remap_audit is not None else {}),
        },
    }
    payload["readiness_sha256"] = _canonical_sha(payload, "readiness_sha256")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--launch", required=True)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--r19-run", required=True)
    parser.add_argument("--multiresolution", required=True)
    parser.add_argument("--cuda-summary", required=True)
    parser.add_argument("--remap-audit")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    revision, clean = _git_revision_and_clean()
    if not clean:
        raise RuntimeError("readiness artifact requires a clean checksum-bound worktree")
    payload = derive(
        args.launch, args.evaluation, args.r19_run, args.multiresolution,
        args.cuda_summary, current_revision=revision,
        current_sources=current_source_manifest(),
        remap_audit_path=args.remap_audit)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(output)
    print(json.dumps({
        "status": payload["status"],
        "blockers": [item["code"] for item in payload["blockers"]],
        "next": payload["next_bounded_sequence"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
