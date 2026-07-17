#!/usr/bin/env python3
"""Bounded same-state audit of Krueger surface-state remap backends.

The remap is part of the physical operator because polymer coverage and
material inventories influence the *next* profile step.  This audit therefore
starts every backend from the same analytic geometry, state, seed epoch, and
timestep schedule.  It advances exactly two short steps: step one proves the
pre-remap transport/profile update is paired, while step two measures the first
downstream consequence of the transferred surface state.

This is an operator-selection gate, not a calibration or validation run.  It
does not read held-out Krueger outcomes and refuses more than 0.05 s of physical
etch time per backend.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import krueger_2024_multiresolution_audit as multires
import krueger_2024_trench_pilot as pilot
from petch.feature_step_3d import make_rectangular_trench_geometry_3d


DEFAULT_OUTPUT = ROOT / "results" / "krueger_2024_remap_backend_audit"
MAXIMUM_STEPS = 2
MAXIMUM_PHYSICAL_TIME_S = 0.05


def _write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            multires._jsonable(payload), indent=2, sort_keys=True,
            allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _array_digest(*arrays):
    digest = sha256()
    for value in arrays:
        array = np.ascontiguousarray(value)
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _geometry_digest(geometry):
    arrays = [np.asarray(geometry.phi), np.asarray(geometry.material_id)]
    arrays.extend(
        np.asarray(field)
        for _, field in sorted(dict(geometry.material_levelsets or {}).items())
    )
    return _array_digest(*arrays)


def _state_receipt(state, area):
    area = np.asarray(area, dtype=float)
    fields = dict(state.conservative_surface_fields())
    modes = dict(state.surface_field_remap_modes())
    digest = sha256()
    receipt = {}
    for name in sorted(fields):
        values = np.asarray(fields[name], dtype=float)
        digest.update(name.encode("utf-8"))
        digest.update(bytes.fromhex(_array_digest(values)))
        record = {
            "mode": str(modes[name]),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            "area_weighted_mean": float(np.dot(values, area) / area.sum()),
        }
        if modes[name] == "conservative":
            record["area_integral_mesh_units2"] = float(np.dot(values, area))
        receipt[name] = record
    return {"sha256": digest.hexdigest(), "fields": receipt}


def _remap_receipt(diagnostics):
    diagnostics = dict(diagnostics)
    materials = {}
    for material_id, value in dict(diagnostics.get("materials", {})).items():
        record = dict(value)
        materials[str(material_id)] = {
            name: record[name]
            for name in (
                "old_area_m2", "new_area_m2", "retained_area_m2",
                "removed_area_m2", "newly_exposed_area_m2",
                "max_relative_conservation_residual",
            )
            if name in record
        }
    geometry = dict(diagnostics.get("geometry_receipt", {}))
    scalar_geometry = {
        name: value for name, value in geometry.items()
        if np.isscalar(value) and not isinstance(value, (str, bytes))
    }
    return {
        "method": diagnostics.get("method"),
        "surface_state_remap_backend": diagnostics.get(
            "surface_state_remap_backend"),
        "old_topology": diagnostics.get("old_topology"),
        "new_topology": diagnostics.get("new_topology"),
        "topology_event": diagnostics.get("topology_event"),
        "fresh_surface_closure": diagnostics.get("fresh_surface_closure"),
        "total_removed_area_m2": diagnostics.get("total_removed_area_m2"),
        "total_newly_exposed_area_m2": diagnostics.get(
            "total_newly_exposed_area_m2"),
        "materials": materials,
        "geometry_receipt": scalar_geometry,
        "candidate_pair_count": diagnostics.get("candidate_pair_count"),
        "aligned_pair_count": diagnostics.get("aligned_pair_count"),
        "positive_pair_image_count": diagnostics.get("positive_pair_image_count"),
        "combined_pair_count": diagnostics.get("combined_pair_count"),
        "minimum_normal_dot": diagnostics.get("minimum_normal_dot"),
    }


def _maximum_remap_residual(diagnostics):
    values = [
        abs(float(record.get("max_relative_conservation_residual", 0.0)))
        for record in dict(diagnostics.get("materials", {})).values()
    ]
    return max(values, default=0.0)


def build_worker_command(args, backend, output):
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        "--output-dir", str(output),
        "--dx-nm", str(float(args.dx_nm)),
        "--steps", str(int(args.steps)),
        "--step-duration-s", str(float(args.step_duration_s)),
        "--seed", str(int(args.seed)),
        "--device", str(args.device),
        "--max-wall-s", str(float(args.max_wall_s)),
        "--worker-backend", str(backend),
    ]


def _initial_geometry(dx_nm):
    return make_rectangular_trench_geometry_3d(
        cell_width=0.13,
        cell_length=0.02,
        domain_height=2.8,
        dx=float(dx_nm) / 1000.0,
        opening_width=0.09,
        mask_thickness=0.85,
        substrate_top=1.8,
        etched_depth=0.0,
    )


def _worker(args):
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    backend = str(args.worker_backend)
    started = perf_counter()
    payload = {
        "schema": "petch.krueger_2024_remap_backend_worker.v1",
        "status": "running",
        "backend": backend,
        "configuration": {
            "dx_nm": float(args.dx_nm),
            "steps": int(args.steps),
            "step_duration_s": float(args.step_duration_s),
            "seed": int(args.seed),
            "device": str(args.device),
            "operator": multires.operator_contract(backend),
            "calibration_parameters": multires.CALIBRATION,
        },
        "scientific_scope": (
            "base boundary only; no held-out outcomes read; operator selection only"),
        "held_out_profile_data_read": False,
        "steps": [],
    }
    try:
        geometry = _initial_geometry(args.dx_nm)
        bootstrap = multires._advance(
            geometry, None, None, duration_s=0.0, seed=int(args.seed),
            device=args.device, topology_policy="refuse", remap_backend=backend)
        state = bootstrap.next_surface_state
        fingerprint = bootstrap.next_surface_state_mesh_fingerprint
        geometry = bootstrap.geometry
        payload["initial_geometry_sha256"] = _geometry_digest(geometry)
        payload["initial_state"] = _state_receipt(
            state, bootstrap.next_active_face_area)
        for index in range(int(args.steps)):
            if perf_counter() - started >= float(args.max_wall_s):
                payload["status"] = "wall_budget_checkpoint"
                break
            step_started = perf_counter()
            result = multires._advance(
                geometry, state, fingerprint,
                duration_s=float(args.step_duration_s),
                seed=int(args.seed) + index,
                device=args.device,
                topology_policy="refuse",
                remap_backend=backend,
            )
            record = {
                "step": index + 1,
                "physical_time_s": (index + 1) * float(args.step_duration_s),
                "wall_time_s": perf_counter() - step_started,
                "geometry_sha256": _geometry_digest(result.geometry),
                "metrics": pilot.measure_krueger_metrics(
                    result.geometry, substrate_top_um=1.8),
                "operator": multires._operator_summary(result),
                "next_state": _state_receipt(
                    result.next_surface_state, result.next_active_face_area),
                "remap": _remap_receipt(result.state_remap_diagnostics),
                "maximum_remap_relative_conservation_residual": (
                    _maximum_remap_residual(result.state_remap_diagnostics)),
                "topology_event": result.diagnostics.get("topology_event"),
            }
            payload["steps"].append(record)
            geometry = result.geometry
            state = result.next_surface_state
            fingerprint = result.next_surface_state_mesh_fingerprint
        if len(payload["steps"]) == int(args.steps):
            payload["status"] = "complete"
    except Exception as error:  # refusal details are part of this operator audit
        payload["status"] = "refused"
        payload["exception"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
    payload["wall_time_s"] = perf_counter() - started
    _write_json(output / "audit.json", payload)
    return payload


def _comparison(cases, requested_backends):
    completed = {
        name: case for name, case in cases.items()
        if case.get("status") == "complete"
    }
    first_geometry = {
        case["steps"][0]["geometry_sha256"] for case in completed.values()
        if case.get("steps")
    }
    initial_geometry = {
        case.get("initial_geometry_sha256") for case in completed.values()
    }
    initial_state = {
        case.get("initial_state", {}).get("sha256") for case in completed.values()
    }
    backend_gates = {}
    for backend in requested_backends:
        case = cases.get(backend, {})
        steps = case.get("steps", ())
        residual_values = [
            float(step["maximum_remap_relative_conservation_residual"])
            for step in steps]
        ledger_values = [
            float(step["operator"]["maximum_material_ledger_residual_units_m2"])
            for step in steps]
        residual = max(residual_values) if residual_values else None
        ledger = max(ledger_values) if ledger_values else None
        backend_gates[backend] = {
            "complete": case.get("status") == "complete",
            "step_count": len(steps),
            "topology_clean": bool(steps) and all(
                step.get("topology_event") is None for step in steps),
            "maximum_remap_relative_conservation_residual": residual,
            "maximum_material_ledger_residual_units_m2": ledger,
            "conservation_pass": bool(
                residual is not None and ledger is not None
                and residual <= 1.0e-12 and ledger <= 1.0e-20),
        }
    indexed = backend_gates.get("indexed_knn", {})
    common = backend_gates.get("common_refinement", {})
    common_pass = bool(
        indexed.get("complete")
        and indexed.get("step_count") == 2
        and indexed.get("topology_clean")
        and indexed.get("conservation_pass")
        and common.get("complete")
        and common.get("step_count") == 2
        and common.get("topology_clean")
        and common.get("conservation_pass")
        and len(initial_geometry) == 1
        and len(initial_state) == 1
        and len(first_geometry) == 1
    )
    return {
        "identical_initial_geometry": len(initial_geometry) == 1,
        "identical_initial_state": len(initial_state) == 1,
        "identical_first_step_geometry": len(first_geometry) == 1,
        "backend_gates": backend_gates,
        "common_refinement_candidate_pass": common_pass,
        "selection": (
            "common_refinement_candidate_for_bounded_5nm_confirmation"
            if common_pass else "blocked; inspect operator refusals before calibration"),
    }


def run(args):
    if args.worker_backend:
        return _worker(args)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    cases = {}
    for backend in args.backends:
        case_output = output / backend
        audit_path = case_output / "audit.json"
        if args.reuse_worker_audits and audit_path.is_file():
            case = json.loads(audit_path.read_text(encoding="utf-8"))
            case["worker"] = {"reused_existing_audit": True}
            cases[backend] = case
            continue
        command = build_worker_command(args, backend, case_output)
        environment = dict(os.environ)
        environment["PETCH_DEVICE"] = str(args.device)
        environment["OMP_NUM_THREADS"] = "1"
        try:
            completed = subprocess.run(
                command, cwd=ROOT, env=environment, capture_output=True,
                text=True, timeout=float(args.max_wall_s), check=False)
            if audit_path.is_file():
                case = json.loads(audit_path.read_text(encoding="utf-8"))
            else:
                case = {
                    "status": "worker_failed_without_audit",
                    "returncode": completed.returncode,
                }
            case["worker"] = {
                "returncode": completed.returncode,
                "stdout_tail": completed.stdout[-2000:],
                "stderr_tail": completed.stderr[-2000:],
            }
        except subprocess.TimeoutExpired as error:
            case = {
                "status": "wall_timeout",
                "timeout_s": float(args.max_wall_s),
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
        cases[backend] = case
    comparison = _comparison(cases, args.backends)
    payload = {
        "schema": "petch.krueger_2024_remap_backend_audit.v1",
        "status": (
            "candidate_selected"
            if comparison["common_refinement_candidate_pass"]
            else "operator_selection_blocked"),
        "configuration": {
            "dx_nm": float(args.dx_nm),
            "steps": int(args.steps),
            "step_duration_s": float(args.step_duration_s),
            "total_physical_time_s": int(args.steps) * float(args.step_duration_s),
            "seed": int(args.seed),
            "device": str(args.device),
            "backends": list(args.backends),
        },
        "scientific_scope": (
            "bounded base-boundary operator selection; no experimental outcomes read"),
        "held_out_profile_data_read": False,
        "cases": cases,
        "comparison": comparison,
    }
    _write_json(output / "audit.json", payload)
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dx-nm", type=float, default=10.0)
    parser.add_argument("--steps", type=int, default=2)
    parser.add_argument("--step-duration-s", type=float, default=0.025)
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-wall-s", type=float, default=300.0)
    parser.add_argument(
        "--backends", nargs="+", choices=multires.REMAP_BACKENDS,
        default=list(multires.REMAP_BACKENDS))
    parser.add_argument(
        "--reuse-worker-audits", action="store_true",
        help="rebuild only the aggregate receipt from existing backend audits")
    parser.add_argument(
        "--worker-backend", choices=multires.REMAP_BACKENDS,
        help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.steps <= 0 or args.steps > MAXIMUM_STEPS:
        parser.error(f"steps must be in [1, {MAXIMUM_STEPS}]")
    if args.step_duration_s <= 0.0:
        parser.error("step-duration-s must be positive")
    if args.steps * args.step_duration_s > MAXIMUM_PHYSICAL_TIME_S + 1.0e-15:
        parser.error(
            f"bounded audit refuses more than {MAXIMUM_PHYSICAL_TIME_S:g} s")
    if args.dx_nm not in (5.0, 10.0):
        parser.error("bounded audit supports only endpoint-compatible 5 or 10 nm grids")
    if args.max_wall_s <= 0.0 or args.max_wall_s > 600.0:
        parser.error("max-wall-s must be in (0, 600]")
    return args


if __name__ == "__main__":
    run(parse_args())
