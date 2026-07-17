#!/usr/bin/env python3
"""Bounded 2x2 frozen-checkpoint diagnostic for the Krueger base case.

This script evaluates two completed profile checkpoints (R17 and R19) with two
mechanism parameter pairs (R17 and R19).  Every cell of the matrix uses the
current production transport/chemistry operator with ``duration_s=0``.  The
geometry and accumulated surface state are therefore frozen:

====================  ====================  ====================
checkpoint state      R17 parameters        R19 parameters
====================  ====================  ====================
R17 geometry/state    instantaneous rate    instantaneous rate
R19 geometry/state    instantaneous rate    instantaneous rate
====================  ====================  ====================

The row contrasts measure accumulated geometry/state effects at fixed
parameters.  The column contrasts measure instantaneous parameter effects at
fixed geometry/state.  A difference-in-differences reports their interaction.
No profile evolution, fitting, transfer-boundary construction, or held-out
observation loading occurs here.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from krueger_2024_endpoint_operator_audit import (  # noqa: E402
    _evaluate,
    _surface_flux_by_species,
)
from krueger_2024_trench_pilot import (  # noqa: E402
    _load_checkpoint,
    _maximum_ledger_residual,
    measure_krueger_metrics,
)


SCHEMA = "petch.krueger-2024.frozen-checkpoint-2x2.v1"
PARAMETER_FIELDS = (
    "effective_mask_crosslinked_growth_fraction",
    "oxide_etch_yield_scale",
)
MATERIAL_NAMES = {1: "sio2", 2: "amorphous_carbon_mask"}
CURRENT_OPERATOR = {
    "boundary_case": "base",
    "boundary_mode": "angular_8x16",
    "ion_energy_bin_eV": 250.0,
    "ion_angle_bin_deg": 0.25,
    "ion_azimuthal_closure": "axisymmetric_uniform",
    "ion_azimuthal_order": 16,
    "ballistic_transport": "face_gather",
    "ballistic_face_quadrature_points": 3,
    "n_position": 16,
    "neutral_radiosity_rays_per_face": 8,
    "transport_device": "cpu",
    "duration_s": 0.0,
}
DIAGNOSTIC_LOW_OPERATOR = {
    **CURRENT_OPERATOR,
    "boundary_mode": "legacy_compressed_tensor",
    "ion_energy_bin_eV": 500.0,
    "ion_angle_bin_deg": 0.5,
    "ballistic_face_quadrature_points": 1,
    "fidelity": "diagnostic_low_screen",
    "authority": False,
}
DIAGNOSTIC_Q3_OPERATOR = {
    **DIAGNOSTIC_LOW_OPERATOR,
    "ballistic_face_quadrature_points": 3,
    "fidelity": "diagnostic_q3_screen",
}
CURRENT_OPERATOR.update(
    fidelity="current_production",
    authority=True,
)
FIDELITY_OPTIONS = {
    "current_production": CURRENT_OPERATOR,
    "diagnostic_low_screen": DIAGNOSTIC_LOW_OPERATOR,
    "diagnostic_q3_screen": DIAGNOSTIC_Q3_OPERATOR,
}
GATES = {
    "maximum_material_ledger_residual_units_m2": 0.0,
    "maximum_radiosity_relative_balance_error": 5.0e-12,
}
SOURCE_PATHS = (
    "scripts/krueger_2024_frozen_checkpoint_2x2.py",
    "scripts/krueger_2024_endpoint_operator_audit.py",
    "scripts/krueger_2024_trench_pilot.py",
    "src/petch/amorphous_carbon_mask.py",
    "src/petch/boundary_transport_3d.py",
    "src/petch/experimental_data.py",
    "src/petch/feature_step_3d.py",
    "src/petch/material_mechanism_3d.py",
    "src/petch/neutral_radiosity_3d.py",
    "src/petch/reactor_boundary.py",
)
BASE_INPUT_PATHS = (
    "data/experimental/krueger_2024/base_case_boundary_fluxes.csv",
    "data/experimental/krueger_2024/digitized_figure4_iead.csv",
    "data/experimental/krueger_2024/digitized_figure4_iead_metadata.json",
)
CELL_NAMES = tuple(
    f"{checkpoint}_checkpoint__{parameters}_parameters"
    for checkpoint in ("r17", "r19")
    for parameters in ("r17", "r19")
)


class EvaluationDeadlineExceeded(TimeoutError):
    """Raised when a native transport evaluation crosses its hard wall budget."""


@contextmanager
def _hard_deadline(seconds):
    """Interrupt Python or an interruptible native Warp launch at a hard deadline."""
    seconds = float(seconds)
    if seconds <= 0.0:
        raise EvaluationDeadlineExceeded("evaluation has no wall-time budget remaining")
    if not hasattr(signal, "SIGALRM"):
        raise RuntimeError("hard evaluation deadlines require SIGALRM on this platform")

    def expired(_signum, _frame):
        raise EvaluationDeadlineExceeded(
            f"zero-duration operator evaluation exceeded {seconds:.3f} s")

    previous_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, expired)
    previous_timer = signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_timer[0] > 0.0:
            signal.setitimer(signal.ITIMER_REAL, *previous_timer)


def _jsonable(value):
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return str(value)


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _hash_manifest(relative_paths):
    entries = {}
    aggregate = sha256()
    for relative in sorted(relative_paths):
        path = ROOT / relative
        digest = _sha256(path)
        entries[relative] = digest
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
    return {"aggregate_sha256": aggregate.hexdigest(), "files": entries}


def _write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _maximum_radiosity_error(result):
    return max(
        (
            float(item["relative_balance_error"])
            for item in result.diagnostics.get("neutral_radiosity", {}).values()
        ),
        default=0.0,
    )


def _integrated_flux(flux, active, area_m2, material):
    active_flux = np.asarray(flux, dtype=float)[active]
    output = {
        "total": {
            "incident_rate_s-1": float(np.sum(active_flux * area_m2)),
            "area_weighted_mean_flux_m-2_s-1": float(
                np.sum(active_flux * area_m2) / max(np.sum(area_m2), np.finfo(float).tiny)
            ),
        },
        "by_material": {},
    }
    for material_id in sorted(set(int(value) for value in material)):
        selected = material == material_id
        selected_area = area_m2[selected]
        name = MATERIAL_NAMES.get(material_id, f"material_{material_id}")
        output["by_material"][name] = {
            "material_id": material_id,
            "incident_rate_s-1": float(
                np.sum(active_flux[selected] * selected_area)
            ),
            "area_weighted_mean_flux_m-2_s-1": float(
                np.sum(active_flux[selected] * selected_area)
                / max(np.sum(selected_area), np.finfo(float).tiny)
            ),
        }
    return output


def _velocity_statistics(velocity_m_s, area_m2, material):
    def summarize(selected):
        values = velocity_m_s[selected]
        area = area_m2[selected]
        total_area = float(np.sum(area))
        tiny = np.finfo(float).tiny
        return {
            "surface_area_m2": total_area,
            "net_signed_volume_rate_m3_s": float(np.sum(values * area)),
            "gross_absolute_volume_rate_m3_s": float(np.sum(np.abs(values) * area)),
            "area_weighted_mean_velocity_m_s": float(
                np.sum(values * area) / max(total_area, tiny)
            ),
            "area_weighted_mean_absolute_velocity_m_s": float(
                np.sum(np.abs(values) * area) / max(total_area, tiny)
            ),
            "area_weighted_rms_velocity_m_s": float(
                np.sqrt(np.sum(values ** 2 * area) / max(total_area, tiny))
            ),
            "minimum_velocity_m_s": float(np.min(values)) if values.size else 0.0,
            "maximum_velocity_m_s": float(np.max(values)) if values.size else 0.0,
            "active_face_count": int(values.size),
        }

    output = {"total": summarize(np.ones(len(material), dtype=bool)), "by_material": {}}
    for material_id in sorted(set(int(value) for value in material)):
        name = MATERIAL_NAMES.get(material_id, f"material_{material_id}")
        output["by_material"][name] = {
            "material_id": material_id,
            **summarize(material == material_id),
        }
    return output


def summarize_result(result, geometry):
    """Reduce one mesh-dependent operator result to comparable physical integrals."""
    active = np.asarray(result.active_face_index, dtype=int)
    area_mesh2 = np.asarray(result.active_face_area, dtype=float)
    if area_mesh2.shape != active.shape:
        raise ValueError("active face area does not match active face index")
    length_unit = float(geometry.mesh_length_unit_m)
    area_m2 = area_mesh2 * length_unit ** 2
    material = np.asarray(result.face_material_id, dtype=int)[active]
    velocity_m_s = (
        np.asarray(result.face_velocity_mesh_units_s, dtype=float)[active] * length_unit
    )
    fluxes = {
        name: _integrated_flux(values, active, area_m2, material)
        for name, values in sorted(_surface_flux_by_species(result).items())
    }
    ledger = float(_maximum_ledger_residual(result.surface.material_exchange))
    radiosity = _maximum_radiosity_error(result)
    validity = result.validity
    gates = {
        "material_ledger_exact": bool(
            ledger == GATES["maximum_material_ledger_residual_units_m2"]
        ),
        "radiosity_balance_within_tolerance": bool(
            radiosity <= GATES["maximum_radiosity_relative_balance_error"]
        ),
        "within_declared_scope": bool(validity.within_declared_scope),
    }
    return {
        "velocity": _velocity_statistics(velocity_m_s, area_m2, material),
        "incident_flux_by_species": fluxes,
        "conservation_and_validity": {
            "maximum_material_ledger_residual_units_m2": ledger,
            "maximum_radiosity_relative_balance_error": radiosity,
            "validity_reasons": list(validity.reasons),
            "gates": gates,
            "all_gates_pass": bool(all(gates.values())),
        },
        "boundary_provenance": _jsonable(result.transport.boundary.provenance)
        if hasattr(result.transport, "boundary")
        else None,
    }


def _snapshot_inputs(geometry, state):
    return {
        "phi": np.asarray(geometry.phi).copy(),
        "material_id": np.asarray(geometry.material_id).copy(),
        "material_levelsets": {
            int(key): np.asarray(value).copy()
            for key, value in dict(geometry.material_levelsets or {}).items()
        },
        "state": {
            str(key): np.asarray(value).copy() for key, value in state.fields.items()
        },
    }


def _inputs_unchanged(snapshot, geometry, state):
    return bool(
        np.array_equal(snapshot["phi"], geometry.phi)
        and np.array_equal(snapshot["material_id"], geometry.material_id)
        and set(snapshot["material_levelsets"]) == set(geometry.material_levelsets or {})
        and all(
            np.array_equal(value, geometry.material_levelsets[key])
            for key, value in snapshot["material_levelsets"].items()
        )
        and set(snapshot["state"]) == set(state.fields)
        and all(
            np.array_equal(value, state.fields[key])
            for key, value in snapshot["state"].items()
        )
    )


def _load_source(label, source):
    source = Path(source)
    audit_path = source / "audit.json"
    checkpoint_path = source / "checkpoint.npz"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "complete":
        raise RuntimeError(f"{label} source is not a completed run")
    config = audit["configuration"]
    if config.get("boundary_case") != "base":
        raise ValueError(f"{label} source is not the sealed Krueger base boundary")
    geometry, state, fingerprint, metadata = _load_checkpoint(checkpoint_path)
    geometry_config = config.get("geometry", {})
    substrate_top = float(geometry_config.get("substrate_top_um", 1.8))
    opening_width = float(geometry_config.get("opening_width_um", 0.09))
    cell_width = float(geometry_config.get("cell_width_um", 0.13))
    metrics = measure_krueger_metrics(
        geometry,
        substrate_top_um=substrate_top,
        opening_center_um=0.5 * cell_width,
        opening_width_um=opening_width,
    )
    parameters = {name: float(config[name]) for name in PARAMETER_FIELDS}
    return {
        "label": label,
        "source": source,
        "audit_path": audit_path,
        "checkpoint_path": checkpoint_path,
        "audit": audit,
        "config": config,
        "geometry": geometry,
        "state": state,
        "fingerprint": fingerprint,
        "checkpoint_metadata": metadata,
        "parameters": parameters,
        "metrics": metrics,
    }


def _operator_config(reference_config, parameters):
    config = deepcopy(reference_config)
    config.update(parameters)
    config.update(
        neutral_direction_polar_order=8,
        neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=250.0,
        ion_angle_bin_deg=0.25,
        ion_azimuthal_closure="axisymmetric_uniform",
        ion_azimuthal_order=16,
        ballistic_transport="face_gather",
        ballistic_face_quadrature_points=3,
        n_position=16,
        radiosity_rays_per_face=8,
        radiosity_enabled=True,
        transport_device="cpu",
        boundary_case="base",
        duration_s=0.0,
    )
    return config


def _flatten_comparable(summary):
    output = {}
    for scope, values in [("total", summary["velocity"]["total"]), *(
            (name, values)
            for name, values in summary["velocity"]["by_material"].items())]:
        for field in (
            "net_signed_volume_rate_m3_s",
            "gross_absolute_volume_rate_m3_s",
            "area_weighted_mean_velocity_m_s",
            "area_weighted_mean_absolute_velocity_m_s",
        ):
            output[f"velocity.{scope}.{field}"] = float(values[field])
    for species, values in summary["incident_flux_by_species"].items():
        output[f"flux.{species}.total.incident_rate_s-1"] = float(
            values["total"]["incident_rate_s-1"]
        )
        for material, material_values in values["by_material"].items():
            output[f"flux.{species}.{material}.incident_rate_s-1"] = float(
                material_values["incident_rate_s-1"]
            )
    return output


def _difference(left, right):
    """Return left minus right for all mesh-independent scalar summaries."""
    left = _flatten_comparable(left)
    right = _flatten_comparable(right)
    if set(left) != set(right):
        raise RuntimeError("2x2 cells expose different comparable scalar fields")
    output = {}
    for name in sorted(left):
        difference = left[name] - right[name]
        output[name] = {
            "difference": difference,
            "fraction_of_right": (
                difference / right[name] if right[name] != 0.0 else None
            ),
        }
    return output


def _interaction(a, b, c, d):
    """Return (a-b)-(c-d), the parameter-by-checkpoint interaction."""
    values = [_flatten_comparable(item) for item in (a, b, c, d)]
    if any(set(item) != set(values[0]) for item in values[1:]):
        raise RuntimeError("2x2 cells expose different comparable scalar fields")
    return {
        name: float((values[0][name] - values[1][name])
                    - (values[2][name] - values[3][name]))
        for name in sorted(values[0])
    }


def _write_report(path, payload):
    cells = payload["evaluations"]
    lines = [
        "# Krueger R17/R19 frozen-checkpoint 2x2 diagnostic",
        "",
        "All four entries are current-operator, base-boundary, zero-duration evaluations. "
        "They diagnose instantaneous rates; they do not evolve or validate a profile.",
        "",
        "| Frozen checkpoint | Mechanism pair | SiO2 net volume rate (m3/s) | "
        "Mask net volume rate (m3/s) | Wall time (s) | Gates |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for checkpoint in ("r17", "r19"):
        for parameters in ("r17", "r19"):
            name = f"{checkpoint}_checkpoint__{parameters}_parameters"
            if name not in cells:
                continue
            cell = cells[name]
            velocity = cell["summary"]["velocity"]["by_material"]
            gates = cell["summary"]["conservation_and_validity"]["all_gates_pass"]
            lines.append(
                f"| {checkpoint.upper()} | {parameters.upper()} | "
                f"{velocity['sio2']['net_signed_volume_rate_m3_s']:.9e} | "
                f"{velocity['amorphous_carbon_mask']['net_signed_volume_rate_m3_s']:.9e} | "
                f"{cell['wall_time_s']:.3f} | {'pass' if gates else 'fail'} |"
            )
    lines.extend([
        "",
        "## Frozen profile observables",
        "",
        "| Checkpoint | Opening (nm) | Depth (nm) | Physical time (s) |",
        "| --- | ---: | ---: | ---: |",
    ])
    for checkpoint in ("r17", "r19"):
        item = payload["checkpoints"][checkpoint]
        lines.append(
            f"| {checkpoint.upper()} | {item['metrics']['mask_opening_nm']:.6f} | "
            f"{item['metrics']['etch_depth_nm']:.6f} | "
            f"{item['checkpoint_metadata']['physical_time_s']:.6f} |"
        )
    lines.extend([
        "",
        "## Interpretation contract",
        "",
        "- Column differences at one checkpoint are instantaneous parameter effects.",
        "- Row differences at one parameter pair are accumulated geometry/state effects.",
        "- The interaction is the change in parameter sensitivity between the two checkpoints.",
        "- None of these rates is an endpoint prediction or held-out validation result.",
        "",
        f"Overall diagnostic status: **{payload['status']}**.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args):
    started = perf_counter()
    sources = {
        "r17": _load_source("r17", args.r17_source),
        "r19": _load_source("r19", args.r19_source),
    }
    reference_config = sources["r19"]["config"]
    requested = list(CELL_NAMES)
    selected = getattr(args, "cell", "all")
    if selected != "all":
        if selected not in CELL_NAMES:
            raise ValueError(f"unknown 2x2 cell {selected!r}")
        requested = [selected]
    evaluations = {}
    boundary_hashes = set()
    timeout = None
    gate_failure = None
    fidelity_name = getattr(args, "fidelity", "current_production")
    operator = dict(FIDELITY_OPTIONS[fidelity_name])
    operator["transport_seed"] = int(args.seed)
    operator["neutral_radiosity_seed"] = int(args.seed) + 10000
    for checkpoint_label in ("r17", "r19"):
        source = sources[checkpoint_label]
        snapshot = _snapshot_inputs(source["geometry"], source["state"])
        for parameter_label in ("r17", "r19"):
            cell_name = f"{checkpoint_label}_checkpoint__{parameter_label}_parameters"
            if cell_name not in requested:
                continue
            config = _operator_config(
                reference_config, sources[parameter_label]["parameters"]
            )
            cell_started = perf_counter()
            remaining_total = (
                float(getattr(args, "maximum_total_wall_s", 480.0))
                - (perf_counter() - started)
            )
            deadline = min(float(args.maximum_evaluation_wall_s), remaining_total)
            try:
                with _hard_deadline(deadline):
                    result, boundary, evaluator_wall = _evaluate(
                        source["geometry"], source["state"], source["fingerprint"],
                        boundary_mode=operator["boundary_mode"],
                        ion_bins=(
                            operator["ion_energy_bin_eV"],
                            operator["ion_angle_bin_deg"],
                        ),
                        face_points=operator["ballistic_face_quadrature_points"],
                        pilot_config=config,
                        radiosity_rays=operator[
                            "neutral_radiosity_rays_per_face"
                        ],
                        seed=int(args.seed),
                        ballistic_transport=operator["ballistic_transport"],
                        n_position=operator["n_position"],
                        transport_device="cpu",
                    )
            except EvaluationDeadlineExceeded as error:
                timeout = {
                    "cell": cell_name,
                    "hard_deadline_s": float(deadline),
                    "elapsed_wall_time_s": float(perf_counter() - cell_started),
                    "reason": str(error),
                }
                break
            wall = perf_counter() - cell_started
            if not _inputs_unchanged(snapshot, source["geometry"], source["state"]):
                raise RuntimeError("zero-duration evaluation mutated checkpoint inputs")
            boundary_provenance = _jsonable(boundary.provenance)
            boundary_hashes.add(json.dumps(boundary_provenance, sort_keys=True))
            summary = summarize_result(result, source["geometry"])
            summary["boundary_provenance"] = boundary_provenance
            evaluations[cell_name] = {
                "checkpoint": checkpoint_label,
                "mechanism_parameters": parameter_label,
                "parameter_values": sources[parameter_label]["parameters"],
                "transport_seed": int(args.seed),
                "neutral_radiosity_seed": int(args.seed) + 10000,
                "wall_time_s": float(wall),
                "reported_evaluator_wall_time_s": float(evaluator_wall),
                "input_checkpoint_unchanged": True,
                "summary": summary,
            }
            _write_json_atomic(Path(args.output), {
                "schema": SCHEMA,
                "status": "running_partial",
                "scientific_scope": "zero-duration diagnostic screen only",
                "current_operator": operator,
                "execution_budget": {
                    "requested_cells": requested,
                    "completed_cells": list(evaluations),
                    "transport_seed": int(args.seed),
                    "neutral_radiosity_seed": int(args.seed) + 10000,
                    "maximum_evaluation_wall_s": float(
                        args.maximum_evaluation_wall_s
                    ),
                },
                "evaluations": evaluations,
            })
            if not summary["conservation_and_validity"]["all_gates_pass"]:
                gate_failure = {
                    "cell": cell_name,
                    "gates": summary["conservation_and_validity"]["gates"],
                    "reason": "matrix stopped on the first numerical gate failure",
                }
                break
        if timeout is not None or gate_failure is not None:
            break
    if boundary_hashes and len(boundary_hashes) != 1:
        raise RuntimeError("2x2 cells did not use one identical base-boundary operator")

    all_gates = all(
        cell["summary"]["conservation_and_validity"]["all_gates_pass"]
        for cell in evaluations.values()
    )
    complete_matrix = all(name in evaluations for name in CELL_NAMES)
    if timeout is not None:
        status = "bounded_timeout"
    elif gate_failure is not None:
        status = "gate_failure"
    elif complete_matrix:
        status = "pass" if all_gates else "fail"
    else:
        status = "single_cell_complete" if all_gates else "single_cell_gate_failure"
    payload = {
        "schema": SCHEMA,
        "status": status,
        "scientific_scope": (
            "instantaneous current-operator decomposition on frozen completed base-case "
            "checkpoints; no profile evolution and no validation claim"
        ),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
            "base_loader_contract": "load_krueger_2024_base_boundary_fluxes only",
        },
        "current_operator": operator,
        "fidelity_contract": {
            "classification": (
                "authority" if operator["authority"] else "diagnostic_screen_only"
            ),
            "may_drive_physics_or_calibration": bool(operator["authority"]),
            "promotion_gate": (
                "Before any low-fidelity contrast can drive physics or calibration, at least "
                "one selected contrast must have the same direction under paired q3 "
                "current-production evaluations on both relevant matrix cells."
            ),
        },
        "gates": GATES,
        "execution_budget": {
            "requested_cells": requested,
            "transport_seed": int(args.seed),
            "neutral_radiosity_seed": int(args.seed) + 10000,
            "maximum_evaluation_wall_s": float(args.maximum_evaluation_wall_s),
            "maximum_total_wall_s": float(
                getattr(args, "maximum_total_wall_s", 480.0)
            ),
            "timeout": timeout,
            "gate_failure": gate_failure,
        },
        "checkpoints": {
            label: {
                "metrics": _jsonable(source["metrics"]),
                "parameters": source["parameters"],
                "checkpoint_metadata": _jsonable(source["checkpoint_metadata"]),
                "audit_sha256": _sha256(source["audit_path"]),
                "checkpoint_sha256": _sha256(source["checkpoint_path"]),
                "archived_config_hash": source["audit"].get("config_hash"),
            }
            for label, source in sources.items()
        },
        "evaluations": evaluations,
        "contrasts": {},
        "provenance": {
            "source": _hash_manifest(SOURCE_PATHS),
            "base_inputs": _hash_manifest(BASE_INPUT_PATHS),
            "identical_boundary_provenance_sha256": (
                sha256(next(iter(boundary_hashes)).encode("utf-8")).hexdigest()
                if boundary_hashes else None
            ),
        },
        "total_wall_time_s": float(perf_counter() - started),
    }
    if complete_matrix:
        rr = evaluations["r17_checkpoint__r17_parameters"]["summary"]
        rn = evaluations["r17_checkpoint__r19_parameters"]["summary"]
        nr = evaluations["r19_checkpoint__r17_parameters"]["summary"]
        nn = evaluations["r19_checkpoint__r19_parameters"]["summary"]
        payload["contrasts"] = {
            "instantaneous_parameter_effect_at_r17_checkpoint": _difference(rn, rr),
            "instantaneous_parameter_effect_at_r19_checkpoint": _difference(nn, nr),
            "accumulated_checkpoint_effect_at_r17_parameters": _difference(nr, rr),
            "accumulated_checkpoint_effect_at_r19_parameters": _difference(nn, rn),
            "parameter_by_checkpoint_interaction": _interaction(nn, nr, rn, rr),
        }
    destination = Path(args.output)
    _write_json_atomic(destination, payload)
    if args.report is not None:
        _write_report(Path(args.report), payload)
    print(json.dumps({
        "status": payload["status"],
        "output": str(destination),
        "report": None if args.report is None else str(args.report),
        "total_wall_time_s": payload["total_wall_time_s"],
    }, indent=2, sort_keys=True))
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--r17-source",
        default=ROOT / "results" / "krueger_2024_base_calibration_r17" / "axisym_candidate",
    )
    parser.add_argument(
        "--r19-source",
        default=ROOT / "results" / "krueger_2024_r19_response_check" / "remote_artifacts",
    )
    output = (
        ROOT / "results" / "krueger_2024_r19_response_check"
        / "frozen_checkpoint_2x2" / "audit.json"
    )
    parser.add_argument("--output", default=output)
    parser.add_argument("--report", default=output.with_name("report.md"))
    parser.add_argument("--seed", type=int, default=241)
    parser.add_argument("--maximum-evaluation-wall-s", type=float, default=120.0)
    parser.add_argument("--maximum-total-wall-s", type=float, default=480.0)
    parser.add_argument("--cell", choices=("all",) + CELL_NAMES, default="all")
    parser.add_argument(
        "--fidelity", choices=tuple(FIDELITY_OPTIONS), default="current_production"
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
