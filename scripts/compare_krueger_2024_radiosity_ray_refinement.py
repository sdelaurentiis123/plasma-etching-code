#!/usr/bin/env python3
"""Compare two nested-ray Krueger frozen radiosity/chemistry receipts.

This is a measurement tool, not a convergence gate.  It requires two completed
audits of the same checkpoint, direct transport, chemistry parameters, and
physical horizon, with distinct nested form-factor ray counts.  It then reports
the observed cross-ray differences without applying a post-hoc tolerance.

Every persisted final surface-state field, per-face exchange inventory,
per-face recession/growth displacement, integrated exchange, and oxide-removal
scalar is compared for R17/R19 and nominal/tight chemistry integration.  Array
errors use symmetric relative Linf and, when a compatible face-area vector is
persisted by both inputs, physical-area-weighted relative L1.  Otherwise the
receipt explicitly labels the second norm as normalized (unweighted) L1.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


SCHEMA = "petch.krueger-2024.radiosity-ray-refinement-comparison.v1"
EXPECTED_AUDIT_SCHEMA = "petch.krueger-2024.frozen-radiosity-chemistry.v1"
PARAMETER_LABELS = ("r17", "r19")
INTEGRATION_LABELS = ("nominal", "tight")
FACE_AREA_PATHS = (
    ("face_area_m2",),
    ("active_face_area_m2",),
    ("full_face_area_m2",),
    ("surface_mesh", "active_face_area_m2"),
    ("surface_mesh", "face_area_m2"),
    ("direct_transport", "active_face_area_m2"),
    ("form_factors", "face_area_m2"),
)


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_atomic(path, payload):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)


def _load_completed_audit(path):
    source = Path(path)
    payload = json.loads(source.read_text(encoding="utf-8"))
    if payload.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError(
            f"{source.name} has unsupported schema {payload.get('schema')!r}")
    if payload.get("status") != "pass":
        raise ValueError(
            f"{source.name} is not a completed passing audit: "
            f"status={payload.get('status')!r}")
    return payload


def _get_path(mapping, path):
    current = mapping
    for name in path:
        if not isinstance(current, dict) or name not in current:
            raise KeyError(".".join(path))
        current = current[name]
    return current


def _optional_path(mapping, path):
    try:
        return _get_path(mapping, path)
    except KeyError:
        return None


def _require_equal(left, right, path):
    left_value = _get_path(left, path)
    right_value = _get_path(right, path)
    if left_value != right_value:
        dotted = ".".join(path)
        raise ValueError(f"incompatible audit provenance at {dotted}")


def _is_power_of_two(value):
    return isinstance(value, int) and value > 0 and value & (value - 1) == 0


def _ray_count(audit):
    actual = int(_get_path(audit, ("form_factors", "rays_per_face")))
    requested = int(_get_path(
        audit, ("form_factors", "requested_rays_per_face")))
    runtime_requested = int(_get_path(
        audit, ("provenance", "runtime_selection", "requested_rays_per_face")))
    if requested != runtime_requested:
        raise ValueError("form-factor and runtime requested ray counts disagree")
    if not _is_power_of_two(actual) or not _is_power_of_two(requested):
        raise ValueError("actual and requested ray counts must be powers of two")
    if actual < requested or actual % requested:
        raise ValueError("actual ray count is not a nested extension of its request")
    return actual


def _validate_audit_pair(left, right):
    """Validate every invariant that must be shared by a ray refinement pair."""
    exact_paths = (
        ("schema",),
        ("scientific_scope",),
        ("data_firewall",),
        ("checkpoint", "audit_sha256"),
        ("checkpoint", "checkpoint_sha256"),
        ("checkpoint", "metadata"),
        ("parameter_provenance",),
        ("transport_operator",),
        ("direct_transport", "direct_surface_flux_sha256"),
        ("direct_transport", "boundary_provenance"),
        ("form_factors", "face_count"),
        ("execution_budget", "next_profile_step_s"),
        ("provenance", "source"),
        ("provenance", "base_inputs"),
    )
    for path in exact_paths:
        _require_equal(left, right, path)

    ignored_radiosity_keys = {"rays_per_face", "maximum_rays_per_face"}
    left_radiosity = {
        key: value for key, value in left["radiosity_operator"].items()
        if key not in ignored_radiosity_keys
    }
    right_radiosity = {
        key: value for key, value in right["radiosity_operator"].items()
        if key not in ignored_radiosity_keys
    }
    if left_radiosity != right_radiosity:
        raise ValueError("radiosity operators differ beyond their ray budgets")

    left_rays = _ray_count(left)
    right_rays = _ray_count(right)
    if left_rays == right_rays:
        raise ValueError("ray-refinement inputs use the same actual ray count")
    lower, higher = sorted((left_rays, right_rays))
    ratio = higher // lower
    if higher % lower or not _is_power_of_two(ratio):
        raise ValueError(
            "actual ray counts are not members of one nested power-of-two sequence")
    return {
        "lower_actual_rays_per_face": lower,
        "higher_actual_rays_per_face": higher,
        "nested_extension_factor": ratio,
        "shared_checkpoint_sha256": left["checkpoint"]["checkpoint_sha256"],
        "shared_direct_surface_flux_sha256": left["direct_transport"][
            "direct_surface_flux_sha256"],
        "shared_face_count": int(left["form_factors"]["face_count"]),
    }


def _passing_horizons(audit):
    output = {}
    for entry in audit.get("horizons", []):
        if not entry.get("common_pass", False):
            continue
        fraction = float(entry["fraction_of_next_profile_step"])
        if fraction in output:
            raise ValueError(f"duplicate completed horizon fraction {fraction}")
        for parameter in PARAMETER_LABELS:
            result = entry.get("parameter_results", {}).get(parameter)
            if not isinstance(result, dict) or not result.get("all_gates_pass", False):
                raise ValueError(
                    f"horizon {fraction} lacks a complete {parameter} result")
        output[fraction] = entry
    if not output:
        raise ValueError("audit has no completed common-passing horizon")
    return output


def _select_horizon(left, right, requested_fraction=None):
    left_entries = _passing_horizons(left)
    right_entries = _passing_horizons(right)
    common = sorted(set(left_entries) & set(right_entries))
    if not common:
        raise ValueError("audits do not share a completed physical horizon")
    if requested_fraction is None:
        fraction = common[-1]
        rule = "largest_common_completed_horizon"
    else:
        fraction = float(requested_fraction)
        if fraction not in common:
            raise ValueError(
                f"requested horizon fraction {fraction} is not complete in both audits")
        rule = "explicit_horizon_fraction"
    left_entry = left_entries[fraction]
    right_entry = right_entries[fraction]
    if left_entry["horizon_s"] != right_entry["horizon_s"]:
        raise ValueError("matching horizon fractions have different physical durations")
    return left_entry, right_entry, {
        "selection_rule": rule,
        "fraction_of_next_profile_step": fraction,
        "horizon_s": float(left_entry["horizon_s"]),
        "all_common_completed_fractions": common,
    }


def _candidate_face_area(audit, face_count):
    found = []
    for path in FACE_AREA_PATHS:
        value = _optional_path(audit, path)
        if value is None:
            continue
        area = np.asarray(value, dtype=float)
        if area.shape != (face_count,):
            raise ValueError(
                f"persisted face area at {'.'.join(path)} has shape {area.shape}, "
                f"expected {(face_count,)}")
        if not np.all(np.isfinite(area)) or np.any(area <= 0.0):
            raise ValueError("persisted face area must be finite and strictly positive")
        found.append((".".join(path), area))
    if not found:
        return None, None
    reference_path, reference = found[0]
    for path, candidate in found[1:]:
        if not np.array_equal(reference, candidate):
            raise ValueError(
                f"audit persists inconsistent face areas at {reference_path} and {path}")
    return reference_path, reference


def _resolve_face_area(left, right, face_count):
    left_path, left_area = _candidate_face_area(left, face_count)
    right_path, right_area = _candidate_face_area(right, face_count)
    if (left_area is None) != (right_area is None):
        raise ValueError("only one audit persists a face-area vector")
    if left_area is None:
        return None, {
            "kind": "normalized_unweighted_relative_l1",
            "reason": "no compatible per-face area vector is persisted in either audit",
        }
    if not np.array_equal(left_area, right_area):
        raise ValueError("persisted face-area vectors differ across ray levels")
    return left_area, {
        "kind": "physical_area_weighted_relative_l1",
        "lower_ray_source": left_path,
        "higher_ray_source": right_path,
    }


def _finite_array(values, name):
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _array_error(lower, higher, *, face_area, name):
    lower = _finite_array(lower, f"lower-ray {name}")
    higher = _finite_array(higher, f"higher-ray {name}")
    if lower.shape != higher.shape:
        raise ValueError(
            f"{name} shape differs across ray levels: {lower.shape} != {higher.shape}")
    difference = np.abs(higher - lower)
    absolute_linf = float(np.max(difference)) if difference.size else 0.0
    linf_scale = max(
        float(np.max(np.abs(lower))) if lower.size else 0.0,
        float(np.max(np.abs(higher))) if higher.size else 0.0,
    )
    relative_linf = absolute_linf / linf_scale if linf_scale > 0.0 else 0.0

    use_area = (
        face_area is not None
        and lower.ndim == 1
        and lower.shape == face_area.shape
    )
    weights = face_area if use_area else np.ones(lower.shape, dtype=float)
    absolute_l1 = float(np.sum(weights * difference))
    l1_scale = max(
        float(np.sum(weights * np.abs(lower))),
        float(np.sum(weights * np.abs(higher))),
    )
    relative_l1 = absolute_l1 / l1_scale if l1_scale > 0.0 else 0.0
    return {
        "shape": list(lower.shape),
        "absolute_linf_error": absolute_linf,
        "symmetric_relative_linf_error": relative_linf,
        "absolute_l1_error": absolute_l1,
        "symmetric_relative_l1_error": relative_l1,
        "l1_norm": (
            "physical_area_weighted_relative_l1"
            if use_area else "normalized_unweighted_relative_l1"),
    }


def _flatten_arrays(mapping, prefix=""):
    if not isinstance(mapping, dict):
        raise ValueError(f"{prefix or 'array collection'} must be a mapping")
    output = {}
    for key, value in sorted(mapping.items()):
        name = f"{prefix}/{key}" if prefix else str(key)
        if isinstance(value, dict):
            output.update(_flatten_arrays(value, name))
        else:
            output[name] = value
    return output


def _compare_array_collection(lower, higher, *, face_area, label):
    lower_flat = _flatten_arrays(lower)
    higher_flat = _flatten_arrays(higher)
    if set(lower_flat) != set(higher_flat):
        missing_lower = sorted(set(higher_flat) - set(lower_flat))
        missing_higher = sorted(set(lower_flat) - set(higher_flat))
        raise ValueError(
            f"{label} keys differ; missing lower={missing_lower}, "
            f"missing higher={missing_higher}")
    by_item = {
        name: _array_error(
            lower_flat[name], higher_flat[name], face_area=face_area,
            name=f"{label}/{name}")
        for name in sorted(lower_flat)
    }
    return {
        "item_count": len(by_item),
        "maximum_symmetric_relative_linf_error": max(
            (item["symmetric_relative_linf_error"] for item in by_item.values()),
            default=0.0),
        "maximum_symmetric_relative_l1_error": max(
            (item["symmetric_relative_l1_error"] for item in by_item.values()),
            default=0.0),
        "by_item": by_item,
    }


def _flatten_scalars(mapping, prefix=""):
    if not isinstance(mapping, dict):
        raise ValueError(f"{prefix or 'scalar collection'} must be a mapping")
    output = {}
    for key, value in sorted(mapping.items()):
        name = f"{prefix}/{key}" if prefix else str(key)
        if isinstance(value, dict):
            output.update(_flatten_scalars(value, name))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            output[name] = float(value)
        else:
            raise ValueError(f"{name} is not a numeric scalar")
    return output


def _scalar_error(lower, higher, name):
    lower = float(lower)
    higher = float(higher)
    if not np.isfinite(lower) or not np.isfinite(higher):
        raise ValueError(f"{name} contains a non-finite scalar")
    absolute = abs(higher - lower)
    scale = max(abs(lower), abs(higher))
    return {
        "lower_ray_value": lower,
        "higher_ray_value": higher,
        "absolute_error": absolute,
        "symmetric_relative_error": absolute / scale if scale > 0.0 else 0.0,
    }


def _compare_scalar_collection(lower, higher, label):
    lower_flat = _flatten_scalars(lower)
    higher_flat = _flatten_scalars(higher)
    if set(lower_flat) != set(higher_flat):
        raise ValueError(f"{label} scalar keys differ across ray levels")
    by_item = {
        name: _scalar_error(lower_flat[name], higher_flat[name], f"{label}/{name}")
        for name in sorted(lower_flat)
    }
    return {
        "item_count": len(by_item),
        "maximum_symmetric_relative_error": max(
            (item["symmetric_relative_error"] for item in by_item.values()),
            default=0.0),
        "by_item": by_item,
    }


def _compare_integration(lower, higher, *, face_area, label):
    required = (
        "final_state_fields",
        "final_state_fields_sha256",
        "per_face_integrated_exchange_units_m2",
        "per_face_integrated_exchange_sha256",
        "integrated_exchange",
        "oxide_removal",
        "displacement",
    )
    for key in required:
        if key not in lower or key not in higher:
            raise ValueError(f"{label} lacks required persisted field {key}")
    displacement_keys = (
        "per_face_integrated_recession_m",
        "per_face_integrated_growth_m",
    )
    displacement_lower = {key: lower["displacement"][key] for key in displacement_keys}
    displacement_higher = {key: higher["displacement"][key] for key in displacement_keys}
    return {
        "exact_hashes": {
            "lower_final_state_fields_sha256": lower[
                "final_state_fields_sha256"],
            "higher_final_state_fields_sha256": higher[
                "final_state_fields_sha256"],
            "final_state_fields_hash_equal": bool(
                lower["final_state_fields_sha256"]
                == higher["final_state_fields_sha256"]),
            "lower_per_face_exchange_sha256": lower[
                "per_face_integrated_exchange_sha256"],
            "higher_per_face_exchange_sha256": higher[
                "per_face_integrated_exchange_sha256"],
            "per_face_exchange_hash_equal": bool(
                lower["per_face_integrated_exchange_sha256"]
                == higher["per_face_integrated_exchange_sha256"]),
        },
        "final_state_fields": _compare_array_collection(
            lower["final_state_fields"], higher["final_state_fields"],
            face_area=face_area, label=f"{label}/final_state_fields"),
        "per_face_integrated_exchange_units_m2": _compare_array_collection(
            lower["per_face_integrated_exchange_units_m2"],
            higher["per_face_integrated_exchange_units_m2"],
            face_area=face_area,
            label=f"{label}/per_face_integrated_exchange_units_m2"),
        "per_face_displacement": _compare_array_collection(
            displacement_lower, displacement_higher, face_area=face_area,
            label=f"{label}/displacement"),
        "integrated_exchange": _compare_scalar_collection(
            lower["integrated_exchange"], higher["integrated_exchange"],
            f"{label}/integrated_exchange"),
        "oxide_removal": _compare_scalar_collection(
            lower["oxide_removal"], higher["oxide_removal"],
            f"{label}/oxide_removal"),
    }


def _paired_direction(entry, integration_label):
    r17 = float(entry["parameter_results"]["r17"][integration_label][
        "oxide_removal"]["integrated_formula_units"])
    r19 = float(entry["parameter_results"]["r19"][integration_label][
        "oxide_removal"]["integrated_formula_units"])
    difference = r19 - r17
    direction = (
        "r19_lower" if difference < 0.0
        else "r19_higher" if difference > 0.0 else "equal")
    return {
        "r19_minus_r17_integrated_formula_units": difference,
        "r19_to_r17_ratio": r19 / r17 if r17 > 0.0 else None,
        "direction": direction,
    }


def _validate_persisted_tight_direction(entry):
    persisted = entry.get("paired_oxide_removal_direction")
    if not isinstance(persisted, dict):
        raise ValueError("horizon lacks its persisted tight paired direction")
    computed = _paired_direction(entry, "tight")
    for key in ("r19_minus_r17_integrated_formula_units", "r19_to_r17_ratio"):
        if persisted.get(key) != computed[key]:
            raise ValueError(f"persisted paired direction disagrees at {key}")
    if persisted.get("direction") != computed["direction"]:
        raise ValueError("persisted paired direction label disagrees with tight results")


def compare_audits(audit_a, audit_b, *, source_a=None, source_b=None,
                   horizon_fraction=None):
    compatibility = _validate_audit_pair(audit_a, audit_b)
    rays_a = _ray_count(audit_a)
    rays_b = _ray_count(audit_b)
    if rays_a < rays_b:
        lower_audit, higher_audit = audit_a, audit_b
        lower_source, higher_source = source_a, source_b
    else:
        lower_audit, higher_audit = audit_b, audit_a
        lower_source, higher_source = source_b, source_a
    lower_entry, higher_entry, selected = _select_horizon(
        lower_audit, higher_audit, horizon_fraction)
    _validate_persisted_tight_direction(lower_entry)
    _validate_persisted_tight_direction(higher_entry)

    face_count = compatibility["shared_face_count"]
    face_area, norm_provenance = _resolve_face_area(
        lower_audit, higher_audit, face_count)
    comparisons = {}
    global_array_linf = 0.0
    global_array_l1 = 0.0
    global_scalar = 0.0
    for parameter in PARAMETER_LABELS:
        comparisons[parameter] = {}
        for integration in INTEGRATION_LABELS:
            label = f"{parameter}/{integration}"
            receipt = _compare_integration(
                lower_entry["parameter_results"][parameter][integration],
                higher_entry["parameter_results"][parameter][integration],
                face_area=face_area, label=label)
            comparisons[parameter][integration] = receipt
            for category in (
                    "final_state_fields",
                    "per_face_integrated_exchange_units_m2",
                    "per_face_displacement"):
                global_array_linf = max(
                    global_array_linf,
                    receipt[category]["maximum_symmetric_relative_linf_error"])
                global_array_l1 = max(
                    global_array_l1,
                    receipt[category]["maximum_symmetric_relative_l1_error"])
            for category in ("integrated_exchange", "oxide_removal"):
                global_scalar = max(
                    global_scalar,
                    receipt[category]["maximum_symmetric_relative_error"])

    paired = {}
    for integration in INTEGRATION_LABELS:
        lower = _paired_direction(lower_entry, integration)
        higher = _paired_direction(higher_entry, integration)
        paired[integration] = {
            "lower_ray": lower,
            "higher_ray": higher,
            "direction_preserved": bool(lower["direction"] == higher["direction"]),
            "difference_comparison": _scalar_error(
                lower["r19_minus_r17_integrated_formula_units"],
                higher["r19_minus_r17_integrated_formula_units"],
                f"paired/{integration}/r19_minus_r17"),
        }

    def source_record(path, audit, rays):
        record = {
            "actual_rays_per_face": rays,
            "audit_payload_schema": audit["schema"],
        }
        if path is not None:
            path = Path(path)
            record.update({
                "input_name": f"{path.parent.name}/{path.name}",
                "input_sha256": _sha256(path),
            })
        return record

    return {
        "schema": SCHEMA,
        "status": "complete",
        "scientific_scope": (
            "bounded cross-ray measurement on one frozen checkpoint and one "
            "identical physical chemistry horizon"),
        "decision_contract": {
            "post_hoc_pass_tolerance_applied": False,
            "interpretation": (
                "This receipt reports observed refinement errors and paired direction only; "
                "it does not declare ray convergence."),
        },
        "sources": {
            "lower_ray": source_record(
                lower_source, lower_audit, compatibility[
                    "lower_actual_rays_per_face"]),
            "higher_ray": source_record(
                higher_source, higher_audit, compatibility[
                    "higher_actual_rays_per_face"]),
        },
        "compatibility": compatibility,
        "selected_horizon": selected,
        "array_norm_contract": {
            "relative_linf": (
                "symmetric: max(abs(high-low)) / "
                "max(max(abs(low)), max(abs(high))); both-zero maps to zero"),
            "relative_l1": norm_provenance,
            "per_array_override": (
                "arrays not shaped as the persisted face-area vector use normalized "
                "unweighted relative L1 and say so in their l1_norm field"),
        },
        "comparisons": comparisons,
        "paired_r19_minus_r17_oxide_removal": paired,
        "observed_global_maxima": {
            "array_symmetric_relative_linf_error": global_array_linf,
            "array_symmetric_relative_l1_error": global_array_l1,
            "scalar_symmetric_relative_error": global_scalar,
        },
    }


def run(args):
    path_a = Path(args.audit_a)
    path_b = Path(args.audit_b)
    audit_a = _load_completed_audit(path_a)
    audit_b = _load_completed_audit(path_b)
    receipt = compare_audits(
        audit_a, audit_b, source_a=path_a, source_b=path_b,
        horizon_fraction=args.horizon_fraction)
    _write_json_atomic(args.output, receipt)
    print(json.dumps({
        "status": receipt["status"],
        "selected_horizon": receipt["selected_horizon"],
        "observed_global_maxima": receipt["observed_global_maxima"],
        "paired_r19_minus_r17_oxide_removal": receipt[
            "paired_r19_minus_r17_oxide_removal"],
        "output": str(args.output),
    }, indent=2, sort_keys=True))
    return receipt


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-a", type=Path, required=True)
    parser.add_argument("--audit-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--horizon-fraction", type=float)
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
