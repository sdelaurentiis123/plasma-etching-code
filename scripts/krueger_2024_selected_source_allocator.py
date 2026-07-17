#!/usr/bin/env python3
"""Audit-only selected-source allocator for the Krueger Stage-A RQMC hold.

This script consumes the completed 8/16-ray Stage-A JSON/NPZ artifacts.  It performs no
transport tracing, chemistry integration, profile motion, fitting, or held-out-data read.

For the linear radiosity equation in particle-rate coordinates,

    q = d + E(F) diag(1-s) q,

the fine-minus-coarse change in any linear patch observable ``w.T @ q`` has the exact
source-row decomposition

    delta y_j = (1-s_j) q_coarse,j [delta E[:,j].T lambda_fine],

where ``A_fine.T @ lambda_fine = w``.  The signed row contributions sum to the directly
computed nested change to roundoff.  We use those exact contributions only to rank rows;
we do not pretend that the rowwise absolute-value ranking is itself an uncertainty bound.

The output proposes appending Sobol points 16:32 only on selected source faces, with the
same per-face shift and replicate seed.  The current production estimator accepts one global
ray count, so this script emits a plan and explicitly records that the row-selective tracer API
is a required, separately reviewed implementation step.  It never launches that refinement.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT / "src", ROOT / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import krueger_2024_replicated_form_factor_closure as closure  # noqa: E402
from petch.amorphous_carbon_mask import (  # noqa: E402
    build_krueger_2024_material_router_3d,
)
from petch.neutral_radiosity_3d import DiffuseFormFactors3D  # noqa: E402
from petch.surface_patch_convergence_3d import (  # noqa: E402
    aggregate_surface_field_on_physical_patches_3d,
)


SCHEMA = "petch.krueger-2024.selected-source-form-factor-allocator.v1"
BASE_RAY_LEVEL = 16
SELECTED_RAY_LEVEL = 32
TARGET_CONCENTRATION = 0.90
MAXIMUM_SELECTED_FRACTION = 0.25
MAXIMUM_WALL_S = 120.0
HORIZON_DIVISOR = 1024


def _sha256(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload):
    return sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode("utf-8")).hexdigest()


def _write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _exchange_matrix(factors):
    if not isinstance(factors, DiffuseFormFactors3D):
        raise TypeError("factors must be DiffuseFormFactors3D")
    return sparse.coo_matrix(
        (factors.transfer_fraction,
         (factors.target_face, factors.source_face)),
        shape=(factors.face_count, factors.face_count)).tocsr()


def exact_radiosity_source_row_decomposition(
        coarse_factors, fine_factors, face_area_m2, direct_flux_m2_s,
        reaction_probability, objective_rate_weights):
    """Return an exact signed source-row decomposition for nested radiosity changes.

    ``objective_rate_weights`` has shape ``(objective, face)`` and acts on incident
    particle rate ``q = area * incident_flux``.  The returned contribution array has
    shape ``(objective, source_face)``.
    """
    if (not isinstance(coarse_factors, DiffuseFormFactors3D)
            or not isinstance(fine_factors, DiffuseFormFactors3D)
            or coarse_factors.face_count != fine_factors.face_count):
        raise ValueError("coarse and fine form factors must share one face mesh")
    n_face = coarse_factors.face_count
    area = np.asarray(face_area_m2, dtype=float)
    direct = np.asarray(direct_flux_m2_s, dtype=float)
    reaction = np.asarray(reaction_probability, dtype=float)
    weights = np.asarray(objective_rate_weights, dtype=float)
    if weights.ndim == 1:
        weights = weights[None, :]
    if (area.shape != (n_face,) or direct.shape != (n_face,)
            or reaction.shape != (n_face,) or weights.ndim != 2
            or weights.shape[1] != n_face or np.any(~np.isfinite(area))
            or np.any(area <= 0.0) or np.any(~np.isfinite(direct))
            or np.any(direct < 0.0) or np.any(~np.isfinite(reaction))
            or np.any((reaction < 0.0) | (reaction > 1.0))
            or np.any(~np.isfinite(weights))):
        raise ValueError("invalid radiosity decomposition inputs")

    reflection = 1.0 - reaction
    coarse_exchange = _exchange_matrix(coarse_factors)
    fine_exchange = _exchange_matrix(fine_factors)
    coarse_operator = (
        sparse.eye(n_face, format="csc")
        - (coarse_exchange @ sparse.diags(reflection)).tocsc())
    fine_operator = (
        sparse.eye(n_face, format="csc")
        - (fine_exchange @ sparse.diags(reflection)).tocsc())
    direct_rate = area * direct
    try:
        coarse_lu = splu(coarse_operator)
        fine_lu = splu(fine_operator)
        coarse_rate = coarse_lu.solve(direct_rate)
        fine_rate = fine_lu.solve(direct_rate)
        adjoint = fine_lu.solve(weights.T, trans="T")
    except RuntimeError as error:
        raise RuntimeError(
            "selected-source allocator requires a nonsingular radiosity operator") from error
    for name, operator, rate in (
            ("coarse", coarse_operator, coarse_rate),
            ("fine", fine_operator, fine_rate)):
        scale = max(float(np.linalg.norm(direct_rate)), np.finfo(float).tiny)
        residual = float(np.linalg.norm(operator @ rate - direct_rate) / scale)
        negative = 1.0e-11 * max(float(np.max(rate, initial=0.0)), 1.0)
        if (not np.isfinite(residual) or residual > 2.0e-10
                or np.any(~np.isfinite(rate)) or np.any(rate < -negative)):
            raise RuntimeError(
                f"{name} radiosity solve is not a certified nonnegative solution")
    delta_exchange = fine_exchange - coarse_exchange
    projected = np.asarray(delta_exchange.T @ adjoint)
    contributions = (
        projected * (reflection * coarse_rate)[:, None]).T
    direct_difference = weights @ (fine_rate - coarse_rate)
    decomposed_difference = np.sum(contributions, axis=1)
    closure_error = np.abs(decomposed_difference - direct_difference)
    closure_scale = np.maximum.reduce((
        np.abs(direct_difference),
        np.sum(np.abs(contributions), axis=1),
        np.full(len(direct_difference), np.finfo(float).tiny)))
    relative_closure = closure_error / closure_scale
    if float(np.max(relative_closure, initial=0.0)) > 2.0e-9:
        raise RuntimeError("exact source-row decomposition failed its closure check")
    return {
        "contributions": contributions,
        "direct_difference": direct_difference,
        "decomposed_difference": decomposed_difference,
        "relative_closure_error": relative_closure,
        "coarse_incident_rate": coarse_rate,
        "fine_incident_rate": fine_rate,
    }


def rank_source_rows(
        authority_contributions, replicate_contributions, objective_scale, *,
        target_concentration=TARGET_CONCENTRATION,
        maximum_selected_fraction=MAXIMUM_SELECTED_FRACTION):
    """Rank source rows using authority magnitude plus replicate standard error.

    The per-row score is an allocation heuristic.  Exact signed decomposition receipts
    remain separate so the score is never misreported as an additive confidence bound.
    """
    authority = np.asarray(authority_contributions, dtype=float)
    replicates = np.asarray(replicate_contributions, dtype=float)
    scale = np.asarray(objective_scale, dtype=float)
    if (authority.ndim != 2 or replicates.ndim != 3
            or replicates.shape[1:] != authority.shape
            or scale.shape != (authority.shape[0],)
            or len(replicates) < 4 or np.any(~np.isfinite(authority))
            or np.any(~np.isfinite(replicates)) or np.any(~np.isfinite(scale))
            or np.any(scale <= 0.0)):
        raise ValueError("invalid source-row ranking inputs")
    target = float(target_concentration)
    cap_fraction = float(maximum_selected_fraction)
    if not 0.0 < target <= 1.0 or not 0.0 < cap_fraction <= 1.0:
        raise ValueError("invalid concentration or selected-fraction contract")
    critical = float(student_t.ppf(0.975, len(replicates) - 1))
    stochastic = critical * np.std(
        replicates, axis=0, ddof=1) / np.sqrt(len(replicates))
    per_objective = (np.abs(authority) + stochastic) / scale[:, None]
    score = np.sum(per_objective, axis=0)
    order = np.argsort(-score, kind="stable")
    total = float(np.sum(score))
    cumulative = (
        np.cumsum(score[order]) / total if total > 0.0
        else np.zeros_like(score))
    required = (
        int(np.searchsorted(cumulative, target, side="left") + 1)
        if total > 0.0 else 0)
    cap = max(1, int(np.floor(cap_fraction * len(score))))
    selected_count = min(required, cap)
    selected = order[:selected_count]
    milestones = sorted(set(
        [1, 5, 10, 25, 50, 100, selected_count, required, cap, len(score)]))
    curve = [{
        "source_face_count": int(count),
        "source_face_fraction": float(count / len(score)),
        "cumulative_ranking_score_fraction": float(cumulative[count - 1]),
    } for count in milestones if 0 < count <= len(score)]
    return {
        "row_score": score,
        "ranked_source_faces": order,
        "selected_source_faces": selected,
        "target_concentration": target,
        "maximum_selected_fraction": cap_fraction,
        "required_face_count_for_target": required,
        "selected_face_count": selected_count,
        "selected_score_fraction": (
            float(cumulative[selected_count - 1]) if selected_count else 0.0),
        "target_reached_within_cap": bool(required <= cap),
        "concentration_curve": curve,
    }


def _validate_stage_a_contract(audit, audit_path):
    if (audit.get("schema") != closure.SCHEMA
            or audit.get("stage") != "stage_a"
            or audit.get("status") != "bounded_precision_hold"):
        raise ValueError("allocator requires the completed bounded Stage-A hold")
    firewall = audit.get("data_firewall", {})
    if (firewall.get("boundary_case") != "base"
            or firewall.get("held_out_observations_loaded") is not False
            or firewall.get("held_out_transfer_boundary_constructed") is not False):
        raise ValueError("Stage-A artifact does not preserve the base-case firewall")
    if (audit.get("operator") != closure.OPERATOR
            or tuple(audit.get("sampling", {}).get("ray_levels", ()))
            != closure.RAY_LEVELS
            or tuple(audit.get("sampling", {}).get("replicate_seeds", ()))
            != closure.REPLICATE_SEEDS
            or audit.get("stage_a", {}).get("all_gates_pass") is not False):
        raise ValueError("Stage-A operator, nesting, seed, or failed-gate identity changed")
    required = {
        "exact_nested_sampling_extension", "exact_replicate_count", "row_closure"}
    gates = audit.get("stage_a", {}).get("gates", {})
    if not required.issubset(gates) or not all(gates[name] for name in required):
        raise ValueError("Stage-A sampling authority did not pass its exact gates")
    return {
        "audit_path_name": Path(audit_path).name,
        "audit_sha256": _sha256(audit_path),
        "checkpoint_sha256": audit["checkpoint"]["checkpoint_sha256"],
    }


def _patch_rate_weights(direct, source, patch_scale_m, normalization, reaction=None):
    scheme = aggregate_surface_field_on_physical_patches_3d(
        np.zeros(len(direct["face_area_m2"])), direct["face_area_m2"],
        direct["verts"], direct["faces"], direct["gas_normals"],
        direct["face_material_id"], float(patch_scale_m),
        mesh_length_unit_m=source["geometry"].mesh_length_unit_m)
    multiplier = np.ones(len(direct["face_area_m2"]))
    if reaction is not None:
        multiplier = np.asarray(reaction, dtype=float)
    values = (
        scheme.contribution_area_m2
        / scheme.patch_area_m2[scheme.contribution_patch_index]
        / direct["face_area_m2"][scheme.contribution_face_index]
        / float(normalization)
        * multiplier[scheme.contribution_face_index])
    weights = sparse.coo_matrix(
        (values, (scheme.contribution_patch_index,
                  scheme.contribution_face_index)),
        shape=(len(scheme.patch_key), len(direct["face_area_m2"]))).toarray()
    return scheme, weights


def _paired_patch_gate(coarse, fine, replicates_coarse, replicates_fine):
    paired = np.asarray(replicates_fine) - np.asarray(replicates_coarse)
    critical = float(student_t.ppf(0.975, len(paired) - 1))
    half = critical * np.std(paired, axis=0, ddof=1) / np.sqrt(len(paired))
    difference = np.asarray(fine) - np.asarray(coarse)
    scale = (
        closure.GATES["patch_absolute_normalized"]
        + closure.GATES["patch_relative_tolerance"]
        * np.maximum(np.abs(coarse), np.abs(fine)))
    score = (np.abs(difference) + half) / scale
    return difference, half, scale, score


def _dominant_axis_projected_support_fraction(
        scheme, face_gas_normals, *, patch_scale_m, periodic_y_extent_m):
    """Projected patch support, with the periodic-y footprint counted only once.

    A 40 nm Cartesian patch applied to this benchmark's 20 nm periodic y-cell does not
    contain two independent copies of the surface.  Whenever y is tangential to the
    patch's dominant normal, its nominal projected footprint is therefore
    ``patch_scale * min(patch_scale, periodic_y_extent)``.
    """
    normal = np.asarray(face_gas_normals, dtype=float)
    scale = float(patch_scale_m)
    periodic_y = float(periodic_y_extent_m)
    if (normal.ndim != 2 or normal.shape[1] != 3
            or not np.isfinite(scale) or scale <= 0.0
            or not np.isfinite(periodic_y) or periodic_y <= 0.0):
        raise ValueError("invalid projected-support inputs")
    contribution_axis = scheme.patch_key[
        scheme.contribution_patch_index, 1]
    projected = np.bincount(
        scheme.contribution_patch_index,
        weights=(scheme.contribution_area_m2 * np.abs(normal[
            scheme.contribution_face_index, contribution_axis])),
        minlength=len(scheme.patch_key))
    dominant_axis = scheme.patch_key[:, 1]
    expected = np.full(len(dominant_axis), scale * scale)
    y_is_tangential = dominant_axis != 1
    expected[y_is_tangential] *= min(scale, periodic_y) / scale
    return projected / expected


def _factor_pairs(coarse, fine):
    yield "mean", coarse.mean_form_factors, fine.mean_form_factors
    for seed, left, right in zip(
            coarse.replicate_seeds, coarse.replicate_form_factors,
            fine.replicate_form_factors):
        yield f"seed:{int(seed)}", left, right


def run(args):
    started = perf_counter()
    audit_path = Path(args.stage_a_audit)
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    authority = _validate_stage_a_contract(audit, audit_path)
    source = closure._load_sealed_base_source(args.r19_source)
    if _sha256(source["checkpoint_path"]) != authority["checkpoint_sha256"]:
        raise ValueError("allocator checkpoint differs from Stage A")

    config = closure._production_config(source["config"], closure.PARAMETERS["r19"])
    direct_identity = closure._direct_transport_cache_identity(
        source, config, seed=int(audit["direct_transport"]["transport_seed"]))
    direct_path = audit_path.parent / "direct_transport_cache.npz"
    direct, direct_meta = closure._read_direct_transport_artifact(
        direct_path, source, identity=direct_identity)
    ensembles = {}
    factor_identity = {}
    for level in closure.RAY_LEVELS:
        record = audit["form_factor_ensembles"][str(level)]["artifact"]
        path = audit_path.parent / record["npz_relative_to_audit"]
        ensemble, metadata = closure._read_form_factor_ensemble(
            path, checkpoint_sha256=authority["checkpoint_sha256"],
            ray_level=level, replicate_seeds=closure.REPLICATE_SEEDS,
            expected_face_area_m2=direct["face_area_m2"])
        if (record["npz_sha256"] != metadata["npz_sha256"]
                or record["metadata_sha256"] != _sha256(path.with_suffix(".json"))):
            raise ValueError("Stage-A audit does not bind a form-factor artifact")
        ensembles[level] = ensemble
        factor_identity[str(level)] = {
            "npz_path_name": path.name,
            "npz_sha256": metadata["npz_sha256"],
            "metadata_sha256": _sha256(path.with_suffix(".json")),
            "ensemble_sha256": ensemble.sha256,
        }

    mechanism = build_krueger_2024_material_router_3d(
        **closure.PARAMETERS["r19"])
    probabilities = mechanism.neutral_reaction_probability_by_material(
        source["state"], direct["face_material_id"])
    species = tuple(sorted(direct["direct_surface_fluxes"].neutral_flux_m2_s))
    patch_weights = {}
    for name in species:
        normalization = float(np.max(np.abs(
            direct["direct_surface_fluxes"].neutral_flux_m2_s[name]), initial=0.0))
        for scale in closure.PATCH_SCALES_M:
            scheme, incident = _patch_rate_weights(
                direct, source, scale, normalization)
            _same, reacted = _patch_rate_weights(
                direct, source, scale, normalization,
                reaction=probabilities[name])
            patch_weights[(name, scale, "incident")] = (scheme, incident)
            patch_weights[(name, scale, "reacted")] = (scheme, reacted)

    modes = {
        "all_patches": {
            "objective_records": [], "authority": [], "replicates": [], "scales": []},
        "physically_supported_patches": {
            "objective_records": [], "authority": [], "replicates": [], "scales": []},
    }
    periodic_y_extent_m = (
        (source["geometry"].phi.shape[1] - 1) * source["geometry"].dx
        * source["geometry"].mesh_length_unit_m)
    maximum_decomposition_error = 0.0
    for name in species:
        if perf_counter() - started > float(args.maximum_wall_s):
            raise TimeoutError("selected-source analysis exceeded its bounded wall time")
        direct_flux = direct["direct_surface_fluxes"].neutral_flux_m2_s[name]
        pair_results = {}
        objective_blocks = []
        for kind in ("incident", "reacted"):
            for scale in closure.PATCH_SCALES_M:
                scheme, weights = patch_weights[(name, scale, kind)]
                objective_blocks.append((kind, scale, scheme, weights))
        all_weights = np.concatenate([item[3] for item in objective_blocks], axis=0)
        for member, coarse, fine in _factor_pairs(ensembles[8], ensembles[16]):
            pair_results[member] = exact_radiosity_source_row_decomposition(
                coarse, fine, direct["face_area_m2"], direct_flux,
                probabilities[name], all_weights)
        offset = 0
        for kind, scale, scheme, weights in objective_blocks:
            count = len(weights)
            authority_result = pair_results["mean"]
            coarse_authority = (
                weights @ authority_result["coarse_incident_rate"])
            fine_authority = weights @ authority_result["fine_incident_rate"]
            replicate_coarse = []
            replicate_fine = []
            for seed in closure.REPLICATE_SEEDS:
                result = pair_results[f"seed:{seed}"]
                replicate_coarse.append(weights @ result["coarse_incident_rate"])
                replicate_fine.append(weights @ result["fine_incident_rate"])
            difference, half, gate_scale, score = _paired_patch_gate(
                coarse_authority, fine_authority,
                replicate_coarse, replicate_fine)
            support = _dominant_axis_projected_support_fraction(
                scheme, direct["gas_normals"], patch_scale_m=scale,
                periodic_y_extent_m=periodic_y_extent_m)
            supported = np.flatnonzero(support >= 0.10)
            chosen = {"all_patches": int(np.argmax(score))}
            if len(supported):
                chosen["physically_supported_patches"] = int(
                    supported[np.argmax(score[supported])])
            for mode, worst in chosen.items():
                row = offset + worst
                replicate_rows = np.stack([
                    pair_results[f"seed:{seed}"]["contributions"][row]
                    for seed in closure.REPLICATE_SEEDS])
                authority_row = authority_result["contributions"][row]
                maximum_decomposition_error = max(
                    maximum_decomposition_error,
                    float(authority_result["relative_closure_error"][row]),
                    max(float(pair_results[f"seed:{seed}"][
                        "relative_closure_error"][row])
                        for seed in closure.REPLICATE_SEEDS))
                if float(score[worst]) <= 1.0:
                    continue
                modes[mode]["authority"].append(authority_row)
                modes[mode]["replicates"].append(replicate_rows)
                modes[mode]["scales"].append(float(gate_scale[worst]))
                modes[mode]["objective_records"].append({
                    "species": name,
                    "field": kind,
                    "patch_scale_m": float(scale),
                    "patch_key": [int(value) for value in scheme.patch_key[worst]],
                    "patch_area_m2": float(scheme.patch_area_m2[worst]),
                    "paired_authority_difference_normalized": float(difference[worst]),
                    "paired_95pct_confidence_half_width_normalized": float(half[worst]),
                    "mixed_gate_scale_normalized": float(gate_scale[worst]),
                    "combined_mixed_normalized": float(score[worst]),
                    "dominant_axis_projected_support_fraction": float(support[worst]),
                    "minimum_supported_fraction": (
                        0.10 if mode == "physically_supported_patches" else None),
                    "failed": True,
                })
            offset += count
    if any(not value["authority"] for value in modes.values()):
        raise RuntimeError("no failed radiosity patch objective remains to allocate")
    rankings = {}
    for mode, values in modes.items():
        rankings[mode] = rank_source_rows(
            np.stack(values["authority"]), np.stack(values["replicates"], axis=1),
            np.asarray(values["scales"]),
            target_concentration=args.target_concentration,
            maximum_selected_fraction=args.maximum_selected_fraction)
    ranking = rankings["physically_supported_patches"]
    selected = ranking["selected_source_faces"]
    ranked = ranking["ranked_source_faces"]
    face_count = len(direct["face_area_m2"])
    extra_selected_rays = int(len(selected) * (SELECTED_RAY_LEVEL - BASE_RAY_LEVEL))
    extra_global_rays = int(face_count * (SELECTED_RAY_LEVEL - BASE_RAY_LEVEL))
    next_step_s = float(audit["checkpoint"]["next_profile_step_s"])
    stage_a_peak_dx = max(
        float(audit["stage_a"]["levels"][str(level)]["parameter_scores"][label][
            "maximum_gross_displacement_dx"])
        for level in closure.RAY_LEVELS for label in ("r17", "r19"))
    horizon_s = next_step_s / HORIZON_DIVISOR
    predicted_dx = stage_a_peak_dx * (
        horizon_s / float(audit["checkpoint"]["shortest_horizon_s"]))
    top_count = min(50, len(ranked))
    row_score = ranking["row_score"]
    allocation_comparison = {}
    for mode, local in rankings.items():
        allocation_comparison[mode] = {
            "objective_count": len(modes[mode]["objective_records"]),
            "required_face_count_for_target": local[
                "required_face_count_for_target"],
            "selected_face_count": local["selected_face_count"],
            "selected_face_fraction": float(local["selected_face_count"] / face_count),
            "selected_score_fraction": local["selected_score_fraction"],
            "target_reached_within_cap": local["target_reached_within_cap"],
            "concentration_curve": local["concentration_curve"],
            "selected_source_faces": [
                int(value) for value in local["selected_source_faces"]],
        }
    payload = {
        "schema": SCHEMA,
        "status": (
            "selected_source_plan_ready"
            if ranking["target_reached_within_cap"]
            else "diffuse_source_error_blocker"),
        "scientific_scope": (
            "audit-only exact linear-radiosity row attribution; no new ray trace, "
            "chemistry integration, profile motion, fitting, or held-out data"),
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
        },
        "authority": {
            **authority,
            "operator": closure.OPERATOR,
            "direct_transport_npz_sha256": direct_meta["npz_sha256"],
            "direct_transport_metadata_sha256": _sha256(
                direct_path.with_suffix(".json")),
            "form_factor_artifacts": factor_identity,
            "replicate_seeds": list(closure.REPLICATE_SEEDS),
            "nested_ray_levels": list(closure.RAY_LEVELS),
        },
        "decomposition": {
            "coordinate_system": "incident_particle_rate_q_equals_area_times_flux",
            "identity": "delta_q=A_fine_inverse_delta_E_diag_reflection_q_coarse",
            "objective_count": len(modes[
                "physically_supported_patches"]["objective_records"]),
            "patch_support_contract": {
                "minimum_dominant_axis_projected_support_fraction": 0.10,
                "periodic_y_extent_m": periodic_y_extent_m,
                "forty_nm_y_footprint_uses_one_twenty_nm_periodic_cell": True,
            },
            "worst_failed_patch_per_species_field_and_scale": modes[
                "physically_supported_patches"]["objective_records"],
            "all_patch_worst_objectives_for_diagnostic_comparison": modes[
                "all_patches"]["objective_records"],
            "maximum_relative_signed_closure_error": maximum_decomposition_error,
            "film_attribution_scope": (
                "film fields are nonlinear downstream consumers and were not reintegrated; "
                "the selected union covers their failed upstream incident/reacted neutral fields"),
        },
        "allocation": {
            "executable_plan_basis": "physically_supported_patches",
            "ranking_metric": (
                "sum over failed objectives of (absolute authority row contribution + "
                "Student-t row standard error) / objective mixed gate scale"),
            "ranking_metric_is_not_an_uncertainty_bound": True,
            "target_concentration": ranking["target_concentration"],
            "maximum_selected_fraction": ranking["maximum_selected_fraction"],
            "required_face_count_for_target": ranking[
                "required_face_count_for_target"],
            "selected_face_count": ranking["selected_face_count"],
            "selected_face_fraction": float(len(selected) / face_count),
            "selected_score_fraction": ranking["selected_score_fraction"],
            "target_reached_within_cap": ranking["target_reached_within_cap"],
            "concentration_curve": ranking["concentration_curve"],
            "selected_source_faces": [int(value) for value in selected],
            "all_vs_physically_supported_comparison": allocation_comparison,
            "top_ranked_source_faces": [{
                "rank": index + 1,
                "face_index": int(face),
                "ranking_score": float(row_score[face]),
                "material_id": int(direct["face_material_id"][face]),
                "centroid_mesh_units": [
                    float(value) for value in direct["centroids"][face]],
            } for index, face in enumerate(ranked[:top_count])],
        },
        "proposed_nested_refinement": {
            "base_rays_per_selected_source_face": BASE_RAY_LEVEL,
            "proposed_rays_per_selected_source_face": SELECTED_RAY_LEVEL,
            "unselected_source_face_rays": BASE_RAY_LEVEL,
            "append_sobol_index_interval": [16, 32],
            "same_replicate_seed": True,
            "same_per_face_cranley_shift": True,
            "same_triangle_area_sampling": True,
            "same_cellwise_certified_visibility": True,
            "same_float64_replay_authority": True,
            "estimated_additional_rays_per_replicate": extra_selected_rays,
            "global_32_additional_rays_per_replicate": extra_global_rays,
            "additional_trace_fraction_of_global_32": (
                extra_selected_rays / extra_global_rays),
            "automatic_execution_authorized": False,
            "row_selective_estimator_implementation_required": True,
            "current_global_only_estimator_is_a_blocker": True,
        },
        "recommended_common_chemistry_horizon": {
            "derivation": "dt_next/1024 from worst dt_next/16 Stage-A displacement",
            "next_step_s": next_step_s,
            "horizon_divisor": HORIZON_DIVISOR,
            "recommended_horizon_s": horizon_s,
            "stage_a_peak_gross_displacement_dx": stage_a_peak_dx,
            "linearly_predicted_gross_displacement_dx": predicted_dx,
            "fixed_geometry_limit_dx": closure.GATES[
                "maximum_gross_displacement_dx"],
            "at_least_two_x_margin": bool(
                predicted_dx <= 0.5 * closure.GATES[
                    "maximum_gross_displacement_dx"]),
            "requires_future_direct_verification": True,
        },
        "blockers": [
            "the current estimator accepts only one global rays_per_face value; a bounded "
            "row-selective nested tracing API must be reviewed before executing this plan",
            "film-specific nonlinear adjoint data were not stored by Stage A; this allocator "
            "uses exact upstream radiosity attribution and does not claim exact film attribution",
        ],
        "wall_time_s": perf_counter() - started,
        "maximum_wall_s": float(args.maximum_wall_s),
    }
    payload["plan_identity_sha256"] = _canonical_sha256({
        "schema": payload["schema"],
        "authority": payload["authority"],
        "allocation": payload["allocation"],
        "proposed_nested_refinement": payload["proposed_nested_refinement"],
        "recommended_common_chemistry_horizon": payload[
            "recommended_common_chemistry_horizon"],
    })
    if payload["wall_time_s"] > float(args.maximum_wall_s):
        raise TimeoutError("selected-source analysis exceeded its bounded wall time")
    _write_json_atomic(args.output, payload)
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage-a-audit",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "replicated_form_factor_closure" / "audit.json"))
    parser.add_argument(
        "--r19-source",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "remote_artifacts"))
    parser.add_argument(
        "--output",
        default=(ROOT / "results" / "krueger_2024_r19_response_check"
                 / "selected_source_allocator" / "audit.json"))
    parser.add_argument(
        "--target-concentration", type=float, default=TARGET_CONCENTRATION)
    parser.add_argument(
        "--maximum-selected-fraction", type=float,
        default=MAXIMUM_SELECTED_FRACTION)
    parser.add_argument(
        "--maximum-wall-s", type=float, default=MAXIMUM_WALL_S)
    args = parser.parse_args(argv)
    if (not 0.0 < args.target_concentration <= 1.0
            or not 0.0 < args.maximum_selected_fraction <= 0.25
            or not 0.0 < args.maximum_wall_s <= MAXIMUM_WALL_S):
        parser.error("invalid bounded allocator controls")
    return args


if __name__ == "__main__":
    result = run(parse_args())
    print(json.dumps({
        "status": result["status"],
        "selected_face_count": result["allocation"]["selected_face_count"],
        "selected_score_fraction": result["allocation"]["selected_score_fraction"],
        "wall_time_s": result["wall_time_s"],
    }, indent=2, sort_keys=True))
