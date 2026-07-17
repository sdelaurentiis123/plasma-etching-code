"""Replicated randomized-QMC control for conservative diffuse form factors.

One scrambled Sobol form-factor estimate is deterministic replay evidence, not an
accuracy estimate.  This module combines at least four independently scrambled,
equal-level estimates into one conservative mean geometric operator and uses the
replicate solves to quantify downstream radiosity uncertainty.

The authority remains the existing radiosity equation solved on the mean form-factor
operator.  Replicate solutions score uncertainty; they are not averaged into a new
surface-chemistry model.  This module contains no benchmark, material, calibration,
profile-motion, or adaptive-stop policy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.stats import t as student_t

from .neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    DiffuseNeutralSolve3D,
    solve_diffuse_neutral_radiosity_3d,
)


MINIMUM_REPLICATE_COUNT = 4


def _identity_jsonable(value, path="identity"):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError(f"{path} contains a non-finite float")
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise ValueError(f"{path} keys must be nonempty strings")
            output[key] = _identity_jsonable(item, f"{path}.{key}")
        return output
    if isinstance(value, (tuple, list)):
        return [
            _identity_jsonable(item, f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise ValueError(f"{path} contains unsupported {type(value).__name__}")


def _freeze_json(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze_json(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze_json(item) for item in value)
    return value


def _copy_factors(factors):
    if not isinstance(factors, DiffuseFormFactors3D):
        raise TypeError("replicates must be DiffuseFormFactors3D")
    return DiffuseFormFactors3D(
        factors.face_count,
        factors.source_face,
        factors.target_face,
        factors.transfer_fraction,
        factors.escape_fraction,
        factors.rays_per_face,
    )


def _digest_array(digest, name, supplied):
    array = np.ascontiguousarray(np.asarray(supplied))
    metadata = json.dumps(
        {"name": name, "dtype": array.dtype.str, "shape": list(array.shape)},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    digest.update(len(metadata).to_bytes(8, "big"))
    digest.update(metadata)
    digest.update(array.tobytes(order="C"))


def _mean_form_factors(replicates):
    face_count = replicates[0].face_count
    replicate_count = len(replicates)
    pair = np.concatenate([
        item.source_face.astype(np.int64) * face_count
        + item.target_face.astype(np.int64)
        for item in replicates
    ])
    weight = np.concatenate([
        item.transfer_fraction / replicate_count for item in replicates
    ])
    if pair.size:
        unique, inverse = np.unique(pair, return_inverse=True)
        transfer = np.bincount(inverse, weights=weight, minlength=unique.size)
        source = (unique // face_count).astype(int)
        target = (unique % face_count).astype(int)
    else:
        source = np.asarray([], dtype=int)
        target = np.asarray([], dtype=int)
        transfer = np.asarray([], dtype=float)
    escape = np.mean(
        np.stack([item.escape_fraction for item in replicates]), axis=0)
    # Each replicate row already closes.  Normalize only accumulated roundoff so the
    # combined sparse row remains exactly conservative to the engine's declared tolerance.
    outgoing = escape + np.bincount(
        source, weights=transfer, minlength=face_count)
    if (np.any(~np.isfinite(outgoing)) or np.any(outgoing <= 0.0)
            or not np.allclose(outgoing, 1.0, rtol=0.0, atol=5e-13)):
        raise ValueError("replicate mean does not preserve form-factor row closure")
    transfer = transfer / outgoing[source]
    escape = escape / outgoing
    return DiffuseFormFactors3D(
        face_count,
        source,
        target,
        transfer,
        escape,
        sum(item.rays_per_face for item in replicates),
    )


@dataclass(frozen=True)
class FormFactorReciprocityDiagnostic3D:
    """Area-reciprocity defect for one sampled directed form-factor operator."""

    unordered_pair_count: int
    one_sided_pair_count: int
    absolute_l1_rate_area_m2: float
    relative_l1_error: float
    absolute_linf_rate_area_m2: float
    relative_linf_error: float


def form_factor_reciprocity_diagnostic_3d(factors, face_area_m2):
    """Measure ``A_i F_ij = A_j F_ji`` without modifying the supplied operator."""
    factors = _copy_factors(factors)
    area = np.asarray(face_area_m2, dtype=float)
    if (area.shape != (factors.face_count,) or np.any(~np.isfinite(area))
            or np.any(area <= 0.0)):
        raise ValueError("face areas must be finite, positive, and match form factors")
    conductance = {}
    for source, target, fraction in zip(
            factors.source_face, factors.target_face, factors.transfer_fraction):
        source = int(source)
        target = int(target)
        if source == target:
            continue
        conductance[(source, target)] = float(area[source] * fraction)
    unordered = sorted({tuple(sorted(pair)) for pair in conductance})
    absolute = []
    scale = []
    one_sided = 0
    for left, right in unordered:
        forward = conductance.get((left, right), 0.0)
        reverse = conductance.get((right, left), 0.0)
        if (forward == 0.0) != (reverse == 0.0):
            one_sided += 1
        absolute.append(abs(forward - reverse))
        scale.append(max(abs(forward), abs(reverse)))
    absolute = np.asarray(absolute, dtype=float)
    scale = np.asarray(scale, dtype=float)
    absolute_l1 = float(np.sum(absolute))
    scale_l1 = float(np.sum(scale))
    absolute_linf = float(np.max(absolute, initial=0.0))
    scale_linf = float(np.max(scale, initial=0.0))
    return FormFactorReciprocityDiagnostic3D(
        unordered_pair_count=len(unordered),
        one_sided_pair_count=one_sided,
        absolute_l1_rate_area_m2=absolute_l1,
        relative_l1_error=absolute_l1 / scale_l1 if scale_l1 > 0.0 else 0.0,
        absolute_linf_rate_area_m2=absolute_linf,
        relative_linf_error=(
            absolute_linf / scale_linf if scale_linf > 0.0 else 0.0),
    )


@dataclass(frozen=True)
class ReplicatedDiffuseFormFactors3D:
    """Immutable conservative mean plus independent sampled operator replicates."""

    replicate_form_factors: tuple[DiffuseFormFactors3D, ...]
    replicate_seeds: tuple[int, ...]
    face_area_m2: np.ndarray
    source_sampling: str = "externally_supplied"
    construction_identity: Mapping[str, object] = field(default_factory=dict)
    mean_form_factors: DiffuseFormFactors3D = field(init=False)
    mean_reciprocity: FormFactorReciprocityDiagnostic3D = field(init=False)
    replicate_reciprocity: tuple[FormFactorReciprocityDiagnostic3D, ...] = field(
        init=False)
    rays_per_replicate: int = field(init=False)
    total_rays_per_face: int = field(init=False)
    sha256: str = field(init=False)

    def __post_init__(self):
        factors = tuple(_copy_factors(item) for item in self.replicate_form_factors)
        seeds = tuple(int(value) for value in self.replicate_seeds)
        if len(factors) < MINIMUM_REPLICATE_COUNT or len(seeds) != len(factors):
            raise ValueError(
                f"at least {MINIMUM_REPLICATE_COUNT} form-factor replicates are required")
        if len(set(seeds)) != len(seeds):
            raise ValueError("form-factor replicate seeds must be distinct")
        order = np.argsort(np.asarray(seeds, dtype=np.int64), kind="stable")
        factors = tuple(factors[int(index)] for index in order)
        seeds = tuple(seeds[int(index)] for index in order)
        face_count = factors[0].face_count
        rays = factors[0].rays_per_face
        if any(item.face_count != face_count for item in factors):
            raise ValueError("form-factor replicate face counts differ")
        if any(item.rays_per_face != rays for item in factors):
            raise ValueError("form-factor replicates must use one common nested level")
        area = np.asarray(self.face_area_m2, dtype=float).copy()
        if (area.shape != (face_count,) or np.any(~np.isfinite(area))
                or np.any(area <= 0.0)):
            raise ValueError("face areas must be finite, positive, and match replicates")
        area.setflags(write=False)
        source_sampling = str(self.source_sampling)
        if source_sampling not in {
                "legacy_centroid", "triangle_area", "externally_supplied"}:
            raise ValueError("invalid replicated source-sampling contract")
        identity = dict(self.construction_identity)
        if not identity:
            identity = {"binding": "externally_supplied_unbound"}
        identity = _identity_jsonable(identity)
        identity_encoding = json.dumps(
            identity, sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode("utf-8")
        mean = _mean_form_factors(factors)
        replicate_reciprocity = tuple(
            form_factor_reciprocity_diagnostic_3d(item, area)
            for item in factors
        )
        mean_reciprocity = form_factor_reciprocity_diagnostic_3d(mean, area)
        digest = sha256(b"petch.replicated-diffuse-form-factors-3d.v1\0")
        digest.update(f"source_sampling={source_sampling}\n".encode("ascii"))
        digest.update(identity_encoding)
        _digest_array(digest, "face_area_m2", area)
        for seed, item in zip(seeds, factors):
            digest.update(f"seed={seed};rays={item.rays_per_face}\n".encode("ascii"))
            for name in (
                    "source_face", "target_face", "transfer_fraction",
                    "escape_fraction"):
                _digest_array(digest, f"seed={seed}:{name}", getattr(item, name))
        object.__setattr__(self, "replicate_form_factors", factors)
        object.__setattr__(self, "replicate_seeds", seeds)
        object.__setattr__(self, "face_area_m2", area)
        object.__setattr__(self, "source_sampling", source_sampling)
        object.__setattr__(self, "construction_identity", _freeze_json(identity))
        object.__setattr__(self, "mean_form_factors", mean)
        object.__setattr__(self, "mean_reciprocity", mean_reciprocity)
        object.__setattr__(self, "replicate_reciprocity", replicate_reciprocity)
        object.__setattr__(self, "rays_per_replicate", int(rays))
        object.__setattr__(self, "total_rays_per_face", int(rays * len(factors)))
        object.__setattr__(self, "sha256", digest.hexdigest())


def estimate_replicated_diffuse_form_factors_3d(
        verts, faces, centroids, gas_normals, face_area_m2, *,
        rays_per_replicate=16, replicate_seeds=(0, 1, 2, 3),
        source_sampling="triangle_area", visibility_mode="cellwise_certified",
        maximum_visibility_wraps=1024, operator_identity=None, **options):
    """Trace independent scrambled-Sobol replicates through the existing estimator."""
    from .boundary_transport_3d import estimate_diffuse_form_factors_3d

    seeds = tuple(int(value) for value in replicate_seeds)
    if len(seeds) < MINIMUM_REPLICATE_COUNT or len(set(seeds)) != len(seeds):
        raise ValueError("at least four distinct replicate seeds are required")
    if any(name in options for name in (
            "seed", "rays_per_face", "source_sampling", "visibility_mode",
            "maximum_visibility_wraps", "return_visibility_receipt")):
        raise ValueError(
            "sampling and visibility controls are controlled by the replicate contract")
    receipts = tuple(
        estimate_diffuse_form_factors_3d(
            verts, faces, centroids, gas_normals,
            rays_per_face=rays_per_replicate, seed=seed,
            source_sampling=source_sampling, visibility_mode=visibility_mode,
            maximum_visibility_wraps=maximum_visibility_wraps,
            return_visibility_receipt=True, **options)
        for seed in seeds
    )
    estimates = tuple(item.form_factors for item in receipts)
    mesh_digest = sha256(b"petch.diffuse-form-factor-mesh.v1\0")
    for name, values in (
            ("verts", verts), ("faces", faces), ("centroids", centroids),
            ("gas_normals", gas_normals), ("face_area_m2", face_area_m2)):
        _digest_array(mesh_digest, name, values)
    verts_array = np.asarray(verts, dtype=float)
    domain = options.get("domain_size")
    if domain is None:
        domain = np.maximum(np.ptp(verts_array, axis=0), 1.0)
    domain = np.asarray(domain, dtype=float)
    identity = {
        "estimator": "petch.diffuse-form-factor-area-direction-rqmc.v1",
        "source_sampling": str(source_sampling),
        "sampling_dimension": 4 if source_sampling == "triangle_area" else 2,
        "mesh_sha256": mesh_digest.hexdigest(),
        "domain_size": domain.tolist(),
        "periodic_lateral": bool(options.get("periodic_lateral", False)),
        "ray_offset": float(options.get("ray_offset", 1.0e-5)),
        "device": str(options.get("device", "engine_default")),
        "replicate_count": len(seeds),
        "rays_per_replicate": int(rays_per_replicate),
        "construction_call_count": len(seeds),
        "visibility_operator": str(visibility_mode),
        "visibility_certification": (
            "full_float64_reference"
            if visibility_mode == "float64_reference" else
            "full_event_float64_parity_14664_ray_real_checkpoint"
            if visibility_mode == "cellwise_certified" else
            "ambiguous_miss_and_invalid_hit_float64_replay_full_event_parity_pending"
            if visibility_mode == "replay_hardened" else
            "legacy_uncertified"),
        "visibility_float64_evaluated_count": sum(
            item.float64_evaluated_count for item in receipts),
        "visibility_recovered_hit_count": sum(
            item.float64_recovered_hit_count for item in receipts),
        "visibility_open_escape_count": sum(
            item.open_escape_count for item in receipts),
        "visibility_maximum_wrap_count": max(
            item.maximum_wrap_count for item in receipts),
        "visibility_launch_inset_count": sum(
            item.launch_inset_count for item in receipts),
        "visibility_centroid_limit_count": sum(
            item.centroid_limit_count for item in receipts),
        "maximum_visibility_wraps": int(maximum_visibility_wraps),
        "caller_operator_identity": (
            {} if operator_identity is None else operator_identity),
    }
    return ReplicatedDiffuseFormFactors3D(
        estimates, seeds, face_area_m2, source_sampling=source_sampling,
        construction_identity=identity)


@dataclass(frozen=True)
class RadiosityReplicateUncertainty3D:
    """Mean-operator authority and replicate uncertainty for one neutral species."""

    authority: DiffuseNeutralSolve3D
    replicate_solutions: tuple[DiffuseNeutralSolve3D, ...]
    confidence_level: float
    student_t_critical: float
    replicate_mean_incident_flux_m2_s: np.ndarray
    incident_standard_error_m2_s: np.ndarray
    incident_confidence_half_width_m2_s: np.ndarray
    replicate_mean_reacted_flux_m2_s: np.ndarray
    reacted_standard_error_m2_s: np.ndarray
    reacted_confidence_half_width_m2_s: np.ndarray
    incident_relative_confidence_linf: float
    incident_area_weighted_relative_confidence_l1: float
    reacted_relative_confidence_linf: float
    reacted_area_weighted_relative_confidence_l1: float
    authority_to_replicate_mean_incident_relative_linf: float
    authority_to_replicate_mean_incident_area_weighted_relative_l1: float

    def __post_init__(self):
        for name in (
                "replicate_mean_incident_flux_m2_s", "incident_standard_error_m2_s",
                "incident_confidence_half_width_m2_s",
                "replicate_mean_reacted_flux_m2_s", "reacted_standard_error_m2_s",
                "reacted_confidence_half_width_m2_s"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("replicate uncertainty arrays must be finite and nonnegative")
            value.setflags(write=False)
            object.__setattr__(self, name, value)


def _relative_array_scores(reference, difference, area):
    reference = np.asarray(reference, dtype=float)
    difference = np.asarray(difference, dtype=float)
    linf_scale = max(float(np.max(np.abs(reference), initial=0.0)), 1.0)
    l1_scale = max(float(np.sum(area * np.abs(reference))), 1.0)
    return (
        float(np.max(np.abs(difference), initial=0.0)) / linf_scale,
        float(np.sum(area * np.abs(difference))) / l1_scale,
    )


def solve_replicated_diffuse_neutral_radiosity_3d(
        direct_flux_m2_s, ensemble, reaction_probability, *,
        relative_tolerance=1.0e-10, maximum_iterations=500,
        confidence_level=0.95):
    """Solve the mean operator and score it with independent replicate operators."""
    if not isinstance(ensemble, ReplicatedDiffuseFormFactors3D):
        raise TypeError("ensemble must be ReplicatedDiffuseFormFactors3D")
    confidence = float(confidence_level)
    if not np.isfinite(confidence) or not 0.0 < confidence < 1.0:
        raise ValueError("confidence_level must lie strictly between zero and one")
    common = dict(
        emitted_flux_m2_s=direct_flux_m2_s,
        face_area_m2=ensemble.face_area_m2,
        reaction_probability=reaction_probability,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
    )

    def solve(factors):
        return solve_diffuse_neutral_radiosity_3d(
            common["emitted_flux_m2_s"], common["face_area_m2"],
            factors.source_face, factors.target_face, factors.transfer_fraction,
            factors.escape_fraction, common["reaction_probability"],
            relative_tolerance=common["relative_tolerance"],
            maximum_iterations=common["maximum_iterations"])

    authority = solve(ensemble.mean_form_factors)
    replicates = tuple(solve(item) for item in ensemble.replicate_form_factors)
    incident = np.stack([item.incident_flux_m2_s for item in replicates])
    reacted = np.stack([item.reacted_flux_m2_s for item in replicates])
    count = len(replicates)
    critical = float(student_t.ppf(0.5 + 0.5 * confidence, count - 1))

    def statistics(values):
        mean = np.mean(values, axis=0)
        standard_error = np.std(values, axis=0, ddof=1) / np.sqrt(count)
        return mean, standard_error, critical * standard_error

    incident_mean, incident_se, incident_half = statistics(incident)
    reacted_mean, reacted_se, reacted_half = statistics(reacted)
    incident_ci = _relative_array_scores(
        authority.incident_flux_m2_s, incident_half, ensemble.face_area_m2)
    reacted_ci = _relative_array_scores(
        authority.reacted_flux_m2_s, reacted_half, ensemble.face_area_m2)
    authority_difference = incident_mean - authority.incident_flux_m2_s
    authority_score = _relative_array_scores(
        authority.incident_flux_m2_s, authority_difference,
        ensemble.face_area_m2)
    return RadiosityReplicateUncertainty3D(
        authority=authority,
        replicate_solutions=replicates,
        confidence_level=confidence,
        student_t_critical=critical,
        replicate_mean_incident_flux_m2_s=incident_mean,
        incident_standard_error_m2_s=incident_se,
        incident_confidence_half_width_m2_s=incident_half,
        replicate_mean_reacted_flux_m2_s=reacted_mean,
        reacted_standard_error_m2_s=reacted_se,
        reacted_confidence_half_width_m2_s=reacted_half,
        incident_relative_confidence_linf=incident_ci[0],
        incident_area_weighted_relative_confidence_l1=incident_ci[1],
        reacted_relative_confidence_linf=reacted_ci[0],
        reacted_area_weighted_relative_confidence_l1=reacted_ci[1],
        authority_to_replicate_mean_incident_relative_linf=authority_score[0],
        authority_to_replicate_mean_incident_area_weighted_relative_l1=(
            authority_score[1]),
    )
