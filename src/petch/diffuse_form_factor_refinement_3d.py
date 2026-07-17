"""Bounded row-selective nested-RQMC refinement for diffuse form factors.

The production form-factor estimator emits one equal-level rule for every source
face.  Accuracy audits can show that only a predeclared subset of source rows is
worth extending.  This module appends the *same* scrambled Sobol sequence indices
that a global refinement would use, traces only those source rows, and merges
integer hit/escape counts before converting back to fractions.

The returned operator is suitable for the existing radiosity equation.  Its
``rays_per_face`` field is deliberately the minimum (unselected-row) level for
backward compatibility; ``row_ray_count`` in the immutable receipt is the
authoritative heterogeneous sampling ledger.  No benchmark policy, tolerance,
patch selection, or automatic refinement decision lives here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_transport_3d import (
    DIFFUSE_VISIBILITY_EMERGENCY_MAXIMUM_WRAPS,
    _diffuse_form_factor_ray_sample_block_3d,
    trace_diffuse_form_factor_events_cellwise_certified_3d,
)
from .neutral_radiosity_3d import DiffuseFormFactors3D


MAXIMUM_SELECTED_SOURCE_FRACTION = 0.25


class DiffuseFormFactorRefinementRefusal(RuntimeError):
    """A selective refinement would violate its declared numerical contract."""


def _jsonable(value, path="identity"):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        result = float(value)
        if not np.isfinite(result):
            raise ValueError(f"{path} contains a non-finite float")
        return result
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise ValueError(f"{path} keys must be nonempty strings")
            output[key] = _jsonable(item, f"{path}.{key}")
        return output
    if isinstance(value, (tuple, list)):
        return [
            _jsonable(item, f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise ValueError(f"{path} contains unsupported {type(value).__name__}")


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _encode_identity(value):
    canonical = _jsonable(value)
    return canonical, json.dumps(
        canonical, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")


def _digest_array(digest, name, supplied):
    array = np.ascontiguousarray(np.asarray(supplied))
    metadata = json.dumps(
        {"name": name, "dtype": array.dtype.str, "shape": list(array.shape)},
        sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    digest.update(len(metadata).to_bytes(8, "big"))
    digest.update(metadata)
    digest.update(array.tobytes(order="C"))


def _power_of_two(value):
    return value > 0 and not value & (value - 1)


def _integer_bincount(index, weight, minlength):
    """Accumulate an integer ledger without NumPy's weighted-float conversion."""
    index = np.asarray(index, dtype=int)
    weight = np.asarray(weight, dtype=np.int64)
    if index.shape != weight.shape:
        raise ValueError("integer bincount index and weight shapes differ")
    result = np.zeros(int(minlength), dtype=np.int64)
    np.add.at(result, index, weight)
    return result


@dataclass(frozen=True)
class DiffuseFormFactorOperatorIdentity3D:
    """Canonical identity of one uniform nested form-factor construction."""

    payload: Mapping[str, object]
    sha256: str = field(init=False)
    _encoding: bytes = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        payload, encoding = _encode_identity(dict(self.payload))
        if payload.get("schema") != "petch.diffuse-form-factor-operator.v1":
            raise ValueError("invalid diffuse form-factor operator identity schema")
        object.__setattr__(self, "payload", _freeze(payload))
        object.__setattr__(self, "_encoding", encoding)
        object.__setattr__(self, "sha256", sha256(encoding).hexdigest())


def diffuse_form_factor_operator_identity_3d(
        verts, faces, centroids, gas_normals, *, rays_per_face, seed,
        domain_size=None, periodic_lateral=False, ray_offset=1.0e-5,
        source_sampling="triangle_area", visibility_mode="cellwise_certified",
        maximum_visibility_wraps=1024, maximum_visibility_replay_wraps=None,
        device=None, operator_identity=None, form_factors=None):
    """Bind a sampled operator to exact mesh, seed, and visibility controls."""
    vertices = np.asarray(verts, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    centers = np.asarray(centroids, dtype=float)
    normals = np.asarray(gas_normals, dtype=float)
    rays = int(rays_per_face)
    wraps = int(maximum_visibility_wraps)
    replay_wraps = wraps if maximum_visibility_replay_wraps is None else int(
        maximum_visibility_replay_wraps)
    if (vertices.ndim != 2 or vertices.shape[1:] != (3,)
            or triangles.ndim != 2 or triangles.shape[1:] != (3,)
            or centers.shape != (len(triangles), 3)
            or normals.shape != centers.shape
            or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(centers))
            or np.any(~np.isfinite(normals)) or np.any(triangles < 0)
            or np.any(triangles >= len(vertices))
            or int(rays_per_face) != rays_per_face or not _power_of_two(rays)
            or int(maximum_visibility_wraps) != maximum_visibility_wraps
            or wraps <= 0 or replay_wraps < wraps
            or not np.isfinite(ray_offset) or ray_offset <= 0.0
            or source_sampling != "triangle_area"
            or visibility_mode != "cellwise_certified"):
        raise ValueError("invalid diffuse form-factor operator identity inputs")
    normal_length = np.linalg.norm(normals, axis=1)
    if not np.allclose(normal_length, 1.0, rtol=0.0, atol=2.0e-6):
        raise ValueError("gas normals must be unit length")
    if domain_size is None:
        domain = np.maximum(np.ptp(vertices, axis=0), 1.0)
    else:
        domain = np.asarray(domain_size, dtype=float)
    if domain.shape != (3,) or np.any(~np.isfinite(domain)) or np.any(domain <= 0.0):
        raise ValueError("domain_size must contain three positive lengths")
    if periodic_lateral and (
            np.min(vertices[:, :2]) < -1.0e-7
            or np.any(np.max(vertices[:, :2], axis=0) > domain[:2] + 1.0e-7)
            or np.max(vertices[:, 2]) > domain[2] + 1.0e-7):
        raise ValueError("periodic mesh must lie inside [0, domain_size]")

    mesh_digest = sha256(b"petch.diffuse-form-factor-refinement-mesh.v1\0")
    for name, value in (
            ("verts", vertices), ("faces", triangles),
            ("centroids", centers), ("gas_normals", normals)):
        _digest_array(mesh_digest, name, value)
    sampled_digest = None
    if form_factors is not None:
        if (not isinstance(form_factors, DiffuseFormFactors3D)
                or form_factors.face_count != len(triangles)
                or form_factors.rays_per_face != rays):
            raise ValueError("bound form factors do not match the operator identity")
        digest = sha256(b"petch.sampled-diffuse-form-factors.v1\0")
        for name in (
                "source_face", "target_face", "transfer_fraction",
                "escape_fraction"):
            _digest_array(digest, name, getattr(form_factors, name))
        digest.update(f"rays_per_face={rays}\n".encode("ascii"))
        sampled_digest = digest.hexdigest()
    caller, _encoding = _encode_identity(
        {} if operator_identity is None else operator_identity)
    return DiffuseFormFactorOperatorIdentity3D({
        "schema": "petch.diffuse-form-factor-operator.v1",
        "estimator": "petch.diffuse-form-factor-area-direction-rqmc.v1",
        "mesh_sha256": mesh_digest.hexdigest(),
        "sampled_form_factors_sha256": sampled_digest,
        "face_count": int(len(triangles)),
        "rays_per_face": rays,
        "seed": int(seed),
        "source_sampling": str(source_sampling),
        "sampling_dimension": 4,
        "visibility_mode": str(visibility_mode),
        "periodic_lateral": bool(periodic_lateral),
        "domain_size": domain.tolist(),
        "ray_offset": float(ray_offset),
        "maximum_visibility_wraps": wraps,
        "maximum_visibility_replay_wraps": replay_wraps,
        "exact_replay_horizon_policy": "configured_then_geometry_derived_open_top_v1",
        "derived_replay_emergency_maximum_wraps": (
            DIFFUSE_VISIBILITY_EMERGENCY_MAXIMUM_WRAPS),
        "device": "engine_default" if device is None else str(device),
        "caller_operator_identity": caller,
    })


def _integer_count_ledger(factors):
    if not isinstance(factors, DiffuseFormFactors3D):
        raise TypeError("base_form_factors must be DiffuseFormFactors3D")
    rays = int(factors.rays_per_face)
    scaled_transfer = np.asarray(factors.transfer_fraction) * rays
    scaled_escape = np.asarray(factors.escape_fraction) * rays
    transfer = np.rint(scaled_transfer).astype(np.int64)
    escape = np.rint(scaled_escape).astype(np.int64)
    tolerance = 32.0 * np.finfo(float).eps * max(rays, 1)
    if (not np.allclose(scaled_transfer, transfer, rtol=0.0, atol=tolerance)
            or not np.allclose(scaled_escape, escape, rtol=0.0, atol=tolerance)):
        raise DiffuseFormFactorRefinementRefusal(
            "base fractions do not reconstruct integer ray counts")
    classified = escape + _integer_bincount(
        factors.source_face, transfer, factors.face_count)
    if not np.array_equal(classified, np.full(factors.face_count, rays, dtype=np.int64)):
        raise DiffuseFormFactorRefinementRefusal(
            "base integer ray-count ledger does not close")
    return transfer, escape


@dataclass(frozen=True)
class NestedRowDiffuseFormFactorReceipt3D:
    """Immutable refined operator, integer count ledger, and provenance receipt."""

    form_factors: DiffuseFormFactors3D
    selected_source_face: np.ndarray
    row_ray_count: np.ndarray
    transfer_ray_count: np.ndarray
    escape_ray_count: np.ndarray
    base_rays_per_face: int
    refined_rays_per_face: int
    seed: int
    selected_source_fraction_cap: float
    base_operator_identity: DiffuseFormFactorOperatorIdentity3D
    refined_operator_identity: DiffuseFormFactorOperatorIdentity3D
    traced_ray_count: int
    float64_evaluated_count: int
    float64_recovered_hit_count: int
    open_escape_count: int
    maximum_wrap_count: int
    launch_inset_count: int
    centroid_limit_count: int
    construction_identity: Mapping[str, object]
    sha256: str = field(init=False)

    def __post_init__(self):
        if not isinstance(self.form_factors, DiffuseFormFactors3D):
            raise TypeError("form_factors must be DiffuseFormFactors3D")
        if (not isinstance(self.base_operator_identity,
                           DiffuseFormFactorOperatorIdentity3D)
                or not isinstance(self.refined_operator_identity,
                                  DiffuseFormFactorOperatorIdentity3D)):
            raise TypeError("refinement identities must be bound operator identities")
        selected = np.asarray(self.selected_source_face, dtype=int).copy()
        row_count = np.asarray(self.row_ray_count, dtype=np.int64).copy()
        transfer_count = np.asarray(self.transfer_ray_count, dtype=np.int64).copy()
        escape_count = np.asarray(self.escape_ray_count, dtype=np.int64).copy()
        n_face = self.form_factors.face_count
        base = int(self.base_rays_per_face)
        refined = int(self.refined_rays_per_face)
        cap = float(self.selected_source_fraction_cap)
        diagnostic_values = (
            self.traced_ray_count, self.float64_evaluated_count,
            self.float64_recovered_hit_count,
            self.open_escape_count, self.maximum_wrap_count,
            self.launch_inset_count, self.centroid_limit_count)
        if (selected.ndim != 1 or selected.size == 0
                or np.any(selected < 0) or np.any(selected >= n_face)
                or not np.array_equal(selected, np.unique(selected))
                or row_count.shape != (n_face,)
                or transfer_count.shape != self.form_factors.transfer_fraction.shape
                or escape_count.shape != (n_face,)
                or np.any(row_count <= 0) or np.any(transfer_count <= 0)
                or np.any(escape_count < 0)
                or not _power_of_two(base) or not _power_of_two(refined)
                or refined <= base or self.form_factors.rays_per_face != base
                or not 0.0 < cap <= MAXIMUM_SELECTED_SOURCE_FRACTION
                or selected.size / n_face > cap + 8.0 * np.finfo(float).eps
                or any(int(value) != value or value < 0 for value in diagnostic_values)
                or self.traced_ray_count != selected.size * (refined - base)
                or self.float64_recovered_hit_count > self.float64_evaluated_count
                or self.float64_evaluated_count > self.traced_ray_count
                or self.open_escape_count > self.traced_ray_count
                or self.centroid_limit_count > self.launch_inset_count):
            raise ValueError("invalid nested row-refinement receipt")
        expected_row_count = np.full(n_face, base, dtype=np.int64)
        expected_row_count[selected] = refined
        if not np.array_equal(row_count, expected_row_count):
            raise ValueError("row ray-count ledger does not match selected refinement")
        classified = escape_count + _integer_bincount(
            self.form_factors.source_face, transfer_count, n_face)
        if not np.array_equal(classified, row_count):
            raise ValueError("refined integer ray-count ledger does not close")
        expected_transfer = transfer_count / row_count[self.form_factors.source_face]
        expected_escape = escape_count / row_count
        if (not np.array_equal(expected_transfer, self.form_factors.transfer_fraction)
                or not np.array_equal(expected_escape, self.form_factors.escape_fraction)):
            raise ValueError("refined fractions do not reproduce the integer count ledger")
        identity, encoding = _encode_identity(dict(self.construction_identity))
        if identity.get("schema") != "petch.nested-row-form-factor-refinement.v1":
            raise ValueError("invalid nested row-refinement construction identity")
        digest = sha256(b"petch.nested-row-form-factor-refinement.v1\0")
        digest.update(encoding)
        for name, value in (
                ("selected_source_face", selected), ("row_ray_count", row_count),
                ("source_face", self.form_factors.source_face),
                ("target_face", self.form_factors.target_face),
                ("transfer_ray_count", transfer_count),
                ("escape_ray_count", escape_count),
                ("transfer_fraction", self.form_factors.transfer_fraction),
                ("escape_fraction", self.form_factors.escape_fraction)):
            _digest_array(digest, name, value)
        for value in (selected, row_count, transfer_count, escape_count):
            value.setflags(write=False)
        object.__setattr__(self, "selected_source_face", selected)
        object.__setattr__(self, "row_ray_count", row_count)
        object.__setattr__(self, "transfer_ray_count", transfer_count)
        object.__setattr__(self, "escape_ray_count", escape_count)
        object.__setattr__(self, "base_rays_per_face", base)
        object.__setattr__(self, "refined_rays_per_face", refined)
        object.__setattr__(self, "seed", int(self.seed))
        object.__setattr__(self, "selected_source_fraction_cap", cap)
        for name in (
                "traced_ray_count", "float64_evaluated_count",
                "float64_recovered_hit_count",
                "open_escape_count", "maximum_wrap_count", "launch_inset_count",
                "centroid_limit_count"):
            object.__setattr__(self, name, int(getattr(self, name)))
        object.__setattr__(self, "construction_identity", _freeze(identity))
        object.__setattr__(self, "sha256", digest.hexdigest())


def refine_diffuse_form_factor_rows_nested_3d(
        base_form_factors, base_operator_identity, verts, faces, centroids,
        gas_normals, *, selected_source_faces, refined_rays_per_face,
        expected_operator_identity=None, selected_source_fraction_cap=0.25):
    """Append nested Sobol indices for selected source rows only.

    ``base_operator_identity`` must have been constructed for the supplied uniform
    base operator.  The function re-derives that identity from the current mesh and
    caller epoch before tracing, so a seed, geometry, visibility, or operator-epoch
    mismatch refuses.  Only ``triangle_area`` + ``cellwise_certified`` is accepted.
    """
    if not isinstance(base_form_factors, DiffuseFormFactors3D):
        raise TypeError("base_form_factors must be DiffuseFormFactors3D")
    if not isinstance(base_operator_identity, DiffuseFormFactorOperatorIdentity3D):
        raise TypeError("base_operator_identity must be a bound operator identity")
    payload = dict(base_operator_identity.payload)
    if payload.get("sampled_form_factors_sha256") is None:
        raise DiffuseFormFactorRefinementRefusal(
            "base operator identity is not bound to sampled form factors")
    caller_expected, _encoding = _encode_identity(
        {} if expected_operator_identity is None else expected_operator_identity)
    if caller_expected != payload["caller_operator_identity"]:
        raise DiffuseFormFactorRefinementRefusal(
            "caller operator identity does not match the base operator")
    device_value = payload["device"]
    device = None if device_value == "engine_default" else device_value
    reproduced_identity = diffuse_form_factor_operator_identity_3d(
        verts, faces, centroids, gas_normals,
        rays_per_face=base_form_factors.rays_per_face,
        seed=payload["seed"], domain_size=payload["domain_size"],
        periodic_lateral=payload["periodic_lateral"],
        ray_offset=payload["ray_offset"],
        source_sampling=payload["source_sampling"],
        visibility_mode=payload["visibility_mode"],
        maximum_visibility_wraps=payload["maximum_visibility_wraps"],
        maximum_visibility_replay_wraps=payload[
            "maximum_visibility_replay_wraps"],
        device=device, operator_identity=caller_expected,
        form_factors=base_form_factors)
    if reproduced_identity.sha256 != base_operator_identity.sha256:
        raise DiffuseFormFactorRefinementRefusal(
            "mesh or sampling operator does not match the base operator identity")
    n_face = base_form_factors.face_count
    if n_face != int(payload["face_count"]):
        raise DiffuseFormFactorRefinementRefusal(
            "base form factors use another face count")
    base_rays = int(base_form_factors.rays_per_face)
    refined_rays = int(refined_rays_per_face)
    if (int(refined_rays_per_face) != refined_rays_per_face
            or not _power_of_two(base_rays) or not _power_of_two(refined_rays)
            or refined_rays <= base_rays):
        raise DiffuseFormFactorRefinementRefusal(
            "base and refined levels must be increasing powers of two")
    selected_supplied = np.asarray(selected_source_faces)
    if (selected_supplied.ndim != 1 or selected_supplied.size == 0
            or not np.issubdtype(selected_supplied.dtype, np.integer)):
        raise DiffuseFormFactorRefinementRefusal(
            "selected source faces must be a nonempty integer vector")
    selected = selected_supplied.astype(int)
    if (np.any(selected < 0) or np.any(selected >= n_face)
            or len(np.unique(selected)) != len(selected)):
        raise DiffuseFormFactorRefinementRefusal(
            "selected source faces are duplicate or out of range")
    selected = np.sort(selected)
    cap = float(selected_source_fraction_cap)
    if (not np.isfinite(cap) or cap <= 0.0
            or cap > MAXIMUM_SELECTED_SOURCE_FRACTION
            or selected.size / n_face > cap + 8.0 * np.finfo(float).eps):
        raise DiffuseFormFactorRefinementRefusal(
            "selected source fraction exceeds the predeclared 25% cap")

    base_transfer_count, base_escape_count = _integer_count_ledger(
        base_form_factors)
    vertices = np.asarray(verts, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    centers = np.asarray(centroids, dtype=float)
    normals = np.asarray(gas_normals, dtype=float)
    authority_normals = normals / np.linalg.norm(normals, axis=1)[:, None]
    source, origin, direction, launch = _diffuse_form_factor_ray_sample_block_3d(
        vertices, triangles, centers, authority_normals,
        source_faces=selected, sobol_index_start=base_rays,
        sobol_index_stop=refined_rays, seed=payload["seed"],
        ray_offset=payload["ray_offset"],
        source_sampling=payload["source_sampling"],
        return_launch_diagnostics=True)
    events = trace_diffuse_form_factor_events_cellwise_certified_3d(
        origin, direction, vertices, triangles, authority_normals,
        domain_size=payload["domain_size"],
        periodic_lateral=payload["periodic_lateral"],
        maximum_wraps=payload["maximum_visibility_wraps"],
        maximum_exact_replay_wraps=payload[
            "maximum_visibility_replay_wraps"],
        device=device)

    escaped = events.hit_face < 0
    delta_escape = np.bincount(
        source[escaped], minlength=n_face).astype(np.int64)
    valid_source = source[~escaped]
    valid_target = events.hit_face[~escaped]
    delta_pair = valid_source.astype(np.int64) * n_face + valid_target
    delta_unique, delta_count = np.unique(delta_pair, return_counts=True)
    base_pair = (
        base_form_factors.source_face.astype(np.int64) * n_face
        + base_form_factors.target_face.astype(np.int64))
    pair = np.concatenate((base_pair, delta_unique))
    count = np.concatenate((base_transfer_count, delta_count.astype(np.int64)))
    unique, inverse = np.unique(pair, return_inverse=True)
    transfer_count = _integer_bincount(inverse, count, len(unique))
    source_face = (unique // n_face).astype(int)
    target_face = (unique % n_face).astype(int)
    escape_count = base_escape_count + delta_escape
    row_count = np.full(n_face, base_rays, dtype=np.int64)
    row_count[selected] = refined_rays
    classified = escape_count + _integer_bincount(
        source_face, transfer_count, n_face)
    if not np.array_equal(classified, row_count):
        raise RuntimeError("selective form-factor integer merge failed row closure")
    factors = DiffuseFormFactors3D(
        n_face, source_face, target_face,
        transfer_count / row_count[source_face], escape_count / row_count,
        base_rays)

    base_untouched = ~np.isin(base_form_factors.source_face, selected)
    refined_untouched = ~np.isin(factors.source_face, selected)
    untouched_face = np.setdiff1d(
        np.arange(n_face), selected, assume_unique=True)
    base_order = np.argsort(base_pair[base_untouched], kind="stable")
    refined_order = np.argsort(unique[refined_untouched], kind="stable")
    if (not np.array_equal(
                base_pair[base_untouched][base_order],
                unique[refined_untouched][refined_order])
            or not np.array_equal(
                base_form_factors.transfer_fraction[base_untouched][base_order],
                factors.transfer_fraction[refined_untouched][refined_order])
            or not np.array_equal(
                base_form_factors.escape_fraction[untouched_face],
                factors.escape_fraction[untouched_face])):
        raise RuntimeError("selective refinement changed an unselected source row")

    refined_identity = diffuse_form_factor_operator_identity_3d(
        vertices, triangles, centers, normals, rays_per_face=refined_rays,
        seed=payload["seed"], domain_size=payload["domain_size"],
        periodic_lateral=payload["periodic_lateral"],
        ray_offset=payload["ray_offset"],
        source_sampling=payload["source_sampling"],
        visibility_mode=payload["visibility_mode"],
        maximum_visibility_wraps=payload["maximum_visibility_wraps"],
        maximum_visibility_replay_wraps=payload[
            "maximum_visibility_replay_wraps"],
        device=device, operator_identity=caller_expected)
    construction = {
        "schema": "petch.nested-row-form-factor-refinement.v1",
        "base_operator_identity_sha256": base_operator_identity.sha256,
        "refined_operator_identity_sha256": refined_identity.sha256,
        "seed": int(payload["seed"]),
        "base_rays_per_face": base_rays,
        "refined_rays_per_face": refined_rays,
        "sobol_index_interval": [base_rays, refined_rays],
        "selected_source_face": selected.tolist(),
        "selected_source_fraction": float(selected.size / n_face),
        "selected_source_fraction_cap": cap,
        "count_merge": "integer_hit_and_escape_counts_before_fraction_conversion",
        "untouched_row_contract": "bitwise_equal_fraction_rows",
        "source_sampling": payload["source_sampling"],
        "visibility_mode": payload["visibility_mode"],
    }
    return NestedRowDiffuseFormFactorReceipt3D(
        form_factors=factors, selected_source_face=selected,
        row_ray_count=row_count, transfer_ray_count=transfer_count,
        escape_ray_count=escape_count, base_rays_per_face=base_rays,
        refined_rays_per_face=refined_rays, seed=payload["seed"],
        selected_source_fraction_cap=cap,
        base_operator_identity=base_operator_identity,
        refined_operator_identity=refined_identity,
        traced_ray_count=len(origin),
        float64_evaluated_count=events.replay_count,
        float64_recovered_hit_count=events.recovered_hit_count,
        open_escape_count=events.open_escape_count,
        maximum_wrap_count=events.maximum_wrap_count,
        launch_inset_count=launch["launch_inset_count"],
        centroid_limit_count=launch["centroid_limit_count"],
        construction_identity=construction)
