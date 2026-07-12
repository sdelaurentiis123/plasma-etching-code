"""Flux-conservative diffuse neutral transport on an arbitrary triangle surface.

The geometry estimator supplies face-to-face diffuse form factors. This module performs the physical
multiple-reflection solve without species, material, benchmark, or aspect-ratio branches. Its unknown
is incident flux density on each face; the source/target area ratio is therefore required by diffuse
form-factor reciprocity. Omitting that ratio conserves neither particles nor the continuum equation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import gmres


@dataclass(frozen=True)
class DiffuseFormFactors3D:
    face_count: int
    source_face: np.ndarray
    target_face: np.ndarray
    transfer_fraction: np.ndarray
    escape_fraction: np.ndarray
    rays_per_face: int

    def __post_init__(self):
        source = np.asarray(self.source_face, dtype=int).copy()
        target = np.asarray(self.target_face, dtype=int).copy()
        fraction = np.asarray(self.transfer_fraction, dtype=float).copy()
        escape = np.asarray(self.escape_fraction, dtype=float).copy()
        n_face = int(self.face_count)
        n_ray = int(self.rays_per_face)
        if (n_face <= 0 or n_ray <= 0 or source.ndim != 1 or target.shape != source.shape
                or fraction.shape != source.shape or escape.shape != (n_face,)
                or np.any(source < 0) or np.any(source >= n_face)
                or np.any(target < 0) or np.any(target >= n_face)
                or np.any(~np.isfinite(fraction)) or np.any(fraction <= 0.0)
                or np.any(~np.isfinite(escape)) or np.any(escape < 0.0)):
            raise ValueError("invalid diffuse form factors")
        outgoing = escape + np.bincount(source, weights=fraction, minlength=n_face)
        if not np.allclose(outgoing, 1.0, rtol=0.0, atol=5e-13):
            raise ValueError("diffuse form factors must classify every emitted ray")
        for value in (source, target, fraction, escape):
            value.setflags(write=False)
        object.__setattr__(self, "face_count", n_face)
        object.__setattr__(self, "rays_per_face", n_ray)
        object.__setattr__(self, "source_face", source)
        object.__setattr__(self, "target_face", target)
        object.__setattr__(self, "transfer_fraction", fraction)
        object.__setattr__(self, "escape_fraction", escape)


@dataclass(frozen=True)
class DiffuseNeutralSolve3D:
    incident_flux_m2_s: np.ndarray
    reacted_flux_m2_s: np.ndarray
    reflected_flux_m2_s: np.ndarray
    source_rate_s: float
    reacted_rate_s: float
    escaped_rate_s: float
    relative_balance_error: float
    relative_linear_residual: float
    iterations_converged: bool


@dataclass(frozen=True)
class DiffuseSurfaceEmissionSolve3D:
    """Transport and reaction balance for a population emitted by the feature surface."""

    emitted_flux_m2_s: np.ndarray
    first_incident_flux_m2_s: np.ndarray
    total_incident_flux_m2_s: np.ndarray
    reacted_flux_m2_s: np.ndarray
    reflected_flux_m2_s: np.ndarray
    emitted_rate_s: float
    reacted_rate_s: float
    escaped_without_impact_rate_s: float
    escaped_after_reflection_rate_s: float
    relative_balance_error: float
    relative_linear_residual: float
    iterations_converged: bool

    def __post_init__(self):
        for name in (
                "emitted_flux_m2_s", "first_incident_flux_m2_s", "total_incident_flux_m2_s",
                "reacted_flux_m2_s", "reflected_flux_m2_s"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            value.setflags(write=False); object.__setattr__(self, name, value)


def solve_diffuse_neutral_radiosity_3d(
        direct_flux_m2_s, face_area_m2, source_face, target_face, transfer_fraction,
        escape_fraction, reaction_probability, *, relative_tolerance=1e-10,
        maximum_iterations=500):
    """Solve ``H = D + B (1-s) H`` and audit the global projectile balance.

    ``transfer_fraction[k]`` is the diffuse form factor from ``source_face[k]`` to
    ``target_face[k]``. Fractions leaving each source face, including ``escape_fraction``, must sum
    to one. ``B[i,j] = A[j] F[j->i] / A[i]`` converts reflected rate on source face ``j`` back to
    incident flux density on target face ``i``.
    """
    direct = np.asarray(direct_flux_m2_s, dtype=float)
    area = np.asarray(face_area_m2, dtype=float)
    source = np.asarray(source_face, dtype=int)
    target = np.asarray(target_face, dtype=int)
    fraction = np.asarray(transfer_fraction, dtype=float)
    escape = np.asarray(escape_fraction, dtype=float)
    reaction = np.asarray(reaction_probability, dtype=float)
    n_face = direct.size
    if (direct.ndim != 1 or area.shape != direct.shape or escape.shape != direct.shape
            or reaction.shape != direct.shape or source.ndim != 1 or target.shape != source.shape
            or fraction.shape != source.shape or np.any(~np.isfinite(direct))
            or np.any(~np.isfinite(area)) or np.any(~np.isfinite(fraction))
            or np.any(~np.isfinite(escape)) or np.any(~np.isfinite(reaction))
            or np.any(direct < 0.0) or np.any(area <= 0.0) or np.any(fraction < 0.0)
            or np.any(escape < 0.0) or np.any((reaction < 0.0) | (reaction > 1.0))
            or np.any(source < 0) or np.any(source >= n_face)
            or np.any(target < 0) or np.any(target >= n_face)
            or relative_tolerance <= 0.0 or int(maximum_iterations) <= 0):
        raise ValueError("invalid diffuse-neutral radiosity inputs")
    outgoing_fraction = escape + np.bincount(source, weights=fraction, minlength=n_face)
    if not np.allclose(outgoing_fraction, 1.0, rtol=0.0, atol=5e-13):
        raise ValueError("each face's transfer and escape fractions must sum to one")

    exchange = sparse.coo_matrix(
        (fraction * area[source] / area[target], (target, source)),
        shape=(n_face, n_face)).tocsr()
    reflection = 1.0 - reaction
    transport = exchange @ sparse.diags(reflection)
    operator = sparse.eye(n_face, format="csr") - transport
    callback_count = [0]

    def count_iteration(_):
        callback_count[0] += 1

    try:
        incident, info = gmres(
            operator, direct, rtol=relative_tolerance, atol=0.0,
            maxiter=int(maximum_iterations), callback=count_iteration,
            callback_type="pr_norm")
    except TypeError:  # scipy before callback_type/rtol
        incident, info = gmres(
            operator, direct, tol=relative_tolerance,
            maxiter=int(maximum_iterations), callback=count_iteration)
    incident = np.asarray(incident, dtype=float)
    scale = max(float(np.linalg.norm(direct)), np.finfo(float).tiny)
    residual = float(np.linalg.norm(operator @ incident - direct) / scale)
    if info != 0 or np.any(incident < -1e-12 * max(float(np.max(incident)), 1.0)):
        raise RuntimeError(
            f"diffuse-neutral radiosity did not converge to a nonnegative solution: info={info}")
    incident = np.maximum(incident, 0.0)
    reacted = reaction * incident
    reflected = reflection * incident
    source_rate = float(np.dot(area, direct))
    reacted_rate = float(np.dot(area, reacted))
    escaped_rate = float(np.dot(area * escape, reflected))
    balance = abs(source_rate - reacted_rate - escaped_rate) / max(
        source_rate, np.finfo(float).tiny)
    if balance > max(20.0 * relative_tolerance, 2e-12):
        raise RuntimeError(f"diffuse-neutral projectile balance failed: {balance:.3e}")
    for value in (incident, reacted, reflected):
        value.setflags(write=False)
    return DiffuseNeutralSolve3D(
        incident, reacted, reflected, source_rate, reacted_rate, escaped_rate,
        balance, residual, info == 0)


def transport_diffuse_surface_emission_3d(
        emitted_flux_m2_s, face_area_m2, form_factors: DiffuseFormFactors3D,
        reaction_probability, *, relative_tolerance=1e-10, maximum_iterations=500):
    """Transport a diffuse population emitted by surface faces and close its global balance.

    The first flight differs from plasma-boundary illumination: the source is an outgoing rate density
    on each face. Form factors convert it to first-incident target flux, while their escape fractions
    account for material that leaves the feature without another impact. Subsequent nonreacting impacts
    use the same diffuse radiosity equation as neutral re-emission.
    """
    emitted = np.asarray(emitted_flux_m2_s, dtype=float)
    area = np.asarray(face_area_m2, dtype=float)
    reaction = np.asarray(reaction_probability, dtype=float)
    if (not isinstance(form_factors, DiffuseFormFactors3D)
            or emitted.ndim != 1 or emitted.shape != (form_factors.face_count,)
            or area.shape != emitted.shape or reaction.shape != emitted.shape
            or np.any(~np.isfinite(emitted)) or np.any(emitted < 0.0)
            or np.any(~np.isfinite(area)) or np.any(area <= 0.0)
            or np.any(~np.isfinite(reaction))
            or np.any((reaction < 0.0) | (reaction > 1.0))):
        raise ValueError("invalid diffuse surface-emission inputs")
    source = form_factors.source_face
    target = form_factors.target_face
    fraction = form_factors.transfer_fraction
    first_incident = np.bincount(
        target, weights=(fraction * area[source] * emitted[source] / area[target]),
        minlength=form_factors.face_count)
    escaped_without_impact = float(np.dot(
        area * form_factors.escape_fraction, emitted))
    emitted_rate = float(np.dot(area, emitted))
    entered_rate = float(np.dot(area, first_incident))
    first_balance = abs(emitted_rate - entered_rate - escaped_without_impact) / max(
        emitted_rate, np.finfo(float).tiny)
    if first_balance > 2e-12:
        raise RuntimeError(f"surface-emission first-flight balance failed: {first_balance:.3e}")
    downstream = solve_diffuse_neutral_radiosity_3d(
        first_incident, area, source, target, fraction,
        form_factors.escape_fraction, reaction,
        relative_tolerance=relative_tolerance, maximum_iterations=maximum_iterations)
    escaped_total = escaped_without_impact + downstream.escaped_rate_s
    balance = abs(emitted_rate - downstream.reacted_rate_s - escaped_total) / max(
        emitted_rate, np.finfo(float).tiny)
    if balance > max(20.0 * relative_tolerance, 2e-12):
        raise RuntimeError(f"surface-emission projectile balance failed: {balance:.3e}")
    return DiffuseSurfaceEmissionSolve3D(
        emitted, first_incident, downstream.incident_flux_m2_s,
        downstream.reacted_flux_m2_s, downstream.reflected_flux_m2_s,
        emitted_rate, downstream.reacted_rate_s, escaped_without_impact,
        downstream.escaped_rate_s, balance, downstream.relative_linear_residual,
        downstream.iterations_converged)
