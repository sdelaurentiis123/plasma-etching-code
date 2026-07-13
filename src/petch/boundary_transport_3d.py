"""Species-resolved 3-D transport from ``PlasmaBoundaryState`` to an arbitrary triangle mesh.

This module is the dimensional bridge between the common reactor/sheath boundary contract and surface
kinetics.  Its forward and reversible-adjoint backends are deliberately limited to collisionless,
absorbing, first-hit transport.  They preserve the declared boundary measure at every hit and report
their limitations; reflection, re-emission, and collisional charged transport remain out of scope.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.stats import qmc
import warp as wp

from .boundary_state import (
    FoldedNormalTangentialDensity, MixtureBoundaryDensity, PlasmaBoundaryState,
    qmc_boundary_proposal,
)
from .surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes
from .neutral_radiosity_3d import DiffuseFormFactors3D
from .threed import DEVICE, _apply_bc
from .warp_runtime import ensure_writable_warp_cache


def _contains_surface_local_density(density):
    if isinstance(density, FoldedNormalTangentialDensity):
        return True
    if isinstance(density, MixtureBoundaryDensity):
        return any(_contains_surface_local_density(item) for item in density.components)
    return False


@wp.kernel
def _first_hit_events_3d(
        mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), direction: wp.array(dtype=wp.vec3),
        max_distance: float, hit_face: wp.array(dtype=int), hit_cosine: wp.array(dtype=float)):
    particle = wp.tid()
    ray = wp.mesh_query_ray(mesh, origin[particle], direction[particle], max_distance)
    if ray.result:
        normal = ray.normal
        if wp.dot(normal, direction[particle]) > 0.0:
            normal = -normal
        cosine = -wp.dot(direction[particle], normal)
        hit_face[particle] = ray.face
        hit_cosine[particle] = wp.clamp(cosine, 0.0, 1.0)


@wp.kernel
def _diffuse_form_factor_events_3d(
        mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), direction: wp.array(dtype=wp.vec3),
        domain_x: float, domain_y: float, domain_z: float, periodic_lateral: int,
        hit_face: wp.array(dtype=int)):
    ray_index = wp.tid()
    boundary_ray = _apply_bc(
        mesh, origin[ray_index], direction[ray_index],
        domain_x, domain_y, domain_z, periodic_lateral)
    ray = wp.mesh_query_ray(mesh, boundary_ray.o, boundary_ray.d, 1.0e6)
    if ray.result:
        hit_face[ray_index] = ray.face


@wp.kernel
def _periodic_first_hit_events_3d(
        mesh: wp.uint64, origin: wp.array(dtype=wp.vec3), direction: wp.array(dtype=wp.vec3),
        domain_x: float, domain_y: float, domain_z: float,
        hit_face: wp.array(dtype=int), hit_cosine: wp.array(dtype=float)):
    particle = wp.tid()
    boundary_ray = _apply_bc(
        mesh, origin[particle], direction[particle], domain_x, domain_y, domain_z, 1)
    ray = wp.mesh_query_ray(mesh, boundary_ray.o, boundary_ray.d, 1.0e6)
    if ray.result:
        normal = ray.normal
        if wp.dot(normal, boundary_ray.d) > 0.0:
            normal = -normal
        hit_face[particle] = ray.face
        hit_cosine[particle] = wp.clamp(-wp.dot(boundary_ray.d, normal), 0.0, 1.0)


@wp.func
def _trilinear_electric_field_3d(
        potential: wp.array3d(dtype=float), position: wp.vec3,
        grid_origin: wp.vec3, spacing: wp.vec3):
    displacement = position - grid_origin
    coordinate = wp.vec3(
        displacement[0] / spacing[0], displacement[1] / spacing[1],
        displacement[2] / spacing[2])
    nx = potential.shape[0]; ny = potential.shape[1]; nz = potential.shape[2]
    i = wp.int32(wp.floor(coordinate[0])); j = wp.int32(wp.floor(coordinate[1]))
    k = wp.int32(wp.floor(coordinate[2]))
    i = wp.clamp(i, 0, nx - 2); j = wp.clamp(j, 0, ny - 2); k = wp.clamp(k, 0, nz - 2)
    fx = wp.clamp(coordinate[0] - float(i), 0.0, 1.0)
    fy = wp.clamp(coordinate[1] - float(j), 0.0, 1.0)
    fz = wp.clamp(coordinate[2] - float(k), 0.0, 1.0)
    one = float(1.0)
    dvx = (
        (one - fy) * (one - fz) * (potential[i + 1, j, k] - potential[i, j, k])
        + fy * (one - fz) * (potential[i + 1, j + 1, k] - potential[i, j + 1, k])
        + (one - fy) * fz * (potential[i + 1, j, k + 1] - potential[i, j, k + 1])
        + fy * fz * (potential[i + 1, j + 1, k + 1] - potential[i, j + 1, k + 1])) / spacing[0]
    dvy = (
        (one - fx) * (one - fz) * (potential[i, j + 1, k] - potential[i, j, k])
        + fx * (one - fz) * (potential[i + 1, j + 1, k] - potential[i + 1, j, k])
        + (one - fx) * fz * (potential[i, j + 1, k + 1] - potential[i, j, k + 1])
        + fx * fz * (potential[i + 1, j + 1, k + 1] - potential[i + 1, j, k + 1])) / spacing[1]
    dvz = (
        (one - fx) * (one - fy) * (potential[i, j, k + 1] - potential[i, j, k])
        + fx * (one - fy) * (potential[i + 1, j, k + 1] - potential[i + 1, j, k])
        + (one - fx) * fy * (potential[i, j + 1, k + 1] - potential[i, j + 1, k])
        + fx * fy * (potential[i + 1, j + 1, k + 1] - potential[i + 1, j + 1, k])) / spacing[2]
    return wp.vec3(-dvx, -dvy, -dvz)


@wp.kernel
def _field_hit_events_3d(
        mesh: wp.uint64, potential: wp.array3d(dtype=float), grid_origin: wp.vec3,
        grid_spacing: wp.vec3, grid_maximum: wp.vec3,
        origin: wp.array(dtype=wp.vec3), velocity: wp.array(dtype=wp.vec3),
        charge_number: float, fixed_dt: float, max_steps: int, periodic_lateral: int,
        hit_face: wp.array(dtype=int), hit_cosine: wp.array(dtype=float),
        hit_energy: wp.array(dtype=float), termination: wp.array(dtype=wp.int8),
        terminal_position: wp.array(dtype=wp.vec3), terminal_velocity: wp.array(dtype=wp.vec3)):
    particle = wp.tid()
    position = origin[particle]
    speed_vector = velocity[particle]
    for _step in range(max_steps):
        field0 = _trilinear_electric_field_3d(
            potential, position, grid_origin, grid_spacing)
        half_velocity = speed_vector + 0.25 * charge_number * fixed_dt * field0
        next_position = position + fixed_dt * half_velocity
        if periodic_lateral != 0:
            domain_x = grid_maximum[0] - grid_origin[0]
            domain_y = grid_maximum[1] - grid_origin[1]
            for _wrap in range(4):
                if next_position[0] < grid_origin[0]:
                    next_position = wp.vec3(
                        next_position[0] + domain_x, next_position[1], next_position[2])
                elif next_position[0] > grid_maximum[0]:
                    next_position = wp.vec3(
                        next_position[0] - domain_x, next_position[1], next_position[2])
                if next_position[1] < grid_origin[1]:
                    next_position = wp.vec3(
                        next_position[0], next_position[1] + domain_y, next_position[2])
                elif next_position[1] > grid_maximum[1]:
                    next_position = wp.vec3(
                        next_position[0], next_position[1] - domain_y, next_position[2])
        field1 = _trilinear_electric_field_3d(
            potential, next_position, grid_origin, grid_spacing)
        next_velocity = half_velocity + 0.25 * charge_number * fixed_dt * field1
        # The physical Verlet segment is straight in the covering space.  Under periodic lateral
        # boundaries it may split into several in-cell ray segments; query each one so a wrap cannot
        # tunnel through a surface near the opposite boundary.
        segment = fixed_dt * half_velocity
        segment_length = wp.length(segment)
        if segment_length > 0.0:
            direction = segment / segment_length
            query_origin = position
            remaining = segment_length
            travelled = float(0.0)
            for _piece in range(8):
                query_length = remaining
                wrap_axis = wp.int32(-1)
                if periodic_lateral != 0:
                    if direction[0] > 1.0e-12:
                        distance = (grid_maximum[0] - query_origin[0]) / direction[0]
                        if distance >= 0.0 and distance < query_length:
                            query_length = distance
                            wrap_axis = wp.int32(0)
                    elif direction[0] < -1.0e-12:
                        distance = (grid_origin[0] - query_origin[0]) / direction[0]
                        if distance >= 0.0 and distance < query_length:
                            query_length = distance
                            wrap_axis = wp.int32(0)
                    if direction[1] > 1.0e-12:
                        distance = (grid_maximum[1] - query_origin[1]) / direction[1]
                        if distance >= 0.0 and distance < query_length:
                            query_length = distance
                            wrap_axis = wp.int32(1)
                    elif direction[1] < -1.0e-12:
                        distance = (grid_origin[1] - query_origin[1]) / direction[1]
                        if distance >= 0.0 and distance < query_length:
                            query_length = distance
                            wrap_axis = wp.int32(1)
                ray = wp.mesh_query_ray(
                    mesh, query_origin, direction, query_length * 1.000001 + 1.0e-9)
                if ray.result and ray.t <= query_length * 1.000001 + 1.0e-9:
                    fraction = wp.clamp((travelled + ray.t) / segment_length, 0.0, 1.0)
                    impact_velocity = speed_vector + fraction * (next_velocity - speed_vector)
                    impact_speed = wp.length(impact_velocity)
                    normal = ray.normal
                    if wp.dot(normal, impact_velocity) > 0.0:
                        normal = -normal
                    cosine = float(0.0)
                    if impact_speed > 0.0:
                        cosine = -wp.dot(impact_velocity / impact_speed, normal)
                    hit_face[particle] = ray.face
                    hit_cosine[particle] = wp.clamp(cosine, 0.0, 1.0)
                    hit_energy[particle] = wp.dot(impact_velocity, impact_velocity)
                    termination[particle] = wp.int8(1)
                    terminal_position[particle] = position + fraction * segment
                    terminal_velocity[particle] = impact_velocity
                    break
                if wrap_axis < 0:
                    break
                remaining = remaining - query_length
                travelled = travelled + query_length
                boundary_point = query_origin + query_length * direction
                if wrap_axis == 0:
                    wrapped_x = grid_origin[0]
                    if direction[0] < 0.0:
                        wrapped_x = grid_maximum[0]
                    query_origin = wp.vec3(wrapped_x, boundary_point[1], boundary_point[2])
                else:
                    wrapped_y = grid_origin[1]
                    if direction[1] < 0.0:
                        wrapped_y = grid_maximum[1]
                    query_origin = wp.vec3(boundary_point[0], wrapped_y, boundary_point[2])
            if termination[particle] == wp.int8(1):
                break
        position = next_position
        speed_vector = next_velocity
        terminal_position[particle] = position
        terminal_velocity[particle] = speed_vector
        outside = (position[2] < grid_origin[2] or position[2] > grid_maximum[2])
        if periodic_lateral == 0:
            outside = (outside or position[0] < grid_origin[0] or position[1] < grid_origin[1]
                       or position[0] > grid_maximum[0] or position[1] > grid_maximum[1])
        if outside:
            termination[particle] = wp.int8(2)
            break


@dataclass(frozen=True)
class BoundaryTransport3DResult:
    surface_fluxes: SurfaceFluxes
    hit_probability: Mapping[str, float]
    escape_probability: Mapping[str, float]
    truncation_probability: Mapping[str, float]
    transport_model: str
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        object.__setattr__(self, "hit_probability", MappingProxyType(dict(self.hit_probability)))
        object.__setattr__(self, "escape_probability", MappingProxyType(dict(self.escape_probability)))
        object.__setattr__(
            self, "truncation_probability", MappingProxyType(dict(self.truncation_probability)))
        object.__setattr__(self, "known_limitations", tuple(self.known_limitations))


@dataclass(frozen=True)
class BidirectionalFaceSelection3D:
    population: FaceResolvedEnergeticFlux
    face_stderr_m2_s: np.ndarray
    method: np.ndarray
    forward_face_mean_m2_s: np.ndarray
    forward_face_stderr_m2_s: np.ndarray
    adjoint_face_mean_m2_s: np.ndarray
    adjoint_face_stderr_m2_s: np.ndarray
    estimator_consistent: np.ndarray
    method_within_tolerance: np.ndarray
    converged: bool

    def __post_init__(self):
        face_count = self.population.face_count
        for name in (
                "face_stderr_m2_s", "forward_face_mean_m2_s", "forward_face_stderr_m2_s",
                "adjoint_face_mean_m2_s", "adjoint_face_stderr_m2_s"):
            value = np.asarray(getattr(self, name), dtype=float).copy()
            if value.shape != (face_count,) or np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise ValueError("bidirectional face statistics must be finite nonnegative arrays")
            value.setflags(write=False); object.__setattr__(self, name, value)
        method = np.asarray(self.method).astype("U7", copy=True)
        consistent = np.asarray(self.estimator_consistent, dtype=bool).copy()
        within = np.asarray(self.method_within_tolerance, dtype=bool).copy()
        if (method.shape != (face_count,) or np.any(~np.isin(method, ("forward", "adjoint")))
                or consistent.shape != (face_count,) or within.shape != (face_count,)):
            raise ValueError("invalid bidirectional face classifications")
        for value in (method, consistent, within): value.setflags(write=False)
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "estimator_consistent", consistent)
        object.__setattr__(self, "method_within_tolerance", within)
        object.__setattr__(self, "converged", bool(self.converged))


@dataclass(frozen=True)
class BidirectionalBoundaryTransport3DResult:
    transport: BoundaryTransport3DResult
    selection_by_species: Mapping[str, BidirectionalFaceSelection3D]

    def __post_init__(self):
        if not isinstance(self.transport, BoundaryTransport3DResult):
            raise TypeError("transport must be a BoundaryTransport3DResult")
        selection = dict(self.selection_by_species)
        if (set(selection) != set(self.transport.hit_probability)
                or any(not isinstance(item, BidirectionalFaceSelection3D)
                       for item in selection.values())):
            raise ValueError("bidirectional selections must classify every transported species")
        object.__setattr__(self, "selection_by_species", MappingProxyType(selection))


def select_bidirectional_face_events_3d(
        forward_populations, adjoint_populations, *, forward_zero_upper_m2_s,
        absolute_tolerance_m2_s, relative_tolerance=0.05, consistency_sigma=5.0,
        support_sigma=2.0, support_ratio=0.5, method_hint=None,
        require_certification=True):
    """Select independently replicated forward or adjoint event measures per triangle.

    Selection uses measured standard error only. The direct forward estimator also audits adjoint
    support: a precise adjoint zero cannot win when forward transport resolves nonzero current. When
    both directions claim the requested precision, they must agree within the combined uncertainty.
    The selected replicate event measures are averaged without energy/angle binning.
    """
    forward = tuple(forward_populations); adjoint = tuple(adjoint_populations)
    if len(forward) < 4 or len(adjoint) != len(forward):
        raise ValueError("bidirectional selection requires matching sets of at least four replicates")
    first = forward[0]
    if (not isinstance(first, FaceResolvedEnergeticFlux)
            or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in forward + adjoint)
            or any(item.name != first.name or item.face_count != first.face_count
                   for item in forward + adjoint)):
        raise ValueError("bidirectional replicates must describe one matching face-resolved species")
    face_count = first.face_count; replicate_count = len(forward)
    if (not np.isfinite(absolute_tolerance_m2_s) or absolute_tolerance_m2_s < 0.0
            or not np.isfinite(relative_tolerance) or relative_tolerance < 0.0
            or not np.isfinite(consistency_sigma) or consistency_sigma <= 0.0
            or not np.isfinite(support_sigma) or support_sigma <= 0.0
            or not np.isfinite(support_ratio) or not 0.0 < support_ratio < 1.0):
        raise ValueError("invalid bidirectional estimator controls")
    zero_upper = np.asarray(forward_zero_upper_m2_s, dtype=float)
    if (zero_upper.shape != (face_count,) or np.any(~np.isfinite(zero_upper))
            or np.any(zero_upper < 0.0)):
        raise ValueError("forward zero-hit bounds must match the surface mesh")
    forward_replicates = np.stack([item.flux_m2_s for item in forward])
    adjoint_replicates = np.stack([item.flux_m2_s for item in adjoint])
    forward_mean = forward_replicates.mean(axis=0)
    adjoint_mean = adjoint_replicates.mean(axis=0)
    forward_stderr = forward_replicates.std(axis=0, ddof=1) / np.sqrt(replicate_count)
    adjoint_stderr = adjoint_replicates.std(axis=0, ddof=1) / np.sqrt(replicate_count)
    forward_stderr = np.where(forward_mean == 0.0, np.maximum(forward_stderr, zero_upper),
                              forward_stderr)
    forward_allowed = absolute_tolerance_m2_s + relative_tolerance * np.abs(forward_mean)
    adjoint_allowed = absolute_tolerance_m2_s + relative_tolerance * np.abs(adjoint_mean)
    forward_ok = forward_stderr <= forward_allowed
    adjoint_support_unresolved = (
        adjoint_mean + support_sigma * adjoint_stderr
        < support_ratio * np.maximum(forward_mean - support_sigma * forward_stderr, 0.0))
    adjoint_ok = (adjoint_stderr <= adjoint_allowed) & ~adjoint_support_unresolved
    combined = np.hypot(forward_stderr, adjoint_stderr)
    discrepancy_sigma = np.divide(
        np.abs(forward_mean - adjoint_mean), combined,
        out=np.where(forward_mean == adjoint_mean, 0.0, np.inf), where=combined > 0.0)
    consistent = ~(forward_ok & adjoint_ok) | (discrepancy_sigma <= consistency_sigma)
    metric_floor = max(float(absolute_tolerance_m2_s), np.finfo(float).tiny)
    forward_metric = forward_stderr / np.maximum(np.abs(forward_mean), metric_floor)
    adjoint_metric = adjoint_stderr / np.maximum(np.abs(adjoint_mean), metric_floor)
    method = np.where(forward_metric < adjoint_metric, "forward", "adjoint")
    method = np.where(forward_ok & ~adjoint_ok, "forward", method)
    method = np.where(adjoint_ok & ~forward_ok, "adjoint", method)
    if method_hint is not None:
        hint = np.asarray(method_hint).astype("U7")
        if hint.shape != (face_count,) or np.any(~np.isin(hint, ("forward", "adjoint"))):
            raise ValueError("method_hint must select forward or adjoint for every face")
        method = hint
    within = np.where(method == "forward", forward_ok, adjoint_ok)
    if method_hint is not None and not require_certification:
        within = np.ones(face_count, dtype=bool)
        consistent = np.ones(face_count, dtype=bool)
    selected_events = []
    for populations, selected_method in ((forward, "forward"), (adjoint, "adjoint")):
        for population in populations:
            keep = method[population.event_face] == selected_method
            selected_events.append((
                population.event_face[keep],
                population.event_flux_m2_s[keep] / replicate_count,
                population.event_energy_eV[keep],
                population.event_cosine_incidence[keep]))
    population = FaceResolvedEnergeticFlux(
        first.name, face_count,
        np.concatenate([item[0] for item in selected_events]),
        np.concatenate([item[1] for item in selected_events]),
        np.concatenate([item[2] for item in selected_events]),
        np.concatenate([item[3] for item in selected_events]))
    selected_stderr = np.where(method == "forward", forward_stderr, adjoint_stderr)
    return BidirectionalFaceSelection3D(
        population, selected_stderr, method,
        forward_mean, forward_stderr, adjoint_mean, adjoint_stderr,
        consistent, within, bool(np.all(consistent & within)))


def _replace_population_face_events_3d(population, replacement, face_indices):
    """Replace selected triangle events while preserving the untouched sparse measure."""
    selected = np.zeros(population.face_count, dtype=bool)
    selected[np.asarray(face_indices, dtype=int)] = True
    keep = ~selected[population.event_face]
    return FaceResolvedEnergeticFlux(
        population.name, population.face_count,
        np.concatenate((population.event_face[keep], replacement.event_face)),
        np.concatenate((population.event_flux_m2_s[keep], replacement.event_flux_m2_s)),
        np.concatenate((population.event_energy_eV[keep], replacement.event_energy_eV)),
        np.concatenate((population.event_cosine_incidence[keep],
                        replacement.event_cosine_incidence)))


def merge_boundary_transport_results_3d(*results):
    """Merge disjoint species transports without changing any event measure."""
    results = tuple(results)
    if not results or any(not isinstance(item, BoundaryTransport3DResult) for item in results):
        raise ValueError("one or more BoundaryTransport3DResult objects are required")
    neutral = {}; energetic = []; hit = {}; escaped = {}; truncated = {}
    species_seen = set(); models = []; limitations = []
    for result in results:
        species = set(result.hit_probability)
        if (species != set(result.escape_probability)
                or species != set(result.truncation_probability)):
            raise ValueError("transport probability maps must classify identical species")
        if species_seen & species:
            raise ValueError("merged transport results must contain disjoint species")
        species_seen |= species
        for name, value in result.surface_fluxes.neutral_flux_m2_s.items():
            if name in neutral:
                raise ValueError("merged neutral species names must be unique")
            neutral[name] = value
        energetic.extend(result.surface_fluxes.energetic_fluxes)
        hit.update(result.hit_probability); escaped.update(result.escape_probability)
        truncated.update(result.truncation_probability)
        models.append(result.transport_model); limitations.extend(result.known_limitations)
    energetic_names = [item.name for item in energetic]
    if len(set(energetic_names)) != len(energetic_names):
        raise ValueError("merged energetic species names must be unique")
    return BoundaryTransport3DResult(
        surface_fluxes=SurfaceFluxes(neutral, tuple(energetic)),
        hit_probability=hit, escape_probability=escaped,
        truncation_probability=truncated,
        transport_model=" + ".join(dict.fromkeys(models)),
        known_limitations=tuple(dict.fromkeys(limitations)))


def estimate_diffuse_form_factors_3d(
        verts, faces, centroids, gas_normals, *, rays_per_face=64, seed=0,
        domain_size=None, periodic_lateral=False, ray_offset=1e-5, device=None):
    """Estimate deterministic diffuse face exchange and classify every emitted ray.

    A scrambled Sobol hemisphere rule is Cranley-shifted per source face. The estimator produces
    geometric form factors only; sticking and chemistry enter the separate conservative radiosity
    solve. ``periodic_lateral`` uses the same periodic-cell ray geometry as the legacy trench engine,
    with the top remaining an open escape boundary.
    """
    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=int)
    centroids = np.asarray(centroids, dtype=float)
    normals = np.asarray(gas_normals, dtype=float)
    if (verts.ndim != 2 or verts.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3
            or centroids.shape != (faces.shape[0], 3) or normals.shape != centroids.shape
            or np.any(~np.isfinite(verts)) or np.any(~np.isfinite(centroids))
            or np.any(~np.isfinite(normals)) or np.any(faces < 0)
            or np.any(faces >= len(verts)) or ray_offset <= 0.0):
        raise ValueError("invalid mesh or gas-normal input for diffuse form factors")
    normal_length = np.linalg.norm(normals, axis=1)
    if not np.allclose(normal_length, 1.0, rtol=0.0, atol=2e-6):
        raise ValueError("gas normals must be unit length")
    if int(rays_per_face) != rays_per_face:
        raise ValueError("rays_per_face must be an integer")
    rays_per_face = int(rays_per_face)
    if rays_per_face <= 0 or rays_per_face & (rays_per_face - 1):
        raise ValueError("rays_per_face must be a positive power of two")
    if domain_size is None:
        domain = np.maximum(np.ptp(verts, axis=0), 1.0)
    else:
        domain = np.asarray(domain_size, dtype=float)
    if domain.shape != (3,) or np.any(~np.isfinite(domain)) or np.any(domain <= 0.0):
        raise ValueError("domain_size must contain three positive lengths")
    if periodic_lateral and (
            np.min(verts[:, :2]) < -1e-7
            or np.any(np.max(verts[:, :2], axis=0) > domain[:2] + 1e-7)
            or np.max(verts[:, 2]) > domain[2] + 1e-7):
        raise ValueError("periodic mesh must lie inside [0, domain_size]")

    base = qmc.Sobol(2, scramble=True, seed=int(seed)).random_base2(
        int(np.log2(rays_per_face)))
    face_index = np.arange(faces.shape[0], dtype=float)[:, None]
    shift = np.mod(face_index * np.array([[0.6180339887498949, 0.4142135623730950]]), 1.0)
    sample = np.mod(base[None, :, :] + shift[:, None, :], 1.0)
    normal = np.repeat(normals, rays_per_face, axis=0)
    sample = sample.reshape(-1, 2)
    tangent_seed = np.where(
        np.abs(normal[:, :1]) > 0.9,
        np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    tangent = np.cross(tangent_seed, normal)
    tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
    bitangent = np.cross(normal, tangent)
    cosine = np.sqrt(sample[:, 0])
    sine = np.sqrt(1.0 - sample[:, 0])
    azimuth = 2.0 * np.pi * sample[:, 1]
    direction = ((sine * np.cos(azimuth))[:, None] * tangent
                 + (sine * np.sin(azimuth))[:, None] * bitangent
                 + cosine[:, None] * normal)
    source = np.repeat(np.arange(faces.shape[0]), rays_per_face)
    origin = centroids[source] + float(ray_offset) * normal

    selected_device = DEVICE if device is None else str(device)
    if selected_device.startswith("warp:"):
        selected_device = selected_device.split(":", 1)[1]
    ensure_writable_warp_cache(wp)
    mesh = wp.Mesh(
        points=wp.array(verts.astype(np.float32), dtype=wp.vec3, device=selected_device),
        indices=wp.array(faces.astype(np.int32).ravel(), dtype=wp.int32, device=selected_device))
    hit_wp = wp.full(source.size, -1, dtype=wp.int32, device=selected_device)
    wp.launch(
        _diffuse_form_factor_events_3d, dim=source.size, device=selected_device,
        inputs=[mesh.id,
                wp.array(origin.astype(np.float32), dtype=wp.vec3, device=selected_device),
                wp.array(direction.astype(np.float32), dtype=wp.vec3, device=selected_device),
                float(domain[0]), float(domain[1]), float(domain[2]),
                int(bool(periodic_lateral)), hit_wp])
    hit = hit_wp.numpy().astype(int)
    escaped = hit < 0
    escape_fraction = np.bincount(
        source[escaped], minlength=faces.shape[0]).astype(float) / rays_per_face
    valid_source = source[~escaped]
    valid_target = hit[~escaped]
    pair = valid_source.astype(np.int64) * faces.shape[0] + valid_target
    unique, count = np.unique(pair, return_counts=True)
    source_face = (unique // faces.shape[0]).astype(int)
    target_face = (unique % faces.shape[0]).astype(int)
    return DiffuseFormFactors3D(
        faces.shape[0], source_face, target_face, count.astype(float) / rays_per_face,
        escape_fraction, rays_per_face)


def gather_boundary_state_ballistic_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str], verts, faces, areas,
        centroids, gas_normals, *, source_bounds, source_z, mesh_length_unit_m=1e-6,
        mesh_origin_m=(0.0, 0.0, 0.0), face_quadrature_points=1,
        periodic_lateral=False, domain_size=None, ray_offset=1e-5, device=None):
    """Deterministically gather collisionless boundary flux onto every visible triangle.

    For boundary direction ``d`` (pointing from the horizontal source plane toward the surface),
    conservation of projected area gives

    ``Gamma_face = Gamma_plane * w * max(0, -d.n) / abs(d.z)``.

    A reverse ray from each triangle quadrature point supplies only the geometric visibility.  This
    is the adjoint form of absorbing first-hit transport: unlike a forward particle tally it gives
    every face a deterministic local flux and cannot imprint zero-count triangles into an evolving
    interface.  The boundary's discrete energy-angle measure is retained exactly.
    """
    verts = np.asarray(verts, dtype=float); faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float); centroids = np.asarray(centroids, dtype=float)
    normals = np.asarray(gas_normals, dtype=float); bounds = np.asarray(source_bounds, dtype=float)
    if (verts.ndim != 2 or verts.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3
            or areas.shape != (faces.shape[0],) or centroids.shape != (faces.shape[0], 3)
            or normals.shape != centroids.shape or np.any(~np.isfinite(verts))
            or np.any(~np.isfinite(centroids)) or np.any(~np.isfinite(normals))
            or np.any(faces < 0) or np.any(faces >= len(verts)) or np.any(areas <= 0.0)
            or bounds.shape != (4,) or np.any(~np.isfinite(bounds))
            or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]
            or not np.isfinite(source_z) or not np.isfinite(ray_offset) or ray_offset <= 0.0):
        raise ValueError("invalid deterministic face-gather geometry")
    geometric_areas = 0.5 * np.linalg.norm(np.cross(
        verts[faces[:, 1]] - verts[faces[:, 0]],
        verts[faces[:, 2]] - verts[faces[:, 0]]), axis=1)
    if not np.allclose(areas, geometric_areas, rtol=1e-7, atol=0.0):
        raise ValueError("triangle areas must match the supplied mesh geometry")
    if not np.allclose(np.linalg.norm(normals, axis=1), 1.0, rtol=0.0, atol=2e-6):
        raise ValueError("gas normals must be unit length")
    role = dict(species_role); names = {item.name for item in boundary.species}
    if set(role) != names:
        raise ValueError("species_role must classify every and only boundary species")
    allowed_roles = {"neutral_reactant", "energetic_bombardment", "charge_carrier"}
    if any(value not in allowed_roles for value in role.values()):
        raise ValueError(f"species roles must be one of {sorted(allowed_roles)}")
    if any(item.position_m is not None for item in boundary.species):
        raise ValueError("face gather currently requires a spatially uniform boundary state")
    if int(face_quadrature_points) != face_quadrature_points or int(face_quadrature_points) not in (1, 3):
        raise ValueError("face_quadrature_points must be 1 or 3")
    face_quadrature_points = int(face_quadrature_points)
    origin_m = np.asarray(mesh_origin_m, dtype=float)
    if (origin_m.shape != (3,) or np.any(~np.isfinite(origin_m))
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("invalid mesh SI coordinate mapping")
    mapped_reference = origin_m[2] + float(source_z) * float(mesh_length_unit_m)
    if not np.isclose(
            mapped_reference, boundary.reference_plane_m, rtol=0.0,
            atol=max(1e-15, 1e-9 * float(mesh_length_unit_m))):
        raise ValueError("mesh source plane does not match PlasmaBoundaryState.reference_plane_m")
    if domain_size is None:
        domain = np.maximum(np.ptp(verts, axis=0), 1.0)
    else:
        domain = np.asarray(domain_size, dtype=float)
    if domain.shape != (3,) or np.any(~np.isfinite(domain)) or np.any(domain <= 0.0):
        raise ValueError("domain_size must contain three positive lengths")

    if face_quadrature_points == 1:
        points = centroids[:, None, :]
    else:
        barycentric = np.array([
            [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
        ])
        points = np.einsum("qv,fvc->fqc", barycentric, verts[faces])
    selected_device = DEVICE if device is None else str(device)
    if selected_device.startswith("warp:"):
        selected_device = selected_device.split(":", 1)[1]
    ensure_writable_warp_cache(wp)
    mesh = wp.Mesh(
        points=wp.array(verts.astype(np.float32), dtype=wp.vec3, device=selected_device),
        indices=wp.array(faces.astype(np.int32).ravel(), dtype=wp.int32, device=selected_device))
    source_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
    neutral_flux = {}; energetic_flux = []; hit_probability = {}; escape_probability = {}
    face_count = faces.shape[0]; point_count = face_quadrature_points
    for species in boundary.species:
        velocity = np.asarray(species.velocity_sqrt_eV, dtype=float).copy()
        velocity[:, 2] *= -1.0
        speed = np.linalg.norm(velocity, axis=1)
        direction = velocity / speed[:, None]
        if np.any(direction[:, 2] >= 0.0):
            raise ValueError("face gather requires boundary velocities directed toward the surface")
        incidence_cosine = np.clip(
            -np.einsum("sd,fd->sf", direction, normals), 0.0, 1.0)
        projection = incidence_cosine / (-direction[:, 2, None])
        normalized_gathered = np.zeros((direction.shape[0], face_count))
        for sample_index, incident_direction in enumerate(direction):
            reverse = -incident_direction
            origin = points.reshape(-1, 3) + np.repeat(
                float(ray_offset) * normals, point_count, axis=0)
            ray_direction = np.broadcast_to(reverse, origin.shape).copy()
            hit_wp = wp.full(origin.shape[0], -1, dtype=wp.int32, device=selected_device)
            wp.launch(
                _diffuse_form_factor_events_3d, dim=origin.shape[0], device=selected_device,
                inputs=[mesh.id,
                        wp.array(origin.astype(np.float32), dtype=wp.vec3, device=selected_device),
                        wp.array(ray_direction.astype(np.float32), dtype=wp.vec3, device=selected_device),
                        float(domain[0]), float(domain[1]), float(domain[2]),
                        int(bool(periodic_lateral)), hit_wp])
            visible = hit_wp.numpy().reshape(face_count, point_count) < 0
            if not periodic_lateral:
                travel = (float(source_z) - points[:, :, 2]) / reverse[2]
                source_point = points[:, :, :2] + travel[:, :, None] * reverse[None, None, :2]
                visible &= ((travel >= 0.0)
                            & (source_point[:, :, 0] >= bounds[0])
                            & (source_point[:, :, 0] <= bounds[1])
                            & (source_point[:, :, 1] >= bounds[2])
                            & (source_point[:, :, 1] <= bounds[3]))
            normalized_gathered[sample_index] = (
                species.weight[sample_index] * projection[sample_index]
                * visible.mean(axis=1))
        sample_probability = np.einsum("sf,f->s", normalized_gathered, areas) / source_area
        probability = float(sample_probability.sum())
        if periodic_lateral:
            # An opaque periodic feature cell has no lateral loss: every downward source ray must
            # land. Enforce this independently for every discrete energy-angle atom so geometric
            # quadrature error cannot distort the boundary distribution while conserving its total.
            expected_probability = np.asarray(species.weight, dtype=float)
            positive_weight = expected_probability > 0.0
            if (np.any(sample_probability[positive_weight] <= 0.0)
                    or np.any(~np.isfinite(sample_probability[positive_weight]))):
                raise RuntimeError(
                    f"periodic face visibility has no landed measure for {species.name!r}")
            scale = np.ones(expected_probability.shape, dtype=float)
            scale[positive_weight] = (
                expected_probability[positive_weight]
                / sample_probability[positive_weight])
            normalized_gathered *= scale[:, None]
            sample_probability = expected_probability.copy()
            probability = float(sample_probability.sum())
        elif not -5e-13 <= probability <= 1.0 + 5e-13:
            raise RuntimeError(
                f"face visibility quadrature violates projected-area conservation for "
                f"{species.name!r}: landed probability={probability:.8g}, "
                f"per-sample={sample_probability.tolist()}; refine triangle visibility")
        else:
            probability = min(max(probability, 0.0), 1.0)
        gathered = species.flux_m2_s * normalized_gathered
        hit_probability[species.name] = probability
        escape_probability[species.name] = 1.0 - probability
        if role[species.name] == "neutral_reactant":
            neutral_flux[species.name] = gathered.sum(axis=0)
        else:
            event_sample, event_face = np.where(gathered > 0.0)
            energetic_flux.append(FaceResolvedEnergeticFlux(
                species.name, face_count, event_face,
                gathered[event_sample, event_face],
                species.kinetic_energy_eV[event_sample],
                incidence_cosine[event_sample, event_face]))
    return BoundaryTransport3DResult(
        surface_fluxes=SurfaceFluxes(neutral_flux, tuple(energetic_flux)),
        hit_probability=hit_probability, escape_probability=escape_probability,
        truncation_probability={name: 0.0 for name in hit_probability},
        transport_model="collisionless_deterministic_face_gather_3d",
        known_limitations=(
            "no intra-feature electric-field trajectory coupling",
            "no surface reflection or neutral re-emission",
            "no spatially varying boundary density",
            "triangle visibility quadrature requires refinement at partial occlusion",
            "float32 triangle-ray intersection",
        ))


def trace_boundary_state_first_hit_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str], verts, faces, areas, *,
        source_bounds, source_z, mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        n_position=256, seed=0, max_distance=None, periodic_lateral=False,
        domain_size=None, device=None):
    """Transport a spatially uniform boundary state to exact triangle-hit events.

    ``species_role`` is a physical input mapping each boundary species to
    ``"neutral_reactant"``, ``"energetic_bombardment"``, or ``"charge_carrier"``; no species name
    selects a formula. Charge carriers contribute charged-particle hit events for current deposition,
    but a downstream surface mechanism must explicitly select which energetic populations drive its
    chemistry. Mesh coordinates may use any length unit declared by ``mesh_length_unit_m``. The mapped
    source plane must equal the boundary's SI reference plane, preventing a silent geometry/boundary
    unit mismatch.

    Boundary velocity quadrature is retained exactly. Scrambled Sobol points integrate only the uniform
    source-plane position, so changing ``n_position`` cannot change the physical energy-angle law.
    """
    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float)
    if (verts.ndim != 2 or verts.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3
            or areas.shape != (faces.shape[0],) or np.any(~np.isfinite(verts))
            or np.any(faces < 0) or np.any(faces >= len(verts))
            or np.any(~np.isfinite(areas)) or np.any(areas <= 0.0)):
        raise ValueError("invalid triangle mesh")
    edge_a = verts[faces[:, 1]] - verts[faces[:, 0]]
    edge_b = verts[faces[:, 2]] - verts[faces[:, 0]]
    geometric_areas = 0.5 * np.linalg.norm(np.cross(edge_a, edge_b), axis=1)
    if (np.any(geometric_areas <= 0.0)
            or not np.allclose(areas, geometric_areas, rtol=1e-7, atol=0.0)):
        raise ValueError("triangle areas must match the supplied mesh geometry")
    bounds = np.asarray(source_bounds, dtype=float)
    if (bounds.shape != (4,) or np.any(~np.isfinite(bounds))
            or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]):
        raise ValueError("source_bounds must be (x_min, x_max, y_min, y_max)")
    if not np.isfinite(source_z) or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0:
        raise ValueError("source coordinate and mesh length unit must be finite and physical")
    origin_m = np.asarray(mesh_origin_m, dtype=float)
    if origin_m.shape != (3,) or np.any(~np.isfinite(origin_m)):
        raise ValueError("mesh_origin_m must contain three finite SI coordinates")
    mapped_reference = origin_m[2] + float(source_z) * float(mesh_length_unit_m)
    reference_tolerance = max(1e-15, 1e-9 * float(mesh_length_unit_m))
    if not np.isclose(
            mapped_reference, boundary.reference_plane_m, rtol=0.0, atol=reference_tolerance):
        raise ValueError(
            "mesh source plane does not match PlasmaBoundaryState.reference_plane_m")
    if int(n_position) != n_position:
        raise ValueError("n_position must be an integer")
    n_position = int(n_position)
    if n_position <= 0 or n_position & (n_position - 1):
        raise ValueError("n_position must be a positive power of two for a balanced Sobol rule")
    role = dict(species_role)
    names = {item.name for item in boundary.species}
    if set(role) != names:
        raise ValueError("species_role must classify every and only boundary species")
    allowed_roles = {"neutral_reactant", "energetic_bombardment", "charge_carrier"}
    if any(value not in allowed_roles for value in role.values()):
        raise ValueError(f"species roles must be one of {sorted(allowed_roles)}")
    if any(item.position_m is not None for item in boundary.species):
        raise ValueError(
            "first-hit 3-D transport currently requires a spatially uniform boundary state")

    source_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
    sampler = qmc.Sobol(2, scramble=True, seed=int(seed))
    positions = sampler.random_base2(int(np.log2(n_position)))
    x_position = bounds[0] + positions[:, 0] * (bounds[1] - bounds[0])
    y_position = bounds[2] + positions[:, 1] * (bounds[3] - bounds[2])

    selected_device = DEVICE if device is None else str(device)
    if selected_device.startswith("warp:"):
        selected_device = selected_device.split(":", 1)[1]
    ensure_writable_warp_cache(wp)
    mesh = wp.Mesh(
        points=wp.array(verts.astype(np.float32), dtype=wp.vec3, device=selected_device),
        indices=wp.array(faces.astype(np.int32).ravel(), dtype=wp.int32, device=selected_device))
    if max_distance is None:
        source_corners = np.array([
            [bounds[0], bounds[2], source_z], [bounds[0], bounds[3], source_z],
            [bounds[1], bounds[2], source_z], [bounds[1], bounds[3], source_z],
        ])
        max_distance = 1.01 * float(np.max(np.linalg.norm(
            source_corners[:, None, :] - verts[None, :, :], axis=2)))
    if not np.isfinite(max_distance) or max_distance <= 0.0:
        raise ValueError("max_distance must be positive and finite")
    if periodic_lateral:
        domain = np.asarray(domain_size, dtype=float)
        if (domain.shape != (3,) or np.any(~np.isfinite(domain)) or np.any(domain <= 0.0)
                or np.min(verts[:, :2]) < -1e-7
                or np.any(np.max(verts[:, :2], axis=0) > domain[:2] + 1e-7)
                or np.max(verts[:, 2]) > domain[2] + 1e-7):
            raise ValueError("periodic first-hit transport requires a containing domain_size")
    else:
        domain = np.ones(3)

    neutral_flux = {}; energetic_flux = []; hit_probability = {}; escape_probability = {}
    for species in boundary.species:
        sample_count = species.velocity_sqrt_eV.shape[0]
        ray_count = sample_count * n_position
        origin = np.column_stack((
            np.tile(x_position, sample_count), np.tile(y_position, sample_count),
            np.full(ray_count, float(source_z))))
        velocity = np.repeat(species.velocity_sqrt_eV, n_position, axis=0)
        speed = np.linalg.norm(velocity, axis=1)
        if np.any(speed <= 0.0):
            raise ValueError(f"species {species.name!r} contains a zero-speed incident sample")
        direction = velocity / speed[:, None]
        direction[:, 2] *= -1.0
        physical_weight = np.repeat(
            species.weight / n_position, n_position)
        energy = np.repeat(species.kinetic_energy_eV, n_position)

        hit_face_wp = wp.full(ray_count, -1, dtype=wp.int32, device=selected_device)
        hit_cosine_wp = wp.zeros(ray_count, dtype=float, device=selected_device)
        origin_wp = wp.array(origin.astype(np.float32), dtype=wp.vec3, device=selected_device)
        direction_wp = wp.array(direction.astype(np.float32), dtype=wp.vec3, device=selected_device)
        if periodic_lateral:
            wp.launch(
                _periodic_first_hit_events_3d, dim=ray_count, device=selected_device,
                inputs=[mesh.id, origin_wp, direction_wp,
                        float(domain[0]), float(domain[1]), float(domain[2]),
                        hit_face_wp, hit_cosine_wp])
        else:
            wp.launch(
                _first_hit_events_3d, dim=ray_count, device=selected_device,
                inputs=[mesh.id, origin_wp, direction_wp, float(max_distance),
                        hit_face_wp, hit_cosine_wp])
        hit_face = hit_face_wp.numpy().astype(int)
        hit_cosine = hit_cosine_wp.numpy().astype(float)
        hit = hit_face >= 0
        hit_probability[species.name] = float(physical_weight[hit].sum())
        escape_probability[species.name] = float(physical_weight[~hit].sum())
        event_flux = (species.flux_m2_s * source_area * physical_weight[hit]
                      / areas[hit_face[hit]])
        if role[species.name] == "neutral_reactant":
            neutral_flux[species.name] = np.bincount(
                hit_face[hit], weights=event_flux, minlength=faces.shape[0])
        else:
            energetic_flux.append(FaceResolvedEnergeticFlux(
                species.name, faces.shape[0], hit_face[hit], event_flux, energy[hit],
                hit_cosine[hit]))

    return BoundaryTransport3DResult(
        surface_fluxes=SurfaceFluxes(neutral_flux, tuple(energetic_flux)),
        hit_probability=hit_probability,
        escape_probability=escape_probability,
        truncation_probability={name: 0.0 for name in hit_probability},
        transport_model=("collisionless_absorbing_first_hit_3d_periodic_cell"
                         if periodic_lateral else "collisionless_absorbing_first_hit_3d"),
        known_limitations=(
            "no intra-feature electric-field trajectory coupling",
            "no surface reflection or neutral re-emission",
            "no spatially varying boundary density",
            "float32 triangle-ray intersection",
        ))


def trace_boundary_state_field_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str], verts, faces, areas, *,
        source_bounds, source_z, nodal_potential_v, potential_origin, potential_spacing,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0), n_position=256,
        seed=0, fixed_dt=0.01, max_steps=10000, allow_truncation=False,
        phase_space_log2_samples=None, periodic_lateral=False, device=None):
    """Trace collisionless species through a supplied 3-D nodal electrostatic potential.

    Velocity coordinates retain the ``sqrt(eV)`` convention of ``PlasmaBoundaryState``. With mesh
    coordinates as the spatial variable, the scaled Hamiltonian equations are ``dx/dtau=v`` and
    ``dv/dtau=qE/2``. A fixed-step velocity-Verlet map is used so timestep refinement and time reversal
    are meaningful. The potential is an input here; this function does not claim self-consistent charging.
    """
    verts = np.asarray(verts, dtype=float); faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float); potential = np.asarray(nodal_potential_v, dtype=float)
    bounds = np.asarray(source_bounds, dtype=float); grid_origin = np.asarray(potential_origin, dtype=float)
    if (verts.ndim != 2 or verts.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3
            or areas.shape != (faces.shape[0],) or np.any(~np.isfinite(verts))
            or np.any(faces < 0) or np.any(faces >= len(verts)) or np.any(areas <= 0.0)):
        raise ValueError("invalid triangle mesh")
    edge_a = verts[faces[:, 1]] - verts[faces[:, 0]]
    edge_b = verts[faces[:, 2]] - verts[faces[:, 0]]
    geometric_areas = 0.5 * np.linalg.norm(np.cross(edge_a, edge_b), axis=1)
    if not np.allclose(areas, geometric_areas, rtol=1e-7, atol=0.0):
        raise ValueError("triangle areas must match the supplied mesh geometry")
    if (bounds.shape != (4,) or bounds[1] <= bounds[0] or bounds[3] <= bounds[2]
            or np.any(~np.isfinite(bounds)) or not np.isfinite(source_z)):
        raise ValueError("invalid source plane")
    grid_spacing = np.asarray(potential_spacing, dtype=float)
    if grid_spacing.ndim == 0:
        grid_spacing = np.full(3, float(grid_spacing))
    if (potential.ndim != 3 or min(potential.shape) < 2 or np.any(~np.isfinite(potential))
            or grid_origin.shape != (3,) or np.any(~np.isfinite(grid_origin))
            or grid_spacing.shape != (3,) or np.any(~np.isfinite(grid_spacing))
            or np.any(grid_spacing <= 0.0)):
        raise ValueError("invalid nodal potential grid")
    grid_maximum = grid_origin + (np.asarray(potential.shape) - 1) * grid_spacing
    source_corners = np.array([
        [bounds[0], bounds[2], source_z], [bounds[0], bounds[3], source_z],
        [bounds[1], bounds[2], source_z], [bounds[1], bounds[3], source_z]])
    points = np.vstack((verts, source_corners))
    tolerance = 1e-7 * max(float(np.max(grid_spacing)), 1.0)
    if np.any(points < grid_origin - tolerance) or np.any(points > grid_maximum + tolerance):
        raise ValueError("mesh and source plane must lie inside the nodal potential grid")
    if periodic_lateral and not np.allclose(
            bounds, (grid_origin[0], grid_maximum[0], grid_origin[1], grid_maximum[1]),
            rtol=0.0, atol=tolerance):
        raise ValueError(
            "periodic field transport requires source bounds equal to the lateral potential domain")
    origin_m = np.asarray(mesh_origin_m, dtype=float)
    if (origin_m.shape != (3,) or np.any(~np.isfinite(origin_m))
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("invalid mesh SI coordinate mapping")
    mapped_reference = origin_m[2] + float(source_z) * float(mesh_length_unit_m)
    if not np.isclose(
            mapped_reference, boundary.reference_plane_m, rtol=0.0,
            atol=max(1e-15, 1e-9 * float(mesh_length_unit_m))):
        raise ValueError("mesh source plane does not match PlasmaBoundaryState.reference_plane_m")
    if phase_space_log2_samples is None:
        if int(n_position) != n_position:
            raise ValueError("n_position must be an integer")
        n_position = int(n_position)
        if n_position <= 0 or n_position & (n_position - 1):
            raise ValueError("n_position must be a positive power of two for a balanced Sobol rule")
    elif (int(phase_space_log2_samples) != phase_space_log2_samples
          or phase_space_log2_samples < 0):
        raise ValueError("phase_space_log2_samples must be a nonnegative integer")
    if not np.isfinite(fixed_dt) or fixed_dt <= 0.0 or int(max_steps) != max_steps or max_steps <= 0:
        raise ValueError("fixed_dt and max_steps must be positive")
    role = dict(species_role); names = {item.name for item in boundary.species}
    if set(role) != names:
        raise ValueError("species_role must classify every and only boundary species")
    allowed_roles = {"neutral_reactant", "energetic_bombardment", "charge_carrier"}
    if any(value not in allowed_roles for value in role.values()):
        raise ValueError(f"species roles must be one of {sorted(allowed_roles)}")
    if any(item.position_m is not None for item in boundary.species):
        raise ValueError("field 3-D transport currently requires a spatially uniform boundary state")

    if phase_space_log2_samples is None:
        sampler = qmc.Sobol(2, scramble=True, seed=int(seed))
        positions = sampler.random_base2(int(np.log2(n_position)))
        x_position = bounds[0] + positions[:, 0] * (bounds[1] - bounds[0])
        y_position = bounds[2] + positions[:, 1] * (bounds[3] - bounds[2])
    source_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
    selected_device = DEVICE if device is None else str(device)
    if selected_device.startswith("warp:"):
        selected_device = selected_device.split(":", 1)[1]
    ensure_writable_warp_cache(wp)
    mesh = wp.Mesh(
        points=wp.array(verts.astype(np.float32), dtype=wp.vec3, device=selected_device),
        indices=wp.array(faces.astype(np.int32).ravel(), dtype=wp.int32, device=selected_device))
    potential_wp = wp.array(
        np.ascontiguousarray(potential.astype(np.float32)), dtype=float, device=selected_device)

    neutral_flux = {}; energetic_flux = []
    hit_probability = {}; escape_probability = {}; truncation_probability = {}
    for species in boundary.species:
        if phase_space_log2_samples is None:
            sample_count = species.velocity_sqrt_eV.shape[0]
            ray_count = sample_count * n_position
            origin = np.column_stack((
                np.tile(x_position, sample_count), np.tile(y_position, sample_count),
                np.full(ray_count, float(source_z))))
            velocity = np.repeat(species.velocity_sqrt_eV, n_position, axis=0)
            physical_weight = np.repeat(species.weight / n_position, n_position)
        else:
            density_dimension = species.flux_sampling_dimension
            sampler = qmc.Sobol(
                2 + density_dimension, scramble=True, seed=int(seed))
            phase_space = sampler.random_base2(int(phase_space_log2_samples))
            ray_count = phase_space.shape[0]
            origin = np.column_stack((
                bounds[0] + phase_space[:, 0] * (bounds[1] - bounds[0]),
                bounds[2] + phase_space[:, 1] * (bounds[3] - bounds[2]),
                np.full(ray_count, float(source_z))))
            velocity = species.sample_flux_velocity(phase_space[:, 2:])
            physical_weight = np.full(ray_count, 1.0 / ray_count)
        velocity[:, 2] *= -1.0
        if np.any(np.linalg.norm(velocity, axis=1) <= 0.0):
            raise ValueError(f"species {species.name!r} contains a zero-speed incident sample")
        hit_face_wp = wp.full(ray_count, -1, dtype=wp.int32, device=selected_device)
        hit_cosine_wp = wp.zeros(ray_count, dtype=float, device=selected_device)
        hit_energy_wp = wp.zeros(ray_count, dtype=float, device=selected_device)
        termination_wp = wp.zeros(ray_count, dtype=wp.int8, device=selected_device)
        terminal_position_wp = wp.zeros(ray_count, dtype=wp.vec3, device=selected_device)
        terminal_velocity_wp = wp.zeros(ray_count, dtype=wp.vec3, device=selected_device)
        wp.launch(
            _field_hit_events_3d, dim=ray_count, device=selected_device,
            inputs=[mesh.id, potential_wp, wp.vec3(*grid_origin), wp.vec3(*grid_spacing),
                    wp.vec3(*grid_maximum),
                    wp.array(origin.astype(np.float32), dtype=wp.vec3, device=selected_device),
                    wp.array(velocity.astype(np.float32), dtype=wp.vec3, device=selected_device),
                    float(species.charge_number), float(fixed_dt), int(max_steps),
                    int(bool(periodic_lateral)),
                    hit_face_wp, hit_cosine_wp, hit_energy_wp, termination_wp,
                    terminal_position_wp, terminal_velocity_wp])
        hit_face = hit_face_wp.numpy().astype(int)
        termination = termination_wp.numpy().astype(int)
        hit = termination == 1; escaped = termination == 2; truncated = termination == 0
        hit_cosine = hit_cosine_wp.numpy().astype(float)
        hit_energy = hit_energy_wp.numpy().astype(float)
        hit_probability[species.name] = float(physical_weight[hit].sum())
        escape_probability[species.name] = float(physical_weight[escaped].sum())
        truncation_probability[species.name] = float(physical_weight[truncated].sum())
        if truncation_probability[species.name] > 0.0 and not allow_truncation:
            raise RuntimeError(
                f"3-D trajectories for {species.name!r} exhausted max_steps; "
                "increase the physical time horizon or explicitly allow diagnostic truncation")
        event_flux = (species.flux_m2_s * source_area * physical_weight[hit]
                      / areas[hit_face[hit]])
        if role[species.name] == "neutral_reactant":
            neutral_flux[species.name] = np.bincount(
                hit_face[hit], weights=event_flux, minlength=faces.shape[0])
        else:
            energetic_flux.append(FaceResolvedEnergeticFlux(
                species.name, faces.shape[0], hit_face[hit], event_flux,
                hit_energy[hit], hit_cosine[hit]))
    return BoundaryTransport3DResult(
        surface_fluxes=SurfaceFluxes(neutral_flux, tuple(energetic_flux)),
        hit_probability=hit_probability, escape_probability=escape_probability,
        truncation_probability=truncation_probability,
        transport_model=(
            ("collisionless_fixed_step_nodal_field_3d"
             if phase_space_log2_samples is None
             else "collisionless_fixed_step_nodal_field_joint_qmc_3d")
            + ("_periodic_cell" if periodic_lateral else "")),
        known_limitations=(
            "nodal potential is supplied rather than self-consistently charged",
            "no surface reflection or neutral re-emission",
            "no spatially varying boundary density",
            "float32 field integration and triangle-ray intersection",
        ) + (() if phase_space_log2_samples is None else (
            "joint scrambled-Sobol phase-space quadrature requires replicate/refinement error control",
        )))


def gather_boundary_state_field_adjoint_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str],
        verts, faces, areas, centroids, gas_normals, *, source_bounds, source_z,
        nodal_potential_v, potential_origin, potential_spacing,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        face_quadrature_points=3, ray_offset=1e-5, fixed_dt=0.01,
        max_steps=10000, periodic_lateral=False, proposal_by_species=None,
        proposal_frame_by_species="surface_local", face_position_seed=None,
        gather_face_indices=None, device=None):
    """Gather charged boundary flux on every triangle with a reversible Liouville adjoint.

    A boundary velocity proposal may be interpreted in the triangle's local two-tangent/inward-normal
    frame or in the global source frame. Its exact time reverse is integrated through the same
    fixed-step nodal Hamiltonian map used by forward transport. Trajectories that reach the plasma
    plane are scored by the physical boundary density and
    ``v_normal_surface / v_normal_boundary`` Jacobian. The proposal changes variance only; it does
    not replace the declared plasma density.
    """
    verts = np.asarray(verts, dtype=float); faces = np.asarray(faces, dtype=int)
    areas = np.asarray(areas, dtype=float); centroids = np.asarray(centroids, dtype=float)
    normals = np.asarray(gas_normals, dtype=float); potential = np.asarray(nodal_potential_v, dtype=float)
    bounds = np.asarray(source_bounds, dtype=float); grid_origin = np.asarray(potential_origin, dtype=float)
    grid_spacing = np.asarray(potential_spacing, dtype=float)
    if grid_spacing.ndim == 0:
        grid_spacing = np.full(3, float(grid_spacing))
    if (verts.ndim != 2 or verts.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3
            or areas.shape != (len(faces),) or centroids.shape != (len(faces), 3)
            or normals.shape != centroids.shape or np.any(~np.isfinite(verts))
            or np.any(~np.isfinite(centroids)) or np.any(~np.isfinite(normals))
            or not np.allclose(np.linalg.norm(normals, axis=1), 1.0, rtol=0.0, atol=2e-6)
            or np.any(~np.isfinite(areas)) or np.any(areas <= 0.0)
            or potential.ndim != 3 or min(potential.shape) < 2
            or grid_origin.shape != (3,) or grid_spacing.shape != (3,)
            or np.any(grid_spacing <= 0.0) or np.any(~np.isfinite(potential))
            or bounds.shape != (4,) or np.any(~np.isfinite(bounds))
            or bounds[0] >= bounds[1] or bounds[2] >= bounds[3] or not np.isfinite(source_z)
            or int(face_quadrature_points) != face_quadrature_points
            or face_quadrature_points <= 0
            or (face_position_seed is None and int(face_quadrature_points) not in (1, 3, 7))
            or (face_position_seed is not None
                and int(face_position_seed) != face_position_seed)
            or not np.isfinite(ray_offset) or ray_offset <= 0.0
            or not np.isfinite(fixed_dt) or fixed_dt <= 0.0
            or int(max_steps) != max_steps or max_steps <= 0):
        raise ValueError("invalid adjoint field-gather inputs")
    geometric_areas = 0.5 * np.linalg.norm(np.cross(
        verts[faces[:, 1]] - verts[faces[:, 0]],
        verts[faces[:, 2]] - verts[faces[:, 0]]), axis=1)
    if not np.allclose(areas, geometric_areas, rtol=1e-7, atol=0.0):
        raise ValueError("triangle areas must match the supplied mesh geometry")
    grid_maximum = grid_origin + (np.asarray(potential.shape) - 1) * grid_spacing
    tolerance = 1e-7 * max(float(np.max(grid_spacing)), 1.0)
    if not np.isclose(source_z, grid_maximum[2], rtol=0.0, atol=tolerance):
        raise ValueError("adjoint field gather currently requires the source plane at the grid top")
    if periodic_lateral and not np.allclose(
            bounds, (grid_origin[0], grid_maximum[0], grid_origin[1], grid_maximum[1]),
            rtol=0.0, atol=tolerance):
        raise ValueError("periodic adjoint gather requires full lateral source bounds")
    origin_m = np.asarray(mesh_origin_m, dtype=float)
    if (origin_m.shape != (3,) or np.any(~np.isfinite(origin_m))
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("mesh origin and length unit must define a finite physical coordinate map")
    mapped_reference = origin_m[2] + float(source_z) * float(mesh_length_unit_m)
    if not np.isclose(mapped_reference, boundary.reference_plane_m, rtol=0.0,
                      atol=max(1e-15, 1e-9 * float(mesh_length_unit_m))):
        raise ValueError("mesh source plane does not match PlasmaBoundaryState.reference_plane_m")
    role = dict(species_role); names = {item.name for item in boundary.species}
    allowed_roles = {"energetic_bombardment", "charge_carrier"}
    if set(role) != names or any(value not in allowed_roles for value in role.values()):
        raise ValueError("adjoint field gather currently supports only charged energetic populations")
    if any(item.charge_number == 0 or item.density_model is None for item in boundary.species):
        raise ValueError("adjoint field gather requires charged continuous boundary densities")
    proposals = {} if proposal_by_species is None else dict(proposal_by_species)
    if set(proposals) - names:
        raise ValueError("adjoint proposal names must belong to the physical boundary")
    proposal_frames = (
        {name: proposal_frame_by_species for name in names}
        if isinstance(proposal_frame_by_species, str) else dict(proposal_frame_by_species))
    if (set(proposal_frames) != names
            or any(value not in {"surface_local", "source_aligned"}
                   for value in proposal_frames.values())):
        raise ValueError(
            "proposal_frame_by_species must select surface_local or source_aligned for every species")
    target_faces = (np.arange(len(faces), dtype=int) if gather_face_indices is None
                    else np.asarray(gather_face_indices, dtype=int))
    if (target_faces.ndim != 1 or target_faces.size == 0
            or np.any(target_faces < 0) or np.any(target_faces >= len(faces))
            or np.unique(target_faces).size != target_faces.size):
        raise ValueError("gather_face_indices must select unique triangles from the collision mesh")

    if face_position_seed is not None:
        position_count = int(face_quadrature_points)
        position_level = int(np.ceil(np.log2(position_count)))
        position_u = qmc.Sobol(
            d=2, scramble=True, seed=int(face_position_seed)).random_base2(position_level)
        position_u = position_u[:position_count]
        root = np.sqrt(position_u[:, 0])
        barycentric = np.column_stack((
            1.0 - root, root * (1.0 - position_u[:, 1]), root * position_u[:, 1]))
        point_weight = np.full(
            int(face_quadrature_points), 1.0 / int(face_quadrature_points))
        points = np.einsum("qv,fvc->fqc", barycentric, verts[faces])
    elif int(face_quadrature_points) == 1:
        points = centroids[:, None, :]
        point_weight = np.ones(1)
    elif int(face_quadrature_points) == 3:
        barycentric = np.array([
            [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
        ])
        points = np.einsum("qv,fvc->fqc", barycentric, verts[faces])
        point_weight = np.full(3, 1.0 / 3.0)
    elif int(face_quadrature_points) == 7:
        barycentric = np.array([
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [0.059715871789770, 0.470142064105115, 0.470142064105115],
            [0.470142064105115, 0.059715871789770, 0.470142064105115],
            [0.470142064105115, 0.470142064105115, 0.059715871789770],
            [0.797426985353087, 0.101286507323456, 0.101286507323456],
            [0.101286507323456, 0.797426985353087, 0.101286507323456],
            [0.101286507323456, 0.101286507323456, 0.797426985353087],
        ])
        point_weight = np.array([
            0.225,
            0.132394152788506, 0.132394152788506, 0.132394152788506,
            0.125939180544827, 0.125939180544827, 0.125939180544827,
        ])
        points = np.einsum("qv,fvc->fqc", barycentric, verts[faces])
    reference = np.zeros_like(normals)
    use_z = np.abs(normals[:, 2]) < 0.9
    reference[use_z, 2] = 1.0; reference[~use_z, 0] = 1.0
    tangent_a = np.cross(reference, normals)
    tangent_a /= np.linalg.norm(tangent_a, axis=1)[:, None]
    tangent_b = np.cross(normals, tangent_a)

    selected_device = DEVICE if device is None else str(device)
    if selected_device.startswith("warp:"):
        selected_device = selected_device.split(":", 1)[1]
    ensure_writable_warp_cache(wp)
    mesh = wp.Mesh(
        points=wp.array(verts.astype(np.float32), dtype=wp.vec3, device=selected_device),
        indices=wp.array(faces.astype(np.int32).ravel(), dtype=wp.int32, device=selected_device))
    potential_wp = wp.array(
        np.ascontiguousarray(potential.astype(np.float32)), dtype=float, device=selected_device)
    source_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
    energetic_flux = []; hit_probability = {}; escape_probability = {}; truncation_probability = {}
    for species in boundary.species:
        proposal = proposals.get(species.name, species)
        if (proposal.name != species.name or proposal.charge_number != species.charge_number
                or proposal.density_model is None):
            raise ValueError("adjoint proposals must preserve species name, charge, and density")
        if (proposal_frames[species.name] == "source_aligned"
                and _contains_surface_local_density(proposal.density_model)):
            raise ValueError(
                "folded grazing proposals use surface-local coordinates and cannot be "
                "combined with a source-aligned proposal frame")
        base_velocity = np.asarray(proposal.velocity_sqrt_eV, dtype=float)
        sample_count = base_velocity.shape[0]; point_count = points.shape[1]
        face_index = np.repeat(target_faces, point_count * sample_count)
        point_index = np.tile(np.repeat(np.arange(point_count), sample_count), len(target_faces))
        velocity_index = np.tile(np.arange(sample_count), len(target_faces) * point_count)
        launch_point = points[face_index, point_index] + float(ray_offset) * normals[face_index]
        sampled_velocity = base_velocity[velocity_index]
        if proposal_frames[species.name] == "source_aligned":
            forward_surface_velocity = sampled_velocity.copy()
            forward_surface_velocity[:, 2] *= -1.0
            surface_normal_speed = -np.einsum(
                "rc,rc->r", forward_surface_velocity, normals[face_index])
        else:
            forward_surface_velocity = (
                sampled_velocity[:, 0, None] * tangent_a[face_index]
                + sampled_velocity[:, 1, None] * tangent_b[face_index]
                - sampled_velocity[:, 2, None] * normals[face_index])
            surface_normal_speed = sampled_velocity[:, 2]
        reverse_velocity = -forward_surface_velocity
        ray_count = launch_point.shape[0]
        hit_face_wp = wp.full(ray_count, -1, dtype=wp.int32, device=selected_device)
        hit_cosine_wp = wp.zeros(ray_count, dtype=float, device=selected_device)
        hit_energy_wp = wp.zeros(ray_count, dtype=float, device=selected_device)
        termination_wp = wp.zeros(ray_count, dtype=wp.int8, device=selected_device)
        terminal_position_wp = wp.zeros(ray_count, dtype=wp.vec3, device=selected_device)
        terminal_velocity_wp = wp.zeros(ray_count, dtype=wp.vec3, device=selected_device)
        wp.launch(
            _field_hit_events_3d, dim=ray_count, device=selected_device,
            inputs=[mesh.id, potential_wp, wp.vec3(*grid_origin), wp.vec3(*grid_spacing),
                    wp.vec3(*grid_maximum),
                    wp.array(launch_point.astype(np.float32), dtype=wp.vec3, device=selected_device),
                    wp.array(reverse_velocity.astype(np.float32), dtype=wp.vec3,
                             device=selected_device),
                    float(species.charge_number), float(fixed_dt), int(max_steps),
                    int(bool(periodic_lateral)), hit_face_wp, hit_cosine_wp, hit_energy_wp,
                    termination_wp, terminal_position_wp, terminal_velocity_wp])
        termination = termination_wp.numpy().astype(int)
        terminal_position = terminal_position_wp.numpy().astype(float)
        terminal_velocity = terminal_velocity_wp.numpy().astype(float)
        reached_source = (termination == 2) & (terminal_position[:, 2] > grid_maximum[2])
        if not periodic_lateral:
            reached_source &= (
                (terminal_position[:, 0] >= bounds[0])
                & (terminal_position[:, 0] <= bounds[1])
                & (terminal_position[:, 1] >= bounds[2])
                & (terminal_position[:, 1] <= bounds[3]))
        truncated = termination == 0
        if np.any(truncated):
            raise RuntimeError(
                f"adjoint trajectories for {species.name!r} exhausted max_steps")
        boundary_velocity = np.column_stack((
            -terminal_velocity[:, 0], -terminal_velocity[:, 1], terminal_velocity[:, 2]))
        log_physical = species.log_flux_density(boundary_velocity)
        log_proposal = proposal.log_flux_density(base_velocity)[velocity_index]
        source_normal_speed = np.maximum(terminal_velocity[:, 2], 1e-300)
        with np.errstate(over="ignore", invalid="ignore"):
            density_ratio = np.exp(log_physical - log_proposal)
        contribution = np.where(
            reached_source & (surface_normal_speed > 0.0) & np.isfinite(log_physical),
            proposal.weight[velocity_index] * point_weight[point_index]
            * surface_normal_speed / source_normal_speed * density_ratio,
            0.0)
        positive = contribution > 0.0
        event_face = face_index[positive]
        event_flux = species.flux_m2_s * contribution[positive]
        event_energy = np.sum(forward_surface_velocity[positive] ** 2, axis=1)
        event_cosine = (
            surface_normal_speed[positive]
            / np.linalg.norm(forward_surface_velocity[positive], axis=1))
        energetic_flux.append(FaceResolvedEnergeticFlux(
            species.name, len(faces), event_face, event_flux, event_energy, event_cosine))
        normalized_landing = float(np.dot(
            np.bincount(event_face, weights=contribution[positive], minlength=len(faces)),
            areas) / source_area)
        hit_probability[species.name] = normalized_landing
        escape_probability[species.name] = max(0.0, 1.0 - normalized_landing)
        truncation_probability[species.name] = 0.0
    return BoundaryTransport3DResult(
        SurfaceFluxes({}, tuple(energetic_flux)), hit_probability, escape_probability,
        truncation_probability, "collisionless_fixed_step_nodal_field_adjoint_gather_3d"
        + ("_periodic_cell" if periodic_lateral else ""),
        ("reversible adjoint uses the supplied finite surface velocity quadrature",
         "triangle position quadrature requires refinement at partial visibility",
         "no surface reflection or re-emission",
         "float32 field integration and triangle-ray intersection"))


def trace_boundary_state_bidirectional_field_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str],
        verts, faces, areas, centroids, gas_normals, *, source_bounds, source_z,
        nodal_potential_v, potential_origin, potential_spacing,
        mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        forward_log2_samples=12, adjoint_log2_samples=10, n_replicates=4, seed=0,
        max_forward_log2_samples=None, max_adjoint_log2_samples=None,
        max_face_quadrature_points=None,
        element_absolute_tolerance=0.01, element_relative_tolerance=0.05,
        consistency_sigma=5.0, support_sigma=2.0, support_ratio=0.5,
        proposal_by_species=None, proposal_frame_by_species="surface_local", method_hint=None,
        require_certification=True,
        face_quadrature_points=3, ray_offset=1e-5, fixed_dt=0.01, max_steps=10000,
        periodic_lateral=False, device=None):
    """Audit and select forward or adjoint charged transport independently on every triangle.

    Four or more independent scrambled-QMC replicates supply per-face standard errors. A direct
    forward zero-hit upper bound audits modes missed by every adjoint replicate. If both estimators
    claim precision they must agree within their combined uncertainty. ``method_hint`` freezes a
    previously certified map during a nonlinear-root epoch without silently changing estimators.
    """
    if (int(forward_log2_samples) != forward_log2_samples or forward_log2_samples < 0
            or int(adjoint_log2_samples) != adjoint_log2_samples or adjoint_log2_samples < 0
            or int(n_replicates) != n_replicates or n_replicates < 4
            or not np.isfinite(element_absolute_tolerance) or element_absolute_tolerance < 0.0):
        raise ValueError("invalid bidirectional sampling controls")
    maximum_forward_level = (int(forward_log2_samples) if max_forward_log2_samples is None
                             else int(max_forward_log2_samples))
    maximum_adjoint_level = (int(adjoint_log2_samples) if max_adjoint_log2_samples is None
                             else int(max_adjoint_log2_samples))
    maximum_face_points = (int(face_quadrature_points) if max_face_quadrature_points is None
                           else int(max_face_quadrature_points))
    if (maximum_forward_level < int(forward_log2_samples)
            or maximum_adjoint_level < int(adjoint_log2_samples)
            or maximum_face_points < int(face_quadrature_points)):
        raise ValueError("bidirectional refinement ceilings cannot be below base levels")
    names = {species.name for species in boundary.species}
    roles = dict(species_role)
    proposals = {} if proposal_by_species is None else dict(proposal_by_species)
    frames = ({name: proposal_frame_by_species for name in names}
              if isinstance(proposal_frame_by_species, str)
              else dict(proposal_frame_by_species))
    hints = {} if method_hint is None else dict(method_hint)
    if set(roles) != names or set(proposals) - names or set(frames) != names or set(hints) - names:
        raise ValueError("bidirectional species controls must match the physical boundary")
    areas_array = np.asarray(areas, dtype=float)
    bounds = np.asarray(source_bounds, dtype=float)
    source_area = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2])
    selections = {}; populations = []; hit = {}; escaped = {}; truncated = {}
    for species_index, species in enumerate(boundary.species):
        species_boundary = PlasmaBoundaryState(
            (species,), boundary.reference_plane_m, provenance=boundary.provenance)
        species_roles = {species.name: roles[species.name]}
        forward_results = []; adjoint_results = []; replicate_seeds = []
        template = proposals.get(species.name, species)
        supplied_method = hints.get(species.name)
        frozen_adjoint_faces = None
        frozen_forward_needed = True
        frozen_adjoint_needed = True
        if supplied_method is not None and not require_certification:
            supplied_method = np.asarray(supplied_method).astype("U7")
            if (supplied_method.shape != (len(faces),)
                    or np.any(~np.isin(supplied_method, ("forward", "adjoint")))):
                raise ValueError("method_hint must select forward or adjoint for every face")
            selected = np.where(supplied_method == "adjoint")[0]
            # Adjoint work is face-local. Once a separately certified method map is frozen, do not
            # retrace faces whose events will be discarded in favor of the global forward estimator.
            if selected.size:
                frozen_adjoint_faces = selected
            else:
                frozen_adjoint_needed = False
            if np.all(supplied_method == "adjoint"):
                frozen_forward_needed = False
        for replicate in range(int(n_replicates)):
            replicate_seed = int(seed) + 104729 * replicate + 15485863 * species_index
            replicate_seeds.append(replicate_seed)
            if frozen_forward_needed:
                forward_results.append(trace_boundary_state_field_3d(
                    species_boundary, species_roles, verts, faces, areas,
                    source_bounds=source_bounds, source_z=source_z,
                    nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
                    potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
                    mesh_origin_m=mesh_origin_m, seed=replicate_seed, fixed_dt=fixed_dt,
                    max_steps=max_steps, phase_space_log2_samples=int(forward_log2_samples),
                    periodic_lateral=periodic_lateral, device=device))
            if frozen_adjoint_needed:
                proposal = qmc_boundary_proposal(
                    template, int(adjoint_log2_samples), replicate_seed,
                    name=species.name)
                adjoint_results.append(gather_boundary_state_field_adjoint_3d(
                    species_boundary, species_roles, verts, faces, areas, centroids, gas_normals,
                    source_bounds=source_bounds, source_z=source_z,
                    nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
                    potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
                    mesh_origin_m=mesh_origin_m, face_quadrature_points=face_quadrature_points,
                    ray_offset=ray_offset, fixed_dt=fixed_dt, max_steps=max_steps,
                    periodic_lateral=periodic_lateral,
                    proposal_by_species={species.name: proposal},
                    proposal_frame_by_species={species.name: frames[species.name]},
                    face_position_seed=replicate_seed,
                    gather_face_indices=frozen_adjoint_faces, device=device))
        empty_population = FaceResolvedEnergeticFlux(
            species.name, len(faces), np.empty(0, dtype=int), np.empty(0),
            np.empty(0), np.empty(0))
        forward_populations = (
            [item.surface_fluxes.energetic_fluxes[0] for item in forward_results]
            if frozen_forward_needed else [empty_population] * int(n_replicates))
        adjoint_populations = (
            [item.surface_fluxes.energetic_fluxes[0] for item in adjoint_results]
            if frozen_adjoint_needed else [empty_population] * int(n_replicates))
        pooled_forward_samples = int(n_replicates) * 2 ** int(forward_log2_samples)
        zero_upper = (3.0 * species.flux_m2_s * source_area
                      / (pooled_forward_samples * areas_array))
        selection = select_bidirectional_face_events_3d(
            forward_populations, adjoint_populations,
            forward_zero_upper_m2_s=zero_upper,
            absolute_tolerance_m2_s=float(element_absolute_tolerance) * species.flux_m2_s,
            relative_tolerance=element_relative_tolerance,
            consistency_sigma=consistency_sigma, support_sigma=support_sigma,
            support_ratio=support_ratio, method_hint=hints.get(species.name),
            require_certification=require_certification)
        forward_level = int(forward_log2_samples)
        adjoint_level = int(adjoint_log2_samples)
        position_count = int(face_quadrature_points)
        while not selection.converged:
            unresolved = np.where(
                ~(selection.estimator_consistent & selection.method_within_tolerance))[0]
            refine_forward = forward_level < maximum_forward_level
            refine_adjoint = (
                adjoint_level < maximum_adjoint_level
                or position_count < maximum_face_points)
            if unresolved.size == 0 or not (refine_forward or refine_adjoint):
                break
            if refine_forward:
                forward_level += 1
                forward_results = [trace_boundary_state_field_3d(
                    species_boundary, species_roles, verts, faces, areas,
                    source_bounds=source_bounds, source_z=source_z,
                    nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
                    potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
                    mesh_origin_m=mesh_origin_m, seed=replicate_seed, fixed_dt=fixed_dt,
                    max_steps=max_steps, phase_space_log2_samples=forward_level,
                    periodic_lateral=periodic_lateral, device=device)
                    for replicate_seed in replicate_seeds]
                forward_populations = [item.surface_fluxes.energetic_fluxes[0]
                                       for item in forward_results]
            if refine_adjoint:
                adjoint_level = min(adjoint_level + 1, maximum_adjoint_level)
                position_count = min(2 * position_count, maximum_face_points)
                for replicate, replicate_seed in enumerate(replicate_seeds):
                    proposal = qmc_boundary_proposal(
                        template, adjoint_level, replicate_seed, name=species.name)
                    partial = gather_boundary_state_field_adjoint_3d(
                        species_boundary, species_roles, verts, faces, areas,
                        centroids, gas_normals, source_bounds=source_bounds, source_z=source_z,
                        nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
                        potential_spacing=potential_spacing,
                        mesh_length_unit_m=mesh_length_unit_m, mesh_origin_m=mesh_origin_m,
                        face_quadrature_points=position_count, ray_offset=ray_offset,
                        fixed_dt=fixed_dt, max_steps=max_steps, periodic_lateral=periodic_lateral,
                        proposal_by_species={species.name: proposal},
                        proposal_frame_by_species={species.name: frames[species.name]},
                        face_position_seed=replicate_seed, gather_face_indices=unresolved,
                        device=device).surface_fluxes.energetic_fluxes[0]
                    adjoint_populations[replicate] = _replace_population_face_events_3d(
                        adjoint_populations[replicate], partial, unresolved)
            pooled_forward_samples = int(n_replicates) * 2 ** forward_level
            zero_upper = (3.0 * species.flux_m2_s * source_area
                          / (pooled_forward_samples * areas_array))
            selection = select_bidirectional_face_events_3d(
                forward_populations, adjoint_populations,
                forward_zero_upper_m2_s=zero_upper,
                absolute_tolerance_m2_s=float(element_absolute_tolerance) * species.flux_m2_s,
                relative_tolerance=element_relative_tolerance,
                consistency_sigma=consistency_sigma, support_sigma=support_sigma,
                support_ratio=support_ratio, method_hint=hints.get(species.name),
                require_certification=require_certification)
        selections[species.name] = selection; populations.append(selection.population)
        landing = float(np.dot(selection.population.flux_m2_s, areas_array)
                        / (species.flux_m2_s * source_area)) if species.flux_m2_s > 0.0 else 0.0
        hit[species.name] = landing
        escaped[species.name] = max(0.0, 1.0 - landing)
        truncation = [
            item.truncation_probability[species.name]
            for item in forward_results + adjoint_results]
        truncated[species.name] = max(truncation) if truncation else 0.0
    transport = BoundaryTransport3DResult(
        SurfaceFluxes({}, tuple(populations)), hit, escaped, truncated,
        "collisionless_fixed_step_nodal_field_bidirectional_3d"
        + ("_periodic_cell" if periodic_lateral else ""),
        ("per-face estimator map must be frozen during a nonlinear-root epoch",
         "scrambled-QMC uncertainty requires independent replicate and level refinement",
         "no surface reflection or re-emission"))
    return BidirectionalBoundaryTransport3DResult(transport, selections)
