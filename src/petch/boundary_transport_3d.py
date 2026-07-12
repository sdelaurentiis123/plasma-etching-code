"""Species-resolved 3-D transport from ``PlasmaBoundaryState`` to an arbitrary triangle mesh.

This module is the dimensional bridge between the common reactor/sheath boundary contract and surface
kinetics.  The first backend is deliberately limited to collisionless, absorbing, first-hit transport.  It
preserves the exact discrete boundary energy-angle measure at every hit and reports its limitations; it is
not a replacement for the charged, reflecting, state-coupled production transport still under construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.stats import qmc
import warp as wp

from .boundary_state import PlasmaBoundaryState
from .surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes
from .threed import DEVICE
from .warp_runtime import ensure_writable_warp_cache


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


@dataclass(frozen=True)
class BoundaryTransport3DResult:
    surface_fluxes: SurfaceFluxes
    hit_probability: Mapping[str, float]
    escape_probability: Mapping[str, float]
    transport_model: str
    known_limitations: tuple[str, ...]

    def __post_init__(self):
        object.__setattr__(self, "hit_probability", MappingProxyType(dict(self.hit_probability)))
        object.__setattr__(self, "escape_probability", MappingProxyType(dict(self.escape_probability)))
        object.__setattr__(self, "known_limitations", tuple(self.known_limitations))


def trace_boundary_state_first_hit_3d(
        boundary: PlasmaBoundaryState, species_role: Mapping[str, str], verts, faces, areas, *,
        source_bounds, source_z, mesh_length_unit_m=1e-6, mesh_origin_m=(0.0, 0.0, 0.0),
        n_position=256, seed=0, max_distance=None, device=None):
    """Transport a spatially uniform boundary state to exact triangle-hit events.

    ``species_role`` is a physical input mapping each boundary species to ``"neutral_reactant"`` or
    ``"energetic_bombardment"``; no species name selects a formula. Mesh coordinates may use any
    length unit declared by ``mesh_length_unit_m``. The mapped source plane must equal the boundary's
    SI reference plane, preventing a silent geometry/boundary unit mismatch.

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
    allowed_roles = {"neutral_reactant", "energetic_bombardment"}
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
        wp.launch(
            _first_hit_events_3d, dim=ray_count, device=selected_device,
            inputs=[mesh.id,
                    wp.array(origin.astype(np.float32), dtype=wp.vec3, device=selected_device),
                    wp.array(direction.astype(np.float32), dtype=wp.vec3, device=selected_device),
                    float(max_distance), hit_face_wp, hit_cosine_wp])
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
        transport_model="collisionless_absorbing_first_hit_3d",
        known_limitations=(
            "no intra-feature electric-field trajectory coupling",
            "no surface reflection or neutral re-emission",
            "no spatially varying boundary density",
            "float32 triangle-ray intersection",
        ))
