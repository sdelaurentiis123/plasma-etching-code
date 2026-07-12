"""One fully dimensional feature-evolution step through the new physical contracts.

The step is intentionally not wrapped in a multi-step loop yet.  Its returned surface state is attached to
the pre-step mesh; iterating it requires a separately verified conservative state-remap operator.  Refusing
to hide nearest-face history copying is part of the solver's validity contract.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.spatial import cKDTree

from .boundary_state import PlasmaBoundaryState
from .boundary_transport_3d import BoundaryTransport3DResult, trace_boundary_state_first_hit_3d
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    ReducedSiO2FluorocarbonMechanism,
    SiO2SurfaceState,
    SurfaceFluxes,
    SurfaceStepResult,
)
from .threed import advect_3d, extend_velocity_3d, extract_mesh_3d, reinit_narrow


@dataclass(frozen=True)
class FeatureGeometry3D:
    """Eulerian material geometry in declared mesh units; material id zero is gas."""

    phi: np.ndarray
    material_id: np.ndarray
    dx: float
    mesh_length_unit_m: float
    mesh_origin_m: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self):
        phi = np.asarray(self.phi, dtype=float).copy()
        material = np.asarray(self.material_id, dtype=int).copy()
        origin = tuple(float(value) for value in self.mesh_origin_m)
        if (phi.ndim != 3 or min(phi.shape) < 2 or material.shape != phi.shape
                or np.any(~np.isfinite(phi)) or np.any(material < 0)
                or not np.isfinite(self.dx) or self.dx <= 0.0
                or not np.isfinite(self.mesh_length_unit_m) or self.mesh_length_unit_m <= 0.0
                or len(origin) != 3 or np.any(~np.isfinite(origin))):
            raise ValueError("invalid 3-D feature geometry")
        if not np.any(phi < 0.0) or not np.any(phi > 0.0):
            raise ValueError("phi must contain both gas and solid")
        phi.setflags(write=False); material.setflags(write=False)
        object.__setattr__(self, "phi", phi)
        object.__setattr__(self, "material_id", material)
        object.__setattr__(self, "dx", float(self.dx))
        object.__setattr__(self, "mesh_length_unit_m", float(self.mesh_length_unit_m))
        object.__setattr__(self, "mesh_origin_m", origin)

    @property
    def coordinate_arrays(self):
        return tuple(np.arange(size) * self.dx for size in self.phi.shape)


@dataclass(frozen=True)
class FeatureStepValidity:
    within_declared_scope: bool
    reasons: tuple[str, ...]
    known_limitations: tuple[str, ...]


@dataclass(frozen=True)
class FeatureStep3DResult:
    geometry: FeatureGeometry3D
    transport: BoundaryTransport3DResult
    surface: SurfaceStepResult
    active_face_index: np.ndarray
    active_face_centroid: np.ndarray
    active_face_area: np.ndarray
    surface_state_mesh_fingerprint: str
    face_material_id: np.ndarray
    face_velocity_mesh_units_s: np.ndarray
    diagnostics: Mapping[str, object]
    validity: FeatureStepValidity

    def __post_init__(self):
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


def _face_material_ids(centroids, geometry):
    """Assign the nearest positive-phi material node to each interface triangle."""
    solid = (geometry.phi > 0.0) & (geometry.material_id > 0)
    index = np.column_stack(np.where(solid))
    if index.size == 0:
        raise ValueError("geometry contains no labeled solid material")
    points = index * geometry.dx
    _, nearest = cKDTree(points).query(centroids)
    chosen = index[np.asarray(nearest, dtype=int)]
    return geometry.material_id[tuple(chosen.T)]


def _surface_mesh_fingerprint(verts, faces, active_face, face_material, geometry):
    digest = sha256()
    for array, dtype in (
            (verts, "<f8"), (faces, "<i8"), (active_face, "<i8"),
            (face_material, "<i8")):
        digest.update(np.ascontiguousarray(array, dtype=dtype).tobytes())
    digest.update(np.asarray(
        [geometry.dx, geometry.mesh_length_unit_m, *geometry.mesh_origin_m],
        dtype="<f8").tobytes())
    return digest.hexdigest()


def _select_surface_fluxes(fluxes, selected_face, face_count):
    selected_face = np.asarray(selected_face, dtype=int)
    old_to_new = np.full(int(face_count), -1, dtype=int)
    old_to_new[selected_face] = np.arange(selected_face.size)
    neutral = {
        name: np.asarray(value)[selected_face]
        for name, value in fluxes.neutral_flux_m2_s.items()}
    energetic = []
    for population in fluxes.energetic_fluxes:
        if isinstance(population, FaceResolvedEnergeticFlux):
            mapped = old_to_new[population.event_face]
            retained = mapped >= 0
            energetic.append(FaceResolvedEnergeticFlux(
                population.name, selected_face.size, mapped[retained],
                population.event_flux_m2_s[retained], population.event_energy_eV[retained],
                population.event_cosine_incidence[retained]))
        elif isinstance(population, EnergeticFlux):
            flux = np.asarray(population.flux_m2_s)
            selected_flux = flux if flux.ndim == 0 else flux[selected_face]
            energetic.append(EnergeticFlux(
                population.name, selected_flux, population.energy_eV,
                population.cosine_incidence, population.weight))
        else:  # pragma: no cover - SurfaceFluxes already validates this
            raise TypeError(type(population).__name__)
    return SurfaceFluxes(neutral, tuple(energetic))


def advance_feature_step_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism: ReducedSiO2FluorocarbonMechanism, *,
        etchable_material_ids, duration_s, source_bounds, source_z,
        surface_state: SiO2SurfaceState | None = None, n_position=256, seed=0,
        surface_state_mesh_fingerprint=None,
        cfl_number=0.3, reinitialize=True, transport_device=None):
    """Advance one stateful, dimensional, collisionless-absorbing feature step.

    The chemistry is evaluated only on triangles whose nearest positive-phi material id is in
    ``etchable_material_ids``. Other labeled solids are pinned. The method refuses a supplied surface
    state whose shape does not match the current active mesh; it never silently remaps history.
    """
    if not np.isfinite(duration_s) or duration_s < 0.0:
        raise ValueError("duration_s must be finite and nonnegative")
    if not np.isfinite(cfl_number) or not 0.0 < cfl_number < 1.0:
        raise ValueError("cfl_number must lie strictly between zero and one")
    etchable = tuple(sorted({int(value) for value in etchable_material_ids}))
    if not etchable or any(value <= 0 for value in etchable):
        raise ValueError("etchable material ids must be positive")

    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    face_material = _face_material_ids(centroids, geometry)
    active_face = np.where(np.isin(face_material, etchable))[0]
    if active_face.size == 0:
        raise ValueError("current interface contains no requested etchable material")
    mesh_fingerprint = _surface_mesh_fingerprint(
        verts, faces, active_face, face_material, geometry)
    transport = trace_boundary_state_first_hit_3d(
        boundary, species_role, verts, faces, areas,
        source_bounds=source_bounds, source_z=source_z,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
        device=transport_device)
    active_flux = _select_surface_fluxes(
        transport.surface_fluxes, active_face, len(faces))
    if surface_state is None:
        if surface_state_mesh_fingerprint is not None:
            raise ValueError("surface_state_mesh_fingerprint requires a supplied surface_state")
        surface_state = SiO2SurfaceState.bare((active_face.size,))
    else:
        if surface_state_mesh_fingerprint != mesh_fingerprint:
            raise ValueError(
                "surface_state mesh fingerprint mismatch; conservative remap is required")
        if surface_state.complex_fraction.shape != (active_face.size,):
            raise ValueError(
                "surface_state does not match the current active mesh; conservative remap is required")
    surface = mechanism.advance(surface_state, active_flux, float(duration_s))

    face_velocity = np.zeros(len(faces))
    face_velocity[active_face] = (
        surface.etch_velocity_m_s / geometry.mesh_length_unit_m)
    maximum_velocity = float(np.max(face_velocity)) if face_velocity.size else 0.0
    displacement = maximum_velocity * float(duration_s)
    substeps = max(1, int(np.ceil(displacement / (float(cfl_number) * geometry.dx))))
    phi = np.array(geometry.phi, copy=True)
    xs, ys, zs = geometry.coordinate_arrays
    extension_geometry = dict(phi=phi, dx=geometry.dx, xs=xs, ys=ys, zs=zs)
    extended_velocity = extend_velocity_3d(
        face_velocity, centroids, extension_geometry, 4.0 * geometry.dx)
    pinned = (geometry.material_id > 0) & ~np.isin(geometry.material_id, etchable)
    for _ in range(substeps):
        phi = advect_3d(
            phi, extended_velocity, geometry.dx, float(duration_s) / substeps)
        phi[pinned] = geometry.phi[pinned]
    if reinitialize and duration_s > 0.0:
        phi = reinit_narrow(phi, geometry.dx, 4.0 * geometry.dx)
        phi[pinned] = geometry.phi[pinned]

    output_geometry = FeatureGeometry3D(
        phi, geometry.material_id, geometry.dx, geometry.mesh_length_unit_m,
        geometry.mesh_origin_m)
    reasons = []
    if not surface.validity.within_declared_scope:
        reasons.extend(surface.validity.reasons)
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(transport.known_limitations) + (
            "surface state is attached to the pre-step mesh; no conservative remap yet",
            "first-order Godunov interface advection",
        ) + tuple(surface.validity.known_model_form_omissions))
    return FeatureStep3DResult(
        geometry=output_geometry, transport=transport, surface=surface,
        active_face_index=active_face, active_face_centroid=centroids[active_face],
        active_face_area=areas[active_face],
        surface_state_mesh_fingerprint=mesh_fingerprint,
        face_material_id=face_material,
        face_velocity_mesh_units_s=face_velocity,
        diagnostics=dict(
            face_count=int(len(faces)), active_face_count=int(active_face.size),
            max_velocity_m_s=maximum_velocity * geometry.mesh_length_unit_m,
            max_displacement_mesh_units=displacement, cfl_substeps=int(substeps),
            cfl_number=float(cfl_number), reinitialized=bool(reinitialize)),
        validity=validity)
