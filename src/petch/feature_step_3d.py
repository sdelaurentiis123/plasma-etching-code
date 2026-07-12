"""Dimensional feature evolution through the new physical contracts.

Each step transfers surface state material-by-material with an area-conservative, bounded remap. Smooth
CFL-limited motion can be iterated; topology changes, material appearance/disappearance, excessive remap
distance, and impossible coverage compression are refused. This makes the multi-step loop explicit about
the domain in which surface history is numerically supported.
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
    next_surface_state: SiO2SurfaceState
    next_active_face_centroid: np.ndarray
    next_active_face_area: np.ndarray
    next_surface_state_mesh_fingerprint: str
    state_remap_diagnostics: Mapping[str, object]
    face_material_id: np.ndarray
    face_velocity_mesh_units_s: np.ndarray
    diagnostics: Mapping[str, object]
    validity: FeatureStepValidity

    def __post_init__(self):
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))
        object.__setattr__(
            self, "state_remap_diagnostics", MappingProxyType(dict(self.state_remap_diagnostics)))


@dataclass(frozen=True)
class FeatureSolve3DResult:
    geometry: FeatureGeometry3D
    surface_state: SiO2SurfaceState
    surface_state_mesh_fingerprint: str
    steps: tuple[FeatureStep3DResult, ...]
    duration_s: float
    validity: FeatureStepValidity


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


def _surface_topology_signature(faces, active_face):
    active = np.asarray(faces, dtype=int)[np.asarray(active_face, dtype=int)]
    if active.size == 0:
        return 0, 0
    edge = np.concatenate((active[:, [0, 1]], active[:, [1, 2]], active[:, [2, 0]]))
    edge.sort(axis=1)
    edge_count = np.unique(edge, axis=0).shape[0]
    vertex_count = np.unique(active).size
    parent = np.arange(active.shape[0])

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left = find(left); right = find(right)
        if left != right:
            parent[right] = left

    owner = {}
    for face_index, vertices in enumerate(active):
        for vertex in vertices:
            vertex = int(vertex)
            if vertex in owner:
                union(face_index, owner[vertex])
            else:
                owner[vertex] = face_index
    components = len({find(index) for index in range(active.shape[0])})
    euler_characteristic = int(vertex_count - edge_count + active.shape[0])
    return int(components), euler_characteristic


def _conserve_nonnegative_surface_field(raw, target_integral, new_area, *, upper=None):
    raw = np.maximum(np.asarray(raw, dtype=float), 0.0)
    area = np.asarray(new_area, dtype=float)
    target = float(target_integral)
    scale = max(abs(target), 1.0)
    if target < -1e-13 * scale:
        raise ValueError("negative conservative-remap target")
    if target <= 1e-15 * scale:
        return np.zeros_like(raw)
    if upper is None:
        raw_integral = float(np.dot(raw, area))
        if raw_integral <= 0.0:
            raw = np.ones_like(raw)
            raw_integral = float(area.sum())
        return raw * (target / raw_integral)
    capacity = float(upper) * float(area.sum())
    if target > capacity * (1.0 + 5e-13):
        raise ValueError("surface contraction exceeds bounded coverage capacity")
    seed = raw if np.any(raw > 0.0) else np.ones_like(raw)

    def integral(multiplier):
        return float(np.dot(np.minimum(multiplier * seed, upper), area))

    lower = 0.0; upper_multiplier = 1.0
    while integral(upper_multiplier) < target:
        upper_multiplier *= 2.0
        if upper_multiplier > 1e300:
            raise RuntimeError("bounded conservative remap failed to bracket target")
    for _ in range(80):
        midpoint = 0.5 * (lower + upper_multiplier)
        if integral(midpoint) < target:
            lower = midpoint
        else:
            upper_multiplier = midpoint
    return np.minimum(upper_multiplier * seed, upper)


def conservative_remap_surface_state(
        state, old_centroid, old_area, old_material, new_centroid, new_area, new_material, *,
        dx, mesh_length_unit_m, neighbor_count=4, maximum_distance=None):
    """First-order material-local remap with exact area-integrated state conservation.

    The interpolation supplies spatial locality; a constrained correction then preserves each material's
    integrated complex sites, polymer units, and cumulative removed formula units. Complex coverage remains
    in [0,1]. This operator does not authorize topology change; the caller must gate topology separately.
    """
    old_centroid = np.asarray(old_centroid, dtype=float)
    new_centroid = np.asarray(new_centroid, dtype=float)
    old_area = np.asarray(old_area, dtype=float); new_area = np.asarray(new_area, dtype=float)
    old_material = np.asarray(old_material, dtype=int)
    new_material = np.asarray(new_material, dtype=int)
    if (old_centroid.ndim != 2 or old_centroid.shape[1] != 3
            or new_centroid.ndim != 2 or new_centroid.shape[1] != 3
            or old_area.shape != (old_centroid.shape[0],)
            or new_area.shape != (new_centroid.shape[0],)
            or old_material.shape != old_area.shape or new_material.shape != new_area.shape
            or state.complex_fraction.shape != old_area.shape
            or np.any(old_area <= 0.0) or np.any(new_area <= 0.0)):
        raise ValueError("invalid surface-state remap geometry")
    if maximum_distance is None:
        maximum_distance = 2.0 * float(dx)
    if not np.isfinite(maximum_distance) or maximum_distance <= 0.0:
        raise ValueError("maximum remap distance must be positive")
    output = [np.zeros(new_area.shape) for _ in range(3)]
    maximum_nearest = 0.0; material_diagnostics = {}
    physical_area_scale = float(mesh_length_unit_m) ** 2
    for material in sorted(set(old_material) | set(new_material)):
        old_index = np.where(old_material == material)[0]
        new_index = np.where(new_material == material)[0]
        if old_index.size == 0 or new_index.size == 0:
            raise ValueError(
                "material surface appeared or disappeared; initialize/retire state explicitly")
        count = min(int(neighbor_count), old_index.size)
        distance, local = cKDTree(old_centroid[old_index]).query(
            new_centroid[new_index], k=count)
        if count == 1:
            distance = np.asarray(distance)[:, None]; local = np.asarray(local)[:, None]
        nearest = float(np.max(distance[:, 0])); maximum_nearest = max(maximum_nearest, nearest)
        if nearest > maximum_distance:
            raise ValueError(
                f"surface remap distance {nearest:g} exceeds {maximum_distance:g}")
        source_index = old_index[np.asarray(local, dtype=int)]
        regularization = (0.25 * float(dx)) ** 2
        weight = old_area[source_index] / (distance * distance + regularization)
        weight /= weight.sum(axis=1, keepdims=True)
        old_values = (
            state.complex_fraction, state.polymer_units_m2,
            state.removed_formula_units_m2)
        targets = []; residuals = []
        for field_index, old_value in enumerate(old_values):
            raw = np.sum(weight * old_value[source_index], axis=1)
            target = float(np.dot(old_value[old_index], old_area[old_index]))
            remapped = _conserve_nonnegative_surface_field(
                raw, target, new_area[new_index], upper=1.0 if field_index == 0 else None)
            output[field_index][new_index] = remapped
            achieved = float(np.dot(remapped, new_area[new_index]))
            residuals.append(abs(achieved - target) / max(abs(target), 1.0))
            targets.append(target * physical_area_scale)
        material_diagnostics[int(material)] = dict(
            old_face_count=int(old_index.size), new_face_count=int(new_index.size),
            old_area_m2=float(old_area[old_index].sum() * physical_area_scale),
            new_area_m2=float(new_area[new_index].sum() * physical_area_scale),
            target_complex_area_m2=float(targets[0]),
            target_polymer_units=float(targets[1]),
            target_removed_formula_units=float(targets[2]),
            max_relative_conservation_residual=float(max(residuals)))
    remapped_state = SiO2SurfaceState(*output)
    return remapped_state, dict(
        method="material_local_area_conservative_knn",
        neighbor_count=int(neighbor_count), maximum_nearest_distance=float(maximum_nearest),
        maximum_allowed_distance=float(maximum_distance),
        materials=material_diagnostics)


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
    next_verts, next_faces, next_centroids, next_areas = extract_mesh_3d(
        output_geometry.phi, output_geometry.dx)
    next_face_material = _face_material_ids(next_centroids, output_geometry)
    next_active_face = np.where(np.isin(next_face_material, etchable))[0]
    if next_active_face.size == 0:
        raise ValueError("etch step removed every requested material surface")
    old_topology = _surface_topology_signature(faces, active_face)
    next_topology = _surface_topology_signature(next_faces, next_active_face)
    if old_topology != next_topology:
        raise ValueError(
            f"surface topology changed from {old_topology} to {next_topology}; "
            "state transfer requires an explicit topology event")
    next_surface_state, remap_diagnostics = conservative_remap_surface_state(
        surface.state, centroids[active_face], areas[active_face], face_material[active_face],
        next_centroids[next_active_face], next_areas[next_active_face],
        next_face_material[next_active_face], dx=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        maximum_distance=displacement + 1.5 * geometry.dx)
    next_mesh_fingerprint = _surface_mesh_fingerprint(
        next_verts, next_faces, next_active_face, next_face_material, output_geometry)
    remap_diagnostics = dict(
        remap_diagnostics, old_topology=old_topology, new_topology=next_topology,
        next_active_face_count=int(next_active_face.size))
    reasons = []
    if not surface.validity.within_declared_scope:
        reasons.extend(surface.validity.reasons)
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(transport.known_limitations) + (
            "first-order material-local conservative surface-state remap",
            "topology-changing surface steps are refused",
            "first-order Godunov interface advection",
        ) + tuple(surface.validity.known_model_form_omissions))
    return FeatureStep3DResult(
        geometry=output_geometry, transport=transport, surface=surface,
        active_face_index=active_face, active_face_centroid=centroids[active_face],
        active_face_area=areas[active_face],
        surface_state_mesh_fingerprint=mesh_fingerprint,
        next_surface_state=next_surface_state,
        next_active_face_centroid=next_centroids[next_active_face],
        next_active_face_area=next_areas[next_active_face],
        next_surface_state_mesh_fingerprint=next_mesh_fingerprint,
        state_remap_diagnostics=remap_diagnostics,
        face_material_id=face_material,
        face_velocity_mesh_units_s=face_velocity,
        diagnostics=dict(
            face_count=int(len(faces)), active_face_count=int(active_face.size),
            max_velocity_m_s=maximum_velocity * geometry.mesh_length_unit_m,
            max_displacement_mesh_units=displacement, cfl_substeps=int(substeps),
            cfl_number=float(cfl_number), reinitialized=bool(reinitialize)),
        validity=validity)


def solve_feature_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism: ReducedSiO2FluorocarbonMechanism, *,
        etchable_material_ids, duration_s, n_steps, source_bounds, source_z,
        n_position=256, seed=0, cfl_number=0.3, reinitialize=True,
        transport_device=None):
    """Run multiple verified feature steps, carrying only conservatively remapped surface state."""
    if int(n_steps) != n_steps or int(n_steps) <= 0:
        raise ValueError("n_steps must be a positive integer")
    if not np.isfinite(duration_s) or duration_s < 0.0:
        raise ValueError("duration_s must be finite and nonnegative")
    step_duration = float(duration_s) / int(n_steps)
    current_geometry = geometry; current_state = None; current_fingerprint = None
    results = []
    for step_index in range(int(n_steps)):
        result = advance_feature_step_3d(
            current_geometry, boundary, species_role, mechanism,
            etchable_material_ids=etchable_material_ids, duration_s=step_duration,
            source_bounds=source_bounds, source_z=source_z,
            surface_state=current_state,
            surface_state_mesh_fingerprint=current_fingerprint,
            n_position=n_position, seed=int(seed) + step_index,
            cfl_number=cfl_number, reinitialize=reinitialize,
            transport_device=transport_device)
        results.append(result)
        current_geometry = result.geometry
        current_state = result.next_surface_state
        current_fingerprint = result.next_surface_state_mesh_fingerprint
    reasons = tuple(reason for result in results for reason in result.validity.reasons)
    limitations = tuple(dict.fromkeys(
        limitation for result in results for limitation in result.validity.known_limitations))
    return FeatureSolve3DResult(
        geometry=current_geometry, surface_state=current_state,
        surface_state_mesh_fingerprint=current_fingerprint,
        steps=tuple(results), duration_s=float(duration_s),
        validity=FeatureStepValidity(not reasons, reasons, limitations))
