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
from .boundary_transport_3d import (
    BoundaryTransport3DResult,
    merge_boundary_transport_results_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from .charging_coupled_3d import (
    SteadyDielectricCharging3DResult, solve_dielectric_charging_steady_3d,
)
from .charging_poisson_3d import NodalPoissonSystem3D
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
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
    charging: SteadyDielectricCharging3DResult | None
    surface: object
    active_face_index: np.ndarray
    active_face_centroid: np.ndarray
    active_face_area: np.ndarray
    surface_state_mesh_fingerprint: str
    next_surface_state: object
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
    surface_state: object
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

    The state declares named nonnegative fields, optional upper bounds, and reconstruction. Interpolation
    supplies spatial locality; a constrained correction then preserves every material/field area integral.
    This operator does not authorize topology change; the caller must gate topology separately.
    """
    old_centroid = np.asarray(old_centroid, dtype=float)
    new_centroid = np.asarray(new_centroid, dtype=float)
    old_area = np.asarray(old_area, dtype=float); new_area = np.asarray(new_area, dtype=float)
    old_material = np.asarray(old_material, dtype=int)
    new_material = np.asarray(new_material, dtype=int)
    if (not hasattr(state, "conservative_surface_fields")
            or not hasattr(state, "conservative_surface_upper_bounds")
            or not hasattr(state, "with_conservative_surface_fields")):
        raise TypeError("surface state does not implement the conservative remap contract")
    old_values = dict(state.conservative_surface_fields())
    upper_bounds = dict(state.conservative_surface_upper_bounds())
    if not old_values or set(upper_bounds) != set(old_values):
        raise ValueError("surface-state remap fields and upper bounds must match")
    old_values = {name: np.asarray(value, dtype=float) for name, value in old_values.items()}
    if (old_centroid.ndim != 2 or old_centroid.shape[1] != 3
            or new_centroid.ndim != 2 or new_centroid.shape[1] != 3
            or old_area.shape != (old_centroid.shape[0],)
            or new_area.shape != (new_centroid.shape[0],)
            or old_material.shape != old_area.shape or new_material.shape != new_area.shape
            or any(value.shape != old_area.shape for value in old_values.values())
            or np.any(old_area <= 0.0) or np.any(new_area <= 0.0)):
        raise ValueError("invalid surface-state remap geometry")
    if maximum_distance is None:
        maximum_distance = 2.0 * float(dx)
    if not np.isfinite(maximum_distance) or maximum_distance <= 0.0:
        raise ValueError("maximum remap distance must be positive")
    output = {name: np.zeros(new_area.shape) for name in old_values}
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
        targets = {}; residuals = []
        for field_name, old_value in old_values.items():
            raw = np.sum(weight * old_value[source_index], axis=1)
            target = float(np.dot(old_value[old_index], old_area[old_index]))
            remapped = _conserve_nonnegative_surface_field(
                raw, target, new_area[new_index], upper=upper_bounds[field_name])
            output[field_name][new_index] = remapped
            achieved = float(np.dot(remapped, new_area[new_index]))
            residuals.append(abs(achieved - target) / max(abs(target), 1.0))
            targets[field_name] = target * physical_area_scale
        material_diagnostics[int(material)] = dict(
            old_face_count=int(old_index.size), new_face_count=int(new_index.size),
            old_area_m2=float(old_area[old_index].sum() * physical_area_scale),
            new_area_m2=float(new_area[new_index].sum() * physical_area_scale),
            target_field_integrals=targets,
            max_relative_conservation_residual=float(max(residuals)))
    remapped_state = state.with_conservative_surface_fields(output)
    return remapped_state, dict(
        method="material_local_area_conservative_knn",
        neighbor_count=int(neighbor_count), maximum_nearest_distance=float(maximum_nearest),
        maximum_allowed_distance=float(maximum_distance),
        materials=material_diagnostics)


def _select_surface_fluxes(fluxes, selected_face, face_count, species_role=None):
    selected_face = np.asarray(selected_face, dtype=int)
    role = None if species_role is None else dict(species_role)
    old_to_new = np.full(int(face_count), -1, dtype=int)
    old_to_new[selected_face] = np.arange(selected_face.size)
    neutral = {
        name: np.asarray(value)[selected_face]
        for name, value in fluxes.neutral_flux_m2_s.items()
        if role is None or role.get(name) == "neutral_reactant"}
    energetic = []
    for population in fluxes.energetic_fluxes:
        if role is not None and role.get(population.name) != "energetic_bombardment":
            continue
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
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, source_bounds, source_z,
        surface_state=None, n_position=256, seed=0,
        surface_state_mesh_fingerprint=None,
        nodal_potential_v=None, potential_origin=None, potential_spacing=None,
        trajectory_fixed_dt=None, trajectory_max_steps=10000,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        initial_charge_node_c=None, charging_options=None,
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
    role = dict(species_role)
    if set(role) != {species.name for species in boundary.species}:
        raise ValueError("species_role must classify every and only boundary species")
    allowed_roles = {"neutral_reactant", "energetic_bombardment", "charge_carrier"}
    if any(value not in allowed_roles for value in role.values()):
        raise ValueError(f"species roles must be one of {sorted(allowed_roles)}")
    if any(species.charge_number != 0 and role[species.name] == "neutral_reactant"
           for species in boundary.species):
        raise ValueError("charged species cannot be classified as neutral_reactant")
    if any(species.charge_number == 0 and role[species.name] == "charge_carrier"
           for species in boundary.species):
        raise ValueError("charge_carrier species must carry nonzero charge")
    if charging_poisson_system is None and (
            initial_charge_node_c is not None or charging_options is not None):
        raise ValueError("charging state/options require charging_poisson_system")

    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    face_material = _face_material_ids(centroids, geometry)
    active_face = np.where(np.isin(face_material, etchable))[0]
    if active_face.size == 0:
        raise ValueError("current interface contains no requested etchable material")
    mesh_fingerprint = _surface_mesh_fingerprint(
        verts, faces, active_face, face_material, geometry)
    common_transport = dict(
        boundary=boundary, species_role=species_role, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
        device=transport_device)
    charging = None
    if charging_poisson_system is not None:
        if nodal_potential_v is not None:
            raise ValueError("self-consistent charging and a supplied nodal potential are exclusive")
        if potential_origin is None or potential_spacing is None or trajectory_fixed_dt is None:
            raise ValueError(
                "self-consistent charging requires potential_origin, potential_spacing, "
                "and trajectory_fixed_dt")
        initial_charge = (np.zeros(charging_poisson_system.shape)
                          if initial_charge_node_c is None
                          else np.asarray(initial_charge_node_c, dtype=float))
        options = {} if charging_options is None else dict(charging_options)
        if options.get("require_converged", True) is not True:
            raise ValueError("feature evolution requires a converged steady charging solve")
        charging = solve_dielectric_charging_steady_3d(
            charging_poisson_system, initial_charge, boundary, verts, faces, areas,
            source_bounds=source_bounds, source_z=source_z,
            potential_origin=potential_origin, potential_spacing=potential_spacing,
            mesh_length_unit_m=geometry.mesh_length_unit_m,
            mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
            trajectory_fixed_dt=trajectory_fixed_dt,
            trajectory_max_steps=trajectory_max_steps,
            transport_device=transport_device, **options)
        if not charging.converged:
            raise RuntimeError("feature evolution cannot consume a nonconverged charging field")
        transport = charging.transport
        uncharged_species = tuple(
            species for species in boundary.species if species.charge_number == 0)
        if uncharged_species:
            uncharged_boundary = PlasmaBoundaryState(
                uncharged_species, boundary.reference_plane_m, provenance=boundary.provenance)
            uncharged_role = {species.name: role[species.name] for species in uncharged_species}
            uncharged_transport = trace_boundary_state_field_3d(
                uncharged_boundary, uncharged_role, verts, faces, areas,
                source_bounds=source_bounds, source_z=source_z,
                nodal_potential_v=charging.potential_v,
                potential_origin=potential_origin, potential_spacing=potential_spacing,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
                fixed_dt=trajectory_fixed_dt, max_steps=trajectory_max_steps,
                device=transport_device)
            transport = merge_boundary_transport_results_3d(
                charging.transport, uncharged_transport)
    elif nodal_potential_v is None:
        if initial_charge_node_c is not None or charging_options is not None:
            raise ValueError("charging state/options require charging_poisson_system")
        if any(value is not None for value in (
                potential_origin, potential_spacing, trajectory_fixed_dt)):
            raise ValueError("field trajectory options require nodal_potential_v")
        transport = trace_boundary_state_first_hit_3d(**common_transport)
    else:
        if potential_origin is None or potential_spacing is None or trajectory_fixed_dt is None:
            raise ValueError(
                "nodal_potential_v requires potential_origin, potential_spacing, and trajectory_fixed_dt")
        transport = trace_boundary_state_field_3d(
            **common_transport, nodal_potential_v=nodal_potential_v,
            potential_origin=potential_origin, potential_spacing=potential_spacing,
            fixed_dt=trajectory_fixed_dt, max_steps=trajectory_max_steps)
    active_flux = _select_surface_fluxes(
        transport.surface_fluxes, active_face, len(faces), role)
    if surface_state is None:
        if surface_state_mesh_fingerprint is not None:
            raise ValueError("surface_state_mesh_fingerprint requires a supplied surface_state")
        if not hasattr(mechanism, "initial_state"):
            raise TypeError("surface mechanism must provide initial_state(shape)")
        surface_state = mechanism.initial_state((active_face.size,))
    else:
        if surface_state_mesh_fingerprint != mesh_fingerprint:
            raise ValueError(
                "surface_state mesh fingerprint mismatch; conservative remap is required")
        if not hasattr(surface_state, "conservative_surface_fields"):
            raise TypeError("surface state does not implement the conservative remap contract")
        state_fields = dict(surface_state.conservative_surface_fields())
        if not state_fields or any(
                np.asarray(value).shape != (active_face.size,) for value in state_fields.values()):
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
    transport_limitations = tuple(transport.known_limitations)
    if charging is not None:
        extra = tuple(
            limitation for limitation in transport.known_limitations
            if limitation not in charging.transport.known_limitations
            and limitation != "nodal potential is supplied rather than self-consistently charged")
        transport_limitations = tuple(charging.known_limitations) + extra
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(dict.fromkeys(transport_limitations)) + (
            "first-order material-local conservative surface-state remap",
            "topology-changing surface steps are refused",
            "first-order Godunov interface advection",
        ) + tuple(surface.validity.known_model_form_omissions))
    return FeatureStep3DResult(
        geometry=output_geometry, transport=transport, charging=charging, surface=surface,
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
            cfl_number=float(cfl_number), reinitialized=bool(reinitialize),
            self_consistent_charging=charging is not None,
            charging_iterations=(0 if charging is None else len(charging.history)),
            charging_converged=(None if charging is None else charging.converged)),
        validity=validity)


def solve_feature_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, n_steps, source_bounds, source_z,
        n_position=256, seed=0, cfl_number=0.3, reinitialize=True,
        transport_device=None, nodal_potential_v=None, potential_origin=None,
        potential_spacing=None, trajectory_fixed_dt=None, trajectory_max_steps=10000,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        charging_system_builder=None, initial_charge_node_c=None, charging_options=None):
    """Run verified feature steps with conserved surface state and optional quasi-static charging.

    A fixed ``charging_poisson_system`` is valid for one geometry only. Repeated charged evolution
    instead requires ``charging_system_builder(geometry)`` to rebuild the physical material operator
    after every interface update. Each geometry is independently converged from zero stored charge;
    this is the quasi-static charging limit, not a claim that transient surface charge was remapped.
    """
    if int(n_steps) != n_steps or int(n_steps) <= 0:
        raise ValueError("n_steps must be a positive integer")
    if not np.isfinite(duration_s) or duration_s < 0.0:
        raise ValueError("duration_s must be finite and nonnegative")
    if charging_poisson_system is not None and charging_system_builder is not None:
        raise ValueError("supply either a fixed charging system or a geometry-dependent builder")
    if charging_poisson_system is not None and int(n_steps) > 1:
        raise ValueError(
            "multi-step charged profile evolution requires a geometry-dependent Poisson builder")
    if charging_system_builder is not None and not callable(charging_system_builder):
        raise TypeError("charging_system_builder must be callable")
    step_duration = float(duration_s) / int(n_steps)
    current_geometry = geometry; current_state = None; current_fingerprint = None
    results = []
    for step_index in range(int(n_steps)):
        step_poisson_system = charging_poisson_system
        step_initial_charge = initial_charge_node_c
        if charging_system_builder is not None:
            step_poisson_system = charging_system_builder(current_geometry)
            if not isinstance(step_poisson_system, NodalPoissonSystem3D):
                raise TypeError("charging_system_builder must return NodalPoissonSystem3D")
            if step_poisson_system.shape != current_geometry.phi.shape:
                raise ValueError("rebuilt Poisson nodal grid must match the feature geometry")
            # A previous nodal charge grid is not a conservative representation on the moved surface.
            # In the quasi-static limit the new geometry owns a new independently converged root.
            if step_index > 0:
                step_initial_charge = np.zeros(step_poisson_system.shape)
        result = advance_feature_step_3d(
            current_geometry, boundary, species_role, mechanism,
            etchable_material_ids=etchable_material_ids, duration_s=step_duration,
            source_bounds=source_bounds, source_z=source_z,
            surface_state=current_state,
            surface_state_mesh_fingerprint=current_fingerprint,
            n_position=n_position, seed=int(seed) + step_index,
            nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
            potential_spacing=potential_spacing, trajectory_fixed_dt=trajectory_fixed_dt,
            trajectory_max_steps=trajectory_max_steps,
            charging_poisson_system=step_poisson_system,
            initial_charge_node_c=step_initial_charge,
            charging_options=charging_options,
            cfl_number=cfl_number, reinitialize=reinitialize,
            transport_device=transport_device)
        results.append(result)
        current_geometry = result.geometry
        current_state = result.next_surface_state
        current_fingerprint = result.next_surface_state_mesh_fingerprint
    reasons = tuple(reason for result in results for reason in result.validity.reasons)
    limitations = tuple(dict.fromkeys(
        limitation for result in results for limitation in result.validity.known_limitations))
    if charging_system_builder is not None:
        limitations += (
            "quasi-static charging re-solves each geometry independently; transient charge memory "
            "requires a conservative moving-surface charge equation",
        )
    return FeatureSolve3DResult(
        geometry=current_geometry, surface_state=current_state,
        surface_state_mesh_fingerprint=current_fingerprint,
        steps=tuple(results), duration_s=float(duration_s),
        validity=FeatureStepValidity(not reasons, reasons, limitations))
