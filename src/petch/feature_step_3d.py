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
from scipy.ndimage import label, map_coordinates
from scipy.spatial import cKDTree
from skimage.measure import euler_number

from .boundary_state import PlasmaBoundaryState
from .boundary_transport_3d import (
    BoundaryTransport3DResult,
    gather_boundary_state_ballistic_3d,
    estimate_diffuse_form_factors_3d,
    merge_boundary_transport_results_3d,
    trace_boundary_state_field_3d,
    trace_boundary_state_first_hit_3d,
)
from .neutral_radiosity_3d import solve_diffuse_neutral_radiosity_3d
from .charging_coupled_3d import (
    SteadyDielectricCharging3DResult, solve_dielectric_charging_steady_3d,
)
from .charging_poisson_3d import NodalPoissonSystem3D
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)
from .threed import advect_3d, extend_velocity_3d, extract_mesh_3d, reinit_cr2, reinit_fsm, reinit_narrow


@dataclass(frozen=True)
class FeatureGeometry3D:
    """Eulerian material geometry in declared mesh units; material id zero is gas."""

    phi: np.ndarray
    material_id: np.ndarray
    dx: float
    mesh_length_unit_m: float
    mesh_origin_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    material_levelsets: Mapping[int, np.ndarray] | None = None

    def __post_init__(self):
        phi = np.asarray(self.phi, dtype=float).copy()
        material = np.asarray(self.material_id, dtype=int).copy()
        origin = tuple(float(value) for value in self.mesh_origin_m)
        layers = (None if self.material_levelsets is None
                  else {int(key): np.asarray(value, dtype=float).copy()
                        for key, value in self.material_levelsets.items()})
        if (phi.ndim != 3 or min(phi.shape) < 2 or material.shape != phi.shape
                or np.any(~np.isfinite(phi)) or np.any(material < 0)
                or not np.isfinite(self.dx) or self.dx <= 0.0
                or not np.isfinite(self.mesh_length_unit_m) or self.mesh_length_unit_m <= 0.0
                or len(origin) != 3 or np.any(~np.isfinite(origin))):
            raise ValueError("invalid 3-D feature geometry")
        if not np.any(phi < 0.0) or not np.any(phi > 0.0):
            raise ValueError("phi must contain both gas and solid")
        if layers is not None:
            material_ids = set(np.unique(material)) - {0}
            if (not layers or set(layers) != material_ids or any(key <= 0 for key in layers)
                    or any(value.shape != phi.shape or np.any(~np.isfinite(value))
                           for value in layers.values())):
                raise ValueError("invalid material level-set fields")
            union = np.maximum.reduce(tuple(layers.values()))
            if np.any((union >= 0.0) != (phi >= 0.0)):
                raise ValueError("material level sets do not reconstruct the combined solid")
            for value in layers.values():
                value.setflags(write=False)
        phi.setflags(write=False); material.setflags(write=False)
        object.__setattr__(self, "phi", phi)
        object.__setattr__(self, "material_id", material)
        object.__setattr__(self, "dx", float(self.dx))
        object.__setattr__(self, "mesh_length_unit_m", float(self.mesh_length_unit_m))
        object.__setattr__(self, "mesh_origin_m", origin)
        object.__setattr__(
            self, "material_levelsets",
            None if layers is None else MappingProxyType(layers))

    @property
    def coordinate_arrays(self):
        return tuple(np.arange(size) * self.dx for size in self.phi.shape)


def make_rectangular_trench_geometry_3d(
        *, cell_width, cell_length, domain_height, dx, opening_width, mask_thickness,
        substrate_top, etched_depth, mesh_length_unit_m=1e-6,
        substrate_material_id=1, mask_material_id=2):
    """Construct a periodic-cell rectangular trench from units-explicit physical geometry.

    The trench is translationally invariant along the cell-length axis. ``etched_depth=0`` gives an
    unetched substrate under an open mask; positive depth creates vertical SiO2 sidewalls and a flat
    floor without a benchmark- or aspect-ratio-specific branch.
    """
    values = np.asarray([
        cell_width, cell_length, domain_height, dx, opening_width, mask_thickness,
        substrate_top, etched_depth, mesh_length_unit_m], dtype=float)
    if (np.any(~np.isfinite(values)) or np.any(values[:7] <= 0.0) or etched_depth < 0.0
            or opening_width >= cell_width or substrate_top <= etched_depth
            or substrate_top + mask_thickness >= domain_height
            or int(substrate_material_id) <= 0 or int(mask_material_id) <= 0
            or int(substrate_material_id) == int(mask_material_id)):
        raise ValueError("invalid rectangular trench geometry inputs")
    # ``phi`` is nodal: N intervals require N+1 nodes so the requested physical endpoint belongs to
    # the mesh.  Using only N nodes silently shortened every domain by one dx while source/wrap bounds
    # still used the requested length, leaving a periodic gap with no surface geometry.
    shape = tuple(max(3, int(round(length / dx)) + 1)
                  for length in (cell_width, cell_length, domain_height))
    x, y, z = (np.arange(size) * dx for size in shape)
    X, _, Z = np.meshgrid(x, y, z, indexing="ij")
    radius = np.abs(X - 0.5 * cell_width)
    floor = substrate_top - etched_depth
    base = floor - Z
    if etched_depth > 0.0:
        substrate_wall_slab = np.minimum(Z - floor, substrate_top - Z)
        substrate_wall = np.minimum(substrate_wall_slab, radius - 0.5 * opening_width)
        substrate_levelset = np.maximum(base, substrate_wall)
    else:
        substrate_levelset = substrate_top - Z
    mask_slab = np.minimum(Z - substrate_top, substrate_top + mask_thickness - Z)
    mask_levelset = np.minimum(mask_slab, radius - 0.5 * opening_width)
    substrate_phi = reinit_narrow(substrate_levelset, dx, domain_height + cell_width)
    mask_phi = reinit_narrow(mask_levelset, dx, domain_height + cell_width)
    analytic = np.maximum(substrate_phi, mask_phi)
    phi = reinit_narrow(analytic, dx, domain_height + cell_width)
    substrate_solid = (Z < substrate_top) & ~(
        (etched_depth > 0.0) & (Z > floor) & (radius < 0.5 * opening_width))
    mask_solid = ((Z >= substrate_top) & (Z < substrate_top + mask_thickness)
                  & (radius >= 0.5 * opening_width))
    material = np.zeros(shape, dtype=int)
    material[substrate_solid] = int(substrate_material_id)
    material[mask_solid] = int(mask_material_id)
    unlabeled_solid = (phi > 0.0) & (material == 0)
    # Reinitialization assigns exact-zero interface nodes to the positive (solid) side.  Material
    # ownership at those nodes must follow the CSG union winner, not a z threshold: at a flat mask
    # opening z==substrate_top belongs to the substrate level set, while the adjacent mask surface
    # belongs to the mask level set.
    substrate_owner = substrate_levelset >= mask_levelset
    material[unlabeled_solid] = np.where(
        substrate_owner[unlabeled_solid],
        int(substrate_material_id), int(mask_material_id))
    return FeatureGeometry3D(
        phi, material, dx, mesh_length_unit_m,
        material_levelsets={
            int(substrate_material_id): substrate_phi,
            int(mask_material_id): mask_phi,
        })


@dataclass(frozen=True)
class FeatureStepValidity:
    within_declared_scope: bool
    reasons: tuple[str, ...]
    known_limitations: tuple[str, ...]
    parameter_evidence_supports_prediction: bool
    nonpredictive_parameters: tuple[str, ...]


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
    """Assign each interface triangle by probing locally into its positive-phi solid.

    A global nearest-solid lookup is ambiguous at a material junction: an unetched substrate face
    inside a mask opening can be closer to the mask corner than to the next substrate grid node.
    The signed-distance gradient gives the physical solid-side normal and therefore the local owner.
    A nearest-solid search remains only as a fallback at degenerate zero-gradient CSG corners.
    """
    centroids = np.asarray(centroids, dtype=float)
    if geometry.material_levelsets is not None:
        material_ids = np.asarray(tuple(geometry.material_levelsets), dtype=int)
        coordinates = (centroids / geometry.dx).T
        values = np.vstack([
            map_coordinates(
                geometry.material_levelsets[int(material_id)], coordinates,
                order=1, mode="nearest", prefilter=False)
            for material_id in material_ids])
        return material_ids[np.argmax(values, axis=0)]
    solid = (geometry.phi > 0.0) & (geometry.material_id > 0)
    index = np.column_stack(np.where(solid))
    if index.size == 0:
        raise ValueError("geometry contains no labeled solid material")
    gradient = np.gradient(geometry.phi, geometry.dx)
    nearest_grid = np.rint(centroids / geometry.dx).astype(int)
    for axis in range(3):
        nearest_grid[:, axis] = np.clip(
            nearest_grid[:, axis], 0, geometry.phi.shape[axis] - 1)
    solid_normal = np.column_stack([
        component[tuple(nearest_grid.T)] for component in gradient])
    magnitude = np.linalg.norm(solid_normal, axis=1)
    valid_normal = magnitude > 1e-12
    solid_normal[valid_normal] /= magnitude[valid_normal, None]

    material = np.zeros(centroids.shape[0], dtype=int)
    unresolved = np.ones(centroids.shape[0], dtype=bool)
    for distance in (0.35, 0.75, 1.25):
        selected = np.where(unresolved & valid_normal)[0]
        if not selected.size:
            break
        probe = centroids[selected] + distance * geometry.dx * solid_normal[selected]
        probe_index = np.rint(probe / geometry.dx).astype(int)
        for axis in range(3):
            probe_index[:, axis] = np.clip(
                probe_index[:, axis], 0, geometry.phi.shape[axis] - 1)
        local_solid = geometry.phi[tuple(probe_index.T)] > 0.0
        local_material = geometry.material_id[tuple(probe_index.T)]
        accepted = local_solid & (local_material > 0)
        material[selected[accepted]] = local_material[accepted]
        unresolved[selected[accepted]] = False

    if np.any(unresolved):
        points = index * geometry.dx
        _, nearest = cKDTree(points).query(centroids[unresolved])
        chosen = index[np.asarray(nearest, dtype=int)]
        material[unresolved] = geometry.material_id[tuple(chosen.T)]
    return material


def _surface_gas_normals(verts, faces, centroids, geometry):
    triangle = np.asarray(verts)[np.asarray(faces)]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    gradient = np.gradient(geometry.phi, geometry.dx)
    index = np.rint(np.asarray(centroids) / geometry.dx).astype(int)
    for axis in range(3):
        index[:, axis] = np.clip(index[:, axis], 0, geometry.phi.shape[axis] - 1)
    # phi is positive in solid, so -grad(phi) points into gas.
    into_gas = -np.column_stack([
        component[tuple(index.T)] for component in gradient])
    flip = np.einsum("ij,ij->i", normal, into_gas) < 0.0
    normal[flip] *= -1.0
    return normal


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


def _physical_volume_topology_signature(geometry, etchable_material_ids):
    solid = (geometry.phi > 0.0) & np.isin(
        geometry.material_id, tuple(etchable_material_ids))
    _, components = label(solid)
    return int(components), int(euler_number(solid, connectivity=1))


def _physical_volume_component_sizes(geometry, etchable_material_ids):
    solid = (geometry.phi > 0.0) & np.isin(
        geometry.material_id, tuple(etchable_material_ids))
    component, count = label(solid)
    if count == 0:
        return ()
    return tuple(sorted(np.bincount(component.ravel())[1:].tolist(), reverse=True))


def _changed_physical_slice_topology(old_geometry, new_geometry, etchable_material_ids):
    """Locate a refused 3-D topology event without weakening the physical-volume gate."""
    materials = tuple(etchable_material_ids)
    old_solid = ((old_geometry.phi > 0.0)
                 & np.isin(old_geometry.material_id, materials))
    new_solid = ((new_geometry.phi > 0.0)
                 & np.isin(new_geometry.material_id, materials))
    changed = {}
    for axis, name in enumerate("xyz"):
        axis_changes = []
        for index in range(old_solid.shape[axis]):
            old_slice = np.take(old_solid, index, axis=axis)
            new_slice = np.take(new_solid, index, axis=axis)
            old_signature = (
                int(label(old_slice)[1]), int(euler_number(old_slice, connectivity=1)))
            new_signature = (
                int(label(new_slice)[1]), int(euler_number(new_slice, connectivity=1)))
            if old_signature != new_signature:
                axis_changes.append((index, old_signature, new_signature))
        if axis_changes:
            changed[name] = tuple(axis_changes[:12])
    return changed


def _remove_unresolved_subcell_solid_components(
        phi, material_id, etchable_material_ids, dx):
    updated = np.array(phi, copy=True)
    solid = (updated > 0.0) & np.isin(material_id, tuple(etchable_material_ids))
    component, count = label(solid)
    if count == 0:
        return updated, 0
    sizes = np.bincount(component.ravel())
    # Eight corner nodes are the minimum support of one resolved hexahedral volume cell.
    unresolved_label = np.flatnonzero(sizes < 8)
    unresolved_label = unresolved_label[unresolved_label != 0]
    if unresolved_label.size == 0:
        return updated, 0
    unresolved = np.isin(component, unresolved_label)
    # A subcell component has no resolved 3-D volume. Give it an unambiguous gas sign,
    # then let the signed-distance reconstruction restore a consistent neighborhood.
    updated[unresolved] = -np.maximum(np.abs(updated[unresolved]), float(dx))
    return updated, int(np.count_nonzero(unresolved))


def _redistance_feature_field(phi, dx, method):
    if method == "fsm":
        return reinit_fsm(phi, dx, 4.0 * dx)
    if method == "cr2":
        return reinit_cr2(phi, dx, 4.0 * dx)
    return reinit_narrow(phi, dx, 4.0 * dx)


def _advect_exposed_material_levelsets(
        material_levelsets, etchable_material_ids, extended_velocity,
        dx, duration_s, substeps):
    """Move each material only where its level set is the exposed union boundary."""
    current = {
        int(material_id): np.asarray(levelset, dtype=float).copy()
        for material_id, levelset in material_levelsets.items()}
    etchable = set(int(value) for value in etchable_material_ids)
    step_duration = float(duration_s) / int(substeps)
    for _ in range(int(substeps)):
        previous = current
        current = {}
        for material_id, levelset in previous.items():
            if material_id not in etchable:
                current[material_id] = levelset
                continue
            other_fields = tuple(
                value for key, value in previous.items() if key != material_id)
            exposed = (np.ones(levelset.shape, dtype=bool) if not other_fields
                       else levelset >= np.maximum.reduce(other_fields))
            current[material_id] = advect_3d(
                levelset, np.where(exposed, extended_velocity, 0.0),
                dx, step_duration)
    return current


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
        # A fixed O(dx) denominator regularization makes every remap smooth neighboring states by
        # O(1), even as the interface displacement tends to zero.  Repeating smaller time steps then
        # increases artificial diffusion and the method has no dt->0 limit.  Use only a roundoff-scale
        # floor: coincident predecessor faces map identically, while inverse-distance interpolation is
        # recovered when marching-cubes connectivity genuinely changes.
        coordinate_scale = max(
            float(dx), float(np.max(np.abs(old_centroid[old_index]))),
            float(np.max(np.abs(new_centroid[new_index]))), 1.0)
        distance_floor = 64.0 * np.finfo(float).eps * coordinate_scale
        exact = distance[:, 0] <= distance_floor
        weight = old_area[source_index] / np.maximum(distance * distance, distance_floor ** 2)
        if np.any(exact):
            weight[exact] = 0.0
            weight[exact, 0] = 1.0
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
                population.event_cosine_incidence[retained],
                event_position=(None if population.event_position is None
                                else population.event_position[retained]),
                event_incident_direction=(
                    None if population.event_incident_direction is None
                    else population.event_incident_direction[retained])))
        elif isinstance(population, EnergeticFlux):
            flux = np.asarray(population.flux_m2_s)
            selected_flux = flux if flux.ndim == 0 else flux[selected_face]
            energetic.append(EnergeticFlux(
                population.name, selected_flux, population.energy_eV,
                population.cosine_incidence, population.weight))
        else:  # pragma: no cover - SurfaceFluxes already validates this
            raise TypeError(type(population).__name__)
    return SurfaceFluxes(neutral, tuple(energetic))


def _apply_diffuse_neutral_transport(
        transport, geometry, verts, faces, centroids, areas, face_material, active_face,
        surface_state, mechanism, species_role, options, transport_device):
    options = dict(options)
    allowed = {
        "rays_per_face", "seed", "periodic_lateral", "domain_size", "ray_offset",
        "nonetchable_reaction_probability_by_material", "relative_tolerance",
        "maximum_iterations",
    }
    unknown = set(options) - allowed
    if unknown:
        raise ValueError("unknown neutral radiosity options: " + ", ".join(sorted(unknown)))
    material_probability = dict(options.pop(
        "nonetchable_reaction_probability_by_material", {}))
    solver_tolerance = float(options.pop("relative_tolerance", 1e-10))
    maximum_iterations = int(options.pop("maximum_iterations", 500))
    if "domain_size" not in options:
        options["domain_size"] = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    if "ray_offset" not in options:
        options["ray_offset"] = 1e-3 * geometry.dx
    factors = estimate_diffuse_form_factors_3d(
        verts, faces, centroids, _surface_gas_normals(
            verts, faces, centroids, geometry),
        device=transport_device, **options)
    if not hasattr(mechanism, "neutral_reaction_probability"):
        raise TypeError("diffuse neutral transport requires a mechanism reaction-probability contract")
    active_probability = dict(mechanism.neutral_reaction_probability(surface_state))
    neutral_names = [
        name for name, value in species_role.items() if value == "neutral_reactant"]
    reaction_probability = {}
    for name in neutral_names:
        probability = np.zeros(len(faces))
        if name in active_probability:
            value = np.asarray(active_probability[name], dtype=float)
            if value.shape != (active_face.size,):
                raise ValueError("mechanism neutral probability does not match active surface")
            probability[active_face] = value
        inactive = np.ones(len(faces), dtype=bool)
        inactive[active_face] = False
        for material in np.unique(face_material[inactive]):
            material_input = dict(material_probability.get(int(material), {}))
            if name not in material_input:
                raise ValueError(
                    f"missing neutral reaction probability for material {int(material)}, {name}")
            probability[inactive & (face_material == material)] = float(material_input[name])
        if np.any((probability < 0.0) | (probability > 1.0)):
            raise ValueError("material neutral reaction probabilities must lie in [0,1]")
        reaction_probability[name] = probability

    neutral_flux = {}
    diagnostics = {}
    physical_area = np.asarray(areas) * geometry.mesh_length_unit_m ** 2
    for name, direct in transport.surface_fluxes.neutral_flux_m2_s.items():
        if name not in reaction_probability:
            neutral_flux[name] = direct
            continue
        solution = solve_diffuse_neutral_radiosity_3d(
            direct, physical_area, factors.source_face, factors.target_face,
            factors.transfer_fraction, factors.escape_fraction,
            reaction_probability[name], relative_tolerance=solver_tolerance,
            maximum_iterations=maximum_iterations)
        neutral_flux[name] = solution.incident_flux_m2_s
        diagnostics[name] = dict(
            source_rate_s=solution.source_rate_s,
            reacted_rate_s=solution.reacted_rate_s,
            escaped_rate_s=solution.escaped_rate_s,
            relative_balance_error=solution.relative_balance_error,
            relative_linear_residual=solution.relative_linear_residual)
    limitations = tuple(
        item for item in transport.known_limitations
        if item != "no surface reflection or neutral re-emission") + (
        "neutral re-emission is diffuse with material/state reaction probabilities",
    )
    updated = BoundaryTransport3DResult(
        SurfaceFluxes(neutral_flux, transport.surface_fluxes.energetic_fluxes),
        transport.hit_probability, transport.escape_probability,
        transport.truncation_probability,
        transport.transport_model + " + flux_conservative_diffuse_radiosity",
        limitations, transport.lineage_replay_count)
    return updated, MappingProxyType(diagnostics)


def advance_feature_step_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, source_bounds, source_z,
        surface_state=None, n_position=256, seed=0,
        surface_state_mesh_fingerprint=None,
        nodal_potential_v=None, potential_origin=None, potential_spacing=None,
        trajectory_fixed_dt=None, trajectory_max_steps=10000,
        field_periodic_lateral=False,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        initial_charge_node_c=None, charging_options=None,
        precomputed_transport: BoundaryTransport3DResult | None = None,
        neutral_radiosity_options=None, ballistic_transport="forward",
        ballistic_face_quadrature_points=1, cfl_number=0.3, reinitialize=True,
        reinitialization_method="skfmm",
        transport_device=None):
    """Advance one stateful, dimensional feature step.

    The chemistry is evaluated only on triangles whose nearest positive-phi material id is in
    ``etchable_material_ids``. Other labeled solids are pinned. The method refuses a supplied surface
    state whose shape does not match the current active mesh; it never silently remaps history.
    ``precomputed_transport`` lets an orchestrating physical-time charging driver reuse its final
    exact charged/re-impact measure for chemistry instead of retracing a second kinetic operator.
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
    if precomputed_transport is not None and not isinstance(
            precomputed_transport, BoundaryTransport3DResult):
        raise TypeError("precomputed_transport must be BoundaryTransport3DResult")
    if precomputed_transport is not None and (
            charging_poisson_system is not None or nodal_potential_v is not None):
        raise ValueError(
            "precomputed transport is exclusive with an internally evaluated charging/field path")
    if ballistic_transport not in ("forward", "face_gather"):
        raise ValueError("ballistic_transport must be 'forward' or 'face_gather'")
    if reinitialization_method not in ("skfmm", "fsm", "cr2"):
        raise ValueError("reinitialization_method must be 'skfmm', 'fsm', or 'cr2'")
    if ballistic_transport == "face_gather" and (
            charging_poisson_system is not None or nodal_potential_v is not None):
        raise ValueError("deterministic ballistic face gather does not yet trace electric fields")

    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    face_material = _face_material_ids(centroids, geometry)
    active_face = np.where(np.isin(face_material, etchable))[0]
    if active_face.size == 0:
        raise ValueError("current interface contains no requested etchable material")
    mesh_fingerprint = _surface_mesh_fingerprint(
        verts, faces, active_face, face_material, geometry)
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
    radiosity_options = (None if neutral_radiosity_options is None
                         else dict(neutral_radiosity_options))
    periodic_neutral = bool(
        radiosity_options is not None and radiosity_options.get("periodic_lateral", False))
    charging_periodic = bool(
        charging_options is not None and charging_options.get("periodic_lateral", False))
    if periodic_neutral and (charging_poisson_system is not None or nodal_potential_v is not None):
        if not (charging_periodic if charging_poisson_system is not None
                else bool(field_periodic_lateral)):
            raise ValueError(
                "periodic neutral radiosity with a field requires periodic charged trajectories")
    common_transport = dict(
        boundary=boundary, species_role=species_role, verts=verts, faces=faces, areas=areas,
        source_bounds=source_bounds, source_z=source_z,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
        device=transport_device)
    face_gas_normals = _surface_gas_normals(verts, faces, centroids, geometry)
    charging = None
    if precomputed_transport is not None:
        transport = precomputed_transport
        available = set(transport.surface_fluxes.neutral_flux_m2_s)
        available.update(
            population.name for population in transport.surface_fluxes.energetic_fluxes)
        expected = {species.name for species in boundary.species}
        if not expected.issubset(available):
            raise ValueError(
                "precomputed transport omits boundary species: "
                + ", ".join(sorted(expected - available)))
        for name, value in transport.surface_fluxes.neutral_flux_m2_s.items():
            if np.asarray(value).shape != (len(faces),):
                raise ValueError(
                    f"precomputed neutral flux {name!r} does not match the current surface mesh")
        for population in transport.surface_fluxes.energetic_fluxes:
            if (isinstance(population, FaceResolvedEnergeticFlux)
                    and population.face_count != len(faces)):
                raise ValueError(
                    f"precomputed energetic flux {population.name!r} uses another surface mesh")
    elif charging_poisson_system is not None:
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
            face_centroids=centroids,
            face_gas_normals=_surface_gas_normals(verts, faces, centroids, geometry),
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
                periodic_lateral=charging_periodic,
                face_gas_normals=face_gas_normals,
                device=transport_device)
            transport = merge_boundary_transport_results_3d(
                charging.transport, uncharged_transport)
    elif nodal_potential_v is None:
        if initial_charge_node_c is not None or charging_options is not None:
            raise ValueError("charging state/options require charging_poisson_system")
        if any(value is not None for value in (
                potential_origin, potential_spacing, trajectory_fixed_dt)):
            raise ValueError("field trajectory options require nodal_potential_v")
        first_hit_options = {}
        if periodic_neutral:
            first_hit_options = dict(
                periodic_lateral=True,
                domain_size=radiosity_options.get(
                    "domain_size", (np.asarray(geometry.phi.shape) - 1) * geometry.dx))
        if ballistic_transport == "face_gather":
            transport = gather_boundary_state_ballistic_3d(
                boundary, species_role, verts, faces, areas, centroids,
                _surface_gas_normals(verts, faces, centroids, geometry),
                source_bounds=source_bounds, source_z=source_z,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                mesh_origin_m=geometry.mesh_origin_m,
                face_quadrature_points=ballistic_face_quadrature_points,
                periodic_lateral=periodic_neutral,
                domain_size=first_hit_options.get("domain_size"),
                ray_offset=1e-3 * geometry.dx, device=transport_device)
        else:
            transport = trace_boundary_state_first_hit_3d(
                **common_transport, **first_hit_options)
    else:
        if potential_origin is None or potential_spacing is None or trajectory_fixed_dt is None:
            raise ValueError(
                "nodal_potential_v requires potential_origin, potential_spacing, and trajectory_fixed_dt")
        transport = trace_boundary_state_field_3d(
            **common_transport, nodal_potential_v=nodal_potential_v,
            potential_origin=potential_origin, potential_spacing=potential_spacing,
            fixed_dt=trajectory_fixed_dt, max_steps=trajectory_max_steps,
            periodic_lateral=bool(field_periodic_lateral),
            face_gas_normals=face_gas_normals)
    neutral_radiosity_diagnostics = MappingProxyType({})
    if radiosity_options is not None:
        transport, neutral_radiosity_diagnostics = _apply_diffuse_neutral_transport(
            transport, geometry, verts, faces, centroids, areas, face_material, active_face,
            surface_state, mechanism, role, radiosity_options, transport_device)
    active_flux = _select_surface_fluxes(
        transport.surface_fluxes, active_face, len(faces), role)
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
    # Extend only from the material surface that is actually evolving.  Including pinned mask
    # triangles with zero velocity lets them win the nearest-face query below a narrow opening and
    # numerically pins a physically bombarded floor after roughly one grid cell of motion.
    extended_velocity = extend_velocity_3d(
        face_velocity[active_face], centroids[active_face],
        extension_geometry, 4.0 * geometry.dx)
    center = (geometry.phi.shape[0] // 2, geometry.phi.shape[1] // 2)
    centerline = geometry.phi[center]
    center_crossing = np.flatnonzero(
        (centerline[:-1] >= 0.0) & (centerline[1:] < 0.0))
    center_diagnostics = {}
    if center_crossing.size == 1:
        lower = int(center_crossing[0])
        fraction = centerline[lower] / (centerline[lower] - centerline[lower + 1])
        center_diagnostics = dict(
            centerline_interface_lower_index=lower,
            centerline_interface_fraction=float(fraction),
            centerline_extended_velocity_mesh_units_s=float(
                (1.0 - fraction) * extended_velocity[center + (lower,)]
                + fraction * extended_velocity[center + (lower + 1,)]),
            centerline_phi_lower_before=float(centerline[lower]),
            centerline_phi_upper_before=float(centerline[lower + 1]))
    pinned = (geometry.material_id > 0) & ~np.isin(geometry.material_id, etchable)
    material_levelsets = None
    if geometry.material_levelsets is None:
        for _ in range(substeps):
            phi = advect_3d(
                phi, extended_velocity, geometry.dx, float(duration_s) / substeps)
            phi[pinned] = geometry.phi[pinned]
    else:
        material_levelsets = _advect_exposed_material_levelsets(
            geometry.material_levelsets, etchable, extended_velocity,
            geometry.dx, duration_s, substeps)
        phi = np.maximum.reduce(tuple(material_levelsets.values()))
    advected_centerline = phi[center]
    advected_crossing = np.flatnonzero(
        (advected_centerline[:-1] >= 0.0) & (advected_centerline[1:] < 0.0))
    if advected_crossing.size == 1:
        lower = int(advected_crossing[0])
        center_diagnostics["centerline_advected_interface_fraction"] = float(
            advected_centerline[lower]
            / (advected_centerline[lower] - advected_centerline[lower + 1]))
    if reinitialize and duration_s > 0.0:
        if material_levelsets is not None:
            material_levelsets = {
                material_id: (
                    _redistance_feature_field(levelset, geometry.dx, reinitialization_method)
                    if material_id in etchable else levelset)
                for material_id, levelset in material_levelsets.items()}
            phi = np.maximum.reduce(tuple(material_levelsets.values()))
        phi = _redistance_feature_field(phi, geometry.dx, reinitialization_method)
        if material_levelsets is None:
            phi[pinned] = geometry.phi[pinned]
    reinitialized_centerline = phi[center]
    reinitialized_crossing = np.flatnonzero(
        (reinitialized_centerline[:-1] >= 0.0) & (reinitialized_centerline[1:] < 0.0))
    if reinitialized_crossing.size == 1:
        lower = int(reinitialized_crossing[0])
        center_diagnostics["centerline_reinitialized_interface_fraction"] = float(
            reinitialized_centerline[lower]
            / (reinitialized_centerline[lower] - reinitialized_centerline[lower + 1]))
    phi, removed_unresolved_solid_cells = _remove_unresolved_subcell_solid_components(
        phi, geometry.material_id, etchable, geometry.dx)
    if removed_unresolved_solid_cells:
        if material_levelsets is not None:
            raise RuntimeError(
                "subcell cleanup requires an explicit material-layer topology update")
        phi = _redistance_feature_field(phi, geometry.dx, reinitialization_method)
        phi[pinned] = geometry.phi[pinned]

    output_geometry = FeatureGeometry3D(
        phi, geometry.material_id, geometry.dx, geometry.mesh_length_unit_m,
        geometry.mesh_origin_m, material_levelsets=material_levelsets)
    next_verts, next_faces, next_centroids, next_areas = extract_mesh_3d(
        output_geometry.phi, output_geometry.dx)
    next_face_material = _face_material_ids(next_centroids, output_geometry)
    next_active_face = np.where(np.isin(next_face_material, etchable))[0]
    if next_active_face.size == 0:
        raise ValueError("etch step removed every requested material surface")
    old_mesh_topology = _surface_topology_signature(faces, active_face)
    next_mesh_topology = _surface_topology_signature(next_faces, next_active_face)
    old_topology = _physical_volume_topology_signature(geometry, etchable)
    next_topology = _physical_volume_topology_signature(output_geometry, etchable)
    if old_topology != next_topology:
        raise ValueError(
            f"surface topology changed from {old_topology} to {next_topology}; "
            f"component sizes changed from "
            f"{_physical_volume_component_sizes(geometry, etchable)} to "
            f"{_physical_volume_component_sizes(output_geometry, etchable)}; "
            f"changed slice topology="
            f"{_changed_physical_slice_topology(geometry, output_geometry, etchable)}; "
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
        old_mesh_topology=old_mesh_topology, new_mesh_topology=next_mesh_topology,
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
    material_exchange = getattr(surface, "material_exchange", None)
    if material_exchange is None:
        exchange_limitations = ("surface mechanism does not expose a material-exchange ledger",)
        product_routing_complete = None
    else:
        exchange_limitations = tuple(material_exchange.known_limitations)
        product_routing_complete = bool(material_exchange.product_routing_complete)
    product_populations = tuple(getattr(surface, "product_populations", ()))
    outgoing_material = bool(
        material_exchange is not None
        and any(np.any(value > 0.0) for value in material_exchange.outgoing_units_m2.values()))
    if not outgoing_material:
        product_transport_ready = None
    elif not product_populations:
        product_transport_ready = False
        exchange_limitations += (
            "outgoing material has no declared surface-product populations",)
    else:
        product_transport_ready = all(item.transport_ready for item in product_populations)
        if not product_transport_ready:
            exchange_limitations += (
                "surface-product populations lack a complete energy/angular launch model",)
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(dict.fromkeys(transport_limitations)) + (
            "first-order material-local conservative surface-state remap",
            "physical volume-topology-changing surface steps are refused",
            "first-order Godunov interface advection",
        ) + tuple(surface.validity.known_model_form_omissions) + exchange_limitations,
        parameter_evidence_supports_prediction=(
            surface.validity.parameter_evidence_supports_prediction),
        nonpredictive_parameters=surface.validity.nonpredictive_parameters)
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
            reinitialization_method=(reinitialization_method if reinitialize else None),
            removed_unresolved_solid_cells=removed_unresolved_solid_cells,
            self_consistent_charging=charging is not None,
            charging_iterations=(0 if charging is None else len(charging.history)),
            charging_converged=(None if charging is None else charging.converged),
            product_routing_complete=product_routing_complete,
            product_population_count=len(product_populations),
            product_transport_ready=product_transport_ready,
            neutral_radiosity=neutral_radiosity_diagnostics,
            **center_diagnostics),
        validity=validity)


def solve_feature_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, n_steps, source_bounds, source_z,
        n_position=256, seed=0, cfl_number=0.3, reinitialize=True,
        transport_device=None, nodal_potential_v=None, potential_origin=None,
        potential_spacing=None, trajectory_fixed_dt=None, trajectory_max_steps=10000,
        field_periodic_lateral=False,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        charging_system_builder=None, initial_charge_node_c=None, charging_options=None,
        neutral_radiosity_options=None, ballistic_transport="forward",
        ballistic_face_quadrature_points=1, reinitialization_method="skfmm"):
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
        try:
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
                field_periodic_lateral=field_periodic_lateral,
                charging_poisson_system=step_poisson_system,
                initial_charge_node_c=step_initial_charge,
                charging_options=charging_options,
                neutral_radiosity_options=neutral_radiosity_options,
                ballistic_transport=ballistic_transport,
                ballistic_face_quadrature_points=ballistic_face_quadrature_points,
                reinitialization_method=reinitialization_method,
                cfl_number=cfl_number, reinitialize=reinitialize,
                transport_device=transport_device)
        except (ValueError, RuntimeError) as error:
            raise type(error)(f"feature step {step_index + 1}/{int(n_steps)}: {error}") from error
        results.append(result)
        current_geometry = result.geometry
        current_state = result.next_surface_state
        current_fingerprint = result.next_surface_state_mesh_fingerprint
    reasons = tuple(reason for result in results for reason in result.validity.reasons)
    limitations = tuple(dict.fromkeys(
        limitation for result in results for limitation in result.validity.known_limitations))
    nonpredictive = tuple(dict.fromkeys(
        name for result in results for name in result.validity.nonpredictive_parameters))
    if charging_system_builder is not None:
        limitations += (
            "quasi-static charging re-solves each geometry independently; transient charge memory "
            "requires a conservative moving-surface charge equation",
        )
    return FeatureSolve3DResult(
        geometry=current_geometry, surface_state=current_state,
        surface_state_mesh_fingerprint=current_fingerprint,
        steps=tuple(results), duration_s=float(duration_s),
        validity=FeatureStepValidity(
            not reasons, reasons, limitations, not nonpredictive, nonpredictive))
