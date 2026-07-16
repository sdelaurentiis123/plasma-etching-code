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
from .charged_surface_cascade_3d import (
    ChargedSurfaceCascade3DResult,
    apply_charged_surface_response_to_transport_3d,
)
from .charged_surface_response_3d import ChargedSurfaceContext3D
from .hwang_giapis_scatter_3d import (
    HwangGiapisForwardScatter3DResult, HwangGiapisSiO2ForwardScatter3D,
    apply_hwang_giapis_forward_scatter_to_transport_3d,
)
from .surface_exchange import SurfaceProductPopulation
from .surface_product_redeposition_3d import (
    SurfaceProductRedepositionContract3D,
    transport_surface_product_redeposition_3d,
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
    charged_surface_cascade: ChargedSurfaceCascade3DResult | None
    neutral_forward_scatter: HwangGiapisForwardScatter3DResult | None
    surface_product_redeposition: object | None
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

    def trilinear(field, point):
        coordinate = np.asarray(point, dtype=float) / geometry.dx
        lower = np.floor(coordinate).astype(int)
        for axis in range(3):
            lower[:, axis] = np.clip(lower[:, axis], 0, field.shape[axis] - 2)
        fraction = np.clip(coordinate - lower, 0.0, 1.0)
        value = np.zeros(len(point))
        for ox in (0, 1):
            wx = fraction[:, 0] if ox else 1.0 - fraction[:, 0]
            for oy in (0, 1):
                wy = fraction[:, 1] if oy else 1.0 - fraction[:, 1]
                for oz in (0, 1):
                    wz = fraction[:, 2] if oz else 1.0 - fraction[:, 2]
                    index = lower + np.array([ox, oy, oz])
                    value += wx * wy * wz * field[tuple(index.T)]
        return value

    centroid = np.asarray(centroids, dtype=float)
    probe_distance = 0.25 * geometry.dx
    domain_maximum = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    plus = np.clip(centroid + probe_distance * normal, 0.0, domain_maximum)
    minus = np.clip(centroid - probe_distance * normal, 0.0, domain_maximum)
    signed_difference = trilinear(geometry.phi, plus) - trilinear(geometry.phi, minus)
    # phi is positive in solid.  A positive signed difference means ``normal`` points into solid.
    flip = signed_difference > 0.0
    ambiguous = np.abs(signed_difference) <= 64.0 * np.finfo(float).eps * geometry.dx
    if np.any(ambiguous):
        gradient = np.gradient(geometry.phi, geometry.dx)
        index = np.rint(centroid[ambiguous] / geometry.dx).astype(int)
        for axis in range(3):
            index[:, axis] = np.clip(index[:, axis], 0, geometry.phi.shape[axis] - 1)
        into_gas = -np.column_stack([
            component[tuple(index.T)] for component in gradient])
        flip[ambiguous] = (
            np.einsum("ij,ij->i", normal[ambiguous], into_gas) < 0.0)
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
        return updated, 0, np.zeros(updated.shape, dtype=bool)
    sizes = np.bincount(component.ravel())
    # Eight corner nodes are the minimum support of one resolved hexahedral volume cell.
    unresolved_label = np.flatnonzero(sizes < 8)
    unresolved_label = unresolved_label[unresolved_label != 0]
    if unresolved_label.size == 0:
        return updated, 0, np.zeros(updated.shape, dtype=bool)
    unresolved = np.isin(component, unresolved_label)
    # A subcell component has no resolved 3-D volume. Give it an unambiguous gas sign,
    # then let the signed-distance reconstruction restore a consistent neighborhood.
    updated[unresolved] = -np.maximum(np.abs(updated[unresolved]), float(dx))
    return updated, int(np.count_nonzero(unresolved)), unresolved


def _apply_subcell_cleanup_to_material_levelsets(
        material_levelsets, removal_mask, owner_material_id, etchable_material_ids,
        dx, reinitialization_method, periodic_lateral):
    """Apply a combined-surface subcell removal to its authoritative material layers."""
    removal = np.asarray(removal_mask, dtype=bool)
    owner = np.asarray(owner_material_id)
    if removal.shape != owner.shape or not np.any(removal):
        raise ValueError("material-layer cleanup requires a nonempty owning removal mask")
    updated = {int(material_id): np.asarray(levelset, dtype=float).copy()
               for material_id, levelset in material_levelsets.items()}
    etchable = tuple(int(value) for value in etchable_material_ids)
    if any(levelset.shape != owner.shape for levelset in updated.values()):
        raise ValueError("material level sets must share the combined surface shape")
    accounted = np.zeros(removal.shape, dtype=bool)
    for material_id in etchable:
        selected = removal & (owner == material_id)
        if not np.any(selected):
            continue
        if material_id not in updated:
            raise ValueError("subcell removal references a missing material level set")
        levelset = updated[material_id]
        levelset[selected] = -np.maximum(np.abs(levelset[selected]), float(dx))
        updated[material_id] = _redistance_feature_field(
            levelset, dx, reinitialization_method,
            periodic_lateral=periodic_lateral)
        accounted |= selected
    if not np.array_equal(accounted, removal):
        raise RuntimeError("subcell removal includes a non-etchable or unowned material node")
    if periodic_lateral:
        updated = {
            material_id: _project_periodic_lateral_endpoints(levelset)[0]
            for material_id, levelset in updated.items()}
    combined = np.maximum.reduce(tuple(updated.values()))
    combined = _redistance_feature_field(
        combined, dx, reinitialization_method,
        periodic_lateral=periodic_lateral)
    material_ids = np.asarray(sorted(updated), dtype=int)
    material_stack = np.stack([updated[int(material_id)] for material_id in material_ids])
    combined_owner = material_ids[np.argmax(material_stack, axis=0)]
    combined_owner = np.where(combined >= 0.0, combined_owner, 0)
    return updated, combined, combined_owner


def _project_periodic_lateral_endpoints(field):
    """Project duplicate x/y endpoint planes onto one nodal-periodic field."""
    output = np.asarray(field, dtype=float).copy()
    maximum_correction = 0.0
    for axis in (0, 1):
        first = [slice(None)] * output.ndim; first[axis] = 0
        last = [slice(None)] * output.ndim; last[axis] = -1
        first = tuple(first); last = tuple(last)
        seam = 0.5 * (output[first] + output[last])
        maximum_correction = max(
            maximum_correction,
            float(np.max(np.abs(output[first] - seam))),
            float(np.max(np.abs(output[last] - seam))))
        output[first] = seam
        output[last] = seam
    return output, maximum_correction


def _periodic_lateral_redistance(phi, dx, method):
    """Redistance a duplicate-endpoint periodic field through wrapped lateral padding."""
    field, projection = _project_periodic_lateral_endpoints(phi)
    core = field[:-1, :-1, :]
    padding = int(np.ceil(4.0 * dx / dx)) + 2
    padded = np.pad(core, ((padding, padding), (padding, padding), (0, 0)), mode="wrap")
    if method == "fsm":
        redistanced = reinit_fsm(padded, dx, 4.0 * dx)
    elif method == "cr2":
        redistanced = reinit_cr2(padded, dx, 4.0 * dx)
    else:
        redistanced = reinit_narrow(padded, dx, 4.0 * dx)
    cropped = redistanced[padding:-padding, padding:-padding, :]
    output = np.empty_like(field)
    output[:-1, :-1, :] = cropped
    output[-1, :-1, :] = cropped[0, :, :]
    output[:-1, -1, :] = cropped[:, 0, :]
    output[-1, -1, :] = cropped[0, 0, :]
    output, final_projection = _project_periodic_lateral_endpoints(output)
    return output, max(projection, final_projection)


def _redistance_feature_field(phi, dx, method, *, periodic_lateral=False):
    if periodic_lateral:
        return _periodic_lateral_redistance(phi, dx, method)[0]
    if method == "fsm":
        return reinit_fsm(phi, dx, 4.0 * dx)
    if method == "cr2":
        return reinit_cr2(phi, dx, 4.0 * dx)
    return reinit_narrow(phi, dx, 4.0 * dx)


def _advect_exposed_material_levelsets(
        material_levelsets, etchable_material_ids, extended_velocity,
        dx, duration_s, substeps, *, periodic_lateral=False):
    """Move each material only where its level set is the exposed union boundary."""
    current = {}
    for material_id, levelset in material_levelsets.items():
        value = np.asarray(levelset, dtype=float).copy()
        if periodic_lateral:
            value = _project_periodic_lateral_endpoints(value)[0]
        current[int(material_id)] = value
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
                dx, step_duration,
                periodic_axes=((0, 1) if periodic_lateral else ()))
            if periodic_lateral:
                current[material_id] = _project_periodic_lateral_endpoints(
                    current[material_id])[0]
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
    """First-order material-local remap with declared intensive/conservative semantics.

    The state declares named nonnegative fields, optional upper bounds, and reconstruction. Interpolation
    supplies spatial locality. A constrained correction preserves each field marked ``conservative``;
    algebraic coverage fractions may explicitly declare ``intensive`` remap and are interpolated without
    inventing an area-integral conservation law. The default remains conservative for legacy states.
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
    remap_modes = (
        dict(state.surface_field_remap_modes())
        if hasattr(state, "surface_field_remap_modes")
        else {name: "conservative" for name in old_values})
    if (not old_values or set(upper_bounds) != set(old_values)
            or set(remap_modes) != set(old_values)
            or any(mode not in {"conservative", "intensive"}
                   for mode in remap_modes.values())):
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
        targets = {}; residuals = []; remapped_integrals = {}
        for field_name, old_value in old_values.items():
            raw = np.sum(weight * old_value[source_index], axis=1)
            target = float(np.dot(old_value[old_index], old_area[old_index]))
            if remap_modes[field_name] == "conservative":
                remapped = _conserve_nonnegative_surface_field(
                    raw, target, new_area[new_index], upper=upper_bounds[field_name])
            else:
                remapped = np.maximum(raw, 0.0)
                if upper_bounds[field_name] is not None:
                    remapped = np.minimum(remapped, upper_bounds[field_name])
            output[field_name][new_index] = remapped
            achieved = float(np.dot(remapped, new_area[new_index]))
            if remap_modes[field_name] == "conservative":
                residuals.append(abs(achieved - target) / max(abs(target), 1.0))
                targets[field_name] = target * physical_area_scale
            remapped_integrals[field_name] = achieved * physical_area_scale
        material_diagnostics[int(material)] = dict(
            old_face_count=int(old_index.size), new_face_count=int(new_index.size),
            old_area_m2=float(old_area[old_index].sum() * physical_area_scale),
            new_area_m2=float(new_area[new_index].sum() * physical_area_scale),
            target_field_integrals=targets,
            remapped_field_integrals=remapped_integrals,
            field_remap_modes=dict(remap_modes),
            max_relative_conservation_residual=float(max(residuals, default=0.0)))
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
    if hasattr(mechanism, "neutral_reaction_probability_by_material"):
        active_probability = dict(mechanism.neutral_reaction_probability_by_material(
            surface_state, face_material[active_face]))
    elif hasattr(mechanism, "neutral_reaction_probability"):
        active_probability = dict(mechanism.neutral_reaction_probability(surface_state))
    else:
        raise TypeError(
            "diffuse neutral transport requires a mechanism reaction-probability contract")
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
        limitations, lineage_replay_count=transport.lineage_replay_count,
        lineage_replay_eligible_count=transport.lineage_replay_eligible_count,
        edge_launch_inset_count=transport.edge_launch_inset_count,
        trajectory_horizon_extension_count=(
            transport.trajectory_horizon_extension_count),
        trajectory_initial_max_steps=transport.trajectory_initial_max_steps,
        trajectory_final_max_steps=transport.trajectory_final_max_steps,
        trajectory_emergency_max_steps=transport.trajectory_emergency_max_steps)
    return updated, MappingProxyType(diagnostics)


def _apply_surface_product_redeposition(
        populations, geometry, verts, faces, centroids, areas, face_material,
        active_face, duration_s, options, transport_device):
    """Run the opt-in same-material product return path on the complete surface mesh."""
    options = dict(options)
    allowed = {
        "contract", "rays_per_face", "seed", "periodic_lateral", "domain_size",
        "ray_offset", "relative_tolerance", "maximum_iterations",
    }
    unknown = set(options) - allowed
    if unknown:
        raise ValueError(
            "unknown surface-product redeposition options: " + ", ".join(sorted(unknown)))
    contract = options.pop("contract", None)
    if not isinstance(contract, SurfaceProductRedepositionContract3D):
        raise TypeError("surface-product redeposition requires an explicit contract")
    relative_tolerance = float(options.pop("relative_tolerance", 1e-10))
    maximum_iterations = int(options.pop("maximum_iterations", 500))
    if "domain_size" not in options:
        options["domain_size"] = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
    if "ray_offset" not in options:
        options["ray_offset"] = 1e-3 * geometry.dx
    factors = estimate_diffuse_form_factors_3d(
        verts, faces, centroids, _surface_gas_normals(
            verts, faces, centroids, geometry),
        device=transport_device, **options)
    full_populations = []
    for population in tuple(populations):
        local_count = np.asarray(population.integrated_particle_count_m2, dtype=float)
        if local_count.shape != (active_face.size,):
            raise ValueError(
                f"surface product {population.name!r} does not match the active surface")
        count = np.zeros(len(faces))
        count[active_face] = local_count
        full_populations.append(SurfaceProductPopulation(
            population.name, population.source_inventory, count,
            population.material_units_per_particle, population.mass_amu,
            angular_model=population.angular_model, energy_model=population.energy_model,
            energy_parameters=population.energy_parameters, provenance=population.provenance,
            relative_standard_uncertainty=population.relative_standard_uncertainty))
    evolving = np.zeros(len(faces), dtype=bool)
    evolving[active_face] = True
    physical_area = np.asarray(areas) * geometry.mesh_length_unit_m ** 2
    return transport_surface_product_redeposition_3d(
        full_populations, float(duration_s), physical_area, factors, face_material, evolving,
        contract, relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations)


def advance_feature_step_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, source_bounds, source_z,
        surface_state=None, n_position=256, seed=0,
        surface_state_mesh_fingerprint=None,
        nodal_potential_v=None, potential_origin=None, potential_spacing=None,
        trajectory_fixed_dt=None, trajectory_max_steps=10000,
        trajectory_adaptive_horizon=False, trajectory_emergency_max_steps=None,
        field_periodic_lateral=False, profile_periodic_lateral=None,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        initial_charge_node_c=None, charging_options=None,
        precomputed_transport: BoundaryTransport3DResult | None = None,
        charged_surface_response=None, charged_surface_response_options=None,
        neutral_forward_scatter=None, neutral_forward_scatter_options=None,
        neutral_radiosity_options=None,
        neutral_surface_fixed_point_tolerance=None,
        neutral_surface_fixed_point_max_iterations=20,
        surface_product_redeposition_options=None,
        ballistic_transport="forward",
        ballistic_face_quadrature_points=1, cfl_number=0.3, reinitialize=True,
        reinitialization_method="skfmm",
        transport_device=None):
    """Advance one stateful, dimensional feature step.

    The chemistry is evaluated only on triangles whose nearest positive-phi material id is in
    ``etchable_material_ids``. Other labeled solids are pinned. Multiple evolving materials require a
    material-resolved mechanism router; one substrate law is never silently applied to a mask. The
    method refuses a supplied surface state whose shape does not match the current active mesh; it never
    silently remaps history.
    ``precomputed_transport`` lets an orchestrating physical-time charging driver reuse its final
    exact charged/re-impact measure for chemistry instead of retracing a second kinetic operator.
    ``charged_surface_response`` applies the same certified common-engine response/re-impact cascade
    to an ordinary supplied-field or explicitly field-free feature step.  It is exclusive with
    precomputed and self-consistent charging transports so the energetic lineage cannot be applied
    twice.
    """
    if not np.isfinite(duration_s) or duration_s < 0.0:
        raise ValueError("duration_s must be finite and nonnegative")
    if not np.isfinite(cfl_number) or not 0.0 < cfl_number < 1.0:
        raise ValueError("cfl_number must lie strictly between zero and one")
    if (int(trajectory_max_steps) != trajectory_max_steps or trajectory_max_steps <= 0
            or not isinstance(trajectory_adaptive_horizon, (bool, np.bool_))
            or (trajectory_emergency_max_steps is not None
                and (int(trajectory_emergency_max_steps) != trajectory_emergency_max_steps
                     or trajectory_emergency_max_steps < trajectory_max_steps))
            or (trajectory_adaptive_horizon and trajectory_emergency_max_steps is None)):
        raise ValueError("invalid feature-step trajectory-horizon controls")
    if (profile_periodic_lateral is not None
            and not isinstance(profile_periodic_lateral, (bool, np.bool_))):
        raise ValueError("profile_periodic_lateral must be boolean or None")
    if neutral_surface_fixed_point_tolerance is not None:
        if (not np.isfinite(neutral_surface_fixed_point_tolerance)
                or not 0.0 < neutral_surface_fixed_point_tolerance < 1.0
                or int(neutral_surface_fixed_point_max_iterations)
                != neutral_surface_fixed_point_max_iterations
                or neutral_surface_fixed_point_max_iterations <= 0):
            raise ValueError("invalid neutral/surface fixed-point controls")
        if neutral_radiosity_options is None:
            raise ValueError(
                "neutral/surface fixed point requires diffuse neutral radiosity")
        if not getattr(mechanism, "quasi_steady_surface_state", False):
            raise ValueError(
                "surface mechanism does not declare a quasi-steady neutral/surface state")
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
    if charged_surface_response is None and charged_surface_response_options is not None:
        raise ValueError(
            "charged_surface_response_options require a charged_surface_response")
    if charged_surface_response is not None and (
            precomputed_transport is not None or charging_poisson_system is not None):
        raise ValueError(
            "ordinary feature-step charged response is exclusive with precomputed or "
            "self-consistent charging transport")
    if ballistic_transport not in ("forward", "face_gather"):
        raise ValueError("ballistic_transport must be 'forward' or 'face_gather'")
    if reinitialization_method not in ("skfmm", "fsm", "cr2"):
        raise ValueError("reinitialization_method must be 'skfmm', 'fsm', or 'cr2'")
    if ballistic_transport == "face_gather" and (
            charging_poisson_system is not None or nodal_potential_v is not None):
        raise ValueError("deterministic ballistic face gather does not yet trace electric fields")
    if charged_surface_response is not None and ballistic_transport == "face_gather":
        raise ValueError(
            "charged surface response requires forward impact-position lineage; "
            "face_gather currently preserves direction but not impact position")

    verts, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    face_material = _face_material_ids(centroids, geometry)
    active_face = np.where(np.isin(face_material, etchable))[0]
    if active_face.size == 0:
        raise ValueError("current interface contains no requested etchable material")
    active_material = face_material[active_face]
    material_resolved_mechanism = hasattr(mechanism, "advance_by_material")
    if len(np.unique(active_material)) > 1 and not material_resolved_mechanism:
        raise ValueError(
            "multiple evolving materials require a material-resolved mechanism router")
    mesh_fingerprint = _surface_mesh_fingerprint(
        verts, faces, active_face, face_material, geometry)
    if surface_state is None:
        if surface_state_mesh_fingerprint is not None:
            raise ValueError("surface_state_mesh_fingerprint requires a supplied surface_state")
        if material_resolved_mechanism:
            if not hasattr(mechanism, "initial_state_by_material"):
                raise TypeError("material-resolved mechanism must initialize by material id")
            surface_state = mechanism.initial_state_by_material(active_material)
        else:
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
    scatter_options = (None if neutral_forward_scatter_options is None
                       else dict(neutral_forward_scatter_options))
    if neutral_forward_scatter is None and scatter_options is not None:
        raise ValueError("neutral_forward_scatter_options require a scatter model")
    if (neutral_forward_scatter is not None
            and not isinstance(neutral_forward_scatter, HwangGiapisSiO2ForwardScatter3D)):
        raise TypeError("neutral_forward_scatter must be HwangGiapisSiO2ForwardScatter3D")
    scatter_options = {} if scatter_options is None else scatter_options
    allowed_scatter_options = {
        "launch_offset", "periodic_lateral", "maximum_periodic_wraps"}
    unknown = set(scatter_options) - allowed_scatter_options
    if unknown:
        raise ValueError(
            "unknown neutral-forward-scatter options: " + ", ".join(sorted(unknown)))
    scatter_periodic = bool(scatter_options.get("periodic_lateral", False))
    periodic_neutral = bool(
        radiosity_options is not None and radiosity_options.get("periodic_lateral", False))
    response_options = None
    response_fixed_dt = None
    response_periodic = False
    if charged_surface_response is not None:
        response_options = ({} if charged_surface_response_options is None
                            else dict(charged_surface_response_options))
        allowed_response_options = {
            "launch_offset", "fixed_dt", "max_steps", "max_bounces",
            "relative_tail_tolerance", "adaptive_bounce_extension",
            "emergency_max_bounces", "trajectory_adaptive_horizon",
            "trajectory_emergency_max_steps", "periodic_lateral",
        }
        unknown = set(response_options) - allowed_response_options
        if unknown:
            raise ValueError(
                "unknown charged-surface response options: "
                + ", ".join(sorted(unknown)))
        response_fixed_dt = response_options.get("fixed_dt", trajectory_fixed_dt)
        if response_fixed_dt is None:
            raise ValueError(
                "charged surface response requires an explicit fixed_dt either in "
                "charged_surface_response_options or trajectory_fixed_dt")
        response_periodic = bool(response_options.get(
            "periodic_lateral",
            bool(field_periodic_lateral) if nodal_potential_v is not None else periodic_neutral))
        if periodic_neutral and not response_periodic:
            raise ValueError(
                "periodic neutral radiosity requires periodic response-enabled trajectories")
    charging_periodic = bool(
        charging_options is not None and charging_options.get("periodic_lateral", False))
    transport_periodic_lateral = bool(
        periodic_neutral or response_periodic or bool(field_periodic_lateral)
        or charging_periodic or scatter_periodic)
    if profile_periodic_lateral is None:
        profile_periodic_lateral = transport_periodic_lateral
    profile_periodic_lateral = bool(profile_periodic_lateral)
    if transport_periodic_lateral and not profile_periodic_lateral:
        raise ValueError(
            "periodic transport requires periodic lateral profile evolution")
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
            trajectory_adaptive_horizon=trajectory_adaptive_horizon,
            trajectory_emergency_max_steps=trajectory_emergency_max_steps,
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
                device=transport_device,
                adaptive_horizon=trajectory_adaptive_horizon,
                emergency_max_steps=trajectory_emergency_max_steps)
            transport = merge_boundary_transport_results_3d(
                charging.transport, uncharged_transport)
    elif nodal_potential_v is None:
        if initial_charge_node_c is not None or charging_options is not None:
            raise ValueError("charging state/options require charging_poisson_system")
        if (potential_origin is not None or potential_spacing is not None
                or (trajectory_fixed_dt is not None
                    and charged_surface_response is None)):
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
        elif charged_surface_response is not None:
            # A reflected/emitted flight must start from the exact primary impact position.
            # Field-free primaries are straight rays, so use the certified one-query hard-hit
            # tracer rather than approximating that ray with a zero-field time integrator.  The
            # latter can place an exact surface crossing on a step boundary and later report the
            # corresponding exit back-face.  Reflected/emitted flights still use the common field
            # cascade below because they launch from an arbitrary surface point.
            charged_species = tuple(
                species for species in boundary.species if species.charge_number != 0)
            if not charged_species:
                raise ValueError(
                    "charged surface response requires at least one charged boundary species")
            charged_boundary = PlasmaBoundaryState(
                charged_species, boundary.reference_plane_m, provenance=boundary.provenance)
            charged_role = {species.name: role[species.name] for species in charged_species}
            charged_first_hit_options = {}
            if response_periodic:
                charged_first_hit_options = {
                    "periodic_lateral": True,
                    "domain_size": (np.asarray(geometry.phi.shape) - 1) * geometry.dx,
                }
            transport = trace_boundary_state_first_hit_3d(
                charged_boundary, charged_role, verts, faces, areas,
                source_bounds=source_bounds, source_z=source_z,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                mesh_origin_m=geometry.mesh_origin_m, n_position=n_position, seed=seed,
                face_gas_normals=face_gas_normals,
                device=transport_device, **charged_first_hit_options)
            uncharged_species = tuple(
                species for species in boundary.species if species.charge_number == 0)
            if uncharged_species:
                uncharged_boundary = PlasmaBoundaryState(
                    uncharged_species, boundary.reference_plane_m,
                    provenance=boundary.provenance)
                uncharged_role = {
                    species.name: role[species.name] for species in uncharged_species}
                uncharged_transport = trace_boundary_state_first_hit_3d(
                    uncharged_boundary, uncharged_role, verts, faces, areas,
                    source_bounds=source_bounds, source_z=source_z,
                    mesh_length_unit_m=geometry.mesh_length_unit_m,
                    mesh_origin_m=geometry.mesh_origin_m,
                    n_position=n_position, seed=seed, device=transport_device,
                    **first_hit_options)
                transport = merge_boundary_transport_results_3d(
                    transport, uncharged_transport)
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
            face_gas_normals=face_gas_normals,
            adaptive_horizon=trajectory_adaptive_horizon,
            emergency_max_steps=trajectory_emergency_max_steps)
    charged_surface_cascade = None
    if charged_surface_response is not None:
        response_options = dict(response_options)
        response_fixed_dt = response_options.pop("fixed_dt", response_fixed_dt)
        response_potential = (
            np.zeros(geometry.phi.shape, dtype=float)
            if nodal_potential_v is None else np.asarray(nodal_potential_v, dtype=float))
        response_origin = (
            np.zeros(3, dtype=float)
            if nodal_potential_v is None else np.asarray(potential_origin, dtype=float))
        response_spacing = (
            float(geometry.dx) if nodal_potential_v is None else potential_spacing)
        charged_names = {
            species.name: int(species.charge_number)
            for species in boundary.species if species.charge_number != 0}
        response_context = ChargedSurfaceContext3D(
            np.asarray(areas, dtype=float) * geometry.mesh_length_unit_m ** 2,
            face_gas_normals, face_material, None)
        transport, charged_surface_cascade = (
            apply_charged_surface_response_to_transport_3d(
                transport, charged_names, charged_surface_response,
                response_context, verts, faces, areas,
                nodal_potential_v=response_potential,
                potential_origin=response_origin,
                potential_spacing=response_spacing,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                launch_offset=response_options.pop("launch_offset", 1e-5),
                fixed_dt=response_fixed_dt,
                max_steps=response_options.pop("max_steps", trajectory_max_steps),
                max_bounces=response_options.pop("max_bounces", 16),
                relative_tail_tolerance=response_options.pop(
                    "relative_tail_tolerance", 0.0),
                adaptive_bounce_extension=response_options.pop(
                    "adaptive_bounce_extension", False),
                emergency_max_bounces=response_options.pop(
                    "emergency_max_bounces", None),
                trajectory_adaptive_horizon=response_options.pop(
                    "trajectory_adaptive_horizon", trajectory_adaptive_horizon),
                trajectory_emergency_max_steps=response_options.pop(
                    "trajectory_emergency_max_steps", trajectory_emergency_max_steps),
                periodic_lateral=response_options.pop(
                    "periodic_lateral", response_periodic),
                device=transport_device))
    neutral_forward_scatter_result = None
    chemistry_role = dict(role)
    if neutral_forward_scatter is not None:
        scatter_context = ChargedSurfaceContext3D(
            np.asarray(areas, dtype=float) * geometry.mesh_length_unit_m ** 2,
            face_gas_normals, face_material, None)
        transport, neutral_forward_scatter_result = (
            apply_hwang_giapis_forward_scatter_to_transport_3d(
                transport, neutral_forward_scatter, scatter_context,
                verts, faces, areas,
                domain_minimum=np.zeros(3),
                domain_maximum=(np.asarray(geometry.phi.shape) - 1) * geometry.dx,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                launch_offset=float(scatter_options.get("launch_offset", 1e-5)),
                periodic_lateral=scatter_periodic,
                maximum_periodic_wraps=int(scatter_options.get(
                    "maximum_periodic_wraps", 10000))))
        chemistry_role[neutral_forward_scatter.neutral_species_name] = (
            "energetic_bombardment")
    base_transport = transport
    neutral_radiosity_diagnostics = MappingProxyType({})
    neutral_surface_iterations = 0
    neutral_surface_residual = None
    if neutral_surface_fixed_point_tolerance is not None:
        if material_resolved_mechanism:
            raise ValueError(
                "neutral/surface fixed point requires a directly inspectable mechanism result")
        working_state = surface_state
        for iteration in range(int(neutral_surface_fixed_point_max_iterations)):
            transport, neutral_radiosity_diagnostics = _apply_diffuse_neutral_transport(
                base_transport, geometry, verts, faces, centroids, areas, face_material,
                active_face, working_state, mechanism, role, radiosity_options,
                transport_device)
            active_flux = _select_surface_fluxes(
                transport.surface_fluxes, active_face, len(faces), chemistry_role)
            trial = mechanism.advance(working_state, active_flux, 0.0)
            change = getattr(trial, "transport_fixed_point_change", None)
            if change is None:
                raise TypeError(
                    "quasi-steady mechanism must report transport_fixed_point_change")
            neutral_surface_residual = float(np.max(np.abs(np.asarray(change, dtype=float))))
            neutral_surface_iterations = iteration + 1
            working_state = trial.state
            if neutral_surface_residual <= float(neutral_surface_fixed_point_tolerance):
                break
        else:
            raise RuntimeError(
                "neutral/surface fixed point did not converge: "
                f"residual={neutral_surface_residual:.6g}, "
                f"tolerance={float(neutral_surface_fixed_point_tolerance):.6g}, "
                f"iterations={int(neutral_surface_fixed_point_max_iterations)}")
        surface_state = working_state
        surface = mechanism.advance(surface_state, active_flux, float(duration_s))
    else:
        if radiosity_options is not None:
            transport, neutral_radiosity_diagnostics = _apply_diffuse_neutral_transport(
                base_transport, geometry, verts, faces, centroids, areas, face_material,
                active_face, surface_state, mechanism, role, radiosity_options,
                transport_device)
        active_flux = _select_surface_fluxes(
            transport.surface_fluxes, active_face, len(faces), chemistry_role)
        surface = (mechanism.advance_by_material(
            surface_state, active_flux, float(duration_s), active_material)
            if material_resolved_mechanism
            else mechanism.advance(surface_state, active_flux, float(duration_s)))

    product_populations = tuple(getattr(surface, "product_populations", ()))
    product_redeposition = None
    if surface_product_redeposition_options is not None:
        if duration_s <= 0.0:
            raise ValueError("surface-product redeposition requires a positive feature duration")
        if not product_populations:
            raise ValueError(
                "surface-product redeposition is enabled but the mechanism emits no populations")
        product_redeposition = _apply_surface_product_redeposition(
            product_populations, geometry, verts, faces, centroids, areas, face_material,
            active_face, duration_s, surface_product_redeposition_options, transport_device)

    surface_etch_velocity = np.asarray(surface.etch_velocity_m_s, dtype=float)
    surface_growth_velocity = np.asarray(
        getattr(surface, "normal_growth_velocity_m_s", 0.0), dtype=float)
    try:
        surface_etch_velocity = np.broadcast_to(
            surface_etch_velocity, (len(active_face),))
        surface_growth_velocity = np.broadcast_to(
            surface_growth_velocity, (len(active_face),))
    except ValueError as error:
        raise ValueError(
            "surface recession/growth velocity does not match the active-face mesh") from error
    if (np.any(~np.isfinite(surface_etch_velocity))
            or np.any(surface_etch_velocity < 0.0)
            or np.any(~np.isfinite(surface_growth_velocity))
            or np.any(surface_growth_velocity < 0.0)):
        raise ValueError("surface recession/growth velocities must be finite and nonnegative")
    face_velocity = np.zeros(len(faces))
    face_velocity[active_face] = (
        (surface_etch_velocity - surface_growth_velocity)
        / geometry.mesh_length_unit_m)
    if product_redeposition is not None:
        face_velocity -= (
            product_redeposition.normal_growth_velocity_m_s
            / geometry.mesh_length_unit_m)
    maximum_speed = float(np.max(np.abs(face_velocity))) if face_velocity.size else 0.0
    maximum_recession = max(
        float(np.max(face_velocity)) if face_velocity.size else 0.0, 0.0)
    maximum_growth = max(
        float(np.max(-face_velocity)) if face_velocity.size else 0.0, 0.0)
    displacement = maximum_speed * float(duration_s)
    substeps = max(1, int(np.ceil(displacement / (float(cfl_number) * geometry.dx))))
    phi = np.array(geometry.phi, copy=True)
    periodic_seam_projection = 0.0
    if profile_periodic_lateral:
        phi, correction = _project_periodic_lateral_endpoints(phi)
        periodic_seam_projection = max(periodic_seam_projection, correction)
    xs, ys, zs = geometry.coordinate_arrays
    extension_geometry = dict(phi=phi, dx=geometry.dx, xs=xs, ys=ys, zs=zs)
    # Extend only from the material surface that is actually evolving.  Including pinned mask
    # triangles with zero velocity lets them win the nearest-face query below a narrow opening and
    # numerically pins a physically bombarded floor after roughly one grid cell of motion.
    extended_velocity = extend_velocity_3d(
        face_velocity[active_face], centroids[active_face],
        extension_geometry, 4.0 * geometry.dx)
    if profile_periodic_lateral:
        extended_velocity, correction = _project_periodic_lateral_endpoints(
            extended_velocity)
        periodic_seam_projection = max(periodic_seam_projection, correction)
        if periodic_seam_projection > 0.25 * geometry.dx:
            raise RuntimeError(
                "periodic profile seam projection exceeds one quarter cell; input geometry or "
                "surface velocity is not a resolved periodic field")
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
    if duration_s == 0.0:
        # A zero-duration transport/chemistry audit is an exact geometry no-op.
        # Reconstructing the union from independently redistanced material fields can
        # otherwise change marching-cubes connectivity even though no material moved.
        # Preserve the authoritative combined level set and material ownership bitwise.
        material_levelsets = (
            None if geometry.material_levelsets is None else {
                material_id: np.array(levelset, copy=True)
                for material_id, levelset in geometry.material_levelsets.items()})
        phi = np.array(geometry.phi, copy=True)
    elif geometry.material_levelsets is None:
        for _ in range(substeps):
            phi = advect_3d(
                phi, extended_velocity, geometry.dx, float(duration_s) / substeps,
                periodic_axes=((0, 1) if profile_periodic_lateral else ()))
            phi[pinned] = geometry.phi[pinned]
            if profile_periodic_lateral:
                phi, correction = _project_periodic_lateral_endpoints(phi)
                periodic_seam_projection = max(periodic_seam_projection, correction)
    else:
        material_levelsets = _advect_exposed_material_levelsets(
            geometry.material_levelsets, etchable, extended_velocity,
            geometry.dx, duration_s, substeps,
            periodic_lateral=profile_periodic_lateral)
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
                    _redistance_feature_field(
                        levelset, geometry.dx, reinitialization_method,
                        periodic_lateral=profile_periodic_lateral)
                    if material_id in etchable else levelset)
                for material_id, levelset in material_levelsets.items()}
            if profile_periodic_lateral:
                material_levelsets = {
                    material_id: _project_periodic_lateral_endpoints(levelset)[0]
                    for material_id, levelset in material_levelsets.items()}
            phi = np.maximum.reduce(tuple(material_levelsets.values()))
        phi = _redistance_feature_field(
            phi, geometry.dx, reinitialization_method,
            periodic_lateral=profile_periodic_lateral)
        if material_levelsets is None:
            phi[pinned] = geometry.phi[pinned]
            if profile_periodic_lateral:
                phi, correction = _project_periodic_lateral_endpoints(phi)
                periodic_seam_projection = max(periodic_seam_projection, correction)
    reinitialized_centerline = phi[center]
    reinitialized_crossing = np.flatnonzero(
        (reinitialized_centerline[:-1] >= 0.0) & (reinitialized_centerline[1:] < 0.0))
    if reinitialized_crossing.size == 1:
        lower = int(reinitialized_crossing[0])
        center_diagnostics["centerline_reinitialized_interface_fraction"] = float(
            reinitialized_centerline[lower]
            / (reinitialized_centerline[lower] - reinitialized_centerline[lower + 1]))
    output_material_id = np.array(geometry.material_id, copy=True)
    if material_levelsets is not None:
        material_ids = np.asarray(sorted(material_levelsets), dtype=int)
        material_stack = np.stack([
            material_levelsets[int(material_id)] for material_id in material_ids])
        owner = material_ids[np.argmax(material_stack, axis=0)]
        output_material_id = np.where(phi >= 0.0, owner, 0)
    if duration_s == 0.0:
        removed_unresolved_solid_cells = 0
        unresolved_solid_mask = np.zeros_like(output_material_id, dtype=bool)
    else:
        phi, removed_unresolved_solid_cells, unresolved_solid_mask = (
            _remove_unresolved_subcell_solid_components(
            phi, output_material_id, etchable, geometry.dx)
        )
    if removed_unresolved_solid_cells:
        if material_levelsets is not None:
            material_levelsets, phi, output_material_id = (
                _apply_subcell_cleanup_to_material_levelsets(
                    material_levelsets, unresolved_solid_mask, output_material_id,
                    etchable, geometry.dx, reinitialization_method,
                    profile_periodic_lateral))
            _, remaining_unresolved_cells, _ = (
                _remove_unresolved_subcell_solid_components(
                    phi, output_material_id, etchable, geometry.dx))
            if remaining_unresolved_cells:
                raise RuntimeError(
                    "material-layer topology update left an unresolved subcell component")
        else:
            phi = _redistance_feature_field(
                phi, geometry.dx, reinitialization_method,
                periodic_lateral=profile_periodic_lateral)
            phi[pinned] = geometry.phi[pinned]
    if profile_periodic_lateral:
        phi, correction = _project_periodic_lateral_endpoints(phi)
        periodic_seam_projection = max(periodic_seam_projection, correction)

    output_geometry = FeatureGeometry3D(
        phi, output_material_id, geometry.dx, geometry.mesh_length_unit_m,
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
    if product_redeposition is not None:
        exchange_limitations = tuple(
            item for item in exchange_limitations
            if item != "outgoing physical-sputter material is not redeposited unless product "
            "transport is enabled") + (
                "redeposition v1 permits same-material growth only; cross-material films are refused",
            )
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(dict.fromkeys(transport_limitations)) + (
            "first-order material-local conservative surface-state remap with declared intensive-field exceptions",
            "physical volume-topology-changing surface steps are refused",
            "first-order Godunov interface advection",
        ) + tuple(surface.validity.known_model_form_omissions) + exchange_limitations,
        parameter_evidence_supports_prediction=(
            surface.validity.parameter_evidence_supports_prediction
            and neutral_forward_scatter_result is None),
        nonpredictive_parameters=(
            surface.validity.nonpredictive_parameters + ((
                "neutral_forward_scatter.critical_angle_deg",
                "neutral_forward_scatter.gas_to_effective_surface_mass_ratio",
            ) if neutral_forward_scatter_result is not None else ())))
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
        charged_surface_cascade=charged_surface_cascade,
        neutral_forward_scatter=neutral_forward_scatter_result,
        surface_product_redeposition=product_redeposition,
        diagnostics=dict(
            face_count=int(len(faces)), active_face_count=int(active_face.size),
            max_velocity_m_s=maximum_speed * geometry.mesh_length_unit_m,
            max_recession_velocity_m_s=maximum_recession * geometry.mesh_length_unit_m,
            max_growth_velocity_m_s=maximum_growth * geometry.mesh_length_unit_m,
            max_surface_mechanism_growth_velocity_m_s=(
                float(np.max(surface_growth_velocity))
                if surface_growth_velocity.size else 0.0),
            max_displacement_mesh_units=displacement, cfl_substeps=int(substeps),
            cfl_number=float(cfl_number), reinitialized=bool(reinitialize),
            reinitialization_method=(reinitialization_method if reinitialize else None),
            profile_periodic_lateral=profile_periodic_lateral,
            periodic_seam_projection_max_mesh_units=periodic_seam_projection,
            removed_unresolved_solid_cells=removed_unresolved_solid_cells,
            self_consistent_charging=charging is not None,
            charging_iterations=(0 if charging is None else len(charging.history)),
            charging_converged=(None if charging is None else charging.converged),
            charged_surface_response_applied=charged_surface_cascade is not None,
            charged_surface_response_field=(
                None if charged_surface_cascade is None
                else ("supplied_nodal_potential" if nodal_potential_v is not None
                      else "explicit_zero_field")),
            charged_surface_response_bounces=(
                0 if charged_surface_cascade is None
                else len(charged_surface_cascade.transfers)),
            charged_surface_response_reimpact_events=(
                0 if charged_surface_cascade is None else sum(
                    flight.incident.event_face.size
                    for bounce in charged_surface_cascade.flights_by_bounce
                    for flight in bounce)),
            charged_surface_response_relative_charge_error=(
                None if charged_surface_cascade is None
                else charged_surface_cascade.relative_charge_balance_error),
            charged_surface_response_maximum_energy_error=(
                None if charged_surface_cascade is None else max(
                    transfer.relative_kinetic_energy_balance_error
                    for transfer in charged_surface_cascade.transfers)),
            charged_surface_response_tail_l1_error_bound=(
                None if charged_surface_cascade is None else
                charged_surface_cascade.tail_closure_l1_current_error_bound_relative),
            charged_surface_response_bounce_budget_extensions=(
                0 if charged_surface_cascade is None else
                charged_surface_cascade.bounce_budget_extension_count),
            neutral_forward_scatter_applied=(
                neutral_forward_scatter_result is not None),
            neutral_forward_scatter_rate_s=(
                0.0 if neutral_forward_scatter_result is None else
                neutral_forward_scatter_result.scattered_rate_s),
            neutral_forward_scatter_landed_rate_s=(
                0.0 if neutral_forward_scatter_result is None else
                neutral_forward_scatter_result.flight.landed_rate_s),
            neutral_forward_scatter_escaped_rate_s=(
                0.0 if neutral_forward_scatter_result is None else
                neutral_forward_scatter_result.flight.escaped_rate_s),
            neutral_forward_scatter_particle_balance_error=(
                None if neutral_forward_scatter_result is None else max(
                    neutral_forward_scatter_result.relative_surface_particle_balance_error,
                    neutral_forward_scatter_result.flight.relative_particle_balance_error)),
            neutral_forward_scatter_energy_balance_error=(
                None if neutral_forward_scatter_result is None else
                neutral_forward_scatter_result.relative_surface_energy_balance_error),
            product_routing_complete=product_routing_complete,
            product_population_count=len(product_populations),
            product_transport_ready=product_transport_ready,
            product_redeposition_enabled=product_redeposition is not None,
            product_redeposition_relative_balance_error=(
                None if product_redeposition is None
                else product_redeposition.maximum_relative_balance_error),
            neutral_radiosity=neutral_radiosity_diagnostics,
            neutral_surface_fixed_point_iterations=neutral_surface_iterations,
            neutral_surface_fixed_point_residual=neutral_surface_residual,
            neutral_surface_fixed_point_tolerance=(
                None if neutral_surface_fixed_point_tolerance is None
                else float(neutral_surface_fixed_point_tolerance)),
            **center_diagnostics),
        validity=validity)


def solve_feature_3d(
        geometry: FeatureGeometry3D, boundary: PlasmaBoundaryState,
        species_role: Mapping[str, str], mechanism, *,
        etchable_material_ids, duration_s, n_steps, source_bounds, source_z,
        n_position=256, seed=0, cfl_number=0.3, reinitialize=True,
        transport_device=None, nodal_potential_v=None, potential_origin=None,
        potential_spacing=None, trajectory_fixed_dt=None, trajectory_max_steps=10000,
        field_periodic_lateral=False, profile_periodic_lateral=None,
        charging_poisson_system: NodalPoissonSystem3D | None = None,
        charging_system_builder=None, initial_charge_node_c=None, charging_options=None,
        charged_surface_response=None, charged_surface_response_options=None,
        neutral_forward_scatter=None, neutral_forward_scatter_options=None,
        neutral_radiosity_options=None,
        neutral_surface_fixed_point_tolerance=None,
        neutral_surface_fixed_point_max_iterations=20,
        surface_product_redeposition_options=None,
        ballistic_transport="forward",
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
    if charged_surface_response is not None and (
            charging_poisson_system is not None or charging_system_builder is not None):
        raise ValueError(
            "ordinary feature response cannot be combined with a self-consistent charging solve; "
            "use the charging co-evolution response path")
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
                profile_periodic_lateral=profile_periodic_lateral,
                charging_poisson_system=step_poisson_system,
                initial_charge_node_c=step_initial_charge,
                charging_options=charging_options,
                charged_surface_response=charged_surface_response,
                charged_surface_response_options=charged_surface_response_options,
                neutral_forward_scatter=neutral_forward_scatter,
                neutral_forward_scatter_options=neutral_forward_scatter_options,
                neutral_radiosity_options=neutral_radiosity_options,
                neutral_surface_fixed_point_tolerance=(
                    neutral_surface_fixed_point_tolerance),
                neutral_surface_fixed_point_max_iterations=(
                    neutral_surface_fixed_point_max_iterations),
                surface_product_redeposition_options=surface_product_redeposition_options,
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
