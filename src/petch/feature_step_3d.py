"""Dimensional feature evolution through the new physical contracts.

Each step transfers surface state material-by-material with an area-conservative, bounded remap. Smooth
CFL-limited motion can be iterated. By default every topology change is refused; an explicit policy may
continue only periodic gas-cavity enclosure/opening while solid/material component counts and domain
breakthrough remain unchanged. Material appearance/disappearance, excessive remap distance, impossible
coverage compression, and every other topology event remain refusals. This makes the multi-step loop
explicit about the domain in which surface history is numerically supported.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from itertools import product
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.ndimage import label
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
from .neutral_radiosity_3d import (
    DiffuseNeutralNoSinkError,
    solve_diffuse_neutral_radiosity_3d,
)
from .extruded_exchange_3d import build_extruded_triangle_exchange_3d
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
from .feature_geometry_state_3d import (
    FeatureGeometry3D,
    face_material_ids_3d as _face_material_ids,
)
from .feature_geometry_backend_3d import UniformFeatureGeometryBackend3D
from .surface_mesh_3d import TriangleSurface3D
from .surface_transfer_3d import build_surface_transfer_3d
from .surface_partitioned_overlap_3d import (
    build_partitioned_surface_overlap_transfer_3d,
)
from .surface_common_refinement_3d import (
    build_surface_common_refinement_transfer_3d,
)
from .threed import advect_3d, extend_velocity_3d, reinit_cr2, reinit_fsm, reinit_narrow


class SurfaceTopologyChangeError(ValueError):
    """Structured refusal when conservative surface-state remap needs a topology event.

    The exception remains a :class:`ValueError` for backwards compatibility, while exposing the
    old/new physical-volume signatures so a campaign driver can distinguish a resolved feature
    closure from computational-domain breakthrough.  It never authorizes an implicit remap.
    """

    def __init__(
            self, message, *, method, old_topology, new_topology,
            old_mesh_topology, new_mesh_topology, changed_slice_topology):
        super().__init__(str(message))
        self.method = str(method)
        self.old_topology = tuple(old_topology)
        self.new_topology = tuple(new_topology)
        self.old_mesh_topology = tuple(old_mesh_topology)
        self.new_mesh_topology = tuple(new_mesh_topology)
        self.changed_slice_topology = MappingProxyType(
            dict(changed_slice_topology))

    @property
    def event_kind(self):
        """Return a geometry-only classification without assigning process semantics."""
        old_solid, old_cavity, old_breakthrough, old_material = self.old_topology
        new_solid, new_cavity, new_breakthrough, new_material = self.new_topology
        if new_solid != old_solid:
            return "solid_component_change"
        if new_material != old_material:
            return "material_component_change"
        if bool(new_breakthrough) != bool(old_breakthrough):
            return (
                "domain_gas_breakthrough" if new_breakthrough
                else "domain_gas_disconnect")
        if new_cavity > old_cavity:
            return "gas_cavity_enclosed"
        if new_cavity < old_cavity:
            return "gas_cavity_opened"
        return "other_topology_change"


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


def _surface_gas_normals(verts, faces, centroids, geometry):
    """Orient geometric face normals toward gas using the local level-set gradient.

    ``geometry.phi`` is positive in solid, so ``-grad(phi)`` is the gas direction at a
    regular interface.  The orientation decision must be local.  A finite probe a sizeable
    fraction of one cell away can cross a second interface in a thin sheet or narrow pocket;
    both probes then sample the same phase and can reverse an otherwise valid marching-cubes
    face.  Evaluating the exact derivative of the trilinear nodal interpolant at the triangle
    centroid avoids that nonlocal ambiguity while retaining the geometric triangle normal used
    by the hard-visibility operator.
    """
    triangle = np.asarray(verts)[np.asarray(faces)]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)

    def trilinear_gradient(field, point):
        coordinate = np.asarray(point, dtype=float) / geometry.dx
        lower = np.floor(coordinate).astype(int)
        for axis in range(3):
            lower[:, axis] = np.clip(lower[:, axis], 0, field.shape[axis] - 2)
        fraction = np.clip(coordinate - lower, 0.0, 1.0)
        gradient = np.zeros((len(point), 3))
        for axis in range(3):
            transverse = tuple(item for item in range(3) if item != axis)
            for first in (0, 1):
                first_weight = (
                    fraction[:, transverse[0]]
                    if first else 1.0 - fraction[:, transverse[0]])
                for second in (0, 1):
                    second_weight = (
                        fraction[:, transverse[1]]
                        if second else 1.0 - fraction[:, transverse[1]])
                    negative = lower.copy()
                    positive = lower.copy()
                    positive[:, axis] += 1
                    negative[:, transverse[0]] += first
                    positive[:, transverse[0]] += first
                    negative[:, transverse[1]] += second
                    positive[:, transverse[1]] += second
                    gradient[:, axis] += (
                        first_weight * second_weight
                        * (field[tuple(positive.T)] - field[tuple(negative.T)])
                        / geometry.dx)
        return gradient

    centroid = np.asarray(centroids, dtype=float)
    gradient = trilinear_gradient(geometry.phi, centroid)
    directional_derivative = np.einsum("ij,ij->i", normal, gradient)
    gradient_scale = np.linalg.norm(gradient, axis=1)
    ambiguous = np.abs(directional_derivative) <= (
        64.0 * np.finfo(float).eps * np.maximum(gradient_scale, 1.0))
    if np.any(ambiguous):
        raise RuntimeError(
            "surface normal is undefined at a level-set critical point for "
            f"{int(np.count_nonzero(ambiguous))} of {len(normal)} faces")
    # phi is positive in solid.  A positive derivative means ``normal`` points into solid.
    flip = directional_derivative > 0.0
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


def _extract_uniform_surface_arrays(geometry):
    """Bridge the read-only backend surface to the legacy writable array contract.

    This is deliberately a behavior-neutral seam.  The existing mesh/state fingerprint remains
    authoritative; backend fingerprints are not substituted into checkpoint or remap semantics.
    """
    surface = UniformFeatureGeometryBackend3D(geometry).extract_surface()
    return tuple(np.array(value, copy=True, order="C") for value in (
        surface.vertices_mesh,
        surface.faces,
        surface.centroids_mesh,
        surface.areas_mesh2,
        surface.face_material_id,
    ))


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


def _periodic_component_roots(field):
    """Return 6-connected component roots after wrapping x/y neighbor pairs."""
    occupied = np.asarray(field, dtype=bool)
    component, count = label(occupied)
    parent = np.arange(int(count) + 1)

    def find(index):
        index = int(index)
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(left, right):
        left = find(left); right = find(right)
        if left and right and left != right:
            parent[right] = left

    for axis in (0, 1):
        first = np.take(component, 0, axis=axis)
        last = np.take(component, -1, axis=axis)
        selected = (first > 0) & (last > 0)
        for left, right in zip(first[selected], last[selected]):
            union(left, right)
    roots = np.zeros(int(count) + 1, dtype=int)
    for index in range(1, int(count) + 1):
        roots[index] = find(index)
    rooted = roots[component]
    return rooted, {int(value) for value in np.unique(rooted) if value > 0}


def _periodic_component_sizes(field):
    """Return periodic 6-connected component sizes, largest first.

    This is diagnostic rather than a topology decision: unlike a bounded-label
    histogram it merges components touching opposite x/y faces before counting.
    Keeping the sizes in refusals distinguishes a resolved material separation
    from a one-cell material-ownership flicker at an interface.
    """
    rooted, roots = _periodic_component_roots(field)
    return tuple(sorted(
        (int(np.count_nonzero(rooted == root)) for root in roots),
        reverse=True,
    ))


def _component_roots_with_resolved_volume_cells(rooted_labels, *, periodic_lateral):
    """Return component roots that own at least one physical hexahedral cell.

    Eight arbitrary nodes are not sufficient support for a 3-D volume cell.  In
    particular, an extruded ``1 x N x 1`` filament can contain many nodes while
    remaining zero cells thick in two directions.  A component is resolved only
    when the eight corners of at least one grid cell all carry the same root.

    Periodic x/y fields contain only their unique nodal core here.  Their last
    physical cells therefore wrap to node zero; z remains bounded.
    """
    rooted = np.asarray(rooted_labels, dtype=int)
    if rooted.ndim != 3:
        raise ValueError("volume-cell support requires a 3-D component-label field")
    nx, ny, nz = rooted.shape
    if nz < 2 or (not periodic_lateral and (nx < 2 or ny < 2)):
        return set()

    if periodic_lateral:
        x_next = np.roll(rooted, -1, axis=0)
        y_next = np.roll(rooted, -1, axis=1)
        xy_next = np.roll(x_next, -1, axis=1)
        corners = (
            rooted[:, :, :-1], x_next[:, :, :-1],
            y_next[:, :, :-1], xy_next[:, :, :-1],
            rooted[:, :, 1:], x_next[:, :, 1:],
            y_next[:, :, 1:], xy_next[:, :, 1:],
        )
    else:
        corners = (
            rooted[:-1, :-1, :-1], rooted[1:, :-1, :-1],
            rooted[:-1, 1:, :-1], rooted[1:, 1:, :-1],
            rooted[:-1, :-1, 1:], rooted[1:, :-1, 1:],
            rooted[:-1, 1:, 1:], rooted[1:, 1:, 1:],
        )
    owner = corners[0]
    resolved = owner > 0
    for corner in corners[1:]:
        resolved &= corner == owner
    return {int(value) for value in np.unique(owner[resolved]) if value > 0}


def _periodic_material_component_sizes(geometry, etchable_material_ids):
    materials = tuple(sorted({int(value) for value in etchable_material_ids}))
    solid = ((geometry.phi > 0.0)
             & np.isin(geometry.material_id, materials))[:-1, :-1, :]
    core_material = np.asarray(geometry.material_id)[:-1, :-1, :]
    return tuple(
        (material, _periodic_component_sizes(
            solid & (core_material == material)))
        for material in materials
    )


def _periodic_physical_volume_topology_signature(
        geometry, etchable_material_ids):
    """Topology gate for duplicate-endpoint x/y-periodic feature cells.

    A bounded-grid Euler number treats the two stored copies of each periodic
    endpoint as distinct boundaries. As an etched groove becomes resolved, that
    representation can manufacture a handle even though the physical periodic
    surface deforms smoothly. The periodic gate instead audits the events that
    invalidate material-local remapping: solid component creation/destruction,
    enclosed gas-cavity creation/destruction, gas breakthrough between the two
    open z boundaries, and per-material component changes.
    """
    materials = tuple(sorted({int(value) for value in etchable_material_ids}))
    if not materials or any(value <= 0 for value in materials):
        raise ValueError("periodic topology requires positive material ids")
    # Drop the duplicate x/y endpoint planes. Adjacency across each remaining
    # first/last pair is then added explicitly by ``_periodic_component_roots``.
    solid = ((geometry.phi > 0.0)
             & np.isin(geometry.material_id, materials))[:-1, :-1, :]
    solid_labels, solid_roots = _periodic_component_roots(solid)
    gas_labels, gas_roots = _periodic_component_roots(~solid)
    open_boundary_roots = {
        int(value) for value in np.unique(
            np.concatenate((gas_labels[:, :, 0].ravel(),
                            gas_labels[:, :, -1].ravel())))
        if value > 0}
    enclosed_gas = gas_roots - open_boundary_roots
    lower = {int(value) for value in np.unique(gas_labels[:, :, 0]) if value > 0}
    upper = {int(value) for value in np.unique(gas_labels[:, :, -1]) if value > 0}
    material_components = []
    core_material = np.asarray(geometry.material_id)[:-1, :-1, :]
    for material in materials:
        _, roots = _periodic_component_roots(solid & (core_material == material))
        material_components.append((material, len(roots)))
    return (
        int(len(solid_roots)),
        int(len(enclosed_gas)),
        bool(lower & upper),
        tuple(material_components),
    )


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
        phi, material_id, etchable_material_ids, dx, *, periodic_lateral=False):
    updated = np.array(phi, copy=True)
    solid = (updated > 0.0) & np.isin(material_id, tuple(etchable_material_ids))
    if periodic_lateral:
        core_solid = solid[:-1, :-1, :]
        rooted, roots = _periodic_component_roots(core_solid)
        resolved_roots = _component_roots_with_resolved_volume_cells(
            rooted, periodic_lateral=True)
        unresolved_roots = tuple(roots - resolved_roots)
        core_mask = np.isin(rooted, unresolved_roots)
        unresolved = np.zeros(updated.shape, dtype=bool)
        unresolved[:-1, :-1, :] = core_mask
        unresolved[-1, :-1, :] = core_mask[0, :, :]
        unresolved[:-1, -1, :] = core_mask[:, 0, :]
        unresolved[-1, -1, :] = core_mask[0, 0, :]
        if not np.any(core_mask):
            return updated, 0, unresolved
        updated[unresolved] = -np.maximum(
            np.abs(updated[unresolved]), float(dx))
        return updated, int(np.count_nonzero(core_mask)), unresolved
    component, count = label(solid)
    if count == 0:
        return updated, 0, np.zeros(updated.shape, dtype=bool)
    resolved_labels = _component_roots_with_resolved_volume_cells(
        component, periodic_lateral=False)
    unresolved_label = tuple(
        index for index in range(1, int(count) + 1)
        if index not in resolved_labels)
    if not unresolved_label:
        return updated, 0, np.zeros(updated.shape, dtype=bool)
    unresolved = np.isin(component, unresolved_label)
    # A subcell component has no resolved 3-D volume. Give it an unambiguous gas sign,
    # then let the signed-distance reconstruction restore a consistent neighborhood.
    updated[unresolved] = -np.maximum(np.abs(updated[unresolved]), float(dx))
    return updated, int(np.count_nonzero(unresolved)), unresolved


def _new_unresolved_subcell_material_component_mask(
        phi, material_id, previous_material_id, etchable_material_ids, *,
        periodic_lateral):
    """Locate newly born material islands too small to own one volume cell.

    A material component supported by fewer than eight corner nodes has no
    resolved hexahedral volume.  It is eligible for cleanup only when *every*
    node changed owner in the candidate step.  A component that existed before
    the step and became disconnected is therefore a real topology event and is
    still refused.
    """
    field = np.asarray(phi, dtype=float)
    owner = np.asarray(material_id)
    previous = np.asarray(previous_material_id)
    if (field.ndim != 3 or owner.shape != field.shape
            or previous.shape != field.shape
            or np.any(~np.isfinite(field))):
        raise ValueError("material-island cleanup requires matching finite 3-D fields")
    materials = tuple(sorted({int(value) for value in etchable_material_ids}))
    if not materials or any(value <= 0 for value in materials):
        raise ValueError("material-island cleanup requires positive material ids")

    core_slice = (slice(None, -1), slice(None, -1), slice(None))
    core_solid = field[core_slice] > 0.0 if periodic_lateral else field > 0.0
    core_owner = owner[core_slice] if periodic_lateral else owner
    core_previous = previous[core_slice] if periodic_lateral else previous
    core_mask = np.zeros(core_solid.shape, dtype=bool)
    for material in materials:
        occupied = core_solid & (core_owner == material)
        if periodic_lateral:
            rooted, roots = _periodic_component_roots(occupied)
            resolved_roots = _component_roots_with_resolved_volume_cells(
                rooted, periodic_lateral=True)
            for root in roots:
                selected = rooted == root
                if (root not in resolved_roots
                        and np.all(core_previous[selected] != material)):
                    core_mask |= selected
        else:
            component, count = label(occupied)
            resolved_labels = _component_roots_with_resolved_volume_cells(
                component, periodic_lateral=False)
            for index in range(1, int(count) + 1):
                selected = component == index
                if (index not in resolved_labels
                        and np.all(core_previous[selected] != material)):
                    core_mask |= selected

    if not periodic_lateral:
        return core_mask, int(np.count_nonzero(core_mask))
    mask = np.zeros(field.shape, dtype=bool)
    mask[:-1, :-1, :] = core_mask
    mask[-1, :-1, :] = core_mask[0, :, :]
    mask[:-1, -1, :] = core_mask[:, 0, :]
    mask[-1, -1, :] = core_mask[0, 0, :]
    return mask, int(np.count_nonzero(core_mask))


def _restore_unresolved_material_ownership(
        material_levelsets, repair_mask, candidate_material_id,
        previous_material_id, dx, reinitialization_method, periodic_lateral):
    """Suppress a newly born subcell material island using prior ownership.

    This changes no resolved material volume: the selected component has fewer
    than the eight nodes required to bound one hexahedral cell.  Previous gas
    ownership retires the spurious layer sign; previous solid ownership restores
    that layer before all affected distance fields are reconstructed.
    """
    repair = np.asarray(repair_mask, dtype=bool)
    candidate = np.asarray(candidate_material_id)
    previous = np.asarray(previous_material_id)
    updated = {int(material_id): np.asarray(levelset, dtype=float).copy()
               for material_id, levelset in material_levelsets.items()}
    if (repair.shape != candidate.shape or previous.shape != candidate.shape
            or not np.any(repair)
            or any(levelset.shape != repair.shape for levelset in updated.values())):
        raise ValueError("material ownership repair requires matching nonempty fields")
    affected = set()
    for material in sorted(int(value) for value in np.unique(candidate[repair])):
        if material <= 0 or material not in updated:
            raise RuntimeError("subcell material island has no authoritative level set")
        selected = repair & (candidate == material)
        layer = updated[material]
        layer[selected] = -np.maximum(np.abs(layer[selected]), float(dx))
        affected.add(material)
        for prior in sorted(int(value) for value in np.unique(previous[selected])):
            if prior <= 0:
                continue
            if prior not in updated:
                raise RuntimeError("prior material owner has no authoritative level set")
            restore = selected & (previous == prior)
            prior_layer = updated[prior]
            prior_layer[restore] = np.maximum(
                np.abs(prior_layer[restore]), float(dx))
            affected.add(prior)
    for material in sorted(affected):
        updated[material] = _redistance_feature_field(
            updated[material], dx, reinitialization_method,
            periodic_lateral=periodic_lateral)
    if periodic_lateral:
        updated = {
            material: _project_periodic_lateral_endpoints(levelset)[0]
            for material, levelset in updated.items()}
    combined = np.maximum.reduce(tuple(updated.values()))
    combined = _redistance_feature_field(
        combined, dx, reinitialization_method,
        periodic_lateral=periodic_lateral)
    material_ids = np.asarray(sorted(updated), dtype=int)
    stack = np.stack([updated[int(material)] for material in material_ids])
    owner = material_ids[np.argmax(stack, axis=0)]
    owner = np.where(combined >= 0.0, owner, 0)
    if np.any(repair & (owner == candidate)):
        raise RuntimeError("subcell material ownership repair did not retire the island")
    return updated, combined, owner


def _unresolved_subcell_gas_cavity_mask(phi, *, periodic_lateral):
    """Find enclosed gas components too small to occupy one resolved volume cell.

    A gas component is resolved only if it owns all eight corners of at least one
    hexahedral cell.  Counting nodes is insufficient: a periodic, translationally
    extruded ``1 x N x 1`` component may contain eight or more nodes but enclose no
    cell volume.  Such pockets arise when redistancing crosses zero by a tiny
    amount near a depositing/re-entrant surface.  The returned mask includes
    duplicate periodic endpoint nodes, while the count is measured on the unique
    periodic core.
    """
    field = np.asarray(phi, dtype=float)
    if field.ndim != 3 or np.any(~np.isfinite(field)):
        raise ValueError("gas-cavity cleanup requires one finite 3-D level-set field")
    gas = field <= 0.0
    if periodic_lateral:
        core_gas = gas[:-1, :-1, :]
        rooted, roots = _periodic_component_roots(core_gas)
        open_roots = {
            int(value) for value in np.unique(np.concatenate((
                rooted[:, :, 0].ravel(), rooted[:, :, -1].ravel())))
            if value > 0}
        resolved_roots = _component_roots_with_resolved_volume_cells(
            rooted, periodic_lateral=True)
        unresolved_roots = {
            root for root in roots - open_roots
            if root not in resolved_roots}
        core_mask = np.isin(rooted, tuple(unresolved_roots))
        mask = np.zeros(field.shape, dtype=bool)
        mask[:-1, :-1, :] = core_mask
        mask[-1, :-1, :] = core_mask[0, :, :]
        mask[:-1, -1, :] = core_mask[:, 0, :]
        mask[-1, -1, :] = core_mask[0, 0, :]
        return mask, int(np.count_nonzero(core_mask))

    component, count = label(gas)
    if count == 0:
        return np.zeros(field.shape, dtype=bool), 0
    boundary = np.concatenate((
        component[0, :, :].ravel(), component[-1, :, :].ravel(),
        component[:, 0, :].ravel(), component[:, -1, :].ravel(),
        component[:, :, 0].ravel(), component[:, :, -1].ravel(),
    ))
    open_component = {int(value) for value in np.unique(boundary) if value > 0}
    resolved_components = _component_roots_with_resolved_volume_cells(
        component, periodic_lateral=False)
    unresolved = tuple(
        index for index in range(1, int(count) + 1)
        if index not in open_component and index not in resolved_components)
    mask = np.isin(component, unresolved)
    return mask, int(np.count_nonzero(mask))


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


def _apply_subcell_gas_fill_to_material_levelsets(
        material_levelsets, fill_mask, etchable_material_ids, dx,
        reinitialization_method, periodic_lateral):
    """Heal an unresolved gas bubble in its nearest authoritative material layer."""
    fill = np.asarray(fill_mask, dtype=bool)
    updated = {int(material_id): np.asarray(levelset, dtype=float).copy()
               for material_id, levelset in material_levelsets.items()}
    if not updated or any(levelset.shape != fill.shape for levelset in updated.values()):
        raise ValueError("gas fill requires matching authoritative material level sets")
    material_ids = np.asarray(sorted(updated), dtype=int)
    material_stack = np.stack([updated[int(material_id)] for material_id in material_ids])
    owner = material_ids[np.argmax(material_stack, axis=0)]
    etchable = {int(value) for value in etchable_material_ids}
    selected_owner = set(int(value) for value in np.unique(owner[fill]))
    if not selected_owner or not selected_owner.issubset(etchable):
        raise RuntimeError("unresolved gas cavity is not bounded by an evolving material")
    accounted = np.zeros(fill.shape, dtype=bool)
    for material_id in sorted(selected_owner):
        selected = fill & (owner == material_id)
        levelset = updated[material_id]
        levelset[selected] = np.maximum(np.abs(levelset[selected]), float(dx))
        updated[material_id] = _redistance_feature_field(
            levelset, dx, reinitialization_method,
            periodic_lateral=periodic_lateral)
        accounted |= selected
    if not np.array_equal(accounted, fill):
        raise RuntimeError("subcell gas fill did not assign every selected node")
    if periodic_lateral:
        updated = {
            material_id: _project_periodic_lateral_endpoints(levelset)[0]
            for material_id, levelset in updated.items()}
    combined = np.maximum.reduce(tuple(updated.values()))
    combined = _redistance_feature_field(
        combined, dx, reinitialization_method,
        periodic_lateral=periodic_lateral)
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


def _periodic_lateral_surface_images(values, centroids, domain_size):
    """Wrap face samples into the eight neighboring lateral periodic cells.

    ``extend_velocity_3d`` is an ordinary nearest-surface extension. Without wrapped
    images, nodes near a periodic seam see only the faces on their stored side of the
    duplicate endpoint and can acquire a nonperiodic velocity even when the physical
    boundary condition is periodic.
    """
    values = np.asarray(values, dtype=float)
    centroids = np.asarray(centroids, dtype=float)
    domain = np.asarray(domain_size, dtype=float)
    if (values.ndim != 1 or centroids.shape != (len(values), 3)
            or domain.shape != (3,) or np.any(~np.isfinite(values))
            or np.any(~np.isfinite(centroids)) or np.any(~np.isfinite(domain))
            or np.any(domain[:2] <= 0.0)):
        raise ValueError("invalid periodic surface-image inputs")
    shifts = np.asarray([
        (ix * domain[0], iy * domain[1], 0.0)
        for ix in (-1, 0, 1) for iy in (-1, 0, 1)], dtype=float)
    return (
        np.tile(values, len(shifts)),
        np.concatenate([centroids + shift for shift in shifts], axis=0),
    )


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


def _point_to_triangles_distance(point, triangles):
    """Return exact Euclidean distance from one point to each triangle.

    The plane projection is used when it lies inside the triangle; otherwise the closest point is
    on one of the three closed edges.  This is deliberately independent of triangle centroids:
    marching-cubes can retriangulate an almost stationary interface and move its centroids by
    multiple cells even though the represented material surface has barely moved.
    """
    point = np.asarray(point, dtype=float)
    triangle = np.asarray(triangles, dtype=float)
    if (point.shape != (3,) or triangle.ndim != 3 or triangle.shape[1:] != (3, 3)
            or np.any(~np.isfinite(point)) or np.any(~np.isfinite(triangle))):
        raise ValueError("point-to-triangle distance requires finite 3-D geometry")
    a = triangle[:, 0]
    b = triangle[:, 1]
    c = triangle[:, 2]
    ab = b - a
    ac = c - a
    normal = np.cross(ab, ac)
    normal_squared = np.einsum("ij,ij->i", normal, normal)
    coordinate_scale = max(
        float(np.max(np.abs(triangle))), float(np.max(np.abs(point))), 1.0)
    degeneracy_floor = (64.0 * np.finfo(float).eps * coordinate_scale) ** 4
    if np.any(normal_squared <= degeneracy_floor):
        raise ValueError("surface remap requires nondegenerate old triangles")

    ap = point[None, :] - a
    plane_parameter = np.einsum("ij,ij->i", ap, normal) / normal_squared
    projection = point[None, :] - plane_parameter[:, None] * normal
    projected = projection - a
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", projected, ab)
    d21 = np.einsum("ij,ij->i", projected, ac)
    denominator = d00 * d11 - d01 * d01
    barycentric_b = (d11 * d20 - d01 * d21) / denominator
    barycentric_c = (d00 * d21 - d01 * d20) / denominator
    barycentric_a = 1.0 - barycentric_b - barycentric_c
    barycentric_tolerance = 256.0 * np.finfo(float).eps
    inside = ((barycentric_a >= -barycentric_tolerance)
              & (barycentric_b >= -barycentric_tolerance)
              & (barycentric_c >= -barycentric_tolerance))
    plane_distance_squared = plane_parameter * plane_parameter * normal_squared

    def segment_distance_squared(start, end):
        edge = end - start
        edge_squared = np.einsum("ij,ij->i", edge, edge)
        parameter = np.einsum(
            "ij,ij->i", point[None, :] - start, edge) / edge_squared
        parameter = np.clip(parameter, 0.0, 1.0)
        delta = point[None, :] - (start + parameter[:, None] * edge)
        return np.einsum("ij,ij->i", delta, delta)

    edge_distance_squared = np.minimum.reduce((
        segment_distance_squared(a, b),
        segment_distance_squared(b, c),
        segment_distance_squared(c, a),
    ))
    return np.sqrt(np.maximum(
        np.where(inside, plane_distance_squared, edge_distance_squared), 0.0))


def _maximum_point_to_surface_distance(points, triangles, maximum_distance):
    """Certify each point against nearby triangles without an all-pairs allocation."""
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3
            or triangles.ndim != 3 or triangles.shape[1:] != (3, 3)
            or triangles.shape[0] == 0 or np.any(~np.isfinite(points))
            or np.any(~np.isfinite(triangles))):
        raise ValueError("surface-distance certification requires finite nonempty geometry")
    center = np.mean(triangles, axis=1)
    radius = np.max(
        np.linalg.norm(triangles - center[:, None, :], axis=2), axis=1)
    maximum_radius = float(np.max(radius))
    coordinate_scale = max(
        float(np.max(np.abs(points))) if points.size else 0.0,
        float(np.max(np.abs(triangles))), 1.0)
    roundoff = 256.0 * np.finfo(float).eps * coordinate_scale
    tree = cKDTree(center)
    candidates = tree.query_ball_point(
        points, r=float(maximum_distance) + maximum_radius + roundoff)
    maximum_nearest = 0.0
    for point, local in zip(points, candidates):
        if not local:
            return np.inf
        distance = _point_to_triangles_distance(
            point, triangles[np.asarray(local, dtype=int)])
        nearest = float(np.min(distance))
        maximum_nearest = max(maximum_nearest, nearest)
        if nearest > float(maximum_distance) + roundoff:
            return nearest
    return maximum_nearest


def _periodic_remap_shifts(periodic_lengths):
    """Return the nearest-image shifts for zero or more declared periodic axes."""
    if periodic_lengths is None:
        return np.zeros((1, 3), dtype=float), (None, None, None)
    if len(periodic_lengths) != 3:
        raise ValueError("surface remap periodic lengths must have three entries")
    normalized = []
    choices = []
    for length in periodic_lengths:
        if length is None or float(length) == 0.0:
            normalized.append(None)
            choices.append((0.0,))
        else:
            value = float(length)
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError("surface remap periodic lengths must be positive")
            normalized.append(value)
            choices.append((-value, 0.0, value))
    return np.asarray(tuple(product(*choices)), dtype=float), tuple(normalized)


def _periodic_images(values, shifts):
    values = np.asarray(values, dtype=float)
    shift_shape = (shifts.shape[0],) + (1,) * (values.ndim - 1) + (3,)
    images = values[None, ...] + shifts.reshape(shift_shape)
    source = np.tile(np.arange(values.shape[0], dtype=int), shifts.shape[0])
    return images.reshape((-1,) + values.shape[1:]), source


def _periodic_unique_neighbors(old_points, new_points, count, shifts):
    """Query wrapped images while returning each physical source face at most once."""
    images, source = _periodic_images(old_points, shifts)
    raw_count = min(images.shape[0], int(count) * shifts.shape[0])
    distance, image_index = cKDTree(images).query(new_points, k=raw_count)
    if raw_count == 1:
        distance = np.asarray(distance)[:, None]
        image_index = np.asarray(image_index)[:, None]
    selected_distance = np.empty((len(new_points), int(count)), dtype=float)
    selected_source = np.empty((len(new_points), int(count)), dtype=int)
    for row in range(len(new_points)):
        used = set(); output = 0
        for candidate_distance, candidate_image in zip(distance[row], image_index[row]):
            candidate_source = int(source[int(candidate_image)])
            if candidate_source in used:
                continue
            selected_distance[row, output] = float(candidate_distance)
            selected_source[row, output] = candidate_source
            used.add(candidate_source)
            output += 1
            if output == int(count):
                break
        if output != int(count):
            raise RuntimeError("periodic remap neighbor query lost physical source faces")
    return selected_distance, selected_source


def conservative_remap_surface_state(
        state, old_centroid, old_area, old_material, new_centroid, new_area, new_material, *,
        dx, mesh_length_unit_m, neighbor_count=4, maximum_distance=None,
        old_triangles=None, periodic_lengths=None):
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
    old_triangles = (
        None if old_triangles is None else np.asarray(old_triangles, dtype=float))
    if (old_centroid.ndim != 2 or old_centroid.shape[1] != 3
            or new_centroid.ndim != 2 or new_centroid.shape[1] != 3
            or old_area.shape != (old_centroid.shape[0],)
            or new_area.shape != (new_centroid.shape[0],)
            or old_material.shape != old_area.shape or new_material.shape != new_area.shape
            or any(value.shape != old_area.shape for value in old_values.values())
            or np.any(old_area <= 0.0) or np.any(new_area <= 0.0)
            or (old_triangles is not None
                and (old_triangles.shape != (old_centroid.shape[0], 3, 3)
                     or np.any(~np.isfinite(old_triangles))))):
        raise ValueError("invalid surface-state remap geometry")
    if maximum_distance is None:
        maximum_distance = 2.0 * float(dx)
    if not np.isfinite(maximum_distance) or maximum_distance <= 0.0:
        raise ValueError("maximum remap distance must be positive")
    periodic_shifts, normalized_periodic_lengths = _periodic_remap_shifts(
        periodic_lengths)
    union_triangle_images = (
        None if old_triangles is None else
        _periodic_images(old_triangles, periodic_shifts)[0])
    output = {name: np.zeros(new_area.shape) for name in old_values}
    maximum_nearest = 0.0; maximum_centroid_nearest = 0.0
    material_diagnostics = {}
    physical_area_scale = float(mesh_length_unit_m) ** 2
    for material in sorted(set(old_material) | set(new_material)):
        old_index = np.where(old_material == material)[0]
        new_index = np.where(new_material == material)[0]
        if old_index.size == 0 or new_index.size == 0:
            raise ValueError(
                "material surface appeared or disappeared; initialize/retire state explicitly")
        count = min(int(neighbor_count), old_index.size)
        distance, local = _periodic_unique_neighbors(
            old_centroid[old_index], new_centroid[new_index], count, periodic_shifts)
        centroid_nearest = float(np.max(distance[:, 0]))
        maximum_centroid_nearest = max(maximum_centroid_nearest, centroid_nearest)
        material_triangle_images = (
            None if old_triangles is None else
            _periodic_images(old_triangles[old_index], periodic_shifts)[0])
        nearest = (
            centroid_nearest if old_triangles is None else
            _maximum_point_to_surface_distance(
                new_centroid[new_index], material_triangle_images, maximum_distance))
        maximum_nearest = max(maximum_nearest, nearest)
        if nearest > maximum_distance:
            union_nearest = (
                nearest if old_triangles is None else
                _maximum_point_to_surface_distance(
                    new_centroid[new_index], union_triangle_images, maximum_distance))
            raise ValueError(
                f"surface remap distance {nearest:g} exceeds {maximum_distance:g} "
                f"using {'centroid' if old_triangles is None else 'point-to-triangle'} metric "
                f"for material {int(material)}; union-surface distance={union_nearest:g}")
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
        maximum_nearest_centroid_distance=float(maximum_centroid_nearest),
        distance_metric=(
            "centroid" if old_triangles is None else "point_to_material_triangle"),
        periodic_lengths=normalized_periodic_lengths,
        maximum_allowed_distance=float(maximum_distance),
        materials=material_diagnostics)


def _remap_surface_state_with_indexed_transfer(
        state, old_surface, new_surface, *, neighbor_count, maximum_distance,
        mesh_length_unit_m):
    """Apply the shared exact-distance/indexed predecessor operator to one state."""
    if (not hasattr(state, "conservative_surface_fields")
            or not hasattr(state, "conservative_surface_upper_bounds")
            or not hasattr(state, "with_conservative_surface_fields")):
        raise TypeError("surface state does not implement the conservative remap contract")
    fields = {
        str(name): np.asarray(value, dtype=float)
        for name, value in state.conservative_surface_fields().items()}
    upper = dict(state.conservative_surface_upper_bounds())
    modes = (
        dict(state.surface_field_remap_modes())
        if hasattr(state, "surface_field_remap_modes")
        else {name: "conservative" for name in fields})
    if (not fields or set(upper) != set(fields) or set(modes) != set(fields)
            or any(mode not in {"conservative", "intensive"} for mode in modes.values())
            or any(value.shape != (len(old_surface.faces),) for value in fields.values())):
        raise ValueError("surface-state remap fields and upper bounds must match")
    transfer = build_surface_transfer_3d(
        old_surface, new_surface, neighbor_count=neighbor_count,
        maximum_distance=maximum_distance)
    output = {}
    applications = {}
    for name, value in fields.items():
        if modes[name] == "conservative":
            application = transfer.apply_extensive(value, upper_bound=upper[name])
        else:
            application = transfer.apply_intensive(
                value, lower_bound=0.0, upper_bound=upper[name])
        output[name] = application.values
        applications[name] = application

    physical_area_scale = float(mesh_length_unit_m) ** 2
    material_diagnostics = {}
    for material in sorted(set(old_surface.face_material_id.tolist())):
        old_selected = old_surface.face_material_id == material
        new_selected = new_surface.face_material_id == material
        targets = {}
        achieved = {}
        residuals = []
        for name, application in applications.items():
            receipt = application.material_integrals[int(material)]
            if modes[name] == "conservative":
                targets[name] = receipt["old_area_integral"] * physical_area_scale
                residuals.append(receipt["relative_difference"])
            achieved[name] = receipt["new_area_integral"] * physical_area_scale
        material_diagnostics[int(material)] = dict(
            old_face_count=int(np.count_nonzero(old_selected)),
            new_face_count=int(np.count_nonzero(new_selected)),
            old_area_m2=float(
                np.sum(old_surface.face_area[old_selected]) * physical_area_scale),
            new_area_m2=float(
                np.sum(new_surface.face_area[new_selected]) * physical_area_scale),
            target_field_integrals=targets,
            remapped_field_integrals=achieved,
            field_remap_modes=dict(modes),
            max_relative_conservation_residual=float(max(residuals, default=0.0)))
    return state.with_conservative_surface_fields(output), dict(
        method="material_local_indexed_exact_surface_knn",
        transfer_fingerprint=transfer.fingerprint,
        neighbor_count=int(neighbor_count),
        maximum_nearest_distance=float(transfer.maximum_exact_surface_distance),
        maximum_nearest_centroid_distance=float(
            transfer.maximum_nearest_centroid_distance),
        distance_metric="indexed_exact_point_to_material_triangle",
        periodic_lengths=old_surface.periodic_lengths,
        maximum_allowed_distance=float(maximum_distance),
        materials=material_diagnostics)


def _remap_surface_state_with_overlap_transfer(
        state, newly_exposed_state, transfer, *, maximum_distance,
        mesh_length_unit_m, method):
    """Apply one sparse overlap authority with explicit fresh-surface closures."""
    old_surface = transfer.old_surface
    new_surface = transfer.new_surface
    required = (
        "conservative_surface_fields", "conservative_surface_upper_bounds",
        "with_conservative_surface_fields")
    if any(not hasattr(state, name) for name in required):
        raise TypeError("surface state does not implement the conservative remap contract")
    if any(not hasattr(newly_exposed_state, name) for name in required):
        raise TypeError("newly exposed surface state does not implement the remap contract")
    fields = {
        str(name): np.asarray(value, dtype=float)
        for name, value in state.conservative_surface_fields().items()}
    fresh = {
        str(name): np.asarray(value, dtype=float)
        for name, value in newly_exposed_state.conservative_surface_fields().items()}
    upper = dict(state.conservative_surface_upper_bounds())
    modes = (
        dict(state.surface_field_remap_modes())
        if hasattr(state, "surface_field_remap_modes")
        else {name: "conservative" for name in fields})
    if (not fields or set(fresh) != set(fields) or set(upper) != set(fields)
            or set(modes) != set(fields)
            or any(mode not in {"conservative", "intensive"} for mode in modes.values())
            or any(value.shape != (len(old_surface.faces),) for value in fields.values())
            or any(value.shape != (len(new_surface.faces),) for value in fresh.values())):
        raise ValueError("surface-state overlap fields and fresh-surface closure must match")
    output = {}
    applications = {}
    for name, value in fields.items():
        if modes[name] == "conservative":
            application = transfer.apply_extensive(
                value, newly_exposed_density=fresh[name])
        else:
            application = transfer.apply_intensive(
                value, uncovered_fill=fresh[name])
        result = application.values
        tolerance = 128.0 * np.finfo(float).eps * max(
            float(np.max(np.abs(result), initial=0.0)), 1.0)
        if np.any(result < -tolerance):
            raise RuntimeError("surface overlap produced a negative state density")
        if upper[name] is not None and np.any(result > float(upper[name]) + tolerance):
            raise ValueError("surface contraction exceeds bounded extensive capacity")
        output[name] = np.maximum(result, 0.0)
        applications[name] = application

    physical_area_scale = float(mesh_length_unit_m) ** 2
    material_diagnostics = {}
    for material in sorted(set(old_surface.face_material_id.tolist())):
        old_selected = old_surface.face_material_id == material
        new_selected = new_surface.face_material_id == material
        targets = {}
        achieved = {}
        removed = {}
        exposed = {}
        residuals = []
        for name, application in applications.items():
            if modes[name] == "conservative":
                ledger = application.material_ledger[int(material)]
                targets[name] = ledger["old_inventory"] * physical_area_scale
                achieved[name] = ledger["new_inventory"] * physical_area_scale
                removed[name] = ledger["removed_inventory"] * physical_area_scale
                exposed[name] = ledger["newly_exposed_inventory"] * physical_area_scale
                residuals.append(ledger["relative_balance_error"])
            else:
                achieved[name] = float(np.dot(
                    output[name][new_selected],
                    new_surface.face_area[new_selected])) * physical_area_scale
        material_diagnostics[int(material)] = dict(
            old_face_count=int(np.count_nonzero(old_selected)),
            new_face_count=int(np.count_nonzero(new_selected)),
            old_area_m2=float(
                np.sum(old_surface.face_area[old_selected]) * physical_area_scale),
            new_area_m2=float(
                np.sum(new_surface.face_area[new_selected]) * physical_area_scale),
            target_field_integrals=targets,
            remapped_field_integrals=achieved,
            removed_field_integrals=removed,
            newly_exposed_field_integrals=exposed,
            removed_area_m2=float(
                np.sum(transfer.old_uncovered_area[old_selected]) * physical_area_scale),
            newly_exposed_area_m2=float(
                np.sum(transfer.new_uncovered_area[new_selected]) * physical_area_scale),
            field_remap_modes=dict(modes),
            max_relative_conservation_residual=float(max(residuals, default=0.0)))
    return state.with_conservative_surface_fields(output), dict(
        method=str(method),
        transfer_fingerprint=transfer.fingerprint,
        total_overlap_area_m2=float(
            np.sum(transfer.overlap_area) * physical_area_scale),
        total_removed_area_m2=float(
            np.sum(transfer.old_uncovered_area) * physical_area_scale),
        total_newly_exposed_area_m2=float(
            np.sum(transfer.new_uncovered_area) * physical_area_scale),
        maximum_allowed_distance=float(maximum_distance),
        periodic_lengths=old_surface.periodic_lengths,
        fresh_surface_closure="mechanism_declared_initial_state",
        geometry_receipt=dict(transfer.geometry_receipt),
        materials=material_diagnostics)


def _remap_surface_state_with_partitioned_overlap(
        state, newly_exposed_state, old_surface, new_surface, *, maximum_distance,
        mesh_length_unit_m):
    """Apply piecewise-planar exact overlap with explicit fresh-surface closures."""
    transfer = build_partitioned_surface_overlap_transfer_3d(
        old_surface, new_surface, maximum_normal_distance=maximum_distance)
    remapped, diagnostics = _remap_surface_state_with_overlap_transfer(
        state, newly_exposed_state, transfer,
        maximum_distance=maximum_distance,
        mesh_length_unit_m=mesh_length_unit_m,
        method="material_local_partitioned_exact_overlap")
    diagnostics.update(
        old_patch_count=transfer.old_patch_count,
        new_patch_count=transfer.new_patch_count,
        candidate_patch_pair_count=transfer.candidate_patch_pair_count,
        positive_patch_pair_count=transfer.positive_patch_pair_count,
        patch_receipts=tuple(dict(item) for item in transfer.patch_receipts))
    return remapped, diagnostics


def _remap_surface_state_with_common_refinement(
        state, newly_exposed_state, old_surface, new_surface, *, maximum_distance,
        mesh_length_unit_m):
    """Apply indexed tangent common refinement to nearby moving surfaces."""
    transfer = build_surface_common_refinement_transfer_3d(
        old_surface, new_surface, maximum_normal_distance=maximum_distance)
    remapped, diagnostics = _remap_surface_state_with_overlap_transfer(
        state, newly_exposed_state, transfer,
        maximum_distance=maximum_distance,
        mesh_length_unit_m=mesh_length_unit_m,
        method="material_local_tangent_common_refinement")
    diagnostics.update(
        candidate_pair_count=transfer.candidate_pair_count,
        aligned_pair_count=transfer.aligned_pair_count,
        positive_pair_image_count=transfer.positive_pair_image_count,
        combined_pair_count=transfer.combined_pair_count,
        minimum_normal_dot=transfer.minimum_normal_dot)
    return remapped, diagnostics


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
            energetic.append(population.remap_faces(
                old_to_new, selected_face.size))
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
        "maximum_iterations", "maximum_rays_per_face", "source_sampling",
        "form_factor_backend", "deterministic_extruded_options",
        "overlap_skip_depth_limit", "launch_surface_distance",
        "unclassified_ray_budget",
    }
    unknown = set(options) - allowed
    if unknown:
        raise ValueError("unknown neutral radiosity options: " + ", ".join(sorted(unknown)))
    material_probability = dict(options.pop(
        "nonetchable_reaction_probability_by_material", {}))
    solver_tolerance = float(options.pop("relative_tolerance", 1e-10))
    maximum_iterations = int(options.pop("maximum_iterations", 500))
    backend = str(options.pop("form_factor_backend", "scrambled_qmc_3d"))
    deterministic_options = dict(options.pop("deterministic_extruded_options", {}))
    if backend not in ("scrambled_qmc_3d", "deterministic_extruded_2d"):
        raise ValueError(
            "neutral radiosity form_factor_backend must be 'scrambled_qmc_3d' or "
            "'deterministic_extruded_2d'")
    deterministic_exchange = None
    deterministic_field_relative_tolerance = None
    deterministic_field_absolute_tolerance = None
    visibility_receipt = None
    refinement = []
    if backend == "scrambled_qmc_3d":
        if deterministic_options:
            raise ValueError(
                "deterministic_extruded_options require the deterministic_extruded_2d backend")
        initial_rays_per_face = int(options.pop("rays_per_face", 64))
        maximum_rays_per_face = int(options.pop(
            "maximum_rays_per_face", 8 * initial_rays_per_face))
        if (initial_rays_per_face <= 0
                or initial_rays_per_face & (initial_rays_per_face - 1)
                or maximum_rays_per_face < initial_rays_per_face
                or maximum_rays_per_face & (maximum_rays_per_face - 1)):
            raise ValueError(
                "neutral radiosity initial/maximum rays must be positive powers of two")
        if "domain_size" not in options:
            options["domain_size"] = (np.asarray(geometry.phi.shape) - 1) * geometry.dx
        if "ray_offset" not in options:
            options["ray_offset"] = 1e-3 * geometry.dx
        if "overlap_skip_depth_limit" not in options:
            # Bounded exit from piecewise-linear sheet interpenetration: marching-cubes facets
            # deviate from the trilinear isosurface by O(dx^2), so an artifact-zone launch lies
            # strictly sub-cell from the authority surface; half a cell bounds it with margin.
            options["overlap_skip_depth_limit"] = 0.5 * geometry.dx
        if "launch_surface_distance" not in options:
            # The trilinear level set is the authority surface (normals are already oriented
            # by its exact gradient); |phi| at a launch point measures how far the
            # piecewise-linear facet strays from it.
            phi_grid = np.asarray(geometry.phi, dtype=float)
            grid_shape = np.asarray(phi_grid.shape)
            dx_local = float(geometry.dx)

            def _trilinear_surface_distance(points, _phi=phi_grid, _shape=grid_shape,
                                            _dx=dx_local):
                q = np.asarray(points, dtype=float) / _dx
                base = np.clip(np.floor(q).astype(int), 0, _shape - 2)
                frac = q - base
                value = np.zeros(len(q))
                for a in (0, 1):
                    for b in (0, 1):
                        for c in (0, 1):
                            weight = (
                                (frac[:, 0] if a else 1.0 - frac[:, 0])
                                * (frac[:, 1] if b else 1.0 - frac[:, 1])
                                * (frac[:, 2] if c else 1.0 - frac[:, 2]))
                            value += weight * _phi[
                                base[:, 0] + a, base[:, 1] + b, base[:, 2] + c]
                return np.abs(value)

            options["launch_surface_distance"] = _trilinear_surface_distance
        if "unclassified_ray_budget" not in options:
            # Bounded disposition of rays that survive every proof-based recovery
            # unclassified: ledger the weight as lost, refuse only when the per-row or
            # global unclassified fraction exceeds these declared budgets (the signature
            # of a structural defect rather than a numerical-tail straggler).
            options["unclassified_ray_budget"] = (0.01, 0.001)
        rays_per_face = initial_rays_per_face
    else:
        incompatible = set(options) & {
            "rays_per_face", "maximum_rays_per_face", "seed", "ray_offset",
            "source_sampling",
        }
        if incompatible:
            raise ValueError(
                "deterministic extruded exchange does not accept sampling controls: "
                + ", ".join(sorted(incompatible)))
        if not bool(options.pop("periodic_lateral", False)):
            raise ValueError(
                "deterministic extruded exchange requires periodic_lateral=True")
        domain_size = np.asarray(options.pop(
            "domain_size", (np.asarray(geometry.phi.shape) - 1) * geometry.dx), dtype=float)
        if domain_size.shape != (3,) or np.any(~np.isfinite(domain_size)):
            raise ValueError("deterministic extruded exchange requires a three-axis domain_size")
        deterministic_allowed = {
            "extrusion_axis", "extrusion_length", "geometry_tolerance",
            "normal_tolerance", "area_relative_tolerance",
            "exchange_relative_tolerance", "exchange_absolute_tolerance",
            "minimum_refinement_level", "maximum_refinement_level",
            "field_relative_tolerance", "field_absolute_tolerance",
        }
        deterministic_unknown = set(deterministic_options) - deterministic_allowed
        if deterministic_unknown:
            raise ValueError(
                "unknown deterministic extruded exchange options: "
                + ", ".join(sorted(deterministic_unknown)))
        extrusion_axis = int(deterministic_options.pop("extrusion_axis", 1))
        deterministic_field_relative_tolerance = float(
            deterministic_options.pop("field_relative_tolerance", 1e-8))
        deterministic_field_absolute_tolerance = float(
            deterministic_options.pop("field_absolute_tolerance", 0.0))
        if (deterministic_field_relative_tolerance < 0.0
                or deterministic_field_absolute_tolerance < 0.0):
            raise ValueError("deterministic extrusion field tolerances must be nonnegative")
        deterministic_options.setdefault("extrusion_axis", extrusion_axis)
        deterministic_options.setdefault("extrusion_length", float(domain_size[extrusion_axis]))
        deterministic_exchange = build_extruded_triangle_exchange_3d(
            verts, faces, _surface_gas_normals(verts, faces, centroids, geometry),
            **deterministic_options)
        factors = deterministic_exchange.form_factors
        rays_per_face = 0
        maximum_rays_per_face = 0
        if options:
            raise ValueError(
                "unused deterministic neutral radiosity options: "
                + ", ".join(sorted(options)))
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

    physical_area = np.asarray(areas) * geometry.mesh_length_unit_m ** 2
    while True:
        if backend == "scrambled_qmc_3d":
            visibility_receipt = estimate_diffuse_form_factors_3d(
                verts, faces, centroids, _surface_gas_normals(
                    verts, faces, centroids, geometry),
                rays_per_face=rays_per_face, device=transport_device,
                return_visibility_receipt=True, **options)
            factors = visibility_receipt.form_factors
        neutral_flux = {}
        diagnostics = {}
        try:
            for name, direct in transport.surface_fluxes.neutral_flux_m2_s.items():
                if name not in reaction_probability:
                    neutral_flux[name] = direct
                    continue
                direct_for_solve = np.asarray(direct, dtype=float)
                probability_for_solve = reaction_probability[name]
                if deterministic_exchange is not None:
                    direct_group = deterministic_exchange.certify_face_field(
                        direct_for_solve,
                        relative_tolerance=deterministic_field_relative_tolerance,
                        absolute_tolerance=deterministic_field_absolute_tolerance)
                    probability_group = deterministic_exchange.certify_face_field(
                        probability_for_solve,
                        relative_tolerance=deterministic_field_relative_tolerance,
                        absolute_tolerance=deterministic_field_absolute_tolerance)
                    direct_for_solve = direct_group[
                        deterministic_exchange.face_group_index]
                    probability_for_solve = probability_group[
                        deterministic_exchange.face_group_index]
                factor_diagnostics = dict(
                    form_factor_backend=backend,
                    form_factor_rays_per_face=int(rays_per_face),
                    form_factor_refinement_count=len(refinement),
                    form_factor_refinement=tuple(refinement))
                if visibility_receipt is not None:
                    factor_diagnostics.update(
                        visibility_mode=visibility_receipt.visibility_mode,
                        visibility_ray_count=visibility_receipt.ray_count,
                        visibility_float64_evaluated_count=(
                            visibility_receipt.float64_evaluated_count),
                        visibility_recovered_hit_count=(
                            visibility_receipt.float64_recovered_hit_count),
                        visibility_open_escape_count=(
                            visibility_receipt.open_escape_count),
                        visibility_maximum_wrap_count=(
                            visibility_receipt.maximum_wrap_count),
                        visibility_derived_horizon_extension_count=(
                            visibility_receipt.derived_horizon_extension_count),
                        visibility_initial_maximum_wraps=(
                            visibility_receipt.initial_maximum_wraps),
                        visibility_final_maximum_wraps=(
                            visibility_receipt.final_maximum_wraps),
                        visibility_launch_inset_count=(
                            visibility_receipt.launch_inset_count),
                        visibility_centroid_limit_count=(
                            visibility_receipt.centroid_limit_count),
                        visibility_source_relaunch_count=(
                            visibility_receipt.source_relaunch_count),
                        visibility_maximum_source_relaunch_distance=(
                            visibility_receipt.maximum_source_relaunch_distance),
                        visibility_overlap_skip_count=(
                            visibility_receipt.overlap_skip_count),
                        visibility_maximum_overlap_skip_depth=(
                            visibility_receipt.maximum_overlap_skip_depth),
                        visibility_unclassified_ray_count=(
                            visibility_receipt.unclassified_ray_count),
                        visibility_source_support_face_count=(
                            visibility_receipt.source_support_face_count),
                        visibility_source_support_area_fraction=(
                            visibility_receipt.source_support_area_fraction),
                        visibility_maximum_source_support_distance=(
                            visibility_receipt.maximum_source_support_distance))
                if deterministic_exchange is not None:
                    factor_diagnostics.update(
                        form_factor_fingerprint=deterministic_exchange.fingerprint,
                        form_factor_group_count=deterministic_exchange.group_count,
                        form_factor_strip_count=deterministic_exchange.strip_count,
                        form_factor_maximum_group_area_relative_error=(
                            deterministic_exchange.maximum_group_area_relative_error),
                        form_factor_maximum_area_reciprocity_error=(
                            deterministic_exchange.maximum_area_reciprocity_error),
                        form_factor_maximum_estimated_exchange_error=(
                            deterministic_exchange.line_exchange.
                            maximum_estimated_absolute_error),
                        form_factor_maximum_refinement_level=int(
                            np.max(deterministic_exchange.line_exchange.refinement_level)),
                        extrusion_field_relative_tolerance=(
                            deterministic_field_relative_tolerance),
                        extrusion_field_absolute_tolerance=(
                            deterministic_field_absolute_tolerance))
                # A globally inert population is causally disconnected from every surface-state
                # and profile equation.  Its repeated diffuse collision count can be arbitrarily
                # large in a high-AR, nearly pinched feature, and a finite form-factor sample can
                # manufacture a singular closed class when the true escape cone is merely tiny.
                # Analytically marginalize that null channel: preserve the first-hit diagnostic,
                # account every launched particle as eventually escaping, and solve radiosity only
                # for species with a nonzero chance to modify the modeled surface.
                if not np.any(probability_for_solve > 0.0):
                    source_rate = float(np.sum(direct_for_solve * physical_area))
                    neutral_flux[name] = direct_for_solve
                    diagnostics[name] = dict(
                        source_rate_s=source_rate,
                        reacted_rate_s=0.0,
                        escaped_rate_s=source_rate,
                        relative_balance_error=0.0,
                        relative_linear_residual=0.0,
                        solver_method="analytic_zero_reaction_elision",
                        iteration_count=0,
                        inactive_face_count=int(len(faces)),
                        repeated_incident_flux_elided=True,
                        **factor_diagnostics,
                    )
                    continue
                solution = solve_diffuse_neutral_radiosity_3d(
                    direct_for_solve, physical_area,
                    factors.source_face, factors.target_face,
                    factors.transfer_fraction, factors.escape_fraction,
                    probability_for_solve, relative_tolerance=solver_tolerance,
                    maximum_iterations=maximum_iterations)
                incident = solution.incident_flux_m2_s
                if deterministic_exchange is not None:
                    incident = deterministic_exchange.area_weighted_group_mean(incident)[
                        deterministic_exchange.face_group_index]
                neutral_flux[name] = incident
                diagnostics[name] = dict(
                    source_rate_s=solution.source_rate_s,
                    reacted_rate_s=solution.reacted_rate_s,
                    escaped_rate_s=solution.escaped_rate_s,
                    relative_balance_error=solution.relative_balance_error,
                    relative_linear_residual=solution.relative_linear_residual,
                    solver_method=solution.solver_method,
                    iteration_count=solution.iteration_count,
                    inactive_face_count=solution.inactive_face_count,
                    repeated_incident_flux_elided=False,
                    **factor_diagnostics)
        except DiffuseNeutralNoSinkError as error:
            if backend == "deterministic_extruded_2d":
                raise RuntimeError(
                    f"deterministic neutral radiosity for {name!r} contains a "
                    "source-fed no-sink class; sampling refinement cannot repair an exact "
                    "closed transport class") from error
            if rays_per_face >= maximum_rays_per_face:
                raise RuntimeError(
                    f"neutral radiosity for {name!r} retained a source-fed no-sink "
                    f"class through {rays_per_face} rays/face") from error
            next_rays = min(2 * rays_per_face, maximum_rays_per_face)
            refinement.append({
                "species": str(name),
                "from_rays_per_face": int(rays_per_face),
                "to_rays_per_face": int(next_rays),
                "closed_class_face_count": int(error.face_count),
                "classification": "nested_form_factor_refinement",
            })
            rays_per_face = next_rays
            continue
        break
    limitations = tuple(
        item for item in transport.known_limitations
        if item != "no surface reflection or neutral re-emission") + (
        "neutral re-emission is diffuse with material/state reaction probabilities",
    ) + ((
        "globally zero-reaction neutral populations retain first-hit flux only; "
        "their chemically irrelevant repeated collision count is analytically elided",
    ) if any(
        item.get("repeated_incident_flux_elided", False)
        for item in diagnostics.values()) else ()) + ((
        "deterministic diffuse exchange assumes exact translational invariance along the "
        "declared periodic extrusion axis",
    ) if backend == "deterministic_extruded_2d" else ())
    updated = BoundaryTransport3DResult(
        SurfaceFluxes(neutral_flux, transport.surface_fluxes.energetic_fluxes),
        transport.hit_probability, transport.escape_probability,
        transport.truncation_probability,
        transport.transport_model + f" + flux_conservative_diffuse_radiosity[{backend}]",
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
        "ray_offset", "relative_tolerance", "maximum_iterations", "source_sampling",
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
        ballistic_transport="forward", ballistic_periodic_lateral=None,
        ballistic_face_quadrature_points=1, cfl_number=0.3, reinitialize=True,
        reinitialization_method="skfmm",
        topology_change_policy="refuse",
        surface_state_remap_backend="legacy_knn",
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
    ``ballistic_periodic_lateral`` declares periodic field-free first-hit transport independently
    of diffuse neutral radiosity.  ``None`` preserves the historical behavior by inheriting the
    radiosity periodic setting; an explicit boolean is the preferred production declaration.
    ``topology_change_policy='continue_gas_cavity'`` is deliberately narrow: only the periodic
    physical-volume signature's gas-cavity count may change. The existing material-local conservative
    remap then transfers surface history; solid components, per-material components, and domain
    breakthrough must remain invariant.
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
    if (ballistic_periodic_lateral is not None
            and not isinstance(ballistic_periodic_lateral, (bool, np.bool_))):
        raise ValueError("ballistic_periodic_lateral must be boolean or None")
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
    if ballistic_periodic_lateral is not None and (
            precomputed_transport is not None or charging_poisson_system is not None
            or nodal_potential_v is not None):
        raise ValueError(
            "ballistic_periodic_lateral controls only field-free first-hit transport; "
            "use the field, charging, or precomputed transport periodic contract")
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
    if topology_change_policy not in ("refuse", "continue_gas_cavity"):
        raise ValueError(
            "topology_change_policy must be 'refuse' or 'continue_gas_cavity'")
    if surface_state_remap_backend not in (
            "legacy_knn", "indexed_knn", "partitioned_overlap", "common_refinement"):
        raise ValueError(
            "surface_state_remap_backend must be 'legacy_knn', 'indexed_knn', or "
            "'partitioned_overlap'/'common_refinement'")
    if ballistic_transport == "face_gather" and (
            charging_poisson_system is not None or nodal_potential_v is not None):
        raise ValueError("deterministic ballistic face gather does not yet trace electric fields")
    if charged_surface_response is not None and ballistic_transport == "face_gather":
        raise ValueError(
            "charged surface response requires forward impact-position lineage; "
            "face_gather currently preserves direction but not impact position")

    verts, faces, centroids, areas, face_material = (
        _extract_uniform_surface_arrays(geometry))
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
    periodic_ballistic = (
        periodic_neutral if ballistic_periodic_lateral is None
        else bool(ballistic_periodic_lateral))
    if periodic_neutral and not periodic_ballistic:
        raise ValueError(
            "periodic neutral radiosity requires periodic ballistic first-hit transport")
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
            bool(field_periodic_lateral) if nodal_potential_v is not None else periodic_ballistic))
        if periodic_neutral and not response_periodic:
            raise ValueError(
                "periodic neutral radiosity requires periodic response-enabled trajectories")
    charging_periodic = bool(
        charging_options is not None and charging_options.get("periodic_lateral", False))
    transport_periodic_lateral = bool(
        periodic_ballistic or response_periodic or bool(field_periodic_lateral)
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
        if periodic_ballistic:
            first_hit_options = dict(
                periodic_lateral=True,
                domain_size=(
                    (np.asarray(geometry.phi.shape) - 1) * geometry.dx
                    if radiosity_options is None
                    else radiosity_options.get(
                        "domain_size",
                        (np.asarray(geometry.phi.shape) - 1) * geometry.dx)))
        if ballistic_transport == "face_gather":
            transport = gather_boundary_state_ballistic_3d(
                boundary, species_role, verts, faces, areas, centroids,
                _surface_gas_normals(verts, faces, centroids, geometry),
                source_bounds=source_bounds, source_z=source_z,
                mesh_length_unit_m=geometry.mesh_length_unit_m,
                mesh_origin_m=geometry.mesh_origin_m,
                face_quadrature_points=ballistic_face_quadrature_points,
                periodic_lateral=periodic_ballistic,
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
    raw_maximum_face_speed = (
        float(np.max(np.abs(face_velocity))) if face_velocity.size else 0.0)
    phi = np.array(geometry.phi, copy=True)
    periodic_seam_projection = 0.0
    periodic_seam_velocity_projection = 0.0
    if profile_periodic_lateral:
        phi, correction = _project_periodic_lateral_endpoints(phi)
        periodic_seam_projection = max(periodic_seam_projection, correction)
    xs, ys, zs = geometry.coordinate_arrays
    extension_geometry = dict(phi=phi, dx=geometry.dx, xs=xs, ys=ys, zs=zs)
    # Extend only from the material surface that is actually evolving.  Including pinned mask
    # triangles with zero velocity lets them win the nearest-face query below a narrow opening and
    # numerically pins a physically bombarded floor after roughly one grid cell of motion.
    extension_velocity = face_velocity[active_face]
    extension_centroid = centroids[active_face]
    if profile_periodic_lateral:
        extension_velocity, extension_centroid = _periodic_lateral_surface_images(
            extension_velocity, extension_centroid,
            (np.asarray(geometry.phi.shape) - 1) * geometry.dx)
    extended_velocity = extend_velocity_3d(
        extension_velocity, extension_centroid,
        extension_geometry, 4.0 * geometry.dx)
    if profile_periodic_lateral:
        extended_velocity, correction = _project_periodic_lateral_endpoints(
            extended_velocity)
        periodic_seam_velocity_projection = max(
            periodic_seam_velocity_projection, correction)
        # The physical boundary condition declares the duplicate endpoints to be
        # one location, so projecting their two finite-sample velocity estimates is
        # the authoritative periodic operation. Report its implied displacement,
        # but do not mix a velocity estimator discrepancy into the geometry gate.
        # CFL substepping independently resolves the accepted projected velocity.
        if periodic_seam_projection > 0.25 * geometry.dx:
            raise RuntimeError(
                "periodic input-geometry seam projection exceeds one quarter cell; "
                f"geometry projection={periodic_seam_projection:.8g} mesh units; "
                "the stored input is not a resolved periodic field")
    # The level set consumes the grid-resolved extended field, not individual marching-cubes face
    # values. A vanishing-area sliver can carry an extreme flux density yet be nearest to no grid
    # node; letting that auxiliary quadrature atom set CFL or the outer coupling step makes the
    # timestep chase subgrid geometry that never enters the evolution operator.
    maximum_speed = (
        float(np.max(np.abs(extended_velocity))) if extended_velocity.size else 0.0)
    maximum_recession = max(
        float(np.max(extended_velocity)) if extended_velocity.size else 0.0, 0.0)
    maximum_growth = max(
        float(np.max(-extended_velocity)) if extended_velocity.size else 0.0, 0.0)
    displacement = maximum_speed * float(duration_s)
    substeps = max(
        1, int(np.ceil(displacement / (float(cfl_number) * geometry.dx))))
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
    if duration_s == 0.0 or material_levelsets is None:
        reassigned_unresolved_material_nodes = 0
    else:
        unresolved_material_mask, reassigned_unresolved_material_nodes = (
            _new_unresolved_subcell_material_component_mask(
                phi, output_material_id, geometry.material_id, etchable,
                periodic_lateral=profile_periodic_lateral))
        if reassigned_unresolved_material_nodes:
            material_levelsets, phi, output_material_id = (
                _restore_unresolved_material_ownership(
                    material_levelsets, unresolved_material_mask,
                    output_material_id, geometry.material_id, geometry.dx,
                    reinitialization_method, profile_periodic_lateral))
    if duration_s == 0.0:
        removed_unresolved_solid_cells = 0
        unresolved_solid_mask = np.zeros_like(output_material_id, dtype=bool)
    else:
        phi, removed_unresolved_solid_cells, unresolved_solid_mask = (
            _remove_unresolved_subcell_solid_components(
                phi, output_material_id, etchable, geometry.dx,
                periodic_lateral=profile_periodic_lateral)
        )
    if removed_unresolved_solid_cells:
        if material_levelsets is not None:
            material_levelsets, phi, output_material_id = (
                _apply_subcell_cleanup_to_material_levelsets(
                    material_levelsets, unresolved_solid_mask, output_material_id,
                    etchable, geometry.dx, reinitialization_method,
                    profile_periodic_lateral))
            _, remaining_unresolved_cells, remaining_unresolved_mask = (
                _remove_unresolved_subcell_solid_components(
                    phi, output_material_id, etchable, geometry.dx,
                    periodic_lateral=profile_periodic_lateral))
            if remaining_unresolved_cells:
                remaining_coordinates = tuple(
                    tuple(int(value) for value in index)
                    for index in np.argwhere(remaining_unresolved_mask)[:12])
                remaining_owners = tuple(sorted(
                    int(value) for value in np.unique(
                        output_material_id[remaining_unresolved_mask])))
                raise RuntimeError(
                    "material-layer topology update left an unresolved subcell component; "
                    f"node_count={remaining_unresolved_cells}; "
                    f"owners={remaining_owners}; coordinates={remaining_coordinates}")
        else:
            phi = _redistance_feature_field(
                phi, geometry.dx, reinitialization_method,
                periodic_lateral=profile_periodic_lateral)
            phi[pinned] = geometry.phi[pinned]
    if duration_s == 0.0:
        unresolved_gas_mask = np.zeros_like(phi, dtype=bool)
        filled_unresolved_gas_cavity_cells = 0
    else:
        unresolved_gas_mask, filled_unresolved_gas_cavity_cells = (
            _unresolved_subcell_gas_cavity_mask(
                phi, periodic_lateral=profile_periodic_lateral))
    if filled_unresolved_gas_cavity_cells:
        if material_levelsets is not None:
            material_levelsets, phi, output_material_id = (
                _apply_subcell_gas_fill_to_material_levelsets(
                    material_levelsets, unresolved_gas_mask, etchable,
                    geometry.dx, reinitialization_method,
                    profile_periodic_lateral))
        else:
            solid_material = tuple(sorted({
                int(value) for value in np.unique(output_material_id[phi > 0.0])
                if int(value) > 0}))
            if len(solid_material) != 1 or solid_material[0] not in etchable:
                raise RuntimeError(
                    "subcell gas-cavity fill without material level sets requires one "
                    "evolving solid material")
            phi = np.asarray(phi, dtype=float).copy()
            phi[unresolved_gas_mask] = np.maximum(
                np.abs(phi[unresolved_gas_mask]), float(geometry.dx))
            phi = _redistance_feature_field(
                phi, geometry.dx, reinitialization_method,
                periodic_lateral=profile_periodic_lateral)
            output_material_id = np.where(phi >= 0.0, solid_material[0], 0)
        _, remaining_unresolved_gas_cells = (
            _unresolved_subcell_gas_cavity_mask(
                phi, periodic_lateral=profile_periodic_lateral))
        if remaining_unresolved_gas_cells:
            raise RuntimeError(
                "subcell gas-cavity cleanup did not restore resolved topology")
    if profile_periodic_lateral:
        phi, correction = _project_periodic_lateral_endpoints(phi)
        periodic_seam_projection = max(periodic_seam_projection, correction)

    output_geometry = FeatureGeometry3D(
        phi, output_material_id, geometry.dx, geometry.mesh_length_unit_m,
        geometry.mesh_origin_m, material_levelsets=material_levelsets)
    (next_verts, next_faces, next_centroids, next_areas,
     next_face_material) = _extract_uniform_surface_arrays(output_geometry)
    next_active_face = np.where(np.isin(next_face_material, etchable))[0]
    if next_active_face.size == 0:
        raise ValueError("etch step removed every requested material surface")
    old_mesh_topology = _surface_topology_signature(faces, active_face)
    next_mesh_topology = _surface_topology_signature(next_faces, next_active_face)
    topology_method = (
        "periodic_xy_component_cavity_breakthrough_v1"
        if profile_periodic_lateral else "bounded_volume_euler_v1")
    topology_operator = (
        _periodic_physical_volume_topology_signature
        if profile_periodic_lateral else _physical_volume_topology_signature)
    old_topology = topology_operator(geometry, etchable)
    next_topology = topology_operator(output_geometry, etchable)
    topology_event = None
    if old_topology != next_topology:
        changed_slice_topology = _changed_physical_slice_topology(
            geometry, output_geometry, etchable)
        message = (
            f"surface topology changed under {topology_method} from "
            f"{old_topology} to {next_topology}; "
            f"marching-cubes topology changed from {old_mesh_topology} "
            f"to {next_mesh_topology}; "
            f"component sizes changed from "
            f"{_physical_volume_component_sizes(geometry, etchable)} to "
            f"{_physical_volume_component_sizes(output_geometry, etchable)}; "
            f"periodic material-component sizes changed from "
            f"{_periodic_material_component_sizes(geometry, etchable)} to "
            f"{_periodic_material_component_sizes(output_geometry, etchable)}; "
            f"changed slice topology="
            f"{changed_slice_topology}; "
            "state transfer requires an explicit topology event")
        topology_error = SurfaceTopologyChangeError(
            message, method=topology_method,
            old_topology=old_topology, new_topology=next_topology,
            old_mesh_topology=old_mesh_topology,
            new_mesh_topology=next_mesh_topology,
            changed_slice_topology=changed_slice_topology)
        permitted_event = (
            topology_change_policy == "continue_gas_cavity"
            and profile_periodic_lateral
            and topology_error.event_kind in (
                "gas_cavity_enclosed", "gas_cavity_opened"))
        if not permitted_event:
            raise topology_error
        topology_event = {
            "accepted": True,
            "policy": str(topology_change_policy),
            "kind": topology_error.event_kind,
            "method": str(topology_method),
            "old_topology": tuple(old_topology),
            "new_topology": tuple(next_topology),
            "old_mesh_topology": tuple(old_mesh_topology),
            "new_mesh_topology": tuple(next_mesh_topology),
            "changed_slice_topology": dict(changed_slice_topology),
            "conservative_surface_state_remap_required": True,
        }
    remap_maximum_distance = displacement + 1.5 * geometry.dx
    remap_periodic_lengths = (
        tuple(((np.asarray(geometry.phi.shape) - 1) * geometry.dx)[:2]) + (None,)
        if profile_periodic_lateral else (None, None, None))
    if surface_state_remap_backend in (
            "indexed_knn", "partitioned_overlap", "common_refinement"):
        old_surface = TriangleSurface3D(
            verts, faces[active_face], face_material[active_face],
            periodic_lengths=remap_periodic_lengths)
        new_surface = TriangleSurface3D(
            next_verts, next_faces[next_active_face], next_face_material[next_active_face],
            periodic_lengths=remap_periodic_lengths)
        if surface_state_remap_backend == "indexed_knn":
            next_surface_state, remap_diagnostics = _remap_surface_state_with_indexed_transfer(
                surface.state, old_surface, new_surface, neighbor_count=4,
                maximum_distance=remap_maximum_distance,
                mesh_length_unit_m=geometry.mesh_length_unit_m)
        else:
            if material_resolved_mechanism:
                newly_exposed_state = mechanism.initial_state_by_material(
                    next_face_material[next_active_face])
            else:
                newly_exposed_state = mechanism.initial_state((next_active_face.size,))
            overlap_remap = (
                _remap_surface_state_with_partitioned_overlap
                if surface_state_remap_backend == "partitioned_overlap"
                else _remap_surface_state_with_common_refinement)
            next_surface_state, remap_diagnostics = (
                overlap_remap(
                    surface.state, newly_exposed_state, old_surface, new_surface,
                    maximum_distance=remap_maximum_distance,
                    mesh_length_unit_m=geometry.mesh_length_unit_m))
    else:
        next_surface_state, remap_diagnostics = conservative_remap_surface_state(
            surface.state, centroids[active_face], areas[active_face],
            face_material[active_face], next_centroids[next_active_face],
            next_areas[next_active_face], next_face_material[next_active_face],
            dx=geometry.dx, mesh_length_unit_m=geometry.mesh_length_unit_m,
            maximum_distance=remap_maximum_distance,
            old_triangles=verts[faces[active_face]],
            periodic_lengths=(
                remap_periodic_lengths if profile_periodic_lateral else None))
    next_mesh_fingerprint = _surface_mesh_fingerprint(
        next_verts, next_faces, next_active_face, next_face_material, output_geometry)
    remap_diagnostics = dict(
        remap_diagnostics, old_topology=old_topology, new_topology=next_topology,
        topology_method=topology_method,
        old_mesh_topology=old_mesh_topology, new_mesh_topology=next_mesh_topology,
        next_active_face_count=int(next_active_face.size),
        surface_state_remap_backend=str(surface_state_remap_backend),
        topology_change_policy=str(topology_change_policy),
        topology_event=topology_event)
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
    topology_limitation = (
        "physical volume-topology-changing surface steps are refused"
        if topology_change_policy == "refuse" else
        "only periodic gas-cavity enclosure/opening may continue through an explicit "
        "conservative state remap; all other physical volume-topology changes are refused")
    validity = FeatureStepValidity(
        within_declared_scope=not reasons,
        reasons=tuple(reasons),
        known_limitations=tuple(dict.fromkeys(transport_limitations)) + (
            "first-order material-local conservative surface-state remap with declared intensive-field exceptions",
            "new subcell material-label components below one resolved volume cell are suppressed and bounded",
            topology_limitation,
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
            raw_maximum_face_velocity_m_s=(
                raw_maximum_face_speed * geometry.mesh_length_unit_m),
            max_recession_velocity_m_s=maximum_recession * geometry.mesh_length_unit_m,
            max_growth_velocity_m_s=maximum_growth * geometry.mesh_length_unit_m,
            max_surface_mechanism_growth_velocity_m_s=(
                float(np.max(surface_growth_velocity))
                if surface_growth_velocity.size else 0.0),
            max_displacement_mesh_units=displacement, cfl_substeps=int(substeps),
            cfl_number=float(cfl_number), reinitialized=bool(reinitialize),
            reinitialization_method=(reinitialization_method if reinitialize else None),
            topology_change_policy=str(topology_change_policy),
            topology_event=topology_event,
            ballistic_periodic_lateral=periodic_ballistic,
            profile_periodic_lateral=profile_periodic_lateral,
            periodic_seam_projection_max_mesh_units=periodic_seam_projection,
            periodic_seam_velocity_projection_max_mesh_units_s=(
                periodic_seam_velocity_projection),
            periodic_seam_velocity_projection_displacement_mesh_units=(
                periodic_seam_velocity_projection * float(duration_s)),
            removed_unresolved_solid_cells=removed_unresolved_solid_cells,
            reassigned_unresolved_material_nodes=(
                reassigned_unresolved_material_nodes),
            unresolved_material_volume_upper_bound_m3=(
                reassigned_unresolved_material_nodes
                * (geometry.dx * geometry.mesh_length_unit_m) ** 3),
            filled_unresolved_gas_cavity_cells=(
                filled_unresolved_gas_cavity_cells),
            unresolved_gas_cavity_volume_upper_bound_m3=(
                filled_unresolved_gas_cavity_cells
                * (geometry.dx * geometry.mesh_length_unit_m) ** 3),
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
        ballistic_transport="forward", ballistic_periodic_lateral=None,
        ballistic_face_quadrature_points=1, reinitialization_method="skfmm",
        topology_change_policy="refuse",
        surface_state_remap_backend="legacy_knn",
        adaptive_timestep_options=None):
    """Run verified feature steps with conserved surface state and optional quasi-static charging.

    A fixed ``charging_poisson_system`` is valid for one geometry only. Repeated charged evolution
    instead requires ``charging_system_builder(geometry)`` to rebuild the physical material operator
    after every interface update. Each geometry is independently converged from zero stored charge;
    this is the quasi-static charging limit, not a claim that transient surface charge was remapped.

    ``adaptive_timestep_options`` controls the outer transport/chemistry/profile coupling step, not
    merely the internal level-set CFL substeps. A trial is rejected without changing state when its
    maximum interface displacement exceeds the declared cell fraction or when a smaller step may
    resolve a topology/remap refusal. The identical state and sampling seed are replayed after
    shrinking ``dt``; accepted steps may then grow toward the displacement target. This preserves the
    physical operator while preventing a stable advection kernel from moving under stale fluxes.
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
    adaptive = (
        None if adaptive_timestep_options is None
        else dict(adaptive_timestep_options))
    if charging_poisson_system is not None and (
            int(n_steps) > 1 or adaptive is not None):
        raise ValueError(
            "multi-step charged profile evolution requires a geometry-dependent Poisson builder")
    if charging_system_builder is not None and not callable(charging_system_builder):
        raise TypeError("charging_system_builder must be callable")
    nominal_step_duration = float(duration_s) / int(n_steps)
    if adaptive is not None:
        allowed = {
            "initial_step_duration_s", "minimum_step_duration_s",
            "maximum_step_duration_s", "target_displacement_cells",
            "maximum_displacement_cells", "shrink_factor", "growth_factor",
            "safety_factor", "maximum_retries_per_step",
            "maximum_accepted_steps",
        }
        unknown = set(adaptive) - allowed
        if unknown:
            raise ValueError(
                "unknown adaptive profile-timestep options: "
                + ", ".join(sorted(unknown)))
        if duration_s <= 0.0:
            raise ValueError("adaptive profile stepping requires positive duration_s")
        initial_step = float(adaptive.get(
            "initial_step_duration_s", nominal_step_duration))
        minimum_step = float(adaptive.get(
            "minimum_step_duration_s", max(
                np.finfo(float).eps * duration_s,
                nominal_step_duration / 128.0)))
        maximum_step = float(adaptive.get(
            "maximum_step_duration_s", nominal_step_duration))
        target_cells = float(adaptive.get("target_displacement_cells", 0.35))
        maximum_cells = float(adaptive.get("maximum_displacement_cells", 0.75))
        shrink_factor = float(adaptive.get("shrink_factor", 0.5))
        growth_factor = float(adaptive.get("growth_factor", 1.5))
        safety_factor = float(adaptive.get("safety_factor", 0.9))
        maximum_retries = int(adaptive.get("maximum_retries_per_step", 20))
        maximum_accepted_steps = int(adaptive.get(
            "maximum_accepted_steps", max(1000, 20 * int(n_steps))))
        numeric = np.asarray([
            initial_step, minimum_step, maximum_step, target_cells,
            maximum_cells, shrink_factor, growth_factor, safety_factor], dtype=float)
        if (np.any(~np.isfinite(numeric)) or np.any(numeric[:5] <= 0.0)
                or minimum_step > initial_step or initial_step > maximum_step
                or target_cells >= maximum_cells
                or not 0.0 < shrink_factor < 1.0 or growth_factor <= 1.0
                or not 0.0 < safety_factor <= 1.0
                or maximum_retries <= 0 or maximum_accepted_steps <= 0):
            raise ValueError("invalid adaptive profile-timestep controls")
        controller = {
            "initial_step_duration_s": initial_step,
            "minimum_step_duration_s": minimum_step,
            "maximum_step_duration_s": maximum_step,
            "target_displacement_cells": target_cells,
            "maximum_displacement_cells": maximum_cells,
            "shrink_factor": shrink_factor,
            "growth_factor": growth_factor,
            "safety_factor": safety_factor,
            "maximum_retries_per_step": maximum_retries,
            "maximum_accepted_steps": maximum_accepted_steps,
        }
    else:
        controller = None
    current_geometry = geometry; current_state = None; current_fingerprint = None
    results = []; physical_time_s = 0.0; step_index = 0
    next_step_duration = (
        nominal_step_duration if controller is None
        else controller["initial_step_duration_s"])
    while (
            step_index < int(n_steps) if controller is None
            else physical_time_s < float(duration_s)):
        if controller is not None and step_index >= controller["maximum_accepted_steps"]:
            raise RuntimeError(
                "adaptive profile stepping exhausted maximum_accepted_steps at "
                f"t={physical_time_s:.8g}/{float(duration_s):.8g} s")
        step_duration = (
            nominal_step_duration if controller is None
            else min(next_step_duration, float(duration_s) - physical_time_s))
        rejected_trials = []
        retry_count = 0
        while True:
            if controller is not None and retry_count > controller["maximum_retries_per_step"]:
                raise RuntimeError(
                    "adaptive profile stepping exhausted its retry budget at "
                    f"t={physical_time_s:.8g} s")
            step_poisson_system = charging_poisson_system
            step_initial_charge = initial_charge_node_c
            if charging_system_builder is not None:
                step_poisson_system = charging_system_builder(current_geometry)
                if not isinstance(step_poisson_system, NodalPoissonSystem3D):
                    raise TypeError("charging_system_builder must return NodalPoissonSystem3D")
                if step_poisson_system.shape != current_geometry.phi.shape:
                    raise ValueError("rebuilt Poisson nodal grid must match the feature geometry")
                # A previous nodal charge grid is not a conservative representation on the moved
                # surface. In the quasi-static limit each geometry owns an independently converged
                # root. Adaptive retries remain on the same geometry and replay this same input.
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
                    ballistic_periodic_lateral=ballistic_periodic_lateral,
                    ballistic_face_quadrature_points=ballistic_face_quadrature_points,
                    reinitialization_method=reinitialization_method,
                    topology_change_policy=topology_change_policy,
                    surface_state_remap_backend=surface_state_remap_backend,
                    cfl_number=cfl_number, reinitialize=reinitialize,
                    transport_device=transport_device)
            except (ValueError, RuntimeError) as error:
                message = str(error)
                retryable = (
                    message.startswith("surface topology changed under ")
                    or message.startswith("surface remap distance ")
                    or message.startswith("material surface appeared or disappeared")
                    or message.startswith("surface contraction exceeds bounded coverage capacity"))
                proposed = step_duration * (
                    controller["shrink_factor"] if controller is not None else 1.0)
                at_minimum = (
                    controller is not None
                    and step_duration <= controller["minimum_step_duration_s"])
                if (controller is None or not retryable
                        or at_minimum):
                    denominator = (
                        int(n_steps) if controller is None
                        else controller["maximum_accepted_steps"])
                    contextual_message = (
                        f"feature step {step_index + 1}/{denominator} at "
                        f"t={physical_time_s:.8g} s, dt={step_duration:.8g} s: "
                        f"{error}")
                    if isinstance(error, SurfaceTopologyChangeError):
                        raise SurfaceTopologyChangeError(
                            contextual_message, method=error.method,
                            old_topology=error.old_topology,
                            new_topology=error.new_topology,
                            old_mesh_topology=error.old_mesh_topology,
                            new_mesh_topology=error.new_mesh_topology,
                            changed_slice_topology=error.changed_slice_topology) from error
                    raise type(error)(
                        contextual_message) from error
                rejected_trials.append({
                    "duration_s": float(step_duration),
                    "reason": message,
                    "classification": "inline_recovery_retry",
                })
                step_duration = max(
                    controller["minimum_step_duration_s"], proposed)
                retry_count += 1
                continue
            if controller is not None:
                displacement = float(
                    result.diagnostics["max_displacement_mesh_units"])
                limit = controller["maximum_displacement_cells"] * current_geometry.dx
                if displacement > limit:
                    proposed = (
                        step_duration * controller["safety_factor"] * limit
                        / displacement)
                    at_minimum = (
                        step_duration <= controller["minimum_step_duration_s"])
                    if at_minimum:
                        raise RuntimeError(
                            "adaptive profile coupling displacement remains unresolved at "
                            f"minimum dt={controller['minimum_step_duration_s']:.8g} s: "
                            f"displacement={displacement:.8g}, limit={limit:.8g} mesh units")
                    rejected_trials.append({
                        "duration_s": float(step_duration),
                        "reason": (
                            f"coupling displacement {displacement:.8g} exceeds "
                            f"{limit:.8g} mesh units"),
                        "classification": "inline_recovery_retry",
                    })
                    step_duration = max(
                        controller["minimum_step_duration_s"],
                        min(step_duration * controller["shrink_factor"], proposed))
                    retry_count += 1
                    continue
            break

        physical_time_s += step_duration
        if controller is not None:
            diagnostics = dict(
                result.diagnostics,
                adaptive_profile_timestep=True,
                accepted_step_duration_s=float(step_duration),
                accepted_physical_time_s=float(physical_time_s),
                adaptive_rejected_trials=tuple(rejected_trials),
                adaptive_retry_count=int(retry_count),
                adaptive_controller=dict(controller),
            )
            result = replace(result, diagnostics=diagnostics)
        results.append(result)
        current_geometry = result.geometry
        current_state = result.next_surface_state
        current_fingerprint = result.next_surface_state_mesh_fingerprint
        step_index += 1
        if controller is not None:
            displacement = float(result.diagnostics["max_displacement_mesh_units"])
            target = controller["target_displacement_cells"] * current_geometry.dx
            if displacement > 0.0:
                factor = controller["safety_factor"] * target / displacement
                factor = float(np.clip(
                    factor, controller["shrink_factor"], controller["growth_factor"]))
            else:
                factor = controller["growth_factor"]
            next_step_duration = float(np.clip(
                step_duration * factor,
                controller["minimum_step_duration_s"],
                controller["maximum_step_duration_s"]))
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
