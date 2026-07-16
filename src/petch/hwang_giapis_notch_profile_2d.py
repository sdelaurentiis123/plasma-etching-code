"""Source-faithful local notch refinement for Hwang--Giapis/Nozawa.

Hwang and Giapis, JVST B 15, 70 (1997), Sec. IV, use a coarse steady
charging solution outside the notch and a 5 nm square-cell local solve inside
it.  The local boundary is the original poly-Si sidewall potential ``Vs(k)``,
the moving poly-Si boundary is the fixed conductor voltage ``Vp``, and the
newly exposed SiO2 potential follows their Eq. (4.4)::

    V(i) = Vp + N(i) / N(1) * (V0 - Vp)
    V0 = (Vs(0) + Vp) / 2.

This module implements that bounded refinement layer.  It deliberately reuses
the common-engine Hwang--Giapis poly-Si yield and SiO2 scattering laws.  Its
profile state is a 2-D square-cell inventory extruded by one cell only for the
shared exact float64 triangle-hit tracer; it is not a second charging model.

The declared profile closure follows the paper:

* charge/field outside the local notch remains the supplied steady state;
* electrons inside the small shadowed notch are neglected;
* the exposed-oxide potential is proportional to cumulative ion count;
* one poly-Si cell is removed after 50 expected reactive collisions;
* poly-Si fragments with no face-connected path to the bulk line detach;
* direct single-bounce neutralized SiO2 scattering is retained;
* poly-Si scattering and spontaneous neutral-Cl etching are omitted.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from .boundary_transport_3d import (
    _first_segment_triangle_hit_float64_3d,
    _trace_field_events_float64_3d,
)
from .chlorine_poly_si import HwangGiapisClSiYield
from .hwang_giapis_scatter_3d import HwangGiapisSiO2ForwardScatter3D


HWANG_GIAPIS_LOCAL_NOTCH_SCHEMA = "petch-hwang-giapis-local-notch-v3"
HWANG_GIAPIS_LOCAL_CELL_UM = 0.005
HWANG_GIAPIS_POLY_THICKNESS_UM = 0.3
HWANG_GIAPIS_LINE_WIDTH_UM = 0.5
HWANG_GIAPIS_REACTIVE_COLLISIONS_PER_CELL = 50.0
HWANG_GIAPIS_FIG13_LAUNCHED_IONS_PER_0P5UM = 18.7e6


def hwang_giapis_exposed_oxide_potential_v(
        cumulative_ion_count, *, poly_potential_v,
        initial_sidewall_floor_potential_v):
    """Evaluate Hwang--Giapis Eq. (4.4) on exposed SiO2 cells.

    Before the reference cell has received an ion, the newly uncovered oxide
    inherits ``Vp`` continuously from the poly-Si/oxide interface.  Once
    ``N(1)`` is positive, the printed equation is applied without clipping.
    """
    count = np.asarray(cumulative_ion_count, dtype=float)
    vp = float(poly_potential_v)
    vs0 = float(initial_sidewall_floor_potential_v)
    if (count.ndim != 1 or not count.size or np.any(~np.isfinite(count))
            or np.any(count < 0.0) or not np.isfinite(vp)
            or not np.isfinite(vs0)):
        raise ValueError("Eq. (4.4) requires finite nonnegative ion counts and potentials")
    v0 = 0.5 * (vs0 + vp)
    if count[0] == 0.0:
        return np.full(count.shape, vp)
    return vp + (count / count[0]) * (v0 - vp)


@dataclass(frozen=True)
class HwangGiapisLocalIonEntry2D:
    """Weighted ions entering the local notch through the original sidewall."""

    height_um: np.ndarray
    velocity_xz_sqrt_eV: np.ndarray
    expected_count: np.ndarray

    def __post_init__(self):
        height = np.asarray(self.height_um, dtype=float).copy()
        velocity = np.asarray(self.velocity_xz_sqrt_eV, dtype=float).copy()
        count = np.asarray(self.expected_count, dtype=float).copy()
        if (height.ndim != 1 or velocity.shape != (height.size, 2)
                or count.shape != height.shape or np.any(~np.isfinite(height))
                or np.any(~np.isfinite(velocity)) or np.any(~np.isfinite(count))
                or np.any(height < 0.0)
                or np.any(height > HWANG_GIAPIS_POLY_THICKNESS_UM)
                or np.any(count < 0.0)
                or np.any(np.linalg.norm(velocity, axis=1) <= 0.0)
                or np.any(velocity[:, 0] <= 0.0)):
            raise ValueError(
                "local ion entries need finite in-range heights, positive inward "
                "velocity, and nonnegative expected counts")
        for value in (height, velocity, count):
            value.setflags(write=False)
        object.__setattr__(self, "height_um", height)
        object.__setattr__(self, "velocity_xz_sqrt_eV", velocity)
        object.__setattr__(self, "expected_count", count)

    @property
    def launched_count(self):
        return float(np.sum(self.expected_count))


@dataclass(frozen=True)
class HwangGiapisLocalBoundary2D:
    """Global charged-transport data reduced onto the 5 nm local notch patch."""

    entries: HwangGiapisLocalIonEntry2D
    sidewall_potential_v: np.ndarray
    poly_potential_v: float
    global_cell_size_um: float
    source_particle_count: int
    target_event_count: int
    published_launched_ions_per_0p5um: float
    provenance: MappingProxyType

    def __post_init__(self):
        sidewall = np.asarray(self.sidewall_potential_v, dtype=float).copy()
        values = np.asarray([
            self.poly_potential_v, self.global_cell_size_um,
            self.published_launched_ions_per_0p5um], dtype=float)
        if (not isinstance(self.entries, HwangGiapisLocalIonEntry2D)
                or sidewall.shape != (
                    int(round(HWANG_GIAPIS_POLY_THICKNESS_UM
                              / HWANG_GIAPIS_LOCAL_CELL_UM)),)
                or np.any(~np.isfinite(sidewall))
                or np.any(~np.isfinite(values))
                or self.global_cell_size_um <= 0.0
                or self.published_launched_ions_per_0p5um <= 0.0
                or int(self.source_particle_count) != self.source_particle_count
                or self.source_particle_count <= 0
                or int(self.target_event_count) != self.target_event_count
                or self.target_event_count <= 0
                or self.target_event_count != len(self.entries.height_um)):
            raise ValueError("invalid global-to-local Hwang--Giapis boundary reduction")
        sidewall.setflags(write=False)
        object.__setattr__(self, "sidewall_potential_v", sidewall)
        object.__setattr__(self, "poly_potential_v", float(self.poly_potential_v))
        object.__setattr__(self, "global_cell_size_um", float(self.global_cell_size_um))
        object.__setattr__(self, "source_particle_count", int(self.source_particle_count))
        object.__setattr__(self, "target_event_count", int(self.target_event_count))
        object.__setattr__(
            self, "published_launched_ions_per_0p5um",
            float(self.published_launched_ions_per_0p5um))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def hwang_giapis_local_boundary_from_edge_array_result(
        result, *,
        launched_ions_per_0p5um=HWANG_GIAPIS_FIG13_LAUNCHED_IONS_PER_0P5UM):
    """Reduce a charged 2-D edge-array solution to the Fig. 12 local inputs.

    The edge-array solver samples the entire lower boundary with equal-weight
    ions.  The paper reports the experimental 200% overetch fluence as
    ``18.7e6 ions / 0.5 um`` of lower-boundary width.  Therefore each retained
    target-wall event represents

    ``18.7e6 * source_width / 0.5 / source_particle_count``

    physical ions.  Coordinates are rotated into the local convention:
    ``x`` points into the edge poly-Si line and ``z`` points upward from SiO2.
    No profile datum is used in this reduction.
    """
    if not isinstance(result, dict):
        raise TypeError("edge-array result must be a dictionary")
    lineage = result.get("final_ion_lineage")
    geometry = result.get("geom")
    if not isinstance(lineage, dict) or not isinstance(geometry, dict):
        raise ValueError(
            "edge-array result needs return_final_ion_lineage=True and geometry metadata")
    required_lineage = {
        "hit_type", "hit_z_grid", "hit_vx_sqrt_eV", "hit_vz_sqrt_eV",
        "impact_energy_eV", "source_particle_count", "cell_size_um",
        "source_width_um", "edge_inner_poly_hit_type"}
    required_geometry = {"nz", "trench0", "poly_cells"}
    if not required_lineage.issubset(lineage) or not required_geometry.issubset(geometry):
        raise ValueError("edge-array lineage is missing required local-reduction fields")
    hit_type = np.asarray(lineage["hit_type"], dtype=int)
    hit_z = np.asarray(lineage["hit_z_grid"], dtype=float)
    hit_vx = np.asarray(lineage["hit_vx_sqrt_eV"], dtype=float)
    hit_vz = np.asarray(lineage["hit_vz_sqrt_eV"], dtype=float)
    impact_energy = np.asarray(lineage["impact_energy_eV"], dtype=float)
    if (hit_type.ndim != 1 or any(
            item.shape != hit_type.shape
            for item in (hit_z, hit_vx, hit_vz, impact_energy))):
        raise ValueError("edge-array final ion lineage has inconsistent shapes")
    selected = hit_type == int(lineage["edge_inner_poly_hit_type"])
    if not np.any(selected):
        raise ValueError("charged edge-array audit delivered no ions to the target inner wall")
    cell_size = float(lineage["cell_size_um"])
    source_count = int(lineage["source_particle_count"])
    source_width = float(lineage["source_width_um"])
    launched = float(launched_ions_per_0p5um)
    if (not np.isfinite(cell_size) or cell_size <= 0.0
            or source_count != len(hit_type) or source_count <= 0
            or not np.isfinite(source_width) or source_width <= 0.0
            or not np.isfinite(launched) or launched <= 0.0):
        raise ValueError("invalid edge-array source measure")
    nz = int(geometry["nz"])
    height = (float(nz - 1) - hit_z[selected]) * cell_size
    velocity = np.column_stack((-hit_vx[selected], -hit_vz[selected]))
    energy_replay = np.einsum("rc,rc->r", velocity, velocity)
    if (np.any(height < -1e-10)
            or np.any(height > HWANG_GIAPIS_POLY_THICKNESS_UM + cell_size)
            or np.any(velocity[:, 0] <= 0.0)
            or not np.allclose(
                energy_replay, impact_energy[selected], rtol=2e-12, atol=2e-12)):
        raise ValueError(
            "target-wall lineage is inconsistent with the local coordinate/energy convention")
    height = np.clip(height, 0.0, HWANG_GIAPIS_POLY_THICKNESS_UM)
    physical_count_per_event = (
        launched * source_width / HWANG_GIAPIS_LINE_WIDTH_UM / source_count)
    entries = HwangGiapisLocalIonEntry2D(
        height, velocity,
        np.full(np.count_nonzero(selected), physical_count_per_event))

    potential = np.asarray(result.get("V"), dtype=float)
    if potential.ndim != 2 or potential.shape[1] != nz:
        raise ValueError("edge-array result does not contain its final potential map")
    gas_x = int(geometry["trench0"])
    if not 0 <= gas_x < potential.shape[0]:
        raise ValueError("edge-array target gas boundary lies outside the potential map")
    global_height = (float(nz - 1) - np.arange(nz)) * cell_size
    order = np.argsort(global_height)
    local_height = (
        np.arange(int(round(
            HWANG_GIAPIS_POLY_THICKNESS_UM / HWANG_GIAPIS_LOCAL_CELL_UM)))
        + 0.5) * HWANG_GIAPIS_LOCAL_CELL_UM
    sidewall = np.interp(
        local_height, global_height[order], potential[gas_x, order])
    poly_potential = float(result["V_poly_edge"])
    provenance = {
        "schema": HWANG_GIAPIS_LOCAL_NOTCH_SCHEMA,
        "global_solver": "petch.charging2d.solve_edge_array_charging",
        "global_resolution_um": cell_size,
        "target_surface": "inner wall of the edge line bordering the ordinary trench",
        "source_fluence": (
            "Hwang & Giapis Fig. 13: 18.7e6 ions per 0.5 um lower-boundary width"),
        "reduction": (
            "equal-weight final ion lineage; 2-D coordinate rotation; gas-adjacent "
            "global potential interpolated to 5 nm cell centers"),
        "target_event_fraction": float(np.mean(selected)),
    }
    return HwangGiapisLocalBoundary2D(
        entries, sidewall, poly_potential, cell_size, source_count,
        int(np.count_nonzero(selected)), launched, MappingProxyType(provenance))


@dataclass(frozen=True)
class HwangGiapisLocalPotential2D:
    """Cell-centred local Laplace solution sampled on a regular tracer grid."""

    potential_v: np.ndarray
    origin_um: tuple[float, float, float]
    spacing_um: tuple[float, float, float]
    gas_cell_potential_v: np.ndarray
    oxide_potential_v: np.ndarray

    def __post_init__(self):
        potential = np.asarray(self.potential_v, dtype=float).copy()
        gas = np.asarray(self.gas_cell_potential_v, dtype=float).copy()
        oxide = np.asarray(self.oxide_potential_v, dtype=float).copy()
        origin = tuple(float(item) for item in self.origin_um)
        spacing = tuple(float(item) for item in self.spacing_um)
        if (potential.ndim != 3 or potential.shape[1] != 2
                or gas.ndim != 2 or oxide.shape != (gas.shape[0],)
                or len(origin) != 3 or len(spacing) != 3
                or np.any(~np.isfinite(potential)) or np.any(~np.isfinite(gas))
                or np.any(~np.isfinite(oxide))
                or any(not np.isfinite(item) for item in (*origin, *spacing))
                or any(item <= 0.0 for item in spacing)):
            raise ValueError("invalid local notch potential")
        for value in (potential, gas, oxide):
            value.setflags(write=False)
        object.__setattr__(self, "potential_v", potential)
        object.__setattr__(self, "gas_cell_potential_v", gas)
        object.__setattr__(self, "oxide_potential_v", oxide)
        object.__setattr__(self, "origin_um", origin)
        object.__setattr__(self, "spacing_um", spacing)


def solve_hwang_giapis_local_laplace_2d(
        poly_cell, sidewall_potential_v, cumulative_oxide_ion_count, *,
        poly_potential_v, cell_size_um=HWANG_GIAPIS_LOCAL_CELL_UM):
    """Solve the Fig. 12 crosshatched local region with a five-point stencil."""
    solid = np.asarray(poly_cell, dtype=bool)
    sidewall = np.asarray(sidewall_potential_v, dtype=float)
    count = np.asarray(cumulative_oxide_ion_count, dtype=float)
    vp = float(poly_potential_v)
    dx = float(cell_size_um)
    if (solid.ndim != 2 or min(solid.shape) < 2
            or sidewall.shape != (solid.shape[1],)
            or count.shape != (solid.shape[0],)
            or np.any(~np.isfinite(sidewall)) or np.any(~np.isfinite(count))
            or np.any(count < 0.0) or not np.isfinite(vp)
            or not np.isfinite(dx) or dx <= 0.0):
        raise ValueError("invalid Hwang--Giapis local Laplace inputs")
    gas = ~solid
    if not np.any(gas):
        raise ValueError("the local Laplace domain needs at least one removed poly-Si cell")
    if np.any(gas[-1]):
        raise ValueError("the local notch reached the right refinement boundary")
    oxide = hwang_giapis_exposed_oxide_potential_v(
        count, poly_potential_v=vp,
        initial_sidewall_floor_potential_v=sidewall[0])
    unknown = -np.ones(solid.shape, dtype=int)
    unknown[gas] = np.arange(np.count_nonzero(gas))
    rows = []
    columns = []
    values = []
    rhs = np.zeros(np.count_nonzero(gas))
    nx, nz = solid.shape
    for x, z in zip(*np.where(gas)):
        row = unknown[x, z]
        rows.append(row); columns.append(row); values.append(4.0)
        for dx_index, dz_index in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            neighbor_x = x + dx_index
            neighbor_z = z + dz_index
            if neighbor_x < 0:
                rhs[row] += sidewall[z]
            elif neighbor_z < 0:
                rhs[row] += oxide[x]
            elif neighbor_x >= nx or neighbor_z >= nz:
                rhs[row] += vp
            elif gas[neighbor_x, neighbor_z]:
                rows.append(row)
                columns.append(unknown[neighbor_x, neighbor_z])
                values.append(-1.0)
            else:
                rhs[row] += vp
    operator = coo_matrix(
        (values, (rows, columns)),
        shape=(len(rhs), len(rhs))).tocsc()
    solution = np.asarray(spsolve(operator, rhs), dtype=float)
    if solution.shape != rhs.shape or np.any(~np.isfinite(solution)):
        raise RuntimeError("local Laplace solve did not produce a finite field")
    cell_potential = np.full(solid.shape, vp)
    cell_potential[gas] = solution
    # The common exact trajectory tracer consumes nodal samples.  Use one
    # ghost-cell-centred layer on every x/z side so all exact cell faces lie
    # strictly inside the interpolation grid.
    tracer = np.full((nx + 2, 2, nz + 2), vp)
    tracer[1:-1, :, 1:-1] = cell_potential[:, None, :]
    tracer[0, :, 1:-1] = sidewall[None, :]
    tracer[1:-1, :, 0] = oxide[:, None]
    tracer[0, :, 0] = 0.5 * (sidewall[0] + oxide[0])
    return HwangGiapisLocalPotential2D(
        tracer, (-0.5 * dx, 0.0, -0.5 * dx), (dx, dx, dx),
        cell_potential, oxide)


@dataclass(frozen=True)
class _LocalSurfaceMesh2D:
    vertices_um: np.ndarray
    faces: np.ndarray
    gas_normal: np.ndarray
    material_id: np.ndarray
    cell_x: np.ndarray
    cell_z: np.ndarray


def _append_extruded_segment(
        vertices, faces, normals, materials, cell_x, cell_z, *,
        start_x, start_z, end_x, end_z, width_y, normal, material, x_index, z_index):
    base = len(vertices)
    vertices.extend((
        (start_x, 0.0, start_z), (end_x, 0.0, end_z),
        (end_x, width_y, end_z), (start_x, width_y, start_z)))
    candidate = ((base, base + 1, base + 2), (base, base + 2, base + 3))
    desired = np.asarray(normal, dtype=float)
    for triangle in candidate:
        triangle = list(triangle)
        points = np.asarray([vertices[index] for index in triangle])
        geometric = np.cross(points[1] - points[0], points[2] - points[0])
        if np.dot(geometric, desired) < 0.0:
            triangle[1], triangle[2] = triangle[2], triangle[1]
        faces.append(tuple(triangle))
        normals.append(tuple(desired))
        materials.append(int(material))
        cell_x.append(int(x_index))
        cell_z.append(int(z_index))


def _local_surface_mesh_2d(poly_cell, cell_size_um):
    solid = np.asarray(poly_cell, dtype=bool)
    nx, nz = solid.shape
    dx = float(cell_size_um)
    vertices = []
    faces = []
    normals = []
    materials = []
    cell_x = []
    cell_z = []
    for x, z in zip(*np.where(solid)):
        x0 = x * dx; x1 = (x + 1) * dx
        z0 = z * dx; z1 = (z + 1) * dx
        if x == 0 or not solid[x - 1, z]:
            _append_extruded_segment(
                vertices, faces, normals, materials, cell_x, cell_z,
                start_x=x0, start_z=z0, end_x=x0, end_z=z1,
                width_y=dx, normal=(-1.0, 0.0, 0.0), material=1,
                x_index=x, z_index=z)
        if x + 1 < nx and not solid[x + 1, z]:
            _append_extruded_segment(
                vertices, faces, normals, materials, cell_x, cell_z,
                start_x=x1, start_z=z1, end_x=x1, end_z=z0,
                width_y=dx, normal=(1.0, 0.0, 0.0), material=1,
                x_index=x, z_index=z)
        if z > 0 and not solid[x, z - 1]:
            _append_extruded_segment(
                vertices, faces, normals, materials, cell_x, cell_z,
                start_x=x1, start_z=z0, end_x=x0, end_z=z0,
                width_y=dx, normal=(0.0, 0.0, -1.0), material=1,
                x_index=x, z_index=z)
        if z + 1 < nz and not solid[x, z + 1]:
            _append_extruded_segment(
                vertices, faces, normals, materials, cell_x, cell_z,
                start_x=x0, start_z=z1, end_x=x1, end_z=z1,
                width_y=dx, normal=(0.0, 0.0, 1.0), material=1,
                x_index=x, z_index=z)
    # Every removed bottom cell exposes one SiO2 floor segment.
    for x in np.flatnonzero(~solid[:, 0]):
        _append_extruded_segment(
            vertices, faces, normals, materials, cell_x, cell_z,
            start_x=x * dx, start_z=0.0, end_x=(x + 1) * dx, end_z=0.0,
            width_y=dx, normal=(0.0, 0.0, 1.0), material=3,
            x_index=x, z_index=-1)
    # Undercutting the uppermost poly-Si cell exposes the photoresist floor
    # drawn in Hwang--Giapis Fig. 12.  Photoresist neither reacts nor scatters
    # in the declared model, but it must remain an absorbing geometric
    # boundary rather than becoming an artificial escape through the top of
    # the local refinement patch.
    for x in np.flatnonzero(~solid[:, -1]):
        _append_extruded_segment(
            vertices, faces, normals, materials, cell_x, cell_z,
            start_x=x * dx, start_z=nz * dx,
            end_x=(x + 1) * dx, end_z=nz * dx,
            width_y=dx, normal=(0.0, 0.0, -1.0), material=2,
            x_index=x, z_index=-2)
    vertex_array = np.asarray(vertices, dtype=float)
    face_array = np.asarray(faces, dtype=int)
    normal_array = np.asarray(normals, dtype=float)
    if (vertex_array.ndim != 2 or vertex_array.shape[1] != 3
            or face_array.ndim != 2 or face_array.shape[1] != 3
            or normal_array.shape != (len(face_array), 3)):
        raise RuntimeError("local surface construction produced an invalid mesh")
    return _LocalSurfaceMesh2D(
        vertex_array, face_array, normal_array,
        np.asarray(materials, dtype=int),
        np.asarray(cell_x, dtype=int), np.asarray(cell_z, dtype=int))


def _local_entry_origins_3d(height_um, cell_size_um, poly_thickness_um):
    """Place incoming ions just outside the original poly-Si sidewall.

    The local ``x=0`` plane is the original sidewall.  Intact cells must stop
    an incoming ion there, while an already opened height lets it enter the
    notch.  A positive-x launch would instead start intact-height ions inside
    solid poly-Si and allow nonphysical tunneling to an internal face.
    """
    height = np.asarray(height_um, dtype=float)
    dx = float(cell_size_um)
    thickness = float(poly_thickness_um)
    if (height.ndim != 1 or np.any(~np.isfinite(height))
            or np.any(height < 0.0) or np.any(height > thickness)
            or not np.isfinite(dx) or dx <= 0.0
            or not np.isfinite(thickness) or thickness <= 0.0):
        raise ValueError("invalid local-notch entry origins")
    return np.column_stack((
        np.full(height.size, -1e-7 * dx),
        np.full(height.size, 0.5 * dx),
        np.clip(height, 1e-7 * dx, thickness - 1e-7 * dx)))


@dataclass(frozen=True)
class HwangGiapisLocalNotchResult2D:
    poly_cell: np.ndarray
    reactive_collision_inventory: np.ndarray
    cumulative_oxide_ion_count: np.ndarray
    notch_depth_by_height_um: np.ndarray
    maximum_notch_depth_um: float
    removed_cell_count: int
    threshold_removed_cell_count: int
    detached_cell_count: int
    threshold_removed_reactive_collisions: float
    detached_reactive_collisions: float
    direct_reactive_collisions: float
    scattered_reactive_collisions: float
    launched_count: float
    landed_poly_count: float
    landed_oxide_count: float
    escaped_count: float
    batches: int
    provenance: MappingProxyType
    landed_photoresist_count: float = 0.0

    def __post_init__(self):
        poly = np.asarray(self.poly_cell, dtype=bool).copy()
        inventory = np.asarray(
            self.reactive_collision_inventory, dtype=float).copy()
        oxide = np.asarray(self.cumulative_oxide_ion_count, dtype=float).copy()
        depth = np.asarray(self.notch_depth_by_height_um, dtype=float).copy()
        values = np.asarray([
            self.maximum_notch_depth_um,
            self.threshold_removed_reactive_collisions,
            self.detached_reactive_collisions,
            self.direct_reactive_collisions,
            self.scattered_reactive_collisions, self.launched_count,
            self.landed_poly_count, self.landed_oxide_count,
            self.escaped_count, self.landed_photoresist_count], dtype=float)
        reactive_total = (
            self.direct_reactive_collisions
            + self.scattered_reactive_collisions)
        reactive_ledger_error = (
            self.threshold_removed_reactive_collisions
            + self.detached_reactive_collisions + float(np.sum(inventory))
            - reactive_total)
        particle_ledger_error = (
            self.landed_poly_count + self.landed_oxide_count
            + self.landed_photoresist_count + self.escaped_count
            - self.launched_count)
        if (poly.ndim != 2 or inventory.shape != poly.shape
                or oxide.shape != (poly.shape[0],)
                or depth.shape != (poly.shape[1],)
                or np.any(~np.isfinite(inventory)) or np.any(inventory < 0.0)
                or np.any(~np.isfinite(oxide)) or np.any(oxide < 0.0)
                or np.any(~np.isfinite(depth)) or np.any(depth < 0.0)
                or np.any(~np.isfinite(values)) or np.any(values < 0.0)
                or int(self.removed_cell_count) != self.removed_cell_count
                or self.removed_cell_count < 0
                or int(self.threshold_removed_cell_count)
                != self.threshold_removed_cell_count
                or self.threshold_removed_cell_count < 0
                or int(self.detached_cell_count) != self.detached_cell_count
                or self.detached_cell_count < 0
                or self.removed_cell_count
                != self.threshold_removed_cell_count + self.detached_cell_count
                or abs(reactive_ledger_error)
                > 1e-10 * max(reactive_total, 1.0)
                or abs(particle_ledger_error)
                > 1e-10 * max(self.launched_count, 1.0)
                or np.any(_detached_poly_cells(poly))
                or int(self.batches) != self.batches or self.batches <= 0):
            raise ValueError("invalid local notch result")
        for value in (poly, inventory, oxide, depth):
            value.setflags(write=False)
        object.__setattr__(self, "poly_cell", poly)
        object.__setattr__(self, "reactive_collision_inventory", inventory)
        object.__setattr__(self, "cumulative_oxide_ion_count", oxide)
        object.__setattr__(self, "notch_depth_by_height_um", depth)
        object.__setattr__(self, "removed_cell_count", int(self.removed_cell_count))
        object.__setattr__(
            self, "threshold_removed_cell_count",
            int(self.threshold_removed_cell_count))
        object.__setattr__(
            self, "detached_cell_count", int(self.detached_cell_count))
        object.__setattr__(self, "batches", int(self.batches))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


@dataclass(frozen=True)
class HwangGiapisLocalNotchCheckpoint2D:
    """Exact resumable state of the event-driven local profile."""

    poly_cell: np.ndarray
    reactive_collision_inventory: np.ndarray
    cumulative_oxide_ion_count: np.ndarray
    remaining_campaign_fraction: float
    direct_reactive_collisions: float
    scattered_reactive_collisions: float
    landed_poly_count: float
    landed_oxide_count: float
    escaped_count: float
    threshold_removed_cell_count: int
    detached_cell_count: int
    threshold_removed_reactive_collisions: float
    detached_reactive_collisions: float
    front_events: int
    landed_photoresist_count: float = 0.0

    def __post_init__(self):
        poly = np.asarray(self.poly_cell, dtype=bool).copy()
        inventory = np.asarray(
            self.reactive_collision_inventory, dtype=float).copy()
        oxide = np.asarray(self.cumulative_oxide_ion_count, dtype=float).copy()
        values = np.asarray([
            self.remaining_campaign_fraction,
            self.direct_reactive_collisions,
            self.scattered_reactive_collisions,
            self.landed_poly_count, self.landed_oxide_count,
            self.escaped_count, self.landed_photoresist_count,
            self.threshold_removed_reactive_collisions,
            self.detached_reactive_collisions], dtype=float)
        reactive_total = (
            self.direct_reactive_collisions
            + self.scattered_reactive_collisions)
        reactive_ledger_error = (
            self.threshold_removed_reactive_collisions
            + self.detached_reactive_collisions + float(np.sum(inventory))
            - reactive_total)
        if (poly.ndim != 2 or inventory.shape != poly.shape
                or oxide.shape != (poly.shape[0],)
                or np.any(~np.isfinite(inventory)) or np.any(inventory < 0.0)
                or np.any(~np.isfinite(oxide)) or np.any(oxide < 0.0)
                or np.any(~np.isfinite(values))
                or not 0.0 <= self.remaining_campaign_fraction <= 1.0
                or np.any(values[1:] < 0.0)
                or int(self.threshold_removed_cell_count)
                != self.threshold_removed_cell_count
                or self.threshold_removed_cell_count < 0
                or int(self.detached_cell_count) != self.detached_cell_count
                or self.detached_cell_count < 0
                or int(self.front_events) != self.front_events
                or self.front_events < 0
                or abs(reactive_ledger_error)
                > 1e-10 * max(reactive_total, 1.0)
                or np.any(_detached_poly_cells(poly))
                or np.any(inventory[~poly] != 0.0)):
            raise ValueError("invalid event-driven local-notch checkpoint")
        for value in (poly, inventory, oxide):
            value.setflags(write=False)
        object.__setattr__(self, "poly_cell", poly)
        object.__setattr__(self, "reactive_collision_inventory", inventory)
        object.__setattr__(self, "cumulative_oxide_ion_count", oxide)
        object.__setattr__(
            self, "threshold_removed_cell_count",
            int(self.threshold_removed_cell_count))
        object.__setattr__(
            self, "detached_cell_count", int(self.detached_cell_count))
        object.__setattr__(self, "front_events", int(self.front_events))


def _profile_depth(poly_cell, cell_size_um):
    solid = np.asarray(poly_cell, dtype=bool)
    depth = np.zeros(solid.shape[1])
    for z in range(solid.shape[1]):
        gas = np.flatnonzero(~solid[:, z])
        if gas.size:
            # Only the connected left-hand removal front is an undercut.
            run = 0
            for x in range(solid.shape[0]):
                if solid[x, z]:
                    break
                run += 1
            depth[z] = run * float(cell_size_um)
    return depth


def _detached_poly_cells(poly_cell):
    """Return cells with no finite-area connection to the bulk poly-Si line.

    The right and top edges of the local refinement patch are the unresolved
    continuation of the material line.  A retained cell must therefore have a
    four-neighbor path to at least one solid cell on either edge.  Diagonal
    corner contact has zero area in the square-cell representation and is not
    treated as a mechanical connection.

    Hwang and Giapis describe one 5 nm cell as a removable cluster of surface
    atoms.  Removing a threshold-crossing surface cell can disconnect a
    subthreshold cluster behind it; that cluster leaves as a detached fragment
    rather than remaining as a floating electrostatic obstacle.  Callers own
    the mass and collision-inventory ledgers for the returned mask.
    """
    solid = np.asarray(poly_cell, dtype=bool)
    if solid.ndim != 2 or not solid.size:
        raise ValueError("detached-fragment audit requires a 2-D material grid")
    connected = np.zeros_like(solid)
    stack = []
    nx, nz = solid.shape
    for z in range(nz):
        if solid[nx - 1, z]:
            connected[nx - 1, z] = True
            stack.append((nx - 1, z))
    for x in range(nx):
        if solid[x, nz - 1] and not connected[x, nz - 1]:
            connected[x, nz - 1] = True
            stack.append((x, nz - 1))
    while stack:
        x, z = stack.pop()
        for next_x, next_z in (
                (x - 1, z), (x + 1, z), (x, z - 1), (x, z + 1)):
            if (0 <= next_x < nx and 0 <= next_z < nz
                    and solid[next_x, next_z]
                    and not connected[next_x, next_z]):
                connected[next_x, next_z] = True
                stack.append((next_x, next_z))
    return solid & ~connected


def _five_cell_surface_reaction_normals(poly_cell, cell_size_um):
    """Reconstruct the Sec. IV C reaction normal from five surface cells.

    Hard staircase faces remain authoritative for collision detection.  Only
    the physical incidence angle used by the yield law is reconstructed from
    a least-squares fit through the impacted height and four neighboring
    heights, as Hwang and Giapis specify.  This keeps a 5 nm grid from
    quantizing every reaction angle to zero or ninety degrees.
    """
    solid = np.asarray(poly_cell, dtype=bool)
    dx = float(cell_size_um)
    if (solid.ndim != 2 or solid.shape[1] < 5
            or not np.isfinite(dx) or dx <= 0.0):
        raise ValueError("five-cell surface normals require a valid local profile")
    depth = _profile_depth(solid, dx)
    height = (np.arange(solid.shape[1], dtype=float) + 0.5) * dx
    normal = np.zeros((solid.shape[1], 3))
    for center in range(solid.shape[1]):
        start = min(max(center - 2, 0), solid.shape[1] - 5)
        selected = slice(start, start + 5)
        z = height[selected]
        x = depth[selected]
        z_centered = z - np.mean(z)
        denominator = float(np.dot(z_centered, z_centered))
        slope = (
            0.0 if denominator <= 0.0
            else float(np.dot(z_centered, x - np.mean(x)) / denominator))
        candidate = np.asarray((-1.0, 0.0, slope))
        normal[center] = candidate / np.linalg.norm(candidate)
    return normal


def evolve_hwang_giapis_local_notch_2d(
        entries: HwangGiapisLocalIonEntry2D, sidewall_potential_v, *,
        poly_potential_v, cell_size_um=HWANG_GIAPIS_LOCAL_CELL_UM,
        line_width_um=HWANG_GIAPIS_LINE_WIDTH_UM,
        poly_thickness_um=HWANG_GIAPIS_POLY_THICKNESS_UM,
        reactive_collisions_per_cell=HWANG_GIAPIS_REACTIVE_COLLISIONS_PER_CELL,
        batches=256, trajectory_fixed_dt=2.5e-4, trajectory_max_steps=8192,
        include_forward_scatter=True, include_exposed_oxide_charging=True):
    """Evolve one expected-value 2-D notch profile from a weighted entry measure."""
    if not isinstance(entries, HwangGiapisLocalIonEntry2D):
        raise TypeError("entries must be HwangGiapisLocalIonEntry2D")
    dx = float(cell_size_um)
    nx = int(round(float(line_width_um) / dx))
    nz = int(round(float(poly_thickness_um) / dx))
    sidewall = np.asarray(sidewall_potential_v, dtype=float)
    if (not np.isclose(nx * dx, line_width_um, atol=1e-12, rtol=0.0)
            or not np.isclose(nz * dx, poly_thickness_um, atol=1e-12, rtol=0.0)
            or sidewall.shape != (nz,) or np.any(~np.isfinite(sidewall))
            or not np.isfinite(poly_potential_v)
            or not np.isfinite(reactive_collisions_per_cell)
            or reactive_collisions_per_cell <= 0.0
            or int(batches) != batches or batches <= 0
            or not np.isfinite(trajectory_fixed_dt) or trajectory_fixed_dt <= 0.0
            or int(trajectory_max_steps) != trajectory_max_steps
            or trajectory_max_steps <= 0
            or not isinstance(include_forward_scatter, (bool, np.bool_))
            or not isinstance(
                include_exposed_oxide_charging, (bool, np.bool_))):
        raise ValueError("invalid local notch evolution controls")
    poly = np.ones((nx, nz), dtype=bool)
    inventory = np.zeros((nx, nz))
    oxide_count = np.zeros(nx)
    yield_law = HwangGiapisClSiYield()
    scatter = HwangGiapisSiO2ForwardScatter3D(3)
    event_count = entries.expected_count / int(batches)
    speed = np.linalg.norm(entries.velocity_xz_sqrt_eV, axis=1)
    direct_reactive = 0.0
    scattered_reactive = 0.0
    landed_poly = 0.0
    landed_oxide = 0.0
    landed_photoresist = 0.0
    escaped = 0.0
    threshold_removed_cells = 0
    detached_cells = 0
    threshold_removed_reactive = 0.0
    detached_reactive = 0.0

    for _batch in range(int(batches)):
        if not np.any(~poly):
            z_cell = np.minimum(
                np.floor(entries.height_um / dx).astype(int), nz - 1)
            cosine = entries.velocity_xz_sqrt_eV[:, 0] / speed
            reactive = event_count * yield_law.evaluate(speed * speed, cosine)
            np.add.at(inventory[0], z_cell, reactive)
            direct_reactive += float(np.sum(reactive))
            landed_poly += float(np.sum(event_count))
        else:
            local = solve_hwang_giapis_local_laplace_2d(
                poly, sidewall, (
                    oxide_count if include_exposed_oxide_charging
                    else np.zeros_like(oxide_count)),
                poly_potential_v=poly_potential_v, cell_size_um=dx)
            mesh = _local_surface_mesh_2d(poly, dx)
            reaction_normal = _five_cell_surface_reaction_normals(poly, dx)
            origin = _local_entry_origins_3d(
                entries.height_um, dx, nz * dx)
            velocity = np.column_stack((
                entries.velocity_xz_sqrt_eV[:, 0],
                np.zeros(len(event_count)),
                entries.velocity_xz_sqrt_eV[:, 1]))
            (hit_face, hit_cosine, hit_energy, termination,
             terminal_position, terminal_velocity) = _trace_field_events_float64_3d(
                origin, velocity, 1, local.potential_v,
                np.asarray(local.origin_um), np.asarray(local.spacing_um),
                mesh.vertices_um, mesh.faces, float(trajectory_fixed_dt),
                int(trajectory_max_steps), False)
            if np.any(termination == 0):
                raise RuntimeError(
                    "local charged trajectory horizon truncated a particle")
            hit = termination == 1
            escaped += float(np.sum(event_count[termination == 2]))
            for event in np.flatnonzero(hit):
                face = hit_face[event]
                material = mesh.material_id[face]
                weight = event_count[event]
                if material == 1:
                    incident_direction = (
                        terminal_velocity[event]
                        / np.linalg.norm(terminal_velocity[event]))
                    cosine = max(-float(np.dot(
                        incident_direction,
                        reaction_normal[mesh.cell_z[face]])), 0.0)
                    reactive = weight * float(yield_law.evaluate(
                        hit_energy[event], cosine))
                    inventory[mesh.cell_x[face], mesh.cell_z[face]] += reactive
                    direct_reactive += reactive
                    landed_poly += weight
                    continue
                if material == 2:
                    landed_photoresist += weight
                    continue
                if material != 3:
                    raise RuntimeError("local trajectory hit an unknown material")
                bottom_cell = mesh.cell_x[face]
                oxide_count[bottom_cell] += weight
                landed_oxide += weight
                if not include_forward_scatter:
                    continue
                probability = float(scatter.scattering_probability(hit_cosine[event]))
                if probability <= 0.0:
                    continue
                retention = float(scatter.energy_retention_fraction(hit_cosine[event]))
                if retention <= 0.0:
                    continue
                incident_direction = terminal_velocity[event] / np.linalg.norm(
                    terminal_velocity[event])
                normal = mesh.gas_normal[face]
                neutral_direction = (
                    incident_direction + 2.0 * hit_cosine[event] * normal)
                neutral_direction /= np.linalg.norm(neutral_direction)
                launch = terminal_position[event] + 1e-7 * dx * normal
                distance = 2.0 * np.hypot(nx * dx, nz * dx)
                neutral_face, _fraction, _impact = (
                    _first_segment_triangle_hit_float64_3d(
                        launch, distance * neutral_direction,
                        mesh.vertices_um, mesh.faces, np.ones(3), False))
                if neutral_face < 0 or mesh.material_id[neutral_face] != 1:
                    continue
                neutral_cosine = max(-float(np.dot(
                    neutral_direction,
                    reaction_normal[mesh.cell_z[neutral_face]])), 0.0)
                neutral_energy = retention * hit_energy[event]
                reactive = weight * probability * float(
                    yield_law.evaluate(neutral_energy, neutral_cosine))
                inventory[
                    mesh.cell_x[neutral_face], mesh.cell_z[neutral_face]] += reactive
                scattered_reactive += reactive
        removable = poly & (inventory >= reactive_collisions_per_cell)
        if np.any(removable):
            threshold_removed_cells += int(np.count_nonzero(removable))
            threshold_removed_reactive += float(np.sum(inventory[removable]))
            poly[removable] = False
            inventory[removable] = 0.0
            detached = _detached_poly_cells(poly)
            if np.any(detached):
                detached_cells += int(np.count_nonzero(detached))
                detached_reactive += float(np.sum(inventory[detached]))
                poly[detached] = False
                inventory[detached] = 0.0
        if not np.any(poly):
            raise RuntimeError("local notch removed the entire poly-Si refinement patch")
    depth = _profile_depth(poly, dx)
    removed = int(np.count_nonzero(~poly))
    provenance = {
        "schema": HWANG_GIAPIS_LOCAL_NOTCH_SCHEMA,
        "source": (
            "Hwang & Giapis, JVST B 15, 70 (1997), "
            "DOI 10.1116/1.589258, Sec. IV, Eqs. (4.1)--(4.4)"),
        "cell_size_um": dx,
        "grid_cells": (nx, nz),
        "reactive_collisions_per_cell": float(reactive_collisions_per_cell),
        "outside_field": "supplied steady global charging state",
        "local_field": "five-point Laplace solve with Eq. (4.4) oxide boundary",
        "entry_plane": (
            "gas-side one-sided limit of the original poly-Si sidewall; "
            "intact heights intercept at x=0 and opened heights enter the notch"),
        "reaction_angle": (
            "five-cell least-squares moving-surface normal, Sec. IV C"),
        "photoresist_ceiling": (
            "absorbing nonreactive boundary exposed by upper undercut; "
            "photoresist reaction and scattering are omitted per Sec. IV A"),
        "detached_fragment_policy": (
            "remove cells lacking a four-neighbor path to the right or top "
            "bulk-poly boundary; diagonal corner contact is zero-area"),
        "forward_scatter": bool(include_forward_scatter),
        "exposed_oxide_charging": bool(include_exposed_oxide_charging),
        "omissions": (
            "spontaneous neutral-Cl etching",
            "poly-Si inelastic scattering after a nonreactive collision",
            "multiple-bounce SiO2 scattering",
            "electrons inside the shadowed local notch",
        ),
    }
    return HwangGiapisLocalNotchResult2D(
        poly, inventory, oxide_count, depth, float(np.max(depth)), removed,
        threshold_removed_cells, detached_cells,
        threshold_removed_reactive, detached_reactive,
        direct_reactive, scattered_reactive, entries.launched_count,
        landed_poly, landed_oxide, escaped, int(batches), provenance,
        landed_photoresist_count=landed_photoresist)


def evolve_hwang_giapis_local_notch_event_driven_2d(
        entries: HwangGiapisLocalIonEntry2D, sidewall_potential_v, *,
        poly_potential_v, cell_size_um=HWANG_GIAPIS_LOCAL_CELL_UM,
        line_width_um=HWANG_GIAPIS_LINE_WIDTH_UM,
        poly_thickness_um=HWANG_GIAPIS_POLY_THICKNESS_UM,
        reactive_collisions_per_cell=HWANG_GIAPIS_REACTIVE_COLLISIONS_PER_CELL,
        trajectory_fixed_dt=2.5e-4, trajectory_max_steps=8192,
        include_forward_scatter=True, include_exposed_oxide_charging=True,
        maximum_front_events=10000,
        initial_checkpoint=None, progress_callback=None,
        checkpoint_callback=None):
    """Advance the local profile from one threshold-crossing event to the next.

    The input weights represent the complete published ion fluence.  At each
    geometry, the weighted trajectories define reactive-collision and oxide-hit
    rates per unit campaign fluence.  The integrator advances by exactly the
    fraction required for the next exposed poly-Si cell to reach 50 reactive
    collisions, removes that cell, and recomputes the field and trajectories.

    This is the source-faithful clock described in Sec. IV C: a new electric
    field is used as the notch boundary advances to the next cell.  It has no
    arbitrary equal-batch timestep and carries every subthreshold collision
    inventory forward.
    """
    if not isinstance(entries, HwangGiapisLocalIonEntry2D):
        raise TypeError("entries must be HwangGiapisLocalIonEntry2D")
    dx = float(cell_size_um)
    nx = int(round(float(line_width_um) / dx))
    nz = int(round(float(poly_thickness_um) / dx))
    sidewall = np.asarray(sidewall_potential_v, dtype=float)
    if (not np.isclose(nx * dx, line_width_um, atol=1e-12, rtol=0.0)
            or not np.isclose(nz * dx, poly_thickness_um, atol=1e-12, rtol=0.0)
            or sidewall.shape != (nz,) or np.any(~np.isfinite(sidewall))
            or not np.isfinite(poly_potential_v)
            or not np.isfinite(reactive_collisions_per_cell)
            or reactive_collisions_per_cell <= 0.0
            or not np.isfinite(trajectory_fixed_dt) or trajectory_fixed_dt <= 0.0
            or int(trajectory_max_steps) != trajectory_max_steps
            or trajectory_max_steps <= 0
            or not isinstance(include_forward_scatter, (bool, np.bool_))
            or not isinstance(
                include_exposed_oxide_charging, (bool, np.bool_))
            or int(maximum_front_events) != maximum_front_events
            or maximum_front_events <= 0
            or (initial_checkpoint is not None
                and not isinstance(
                    initial_checkpoint, HwangGiapisLocalNotchCheckpoint2D))
            or (progress_callback is not None and not callable(progress_callback))
            or (checkpoint_callback is not None
                and not callable(checkpoint_callback))):
        raise ValueError("invalid event-driven local notch evolution controls")
    if initial_checkpoint is None:
        poly = np.ones((nx, nz), dtype=bool)
        inventory = np.zeros((nx, nz))
        oxide_count = np.zeros(nx)
        remaining_fraction = 1.0
        direct_reactive = 0.0
        scattered_reactive = 0.0
        landed_poly = 0.0
        landed_oxide = 0.0
        landed_photoresist = 0.0
        escaped = 0.0
        threshold_removed_cells = 0
        detached_cells = 0
        threshold_removed_reactive = 0.0
        detached_reactive = 0.0
        front_events = 0
    else:
        if (initial_checkpoint.poly_cell.shape != (nx, nz)
                or initial_checkpoint.cumulative_oxide_ion_count.shape
                != (nx,)):
            raise ValueError(
                "event-driven checkpoint shape does not match the local grid")
        poly = initial_checkpoint.poly_cell.copy()
        inventory = initial_checkpoint.reactive_collision_inventory.copy()
        oxide_count = initial_checkpoint.cumulative_oxide_ion_count.copy()
        remaining_fraction = float(
            initial_checkpoint.remaining_campaign_fraction)
        direct_reactive = float(
            initial_checkpoint.direct_reactive_collisions)
        scattered_reactive = float(
            initial_checkpoint.scattered_reactive_collisions)
        landed_poly = float(initial_checkpoint.landed_poly_count)
        landed_oxide = float(initial_checkpoint.landed_oxide_count)
        landed_photoresist = float(
            initial_checkpoint.landed_photoresist_count)
        escaped = float(initial_checkpoint.escaped_count)
        threshold_removed_cells = int(
            initial_checkpoint.threshold_removed_cell_count)
        detached_cells = int(initial_checkpoint.detached_cell_count)
        threshold_removed_reactive = float(
            initial_checkpoint.threshold_removed_reactive_collisions)
        detached_reactive = float(
            initial_checkpoint.detached_reactive_collisions)
        front_events = int(initial_checkpoint.front_events)
    yield_law = HwangGiapisClSiYield()
    scatter = HwangGiapisSiO2ForwardScatter3D(3)
    campaign_weight = entries.expected_count
    speed = np.linalg.norm(entries.velocity_xz_sqrt_eV, axis=1)
    fraction_tolerance = 2e-14

    while remaining_fraction > fraction_tolerance:
        reactive_rate = np.zeros_like(inventory)
        oxide_rate = np.zeros_like(oxide_count)
        direct_rate = 0.0
        scattered_rate = 0.0
        landed_poly_rate = 0.0
        landed_oxide_rate = 0.0
        landed_photoresist_rate = 0.0
        escaped_rate = 0.0
        if not np.any(~poly):
            z_cell = np.minimum(
                np.floor(entries.height_um / dx).astype(int), nz - 1)
            cosine = entries.velocity_xz_sqrt_eV[:, 0] / speed
            reactive = campaign_weight * yield_law.evaluate(
                speed * speed, cosine)
            np.add.at(reactive_rate[0], z_cell, reactive)
            direct_rate = float(np.sum(reactive))
            landed_poly_rate = float(np.sum(campaign_weight))
        else:
            local = solve_hwang_giapis_local_laplace_2d(
                poly, sidewall, (
                    oxide_count if include_exposed_oxide_charging
                    else np.zeros_like(oxide_count)),
                poly_potential_v=poly_potential_v, cell_size_um=dx)
            mesh = _local_surface_mesh_2d(poly, dx)
            reaction_normal = _five_cell_surface_reaction_normals(poly, dx)
            origin = _local_entry_origins_3d(
                entries.height_um, dx, nz * dx)
            velocity = np.column_stack((
                entries.velocity_xz_sqrt_eV[:, 0],
                np.zeros(len(campaign_weight)),
                entries.velocity_xz_sqrt_eV[:, 1]))
            (hit_face, hit_cosine, hit_energy, termination,
             terminal_position, terminal_velocity) = (
                _trace_field_events_float64_3d(
                    origin, velocity, 1, local.potential_v,
                    np.asarray(local.origin_um), np.asarray(local.spacing_um),
                    mesh.vertices_um, mesh.faces, float(trajectory_fixed_dt),
                    int(trajectory_max_steps), False))
            if np.any(termination == 0):
                raise RuntimeError(
                    "local charged trajectory horizon truncated a particle")
            escaped_rate = float(np.sum(
                campaign_weight[termination == 2]))
            for event in np.flatnonzero(termination == 1):
                face = hit_face[event]
                material = mesh.material_id[face]
                weight = campaign_weight[event]
                if material == 1:
                    incident_direction = (
                        terminal_velocity[event]
                        / np.linalg.norm(terminal_velocity[event]))
                    cosine = max(-float(np.dot(
                        incident_direction,
                        reaction_normal[mesh.cell_z[face]])), 0.0)
                    reactive = weight * float(yield_law.evaluate(
                        hit_energy[event], cosine))
                    reactive_rate[
                        mesh.cell_x[face], mesh.cell_z[face]] += reactive
                    direct_rate += reactive
                    landed_poly_rate += weight
                    continue
                if material == 2:
                    landed_photoresist_rate += weight
                    continue
                if material != 3:
                    raise RuntimeError("local trajectory hit an unknown material")
                bottom_cell = mesh.cell_x[face]
                oxide_rate[bottom_cell] += weight
                landed_oxide_rate += weight
                if not include_forward_scatter:
                    continue
                probability = float(
                    scatter.scattering_probability(hit_cosine[event]))
                retention = float(
                    scatter.energy_retention_fraction(hit_cosine[event]))
                if probability <= 0.0 or retention <= 0.0:
                    continue
                incident_direction = terminal_velocity[event] / np.linalg.norm(
                    terminal_velocity[event])
                normal = mesh.gas_normal[face]
                neutral_direction = (
                    incident_direction + 2.0 * hit_cosine[event] * normal)
                neutral_direction /= np.linalg.norm(neutral_direction)
                launch = terminal_position[event] + 1e-7 * dx * normal
                distance = 2.0 * np.hypot(nx * dx, nz * dx)
                neutral_face, _fraction, _impact = (
                    _first_segment_triangle_hit_float64_3d(
                        launch, distance * neutral_direction,
                        mesh.vertices_um, mesh.faces, np.ones(3), False))
                if neutral_face < 0 or mesh.material_id[neutral_face] != 1:
                    continue
                neutral_cosine = max(-float(np.dot(
                    neutral_direction,
                    reaction_normal[mesh.cell_z[neutral_face]])), 0.0)
                neutral_energy = retention * hit_energy[event]
                reactive = weight * probability * float(
                    yield_law.evaluate(neutral_energy, neutral_cosine))
                reactive_rate[
                    mesh.cell_x[neutral_face],
                    mesh.cell_z[neutral_face]] += reactive
                scattered_rate += reactive

        active = poly & (reactive_rate > 0.0)
        if np.any(active):
            needed = np.maximum(
                reactive_collisions_per_cell - inventory[active], 0.0)
            fraction_to_event = float(np.min(
                needed / reactive_rate[active]))
            fraction_to_event = max(fraction_to_event, 0.0)
            advance = min(remaining_fraction, fraction_to_event)
        else:
            advance = remaining_fraction
        if advance <= fraction_tolerance and np.any(active):
            # Roundoff can leave the limiting cell an ulp below threshold.
            limiting = np.full(poly.shape, np.inf)
            limiting[active] = (
                reactive_collisions_per_cell - inventory[active]
            ) / reactive_rate[active]
            cell = np.unravel_index(np.argmin(limiting), limiting.shape)
            inventory[cell] = reactive_collisions_per_cell
            advance = 0.0
        inventory += advance * reactive_rate
        oxide_count += advance * oxide_rate
        direct_reactive += advance * direct_rate
        scattered_reactive += advance * scattered_rate
        landed_poly += advance * landed_poly_rate
        landed_oxide += advance * landed_oxide_rate
        landed_photoresist += advance * landed_photoresist_rate
        escaped += advance * escaped_rate
        remaining_fraction = max(remaining_fraction - advance, 0.0)

        removable = poly & (
            inventory >= reactive_collisions_per_cell * (1.0 - 2e-12))
        if np.any(removable):
            threshold_removed_cells += int(np.count_nonzero(removable))
            threshold_removed_reactive += float(np.sum(inventory[removable]))
            poly[removable] = False
            inventory[removable] = 0.0
            detached = _detached_poly_cells(poly)
            if np.any(detached):
                detached_cells += int(np.count_nonzero(detached))
                detached_reactive += float(np.sum(inventory[detached]))
                poly[detached] = False
                inventory[detached] = 0.0
            front_events += 1
            if front_events > int(maximum_front_events):
                raise RuntimeError(
                    "event-driven local notch exceeded its front-event budget")
            if np.any(~poly[-1]):
                raise RuntimeError(
                    "event-driven local notch reached the right refinement boundary")
        elif advance >= remaining_fraction + advance - fraction_tolerance:
            # The complete remaining campaign fluence was consumed without a
            # threshold crossing.
            break
        if progress_callback is not None:
            progress_callback(
                int(front_events), float(1.0 - remaining_fraction),
                float(np.max(_profile_depth(poly, dx))))
        if checkpoint_callback is not None:
            checkpoint_callback(HwangGiapisLocalNotchCheckpoint2D(
                poly, inventory, oxide_count, remaining_fraction,
                direct_reactive, scattered_reactive, landed_poly,
                landed_oxide, escaped, threshold_removed_cells,
                detached_cells, threshold_removed_reactive,
                detached_reactive, front_events,
                landed_photoresist_count=landed_photoresist))
        if not np.any(poly):
            raise RuntimeError(
                "event-driven local notch removed the entire poly-Si patch")

    depth = _profile_depth(poly, dx)
    removed = int(np.count_nonzero(~poly))
    if checkpoint_callback is not None:
        checkpoint_callback(HwangGiapisLocalNotchCheckpoint2D(
            poly, inventory, oxide_count, remaining_fraction,
            direct_reactive, scattered_reactive, landed_poly,
            landed_oxide, escaped, threshold_removed_cells,
            detached_cells, threshold_removed_reactive,
            detached_reactive, front_events,
            landed_photoresist_count=landed_photoresist))
    provenance = {
        "schema": HWANG_GIAPIS_LOCAL_NOTCH_SCHEMA,
        "source": (
            "Hwang & Giapis, JVST B 15, 70 (1997), "
            "DOI 10.1116/1.589258, Sec. IV, Eqs. (4.1)--(4.4)"),
        "cell_size_um": dx,
        "grid_cells": (nx, nz),
        "reactive_collisions_per_cell": float(
            reactive_collisions_per_cell),
        "outside_field": "supplied steady global charging state",
        "local_field": "five-point Laplace solve with Eq. (4.4) oxide boundary",
        "entry_plane": (
            "gas-side one-sided limit of the original poly-Si sidewall; "
            "intact heights intercept at x=0 and opened heights enter the notch"),
        "reaction_angle": (
            "five-cell least-squares moving-surface normal, Sec. IV C"),
        "photoresist_ceiling": (
            "absorbing nonreactive boundary exposed by upper undercut; "
            "photoresist reaction and scattering are omitted per Sec. IV A"),
        "detached_fragment_policy": (
            "remove cells lacking a four-neighbor path to the right or top "
            "bulk-poly boundary; diagonal corner contact is zero-area"),
        "front_integrator": (
            "event-driven campaign-fluence advance to the next 50-collision "
            "surface-cell threshold; no equal-batch timestep"),
        "front_events": int(front_events),
        "forward_scatter": bool(include_forward_scatter),
        "exposed_oxide_charging": bool(include_exposed_oxide_charging),
        "omissions": (
            "spontaneous neutral-Cl etching",
            "poly-Si inelastic scattering after a nonreactive collision",
            "multiple-bounce SiO2 scattering",
            "electrons inside the shadowed local notch",
        ),
    }
    return HwangGiapisLocalNotchResult2D(
        poly, inventory, oxide_count, depth, float(np.max(depth)), removed,
        threshold_removed_cells, detached_cells,
        threshold_removed_reactive, detached_reactive,
        direct_reactive, scattered_reactive, entries.launched_count,
        landed_poly, landed_oxide, escaped, max(int(front_events), 1),
        provenance, landed_photoresist_count=landed_photoresist)
