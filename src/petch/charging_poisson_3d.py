"""Compatible 3-D electrostatics for feature-scale dielectric charging.

The unknown is nodal voltage and the source is physical free charge in coulombs.  Trilinear Q1
elements assemble ``integral epsilon_r grad(Na).grad(Nb) dV`` on a uniform rectilinear grid, so
material-interface displacement continuity follows from the weak form rather than a fitted field law.
The resulting voltage array is directly consumable by :mod:`petch.boundary_transport_3d`.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import splu

from .charging_poisson import EPS0


_LOCAL_OFFSETS = np.asarray(tuple(product((0, 1), repeat=3)), dtype=int)


def _spacing_vector(spacing_m):
    spacing = np.asarray(spacing_m, dtype=float)
    if spacing.ndim == 0:
        spacing = np.full(3, float(spacing))
    if spacing.shape != (3,) or np.any(~np.isfinite(spacing)) or np.any(spacing <= 0.0):
        raise ValueError("spacing_m must be one positive finite value or a length-three vector")
    return spacing


def _local_q1_stiffness_3d(spacing_m):
    """Return the physical-cell trilinear Laplacian matrix for unit relative permittivity."""
    spacing = _spacing_vector(spacing_m)
    gauss = 0.5 + np.asarray((-1.0, 1.0)) / (2.0 * np.sqrt(3.0))
    local = np.zeros((8, 8))
    for coordinate in product(gauss, repeat=3):
        coordinate = np.asarray(coordinate)
        gradient = np.empty((8, 3))
        for node_index, offset in enumerate(_LOCAL_OFFSETS):
            for axis in range(3):
                derivative = 1.0 if offset[axis] else -1.0
                for other_axis in range(3):
                    if other_axis == axis:
                        continue
                    derivative *= (coordinate[other_axis] if offset[other_axis]
                                   else 1.0 - coordinate[other_axis])
                gradient[node_index, axis] = derivative / spacing[axis]
        # Two-point Gauss weights on [0,1] are 1/2 in each direction.
        local += gradient @ gradient.T * np.prod(spacing) / 8.0
    return 0.5 * (local + local.T)


def assemble_q1_stiffness_3d(epsilon_r, spacing_m):
    """Assemble ``integral epsilon_r grad(Na).grad(Nb) dV`` on a 3-D cell grid."""
    epsilon_r = np.asarray(epsilon_r, dtype=float)
    if (epsilon_r.ndim != 3 or epsilon_r.size == 0 or not np.all(np.isfinite(epsilon_r))
            or np.any(epsilon_r <= 0.0)):
        raise ValueError("epsilon_r must be a nonempty finite positive 3-D cell grid")
    local_unit = _local_q1_stiffness_3d(spacing_m)
    node_shape = tuple(np.asarray(epsilon_r.shape) + 1)
    rows = []; columns = []; values = []
    for cell in np.ndindex(epsilon_r.shape):
        node_coordinates = _LOCAL_OFFSETS + np.asarray(cell)
        nodes = np.ravel_multi_index(node_coordinates.T, node_shape)
        local = epsilon_r[cell] * local_unit
        rows.extend(np.repeat(nodes, 8))
        columns.extend(np.tile(nodes, 8))
        values.extend(local.ravel())
    node_count = int(np.prod(node_shape))
    return coo_matrix(
        (values, (rows, columns)), shape=(node_count, node_count)).tocsc()


@dataclass(frozen=True)
class PoissonDiagnostics3D:
    max_abs_free_charge_residual_c: float
    rms_free_charge_residual_c: float
    free_nodes: int
    electrostatic_energy_j: float
    specified_charge_c: float
    dirichlet_reaction_charge_c: float
    charge_balance_c: float


class NodalPoissonSystem3D:
    """Reusable SI-unit Poisson factorization for one fixed 3-D material grid.

    Cell-centered ``epsilon_r`` may vary arbitrarily. ``charge_node_c`` is already projected onto the
    Q1 nodal basis. Boundaries not present in ``dirichlet_mask`` carry the natural zero-normal-
    displacement condition. At least one Dirichlet node is required to fix the voltage gauge.
    """

    def __init__(self, epsilon_r, spacing_m, dirichlet_mask, dirichlet_voltage=None):
        self.epsilon_r = np.asarray(epsilon_r, dtype=float).copy()
        if self.epsilon_r.ndim != 3:
            raise ValueError("epsilon_r must be a 3-D cell grid")
        self.spacing_m = _spacing_vector(spacing_m).copy()
        self.shape = tuple(np.asarray(self.epsilon_r.shape) + 1)
        mask = np.asarray(dirichlet_mask, dtype=bool)
        if mask.shape != self.shape or not np.any(mask):
            raise ValueError("dirichlet_mask must match the nodal grid and fix at least one node")
        if dirichlet_voltage is None:
            voltage = np.zeros(self.shape)
        else:
            voltage = np.asarray(dirichlet_voltage, dtype=float)
            if voltage.shape != self.shape or not np.all(np.isfinite(voltage)):
                raise ValueError("dirichlet_voltage must be a finite nodal grid")
        self.dirichlet_mask = mask.copy()
        self.dirichlet_voltage = voltage.copy()
        self.stiffness = assemble_q1_stiffness_3d(self.epsilon_r, self.spacing_m)
        flat_mask = self.dirichlet_mask.ravel()
        self.fixed = np.flatnonzero(flat_mask)
        self.free = np.flatnonzero(~flat_mask)
        if self.free.size == 0:
            raise ValueError("dirichlet_mask must leave at least one free node")
        self._free_lookup = np.full(flat_mask.size, -1, dtype=int)
        self._free_lookup[self.free] = np.arange(self.free.size)
        self._kff = csc_matrix(self.stiffness[self.free][:, self.free])
        self._kfc = csc_matrix(self.stiffness[self.free][:, self.fixed])
        self._factor = splu(self._kff)

    def solve(self, charge_node_c=None):
        if charge_node_c is None:
            charge = np.zeros(self.shape)
        else:
            charge = np.asarray(charge_node_c, dtype=float)
            if charge.shape != self.shape or not np.all(np.isfinite(charge)):
                raise ValueError("charge_node_c must be a finite nodal grid")
        fixed_voltage = self.dirichlet_voltage.ravel()[self.fixed]
        load = charge.ravel() / EPS0
        right_hand_side = load[self.free] - self._kfc @ fixed_voltage
        flat_voltage = self.dirichlet_voltage.ravel().copy()
        flat_voltage[self.free] = self._factor.solve(right_hand_side)
        voltage = flat_voltage.reshape(self.shape)
        residual_charge = EPS0 * (self.stiffness @ flat_voltage) - charge.ravel()
        free_residual = residual_charge[self.free]
        specified_charge = float(np.sum(charge))
        dirichlet_reaction = float(np.sum(residual_charge[self.fixed]))
        diagnostics = PoissonDiagnostics3D(
            max_abs_free_charge_residual_c=(float(np.max(np.abs(free_residual)))
                                            if free_residual.size else 0.0),
            rms_free_charge_residual_c=(float(np.sqrt(np.mean(free_residual ** 2)))
                                        if free_residual.size else 0.0),
            free_nodes=int(self.free.size),
            electrostatic_energy_j=(
                0.5 * EPS0 * float(flat_voltage @ (self.stiffness @ flat_voltage))),
            specified_charge_c=specified_charge,
            dirichlet_reaction_charge_c=dirichlet_reaction,
            charge_balance_c=specified_charge + dirichlet_reaction)
        return voltage, diagnostics

    def diagonal_capacitance(self, nodes):
        """Return exact diagonal nodal charge-to-voltage capacitances in farads."""
        nodes = np.asarray(nodes, dtype=int)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or nodes.size == 0:
            raise ValueError("nodes must have shape (n,3)")
        if np.any(nodes < 0) or np.any(nodes >= np.asarray(self.shape)):
            raise ValueError("response nodes lie outside the nodal grid")
        flat = np.ravel_multi_index(nodes.T, self.shape)
        free_rows = self._free_lookup[flat]
        if np.any(free_rows < 0):
            raise ValueError("response is defined only for free nodes")
        right_hand_side = np.zeros((self.free.size, len(flat)))
        right_hand_side[free_rows, np.arange(len(flat))] = 1.0 / EPS0
        response = self._factor.solve(right_hand_side)
        diagonal = response[free_rows, np.arange(len(flat))]
        if np.any(diagonal <= 0.0) or not np.all(np.isfinite(diagonal)):
            raise RuntimeError("Poisson response diagonal must be finite and positive")
        return 1.0 / diagonal


def lump_triangle_sheet_charge_3d(
        shape, vertices, faces, sigma_c_per_m2, *, grid_origin=(0.0, 0.0, 0.0),
        grid_spacing=1.0, coordinate_length_unit_m=1.0):
    """Project piecewise-constant triangle sheet charge onto a matching Q1 volume grid.

    Each triangle must lie inside one grid cell, as marching-cubes triangles do. A degree-three exact
    triangle rule integrates the trilinear basis; total charge and its first spatial moments are thus
    conserved without nearest-node assignment. Coordinates may use any declared length unit.
    """
    shape = tuple(int(value) for value in shape)
    vertices = np.asarray(vertices, dtype=float); faces = np.asarray(faces, dtype=int)
    sigma = np.asarray(sigma_c_per_m2, dtype=float)
    origin = np.asarray(grid_origin, dtype=float)
    spacing = np.asarray(grid_spacing, dtype=float)
    if spacing.ndim == 0:
        spacing = np.full(3, float(spacing))
    if (len(shape) != 3 or min(shape) < 2 or vertices.ndim != 2 or vertices.shape[1] != 3
            or faces.ndim != 2 or faces.shape[1] != 3 or sigma.shape != (len(faces),)
            or origin.shape != (3,) or spacing.shape != (3,)
            or np.any(~np.isfinite(vertices)) or np.any(~np.isfinite(sigma))
            or np.any(~np.isfinite(origin)) or np.any(~np.isfinite(spacing))
            or np.any(spacing <= 0.0) or np.any(faces < 0) or np.any(faces >= len(vertices))
            or not np.isfinite(coordinate_length_unit_m) or coordinate_length_unit_m <= 0.0):
        raise ValueError("invalid triangle sheet-charge projection inputs")
    cell_shape = np.asarray(shape) - 1
    normalized_vertices = (vertices - origin) / spacing
    tolerance = 1e-10
    if (np.any(normalized_vertices < -tolerance)
            or np.any(normalized_vertices > cell_shape + tolerance)):
        raise ValueError("triangle vertices lie outside the nodal grid")

    barycentric = np.asarray([
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.6, 0.2, 0.2], [0.2, 0.6, 0.2], [0.2, 0.2, 0.6],
    ])
    quadrature_weight = np.asarray([-27.0 / 48.0, 25.0 / 48.0,
                                    25.0 / 48.0, 25.0 / 48.0])
    charge = np.zeros(shape)
    for face_index, face in enumerate(faces):
        triangle = vertices[face]
        normalized_triangle = normalized_vertices[face]
        lower = np.floor(np.min(normalized_triangle, axis=0) + tolerance).astype(int)
        lower = np.minimum(lower, cell_shape - 1)
        if np.any(np.max(normalized_triangle, axis=0) > lower + 1.0 + tolerance):
            raise ValueError("each triangle must lie within one potential-grid cell")
        physical_edge_a = (triangle[1] - triangle[0]) * float(coordinate_length_unit_m)
        physical_edge_b = (triangle[2] - triangle[0]) * float(coordinate_length_unit_m)
        area_m2 = 0.5 * float(np.linalg.norm(np.cross(physical_edge_a, physical_edge_b)))
        if area_m2 <= 0.0:
            raise ValueError("sheet-charge triangles must have positive area")
        quadrature_point = barycentric @ normalized_triangle
        for point, weight in zip(quadrature_point, quadrature_weight):
            fraction = np.clip(point - lower, 0.0, 1.0)
            nodal_weight = np.prod(np.where(
                _LOCAL_OFFSETS == 1, fraction[None, :], 1.0 - fraction[None, :]), axis=1)
            contribution = float(sigma[face_index]) * area_m2 * float(weight) * nodal_weight
            nodes = _LOCAL_OFFSETS + lower
            np.add.at(charge, tuple(nodes.T), contribution)
    return charge
