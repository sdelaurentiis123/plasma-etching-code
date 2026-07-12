"""Compatible variable-permittivity electrostatics for physical surface charging.

The charged-particle mover and this solver use the same nodal potential. Dielectric state is free
surface charge (C/m in the 2-D per-unit-depth model), not an independently prescribed surface voltage.
The Q1 weak form enforces dielectric-interface flux continuity without a fitted charge-to-voltage map.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import splu


EPS0 = 8.8541878128e-12

# Bilinear Q1 Laplacian stiffness on a unit square. The cell-size factors cancel in 2-D per unit
# out-of-plane depth. Node order is lower-left, lower-right, upper-right, upper-left.
_Q1_STIFFNESS = np.array([
    [4.0, -1.0, -2.0, -1.0],
    [-1.0, 4.0, -1.0, -2.0],
    [-2.0, -1.0, 4.0, -1.0],
    [-1.0, -2.0, -1.0, 4.0],
]) / 6.0


def assemble_q1_stiffness(epsilon_r):
    """Assemble ``integral epsilon_r grad(Na).grad(Nb) dA`` on a rectilinear cell grid."""
    epsilon_r = np.asarray(epsilon_r, dtype=float)
    if (epsilon_r.ndim != 2 or epsilon_r.size == 0
            or not np.all(np.isfinite(epsilon_r)) or np.any(epsilon_r <= 0.0)):
        raise ValueError("epsilon_r must be a nonempty finite positive 2-D cell grid")
    nx, nz = epsilon_r.shape
    node_nz = nz + 1
    rows = []; cols = []; values = []
    for i in range(nx):
        for j in range(nz):
            nodes = np.array([
                i * node_nz + j,
                (i + 1) * node_nz + j,
                (i + 1) * node_nz + j + 1,
                i * node_nz + j + 1,
            ], dtype=int)
            local = epsilon_r[i, j] * _Q1_STIFFNESS
            rows.extend(np.repeat(nodes, 4)); cols.extend(np.tile(nodes, 4))
            values.extend(local.ravel())
    node_count = (nx + 1) * (nz + 1)
    return coo_matrix((values, (rows, cols)), shape=(node_count, node_count)).tocsc()


@dataclass
class PoissonDiagnostics:
    max_abs_residual_v: float
    rms_residual_v: float
    free_nodes: int
    electrostatic_energy_j_per_m: float
    specified_charge_c_per_m: float
    dirichlet_reaction_charge_c_per_m: float
    charge_balance_c_per_m: float


class NodalPoissonSystem:
    """Reusable physical-units Poisson factorization on a fixed material geometry.

    ``charge_node_c_per_m`` is free line charge per unit out-of-plane depth, already lumped to Q1
    nodes. In a surface edge of length ``h``, a uniform sheet charge ``sigma`` contributes
    ``sigma*h/2`` to each endpoint. Natural outer boundaries are zero normal displacement unless a
    node is included in ``dirichlet_mask``.
    """

    def __init__(self, epsilon_r, dirichlet_mask, dirichlet_voltage=None):
        self.epsilon_r = np.asarray(epsilon_r, dtype=float).copy()
        self.shape = (self.epsilon_r.shape[0] + 1, self.epsilon_r.shape[1] + 1)
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
        self.stiffness = assemble_q1_stiffness(self.epsilon_r)
        flat_mask = self.dirichlet_mask.ravel()
        self.fixed = np.flatnonzero(flat_mask)
        self.free = np.flatnonzero(~flat_mask)
        self._free_lookup = np.full(flat_mask.size, -1, dtype=int)
        self._free_lookup[self.free] = np.arange(self.free.size)
        self._kff = csc_matrix(self.stiffness[self.free][:, self.free])
        self._kfc = csc_matrix(self.stiffness[self.free][:, self.fixed])
        self._factor = splu(self._kff)

    def solve(self, charge_node_c_per_m=None):
        if charge_node_c_per_m is None:
            charge = np.zeros(self.shape)
        else:
            charge = np.asarray(charge_node_c_per_m, dtype=float)
            if charge.shape != self.shape or not np.all(np.isfinite(charge)):
                raise ValueError("charge_node_c_per_m must be a finite nodal grid")
        fixed_voltage = self.dirichlet_voltage.ravel()[self.fixed]
        load = charge.ravel() / EPS0
        rhs = load[self.free] - self._kfc @ fixed_voltage
        flat_voltage = self.dirichlet_voltage.ravel().copy()
        flat_voltage[self.free] = self._factor.solve(rhs)
        voltage = flat_voltage.reshape(self.shape)
        residual = self.stiffness @ flat_voltage - load
        free_residual = residual[self.free]
        specified_charge = float(np.sum(charge))
        dirichlet_reaction_charge = EPS0 * float(np.sum(residual[self.fixed]))
        energy = 0.5 * EPS0 * float(flat_voltage @ (self.stiffness @ flat_voltage))
        diagnostics = PoissonDiagnostics(
            max_abs_residual_v=(float(np.max(np.abs(free_residual)))
                                if free_residual.size else 0.0),
            rms_residual_v=(float(np.sqrt(np.mean(free_residual ** 2)))
                            if free_residual.size else 0.0),
            free_nodes=int(self.free.size),
            electrostatic_energy_j_per_m=energy,
            specified_charge_c_per_m=specified_charge,
            dirichlet_reaction_charge_c_per_m=dirichlet_reaction_charge,
            charge_balance_c_per_m=specified_charge + dirichlet_reaction_charge)
        return voltage, diagnostics

    def diagonal_surface_capacitance(self, surface_nodes):
        """Return exact diagonal charge-to-voltage capacitance for selected free nodes [F/m]."""
        nodes = np.asarray(surface_nodes, dtype=int)
        if nodes.ndim != 2 or nodes.shape[1] != 2 or nodes.size == 0:
            raise ValueError("surface_nodes must have shape (n,2)")
        if (np.any(nodes < 0) or np.any(nodes[:, 0] >= self.shape[0])
                or np.any(nodes[:, 1] >= self.shape[1])):
            raise ValueError("surface_nodes lie outside the nodal grid")
        flat = np.ravel_multi_index((nodes[:, 0], nodes[:, 1]), self.shape)
        free_rows = self._free_lookup[flat]
        if np.any(free_rows < 0):
            raise ValueError("surface response is defined only for free nodes")
        rhs = np.zeros((self.free.size, len(flat)))
        rhs[free_rows, np.arange(len(flat))] = 1.0 / EPS0
        response = self._factor.solve(rhs)
        diagonal = response[free_rows, np.arange(len(flat))]
        if np.any(diagonal <= 0.0) or not np.all(np.isfinite(diagonal)):
            raise RuntimeError("Poisson response diagonal must be finite and positive")
        return 1.0 / diagonal


def lump_edge_sheet_charge(shape, face_nodes, sigma_c_per_m2, cell_size_m):
    """Mass-lump uniform edge sheet charges to physical nodal line charge [C/m]."""
    charge = np.zeros(tuple(shape), dtype=float)
    sigma = np.asarray(sigma_c_per_m2, dtype=float)
    if sigma.shape != (len(face_nodes),) or not np.all(np.isfinite(sigma)):
        raise ValueError("one finite sheet-charge value is required per face")
    h = float(cell_size_m)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("cell_size_m must be finite and positive")
    for value, endpoints in zip(sigma, face_nodes):
        if len(endpoints) != 2:
            raise ValueError("each face must contain two endpoint nodes")
        for node in endpoints:
            charge[tuple(node)] += 0.5 * value * h
    return charge
