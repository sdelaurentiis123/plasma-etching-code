"""Boundary-fitted nodal electrostatics for feature charging (experimental).

The legacy charging field stores voltage at material cell centres while particles are absorbed at cell
faces.  This module instead places voltage on the vertices of the gas-cell complex.  Every axis-aligned
gas/solid interface is therefore an exact mesh edge and its Dirichlet value is imposed at that edge's
nodes.  Covered solid cells never participate in particle interpolation.

This is deliberately standalone until its analytic and reciprocity gates close.
"""
from __future__ import annotations

import numpy as np


def nodal_domain(solid, surface_voltage):
    """Construct active and Dirichlet node data from a cell-centred material grid.

    Returns ``active, dirichlet, boundary_voltage`` on shape ``(nx+1,nz+1)``. The top plasma boundary
    is grounded. Material-boundary node values are the arithmetic trace of incident face values; at a
    dielectric corner this is the continuous piecewise-linear nodal representation of adjacent face data.
    """
    solid = np.asarray(solid, dtype=bool)
    surface_voltage = np.asarray(surface_voltage, dtype=float)
    if solid.shape != surface_voltage.shape:
        raise ValueError("solid and surface_voltage must have identical shape")
    nx, nz = solid.shape
    gas = ~solid
    active = np.zeros((nx + 1, nz + 1), dtype=bool)
    for di in (0, 1):
        for dj in (0, 1):
            active[di:di + nx, dj:dj + nz] |= gas

    value_sum = np.zeros_like(active, dtype=float)
    value_count = np.zeros_like(active, dtype=np.int32)

    # A solid cell contributes its voltage to both vertices of each face whose neighbour is gas.
    left = np.zeros_like(solid); left[1:] = solid[1:] & gas[:-1]
    right = np.zeros_like(solid); right[:-1] = solid[:-1] & gas[1:]
    top = np.zeros_like(solid); top[:, 1:] = solid[:, 1:] & gas[:, :-1]
    bottom = np.zeros_like(solid); bottom[:, :-1] = solid[:, :-1] & gas[:, 1:]
    ii, jj = np.where(left)
    for i, j in zip(ii, jj):
        value_sum[i, j] += surface_voltage[i, j]; value_count[i, j] += 1
        value_sum[i, j + 1] += surface_voltage[i, j]; value_count[i, j + 1] += 1
    ii, jj = np.where(right)
    for i, j in zip(ii, jj):
        value_sum[i + 1, j] += surface_voltage[i, j]; value_count[i + 1, j] += 1
        value_sum[i + 1, j + 1] += surface_voltage[i, j]; value_count[i + 1, j + 1] += 1
    ii, jj = np.where(top)
    for i, j in zip(ii, jj):
        value_sum[i, j] += surface_voltage[i, j]; value_count[i, j] += 1
        value_sum[i + 1, j] += surface_voltage[i, j]; value_count[i + 1, j] += 1
    ii, jj = np.where(bottom)
    for i, j in zip(ii, jj):
        value_sum[i, j + 1] += surface_voltage[i, j]; value_count[i, j + 1] += 1
        value_sum[i + 1, j + 1] += surface_voltage[i, j]; value_count[i + 1, j + 1] += 1

    dirichlet = (value_count > 0) & active
    boundary_voltage = np.zeros_like(value_sum)
    boundary_voltage[dirichlet] = value_sum[dirichlet] / value_count[dirichlet]
    # Plasma plane has priority at any degenerate material contact with z=0.
    dirichlet[:, 0] |= active[:, 0]
    boundary_voltage[:, 0] = 0.0
    return active, dirichlet, boundary_voltage


def solve_nodal_laplace(solid, surface_voltage, sweeps=1000, omega=1.7, tolerance=None, initial=None):
    """Solve Laplace's equation on the gas-node complex with exact material-edge Dirichlet data.

    Missing neighbours at lateral/bottom exterior boundaries implement zero normal flux. Returns
    ``(V, diagnostics)``. The residual is the graph finite-volume balance on free active nodes.
    """
    active, fixed, fixed_value = nodal_domain(solid, surface_voltage)
    gas = ~np.asarray(solid, dtype=bool)
    # A nodal edge participates iff it bounds at least one gas cell. Endpoint activity alone is
    # insufficient near covered corners and can create a graph connection through solid material.
    edge_x = np.zeros((gas.shape[0], gas.shape[1] + 1), dtype=bool)
    edge_x[:, :-1] |= gas; edge_x[:, 1:] |= gas
    edge_z = np.zeros((gas.shape[0] + 1, gas.shape[1]), dtype=bool)
    edge_z[:-1, :] |= gas; edge_z[1:, :] |= gas
    if initial is None:
        V = np.zeros(active.shape, dtype=float)
    else:
        V = np.asarray(initial, dtype=float).copy()
        if V.shape != active.shape:
            raise ValueError("initial nodal voltage has wrong shape")
    V[fixed] = fixed_value[fixed]
    ni, nj = np.indices(V.shape)
    red = (ni + nj) % 2 == 0
    free = active & ~fixed
    performed = 0
    max_residual = 0.0
    rms_residual = 0.0
    for sweep in range(int(sweeps)):
        for color in (red, ~red):
            num = np.zeros_like(V); degree = np.zeros_like(V)
            mask = edge_x
            num[1:] += np.where(mask, V[:-1], 0.0); degree[1:] += mask
            num[:-1] += np.where(mask, V[1:], 0.0); degree[:-1] += mask
            mask = edge_z
            num[:, 1:] += np.where(mask, V[:, :-1], 0.0); degree[:, 1:] += mask
            num[:, :-1] += np.where(mask, V[:, 1:], 0.0); degree[:, :-1] += mask
            update = num / np.maximum(degree, 1.0)
            selected = free & color
            V[selected] = (1.0 - omega) * V[selected] + omega * update[selected]
            V[fixed] = fixed_value[fixed]
        performed = sweep + 1
        if tolerance is not None or sweep == int(sweeps) - 1:
            residual = nodal_laplace_residual(V, active, fixed, edge_x=edge_x, edge_z=edge_z)
            max_residual = residual["max_abs"]
            rms_residual = residual["rms"]
            if tolerance is not None and max_residual <= tolerance:
                break
    return V, dict(sweeps=performed, max_abs=max_residual, rms=rms_residual,
                   active_nodes=int(active.sum()), free_nodes=int(free.sum()))


def nodal_laplace_residual(V, active, fixed, edge_x=None, edge_z=None):
    """Graph-Laplacian residual on free active nodes, normalized by local degree."""
    V = np.asarray(V, dtype=float); active = np.asarray(active, dtype=bool); fixed = np.asarray(fixed, dtype=bool)
    num = np.zeros_like(V); degree = np.zeros_like(V)
    if edge_x is None:
        edge_x = active[1:] & active[:-1]
    if edge_z is None:
        edge_z = active[:, 1:] & active[:, :-1]
    mask = edge_x
    num[1:] += np.where(mask, V[:-1], 0.0); degree[1:] += mask
    num[:-1] += np.where(mask, V[1:], 0.0); degree[:-1] += mask
    mask = edge_z
    num[:, 1:] += np.where(mask, V[:, :-1], 0.0); degree[:, 1:] += mask
    num[:, :-1] += np.where(mask, V[:, 1:], 0.0); degree[:, :-1] += mask
    free = active & ~fixed
    values = (V - num / np.maximum(degree, 1.0))[free]
    return dict(max_abs=float(np.max(np.abs(values))) if values.size else 0.0,
                rms=float(np.sqrt(np.mean(values * values))) if values.size else 0.0)


def nodal_field_at(V, x, z):
    """Q1 element-local electric field at particle coordinates in a gas cell."""
    nx, nz = V.shape[0] - 1, V.shape[1] - 1
    i = min(max(int(np.floor(x)), 0), nx - 1)
    j = min(max(int(np.floor(z)), 0), nz - 1)
    fx = min(max(x - i, 0.0), 1.0); fz = min(max(z - j, 0.0), 1.0)
    ex = -((1.0 - fz) * (V[i + 1, j] - V[i, j])
           + fz * (V[i + 1, j + 1] - V[i, j + 1]))
    ez = -((1.0 - fx) * (V[i, j + 1] - V[i, j])
           + fx * (V[i + 1, j + 1] - V[i + 1, j]))
    return float(ex), float(ez)
