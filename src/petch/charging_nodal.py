"""Boundary-fitted nodal electrostatics for feature charging (experimental).

The legacy charging field stores voltage at material cell centres while particles are absorbed at cell
faces.  This module instead places voltage on the vertices of the gas-cell complex.  Every axis-aligned
gas/solid interface is therefore an exact mesh edge and its Dirichlet value is imposed at that edge's
nodes.  Covered solid cells never participate in particle interpolation.

This is deliberately standalone until its analytic and reciprocity gates close.
"""
from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange
except Exception:  # pragma: no cover
    njit = None
    prange = range


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


def material_face_nodes(cell, normal):
    """Return the two grid vertices of one axis-aligned gas-facing material face."""
    i, j = map(int, cell); nx, nz = map(int, normal)
    if (nx, nz) == (-1, 0):
        return ((i, j), (i, j + 1))
    if (nx, nz) == (1, 0):
        return ((i + 1, j), (i + 1, j + 1))
    if (nx, nz) == (0, -1):
        return ((i, j), (i + 1, j))
    if (nx, nz) == (0, 1):
        return ((i, j + 1), (i + 1, j + 1))
    raise ValueError("material face normal must be axis aligned")


def solve_nodal_laplace(
        solid, surface_voltage=None, sweeps=1000, omega=1.7, tolerance=None, initial=None,
        boundary_nodal_voltage=None):
    """Solve Laplace's equation on the gas-node complex with exact material-edge Dirichlet data.

    Missing neighbours at lateral/bottom exterior boundaries implement zero normal flux. Returns
    ``(V, diagnostics)``. The residual is the graph finite-volume balance on free active nodes.
    """
    if boundary_nodal_voltage is None:
        if surface_voltage is None:
            raise ValueError("surface_voltage or boundary_nodal_voltage is required")
        active, fixed, fixed_value = nodal_domain(solid, surface_voltage)
    else:
        boundary_nodal_voltage = np.asarray(boundary_nodal_voltage, dtype=float)
        expected = (np.asarray(solid).shape[0] + 1, np.asarray(solid).shape[1] + 1)
        if (boundary_nodal_voltage.shape != expected
                or not np.all(np.isfinite(boundary_nodal_voltage))):
            raise ValueError("boundary_nodal_voltage must be a finite nodal grid")
        active, fixed, _ = nodal_domain(
            solid, np.zeros(np.asarray(solid).shape, dtype=float))
        fixed_value = boundary_nodal_voltage.copy()
        fixed_value[:, 0] = 0.0
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


def _nodal_field_scalar(V, x, z, nx, nz):
    i = int(np.floor(x)); j = int(np.floor(z))
    if i < 0: i = 0
    elif i > nx - 1: i = nx - 1
    if j < 0: j = 0
    elif j > nz - 1: j = nz - 1
    fx = x - i; fz = z - j
    if fx < 0.0: fx = 0.0
    elif fx > 1.0: fx = 1.0
    if fz < 0.0: fz = 0.0
    elif fz > 1.0: fz = 1.0
    ex = -((1.0 - fz) * (V[i + 1, j] - V[i, j])
           + fz * (V[i + 1, j + 1] - V[i, j + 1]))
    ez = -((1.0 - fx) * (V[i, j + 1] - V[i, j])
           + fx * (V[i + 1, j + 1] - V[i + 1, j]))
    return ex, ez


if njit is not None:
    _nodal_field_scalar = njit(cache=True, fastmath=True)(_nodal_field_scalar)


def _trace_nodal_py(V, solid, x0, z0, vx0, vz0, q, nx, nz, max_steps, dt_cap, dt_field,
                    fixed_dt=0.0):
    """Trace particles in the boundary-fitted nodal field with exact first-face absorption.

    Midpoint field iteration is second order in smooth elements. Adaptive substeps move less than
    ``dt_cap`` cells per coordinate. Before accepting a cell change, a DDA face test identifies the
    first crossed face and its adjacent cell, preventing particles from tunnelling through thin solids.
    """
    n = x0.shape[0]
    hit_ix = np.full(n, -1, np.int64); hit_iz = np.full(n, -1, np.int64)
    impact_E = np.zeros(n); hit_vx = np.zeros(n); survivor = np.zeros(n, np.uint8)
    exit_vx = np.zeros(n); exit_vz = np.zeros(n)
    hit_nx = np.zeros(n, np.int8); hit_nz = np.zeros(n, np.int8)
    hit_x_position = np.zeros(n); hit_z_position = np.zeros(n)
    for p in prange(n):
        x = x0[p]; z = z0[p]; vx = vx0[p]; vz = vz0[p]
        alive = True
        for _ in range(max_steps):
            ex, ez = _nodal_field_scalar(V, x, z, nx, nz)
            vmax = max(abs(vx), abs(vz), 0.8)
            field = max((ex * ex + ez * ez) ** 0.5, 1.0e-9)
            dt = fixed_dt if fixed_dt > 0.0 else min(dt_cap / vmax, dt_field / field ** 0.5)
            vxn = vx + 0.5 * q * ex * dt
            vzn = vz + 0.5 * q * ez * dt
            xa = x + 0.5 * (vx + vxn) * dt
            za = z + 0.5 * (vz + vzn) * dt
            for _mid in range(4):
                xm = 0.5 * (x + xa); zm = 0.5 * (z + za)
                emx, emz = _nodal_field_scalar(V, xm, zm, nx, nz)
                vxn = vx + 0.5 * q * emx * dt
                vzn = vz + 0.5 * q * emz * dt
                xa = x + 0.5 * (vx + vxn) * dt
                za = z + 0.5 * (vz + vzn) * dt

            dx = xa - x; dz = za - z
            ci = int(np.floor(x)); cj = int(np.floor(z))
            tx = 2.0; tz = 2.0; ni = ci; nj = cj
            if dx > 0.0:
                tx = (ci + 1.0 - x) / dx; ni = ci + 1
            elif dx < 0.0:
                tx = (ci - x) / dx; ni = ci - 1
            if dz > 0.0:
                tz = (cj + 1.0 - z) / dz; nj = cj + 1
            elif dz < 0.0:
                tz = (cj - z) / dz; nj = cj - 1

            crossed = False; hi = -1; hj = -1; hnx = 0; hnz = 0; hit_fraction = 2.0
            # Test every grid face crossed by this short segment. The orthogonal cell index is
            # evaluated AT the crossing, so a gas-gas crossing followed by a solid crossing cannot
            # tunnel through a corner.
            if tx >= 0.0 and tx <= 1.0:
                test_i = ni; test_j = int(np.floor(z + dz * tx))
                if test_i >= 0 and test_i < nx and test_j >= 0 and test_j < nz:
                    if solid[test_i, test_j]:
                        crossed = True; hi = test_i; hj = test_j
                        hnx = -1 if dx > 0.0 else 1; hnz = 0; hit_fraction = tx
            if tz >= 0.0 and tz <= 1.0:
                test_i = int(np.floor(x + dx * tz)); test_j = nj
                if test_j < 0 and tz < hit_fraction:
                    exit_vx[p] = vx + tz * (vxn - vx)
                    exit_vz[p] = vz + tz * (vzn - vz)
                    alive = False; break
                if (test_i >= 0 and test_i < nx and test_j >= 0 and test_j < nz
                        and solid[test_i, test_j] and tz < hit_fraction):
                    crossed = True; hi = test_i; hj = test_j; hit_fraction = tz
                    hnx = 0; hnz = -1 if dz > 0.0 else 1
            if not crossed and tx >= 0.0 and tx <= 1.0 and (ni < 0 or ni >= nx):
                # Specular symmetry reflection keeps the REMAINDER of the proposed step. Clipping
                # at the wall shortens the orbit and is not time reversible.
                xa = -xa if ni < 0 else 2.0 * nx - xa
                vxn = -vxn
            if not crossed and tz >= 0.0 and tz <= 1.0 and nj >= nz:
                za = 2.0 * nz - za; vzn = -vzn
            if crossed:
                xa = x + dx * hit_fraction; za = z + dz * hit_fraction
                hit_vxn = vx + hit_fraction * (vxn - vx)
                hit_vzn = vz + hit_fraction * (vzn - vz)
                hit_ix[p] = hi; hit_iz[p] = hj
                impact_E[p] = hit_vxn * hit_vxn + hit_vzn * hit_vzn; hit_vx[p] = hit_vxn
                hit_nx[p] = hnx; hit_nz[p] = hnz
                hit_x_position[p] = xa; hit_z_position[p] = za
                alive = False; break
            x = xa; z = za; vx = vxn; vz = vzn
            if z <= 0.0:
                exit_vx[p] = vx; exit_vz[p] = vz; alive = False; break
        if alive:
            survivor[p] = 1
    return (hit_ix, hit_iz, impact_E, hit_vx, survivor, exit_vx, exit_vz,
            hit_nx, hit_nz, hit_x_position, hit_z_position)


trace_nodal = (njit(cache=True, parallel=True, fastmath=True)(_trace_nodal_py)
               if njit is not None else _trace_nodal_py)
