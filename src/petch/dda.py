"""Deterministic discrete-ordinates (DDA) neutral transport for petch.

Ported from Craig Xu Chen's plasma_sim DDA gather (gpu_flux.py / solver3d._neutral_flux),
adapted to petch's grid (phi > 0 = solid, phi < 0 = gas; grid origin at 0; z = depth axis,
z increases upward). The benchmark (cross_validate_dda) showed this deterministic
discrete-ordinates method gives a clean high-AR ARDE rolloff with NO Monte-Carlo
floor-starvation — unlike petch's MC at low ray counts. This is petch's noise-free deep-AR
neutral option (`neutral_transport="dda"`).

Method: for each surface point, sum the arriving neutral flux over a FIXED set of
upper-hemisphere directions (cosine-weighted). Each direction-ray marches cell-by-cell
through the occupancy grid; a ray that escapes to the open field (top) carries the source
flux (1.0), a ray that hits a wall carries that wall cell's re-emitted flux from the previous
iteration. The re-emission fixed point (Gamma = direct + K(1-S)Gamma) is iterated like
plasma_sim's "source" solver. No mesh BVH, no random sampling -> deterministic, noise-free.
"""
from __future__ import annotations

import numpy as np


def fib_hemisphere(n):
    """n directions on the upper hemisphere (dz>0), ~uniform, with cosine (dz) weights."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - i / n)              # 0..pi/2 -> upper hemisphere
    gold = np.pi * (1.0 + 5.0 ** 0.5)
    theta = gold * i
    dz = np.cos(phi)
    r = np.sin(phi)
    d = np.stack([r * np.cos(theta), r * np.sin(theta), dz], axis=1)
    w = np.clip(dz, 0.0, None)               # cosine emission weight
    return d.astype(np.float64), w.astype(np.float64)


def _gather(phi, dx, origins, normals, dirs, weights, remit_flat, z_top, max_steps):
    """Cosine-normalized arriving flux per origin from a DDA march over the occupancy grid.

    origins (N,3) just-into-gas points; normals (N,3) into-gas surface normals; dirs (D,3),
    weights (D,) the fixed quadrature. remit_flat (nx*ny*nz,) per-cell re-emitted flux (0 for a
    pure direct gather); a ray reaching z >= z_top carries src=1.0. Returns flux (N,)."""
    nx, ny, nz = phi.shape
    N = len(origins)
    D = len(dirs)
    solid = phi > 0.0                         # petch: phi > 0 = solid
    num = np.zeros(N)
    den = np.zeros(N)
    inv = 1.0 / dx
    for k in range(D):
        dk = dirs[k]
        acc = np.maximum(normals @ dk, 0.0)   # cos(theta) of this direction vs each normal
        den += weights[k] * acc
        active = acc > 1.0e-6
        if not active.any():
            continue
        idx = np.flatnonzero(active)
        pos = origins[idx] + dk * dx * 1.01    # step off the surface into gas
        contrib = np.zeros(len(idx))
        done = np.zeros(len(idx), dtype=bool)
        for _ in range(max_steps):
            live = ~done
            if not live.any():
                break
            p = pos[live]
            reached = p[:, 2] >= z_top
            ix = np.clip((p[:, 0] * inv + 0.5).astype(np.int64), 0, nx - 1)
            iy = np.clip((p[:, 1] * inv + 0.5).astype(np.int64), 0, ny - 1)
            iz = np.clip((p[:, 2] * inv + 0.5).astype(np.int64), 0, nz - 1)
            hit = solid[ix, iy, iz] & (~reached)
            live_idx = np.flatnonzero(live)
            r_idx = live_idx[reached]
            h_idx = live_idx[hit]
            contrib[r_idx] = 1.0                                   # escaped to source/open field
            contrib[h_idx] = remit_flat[(ix[hit] * ny + iy[hit]) * nz + iz[hit]]
            done[r_idx] = True
            done[h_idx] = True
            pos[live] = p + dk * dx
        num[idx] += weights[k] * acc[idx] * contrib
    return num / np.maximum(den, 1.0e-12)


def dda_neutral_flux(phi, dx, zs, centroids, into_gas_normals, s_face,
                     n_dir=64, n_reemit=12, z_top=None):
    """Per-face neutral arrival multiplier via deterministic DDA + diffuse re-emission.

    s_face (F,) per-face sticking (= bare*beta). Re-emission deposited on the SOLID wall cell
    adjacent to each face (one step along -into_gas_normal) so a marching ray that stops at the
    wall picks it up. Iterates Gamma = direct + gather((1-s)Gamma). Returns flux (F,)."""
    nx, ny, nz = phi.shape
    if z_top is None:
        z_top = zs[-1] - 0.5 * dx
    dirs, weights = fib_hemisphere(n_dir)
    origins = centroids
    # solid-side cell index for each face (where its re-emitted flux lives)
    spos = centroids - 0.7 * dx * into_gas_normals
    inv = 1.0 / dx
    six = np.clip((spos[:, 0] * inv + 0.5).astype(np.int64), 0, nx - 1)
    siy = np.clip((spos[:, 1] * inv + 0.5).astype(np.int64), 0, ny - 1)
    siz = np.clip((spos[:, 2] * inv + 0.5).astype(np.int64), 0, nz - 1)
    wall_cell = (six * ny + siy) * nz + siz
    max_steps = int((z_top - zs[0]) / dx) + 4

    zero = np.zeros(nx * ny * nz)
    direct = _gather(phi, dx, origins, into_gas_normals, dirs, weights, zero, z_top, max_steps)
    vf = direct.copy()
    alpha = np.clip(1.0 - np.asarray(s_face, float), 0.0, 1.0)
    for _ in range(max(1, n_reemit)):
        remit = np.zeros(nx * ny * nz)
        np.maximum.at(remit, wall_cell, alpha * np.maximum(vf, 0.0))   # face -> solid wall cell
        vf = direct + _gather(phi, dx, origins, into_gas_normals, dirs, weights, remit, z_top, max_steps)
    return np.clip(np.maximum(vf, 0.0), 0.0, 8.0)
