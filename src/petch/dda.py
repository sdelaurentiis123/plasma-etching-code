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

import os
import numpy as np
import warp as wp

DEVICE = os.environ.get("PETCH_DEVICE", "cpu")


@wp.kernel
def _dda_gather_kernel(phi: wp.array3d(dtype=wp.float32), origins: wp.array(dtype=wp.vec3),
                       normals: wp.array(dtype=wp.vec3), dirs: wp.array(dtype=wp.vec3),
                       weights: wp.array(dtype=wp.float32), remit: wp.array(dtype=wp.float32),
                       dx: float, z_top: float, max_steps: int, nx: int, ny: int, nz: int,
                       out: wp.array(dtype=wp.float32)):
    """One thread per surface point: march each direction-ray through the occupancy grid,
    escaping to source (1.0) or hitting a wall (its remit). Cosine-weighted -> arriving flux."""
    i = wp.tid()
    o = origins[i]
    n = normals[i]
    num = float(0.0)
    den = float(0.0)
    D = dirs.shape[0]
    for k in range(D):
        dk = dirs[k]
        acc = wp.max(wp.dot(n, dk), 0.0)
        den += weights[k] * acc
        if acc > 1.0e-6:
            px = o[0] + dk[0] * dx * 1.01
            py = o[1] + dk[1] * dx * 1.01
            pz = o[2] + dk[2] * dx * 1.01
            contrib = float(0.0)
            done = int(0)
            for s in range(max_steps):
                if done == 0:
                    if pz >= z_top:
                        contrib = 1.0
                        done = 1
                    else:
                        ix = wp.clamp(int(px / dx + 0.5), 0, nx - 1)
                        iy = wp.clamp(int(py / dx + 0.5), 0, ny - 1)
                        iz = wp.clamp(int(pz / dx + 0.5), 0, nz - 1)
                        if phi[ix, iy, iz] > 0.0:
                            contrib = remit[(ix * ny + iy) * nz + iz]
                            done = 1
                    px += dk[0] * dx
                    py += dk[1] * dx
                    pz += dk[2] * dx
            num += weights[k] * acc * contrib
    out[i] = num / wp.max(den, 1.0e-12)


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
    F = len(centroids)
    if z_top is None:
        z_top = zs[-1] - 0.5 * dx
    dirs, weights = fib_hemisphere(n_dir)
    # solid-side cell index for each face (where its re-emitted flux lives)
    spos = centroids - 0.7 * dx * into_gas_normals
    inv = 1.0 / dx
    six = np.clip((spos[:, 0] * inv + 0.5).astype(np.int64), 0, nx - 1)
    siy = np.clip((spos[:, 1] * inv + 0.5).astype(np.int64), 0, ny - 1)
    siz = np.clip((spos[:, 2] * inv + 0.5).astype(np.int64), 0, nz - 1)
    wall_cell = (six * ny + siy) * nz + siz
    max_steps = int((z_top - zs[0]) / dx) + 4

    # Warp gather (CUDA or CPU): upload static arrays once, re-launch per re-emission iteration.
    phi_wp = wp.array(np.ascontiguousarray(phi, dtype=np.float32), dtype=wp.float32, device=DEVICE)
    orig_wp = wp.array(centroids.astype(np.float32), dtype=wp.vec3, device=DEVICE)
    nrm_wp = wp.array(into_gas_normals.astype(np.float32), dtype=wp.vec3, device=DEVICE)
    dirs_wp = wp.array(dirs.astype(np.float32), dtype=wp.vec3, device=DEVICE)
    w_wp = wp.array(weights.astype(np.float32), dtype=wp.float32, device=DEVICE)
    out_wp = wp.zeros(F, dtype=wp.float32, device=DEVICE)

    def gather(remit_flat):
        remit_wp = wp.array(remit_flat.astype(np.float32), dtype=wp.float32, device=DEVICE)
        wp.launch(_dda_gather_kernel, dim=F, device=DEVICE,
                  inputs=[phi_wp, orig_wp, nrm_wp, dirs_wp, w_wp, remit_wp,
                          float(dx), float(z_top), int(max_steps), int(nx), int(ny), int(nz), out_wp])
        return out_wp.numpy().astype(np.float64)

    direct = gather(np.zeros(nx * ny * nz, np.float32))
    vf = direct.copy()
    alpha = np.clip(1.0 - np.asarray(s_face, float), 0.0, 1.0)
    for _ in range(max(1, n_reemit)):
        remit = np.zeros(nx * ny * nz)
        np.maximum.at(remit, wall_cell, alpha * np.maximum(vf, 0.0))   # face -> solid wall cell
        vf = direct + gather(remit)
    return np.clip(np.maximum(vf, 0.0), 0.0, 8.0)
