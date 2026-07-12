"""Warp implementation of the compatible nodal charged-particle tracer."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile

import numpy as np

try:
    import warp as wp
except Exception:  # pragma: no cover
    wp = None


if wp is not None:
    @wp.func
    def _field(V: wp.array2d(dtype=wp.float64), x: wp.float64, z: wp.float64,
               nx: wp.int32, nz: wp.int32):
        one = wp.float64(1.0)
        i = wp.int32(wp.floor(x)); j = wp.int32(wp.floor(z))
        if i < 0:
            i = 0
        elif i > nx - 1:
            i = nx - 1
        if j < 0:
            j = 0
        elif j > nz - 1:
            j = nz - 1
        fx = x - wp.float64(i); fz = z - wp.float64(j)
        fx = wp.clamp(fx, wp.float64(0.0), wp.float64(1.0))
        fz = wp.clamp(fz, wp.float64(0.0), wp.float64(1.0))
        ex = -((one - fz) * (V[i + 1, j] - V[i, j])
               + fz * (V[i + 1, j + 1] - V[i, j + 1]))
        ez = -((one - fx) * (V[i, j + 1] - V[i, j])
               + fx * (V[i + 1, j + 1] - V[i + 1, j]))
        return wp.vec2d(ex, ez)


    @wp.kernel
    def _trace_nodal_kernel(
            V: wp.array2d(dtype=wp.float64), solid: wp.array2d(dtype=wp.int8),
            x0: wp.array(dtype=wp.float64), z0: wp.array(dtype=wp.float64),
            vx0: wp.array(dtype=wp.float64), vz0: wp.array(dtype=wp.float64),
            q: wp.float64, nx: wp.int32, nz: wp.int32, max_steps: wp.int32,
            dt_cap: wp.float64, dt_field: wp.float64, fixed_dt: wp.float64,
            hit_ix: wp.array(dtype=wp.int32), hit_iz: wp.array(dtype=wp.int32),
            impact_E: wp.array(dtype=wp.float64), hit_vx: wp.array(dtype=wp.float64),
            survivor: wp.array(dtype=wp.uint8), exit_vx: wp.array(dtype=wp.float64),
            exit_vz: wp.array(dtype=wp.float64), hit_nx: wp.array(dtype=wp.int8),
            hit_nz: wp.array(dtype=wp.int8), hit_x_position: wp.array(dtype=wp.float64),
            hit_z_position: wp.array(dtype=wp.float64)):
        p = wp.tid()
        zero = wp.float64(0.0); half = wp.float64(0.5)
        one = wp.float64(1.0); two = wp.float64(2.0)
        x = x0[p]; z = z0[p]; vx = vx0[p]; vz = vz0[p]
        alive = wp.int32(1)
        hix = wp.int32(-1); hiz = wp.int32(-1)
        energy = wp.float64(0.0); hvx = wp.float64(0.0)
        sv = wp.uint8(0); evx = wp.float64(0.0); evz = wp.float64(0.0)
        hnx = wp.int8(0); hnz = wp.int8(0)
        hxp = wp.float64(0.0); hzp = wp.float64(0.0)
        for _step in range(max_steps):
            if alive == 1:
                electric = _field(V, x, z, nx, nz)
                vmax = wp.max(wp.max(wp.abs(vx), wp.abs(vz)), wp.float64(0.8))
                field_norm = wp.max(wp.sqrt(electric[0] * electric[0]
                                            + electric[1] * electric[1]), wp.float64(1.0e-9))
                dt = fixed_dt
                if fixed_dt <= 0.0:
                    dt = wp.min(dt_cap / vmax, dt_field / wp.sqrt(field_norm))
                vxn = vx + half * q * electric[0] * dt
                vzn = vz + half * q * electric[1] * dt
                xa = x + half * (vx + vxn) * dt
                za = z + half * (vz + vzn) * dt
                for _mid in range(4):
                    midpoint_field = _field(V, half * (x + xa), half * (z + za), nx, nz)
                    vxn = vx + half * q * midpoint_field[0] * dt
                    vzn = vz + half * q * midpoint_field[1] * dt
                    xa = x + half * (vx + vxn) * dt
                    za = z + half * (vz + vzn) * dt

                dx = xa - x; dz = za - z
                ci = wp.int32(wp.floor(x)); cj = wp.int32(wp.floor(z))
                tx = two; tz = two
                ni = ci; nj = cj
                if dx > 0.0:
                    tx = (wp.float64(ci) + one - x) / dx; ni = ci + 1
                elif dx < 0.0:
                    tx = (wp.float64(ci) - x) / dx; ni = ci - 1
                if dz > 0.0:
                    tz = (wp.float64(cj) + one - z) / dz; nj = cj + 1
                elif dz < 0.0:
                    tz = (wp.float64(cj) - z) / dz; nj = cj - 1

                crossed = wp.int32(0); hi = wp.int32(-1); hj = wp.int32(-1)
                step_nx = wp.int8(0); step_nz = wp.int8(0); hit_fraction = two
                if tx >= 0.0 and tx <= 1.0:
                    test_i = ni; test_j = wp.int32(wp.floor(z + dz * tx))
                    if test_i >= 0 and test_i < nx and test_j >= 0 and test_j < nz:
                        if solid[test_i, test_j] != wp.int8(0):
                            crossed = 1; hi = test_i; hj = test_j; hit_fraction = tx
                            step_nx = wp.int8(-1) if dx > 0.0 else wp.int8(1)
                if tz >= 0.0 and tz <= 1.0:
                    test_i = wp.int32(wp.floor(x + dx * tz)); test_j = nj
                    if test_j < 0 and tz < hit_fraction:
                        evx = vx + tz * (vxn - vx); evz = vz + tz * (vzn - vz)
                        alive = 0
                    elif (test_i >= 0 and test_i < nx and test_j >= 0 and test_j < nz
                          and solid[test_i, test_j] != wp.int8(0) and tz < hit_fraction):
                        crossed = 1; hi = test_i; hj = test_j; hit_fraction = tz
                        step_nx = wp.int8(0)
                        step_nz = wp.int8(-1) if dz > 0.0 else wp.int8(1)
                if alive == 1:
                    if crossed == 0 and tx >= 0.0 and tx <= 1.0 and (ni < 0 or ni >= nx):
                        xa = -xa if ni < 0 else two * wp.float64(nx) - xa
                        vxn = -vxn
                    if crossed == 0 and tz >= 0.0 and tz <= 1.0 and nj >= nz:
                        za = two * wp.float64(nz) - za; vzn = -vzn
                    if crossed == 1:
                        xa = x + dx * hit_fraction; za = z + dz * hit_fraction
                        hit_vxn = vx + hit_fraction * (vxn - vx)
                        hit_vzn = vz + hit_fraction * (vzn - vz)
                        hix = hi; hiz = hj; energy = hit_vxn * hit_vxn + hit_vzn * hit_vzn
                        hvx = hit_vxn; hnx = step_nx; hnz = step_nz; hxp = xa; hzp = za
                        alive = 0
                    else:
                        x = xa; z = za; vx = vxn; vz = vzn
                        if z <= 0.0:
                            evx = vx; evz = vz; alive = 0
        if alive == 1:
            sv = wp.uint8(1)
        hit_ix[p] = hix; hit_iz[p] = hiz; impact_E[p] = energy; hit_vx[p] = hvx
        survivor[p] = sv; exit_vx[p] = evx; exit_vz[p] = evz
        hit_nx[p] = hnx; hit_nz[p] = hnz
        hit_x_position[p] = hxp; hit_z_position[p] = hzp


def trace_nodal_warp(
        V, solid, x0, z0, vx0, vz0, q, nx, nz, max_steps, dt_cap, dt_field,
        fixed_dt=0.0, *, device="cpu"):
    """Run the nodal tracer on a Warp CPU/CUDA device and return NumPy arrays."""
    if wp is None:
        raise RuntimeError("warp unavailable")
    # Sandboxed services and read-only home mounts can make Warp's default user-cache path unusable.
    # Select a process-local temporary cache only after proving the configured directory unwritable.
    # This changes compilation storage, never device choice or numerical behavior.
    cache = Path(wp.config.kernel_cache_dir)
    try:
        probe = Path(tempfile.mkdtemp(prefix="petch-probe-", dir=cache))
        probe.rmdir()
    except OSError:
        fallback = Path(tempfile.gettempdir()) / "petch-warp-cache" / wp.config.version
        fallback.mkdir(parents=True, exist_ok=True)
        wp.config.kernel_cache_dir = os.fspath(fallback)
    inputs = [np.asarray(value) for value in (x0, z0, vx0, vz0)]
    if len({value.shape for value in inputs}) != 1:
        raise ValueError("particle input arrays must have identical shapes")
    count = inputs[0].size
    with wp.ScopedDevice(device):
        V_d = wp.array(np.ascontiguousarray(V, dtype=np.float64), dtype=wp.float64)
        solid_d = wp.array(np.ascontiguousarray(solid, dtype=np.int8), dtype=wp.int8)
        particle = [wp.array(np.ascontiguousarray(value, dtype=np.float64), dtype=wp.float64)
                    for value in inputs]
        hit_ix = wp.empty(count, dtype=wp.int32); hit_iz = wp.empty(count, dtype=wp.int32)
        impact = wp.empty(count, dtype=wp.float64); hit_vx = wp.empty(count, dtype=wp.float64)
        survivor = wp.empty(count, dtype=wp.uint8)
        exit_vx = wp.empty(count, dtype=wp.float64); exit_vz = wp.empty(count, dtype=wp.float64)
        hit_nx = wp.empty(count, dtype=wp.int8); hit_nz = wp.empty(count, dtype=wp.int8)
        hit_x = wp.empty(count, dtype=wp.float64); hit_z = wp.empty(count, dtype=wp.float64)
        wp.launch(
            _trace_nodal_kernel, dim=count, device=device,
            inputs=[V_d, solid_d, *particle, float(q), int(nx), int(nz), int(max_steps),
                    float(dt_cap), float(dt_field), float(fixed_dt), hit_ix, hit_iz, impact,
                    hit_vx, survivor, exit_vx, exit_vz, hit_nx, hit_nz, hit_x, hit_z])
        wp.synchronize_device(device)
        return tuple(array.numpy() for array in (
            hit_ix, hit_iz, impact, hit_vx, survivor, exit_vx, exit_vz,
            hit_nx, hit_nz, hit_x, hit_z))
