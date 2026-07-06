"""GPU (Warp/CUDA) particle tracer for feature charging -- the kinetic core.

Mirrors the numba `_trace_general` push+hit loop as a Warp kernel so millions of ions/electrons
can be pushed through the self-consistent field ON THE GPU per iteration (like the 3D `_trace3d`
kernel that hit 61x). Same leapfrog with adaptive dt (<=0.45 cell/step), periodic x, z<0.5 = escaped
to plasma. Runs on CPU (parity) or CUDA (production) via the `device` arg. This is step 1 of the
GPU-native kinetic charging engine: field solve + charge deposit stay in numpy/Warp around it.
"""
from __future__ import annotations

import numpy as np

try:
    import warp as wp
    _WARP = True
except Exception:  # pragma: no cover
    wp = None
    _WARP = False


if _WARP:
    @wp.kernel
    def _trace_gpu_kernel(Ex: wp.array2d(dtype=wp.float32), Ez: wp.array2d(dtype=wp.float32),
                          solid: wp.array2d(dtype=wp.int8),
                          x0: wp.array(dtype=wp.float32), z0: wp.array(dtype=wp.float32),
                          vx0: wp.array(dtype=wp.float32), vz0: wp.array(dtype=wp.float32),
                          q: wp.float32, nx: wp.int32, nz: wp.int32, max_steps: wp.int32,
                          dt_cap: wp.float32, dt_field: wp.float32,
                          hit_ix: wp.array(dtype=wp.int32), hit_iz: wp.array(dtype=wp.int32),
                          impact_E: wp.array(dtype=wp.float32)):
        p = wp.tid()
        x = x0[p]; z = z0[p]; vx = vx0[p]; vz = vz0[p]
        hix = wp.int32(-1); hiz = wp.int32(-1); E = wp.float32(0.0)
        xmax = wp.float32(nx)
        for _ in range(max_steps):
            ix = wp.int32(x)
            if ix < 0:
                ix = wp.int32(0)
            elif ix > nx - 2:
                ix = nx - 2
            iz = wp.int32(z)
            if iz < 0:
                iz = wp.int32(0)
            elif iz > nz - 2:
                iz = nz - 2
            fx = Ex[ix, iz]; fz = Ez[ix, iz]
            ax = q * fx * 0.5; az = q * fz * 0.5
            avx = wp.abs(vx); avz = wp.abs(vz)
            vmax = wp.max(avx, avz)
            if vmax < 0.8:
                vmax = wp.float32(0.8)
            dt_v = dt_cap / vmax
            field = wp.sqrt(fx * fx + fz * fz)
            if field < 1.0e-9:
                field = wp.float32(1.0e-9)
            dt_e = dt_field / wp.sqrt(field)
            dt = wp.min(dt_v, dt_e)
            vx_half = vx + 0.5 * ax * dt
            vz_half = vz + 0.5 * az * dt
            xa = x + vx_half * dt
            za = z + vz_half * dt
            ix2 = wp.int32(xa)
            if ix2 < 0:
                ix2 = wp.int32(0)
            elif ix2 > nx - 2:
                ix2 = nx - 2
            iz2 = wp.int32(za)
            if iz2 < 0:
                iz2 = wp.int32(0)
            elif iz2 > nz - 2:
                iz2 = nz - 2
            vx = vx_half + 0.25 * q * Ex[ix2, iz2] * dt
            vz = vz_half + 0.25 * q * Ez[ix2, iz2] * dt
            x = wp.mod(xa, xmax)
            if x < 0.0:
                x = x + xmax
            z = za
            if z < 0.5:
                break
            ixh = wp.int32(x)
            if ixh < 0:
                ixh = wp.int32(0)
            elif ixh > nx - 1:
                ixh = nx - 1
            izh = wp.int32(z)
            if izh < 0:
                izh = wp.int32(0)
            elif izh > nz - 1:
                izh = nz - 1
            if solid[ixh, izh] != wp.int8(0):
                hix = ixh; hiz = izh
                E = vx * vx + vz * vz
                break
        hit_ix[p] = hix
        hit_iz[p] = hiz
        impact_E[p] = E


def trace_gpu(Ex, Ez, solid, x0, z0, vx0, vz0, q, max_steps, device="cpu", dt_cap=0.45, dt_field=0.3):
    """Push particles on `device` ("cpu" or "cuda"); returns (hit_ix, hit_iz, impact_E) numpy arrays.
    hit_ix<0 = escaped (survivor). Same semantics as numba `_trace_general`."""
    if not _WARP:
        raise RuntimeError("warp unavailable")
    nx, nz = solid.shape
    n = x0.shape[0]
    with wp.ScopedDevice(device):
        Ex_d = wp.array(np.ascontiguousarray(Ex, dtype=np.float32), dtype=wp.float32)
        Ez_d = wp.array(np.ascontiguousarray(Ez, dtype=np.float32), dtype=wp.float32)
        solid_d = wp.array(np.ascontiguousarray(solid.astype(np.int8)), dtype=wp.int8)
        x_d = wp.array(x0.astype(np.float32), dtype=wp.float32)
        z_d = wp.array(z0.astype(np.float32), dtype=wp.float32)
        vx_d = wp.array(vx0.astype(np.float32), dtype=wp.float32)
        vz_d = wp.array(vz0.astype(np.float32), dtype=wp.float32)
        hix = wp.empty(n, dtype=wp.int32)
        hiz = wp.empty(n, dtype=wp.int32)
        E = wp.empty(n, dtype=wp.float32)
        wp.launch(_trace_gpu_kernel, dim=n,
                  inputs=[Ex_d, Ez_d, solid_d, x_d, z_d, vx_d, vz_d,
                          float(q), int(nx), int(nz), int(max_steps),
                          float(dt_cap), float(dt_field), hix, hiz, E])
        wp.synchronize()
        return hix.numpy(), hiz.numpy(), E.numpy()
