#!/usr/bin/env python3
"""Step 5 spike: prove the Warp bridge works on this M1 (CPU).

Two checks:
  (1) Parity — a ray-vs-segment nearest-hit flux kernel as a `wp.kernel` on device='cpu'
      reproduces a numpy reference exactly. This is the 2D analogue of the flux kernel that
      becomes a DDA/BVH traversal on the GPU in Phase 1.
  (2) Autodiff — `wp.Tape` produces a correct gradient on the CPU backend (the differentiable
      substrate the inverse-design thesis needs).

If both pass, the same kernels run unchanged on NVIDIA RT-core GPUs later.
"""
import numpy as np
import warp as wp

wp.init()
DEVICE = "cpu"


# ---- (1) ray-segment nearest-hit flux kernel (forward parity) ----
@wp.kernel
def trace_direct(seg: wp.array(dtype=wp.vec4),
                 ox: wp.array(dtype=float), oy: wp.array(dtype=float),
                 ux: wp.array(dtype=float), uy: wp.array(dtype=float),
                 flux: wp.array(dtype=float)):
    p = wp.tid()
    px = ox[p]; py = oy[p]; dx = ux[p]; dy = uy[p]
    best_t = float(1.0e18)
    best_s = int(-1)
    for s in range(seg.shape[0]):
        ax = seg[s][0]; ay = seg[s][1]; bx = seg[s][2]; by = seg[s][3]
        ex = bx - ax; ey = by - ay
        den = dx * ey - dy * ex
        if wp.abs(den) > 1.0e-14:
            t = ((ax - px) * ey - (ay - py) * ex) / den
            u = ((ax - px) * dy - (ay - py) * dx) / den
            if t > 1.0e-6 and u >= -1.0e-6 and u <= 1.0 + 1.0e-6:
                if t < best_t:
                    best_t = t
                    best_s = s
    if best_s >= 0:
        wp.atomic_add(flux, best_s, 1.0)


def numpy_reference(seg, ox, oy, ux, uy):
    M = seg.shape[0]
    flux = np.zeros(M)
    for p in range(len(ox)):
        best_t = 1e18; best_s = -1
        px, py, dx, dy = ox[p], oy[p], ux[p], uy[p]
        for s in range(M):
            ax, ay, bx, by = seg[s]
            ex, ey = bx - ax, by - ay
            den = dx * ey - dy * ex
            if abs(den) < 1e-14:
                continue
            t = ((ax - px) * ey - (ay - py) * ex) / den
            u = ((ax - px) * dy - (ay - py) * dx) / den
            if t > 1e-6 and -1e-6 <= u <= 1 + 1e-6 and t < best_t:
                best_t, best_s = t, s
        if best_s >= 0:
            flux[best_s] += 1.0
    return flux


def check_parity():
    rng = np.random.default_rng(0)
    # a few horizontal "floor" segments + random downward rays
    M = 12
    xs = np.linspace(0, 10, M + 1)
    seg = np.stack([xs[:-1], np.zeros(M), xs[1:], np.zeros(M)], axis=1).astype(np.float32)
    N = 4000
    ox = rng.uniform(0, 10, N).astype(np.float32)
    oy = np.full(N, 5.0, np.float32)
    ang = rng.normal(0, 0.3, N).astype(np.float32)
    ux = np.sin(ang).astype(np.float32)
    uy = (-np.cos(ang)).astype(np.float32)

    ref = numpy_reference(seg, ox, oy, ux, uy)

    seg_w = wp.array(seg, dtype=wp.vec4, device=DEVICE)
    flux_w = wp.zeros(M, dtype=float, device=DEVICE)
    wp.launch(trace_direct, dim=N, device=DEVICE,
              inputs=[seg_w,
                      wp.array(ox, dtype=float, device=DEVICE),
                      wp.array(oy, dtype=float, device=DEVICE),
                      wp.array(ux, dtype=float, device=DEVICE),
                      wp.array(uy, dtype=float, device=DEVICE),
                      flux_w])
    got = flux_w.numpy()
    ok = np.array_equal(got, ref)
    print(f"(1) parity: warp kernel vs numpy reference -> {'MATCH' if ok else 'MISMATCH'} "
          f"(total hits warp={int(got.sum())} ref={int(ref.sum())})")
    return ok


# ---- (2) autodiff on CPU ----
@wp.kernel
def sq_loss(x: wp.array(dtype=float), loss: wp.array(dtype=float)):
    i = wp.tid()
    wp.atomic_add(loss, 0, x[i] * x[i])


def check_autodiff():
    x_np = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    x = wp.array(x_np, dtype=float, device=DEVICE, requires_grad=True)
    loss = wp.zeros(1, dtype=float, device=DEVICE, requires_grad=True)
    tape = wp.Tape()
    with tape:
        wp.launch(sq_loss, dim=len(x_np), device=DEVICE, inputs=[x, loss])
    tape.backward(loss=loss)
    grad = x.grad.numpy()
    expected = 2.0 * x_np          # d(sum x^2)/dx = 2x
    ok = np.allclose(grad, expected, atol=1e-5)
    print(f"(2) autodiff: wp.Tape grad={grad.tolist()} expected={expected.tolist()} "
          f"-> {'OK' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    print(f"Warp {wp.config.version} on device='{DEVICE}'")
    p = check_parity()
    a = check_autodiff()
    print("\nWARP BRIDGE:", "READY" if (p and a) else "ISSUES",
          "— same kernels run on NVIDIA GPU in Phase 1." if (p and a) else "")
