#!/usr/bin/env python3
"""GPU benchmark: the 3D Warp flux kernel on RT cores vs CPU, and a full 3D etch on GPU.

Run on a Linux+NVIDIA box. Measures primary-ray throughput (rays/sec) of the _trace3d kernel
(wp.mesh_query_ray BVH traversal -> RT cores on GPU) on CPU vs CUDA, on a realistic etched-trench
mesh, then times a full 3D etch loop on the GPU.
"""
import time
import numpy as np
import warp as wp
import petch
from petch import threed as t3

wp.init()
HAS_CUDA = wp.get_cuda_device_count() > 0
print(f"Warp {wp.config.version} | CUDA devices: {wp.get_cuda_device_count()} | "
      f"{wp.get_cuda_device_count() and wp.get_cuda_device(0).name or 'cpu-only'}")


def get_mesh():
    """A realistic etched-trench mesh (run a few CPU steps to deepen it)."""
    p = dict(petch.PAR); p['rate_scale'] = 0.3
    geo = t3.run_etch_3d(Lx=10, Ly=4, Lz=14, dx=0.3, trench_width=4, mask_th=2, sub_top=10,
                         t_end=2.0, n_steps=8, par=p, flags=petch.Flags(),
                         n_ion=6000, n_neu=6000, verbose=False)
    return t3.extract_mesh_3d(geo['phi'], 0.3)


def bench(dev, verts, faces, N=300000, n_re=12, s=0.2):
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=dev),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=dev))
    rng = np.random.default_rng(0)
    origin = np.stack([rng.uniform(0, 10, N), rng.uniform(0, 4, N),
                       np.full(N, 13.5)], axis=1).astype(np.float32)
    ct = np.sqrt(rng.uniform(0, 1, N)); st = np.sqrt(1 - ct**2); ph = rng.uniform(0, 2*np.pi, N)
    dirs = np.stack([st*np.cos(ph), st*np.sin(ph), -ct], axis=1).astype(np.float32)
    o = wp.array(origin, dtype=wp.vec3, device=dev)
    d = wp.array(dirs, dtype=wp.vec3, device=dev)
    flux = wp.zeros(len(faces), dtype=float, device=dev)
    ang = wp.zeros(len(faces), dtype=float, device=dev)
    # _trace3d(mesh, origin, dir0, sticking, n_reemit, seed, specular, cos_thr, eta, flux, angacc)
    args = [mesh.id, o, d, float(s), int(n_re), 1, int(0), 0.34, 0.8, flux, ang]
    wp.launch(t3._trace3d, dim=N, device=dev, inputs=args); wp.synchronize()  # warmup/compile
    best = 1e18
    for _ in range(3):
        flux.zero_(); ang.zero_()
        t0 = time.time()
        wp.launch(t3._trace3d, dim=N, device=dev, inputs=args); wp.synchronize()
        best = min(best, time.time() - t0)
    return N / best, best


def main():
    verts, faces, _, _ = get_mesh()
    print(f"\nbench mesh: {len(verts)} verts, {len(faces)} faces; N=300k neutral rays, 12 bounces")
    cpu_rps, cpu_t = bench("cpu", verts, faces)
    print(f"  CPU : {cpu_rps/1e6:6.2f} M primary-rays/s  ({cpu_t*1e3:.1f} ms)")
    if HAS_CUDA:
        gpu_rps, gpu_t = bench("cuda", verts, faces)
        print(f"  GPU : {gpu_rps/1e6:6.2f} M primary-rays/s  ({gpu_t*1e3:.1f} ms)  "
              f"-> {gpu_rps/cpu_rps:.1f}x faster than CPU")

        # full 3D etch on GPU
        import os
        os.environ["PETCH_DEVICE"] = "cuda"
        import importlib; importlib.reload(t3)
        p = dict(petch.PAR); p['rate_scale'] = 0.3
        t0 = time.time()
        geo = t3.run_etch_3d(Lx=10, Ly=4, Lz=14, dx=0.3, trench_width=4, mask_th=2, sub_top=10,
                             t_end=2.0, n_steps=12, par=p, flags=petch.Flags(sampling="sobol"),
                             n_ion=20000, n_neu=20000, verbose=False)
        print(f"\n  full 3D etch on GPU (12 steps, 60k rays/step): {time.time()-t0:.1f}s, "
              f"flux {geo['timings']['flux']:.1f}s, depth {t3.center_depth_3d(geo):.2f} um")
    else:
        print("  (no CUDA device — CPU only)")


if __name__ == "__main__":
    main()
