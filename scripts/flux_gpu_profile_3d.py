#!/usr/bin/env python3
"""GPU-side breakdown of ONE coupled-flux call: where do the ~90ms/step actually go after FSM+warm?
Times each sub-op with wp.synchronize_device barriers (so GPU work is attributed correctly, not hidden
behind async launches). Tells us the real next lever: ray launch (-> cuBQL), .numpy() syncs / array
churn (-> device-resident loop), smooth (-> within-step cache), or source-gen. Run: PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DX, DIAM, N = 0.25, 6.0, 30000
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=0.7)
DEV = t3.DEVICE


def sync():
    wp.synchronize_device(DEV)


def tmed(fn, reps=40):
    fn(); sync()
    t0 = time.time()
    for _ in range(reps):
        fn()
    sync()
    return 1000.0 * (time.time() - t0) / reps


# deep mesh
g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=12, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=10000, n_neu=10000, reinit_method="fsm", verbose=False, **GEO)
verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
F = len(faces)
mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV),
               indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
vv = verts[faces]; fn = np.cross(vv[:, 1]-vv[:, 0], vv[:, 2]-vv[:, 0])
fn = fn/(np.linalg.norm(fn, axis=1, keepdims=True)+1e-12)
pairs = t3._edge_adjacency(faces)
rng = np.random.default_rng(0)
o, d = t3._source3d('neutral', N, 14, 14, GEO['Lz']-DX, 0.04, 'sobol', rng, 7)
bare = np.ones(F, np.float32)
flux = np.abs(rng.normal(1.0, 0.3, F))

print(f"GPU flux breakdown  F={F} faces  pairs={len(pairs)}  rays={N}\n", flush=True)

t_src = tmed(lambda: t3._source3d('neutral', N, 14, 14, GEO['Lz']-DX, 0.04, 'sobol', rng, 7))


def upload():
    return (wp.array(o, dtype=wp.vec3, device=DEV), wp.array(d, dtype=wp.vec3, device=DEV),
            wp.array(bare, dtype=float, device=DEV), wp.zeros(F, dtype=float, device=DEV))
t_upload = tmed(upload)

ow, dw, bw, flw = upload()


def launch():
    flw.zero_()
    wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, ow, dw, bw, 0.7, 7, flw])
t_launch = tmed(launch)


def launch_sync():
    flw.zero_()
    wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, ow, dw, bw, 0.7, 7, flw])
    _ = flw.numpy()
t_launch_sync = tmed(launch_sync)

t_smooth_gpu = tmed(lambda: t3.smooth_flux_gpu(flux, fn, pairs, 1, 1.0))
t_smooth_cpu = tmed(lambda: t3.smooth_flux(flux, fn, pairs, 1, 1.0))

print(f"  _source3d (sobol gen, host) : {t_src:6.2f} ms")
print(f"  upload o/d/bare/zeros (H2D)  : {t_upload:6.2f} ms")
print(f"  ray launch ALONE (async)     : {t_launch:6.2f} ms")
print(f"  ray launch + .numpy() (D2H)  : {t_launch_sync:6.2f} ms   (sync cost = {t_launch_sync-t_launch:.2f} ms)")
print(f"  smooth_flux_gpu              : {t_smooth_gpu:6.2f} ms")
print(f"  smooth_flux (numpy)          : {t_smooth_cpu:6.2f} ms")
per_neu = t_src + t_upload + t_launch_sync + t_smooth_gpu
print(f"\n  per neutral() ~ {per_neu:.1f} ms; warm n_fp=1 = 1 ion + 2 neutral ~ {3*per_neu:.0f} ms/step flux")
print("  -> dominant term names the next lever.")
