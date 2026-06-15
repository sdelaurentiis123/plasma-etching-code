#!/usr/bin/env python3
"""Pin where the per-step flux time goes in the CURRENT config (KDTree gone). Time, on a DEEP mesh with
sync barriers: GPU warm-start query, ion trace, neutral RR trace (the suspect), device smooth. If the
neutral trace dominates, the flux is RAY-bound (RR bounces in deep holes), not host-bound."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time
import numpy as np
import warp as wp
import petch
from petch import threed as t3

DEV = t3.DEVICE
DX, W, N = 0.25, 6.0, 30000
GEO = dict(Lx=14, Ly=14, Lz=40, mask_th=2, sub_top=34, hole=True, t_end=2.4)


def sync():
    wp.synchronize_device(DEV)


def tmed(fn, reps=30):
    fn(); sync()
    t0 = time.time()
    for _ in range(reps):
        fn()
    sync()
    return 1000.0 * (time.time() - t0) / reps


# DEEP hole
g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=36, par=dict(petch.PAR, n_fp=1, gpu_source=True, gpu_mesh=True, gpu_warmstart=True, flux_smooth_gpu=True),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=20000, n_neu=20000, reinit_method="fsm", verbose=False, **GEO)
verts, faces, centroids, areas = t3.extract_mesh_3d_gpu(g['phi'], DX)
F = len(faces); A = np.maximum(areas, 0.3 * np.median(areas))
mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV), indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
zc = centroids[:, 2]; ar = (GEO['sub_top'] - zc[zc < GEO['sub_top']].min()) / W
print(f"deep mesh F={F} faces, AR~{ar:.1f}\n", flush=True)

o, d = t3.gen_source_gpu('neutral', N, GEO['Lx'], GEO['Ly'], GEO['Lz']-DX, 0.04, 7)
oi, di = t3.gen_source_gpu('ion', N, GEO['Lx'], GEO['Ly'], GEO['Lz']-DX, petch.PAR['ion_ang_sigma'], 1)
bare = wp.array(np.ones(F, np.float32), dtype=float, device=DEV)
fl = wp.zeros(F, dtype=float, device=DEV); fi = wp.zeros(F, dtype=float, device=DEV); ai = wp.zeros(F, dtype=float, device=DEV)
cpts = wp.array(centroids.astype(np.float32), dtype=wp.vec3, device=DEV)
outb = wp.zeros(F, dtype=float, device=DEV)


def ion_trace():
    wp.launch(t3._trace3d, dim=N, device=DEV, inputs=[mesh.id, oi, di, 1.0, 0, 1, 0, 0.34, 0.8, fi, ai])


def neu_trace():
    wp.launch(t3._trace3d_cov_rr, dim=N, device=DEV, inputs=[mesh.id, o, d, bare, 0.7, 2, fl])


def warmstart():
    wp.launch(t3._nearest_face_cov, dim=F, device=DEV, inputs=[mesh.id, cpts, bare, outb])


t_ion = tmed(ion_trace)
t_neu = tmed(neu_trace)
t_ws = tmed(warmstart)
print(f"  ion trace (directional)        : {t_ion:6.2f} ms")
print(f"  neutral RR trace (the suspect) : {t_neu:6.2f} ms")
print(f"  GPU warm-start query           : {t_ws:6.2f} ms")
print(f"\n  per step ~ 1 ion + 2 neutral + 1 warmstart = {t_ion + 2*t_neu + t_ws:.1f} ms")
print("  -> if neutral RR dominates, the flux is RAY-bound -> wavefront / bounce-control, not host.")
