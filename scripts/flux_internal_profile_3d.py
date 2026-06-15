#!/usr/bin/env python3
"""Where does the per-step flux time go? The flux line is 72% of the loop after FSM+warm-start, but at
30k rays on a 3090 the GPU traversal should be sub-ms -- so the cost is likely HOST post-processing:
Sobol source gen (rebuilt every call!), flux smoothing, coverage algebra, .numpy() syncs.

These host ops are pure numpy -> identical cost on the M1 and the box, so we measure them locally.
Build a realistic DEEP mesh, then time each host sub-op at 30k rays. PETCH_DEVICE=cpu is fine.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cpu")
import time
import numpy as np
import petch
from petch import threed as t3

DX, DIAM, N = 0.25, 6.0, 30000
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=0.6)

# build a deep mesh by etching a few steps
g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=10, par=dict(petch.PAR, n_fp=1),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=8000, n_neu=8000, reinit_method="fsm", verbose=False, **GEO)
verts, faces, centroids, areas = t3.extract_mesh_3d(g['phi'], DX)
F = len(faces)
vv = verts[faces]
fn = np.cross(vv[:, 1] - vv[:, 0], vv[:, 2] - vv[:, 0])
fn = fn / (np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12)
pairs = t3._edge_adjacency(faces)
rng = np.random.default_rng(0)
flux = np.abs(rng.normal(1.0, 0.3, F))
m_i = np.clip(flux, 0, 1.5); cos_i = np.clip(flux, 0, 1)
par = dict(petch.PAR); fl = petch.Flags(coverage_sticking=True)
Lz = GEO['Lz']; z_src = Lz - DX


def tm(fn_call, reps=30):
    fn_call()  # warm
    t0 = time.time()
    for _ in range(reps):
        fn_call()
    return 1000.0 * (time.time() - t0) / reps


print(f"deep mesh: F={F} faces, edge-pairs={len(pairs)}, rays={N}\n", flush=True)
t_src_sobol = tm(lambda: t3._source3d('neutral', N, 14, 14, z_src, 0.04, 'sobol', rng, 7))
t_src_pseudo = tm(lambda: t3._source3d('neutral', N, 14, 14, z_src, 0.04, 'pseudo', rng, 7))
t_smooth = tm(lambda: t3.smooth_flux(flux, fn, pairs, 1, 1.0))
t_belen = tm(lambda: t3._belen_coverages(m_i, flux, flux, cos_i, par, fl))

print(f"  _source3d  sobol  : {t_src_sobol:6.2f} ms   (rebuilds the Sobol engine every call)")
print(f"  _source3d  pseudo : {t_src_pseudo:6.2f} ms")
print(f"  smooth_flux       : {t_smooth:6.2f} ms")
print(f"  _belen_coverages  : {t_belen:6.2f} ms")
# per step (warm n_fp=1): 1 ion + 2 neutral source-gen + 3 smooth + 1 belen
per_step_host = 3 * t_src_sobol + 3 * t_smooth + t_belen
print(f"\n  host per step (3 source-gen + 3 smooth + 1 belen) ~ {per_step_host:.1f} ms", flush=True)
print(f"  measured flux line was ~110 ms/step -> host fraction ~ {100*per_step_host/110:.0f}%")
print(f"  -> if host-dominated, GPU the smoothing/source (or cache Sobol), NOT faster BVH (cuBQL).")
