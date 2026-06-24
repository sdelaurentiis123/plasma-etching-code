#!/usr/bin/env python3
"""Clean isolation of the sampling question: on a STATIC clean trench (fixed AR), measure the normalized
deep-floor neutral flux at dx=0.25/0.15/0.10 with (a) FIXED ray budget and (b) rays scaled with #facets
(ViennaPS style). No etch, no depth metric. A STABLE fix => floor flux ~constant across dx; the fixed
budget should get noisier/biased as dx refines (fewer rays per facet on the deep floor)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np, warp as wp
from petch import threed as t3
DEV = t3.DEVICE
W, DEPTH, XE, YE = 0.5, 3.0, 1.5, 1.0       # AR=6 trench
SUBTOP = DEPTH + 0.3


def floor_flux(dx, n_rays):
    LZ = SUBTOP + 0.4
    geo = t3.make_trench_3d(XE, YE, LZ, dx, W, 2*dx, SUBTOP, hole=False)
    verts, faces, cen, areas = t3.extract_mesh_3d(geo['phi'], dx)
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
    A_src = XE*YE; z_src = LZ - dx
    o, d = t3.gen_source_gpu('neutral', n_rays, XE, YE, z_src, 0.0, 1)
    fl = wp.zeros(len(faces), dtype=float, device=DEV)
    wp.launch(t3._trace3d_cov_rr, dim=n_rays, device=DEV,
              inputs=[mesh.id, o, d, wp.array(np.ones(len(faces), np.float32), dtype=float, device=DEV),
                      0.1, 7, fl, float(XE), float(YE), float(LZ), 1])
    mF = (fl.numpy()/np.maximum(areas, 1e-9)) / (n_rays/A_src)
    cx = cen[:, 0]-XE/2; cz = cen[:, 2]
    floor = (np.abs(cx) < 0.12) & (cz < cz.min()+3*dx)
    field = cz > SUBTOP                                  # open field above substrate (normalization ref)
    fmF = mF[floor].mean() if floor.sum() else np.nan
    return fmF, int(floor.sum()), len(faces)


print(f"device={DEV}\n  trench W={W} AR={DEPTH/W:.0f}\n", flush=True)
print("dx   | FIXED 40k: floorFlux (nfloor/nfaces)  | rays_per_point=80: floorFlux  (rays)", flush=True)
for dx in [0.25, 0.15, 0.10]:
    a, naf, nf = floor_flux(dx, 40000)
    rpp = 80*nf
    b, nbf, _ = floor_flux(dx, rpp)
    print(f"{dx:.2f} | {a:.3f}  ({naf}/{nf})            | {b:.3f}  ({rpp} rays)", flush=True)
print("\n  STABLE => rays_per_point floorFlux ~constant across dx (well-sampled deep floor at every dx).", flush=True)
print("  fixed budget should drift/noisy as dx shrinks (deep-floor facets starved of rays).", flush=True)
