#!/usr/bin/env python3
"""High-quality SMOOTH etch movies: petch trench + hole at fine dx and high per-facet ray sampling
(ViennaPS-style rays_per_point) so the walls come out clean, not MC-noisy. Captures frames -> /root/
smooth_<feat>.pkl for rendering. GPU. Run: PETCH_DEVICE=cuda PYTHONPATH=/root/petch/src python smooth_movies.py"""
import os, pickle
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np, petch
from petch import threed as t3

DX, W, SUB = 0.04, 1.0, 9.0      # fine grid, deep substrate
DUR, STEPS = 2.6, 52


def run(hole, rate=0.12):
    Ly = (3.2 if hole else 0.32)
    GEO = dict(Lx=3.2, Ly=Ly, Lz=2*DX+SUB+0.4, dx=DX, trench_width=W, mask_th=2*DX, sub_top=SUB, hole=hole)
    p_ = dict(petch.PAR); p_['rate_scale'] = rate
    if not hole:
        p_['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", ion_reflection=True)
    # rays_per_point=70 -> total rays scale with #facets -> constant sampling per facet -> smooth walls
    g = t3.run_etch_3d(t_end=DUR, n_steps=STEPS, par=p_, flags=fl, n_ion=40000, n_neu=40000,
                       rays_per_point=70, reinit_method="fsm", verbose=True,
                       record_depth_every=1, record_frames=True, **GEO)
    return g


for hole in [False, True]:
    feat = "hole" if hole else "trench"
    print(f"=== {feat} (dx={DX}, rays_per_point=70) ===", flush=True)
    g = run(hole)
    pickle.dump(dict(W=W, DX=DX, sub_top=SUB, Lx=g['Lx'], Lz=g['Lz'], xs=g['xs'], zs=g['zs'], hole=hole,
                     frames=[dict(t=f['t'], depth=f['depth'], phi=f['phi_xz'].astype(np.float32)) for f in g['frames']]),
                open(f"/root/smooth_{feat}.pkl", "wb"))
    print(f"saved /root/smooth_{feat}.pkl  final depth {t3.max_depth_3d(g):.2f}  frames {len(g['frames'])}", flush=True)
print("DONE", flush=True)
