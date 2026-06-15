#!/usr/bin/env python3
"""Per-step SPEED breakdown for the 3D etch loop + the radiosity-vs-8xMC neutral question.

Two things we need to know before chasing "lightning fast":
  1. Where does each step's wall-clock actually go? (timings dict: mesh/flux/extend/advect/reinit)
     -> tells us the NEXT bottleneck to attack on the GPU.
  2. Is deterministic RADIOSITY (1 neutral ray launch + a cheap host sparse solve) actually faster
     than the coverage-coupled MC path (n_fp=4 -> 8 neutral ray launches + 8 host syncs per step),
     AND does it land at the same depth (accuracy-neutral)?  If yes, radiosity is a free speed win.

Matched 3D hole (same geometry as head_to_head_3d.py). Run on a GPU box: PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
NS, NRAY = 30, 30000


def run(label, flags, par):
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=par, flags=flags,
                       n_ion=NRAY, n_neu=NRAY, reinit_method="skfmm", verbose=False, **GEO)
    wall = time.time() - t0
    tm = g['timings']
    depth = t3.max_depth_3d(g)
    other = tm['total'] - (tm['mesh'] + tm['flux'] + tm['extend'] + tm['advect'] + tm['reinit'])
    print(f"\n=== {label} ===", flush=True)
    print(f"  wall {wall:6.2f}s   depth {depth:5.2f}um   nsub_max {tm['nsub_max']}")
    tot = max(tm['total'], 1e-9)
    for k in ['mesh', 'flux', 'extend', 'advect', 'reinit']:
        print(f"    {k:7s} {tm[k]:6.2f}s  {100*tm[k]/tot:4.1f}%")
    print(f"    {'host':7s} {other:6.2f}s  {100*other/tot:4.1f}%   (surface_rate / coverage / numpy)")
    return dict(label=label, wall=wall, depth=depth, timings=tm)


print(f"device={t3.DEVICE}  matched hole d={DIAM}um dx={DX} steps={NS} rays={NRAY}\n", flush=True)
res = []

# 1) coverage-coupled MC neutrals (the accurate default): 8 neutral ray launches/step
p = dict(petch.PAR)
res.append(run("MC neutrals (coverage_sticking, n_fp=4 -> 8 launches/step)",
               petch.Flags(coverage_sticking=True, sampling="sobol"), p))

# 2) deterministic radiosity neutrals: 1 ray launch/step + host sparse solve
res.append(run("RADIOSITY neutrals (1 launch/step + sparse solve)",
               petch.Flags(neutral_transport="radiosity", sampling="sobol"), p))

mc, rad = res[0], res[1]
sx = mc['wall'] / max(rad['wall'], 1e-3)
dd = abs(mc['depth'] - rad['depth'])
print(f"\nradiosity is {sx:.2f}x {'faster' if sx>1 else 'SLOWER'} than 8x-MC neutrals; "
      f"depth delta {dd:.2f}um ({'accuracy-neutral' if dd < 0.5 else 'CHECK accuracy'})", flush=True)
json.dump(dict(geo=GEO, dx=DX, diam=DIAM, ns=NS, nray=NRAY, result=res, radiosity_speedup=sx,
               depth_delta=dd), open("profile_steps_result.json", "w"), indent=2)
print("wrote profile_steps_result.json")
