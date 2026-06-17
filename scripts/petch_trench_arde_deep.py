#!/usr/bin/env python3
"""petch trench ARDE, DEEPER domain so AR can reach ~12 without the domain-floor clamp that
contaminated the AR>9 points in petch_trench_arde.py. Same DX/W/XE/YE as vps_trench_arde.py.
ViennaPS ref points: depth 1.305/2.460/3.705/4.969 at dur 0.4/0.8/1.3/1.9. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE = 0.03, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9, 2.6, 3.4]
SUB = 8.0
GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)


def depth(dur, nt):
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport=nt)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 25)), par=p, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


print(f"device={t3.DEVICE}  petch trench ARDE DEEP (sub={SUB}, same geom as ViennaPS)\n", flush=True)
for nt in ["radiosity", "mc"]:
    depth(0.4, nt)  # warm
    dep = np.array([depth(dr, nt) for dr in DURS])
    dr = np.array(DURS)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dr)
    nr = rate / rate[0]
    print(f"  petch {nt}: depths {np.round(dep,3).tolist()}", flush=True)
    print(f"    AR  {np.round(armid,2).tolist()}", flush=True)
    print(f"    nr  {np.round(nr,3).tolist()}", flush=True)
