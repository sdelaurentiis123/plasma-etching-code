#!/usr/bin/env python3
"""petch trench ARDE on the SAME geometry as vps_trench_arde.py, to compare vs the ViennaPS curve
(ViennaPS: norm rate 1.0/0.89/0.76 at AR 1.3/3.8/6.2). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE = 0.03, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9, 2.6]
GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + 5.0 + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=5.0 + 0.3, hole=False)


def depth(dur, periodic, nt="radiosity"):
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = periodic
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport=nt)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 25)), par=p, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


def arde(nt, periodic):
    depth(0.4, periodic, nt)  # warm
    pts = [(dr, depth(dr, periodic, nt)) for dr in DURS]
    dr = np.array([x[0] for x in pts]); dep = np.array([x[1] for x in pts])
    armid = 0.5 * (dep[1:] + dep[:-1]) / W; rate = np.diff(dep) / np.diff(dr)
    return armid, rate / rate[0], dep


print(f"device={t3.DEVICE}  petch trench ARDE (same geom as ViennaPS)\n", flush=True)
print(f"  ViennaPS reference: norm rate 1.00/0.89/0.76 at AR 1.3/3.8/6.2\n", flush=True)
for nt in ["radiosity", "mc"]:
    ar, nr, dep = arde(nt, 1)
    print(f"  petch {nt} +periodic: depths {np.round(dep,2)}", flush=True)
    print(f"    AR  {np.round(ar,1)}", flush=True)
    print(f"    nr  {np.round(nr,2)}", flush=True)
