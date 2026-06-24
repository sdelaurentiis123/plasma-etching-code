"""Dump petch trench cross-section (half-width vs depth) at matched deep etch -> /root/petch_prof.json.
Compare to ViennaPS (vps_prof.json): is petch's trench more tapered (narrower floor -> starved -> steeper
ARDE)? PETCH_DEVICE=cuda."""
import os, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE, SUB = 0.04, 0.5, 1.5, 0.3, 7.0
TARGET = 4.9   # match ViennaPS deep depth (~AR 10); fudge-free petch is slower so use a longer dur


def run(dur):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
    return t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 22)), par=p, flags=fl, n_ion=60000, n_neu=60000,
                          reinit_method="fsm", verbose=False, **GEO)


# find a duration reaching ~TARGET depth
g = run(3.0); dep = t3.center_depth_3d(g)
print(f"petch dur3.0 center depth {dep:.2f} um (target ~{TARGET})", flush=True)
phi, xs, ys, zs = g['phi'], g['xs'], g['ys'], g['zs']
sub_top = g['sub_top']; cx = XE / 2.0
jmid = len(ys) // 2
prof = {}
for z in np.arange(0.0, dep, 0.2):
    kz = int(np.argmin(np.abs(zs - (sub_top - z))))      # grid index at depth z below substrate top
    gas = phi[:, jmid, kz] < 0                            # gas (etched) at this depth, mid-y
    if gas.any():
        xg = xs[gas]
        inside = xg[np.abs(xg - cx) < W]
        if len(inside):
            prof[round(float(z), 2)] = float(np.abs(inside - cx).max())
print(f"petch trench depth {dep:.2f} um", flush=True)
for z in sorted(prof): print(f"  z={z:.1f}: half-width {prof[z]:.3f}", flush=True)
json.dump({'depth': float(dep), 'halfwidth': prof, 'W': W}, open("/root/petch_prof.json", "w"))
print("wrote /root/petch_prof.json", flush=True)
