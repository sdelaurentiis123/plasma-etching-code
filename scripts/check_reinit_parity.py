#!/usr/bin/env python3
"""Validate the fixed GPU reinit: (1) depth parity with skfmm, (2) |grad phi| ~ 1 near the front,
(3) no NaN. Uses a deeper cov-ON etch so the front is well-developed (not a 1-cell starved floor)."""
import numpy as np, petch
from petch import threed as t3

GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)


def grad_band_stats(phi, dx, band=4):
    g = np.gradient(phi, dx)
    G = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2)
    m = np.abs(phi) < band * dx
    Gb = G[m]
    return float(np.mean(Gb)), float(np.std(Gb)), float(np.percentile(Gb, 99))


out = {}
for meth in ["skfmm", "gpu"]:
    p = dict(petch.PAR); p['rate_scale'] = 0.15
    g = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=40, par=p,
                       flags=petch.Flags(coverage_sticking=True),
                       n_ion=20000, n_neu=20000, reinit_method=meth, verbose=False, **GEO)
    phi = g['phi']
    d = t3.center_depth_3d(g)
    gm, gs, g99 = grad_band_stats(phi, 0.25)
    out[meth] = (d, gm, gs, g99, bool(np.isnan(phi).any()))
    print(f"reinit={meth:6s}  depth={d:5.2f}  |grad|band mean={gm:.3f} std={gs:.3f} p99={g99:.3f}  nan={out[meth][4]}")

dd = abs(out['gpu'][0] - out['skfmm'][0])
print(f"\ndepth |gpu - skfmm| = {dd:.2f} um ({dd/0.25:.1f} cells)")
print("PARITY OK" if dd <= 0.5 and not out['gpu'][4] and abs(out['gpu'][1]-1.0) < 0.05 else "CHECK")
