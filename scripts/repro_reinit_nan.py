#!/usr/bin/env python3
"""Local reproduction of the GPU-reinit NaN blowup (Warp runs on CPU on the M1).

Deep d=6 hole, cov-OFF, rate 0.05 -> on the A4000 this NaN'd at step ~16-20 with reinit='gpu'
and was rock-stable with reinit='skfmm'. Reproduce here, then fix the kernel, then re-verify.
"""
import numpy as np, petch
from petch import threed as t3

GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
for meth in ["gpu", "skfmm"]:
    p = dict(petch.PAR); p['rate_scale'] = 0.05
    g = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=40, par=p,
                       flags=petch.Flags(coverage_sticking=False),
                       n_ion=8000, n_neu=8000, reinit_method=meth, verbose=True, **GEO)
    d = t3.center_depth_3d(g)
    phi = g['phi']
    print(f"reinit={meth}  FINAL depth={d:.2f}  any_nan={np.isnan(phi).any()}  "
          f"phi_range=[{np.nanmin(phi):.1f},{np.nanmax(phi):.1f}]")
    print("----")
