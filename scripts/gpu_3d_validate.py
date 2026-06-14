#!/usr/bin/env python3
"""GPU clean validation of the 3D ARDE gap. Run on a Linux+NVIDIA box (PETCH_DEVICE=cuda).

For coverage-sticking OFF vs ON: rate-match each so the d=6 hole reaches ViennaPS depth (~9.24 um),
then measure the d=3/4/6 ARDE at MATCHED depth and compare to ViennaPS-3D [0.832, 0.916, 1.0].
GPU makes deep holes fast, so we can rate-match + use fine resolution + many rays (clean signal).
"""
import os
os.environ["PETCH_DEVICE"] = "cuda"
import numpy as np
import petch
from petch import threed as t3

VPS = {3.0: 7.684, 4.0: 8.469, 6.0: 9.241}      # ViennaPS-3D hole depths
DIA = [3.0, 4.0, 6.0]
VPSN = np.array([VPS[d] for d in DIA]); VPSN = VPSN / VPSN[-1]

GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 40, 30000


def depth(dd, cov, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=cov),
                       n_ion=NR, n_neu=NR, reinit_method="gpu", verbose=False, **GEO)
    return t3.center_depth_3d(g)


def rate_for_d6(cov, target=9.241):
    # bracket + a couple bisections on rate so d6 ~ target
    lo, hi = 0.05, 2.5
    dlo, dhi = depth(6.0, cov, lo), depth(6.0, cov, hi)
    for _ in range(5):
        mid = 0.5 * (lo + hi)
        dm = depth(6.0, cov, mid)
        if dm < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print(f"Warp device check: {t3.DEVICE}")
print(f"ViennaPS-3D hole ARDE target: {np.round(VPSN, 3)}\n")
for cov in [False, True]:
    rate = rate_for_d6(cov)
    d = np.array([depth(dd, cov, rate) for dd in DIA])
    na = d / d[-1]
    rmse = np.sqrt(np.mean((na - VPSN) ** 2))
    tag = "coverage-sticking ON " if cov else "constant sticking OFF"
    print(f"  {tag} (rate={rate:.3f}): depths {np.round(d,2)}  normARDE {np.round(na,3)}  rmse={rmse:.4f}")
