#!/usr/bin/env python3
"""Can we MATCH the de Boer/Blauw real-wafer ARDE by using the EXPERIMENT'S OWN F reaction probability?
We over-starved at betaE=0.7 (the ViennaPS value, too sticky). The de Boer/Blauw experiment is
F-transport-limited with S_F ~ 0.47. Sweep betaE around S_F: lower betaE = more wall reflection =
deeper penetration = gentler ARDE. Hypothesis: betaE ~ 0.47 (the wafer's actual F reaction probability)
reproduces the measured normalized-rate-vs-aspect-ratio curve. Run on a box (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0])           # de Boer/Blauw (drop AR40 -- hard to reach cleanly)
EXP_R = np.array([1.0, 0.43, 0.29])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS, NR = 0.25, 140, 30000
BETAS = [0.7, 0.47, 0.30, 0.15]


def arde_curve(beta):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False,
                       record_depth_every=7, **GEO)
    hist = g['depth_history']
    steps = np.array([h[0] for h in hist]); dd = np.array([h[1] for h in hist])
    t = steps / NS * GEO['t_end']
    rate = np.gradient(dd, t); ar = dd / W
    ref = rate[ar < 2.0]; r0 = float(ref.max()) if len(ref) else float(rate.max())
    nrate = np.clip(rate / max(r0, 1e-9), 0, 2)
    ours = np.interp(EXP_AR, ar, nrate)
    rmse = float(np.sqrt(np.mean((ours - EXP_R) ** 2)))
    return float(ar.max()), ours.tolist(), rmse


print(f"device={t3.DEVICE}  de Boer/Blauw experiment: normARDE {EXP_R} at AR {EXP_AR}")
print(f"  (experiment F reaction prob S_F ~ 0.47)\n")
res = {}
for be in BETAS:
    armax, ours, rmse = arde_curve(be)
    res[be] = dict(ar_max=armax, ours=ours, rmse=rmse)
    print(f"  betaE={be:.2f}: AR_max={armax:5.1f}  ours@AR{list(EXP_AR.astype(int))}={np.round(ours,3)}  "
          f"RMSE-vs-wafer={rmse:.3f}", flush=True)
best = min(BETAS, key=lambda b: res[b]['rmse'])
print(f"\n  BEST betaE={best}  RMSE={res[best]['rmse']:.3f}  (experiment S_F~0.47)")
print("  -> if best betaE ~ 0.47 = the wafer's own F reaction probability, our transport is physically right.")
json.dump(dict(exp_ar=EXP_AR.tolist(), exp_r=EXP_R.tolist(),
               result={str(k): v for k, v in res.items()}, best=best),
          open("validate_deBoer_betaE_result.json", "w"), indent=2)
print("\nwrote validate_deBoer_betaE_result.json")
