#!/usr/bin/env python3
"""Does deterministic RADIOSITY (no deep-floor under-sampling) match the de Boer/Blauw real-wafer ARDE
where MC over-starved? MC gave ~0.13 at AR10 vs the wafer's 0.43. Radiosity solves the multi-bounce
conductance exactly -> should give the gentle experimental ARDE. Compare MC vs radiosity on the same
de Boer trench, at betaE near the wafer's S_F~0.47. Run on a box (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json, time
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0])
EXP_R = np.array([1.0, 0.43, 0.29])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 140


def curve(nt, beta):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta
    fl = petch.Flags(coverage_sticking=True, neutral_transport=nt, sampling="sobol")
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=30000, n_neu=30000, reinit_method="skfmm", verbose=False,
                       record_depth_every=7, **GEO)
    wall = time.time() - t0
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    t = st / NS * GEO['t_end']; r = np.gradient(dd, t); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    nr = np.clip(r / max(r0, 1e-9), 0, 2)
    ours = np.interp(EXP_AR, ar, nr)
    rmse = float(np.sqrt(np.mean((ours - EXP_R) ** 2)))
    return float(ar.max()), ours.tolist(), rmse, wall


print(f"device={t3.DEVICE}  de Boer/Blauw wafer: normARDE {EXP_R} at AR {EXP_AR}\n")
res = {}
for nt, beta in [("mc", 0.47), ("radiosity", 0.47), ("radiosity", 0.30)]:
    armax, ours, rmse, wall = curve(nt, beta)
    key = f"{nt}_b{beta}"
    res[key] = dict(ar_max=armax, ours=ours, rmse=rmse, wall=wall)
    print(f"  {nt:9s} betaE={beta}: AR_max={armax:5.1f}  ours@AR{list(EXP_AR.astype(int))}={np.round(ours,3)}  "
          f"RMSE-vs-wafer={rmse:.3f}  ({wall:.0f}s)", flush=True)
print(f"\n  MC over-starves (~0.13 at AR10); radiosity SHOULD approach the wafer's 0.43 if under-sampling was the cause.")
json.dump(dict(exp_ar=EXP_AR.tolist(), exp_r=EXP_R.tolist(), result=res),
          open("validate_deBoer_radiosity_result.json", "w"), indent=2)
print("wrote validate_deBoer_radiosity_result.json")
