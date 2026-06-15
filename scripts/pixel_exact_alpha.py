#!/usr/bin/env python3
"""Pixel-exact calibration: sweep the flux-smoothing STRENGTH alpha to hit ViennaPS exactly. Without
smoothing our ARDE is too STEEP (rmse 0.155); with full smoothing slightly too FLAT (~0.08) -- so an
intermediate alpha brackets ViennaPS. Find alpha minimizing the depth-resolved ARDE rmse. betaE=0.7
(ViennaPS), RR transport, NS=80 (CFL), center_depth. Run on the box (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

GT = json.load(open("viennaps_3d_depth_resolved.json"))
DIAMS = [3.0, 4.0, 6.0]
vps = []
for dur in GT["durs"]:
    dd = np.array([GT["depth_grid"][str(dur)][str(x)] for x in DIAMS])
    vps.append(dict(d6=float(dd[-1]), norm=(dd / dd[-1]).tolist()))

GEO = dict(Lx=14, Ly=14, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 80, 30000
ALPHAS = [0.0, 0.25, 0.45, 0.65, 0.85, 1.0]
RATES = [0.02, 0.035, 0.055, 0.08, 0.11, 0.15, 0.20]


def run(dd, alpha, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = 0.7; p['cal_F'] = 12.0
    p['flux_smooth'] = 1; p['flux_smooth_alpha'] = alpha
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


print(f"device={t3.DEVICE}  flux-smoothing strength (alpha) calibration vs ViennaPS")
print("ViennaPS:", [(round(v['d6'],2), [round(x,3) for x in v['norm']]) for v in vps], "\n")
result = {}
for al in ALPHAS:
    rows = []
    for r in RATES:
        d = np.array([run(dd, al, r) for dd in DIAMS])
        if np.any(d <= 0):
            continue
        rows.append(dict(rate=r, norm=(d/d[-1]).tolist(), d6=float(d[-1])))
    per = []
    for v in vps:
        if not rows:
            continue
        x = min(rows, key=lambda x: abs(x['d6'] - v['d6']))
        per.append(float(np.sqrt(np.mean((np.array(x['norm']) - np.array(v['norm'])) ** 2))))
    mr = float(np.mean(per)) if per else None
    result[al] = dict(rows=rows, mean_rmse=mr)
    print(f"  alpha={al:.2f}: mean ARDE rmse vs ViennaPS = {mr}", flush=True)

print("=" * 56)
valid = {a: result[a]['mean_rmse'] for a in ALPHAS if result[a]['mean_rmse'] is not None}
best = min(valid, key=valid.get) if valid else None
if best is not None:
    print(f"  BEST alpha={best}  mean ARDE rmse={valid[best]:.4f}")
    b = min(result[best]['rows'], key=lambda x: abs(x['d6'] - 9.24)) if result[best]['rows'] else None
    if b:
        print(f"  (at d6~9.24: ours norm={np.round(b['norm'],3)} vs ViennaPS [0.831,0.917,1.0])")
json.dump(dict(vps=vps, result={str(k): v for k, v in result.items()}, best=best),
          open("pixel_exact_alpha_result.json", "w"), indent=2)
print("\nwrote pixel_exact_alpha_result.json")
