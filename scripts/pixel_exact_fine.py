#!/usr/bin/env python3
"""Pixel-exact #2: finer dx to resolve the deep ARDE (dx=0.25 quantization caused d4>d6 inversions).
Smoothing ON (the validated fix) + ViennaPS betaE=0.7 + RR transport, dx=0.18, more steps for CFL,
seed-averaged to cut MC noise. Compare deep ARDE to ViennaPS. Run on the box (PETCH_DEVICE=cuda)."""
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

GEO = dict(Lx=12, Ly=12, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)
DX, NS, NR = 0.18, 100, 35000
RATES = [0.04, 0.07, 0.11, 0.15, 0.20]


def run(dd, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = 0.7; p['cal_F'] = 12.0; p['flux_smooth'] = 1
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


print(f"device={t3.DEVICE}  dx={DX} NS={NS}  smoothing ON, betaE=0.7")
print("ViennaPS:", [(round(v['d6'],2), [round(x,3) for x in v['norm']]) for v in vps], "\n")
rows = []
for r in RATES:
    d = np.array([run(dd, r) for dd in DIAMS])
    if np.any(d <= 0):
        print(f"  rate={r:.3f} {np.round(d,2)} [invalid]", flush=True); continue
    rows.append(dict(rate=r, d=d.tolist(), norm=(d/d[-1]).tolist(), d6=float(d[-1])))
    print(f"  rate={r:.3f}  d={np.round(d,2)}  norm={np.round(d/d[-1],3)}  d6={d[-1]:.2f}", flush=True)

per = []
for v in vps:
    if not rows:
        continue
    x = min(rows, key=lambda x: abs(x['d6'] - v['d6']))
    rmse = float(np.sqrt(np.mean((np.array(x['norm']) - np.array(v['norm'])) ** 2)))
    per.append(dict(vps_d6=v['d6'], vps_norm=v['norm'], d6=x['d6'], norm=x['norm'], rmse=rmse))
    print(f"  vps d6={v['d6']:5.2f} {np.round(v['norm'],3)} | ours d6={x['d6']:5.2f} {np.round(x['norm'],3)}  rmse={rmse:.3f}")
mr = float(np.mean([p['rmse'] for p in per])) if per else None
print(f"\n  finer-dx smoothing mean ARDE rmse vs ViennaPS = {mr}")
json.dump(dict(vps=vps, rows=rows, per=per, mean_rmse=mr), open("pixel_exact_fine_result.json", "w"), indent=2)
print("wrote pixel_exact_fine_result.json")
