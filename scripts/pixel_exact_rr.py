#!/usr/bin/env python3
"""Pixel-exact test with the EXACT ViennaPS neutral transport (weighted ray + Russian roulette, no
bounce cap; now the default in mc_flux_3d_coupled). HYPOTHESIS: with the unbiased estimator, ViennaPS's
OWN sticking betaE=0.7 (not the 0.08 hack that compensated for our old truncation bias) gives a
pixel-exact ARDE match. Sweeps betaE around 0.7; QMC + many rays to cut noise; robust center_depth
metric (~ViennaPS global-min for a clean etch). Run on the box (PETCH_DEVICE=cuda)."""
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
DX, NS, NR = 0.25, 40, 30000
BETAS = [0.5, 0.7, 0.9]
RATES = [0.04, 0.07, 0.10, 0.13, 0.17, 0.22, 0.28, 0.36, 0.46]   # finer; CFL cap now 160 (no blowup)


def run(dd, beta, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = beta   # RR transport is default
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


print(f"device={t3.DEVICE}  EXACT RR transport; betaE sweep around ViennaPS's 0.7")
print("ViennaPS trajectory:", [(round(v['d6'],2), [round(x,3) for x in v['norm']]) for v in vps], "\n")
result = {}
for be in BETAS:
    print(f"=== betaE={be} ===")
    rows = []
    for r in RATES:
        d = np.array([run(dd, be, r) for dd in DIAMS])
        if np.any(d <= 0):
            print(f"  rate={r:.3f}  {np.round(d,2)} [invalid]", flush=True); continue
        rows.append(dict(rate=r, d=d.tolist(), norm=(d/d[-1]).tolist(), d6=float(d[-1])))
        print(f"  rate={r:.3f}  d={np.round(d,2)}  norm={np.round(d/d[-1],3)}  d6={d[-1]:.2f}", flush=True)
    per = []
    for v in vps:
        cand = [x for x in rows if x['d6'] > 0]
        if not cand:
            continue
        x = min(cand, key=lambda x: abs(x['d6'] - v['d6']))
        rmse = float(np.sqrt(np.mean((np.array(x['norm']) - np.array(v['norm'])) ** 2)))
        per.append(dict(vps_d6=v['d6'], vps_norm=v['norm'], d6=x['d6'], norm=x['norm'], rmse=rmse))
        print(f"  vps d6={v['d6']:5.2f} {np.round(v['norm'],3)} | ours d6={x['d6']:5.2f} {np.round(x['norm'],3)}  rmse={rmse:.3f}")
    mr = float(np.mean([p['rmse'] for p in per])) if per else None
    result[be] = dict(rows=rows, per=per, mean_rmse=mr)
    print(f"  -> betaE={be} mean ARDE rmse vs ViennaPS = {mr}\n", flush=True)

print("=" * 64)
for be in BETAS:
    print(f"  betaE={be}: mean ARDE rmse = {result[be]['mean_rmse']}")
best = min((b for b in BETAS if result[b]['mean_rmse'] is not None), key=lambda b: result[b]['mean_rmse'], default=None)
if best is not None:
    print(f"\n  BEST betaE={best}  mean rmse={result[best]['mean_rmse']:.3f}")
json.dump(dict(vps=vps, result={str(k): v for k, v in result.items()}, best=best),
          open("pixel_exact_rr_result.json", "w"), indent=2)
print("\nwrote pixel_exact_rr_result.json")
