#!/usr/bin/env python3
"""Pixel-exact push vs depth-resolved ViennaPS-3D. Tests the two no-new-physics alignments:
(1) adaptive re-emission bounce cap (now default; feeds the narrow HARC floor at low betaE),
(2) max_depth metric matching ViennaPS's -surfaceNode.z.min() (vs our median-center center_depth).

Sweeps betaE x rate with coverage_sticking; for EACH config records d3/d4/d6 under BOTH metrics; for
each ViennaPS depth point finds our nearest-d6 config and reports normalized-ARDE rmse under each metric.
Goal: a single betaE whose ARDE tracks ViennaPS's [d3/d6,d4/d6] across all four depths to rmse -> ~0.
Run on the box (PETCH_DEVICE=cuda)."""
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
    vps.append(dict(dur=dur, d=dd.tolist(), norm=(dd / dd[-1]).tolist(), d6=float(dd[-1])))

GEO = dict(Lx=14, Ly=14, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 40, 30000
BETAS = [0.08, 0.15, 0.25, 0.40]
RATES = [0.05, 0.09, 0.14, 0.20, 0.28, 0.40]


def run(dd, beta, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = beta
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    return t3.center_depth_3d(g), t3.max_depth_3d(g)


def norm_rmse(d, vnorm):
    d = np.array(d)
    if np.any(d <= 0):
        return None, None
    na = d / d[-1]
    return na.tolist(), float(np.sqrt(np.mean((na - np.array(vnorm)) ** 2)))


print(f"device={t3.DEVICE}")
print("ViennaPS-3D trajectory (global-min depth):")
for v in vps:
    print(f"  d6={v['d6']:5.2f}  norm={np.round(v['norm'],3)}")
print()

result = {}
for be in BETAS:
    print(f"=== betaE={be} (adaptive bounces) ===")
    rows = []
    for r in RATES:
        cc = []; mm = []
        for dd in DIAMS:
            c, m = run(dd, be, r); cc.append(c); mm.append(m)
        rows.append(dict(rate=r, cen=cc, mx=mm))
        print(f"  rate={r:.3f}  center d={np.round(cc,2)}  max d={np.round(mm,2)}", flush=True)
    # for each VPS depth, nearest-d6 under the max metric (apples-to-apples), report rmse both metrics
    per = []
    for v in vps:
        cand = [x for x in rows if x['mx'][-1] > 0]
        if not cand:
            continue
        x = min(cand, key=lambda x: abs(x['mx'][-1] - v['d6']))
        cn, cr = norm_rmse(x['cen'], v['norm'])
        mn, mr = norm_rmse(x['mx'], v['norm'])
        per.append(dict(vps_d6=v['d6'], vps_norm=v['norm'], mx_d6=x['mx'][-1],
                        cen_norm=cn, cen_rmse=cr, mx_norm=mn, mx_rmse=mr))
        print(f"  vps d6={v['d6']:5.2f} {np.round(v['norm'],3)} | ours(max) d6={x['mx'][-1]:5.2f} "
              f"{np.round(mn,3) if mn else mn} rmse_max={mr if mr is None else round(mr,3)} "
              f"rmse_cen={cr if cr is None else round(cr,3)}")
    rmm = [p['mx_rmse'] for p in per if p['mx_rmse'] is not None]
    result[be] = dict(rows=rows, per=per, mean_rmse_max=float(np.mean(rmm)) if rmm else None)
    print(f"  -> betaE={be} mean ARDE rmse (max metric) across depths = {result[be]['mean_rmse_max']}\n", flush=True)

print("=" * 70)
best = min((b for b in BETAS if result[b]['mean_rmse_max'] is not None),
          key=lambda b: result[b]['mean_rmse_max'], default=None)
for be in BETAS:
    print(f"  betaE={be}: mean ARDE rmse (vs ViennaPS, max metric) = {result[be]['mean_rmse_max']}")
if best is not None:
    print(f"\n  BEST betaE={best}  mean rmse={result[best]['mean_rmse_max']:.3f}")
json.dump(dict(vps=vps, result={str(k): v for k, v in result.items()}, best=best),
          open("pixel_exact_3d_result.json", "w"), indent=2)
print("\nwrote pixel_exact_3d_result.json")
