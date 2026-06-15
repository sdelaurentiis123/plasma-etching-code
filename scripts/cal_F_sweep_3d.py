#!/usr/bin/env python3
"""3D flux-saturation reconciliation: sweep cal_F (the F-flux normalization in the coverage balance)
to test whether pushing the floor toward F-SATURATION un-starves our deep holes and flattens ARDE to
match ViennaPS-3D across depths. rate_scale only scales speed; cal_F sets floor theta_F -> the actual
ARDE/starvation lever (Coburn-Winters: as the floor stays F-fed, R_b/R_t -> 1/K, ARDE flattens).

For each cal_F, sweep rate; record d3/d4/d6 + normalized ARDE; find the config whose d6 is closest to
ViennaPS d6=9.24 um (dur=3) and report its ARDE rmse vs [0.831,0.917,1.0]. The win = a cal_F that both
REACHES 9.24 um AND flattens ARDE toward ViennaPS. Run on the box (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

VPS_D6 = 9.241
VPS_NORM = np.array([0.831, 0.917, 1.0])     # ViennaPS dur=3 hole ARDE
DIAMS = [3.0, 4.0, 6.0]
GEO = dict(Lx=14, Ly=14, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 40, 30000
CALFS = [12.0, 35.0, 70.0, 130.0]
RATES = [0.05, 0.10, 0.18, 0.30, 0.45]


def depth(dd, cal_F, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['cal_F'] = cal_F
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    d = t3.center_depth_3d(g)
    return float(d) if np.isfinite(d) else -1.0


print(f"device={t3.DEVICE}  target ViennaPS d6={VPS_D6} norm={VPS_NORM}\n")
best_overall = None
summary = {}
for cf in CALFS:
    print(f"=== cal_F={cf} ===")
    rows = []
    for r in RATES:
        d = np.array([depth(dd, cf, r) for dd in DIAMS])
        if np.any(d <= 0):
            print(f"  rate={r:.3f}  {np.round(d,2)}  [invalid]", flush=True); continue
        na = d / d[-1]
        rows.append(dict(rate=r, d=d.tolist(), norm=na.tolist(), d6=float(d[-1])))
        print(f"  rate={r:.3f}  d={np.round(d,2)}  normARDE={np.round(na,3)}  d6={d[-1]:.2f}", flush=True)
    if not rows:
        summary[cf] = None; continue
    # config closest to ViennaPS d6
    o = min(rows, key=lambda x: abs(x['d6'] - VPS_D6))
    rmse = float(np.sqrt(np.mean((np.array(o['norm']) - VPS_NORM) ** 2)))
    reached = o['d6'] >= 0.9 * VPS_D6
    summary[cf] = dict(best=o, rmse=rmse, reached=reached, max_d6=max(x['d6'] for x in rows))
    print(f"  -> nearest d6={o['d6']:.2f} (reached={reached}) normARDE={np.round(o['norm'],3)} rmse-vs-VPS={rmse:.3f}  max_d6={summary[cf]['max_d6']:.2f}\n", flush=True)
    if best_overall is None or (reached and rmse < best_overall[1]):
        best_overall = (cf, rmse, o)

print("=" * 70)
print(f"ViennaPS dur=3:  d6={VPS_D6}  norm={VPS_NORM}")
for cf in CALFS:
    s = summary.get(cf)
    if s:
        print(f"  cal_F={cf:6.1f}: max_d6={s['max_d6']:5.2f}  nearest-d6 norm={np.round(s['best']['norm'],3)}  rmse={s['rmse']:.3f}  reached={s['reached']}")
if best_overall:
    cf, rmse, o = best_overall
    print(f"\n  BEST (reaches depth + flattest): cal_F={cf}  d6={o['d6']:.2f}  norm={np.round(o['norm'],3)}  rmse={rmse:.3f}")
else:
    print("\n  NO cal_F both reached 9.24um AND matched ARDE -> deeper model limit (charging/transport).")

json.dump(dict(VPS_D6=VPS_D6, VPS_NORM=VPS_NORM.tolist(),
               summary={str(k): v for k, v in summary.items()}), open("cal_F_sweep_3d_result.json", "w"), indent=2)
print("\nwrote cal_F_sweep_3d_result.json")
