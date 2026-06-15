#!/usr/bin/env python3
"""3D ARDE via the STICKING coefficient (the Coburn-Winters lever). Lower betaE -> radicals reflect
off the (saturated) walls and penetrate deeper to the under-fed floor -> R_b/R_t -> 1/K (floor keeps
etching) AND flatter ARDE. Literature backs low F sticking (Donnelly 2017: 0.001-0.03 at high flux).

For each betaE, sweep rate; find the config whose d6 is nearest ViennaPS d6=9.24 um; report whether it
REACHES depth and its ARDE rmse vs [0.831,0.917,1.0]. The win = a betaE that reaches 9.24um AND flattens
ARDE. If even betaE->0.03 can't -> the floor-starvation is a genuine transport model-form limit (bounce
model / resolution), not a calibration gap. Run on the box (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

VPS_D6 = 9.241
VPS_NORM = np.array([0.831, 0.917, 1.0])
DIAMS = [3.0, 4.0, 6.0]
GEO = dict(Lx=14, Ly=14, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 40, 30000
BETAS = [0.7, 0.4, 0.2, 0.08, 0.03]
RATES = [0.05, 0.10, 0.18, 0.30, 0.45]


def depth(dd, beta, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = beta
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    d = t3.center_depth_3d(g)
    return float(d) if np.isfinite(d) else -1.0


print(f"device={t3.DEVICE}  target ViennaPS d6={VPS_D6} norm={VPS_NORM}\n")
best = None; summary = {}
for be in BETAS:
    print(f"=== betaE={be} ===")
    rows = []
    for r in RATES:
        d = np.array([depth(dd, be, r) for dd in DIAMS])
        if np.any(d <= 0):
            print(f"  rate={r:.3f}  {np.round(d,2)}  [invalid]", flush=True); continue
        na = d / d[-1]
        rows.append(dict(rate=r, d=d.tolist(), norm=na.tolist(), d6=float(d[-1])))
        print(f"  rate={r:.3f}  d={np.round(d,2)}  normARDE={np.round(na,3)}  d6={d[-1]:.2f}", flush=True)
    if not rows:
        summary[be] = None; continue
    o = min(rows, key=lambda x: abs(x['d6'] - VPS_D6))
    rmse = float(np.sqrt(np.mean((np.array(o['norm']) - VPS_NORM) ** 2)))
    reached = o['d6'] >= 0.9 * VPS_D6
    summary[be] = dict(best=o, rmse=rmse, reached=reached, max_d6=max(x['d6'] for x in rows))
    print(f"  -> nearest d6={o['d6']:.2f} (reached={reached}) norm={np.round(o['norm'],3)} rmse={rmse:.3f}  max_d6={summary[be]['max_d6']:.2f}\n", flush=True)
    if reached and (best is None or rmse < best[1]):
        best = (be, rmse, o)

print("=" * 70)
print(f"ViennaPS dur=3:  d6={VPS_D6}  norm={VPS_NORM}")
for be in BETAS:
    s = summary.get(be)
    if s:
        print(f"  betaE={be:5.2f}: max_d6={s['max_d6']:5.2f}  nearest-d6 norm={np.round(s['best']['norm'],3)}  rmse={s['rmse']:.3f}  reached={s['reached']}")
if best:
    be, rmse, o = best
    print(f"\n  BEST (reaches depth + flattest): betaE={be}  d6={o['d6']:.2f}  norm={np.round(o['norm'],3)}  rmse={rmse:.3f}")
else:
    md = max((s['max_d6'] for s in summary.values() if s), default=0)
    print(f"\n  NO betaE reached 9.24um (max reached={md:.2f}um) -> floor-starvation is a TRANSPORT MODEL-FORM")
    print("  limit (bounce model / conductance), NOT a calibration gap. Matches the literature: deep-HARC")
    print("  needs Knudsen molecular-flow transport (+ charging), the physics ViennaPS approximates too.")

json.dump(dict(VPS_D6=VPS_D6, VPS_NORM=VPS_NORM.tolist(),
               summary={str(k): v for k, v in summary.items()}), open("beta_sweep_3d_result.json", "w"), indent=2)
print("\nwrote beta_sweep_3d_result.json")
