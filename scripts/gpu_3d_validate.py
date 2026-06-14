#!/usr/bin/env python3
"""GPU clean validation of the 3D ARDE gap. Run on a Linux+NVIDIA box (PETCH_DEVICE=cuda).

Question: does coverage-dependent (Langmuir) neutral sticking close the 3D deep-HARC ARDE gap?

Flux transport runs on the GPU (Warp BVH). Level-set REINIT uses CPU skfmm: the GPU Russo-Smereka
reinit develops a NaN instability in deep holes (~step 16-20) that corrupts phi -> negative depths;
skfmm is rock-stable and the per-step cost is negligible vs the GPU flux solve. (GPU-reinit NaN bug
is logged separately as a production follow-up.)

Method (rate-MATCHED, compare shape not absolute rate):
  For coverage-sticking OFF and ON, sweep rate_scale and record the full d3/d4/d6 ARDE at each rate
  (rejecting any blowup: a depth <= 0). Then at MATCHED d6 depth compare the normalized ARDE curve
  (/d6) to ViennaPS-3D [0.832, 0.916, 1.0]. If coverage-sticking is the missing mechanism, at the
  SAME d6 its normalized ARDE is FLATTER (deep holes keep up) -> closer to ViennaPS.
"""
import os
os.environ["PETCH_DEVICE"] = "cuda"
import json
import numpy as np
import petch
from petch import threed as t3

VPS = {3.0: 7.684, 4.0: 8.469, 6.0: 9.241}      # ViennaPS-3D hole depths (um)
DIA = [3.0, 4.0, 6.0]
VPSN = np.array([VPS[d] for d in DIA]); VPSN = VPSN / VPSN[-1]

GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
DX, NS, NR = 0.25, 40, 30000
RATES = [0.05, 0.09, 0.15, 0.24, 0.36, 0.50]


def depth(dd, cov, rate):
    p = dict(petch.PAR); p['rate_scale'] = rate
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=cov),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    d = t3.center_depth_3d(g)
    return d if np.isfinite(d) else -1.0


print(f"Warp device check: {t3.DEVICE}")
print(f"ViennaPS-3D hole ARDE target: {np.round(VPSN, 3)}  (depths {[VPS[d] for d in DIA]})\n")

table = {False: [], True: []}   # cov -> list of dict(rate,d3,d4,d6,normARDE,rmse)
for cov in [False, True]:
    tag = "coverage ON " if cov else "coverage OFF"
    print(f"=== {tag} : rate sweep (d3/d4/d6, normalized ARDE vs VPS) ===")
    for r in RATES:
        d = np.array([depth(dd, cov, r) for dd in DIA])
        if np.any(d <= 0):
            print(f"  rate={r:.3f}  depths {np.round(d,2)}  [BLOWUP/invalid - skip]", flush=True)
            continue
        na = d / d[-1]
        rmse = float(np.sqrt(np.mean((na - VPSN) ** 2)))
        table[cov].append(dict(rate=r, d=d.tolist(), normARDE=na.tolist(), rmse=rmse, d6=float(d[-1])))
        print(f"  rate={r:.3f}  depths {np.round(d,2)}  normARDE {np.round(na,3)}  rmse={rmse:.4f}", flush=True)
    print()

# Matched-depth comparison: pick, per config, the row whose d6 is closest to a COMMON target =
# the deepest d6 BOTH configs can reach (so we compare ARDE shape at equal hole depth).
def maxd6(rows):
    return max((row['d6'] for row in rows), default=0.0)
common = min(maxd6(table[False]), maxd6(table[True]))
print("=" * 70)
print(f"MATCHED-DEPTH COMPARISON  (common reachable d6 ~ {common:.2f} um)")
print(f"  ViennaPS-3D      {np.round(VPSN,3)}")
picked = {}
for cov in [False, True]:
    rows = table[cov]
    if not rows:
        print(f"  {'coverage ON ' if cov else 'coverage OFF'}: no valid rows"); continue
    row = min(rows, key=lambda r: abs(r['d6'] - common))
    picked[cov] = row
    tag = "coverage ON " if cov else "coverage OFF"
    print(f"  {tag}     {np.round(row['normARDE'],3)}  rmse={row['rmse']:.4f}  "
          f"(rate={row['rate']:.3f}, d6={row['d6']:.2f})")
if len(picked) == 2:
    better = "coverage ON" if picked[True]['rmse'] < picked[False]['rmse'] else "coverage OFF"
    print(f"\n  -> {better} matches ViennaPS-3D ARDE shape better at matched d6 depth.")

with open("gpu_3d_validate_result.json", "w") as f:
    json.dump(dict(VPSN=VPSN.tolist(), DIA=DIA,
                   table={str(k): v for k, v in table.items()},
                   picked={str(k): v for k, v in picked.items()}), f, indent=2)
print("\nwrote gpu_3d_validate_result.json")
