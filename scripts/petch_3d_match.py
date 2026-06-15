#!/usr/bin/env python3
"""Match OUR 3D model to the depth-RESOLVED ViennaPS-3D ground truth (run on the box, PETCH_DEVICE=cuda).

Loads viennaps_3d_depth_resolved.json (ViennaPS holes d=3/4/6 at several durations -> depth + wall).
For our model (coverage_sticking ON, the validated config; skfmm reinit = trusted), sweeps rate_scale
to build OUR depth-resolved ARDE trajectory and per-etch wall-clock. Then:
  (1) ACCURACY: for each ViennaPS d6 depth, pick our nearest-d6 config and compare normalized ARDE
      [d3/d6, d4/d6] -> does our ARDE shape track ViennaPS ACROSS depths? Where our floor starves
      (can't reach ViennaPS's deep d6) is the quantified residual gap.
  (2) SPEED: our GPU per-etch wall-clock vs ViennaPS CPU wall-clock (same box) + our CPU time.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json, time
import numpy as np
import petch
from petch import threed as t3

GT = json.load(open("viennaps_3d_depth_resolved.json"))
DIAMS = [3.0, 4.0, 6.0]
# ViennaPS depth-resolved: dur -> [d3,d4,d6], and its normalized ARDE + wall per (dur,diam)
vps_traj = []  # list of dict(dur, d=[..], norm=[..], d6, wall6)
for dur in GT["durs"]:
    dd = np.array([GT["depth_grid"][str(dur)][str(x)] for x in DIAMS])
    w6 = GT["wall_grid"][str(dur)][str(6.0)]
    vps_traj.append(dict(dur=dur, d=dd.tolist(), norm=(dd/dd[-1]).tolist(), d6=float(dd[-1]), wall6=float(w6)))

GEO = dict(Lx=14, Ly=14, Lz=28, mask_th=2, sub_top=22, hole=True, t_end=3.0)  # deeper Lz/sub_top for deep holes
DX, NS, NR = 0.25, 40, 30000
RATES = [0.06, 0.10, 0.15, 0.22, 0.30, 0.40, 0.55, 0.75]


CALIB_BETA = float(os.environ.get("PETCH_BETAE", "0.08"))   # 3D-calibrated F sticking (beta sweep optimum)


def run_ours(dd, rate, cov, timeit=False):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = CALIB_BETA
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=dd, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=cov),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False, **GEO)
    wall = time.time() - t0
    d = t3.center_depth_3d(g)
    return (float(d) if np.isfinite(d) else -1.0, wall)


print(f"device={t3.DEVICE}")
print("ViennaPS-3D depth-resolved ARDE trajectory:")
for v in vps_traj:
    print(f"  dur={v['dur']:4.1f}  d={np.round(v['d'],2)}  normARDE={np.round(v['norm'],3)}  d6={v['d6']:.2f}  (vps wall {v['wall6']:.1f}s)")
print()

# Our trajectory: sweep rate, full d3/d4/d6 + wall, cov ON
ours = []   # dict(rate, d=[..], norm, d6, wall_gpu)
for r in RATES:
    drow = []; wall6 = 0.0
    for dd in DIAMS:
        dep, wall = run_ours(dd, r, cov=True, timeit=(dd == 6.0))
        drow.append(dep)
        if dd == 6.0:
            wall6 = wall
    d = np.array(drow)
    if np.any(d <= 0):
        print(f"  [ours cov-ON] rate={r:.3f}  depths {np.round(d,2)}  [invalid - skip]", flush=True); continue
    ours.append(dict(rate=r, d=d.tolist(), norm=(d/d[-1]).tolist(), d6=float(d[-1]), wall_gpu=wall6))
    print(f"  [ours cov-ON] rate={r:.3f}  d={np.round(d,2)}  normARDE={np.round(d/d[-1],3)}  d6={d[-1]:.2f}  (gpu {wall6:.2f}s)", flush=True)
print()

# Accuracy: match each ViennaPS d6 to our nearest-d6 config; compare normalized ARDE
print("=" * 72)
print("ACCURACY: our ARDE vs ViennaPS at matched d6 depth")
max_ours_d6 = max((o['d6'] for o in ours), default=0.0)
rows = []
for v in vps_traj:
    cand = [o for o in ours]
    if not cand:
        continue
    o = min(cand, key=lambda o: abs(o['d6'] - v['d6']))
    rmse = float(np.sqrt(np.mean((np.array(o['norm']) - np.array(v['norm']))**2)))
    reached = o['d6'] >= 0.85 * v['d6']
    rows.append(dict(vps_d6=v['d6'], ours_d6=o['d6'], vps_norm=v['norm'], ours_norm=o['norm'], rmse=rmse, reached=reached))
    flag = "" if reached else "  <-- our floor STARVES (cannot reach this depth)"
    print(f"  vps d6={v['d6']:5.2f} norm={np.round(v['norm'],3)} | ours d6={o['d6']:5.2f} norm={np.round(o['norm'],3)} rmse={rmse:.3f}{flag}")
print(f"\n  our max reachable d6 = {max_ours_d6:.2f} um  (ViennaPS reaches {vps_traj[-1]['d6']:.2f} um)")

# Speed: our GPU vs ViennaPS CPU (same box), at comparable d6
print("=" * 72)
print("SPEED: per-etch wall-clock (same box; ViennaPS=CPU_TRIANGLE, ours=GPU Warp)")
if ours:
    om = ours[len(ours)//2]
    vmatch = min(vps_traj, key=lambda v: abs(v['d6'] - om['d6']))
    print(f"  ours  (d6~{om['d6']:.1f}): {om['wall_gpu']:.2f}s GPU")
    print(f"  vps   (d6~{vmatch['d6']:.1f}): {vmatch['wall6']:.1f}s CPU")
    if om['wall_gpu'] > 0:
        print(f"  -> ours ~{vmatch['wall6']/max(om['wall_gpu'],1e-3):.1f}x faster (GPU vs CPU; not same-engine)")

json.dump(dict(vps_traj=vps_traj, ours=ours, accuracy_rows=rows), open("petch_3d_match_result.json", "w"), indent=2)
print("\nwrote petch_3d_match_result.json")
