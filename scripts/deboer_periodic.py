#!/usr/bin/env python3
"""Does PERIODIC-Y fix the trench ARDE Craig plotted? de Boer/Blauw cryo SF6/O2 RIE-lag: normalized
bottom etch rate vs aspect ratio = 1.0/0.43/0.29/0.20 at AR 0/10/20/40. A trench is invariant in y;
without periodic-y rays leak out the open ends -> radiosity flat-then-cliff, MC erratic. Compare
periodic_y 0 vs 1 for MC and radiosity. Run on a box: PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0])
EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 140


def curve(nt, periodic, beta=0.47):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta; p['periodic_y'] = periodic
    fl = petch.Flags(coverage_sticking=True, neutral_transport=nt, sampling="sobol", warm_start_coverage=True)
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False,
                       record_depth_every=7, **GEO)
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    t = st / NS * GEO['t_end']; r = np.gradient(dd, t); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    nr = np.clip(r / max(r0, 1e-9), 0, 2)
    ours = np.interp(EXP_AR, ar, nr)
    rmse = float(np.sqrt(np.mean((ours - EXP_R) ** 2)))
    return float(ar.max()), ours, rmse


print(f"device={t3.DEVICE}  de Boer trench: normARDE {EXP_R} at AR {EXP_AR}\n", flush=True)
print(f"  {'config':28s} {'AR_max':>6} {'@AR 0/10/20/40':>26} {'RMSE':>6}", flush=True)
res = {}
for nt in ["mc", "radiosity"]:
    for periodic in [0, 1]:
        armax, ours, rmse = curve(nt, periodic)
        key = f"{nt}_periodic{periodic}"
        res[key] = dict(ar_max=armax, ours=ours.tolist(), rmse=rmse)
        print(f"  {nt+' periodic_y='+str(periodic):28s} {armax:6.1f} {np.round(ours,3)!s:>26} {rmse:6.3f}", flush=True)
print("\n  periodic_y=1 should: radiosity roll off gradually (not flat-then-cliff), MC less erratic,", flush=True)
print("  both lower RMSE vs the de Boer experiment.", flush=True)
json.dump(dict(exp_ar=EXP_AR.tolist(), exp_r=EXP_R.tolist(), result=res),
          open("deboer_periodic_result.json", "w"), indent=2)
print("wrote deboer_periodic_result.json")
