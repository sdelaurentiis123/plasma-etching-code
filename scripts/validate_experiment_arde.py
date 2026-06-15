#!/usr/bin/env python3
"""REAL-WAFER validation: our model's ARDE vs the de Boer 2002 / Blauw cryo SF6/O2 RIE-lag EXPERIMENT
(not vs ViennaPS). Blauw Knudsen/Clausing F-transport-limited model fit to de Boer's wafers gives the
normalized etch-rate-vs-aspect-ratio curve. Experimental anchor points (digitized): normalized rate
~1.0 / 0.43 / 0.29 / 0.20 at AR ~0 / 10 / 20 / 40 (F reaction prob S_F~0.47).

Method: etch ONE deep trench (the de Boer feature is a trench), record depth per step -> instantaneous
AR=depth/width and rate=d(depth)/dt; normalize rate to the shallow (AR~1) value; compare to the
experimental points. Run on a box (PETCH_DEVICE=cuda). Honest: our model is calibrated to ViennaPS
(which is calibrated to a DIFFERENT experiment, Belen), so this is a genuine cross-experiment test."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import json
import numpy as np
import petch
from petch import threed as t3

# Blauw F-transport-limited model R_b/R_t = K^AR/(K^AR + S_F - K^AR*S_F), Clausing K ~ 1/(1+AR/2)
S_F = 0.47
def blauw(ar):
    K = 1.0 / (1.0 + 0.5 * np.asarray(ar, float))   # simple Clausing-like transmission
    return K / (K + S_F - K * S_F)
EXP_AR = np.array([0.0, 10.0, 20.0, 40.0])
EXP_R = np.array([1.0, 0.43, 0.29, 0.20])           # de Boer/Blauw digitized

W = 2.0                                              # trench width (um) -> AR=depth/W
GEO = dict(Lx=10, Ly=5, Lz=60, mask_th=2, sub_top=54, hole=False, t_end=8.0)
DX, NS, NR = 0.25, 120, 30000


def run():
    p = dict(petch.PAR); p['rate_scale'] = 0.18; p['betaE'] = 0.7   # validated 3D config
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=NR, n_neu=NR, reinit_method="skfmm", verbose=False,
                       record_depth_every=6, **GEO)
    hist = g['depth_history']
    for k, d in hist:
        print(f"  step {k:3d}  depth={d:6.2f}  AR={d/W:5.2f}", flush=True)
    return hist


print(f"device={t3.DEVICE}  REAL-WAFER ARDE validation vs de Boer/Blauw (trench W={W}um)")
print(f"  experiment: normARDE {EXP_R} at AR {EXP_AR}\n")
dep = run()
steps = np.array([d[0] for d in dep]); dd = np.array([d[1] for d in dep])
t = steps / NS * GEO['t_end']
# instantaneous rate from finite differences; AR at each point
rate = np.gradient(dd, t)
ar = dd / W
# normalize rate to the shallowest (AR<2) point
ref = rate[ar < 2.0]
r0 = float(ref.max()) if len(ref) else float(rate.max())
nrate = np.clip(rate / max(r0, 1e-9), 0, 2)
print("\n  ours (AR, normRate):")
for a, r in zip(ar, nrate):
    print(f"    AR={a:5.2f}  normRate={r:.3f}  (Blauw {blauw(a):.3f})")
# compare at the experimental AR points (interpolate ours)
ours_at = np.interp(EXP_AR, ar, nrate)
rmse = float(np.sqrt(np.mean((ours_at - EXP_R) ** 2)))
print(f"\n  ours at exp AR {EXP_AR}: {np.round(ours_at,3)}")
print(f"  experiment           : {EXP_R}")
print(f"  RMSE vs de Boer/Blauw experiment = {rmse:.3f}")
json.dump(dict(exp_ar=EXP_AR.tolist(), exp_r=EXP_R.tolist(), ar=ar.tolist(),
               nrate=nrate.tolist(), ours_at=ours_at.tolist(), rmse=rmse),
          open("validate_experiment_arde_result.json", "w"), indent=2)
print("wrote validate_experiment_arde_result.json")
