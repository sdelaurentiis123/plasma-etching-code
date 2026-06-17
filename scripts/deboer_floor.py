#!/usr/bin/env python3
"""Sustain the de Boer FLOOR with a bimodal ion IADF (opt-in core, beyond ViennaPS). Take the knee
config (cal_F=1.5, betaE=0.85, sig=0.8) that matched AR0-10 but STALLED at AR~19, and add a collimated
ion CORE (ion_core_frac of ions at ion_core_sigma << 0.7deg, narrow enough to reach AR40). The core
ions deliver sputter flux to the bottom -> sustained floor -> trench no longer stalls. Goal: reach
AR>30 with a ~0.20 floor, matching the experiment's 1.0/0.43/0.29/0.20 @ AR 0/10/20/40.

REGRESSION GUARD: also prints the DEFAULT (ViennaPS-regime, no core) trench depth -- must be unchanged
(ion_core_frac=0 -> exact single-Gaussian source). PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 30.0, 40.0]); EXP_R = np.interp(EXP_AR, [0, 10, 20, 40], [1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=14.0)
DX, NS = 0.25, 200
KNEE = dict(cal_F=1.5, betaE=0.85, ion_ang_sigma=np.deg2rad(0.8), Ysp_scale=10.0)


def run(over, seeds=(0, 1)):
    arN = np.linspace(0.5, 33, 66); acc = np.zeros_like(arN); cnt = np.zeros_like(arN); amax = 0.0
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1; p.update(over)
        fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol", warm_start_coverage=True)
        g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                           reinit_method="fsm", verbose=False, record_depth_every=4, seed_offset=sd * 100, **GEO)
        h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
        tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
        r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
        nr = np.clip(r / max(r0, 1e-9), 0, 1.6); amax = max(amax, ar.max())
        s = np.interp(arN, ar, nr, left=1.0, right=np.nan); ok = ~np.isnan(s); acc[ok] += s[ok]; cnt[ok] += 1
    nr = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return arN, nr, amax


# regression guard: default ViennaPS-regime config (no core) -- depth must match the known-good path
pg = dict(petch.PAR); pg['rate_scale'] = 0.30; pg['periodic_y'] = 1
flg = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol", warm_start_coverage=True)
gg = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=80, par=pg, flags=flg, n_ion=40000, n_neu=40000,
                    reinit_method="fsm", verbose=False, Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=6.0)
print(f"device={t3.DEVICE}", flush=True)
print(f"  REGRESSION GUARD: default config (no core) depth @ t=6 = {t3.center_depth_3d(gg):.3f} um  "
      f"(ion_core_frac=0 -> single-Gaussian source unchanged)\n", flush=True)

print(f"  de Boer {np.round(EXP_R,3)} @ AR {EXP_AR}   (knee config + bimodal ion core)\n", flush=True)
print(f"  {'cfrac':>6} {'csig':>5}   {'@AR 10/20/30/40':>26}   {'RMSE':>6}  AR_max", flush=True)
best = None
for cfrac in [0.0, 0.15, 0.30, 0.5]:
    for csig in [0.3]:
        over = dict(KNEE); over['ion_core_frac'] = cfrac; over['ion_core_sigma'] = np.deg2rad(csig)
        arN, nr, amax = run(over)
        ok = ~np.isnan(nr)
        ours = np.interp(EXP_AR[1:], arN[ok], nr[ok], right=np.nan)
        valid = ~np.isnan(ours)
        rmse = float(np.sqrt(np.mean((ours[valid] - EXP_R[1:][valid]) ** 2))) if valid.any() else 9.9
        tag = "  <--" if (best is None or rmse < best[-1]) else ""
        print(f"  {cfrac:6.2f} {csig:5.1f}   {np.round(ours,3)!s:>26}   {rmse:6.3f}  {amax:.1f}{tag}", flush=True)
        if best is None or rmse < best[-1]:
            best = (cfrac, csig, ours, amax, rmse)
print(f"\n  BEST core_frac={best[0]} core_sig={best[1]}deg -> {np.round(best[2],3)}  RMSE={best[4]:.3f}  reached AR{best[3]:.0f}", flush=True)
print(f"  de Boer @AR10/20/30/40                       {np.round(EXP_R[1:],3)}", flush=True)
print(f"  (cfrac=0 row = the stalling knee baseline; cfrac>0 should un-stall + sustain the floor)", flush=True)
