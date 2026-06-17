#!/usr/bin/env python3
"""Final de Boer match: dump full ARDE curves for (a) petch with VIENNAPS-DEFAULT params (etchant-rich,
broad IADF -> the ballistic gap) and (b) petch with DE-BOER PROCESS params (etchant-starved cal_F,
sub-degree IADF, lifted sputter floor -> matches the experiment). Seed-averaged to clean the MC floor
noise. Saves /root/deboer_final.npz for local plotting. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 160


def run_curve(par_over, seeds=(0, 1, 2)):
    arN = np.linspace(0.5, 30, 60)
    acc = np.zeros_like(arN); cnt = np.zeros_like(arN)
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1; p.update(par_over)
        fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol",
                         warm_start_coverage=True)
        g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                           reinit_method="fsm", verbose=False, record_depth_every=4, seed_offset=sd * 100, **GEO)
        h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
        tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
        r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
        nr = np.clip(r / max(r0, 1e-9), 0, 1.6)
        samp = np.interp(arN, ar, nr, left=1.0, right=np.nan)
        ok = ~np.isnan(samp); acc[ok] += samp[ok]; cnt[ok] += 1
    nr_avg = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return arN, nr_avg


DEFAULT = dict(cal_F=12.0, ion_ang_sigma=np.deg2rad(2.5))                              # ViennaPS regime
DEBOER = dict(cal_F=3.5, ion_ang_sigma=np.deg2rad(0.8), betaE=0.8, Ysp_scale=10.0)     # de Boer regime

print(f"device={t3.DEVICE}  computing seed-averaged curves...\n", flush=True)
ar_d, nr_d = run_curve(DEFAULT)
ar_b, nr_b = run_curve(DEBOER)
rb = np.interp(EXP_AR, ar_b, nr_b); rd = np.interp(EXP_AR, ar_d, nr_d)
rmse_b = float(np.sqrt(np.nanmean((rb - EXP_R) ** 2))); rmse_d = float(np.sqrt(np.nanmean((rd - EXP_R) ** 2)))
np.savez("/root/deboer_final.npz", ar_d=ar_d, nr_d=nr_d, ar_b=ar_b, nr_b=nr_b,
         exp_ar=EXP_AR, exp_r=EXP_R, rmse_b=rmse_b, rmse_d=rmse_d)
print(f"  petch DEFAULT (ViennaPS params)  @AR0/10/20/40 = {np.round(rd,3)}  RMSE {rmse_d:.3f}", flush=True)
print(f"  petch DE-BOER process params     @AR0/10/20/40 = {np.round(rb,3)}  RMSE {rmse_b:.3f}", flush=True)
print(f"  de Boer experiment                             = {EXP_R}", flush=True)
print("  saved /root/deboer_final.npz", flush=True)
