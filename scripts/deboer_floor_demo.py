"""Demonstrate the bimodal-core floor: run the knee config with vs without the collimated ion core,
long enough (t_end=22) to reach high AR. WITHOUT core -> floor collapses, trench creeps and stalls
low. WITH core -> sustained ion-sputter floor, trench reaches high AR holding ~0.2. Saves
/root/deboer_floor.npz for local plotting. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0
GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=22.0)
DX, NS = 0.25, 300
BASE = dict(cal_F=1.5, betaE=0.85, ion_ang_sigma=np.deg2rad(0.8), Ysp_scale=15.0)


def run(over, seeds=(0, 1)):
    arN = np.linspace(0.5, 33, 66); acc = np.zeros_like(arN); cnt = np.zeros_like(arN); amax = 0.0
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1; p.update(over)
        fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol", warm_start_coverage=True)
        g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                           reinit_method="fsm", verbose=False, record_depth_every=5, seed_offset=sd * 100, **GEO)
        h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
        tm = st / NS * GEO['t_end']; r = np.gradient(dd, tm); ar = dd / W
        r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
        nr = np.clip(r / max(r0, 1e-9), 0, 1.6); amax = max(amax, ar.max())
        s = np.interp(arN, ar, nr, left=1.0, right=np.nan); ok = ~np.isnan(s); acc[ok] += s[ok]; cnt[ok] += 1
    return arN, np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan), amax


nocore = dict(BASE); nocore['ion_core_frac'] = 0.0
core = dict(BASE); core['ion_core_frac'] = 0.3; core['ion_core_sigma'] = np.deg2rad(0.3)
print(f"device={t3.DEVICE}  floor demo (knee config +/- bimodal ion core), t_end=22\n", flush=True)
ar_n, nr_n, am_n = run(nocore)
ar_c, nr_c, am_c = run(core)
np.savez("/root/deboer_floor.npz", ar_n=ar_n, nr_n=nr_n, ar_c=ar_c, nr_c=nr_c,
         exp_ar=EXP_AR, exp_r=EXP_R, am_n=am_n, am_c=am_c)
for lab, ar, nr, am in [("no core", ar_n, nr_n, am_n), ("+core", ar_c, nr_c, am_c)]:
    ok = ~np.isnan(nr)
    vals = np.interp([10, 20, 30], ar[ok], nr[ok], right=np.nan)
    print(f"  {lab:8s}: @AR10/20/30 = {np.round(vals,3)}   reached AR{am:.0f}", flush=True)
print(f"  de Boer  : @AR10/20/30 = {np.round(np.interp([10,20,30],[0,10,20,40],[1,.43,.29,.20]),3)}", flush=True)
print("  saved /root/deboer_floor.npz", flush=True)
