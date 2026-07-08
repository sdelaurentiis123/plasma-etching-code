#!/usr/bin/env python3
"""One arm of the de Boer AR>20 divergence campaign (2026-07-07): config x throttle x seed ->
raw depth-history npz. Curves are assembled by scripts/deboer_divergence_plot.py, which fixes the
r0 normalization bias (coarse-record long runs understate the AR<2 reference rate -- use a
FINE-cadence short run of the same config for r0).

Geometries: std (Lz=72, sub_top=66 -> AR ceiling 32) | tall (Lz=100, sub_top=94 -> AR ceiling 46;
GPU recommended). Physics identical across geometries.

Usage: deboer_divergence_arm.py <knee|default> <0|1 throttle> <seed> <t_end> <n_steps> <record_every> <outtag> [tall]
Campaign settings used: knee long = t_end 90, NS 360, rec 2 (tall); fine r0 runs = t_end 5 (knee) /
2 (default), rec 1. de Boer experiment: nr 1.0/0.43/0.29/0.20 @ AR 0/10/20/40."""
import os, sys, time
os.environ.setdefault("PETCH_DEVICE", "cpu")
import numpy as np, petch
from petch import threed as t3

cfg, thr, sd, t_end, NS, rec, tag = (sys.argv[1], int(sys.argv[2]), int(sys.argv[3]),
                                     float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
tall = len(sys.argv) > 8 and sys.argv[8] == "tall"
W = 2.0
GEO = (dict(Lx=10, Ly=5, Lz=100, mask_th=2, sub_top=94, hole=False) if tall
       else dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False))
OVER = dict(knee=dict(cal_F=1.5, ion_ang_sigma=np.deg2rad(0.8), betaE=0.85, Ysp_scale=10.0),
            default=dict(cal_F=12.0, ion_ang_sigma=np.deg2rad(2.5)))[cfg]
p = dict(petch.PAR); p['rate_scale'] = 0.30; p['periodic_y'] = 1; p.update(OVER)
fl = petch.Flags(coverage_sticking=True, neutral_transport="mc", sampling="sobol",
                 warm_start_coverage=True, floor_charge_throttle=bool(thr))
t0 = time.time()
g = t3.run_etch_3d(trench_width=W, dx=0.25, n_steps=NS, par=p, flags=fl, n_ion=40000, n_neu=40000,
                   reinit_method="fsm", verbose=False, record_depth_every=rec, seed_offset=sd*100,
                   t_end=t_end, **GEO)
h = g['depth_history']; st = np.array([x[0] for x in h], float); dd = np.array([x[1] for x in h], float)
np.savez(f"arm_{tag}.npz", t=st/NS*t_end, depth=dd, W=W, t_end=t_end)
print(f"ARMDONE {tag} maxdepth {dd.max():.2f} maxAR {dd.max()/W:.1f} wall {time.time()-t0:.0f}s", flush=True)
