#!/usr/bin/env python3
"""Validate on-device GPU source-ray generation (par['gpu_source']): does pseudorandom GPU gen hold
ACCURACY vs the Sobol-host source, and how much SPEED does killing the host gen + upload buy?

All configs: warm n_fp=1 + fsm reinit + GPU prep-cached smoothing (the current best). Only the source
differs. Compares the full depth TRAJECTORY (accuracy) + wall (speed). PETCH_DEVICE=cuda.
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np
import petch
from petch import threed as t3

DX, DIAM = 0.25, 6.0
GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)
NS, NRAY = 20, 30000


def run(label, gpu_source, sampling):
    p = dict(petch.PAR); p['n_fp'] = 1; p['gpu_source'] = gpu_source
    fl = petch.Flags(coverage_sticking=True, sampling=sampling, warm_start_coverage=True)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=NS, par=p, flags=fl,
                       n_ion=NRAY, n_neu=NRAY, reinit_method="fsm", verbose=False,
                       record_depth_every=2, **GEO)
    wall = time.time() - t0
    h = g['depth_history']
    st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    print(f"  {label:28s} wall {wall:5.2f}s  depth {dd[-1]:5.2f}um", flush=True)
    return st, dd, wall


print(f"device={t3.DEVICE}  hole d={DIAM} dx={DX} steps={NS} rays={NRAY}\n", flush=True)
# warm up JIT (both source paths + fsm + smooth)
_ = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=2, par=dict(petch.PAR, n_fp=1, gpu_source=True),
                   flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                   n_ion=3000, n_neu=3000, reinit_method="fsm", verbose=False, **GEO)

st0, base, w_sobol = run("Sobol host (current best)", False, "sobol")
st1, gpu, w_gpu = run("GPU source (pseudo)", True, "sobol")
st2, pse, w_pse = run("pseudo host (control)", False, "pseudo")

g_err = np.abs(np.interp(st0, st1, gpu) - base)
p_err = np.abs(np.interp(st0, st2, pse) - base)
print(f"\n  depth trajectory vs Sobol-host:")
print(f"    GPU source   : max {g_err.max():.2f}um  mean {g_err.mean():.2f}um   ({'accuracy-neutral' if g_err.max()<0.6 else 'CHECK'})")
print(f"    pseudo host  : max {p_err.max():.2f}um  mean {p_err.mean():.2f}um   (isolates QMC-vs-pseudo effect)")
print(f"  speed: GPU source {w_sobol/max(w_gpu,1e-3):.2f}x vs Sobol-host ({w_sobol:.2f}s -> {w_gpu:.2f}s)", flush=True)
json.dump(dict(sobol_wall=w_sobol, gpu_wall=w_gpu, pseudo_wall=w_pse,
               gpu_err_max=float(g_err.max()), gpu_err_mean=float(g_err.mean()),
               pseudo_err_max=float(p_err.max())), open("gpu_source_result.json", "w"), indent=2)
print("wrote gpu_source_result.json")
