#!/usr/bin/env python3
"""GPU warm-start (mesh_query_point vs scipy KDTree): depth parity (within-noise) + speed, full stack."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time
import petch
from petch import threed as t3

GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)


def run(gws):
    p = dict(petch.PAR); p["n_fp"] = 1
    p["flux_smooth_gpu"] = True; p["gpu_source"] = True; p["gpu_mesh"] = True; p["gpu_warmstart"] = gws
    fl = petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=20, par=p, flags=fl,
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.max_depth_3d(g), g["timings"]


run(True)  # warmup
wk, dk, tmk = run(False)
wg, dg, tmg = run(True)
print(f"KDTree warm-start : wall {wk:.2f}s  flux {tmk['flux']:.2f}s  depth {dk:.2f}um")
print(f"GPU    warm-start : wall {wg:.2f}s  flux {tmg['flux']:.2f}s  depth {dg:.2f}um")
tag = "within-noise" if abs(dk - dg) < 0.6 else "CHECK"
print(f"-> {wk/max(wg,1e-3):.2f}x faster overall, flux {tmk['flux']/max(tmg['flux'],1e-3):.2f}x; depth delta {abs(dk-dg):.2f}um ({tag})")
