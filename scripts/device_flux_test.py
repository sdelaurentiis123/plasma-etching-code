#!/usr/bin/env python3
"""device_flux (normalize+smooth on GPU, 1 readback/flux): depth parity (within-noise) + speed, full stack."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time
import petch
from petch import threed as t3

GEO = dict(Lx=14, Ly=14, Lz=34, mask_th=2, sub_top=28, hole=True, t_end=1.2)


def run(devflux):
    p = dict(petch.PAR); p["n_fp"] = 1
    p["flux_smooth_gpu"] = True; p["gpu_source"] = True; p["gpu_mesh"] = True
    p["gpu_warmstart"] = True; p["device_flux"] = devflux
    fl = petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=20, par=p, flags=fl,
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.max_depth_3d(g), g["timings"]


run(True)  # warmup
wh, dh, tmh = run(False)
wd, dd, tmd = run(True)
print(f"host-flux   : wall {wh:.2f}s  flux {tmh['flux']:.2f}s  depth {dh:.2f}um")
print(f"device-flux : wall {wd:.2f}s  flux {tmd['flux']:.2f}s  depth {dd:.2f}um")
tag = "within-noise" if abs(dh - dd) < 0.6 else "CHECK"
print(f"-> {wh/max(wd,1e-3):.2f}x faster overall, flux {tmh['flux']/max(tmd['flux'],1e-3):.2f}x; depth delta {abs(dh-dd):.2f}um ({tag})")
