#!/usr/bin/env python3
"""FAIR head-to-head SWEEP: ours (full GPU stack) vs ViennaPS-GPU (OptiX RT-core) across hole diameters
(aspect ratios), MATCHED depth. NOTE: ViennaPS-OptiX and Warp-CUDA conflict if interleaved, so we run
ALL ViennaPS first (before importing petch/warp), then all ours. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np

DX, DUR, XE = 0.25, 3.0, 14.0
DIAMS = [4.0, 6.0, 8.0]

# ---------- PHASE 1: all ViennaPS (no warp yet) ----------
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def vps_run(diam):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=diam/2,
                holeDepth=2.0, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE)
    t0 = time.time(); p.apply(); wall = time.time() - t0
    n = np.array(d.getSurfaceMesh().getNodes())
    return wall, float(-n[:, 2].min())


vps = {}
for diam in DIAMS:
    w, dep = vps_run(diam)                    # ViennaPS-GPU is steady (no cold-start), one measure
    vps[diam] = (w, dep)
    print(f"  ViennaPS-GPU d={diam}: {w:.2f}s depth {dep:.2f}um", flush=True)

# ---------- PHASE 2: all ours (warp) ----------
import petch
from petch import threed as t3


def ours_run(diam, rate, ns=40):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = 0.7
    p.update(flux_smooth_gpu=True, gpu_source=True, gpu_mesh=True, gpu_warmstart=True, device_flux=True, n_fp=1)
    GEO = dict(Lx=XE, Ly=XE, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=diam, dx=DX, n_steps=ns, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.center_depth_3d(g)


ours_run(6.0, 0.08, ns=4)                     # warm ours (JIT)
print(f"\n{'d(um)':>6} {'AR':>5} | {'ViennaPS-GPU':>13} {'vdep':>6} | {'ours':>7} {'odep':>6} | {'speedup':>8}", flush=True)
res = []
for diam in DIAMS:
    vw, vdep = vps[diam]
    best = None
    for r in [0.05, 0.08, 0.12]:
        ow, odep = ours_run(diam, r)
        if best is None or abs(odep - vdep) < abs(best[1] - vdep):
            best = (ow, odep, r)
    ow, odep, r = best
    sx = vw / max(ow, 1e-3)
    res.append(dict(diam=diam, ar=round(vdep/diam, 1), vps_wall=vw, vps_depth=vdep,
                    ours_wall=ow, ours_depth=odep, speedup=sx))
    print(f"{diam:6.1f} {vdep/diam:5.1f} | {vw:11.2f}s {vdep:6.2f} | {ow:6.2f}s {odep:6.2f} | {sx:7.1f}x", flush=True)

sp = sorted(r['speedup'] for r in res)
print(f"\n  ours vs ViennaPS-GPU (OptiX RT-core), matched depth: {sp[0]:.1f}x - {sp[-1]:.1f}x (median {sp[len(sp)//2]:.1f}x)", flush=True)
json.dump(res, open("vps_sweep_result.json", "w"), indent=2)
print("wrote vps_sweep_result.json")
