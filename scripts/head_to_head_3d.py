#!/usr/bin/env python3
"""Honest SAME-ENGINE speed benchmark: our differentiable Warp etcher (GPU) vs ViennaPS-3D, on a
MATCHED 3D hole etch, measuring wall-clock. Tries ViennaPS GPU_TRIANGLE (OptiX) for a fair GPU-vs-GPU
comparison; falls back to CPU_TRIANGLE and LABELS which engine ran. Run on a box with a recent driver
so OptiX is available (PETCH_DEVICE=cuda)."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import time, json
import numpy as np

# ---- ViennaPS side ----
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
DX, XE, MASK, DUR, DIAM = 0.25, 14.0, 2.0, 3.0, 6.0

def vps_run(engine):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=DIAM/2,
                holeDepth=MASK, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
    p.setFluxEngineType(engine)
    t0 = time.time(); p.apply(); wall = time.time() - t0
    n = np.array(d.getSurfaceMesh().getNodes())
    return wall, float(-n[:, 2].min())

vps_engine = None
for name in ["GPU_TRIANGLE", "CPU_TRIANGLE"]:
    eng = getattr(ps.FluxEngineType, name, None)
    if eng is None:
        continue
    try:
        w, dep = vps_run(eng)
        vps_engine, vps_wall, vps_depth = name, w, dep
        print(f"ViennaPS {name}: {w:.2f}s  depth={dep:.2f}um", flush=True)
        break
    except Exception as e:
        print(f"ViennaPS {name} failed: {str(e)[:70]}", flush=True)

# ---- our side (GPU Warp), matched depth ----
import petch
from petch import threed as t3
# rate-match our d6 to ViennaPS depth
GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
def ours_run(rate, ns=40):
    p = dict(petch.PAR); p['rate_scale'] = rate; p['betaE'] = 0.7
    p.update(flux_smooth_gpu=True, gpu_source=True, gpu_mesh=True, gpu_warmstart=True, device_flux=True, n_fp=1)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=DIAM, dx=DX, n_steps=ns, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol", warm_start_coverage=True),
                       n_ion=30000, n_neu=30000, reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.center_depth_3d(g)

ours_run(0.10, ns=4)   # warmup (JIT compile all GPU kernels)
# quick rate find for ~vps_depth
best = None
for r in [0.04, 0.07, 0.10, 0.14]:
    w, dep = ours_run(r)
    if best is None or abs(dep - vps_depth) < abs(best[2] - vps_depth):
        best = (r, w, dep)
ours_wall, ours_depth = best[1], best[2]
print(f"ours (Warp GPU, RR+smoothing): {ours_wall:.2f}s  depth={ours_depth:.2f}um  (rate={best[0]})", flush=True)

print("=" * 56)
print(f"ViennaPS ({vps_engine}): {vps_wall:.2f}s @ {vps_depth:.2f}um")
print(f"ours     (Warp GPU)   : {ours_wall:.2f}s @ {ours_depth:.2f}um")
print(f"-> ours {vps_wall/max(ours_wall,1e-3):.1f}x faster than ViennaPS-{vps_engine}"
      + ("  (SAME-ENGINE GPU-vs-GPU)" if vps_engine == "GPU_TRIANGLE" else "  (ours-GPU vs ViennaPS-CPU)"))
json.dump(dict(vps_engine=vps_engine, vps_wall=vps_wall, vps_depth=vps_depth,
               ours_wall=ours_wall, ours_depth=ours_depth,
               speedup=vps_wall/max(ours_wall,1e-3)), open("head_to_head_result.json", "w"), indent=2)
print("wrote head_to_head_result.json")
