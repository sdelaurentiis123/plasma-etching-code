#!/usr/bin/env python3
"""Depth-RESOLVED ViennaPS-3D ground truth: SF6O2 contact holes d=3/4/6 run to SEVERAL durations,
so we get the ARDE-vs-depth TRAJECTORY (not a single point). The real parity test for our model is
matching how the d3/d6 ratio evolves as the holes deepen. Also records wall-clock per etch (for the
head-to-head speed benchmark) and tries the GPU engine if the driver supports it. Run on a Linux+NVIDIA box.
"""
import json, time, os
import numpy as np
import viennaps as ps
import viennaps.d3 as v3

ps.Logger.setLogLevel(ps.LogLevel.ERROR)
ps.Length.setUnit("micrometer")
ps.Time.setUnit("min")
Mat = ps.Material

DX, XE, MASK = 0.3, 14.0, 2.0
DIAMS = [3.0, 4.0, 6.0]
DURS = [1.0, 2.0, 3.0, 4.5]          # depth-resolved: shallow -> deep

# pick the fastest working flux engine. NOTE: GPU_TRIANGLE is skipped by default because OptiX
# (libnvoptix.so) is frequently absent in vast.ai containers and a failed optixInit() HANGS apply()
# rather than raising — set PETCH_TRY_GPU=1 only on a box where ps.gpuAvailable() is truly True.
def pick_engine():
    names = (["GPU_TRIANGLE"] if os.environ.get("PETCH_TRY_GPU") == "1" else []) + ["CPU_TRIANGLE", "CPU_DISK"]
    for name in names:
        eng = getattr(ps.FluxEngineType, name, None)
        if eng is None:
            continue
        try:
            d = v3.Domain()
            v3.MakeHole(domain=d, gridDelta=0.6, xExtent=6.0, yExtent=6.0, holeRadius=1.0,
                        holeDepth=MASK, makeMask=True, material=Mat.Si).apply()
            m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
            p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(0.5)
            p.setFluxEngineType(eng); p.apply()
            print(f"engine OK: {name}")
            return eng, name
        except Exception as e:
            print(f"engine {name} failed: {str(e)[:70]}")
    raise RuntimeError("no working flux engine")


def etch_hole(diam, dur, engine):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=diam / 2.0,
                holeDepth=MASK, makeMask=True, material=Mat.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(engine)
    t0 = time.time(); p.apply(); wall = time.time() - t0
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min()), wall


engine, ename = pick_engine()
print(f"using engine={ename}\n")   # NB: do NOT call ps.gpuAvailable() — it triggers a fatal OptiX probe
grid = {}     # dur -> {diam -> depth}
walls = {}    # dur -> {diam -> seconds}
for dur in DURS:
    grid[dur] = {}; walls[dur] = {}
    for dd in DIAMS:
        dep, wall = etch_hole(dd, dur, engine)
        grid[dur][dd] = dep; walls[dur][dd] = wall
        print(f"  dur={dur:4.1f} d={dd}: depth={dep:6.3f} um  ({wall:5.1f}s)", flush=True)
    dvec = np.array([grid[dur][d] for d in DIAMS])
    print(f"    -> normARDE {np.round(dvec/dvec[-1],3)}\n", flush=True)

out = dict(engine=ename, DX=DX, XE=XE, MASK=MASK, diams=DIAMS, durs=DURS,
           depth_grid={str(k): {str(kk): vv for kk, vv in v.items()} for k, v in grid.items()},
           wall_grid={str(k): {str(kk): vv for kk, vv in v.items()} for k, v in walls.items()})
json.dump(out, open("viennaps_3d_depth_resolved.json", "w"), indent=2)
print("wrote viennaps_3d_depth_resolved.json")
