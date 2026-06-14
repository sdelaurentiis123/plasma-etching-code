#!/usr/bin/env python3
"""ViennaPS-3D ground truth for the 3D calibration. Run on a Linux+NVIDIA box.

Runs SF6O2 3D trenches (reliable in d3) at several widths -> depth + ARDE, with the CPU_TRIANGLE
(Embree) engine. Also attempts MakeHole for a genuine 3D contact-hole reference. Saves JSON.
"""
import json
import time
import numpy as np
import viennaps as ps
import viennaps.d3 as v3

ps.Logger.setLogLevel(ps.LogLevel.ERROR)
ps.Length.setUnit("micrometer")
ps.Time.setUnit("min")
Mat = ps.Material

DX, XE, YE, MASK, DUR = 0.3, 12.0, 6.0, 2.0, 3.0
WIDTHS = [3.0, 4.0, 6.0, 8.0]


def etch(make_fn):
    d = v3.Domain()
    make_fn(d)
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE)
    p.apply()
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min())


print("gpuAvailable:", ps.gpuAvailable(), "| version:", getattr(ps, "__version__", "?"))
print("d3 makers:", [x for x in dir(v3) if x.startswith("Make")])

# --- trenches (reliable) ---
trench_depths = {}
for w in WIDTHS:
    def mk(d, w=w):
        v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=w,
                      trenchDepth=MASK, taperingAngle=0.0, baseHeight=0.0,
                      periodicBoundary=False, makeMask=True, material=Mat.Si).apply()
    t0 = time.time(); dep = etch(mk)
    trench_depths[str(w)] = dep
    print(f"  trench w={w}: depth={dep:.3f} um  ({time.time()-t0:.1f}s)")

# --- holes (genuine 3D; attempt, may not match this API) ---
hole_depths = {}
for r in [1.5, 2.0, 3.0]:
    try:
        def mk(d, r=r):
            v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=r,
                        holeDepth=MASK, makeMask=True, material=Mat.Si).apply()
        dep = etch(mk)
        hole_depths[str(2 * r)] = dep
        print(f"  hole d={2*r}: depth={dep:.3f} um")
    except Exception as e:
        print(f"  hole r={r}: skipped ({str(e)[:80]})")
        break

out = dict(DX=DX, XE=XE, YE=YE, MASK=MASK, DUR=DUR, widths=WIDTHS,
           trench_depths=trench_depths, hole_depths=hole_depths)
json.dump(out, open("viennaps_3d_groundtruth.json", "w"), indent=2)
print("wrote viennaps_3d_groundtruth.json")
