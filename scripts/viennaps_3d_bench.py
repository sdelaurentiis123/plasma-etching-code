#!/usr/bin/env python3
"""ViennaPS 3D timing: does the GPU/OptiX engine work on this driver, and how fast vs CPU?
Run on a Linux+NVIDIA box. Times an SF6O2 3D trench etch with GPU_TRIANGLE then CPU_TRIANGLE."""
import time
import numpy as np
import viennaps as ps
import viennaps.d3 as v3

ps.Logger.setLogLevel(ps.LogLevel.ERROR)
ps.Length.setUnit("micrometer")
ps.Time.setUnit("min")
Mat = ps.Material


def build():
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=0.3, xExtent=10.0, yExtent=4.0, trenchWidth=4.0,
                  trenchDepth=2.0, taperingAngle=0.0, baseHeight=0.0,
                  periodicBoundary=False, makeMask=True, material=Mat.Si).apply()
    return d


def run(engine):
    d = build()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process()
    p.setDomain(d)
    p.setProcessModel(m)
    p.setProcessDuration(2.0)
    p.setFluxEngineType(getattr(ps.FluxEngineType, engine))
    t0 = time.time()
    p.apply()
    dt = time.time() - t0
    n = np.array(d.getSurfaceMesh().getNodes())
    return dt, float(-n[:, 2].min()), len(n)


print("gpuAvailable:", ps.gpuAvailable())
for eng in ["GPU_TRIANGLE", "CPU_TRIANGLE"]:
    try:
        dt, depth, nn = run(eng)
        print(f"  {eng:13s}: {dt:6.1f}s   max-depth~{depth:.2f} um   nodes {nn}")
    except Exception as e:
        print(f"  {eng:13s}: FAILED -> {str(e)[:160]}")
