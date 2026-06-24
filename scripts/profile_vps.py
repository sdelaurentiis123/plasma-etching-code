#!/usr/bin/env python3
"""Profile upstream ViennaPS-GPU: WHERE does an SF6O2 etch spend time? Decides if a fork has
capturable headroom (coverage re-solves + level-set on CPU = yes; all OptiX ray tracing = no).
Uses ViennaPS TIMING log + a coarse step-scaling probe. cuda + libnvoptix in LD_LIBRARY_PATH."""
import time, numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def build(dur):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=0.04, xExtent=1.2, yExtent=0.24, trenchWidth=0.5,
                  trenchDepth=0.08, taperingAngle=0.0, baseHeight=0.0,
                  periodicBoundary=True, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE)
    return p


# 1) TIMING-level log: ViennaPS prints per-phase (ray trace / advection) timings if available
print("===== ViennaPS TIMING log (1.0 min trench etch, GPU) =====", flush=True)
ps.Logger.setLogLevel(ps.LogLevel.TIMING)
t0 = time.time(); build(1.0).apply(); print(f"[total apply] {time.time()-t0:.2f}s", flush=True)

# 2) step-scaling probe: total time vs process duration -> per-step cost is ~linear; the
#    intercept ~ fixed setup (OptiX/mesh build), slope ~ per-advection-step (flux+level-set).
ps.Logger.setLogLevel(ps.LogLevel.ERROR)
print("\n===== step-scaling probe (total wall vs duration) =====", flush=True)
for dur in [0.25, 0.5, 1.0, 2.0]:
    t0 = time.time(); p = build(dur); p.apply(); wall = time.time() - t0
    d = p.getDomain(); n = np.array(d.getSurfaceMesh().getNodes())
    print(f"  dur {dur:.2f}min  wall {wall:.2f}s  depth {-n[:,2].min():.2f}um", flush=True)
print("done", flush=True)
