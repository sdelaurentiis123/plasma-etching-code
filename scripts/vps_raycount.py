#!/usr/bin/env python3
"""Is ViennaPS's trench ARDE (0.73) converged, or is ViennaPS itself under-sampled? Sweep ViennaPS
raysPerPoint (default 1000) on the SAME trench and see if its bottom-rate ARDE drifts GENTLER (toward
petch's converged ~0.79). If yes -> petch and ViennaPS AGREE at convergence; the apparent gap was a
ray-count mismatch. If ViennaPS stays ~0.73 -> petch is genuinely gentler (real model diff). CPU."""
import time
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9]


def vps_depth(dur, rpp):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W,
                  trenchDepth=2 * DX, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE)
    rt = ps.RayTracingParameters(); rt.raysPerPoint = rpp
    p.setParameters(rt)
    p.apply()
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min())


def nr_of(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    r = np.diff(dep) / np.diff(DURS); return armid, r / r[0]


# default 1000-ray reference (already measured): 1.297/2.437/3.666/4.918
ref = np.array([1.297, 2.437, 3.666, 4.918])
ar0, nr0 = nr_of(ref)
print(f"  ViennaPS raysPerPoint=1000 (cached): AR {np.round(ar0,2)} nr {np.round(nr0,3)}", flush=True)
for rpp in [4000]:
    dep = []
    for dr in DURS:
        t0 = time.time(); dd = vps_depth(dr, rpp); dep.append(dd)
        print(f"    rpp{rpp} dur{dr}: {dd:.3f} (AR {dd/W:.1f}) [{time.time()-t0:.0f}s]", flush=True)
    ar, nr = nr_of(dep)
    print(f"  ViennaPS raysPerPoint={rpp}: AR {np.round(ar,2)} nr {np.round(nr,3)}", flush=True)
print("\n  petch converged (400k rays): nr ~0.92/0.79. If ViennaPS drifts toward that -> they AGREE.", flush=True)
