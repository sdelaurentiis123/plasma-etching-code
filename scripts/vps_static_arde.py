"""FAIR ViennaPS reference: static-equivalent nr(AR) on a W=0.5um trench.
Pre-carve the trench at depth D (MakeTrench trenchDepth), etch a SHORT duration with the full
SF6O2 model (reflecting ions included -- ViennaPS default), and take the instantaneous floor
rate = (depth_after - D)/dur. Normalize by the flat-field rate (MakePlane). GPU_TRIANGLE."""
import time
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
W, DX, XE, YE, DUR = 0.5, 0.05, 2.0, 1.0, 0.15

def etch(domain, dur):
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(domain); p.setProcessModel(m); p.setProcessDuration(float(dur))
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()

def field_rate():
    d = v3.Domain()
    v3.MakePlane(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, height=0.0,
                 periodicBoundary=False, material=ps.Material.Si).apply()
    etch(d, DUR)
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min()) / DUR

def floor_rate(D):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W,
                  trenchDepth=D, taperingAngle=0.0, baseHeight=0.0, periodicBoundary=True,
                  makeMask=True, material=ps.Material.Si).apply()
    etch(d, DUR)
    n = np.array(d.getSurfaceMesh().getNodes())
    return (float(-n[:, 2].min()) - D) / DUR

r0 = field_rate()
print(f"field rate {r0:.3f} um/min", flush=True)
etch_warm = floor_rate(0.5)                       # warm OptiX once more on trench geometry
for AR in [2, 4, 6, 8, 10]:
    t0 = time.time(); r = floor_rate(AR * W)
    print(f"VPS_STATIC AR{AR:3d}: nr={r/r0:.3f}   ({time.time()-t0:.0f}s)", flush=True)
print("VPS DONE", flush=True)
