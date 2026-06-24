"""Dump ViennaPS trench cross-section (half-width vs depth) at a deep etch -> /root/vps_prof.json.
For the petch-vs-ViennaPS PROFILE comparison (is petch's trench more tapered -> starved floor?)."""
import json
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
DX, W, XE, YE, DUR = 0.04, 0.5, 1.5, 0.3, 1.9

d = v3.Domain()
v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W, trenchDepth=2 * DX,
              makeMask=True, material=ps.Material.Si).apply()
m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
n = np.array(d.getSurfaceMesh().getNodes())
cx = XE / 2.0
# mid-y slice of the surface; for each depth bin, the trench half-width = max |x-cx| of interface points
mid = n[np.abs(n[:, 1] - YE / 2.0) < DX]
zmin = n[:, 2].min()
prof = {}
for z in np.arange(0.0, -zmin, 0.2):           # depth below the mask plane (z=0 ~ substrate top)
    band = mid[np.abs(mid[:, 2] - (-z)) < 0.12]
    if len(band):
        inside = band[np.abs(band[:, 0] - cx) < W]   # within the trench opening +- margin
        if len(inside):
            prof[round(z, 2)] = float(np.abs(inside[:, 0] - cx).max())   # half-width
print("ViennaPS trench depth", round(-zmin, 2), "um", flush=True)
for z in sorted(prof): print(f"  z={z:.1f}: half-width {prof[z]:.3f}", flush=True)
json.dump({'depth': float(-zmin), 'halfwidth': prof, 'W': W}, open("/root/vps_prof.json", "w"))
print("wrote /root/vps_prof.json", flush=True)
