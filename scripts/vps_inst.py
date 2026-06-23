"""ISOLATION TEST (ViennaPS side): instantaneous bottom etch rate on clean vertical trenches at fixed
depths (MakeTrench trenchDepth=D), via a tiny etch step -> bottom velocity. Compare to petch_inst.py on
the identical geometry to separate flux-physics from etch-dynamics. CPU flux engine (no OptiX needed)."""
import json
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
DS = [1.0, 1.85, 2.7, 3.55, 4.3]
DT = 0.05                              # small step: instantaneous bottom velocity = (depth change)/dt


def vrate(D):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W, trenchDepth=D,
                  makeMask=True, material=ps.Material.Si).apply()
    floor0 = float(np.array(d.getSurfaceMesh().getNodes())[:, 2].min())   # trench floor (~0); field at +D
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DT)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE); p.apply()
    floor1 = float(np.array(d.getSurfaceMesh().getNodes())[:, 2].min())   # floor after (moves down)
    vfloor = (floor0 - floor1) / DT                                       # downward bottom velocity
    print(f"    D={D} AR={D/W:.1f} floor {floor0:.3f}->{floor1:.3f} disp={floor0-floor1:.3f}um "
          f"({(floor0-floor1)/DX:.1f} cells)", flush=True)
    return vfloor, D / W


print("ViennaPS INSTANTANEOUS bottom rate on clean vertical trenches (CPU flux)\n", flush=True)
rows = [vrate(D) for D in DS]
v0 = rows[0][0]
print("  AR    nr(=v/v0)", flush=True)
for v, ar in rows:
    print(f"  {ar:4.1f}   {v / v0:.3f}", flush=True)
json.dump([(ar, v / v0) for v, ar in rows], open("/root/vps_inst.json", "w"))
print("\n  ViennaPS ETCHED nr @AR8.6 was 0.731 (cross-check vs instantaneous).", flush=True)
