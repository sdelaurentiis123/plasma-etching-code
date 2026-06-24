"""ViennaPS-GPU trench ARDE at dx=0.025 (finer) -> /root/vps_dx025.json. For the dx-convergence proof.
Run with LD_LIBRARY_PATH = nvidia cuda libs."""
import json
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
DX, W, XE, YE = 0.025, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9]


def depth(dur):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W, trenchDepth=2 * DX,
                  makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    return float(-np.array(d.getSurfaceMesh().getNodes())[:, 2].min())


out = {}
for dr in DURS:
    out[dr] = depth(dr)
    print(f"  VPS-GPU dx0.025 trench dur{dr}: {out[dr]:.3f} (AR {out[dr]/W:.1f})", flush=True)
json.dump(out, open("/root/vps_dx025.json", "w"))
print("wrote /root/vps_dx025.json", flush=True)
