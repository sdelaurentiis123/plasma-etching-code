"""ViennaPS-GPU (OptiX) speed on holes d=4/6/8, matched-depth benchmark geometry. Warmed (run twice,
take 2nd). Writes /root/vps_gpu.json. Run with LD_LIBRARY_PATH = nvidia cuda libs."""
import time, json
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
DX, DUR, XE = 0.25, 3.0, 14.0


def run(diam):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=diam / 2, holeDepth=2.0,
                makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE)
    t0 = time.time(); p.apply(); wall = time.time() - t0
    return wall, float(-np.array(d.getSurfaceMesh().getNodes())[:, 2].min())


out = {}
for diam in [4.0, 6.0, 8.0]:
    run(diam)                      # warm (OptiX context + kernel cache)
    w, dep = run(diam)             # timed
    out[diam] = {'wall': w, 'depth': dep}
    print(f"  ViennaPS-GPU d={diam}: {w:.2f}s  depth {dep:.2f}um", flush=True)
json.dump(out, open("/root/vps_gpu.json", "w"))
print("wrote /root/vps_gpu.json", flush=True)
