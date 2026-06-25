"""Run PURE ViennaPS on the de Boer trench geometry (W=2 µm, the experimental regime) and sweep duration
to build its ARDE curve, for direct comparison to the measured de Boer wafer points. Saves
/root/vps_deboer.npz. cuda + libnvoptix in LD_LIBRARY_PATH."""
import time, numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
W, DX, MASK = 2.0, 0.2, 2.0          # de Boer-scale trench; dx=0.2 + capped AR so ViennaPS (CPU level-set) finishes

def vps_depth(dur):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=12, yExtent=1.0, trenchWidth=W, trenchDepth=MASK,
                  taperingAngle=0.0, baseHeight=0.0, periodicBoundary=True, makeMask=True,
                  material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(float(dur))
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    n = np.array(d.getSurfaceMesh().getNodes()); return float(-n[:, 2].min())

durs = [1, 2, 4, 6, 9, 13, 18]       # reach ~AR10-12; capped so the deep CPU level-set stays feasible
deps = []
for du in durs:
    t0 = time.time(); z = vps_depth(du); deps.append(z)
    print(f"ViennaPS deBoer-trench dur {du:>2} depth {z:5.1f} AR {z/W:4.1f} ({time.time()-t0:.0f}s)", flush=True)
np.savez("/root/vps_deboer.npz", dur=np.array(durs, float), dep=np.array(deps), W=W)
print("saved /root/vps_deboer.npz  DONE", flush=True)
