"""How reproducible is ViennaPS itself? Run the trench ARDE with different random seeds + ray counts and
measure the spread in the deep-AR normalized rate. If ViennaPS's own deep-AR varies ~0.10, then petch's
~0.10 gap is WITHIN ViennaPS's noise (= as accurate as ViennaPS). LD_LIBRARY_PATH = nvidia cuda libs."""
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9]


def depth(dur, rpp, seed):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W, trenchDepth=2 * DX,
                  makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE)
    rt = ps.RayTracingParameters(); rt.raysPerPoint = rpp; rt.rngSeed = seed; rt.useRandomSeeds = False
    p.setParameters(rt); p.apply()
    return float(-np.array(d.getSurfaceMesh().getNodes())[:, 2].min())


def arde(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    nr = (np.diff(dep) / np.diff(DURS)); nr = nr / nr[0]
    return armid, nr


print("ViennaPS reproducibility: trench ARDE deep-rate vs seed/rays\n", flush=True)
runs = []
for rpp, seed, lab in [(1000, 15, "rpp1000 seedA"), (1000, 92, "rpp1000 seedB"),
                       (1000, 7, "rpp1000 seedC"), (4000, 15, "rpp4000 seedA")]:
    dep = [depth(dr, rpp, seed) for dr in DURS]
    ar, nr = arde(dep)
    deep = float(nr[-1])      # normalized rate at the deepest AR point
    runs.append(deep)
    print(f"  {lab}: AR {np.round(ar,1)}  nr {np.round(nr,3)}  deep={deep:.3f}", flush=True)
print(f"\n  ViennaPS deep-AR nr: mean {np.mean(runs):.3f}  spread (max-min) {max(runs)-min(runs):.3f}  std {np.std(runs):.3f}", flush=True)
print("  if spread ~0.10 -> petch's 0.10 gap is within ViennaPS's own reproducibility.", flush=True)
