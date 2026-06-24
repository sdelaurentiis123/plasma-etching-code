"""ViennaPS reference (5.8.1 working): trench + hole ARDE + speed timing -> /root/vps.json. Run with
LD_LIBRARY_PATH set to the nvidia cuda libs. Separate process from petch (avoid CUDA/Warp conflict)."""
import time, json
import numpy as np
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")

TR = dict(dx=0.04, W=0.5, XE=1.5, YE=0.3, durs=[0.4, 0.8, 1.3, 1.9])
HO = dict(dx=0.05, RAD=0.25, XE=1.5, durs=[0.4, 1.0, 1.8, 2.8])


def trench_depth(dur):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=TR['dx'], xExtent=TR['XE'], yExtent=TR['YE'], trenchWidth=TR['W'],
                  trenchDepth=2 * TR['dx'], makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    return float(-np.array(d.getSurfaceMesh().getNodes())[:, 2].min())


def hole_depth(dur):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=HO['dx'], xExtent=HO['XE'], yExtent=HO['XE'], holeRadius=HO['RAD'],
                holeDepth=2 * HO['dx'], makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    return float(-np.array(d.getSurfaceMesh().getNodes())[:, 2].min())


out = {'trench': {}, 'hole': {}, 'speed': {}}
for dr in TR['durs']:
    t0 = time.time(); dep = trench_depth(dr); out['trench'][dr] = dep
    print(f"  VPS trench dur{dr}: {dep:.3f} (AR {dep/TR['W']:.1f}) [{time.time()-t0:.0f}s]", flush=True)
out['speed']['trench_dur1.9_s'] = time.time() - t0   # last trench wall time (the deep one)
for dr in HO['durs']:
    t0 = time.time(); dep = hole_depth(dr); out['hole'][dr] = dep
    print(f"  VPS hole dur{dr}: {dep:.3f} (AR {dep/(2*HO['RAD']):.1f}) [{time.time()-t0:.0f}s]", flush=True)
out['speed']['hole_dur2.8_s'] = time.time() - t0

json.dump(out, open("/root/vps.json", "w"))
print("wrote /root/vps.json", flush=True)
