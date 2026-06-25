"""Race v2: petch (SMOOTH, high per-facet ray sampling) vs ViennaPS, same hole, depth-matched, both warm.
Smooth petch (rays_per_point) also brings the speed ratio into the documented ~14x regime. ViennaPS first.
Saves /root/race2.pkl. cuda + libnvoptix in LD_LIBRARY_PATH."""
import time, pickle, numpy as np
W, DX, MASK, DUR, EXT, SUB = 0.5, 0.06, 0.08, 1.6, 1.2, 4.0
NFR = 9

import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")

def vps_run(dur):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=EXT, yExtent=EXT, holeRadius=W/2,
                holeDepth=MASK, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(float(dur))
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    return np.array(d.getSurfaceMesh().getNodes())

_ = vps_run(0.3)                                            # warm (OptiX compile)
t0 = time.time(); _ = vps_run(DUR); vps_wall = time.time() - t0
vps_depth = -_[:, 2].min()
print(f"[ViennaPS] warm wall {vps_wall:.2f}s depth {vps_depth:.2f}", flush=True)
vps_frames = []
for fr in np.linspace(0.07, 1.0, NFR):
    n = vps_run(DUR * fr); sl = n[np.abs(n[:, 1]) < DX]
    vps_frames.append((float(fr), sl[:, 0].copy(), sl[:, 2].copy()))
print("vps snapshots done", flush=True)

import petch
from petch import threed as t3
GEO = dict(Lx=EXT, Ly=EXT, Lz=2*DX+SUB+0.5, dx=DX, trench_width=W, mask_th=0.3, sub_top=SUB, hole=True)
fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                 warm_start_coverage=True, sampling="sobol", ion_reflection=True)
def prun(rate, ns, cap=False):
    p_ = dict(petch.PAR); p_['rate_scale'] = rate
    kw = dict(record_depth_every=1, record_frames=True) if cap else {}
    t0 = time.time()
    g = t3.run_etch_3d(t_end=DUR, n_steps=ns, par=p_, flags=fl, n_ion=40000, n_neu=40000,
                       surf_smooth=0.7, reinit_method="fsm", verbose=False, **GEO, **kw)   # SMOOTH (regularized)
    return time.time()-t0, t3.max_depth_3d(g), g
prun(0.1, 2)                                               # warm
best = min([0.1, 0.13, 0.16, 0.2], key=lambda r: abs(prun(r, 40)[1] - vps_depth))
w, dep, g = prun(best, 40, cap=True)
print(f"[petch] SMOOTH warm wall {w:.2f}s depth {dep:.2f} (vps {vps_depth:.2f}) rate {best} -> {vps_wall/w:.0f}x", flush=True)
pickle.dump(dict(W=W, DX=DX, EXT=EXT, DUR=DUR, vps_wall=vps_wall, vps_frames=vps_frames,
                 petch_wall=w, petch_sub=SUB, xs=g['xs'], zs=g['zs'],
                 petch_frames=[dict(t=f['t'], depth=f['depth'], phi=f['phi_xz'].astype(np.float32)) for f in g['frames']]),
            open("/root/race2.pkl", "wb"))
print("saved /root/race2.pkl   DONE", flush=True)
