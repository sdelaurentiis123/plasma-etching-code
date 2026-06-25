#!/usr/bin/env python3
"""Capture data for a petch-vs-ViennaPS real-time etch race: same hole, both GPU, record each engine's
centre cross-section over the etch + the TOTAL wall-clock of a single full run. ViennaPS first (OptiX<->
Warp conflict). Saves /root/race.npz for local rendering. cuda + libnvoptix in LD_LIBRARY_PATH."""
import time, pickle, numpy as np
W, DX, MASK, DUR = 0.5, 0.04, 0.08, 2.2          # hole diameter / grid / mask / etch duration (min)
EXT, SUB = 1.2, 0.9
NFR = 14

# ---------------- ViennaPS (before importing warp) ----------------
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

# WARM the engine (first GPU_TRIANGLE call compiles OptiX PTX) -> then time a warm full etch (fair)
_ = vps_run(0.3)
t0 = time.time(); _ = vps_run(DUR); vps_wall = time.time() - t0
print(f"[ViennaPS] WARM full etch wall {vps_wall:.2f}s", flush=True)
# shape snapshots at increasing fractions (centre y-slice profile)
vps_frames = []
for fr in np.linspace(0.07, 1.0, NFR):
    n = vps_run(DUR * fr); sl = n[np.abs(n[:, 1]) < DX]
    vps_frames.append((float(fr), sl[:, 0].copy(), sl[:, 2].copy()))
    print(f"  vps frac {fr:.2f}  depth {-(sl[:,2].min()):.2f}", flush=True)

# ---------------- petch (same hole, faithful config) ----------------
import petch
from petch import threed as t3
GEO = dict(Lx=EXT, Ly=EXT, Lz=2*DX+SUB+0.4, dx=DX, trench_width=W, mask_th=2*DX, sub_top=SUB, hole=True)
p_ = dict(petch.PAR); p_['rate_scale'] = 0.1
fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                 warm_start_coverage=True, sampling="sobol", ion_reflection=True)
# WARM petch (first launch JIT-compiles the Warp kernels ~9s) on a tiny run -> then time the real etch
_ = t3.run_etch_3d(t_end=0.2, n_steps=2, par=p_, flags=fl, n_ion=4000, n_neu=4000,
                   reinit_method="fsm", verbose=False, **GEO)
t0 = time.time()
g = t3.run_etch_3d(t_end=DUR, n_steps=NFR, par=p_, flags=fl, n_ion=40000, n_neu=40000,
                   reinit_method="fsm", verbose=False, record_depth_every=1, record_frames=True, **GEO)
petch_wall = time.time() - t0
print(f"[petch] WARM full etch wall {petch_wall:.2f}s  -> {vps_wall/petch_wall:.1f}x faster", flush=True)

pframes = [dict(t=f['t'], depth=f['depth'], phi=f['phi_xz'].astype(np.float32)) for f in g['frames']]
with open("/root/race.pkl", "wb") as fh:
    pickle.dump(dict(W=W, DX=DX, SUB=SUB, EXT=EXT, DUR=DUR,
                     vps_wall=vps_wall, petch_wall=petch_wall, vps_frames=vps_frames,
                     petch_frames=pframes, xs=g['xs'], zs=g['zs'], sub_top=g['sub_top']), fh)
print("saved /root/race.pkl", flush=True)
print("DONE", flush=True)
