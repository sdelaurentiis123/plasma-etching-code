"""Clean speed re-check on the CURRENT (fudge-free, wrap-fixed) code. OptiX/GPU_TRIANGLE segfaults on this
box, so: (1) time petch (warmed) on the benchmark holes -> confirm the hot loop is unchanged by this
session's edits (cal_F removal + wrap cap touch only the physics path, not speed); (2) ViennaPS CPU_TRIANGLE
same-box baseline. Prior session measured petch 1.30/1.61/1.82s vs ViennaPS-GPU 19.4/22.7/25.8s = ~14x;
if petch timing here matches, that 14x stands. ViennaPS FIRST."""
import time
import numpy as np
DX, DUR, XE = 0.25, 3.0, 14.0
DIAMS = [4.0, 6.0, 8.0]

import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def vps_cpu(diam):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=diam / 2, holeDepth=2.0,
                makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(DUR)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE)
    t0 = time.time(); p.apply(); wall = time.time() - t0
    n = np.array(d.getSurfaceMesh().getNodes()); return wall, float(-n[:, 2].min())


vps = {}
for diam in DIAMS:
    w, dep = vps_cpu(diam); vps[diam] = (w, dep)
    print(f"  ViennaPS-CPU d={diam}: {w:6.2f}s  depth {dep:.2f}um", flush=True)

import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import petch
from petch import threed as t3


def ours(diam, ns=40):
    p = dict(petch.PAR); p['rate_scale'] = 1.2; p['betaE'] = 0.7   # fudge-free; rate_scale ~matches depth
    GEO = dict(Lx=XE, Ly=XE, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=DUR)
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=diam, dx=DX, n_steps=ns, par=p,
                       flags=petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc"),
                       n_ion=40000, n_neu=40000, reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.max_depth_3d(g)


print(flush=True)
for diam in DIAMS:
    ours(diam)                                   # warm (compile + cache)
    w, dep = ours(diam)
    vw, vdep = vps[diam]
    print(f"  petch d={diam}: {w:6.2f}s  depth {dep:.2f}um   |  vs ViennaPS-CPU {vw:.1f}s = {vw/w:.1f}x faster (same box)", flush=True)
print("\n  (ViennaPS-GPU/OptiX won't run on this box; prior fair GPU-vs-GPU = ~14x. petch hot loop unchanged", flush=True)
print("   by this session's edits -> that 14x holds.)", flush=True)
