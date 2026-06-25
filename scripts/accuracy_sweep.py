#!/usr/bin/env python3
"""Full accuracy sweep: petch vs ViennaPS, TRENCH and HOLE, normalized etch rate nr vs aspect ratio.
ViennaPS-GPU first (OptiX<->Warp single-process conflict), then petch (current faithful config). Each
engine/feature etched to increasing durations -> depth(t) -> nr=rate/rate_lowAR, ar=depth/W; interpolated
to a common AR grid. Saves /root/accuracy_sweep.npz for plotting. PETCH_DEVICE=cuda + libnvoptix in path."""
import os, time, numpy as np
os.environ.setdefault("PETCH_DEVICE", "cuda")
W, DX, MASK = 0.5, 0.04, 0.08
R = {}

# ---------------- ViennaPS (before importing warp) ----------------
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")

def vps_depth(hole, dur):
    d = v3.Domain()
    if hole:
        v3.MakeHole(domain=d, gridDelta=DX, xExtent=1.2, yExtent=1.2, holeRadius=W/2,
                    holeDepth=MASK, makeMask=True, material=ps.Material.Si).apply()
    else:
        v3.MakeTrench(domain=d, gridDelta=DX, xExtent=1.2, yExtent=0.24, trenchWidth=W,
                      trenchDepth=MASK, taperingAngle=0.0, baseHeight=0.0,
                      periodicBoundary=True, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(float(dur))
    p.setFluxEngineType(ps.FluxEngineType.GPU_TRIANGLE); p.apply()
    n = np.array(d.getSurfaceMesh().getNodes()); return float(-n[:, 2].min())

VDUR = {False: [0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2],          # trench
        True:  [0.2, 0.4, 0.6, 0.9, 1.2, 1.6, 2.0]}            # hole (steeper -> longer)
for hole in [False, True]:
    deps = []
    for du in VDUR[hole]:
        z = vps_depth(hole, du); deps.append(z)
        print(f"ViennaPS {'hole ' if hole else 'trench'} dur {du:.2f} depth {z:.3f} AR {z/W:.1f}", flush=True)
    R[f"vps_{'hole' if hole else 'trench'}_dur"] = np.array(VDUR[hole])
    R[f"vps_{'hole' if hole else 'trench'}_dep"] = np.array(deps)
    np.savez("/root/accuracy_sweep.npz", **R)

# ---------------- petch (current faithful config) ----------------
import petch
from petch import threed as t3

def petch_depth(hole, dur, seed=0):
    sub = 6.0
    GEO = dict(Lx=1.2, Ly=(1.2 if hole else 0.24), Lz=2*DX+sub+0.3, dx=DX, trench_width=W,
               mask_th=2*DX, sub_top=sub+0.3, hole=hole)
    p = dict(petch.PAR); p['rate_scale'] = 0.1
    if not hole:
        p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", ion_reflection=True)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur*22)), par=p, flags=fl, n_ion=60000, n_neu=60000,
                       reinit_method="fsm", verbose=False, seed_offset=seed, **GEO)
    return t3.max_depth_3d(g)

PDUR = {False: [0.5, 0.8, 1.1, 1.5, 1.9, 2.3, 2.7],
        True:  [0.7, 1.2, 1.8, 2.5, 3.3, 4.2, 5.2]}
for hole in [False, True]:
    deps = np.mean([[petch_depth(hole, du, s) for du in PDUR[hole]] for s in (0, 100)], axis=0)
    for du, z in zip(PDUR[hole], deps):
        print(f"petch    {'hole ' if hole else 'trench'} dur {du:.2f} depth {z:.3f} AR {z/W:.1f}", flush=True)
    R[f"petch_{'hole' if hole else 'trench'}_dur"] = np.array(PDUR[hole])
    R[f"petch_{'hole' if hole else 'trench'}_dep"] = np.array(deps)
    np.savez("/root/accuracy_sweep.npz", **R)
print("ALL DONE", flush=True)
