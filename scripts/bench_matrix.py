#!/usr/bin/env python3
"""Rigorous petch-vs-ViennaPS validation matrix across feature TYPES (trench + hole) and
aspect ratios. ViennaPS GPU runs FIRST (OptiX<->Warp single-process conflict), then petch.
Metric: normalized bottom-rate nr vs AR (rate_scale-independent), per feature type. cuda."""
import os, time, numpy as np
os.environ.setdefault("PETCH_DEVICE", "cuda")
W, DX, MASK = 0.5, 0.04, 0.08

# ---------------- ViennaPS (before importing petch/warp) ----------------
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")
ENG = "GPU_TRIANGLE"

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
    p.setFluxEngineType(getattr(ps.FluxEngineType, ENG))
    p.apply()
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min())

t0 = time.time()
try:
    d0 = vps_depth(False, 0.2)
    print(f"[ViennaPS GPU OK] trench@0.2min depth {d0:.3f}um ({time.time()-t0:.0f}s incl OptiX compile)", flush=True)
except Exception as e:
    ENG = "CPU_TRIANGLE"; print(f"[ViennaPS GPU failed ({e}); using {ENG}]", flush=True)

VDUR = {False: [0.2, 0.4, 0.6, 0.8, 1.0, 1.2], True: [0.3, 0.6, 0.9, 1.2, 1.6, 2.0]}
vps = {}
for hole in [False, True]:
    deps = []
    for du in VDUR[hole]:
        z = vps_depth(hole, du); deps.append(z)
        print(f"  ViennaPS {'hole ' if hole else 'trench'} dur {du:.1f}  depth {z:.3f}  AR {z/W:.1f}", flush=True)
    vps[hole] = (np.array(VDUR[hole]), np.array(deps))

# ---------------- petch (after) ----------------
import petch
from petch import threed as t3

def petch_depth(hole, dur, seed=0):
    sub = 8.0
    GEO = dict(Lx=1.2, Ly=(1.2 if hole else 0.24), Lz=2*DX+sub+0.3, dx=DX, trench_width=W,
               mask_th=2*DX, sub_top=sub+0.3, hole=hole)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['n_fp'] = 3
    if not hole:
        p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc", ion_reflection=True)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur*22)), par=p, flags=fl,
                       n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False, seed_offset=seed, **GEO)
    return t3.max_depth_3d(g)

PDUR = {False: [0.7, 1.1, 1.5, 1.9, 2.3, 2.7], True: [1.0, 1.6, 2.2, 2.8, 3.4, 4.0]}
petchd = {}
for hole in [False, True]:
    deps = []
    for du in PDUR[hole]:
        z = petch_depth(hole, du); deps.append(z)
        print(f"  petch    {'hole ' if hole else 'trench'} dur {du:.1f}  depth {z:.3f}  AR {z/W:.1f}", flush=True)
    petchd[hole] = (np.array(PDUR[hole]), np.array(deps))

# ---------------- scorecard ----------------
def nr_curve(durs, deps):
    ar = 0.5 * (deps[1:] + deps[:-1]) / W
    nr = np.diff(deps) / np.diff(durs)
    return ar, nr / nr[0]

print("\n================  SCORECARD: normalized etch rate nr(AR), petch vs ViennaPS  ================", flush=True)
out = {}
for hole in [False, True]:
    av, nv = nr_curve(*vps[hole]); ap, npc = nr_curve(*petchd[hole])
    lo, hi = max(av.min(), ap.min()), min(av.max(), ap.max())
    arc = np.linspace(lo, hi, 4)
    vi, pi = np.interp(arc, av, nv), np.interp(arc, ap, npc)
    rmse = float(np.sqrt(np.mean((vi - pi) ** 2)))
    tag = "HOLE  " if hole else "TRENCH"
    print(f"  {tag}  AR {np.round(arc,1)}", flush=True)
    print(f"          ViennaPS nr {np.round(vi,3)}", flush=True)
    print(f"          petch    nr {np.round(pi,3)}    RMSE {rmse:.3f}    {'PASS' if rmse<0.15 else 'CHECK'}", flush=True)
    out[tag] = rmse
print(f"\n  overlap-AR nr RMSE: trench {out['TRENCH']:.3f}  hole {out['HOLE  ']:.3f}", flush=True)
print("done", flush=True)
