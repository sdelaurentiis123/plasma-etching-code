"""Why is fudge-free petch ~0.08-0.10 too steep at the deepest AR (consistent on both geometries)? Test
if it's DISCRETIZATION: run ViennaPS + petch TRENCH ARDE at dx=0.04 AND dx=0.025; if the petch-vs-ViennaPS
gap SHRINKS as dx->0, it's resolution (the narrow deep floor is under-resolved). If it persists, structural.
Fudge-free (no cal_F), wrap-fix in place. ViennaPS FIRST (OptiX<->Warp)."""
import time
import numpy as np
W, XE, YE, SUB = 0.5, 1.5, 0.3, 7.0
DURS = [0.4, 0.8, 1.3, 1.9]

import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def vps_trench(dur, dx):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=dx, xExtent=XE, yExtent=YE, trenchWidth=W, trenchDepth=2 * dx,
                  makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE); p.apply()
    n = np.array(d.getSurfaceMesh().getNodes()); return float(-n[:, 2].min())


def arde(dep):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(DURS)) / (np.diff(dep) / np.diff(DURS))[0]


vps_ref = {}
for dx in [0.04, 0.025]:
    dep = []
    for dr in DURS:
        t0 = time.time(); dd = vps_trench(dr, dx); dep.append(dd)
        print(f"  ViennaPS dx={dx} dur{dr}: {dd:.3f} (AR {dd/W:.1f}) [{time.time()-t0:.0f}s]", flush=True)
    vps_ref[dx] = arde(dep)
    print(f"  ViennaPS dx={dx} ARDE: AR {np.round(vps_ref[dx][0],2)} nr {np.round(vps_ref[dx][1],3)}\n", flush=True)

import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import petch
from petch import threed as t3


def petch_trench(dx, seeds=(0, 1, 2)):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * dx + SUB + 0.3, dx=dx, trench_width=W, mask_th=2 * dx, sub_top=SUB + 0.3, hole=False)
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1   # fudge-free
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


petch_trench(0.04, seeds=(0,))  # warm
print("  petch fudge-free vs ViennaPS, per dx:", flush=True)
for dx in [0.04, 0.025]:
    v_ar, v_nr = vps_ref[dx]
    ar, nr = petch_trench(dx)
    pp = np.interp(v_ar, ar, nr); gap = float(np.sqrt(np.mean((pp - v_nr) ** 2)))
    print(f"    dx={dx}: petch nr@vpsAR {np.round(pp,3)} vs {np.round(v_nr,3)}  gap {gap:.3f}", flush=True)
print("\n  gap shrinks 0.04->0.025 -> residual is DISCRETIZATION (under-resolved floor). persists -> structural.", flush=True)
