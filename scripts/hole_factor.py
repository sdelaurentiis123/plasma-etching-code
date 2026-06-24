#!/usr/bin/env python3
"""Decisive: is petch's ~1.8x deep-floor flux deficit a SINGLE geometry-independent constant (clean
calibration) or geometry-dependent (real bug)? Trench needs Fflux x4 to match ViennaPS. Test the HOLE:
run ViennaPS hole ARDE (reference) + petch hole at Fflux x1/x4/x8. If x4 matches the hole too -> ONE
constant fixes both = defensible calibration. ViennaPS FIRST (OptiX<->Warp). dx=0.05 for tractable
ViennaPS 3D holes."""
import time
import numpy as np
DX, RAD, XE = 0.05, 0.25, 1.5     # hole radius 0.25 -> diameter W=0.5; coarser dx for ViennaPS 3D speed
DURS = [0.4, 1.0, 1.8, 2.8]
W = 2 * RAD

import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def vps_hole(dur):
    d = v3.Domain()
    v3.MakeHole(domain=d, gridDelta=DX, xExtent=XE, yExtent=XE, holeRadius=RAD, holeDepth=2 * DX,
                makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE); p.apply()
    n = np.array(d.getSurfaceMesh().getNodes()); return float(-n[:, 2].min())


def arde(dep):
    dep = np.asarray(dep, float)
    return (dep[1:] + dep[:-1]) * 0.5 / W, (np.diff(dep) / np.diff(DURS)) / (np.diff(dep) / np.diff(DURS))[0]


vps = []
for dr in DURS:
    t0 = time.time(); dd = vps_hole(dr); vps.append(dd)
    print(f"  ViennaPS hole dur {dr}: {dd:.3f} (AR {dd/W:.1f}) [{time.time()-t0:.0f}s]", flush=True)
v_ar, v_nr = arde(vps)
print(f"  ViennaPS hole ARDE: AR {np.round(v_ar,2)} nr {np.round(v_nr,3)}\n", flush=True)

import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import petch
from petch import threed as t3
GEO = dict(Lx=XE, Ly=XE, Lz=2 * DX + 7.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=7.3, hole=True)


def petch_hole(fflux):
    accd = None
    for sd in (0, 1):
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['Fflux'] = fflux
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.max_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / 2)


petch_hole(1800.0)  # warm-ish
print(f"  petch hole vs ViennaPS hole nr {np.round(v_nr,3)}:", flush=True)
for k in [1, 4, 8]:
    ar, nr = petch_hole(1800.0 * k)
    pp = np.interp(v_ar, ar, nr); gap = float(np.sqrt(np.mean((pp - v_nr) ** 2)))
    print(f"    Fflux x{k}: nr@vpsAR {np.round(pp,3)}  gapRMSE {gap:.3f}", flush=True)
print(f"\n  TRENCH matched at x4. If the HOLE also matches near x4 -> ONE constant fixes both (clean", flush=True)
print(f"  calibration). If the hole wants a very different x -> geometry-dependent bug.", flush=True)
