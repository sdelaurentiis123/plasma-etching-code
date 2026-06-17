#!/usr/bin/env python3
"""Decisive test: does ViennaPS's TRENCH ARDE match de Boer (steep, ~0.43 at AR10) or petch (gentle,
~0.9)? If ViennaPS is gentle too, the de Boer gap is a SHARED ballistic limit (petch is fine vs the
tool). If ViennaPS is steep, petch has a real bug. Both run at increasing depth; compare the NORMALIZED
bottom-rate vs AR (shape, no rate-matching needed). ViennaPS first (OptiX<->Warp). PETCH_DEVICE=cuda."""
import time
import numpy as np

DX, W, XE, YE = 0.03, 0.5, 1.5, 0.3       # sub-micron trench, thin-y (quasi-2D)
DURS = [0.4, 0.8, 1.3, 1.9, 2.6]          # min; increasing depth

# ---------- ViennaPS (CPU; profile engine-independent) ----------
import viennaps as ps
import viennaps.d3 as v3
ps.Logger.setLogLevel(ps.LogLevel.ERROR); ps.Length.setUnit("micrometer"); ps.Time.setUnit("min")


def vps_depth(dur):
    d = v3.Domain()
    v3.MakeTrench(domain=d, gridDelta=DX, xExtent=XE, yExtent=YE, trenchWidth=W,
                  trenchDepth=2 * DX, makeMask=True, material=ps.Material.Si).apply()
    m = v3.SF6O2Etching(v3.SF6O2Etching.defaultParameters())
    p = v3.Process(); p.setDomain(d); p.setProcessModel(m); p.setProcessDuration(dur)
    p.setFluxEngineType(ps.FluxEngineType.CPU_TRIANGLE)
    p.apply()
    n = np.array(d.getSurfaceMesh().getNodes())
    return float(-n[:, 2].min())


vps = []
for dr in DURS:
    dep = vps_depth(dr); vps.append((dr, dep))
    print(f"  ViennaPS dur {dr}: depth {dep:.3f} (AR {dep/W:.1f})", flush=True)


def arde(pts):
    dr = np.array([x[0] for x in pts]); dep = np.array([x[1] for x in pts])
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dr)
    nr = rate / rate[0]
    return armid, nr


vps_ar, vps_nr = arde(vps)

# ---------- petch (periodic-y; same trench) ----------
import petch
from petch import threed as t3
GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + 5.0 + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=5.0 + 0.3, hole=False)


def petch_depth(dur):
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol")
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 20)), par=p, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


petch_depth(0.4)  # warm
pe = []
for dr in DURS:
    dep = petch_depth(dr); pe.append((dr, dep))
    print(f"  petch dur {dr}: depth {dep:.3f} (AR {dep/W:.1f})", flush=True)
pe_ar, pe_nr = arde(pe)

print("\n  normalized bottom rate vs AR:", flush=True)
print(f"    {'AR':>5} {'ViennaPS':>9} {'petch':>7}", flush=True)
allar = np.union1d(np.round(vps_ar, 1), np.round(pe_ar, 1))
for a in allar:
    v = np.interp(a, vps_ar, vps_nr); pp = np.interp(a, pe_ar, pe_nr)
    print(f"    {a:5.1f} {v:9.2f} {pp:7.2f}", flush=True)
print("\n  if ViennaPS ~ petch (both gentle) -> de Boer gap is a SHARED ballistic limit (petch matches the tool).", flush=True)
print("  if ViennaPS is much steeper -> petch has a real discrepancy vs ViennaPS to fix.", flush=True)
