#!/usr/bin/env python3
"""Find the REAL petch-vs-ViennaPS trench-ARDE gap across the full AR range (not 3 shallow points) and
FIX petch's over-gentle bias. Hypothesis (from FINDINGS): petch's flux smoothing is stronger than
ViennaPS's -> over-feeds the HARC floor -> gentler ARDE. Sweep flux_smooth_alpha to steepen petch onto
ViennaPS.

ViennaPS FIRST (all of it; OptiX<->Warp conflict), writes /root/vps_curve.txt as it goes, then petch.
PETCH_DEVICE=cuda."""
import time
import numpy as np

DX, W, XE, YE = 0.04, 0.5, 1.5, 0.3
DURS = [0.4, 0.8, 1.3, 1.9, 2.6, 3.4, 4.3]      # deeper -> reach AR~12-14
SUB = 7.0
LZ = 2 * DX + SUB + 0.3

# ---------- ViennaPS (CPU_TRIANGLE) ----------
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
with open("/root/vps_curve.txt", "w") as f:
    for dr in DURS:
        t0 = time.time(); dep = vps_depth(dr); vps.append((dr, dep))
        line = f"  ViennaPS dur {dr}: depth {dep:.3f} (AR {dep/W:.1f})  [{time.time()-t0:.0f}s]"
        print(line, flush=True); f.write(line + "\n"); f.flush()


def arde(pts):
    dr = np.array([x[0] for x in pts]); dep = np.array([x[1] for x in pts])
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dr)
    return armid, rate / rate[0]


vps_ar, vps_nr = arde(vps)

# ---------- petch: sweep flux_smooth_alpha (steepen onto ViennaPS) ----------
import petch
from petch import threed as t3
GEO = dict(Lx=XE, Ly=YE, Lz=LZ, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)


def petch_curve(alpha, n_smooth, seeds=(0, 1)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
        p['flux_smooth'] = n_smooth; p['flux_smooth_alpha'] = alpha
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = []
        for dr in DURS:
            g = t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl,
                               n_ion=40000, n_neu=40000, reinit_method="fsm", verbose=False,
                               seed_offset=sd * 100, **GEO)
            deps.append(t3.center_depth_3d(g))
        deps = np.array(deps)
        accd = deps if accd is None else accd + deps
    dep = accd / len(seeds)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(DURS)
    return armid, rate / rate[0], dep


# warm
petch_curve(1.0, 1, seeds=(0,))
print(f"\n  ViennaPS curve: AR {np.round(vps_ar,2)}  nr {np.round(vps_nr,3)}\n", flush=True)
print(f"  {'alpha':>5} {'nsm':>3}   {'normalized rate vs ViennaPS AR':>34}   gapRMSE", flush=True)
for alpha, nsm in [(1.0, 1), (0.5, 1), (0.25, 1), (0.0, 0)]:
    ar, nr, dep = petch_curve(alpha, nsm)
    # compare petch nr interpolated to ViennaPS AR points
    pp = np.interp(vps_ar, ar, nr)
    gap = float(np.sqrt(np.mean((pp - vps_nr) ** 2)))
    print(f"  {alpha:5.2f} {nsm:3d}   {np.round(pp,3)!s:>34}   {gap:.3f}", flush=True)
    print(f"             petch own AR {np.round(ar,1)}  maxAR {ar.max():.1f}", flush=True)
print(f"\n  target (ViennaPS)              {np.round(vps_nr,3)}", flush=True)
print("  lower gapRMSE = petch steepened onto ViennaPS. Pick the alpha that minimizes it.", flush=True)
