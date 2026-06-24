"""Holes match ViennaPS fudge-free; only the TRENCH is too steep -> suspect the periodic-y handling
(_wrap_y), not the core transport. petch trench over-shadows (acts confined in y). Test fudge-free trench:
(a) thin Ly=0.3 + periodic-y [current], (b) wide Ly=2.0 + periodic-y, (c) wide Ly=2.0 NO periodic (truly
infinite-y by extent). If (c)/(b) get GENTLER toward ViennaPS (0.862/0.732) -> the thin periodic-y is the
bug. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.6, 1.4, 2.6, 4.2]
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def arde(dep):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(DURS)) / (np.diff(dep) / np.diff(DURS))[0]


def curve(Ly, periodic, seeds=(0, 1)):
    GEO = dict(Lx=1.5, Ly=Ly, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = periodic   # fudge-free Fflux=1800
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=50000, n_neu=50000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


curve(0.3, 1, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  fudge-free trench: is the periodic-y the bug? (vs ViennaPS 0.862/0.732)\n", flush=True)
for Ly, per, lab in [(0.3, 1, "Ly=0.3 periodic (current)"), (2.0, 1, "Ly=2.0 periodic"),
                     (2.0, 0, "Ly=2.0 NO periodic (wide=infinite)")]:
    ar, nr = curve(Ly, per)
    pp = np.interp(VPS_AR, ar, nr); gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
    print(f"  {lab:36s}: nr@vpsAR {np.round(pp,3)}  gap {gap:.3f}", flush=True)
print(f"  ViennaPS target                       {np.round(VPS_NR,3)}", flush=True)
print("\n  if wide/non-periodic gets gentler toward ViennaPS -> thin periodic-y is the trench bug.", flush=True)
