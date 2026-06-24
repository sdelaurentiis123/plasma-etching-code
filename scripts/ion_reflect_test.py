#!/usr/bin/env python3
"""The deep-floor bias hypothesis: ViennaPS ions REFLECT off sidewalls (ReflectionConedCosine, grazing)
-> reach the deep floor -> gentler ARDE. petch ions are single-bounce (stick on sidewall) -> under-deliver
ions to the floor -> too steep. Test fudge-free petch with ion_reflection ON vs OFF; if ON climbs toward
ViennaPS (0.862/0.732) the missing ion reflection IS the bias. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, SUB = 0.04, 0.5, 7.0
DURS = [0.6, 1.4, 2.6, 4.2]
GEO = dict(Lx=1.5, Ly=0.3, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W, mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
VPS_AR = np.array([3.73, 6.10, 8.58]); VPS_NR = np.array([1.0, 0.862, 0.732])


def arde(dep):
    dep = np.asarray(dep, float)
    armid = 0.5 * (dep[1:] + dep[:-1]) / W
    r = np.diff(dep) / np.diff(DURS); return armid, r / r[0]


def curve(ion_refl, seeds=(0, 1)):
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1   # fudge-free Fflux=1800
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc",
                         ion_reflection=ion_refl)
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False,
                         seed_offset=sd * 100, **GEO)) for dr in DURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds))


curve(False, seeds=(0,))  # warm
print(f"device={t3.DEVICE}  fudge-free trench: ion reflection ON vs OFF (vs ViennaPS 0.862/0.732)\n", flush=True)
for refl, lab in [(False, "ion_reflection OFF (current default)"), (True, "ion_reflection ON  (like ViennaPS)")]:
    ar, nr = curve(refl)
    pp = np.interp(VPS_AR, ar, nr)
    gap = float(np.sqrt(np.mean((pp - VPS_NR) ** 2)))
    print(f"  {lab}: nr@vpsAR {np.round(pp,3)}  gapRMSE {gap:.3f}", flush=True)
print(f"  ViennaPS target                     {np.round(VPS_NR,3)}", flush=True)
print("\n  if ON closes the gap -> missing ion sidewall reflection was the bias.", flush=True)
