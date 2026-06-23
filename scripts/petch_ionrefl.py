"""Mechanism test: does ion reflection close petch's deep-AR ARDE gap vs ViennaPS (0.731 at AR 8.6)?
Run petch trench ARDE with ion_reflection OFF vs ON, interpolate normalized bottom rate onto
ViennaPS AR points [3.7, 6.1, 8.6]. PETCH_DEVICE=cuda."""
import os, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE, SUB = 0.04, 0.5, 1.5, 0.3, 6.0
DURS = [0.7, 1.1, 1.5, 1.9, 2.3, 2.7]      # span AR ~3 to ~9
VPS_AR = np.array([3.7, 6.1, 8.6])
VPS_NR = np.array([1.0, 0.861, 0.731])


def depth(dur, ion_refl):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc",
                     ion_reflection=ion_refl)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 22)), par=p, flags=fl,
                       n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


def arde(deps):
    deps = np.asarray(deps, float)
    armid = 0.5 * (deps[1:] + deps[:-1]) / W
    nr = np.diff(deps) / np.diff(DURS); nr = nr / nr[0]
    return armid, nr


for refl in [False, True]:
    deps = [depth(d, refl) for d in DURS]
    ar, nr = arde(deps)
    nr_at_vps = np.interp(VPS_AR, ar, nr)            # petch nr at ViennaPS AR points
    gap = float(np.sqrt(np.mean((nr_at_vps - VPS_NR) ** 2)))
    print(f"\n  ion_reflection={refl}", flush=True)
    print(f"    depths {np.round(deps,2)}", flush=True)
    print(f"    petch AR   {np.round(ar,1)}  nr {np.round(nr,3)}", flush=True)
    print(f"    @ ViennaPS AR {VPS_AR}: petch {np.round(nr_at_vps,3)} vs VPS {VPS_NR}", flush=True)
    print(f"    deep(AR8.6): petch {nr_at_vps[-1]:.3f} vs VPS 0.731   RMSE {gap:.3f}", flush=True)
print("\ndone", flush=True)
