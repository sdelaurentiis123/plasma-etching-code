"""Does coverage-convergence (n_fp flux<->coverage iterations/step) close petch's residual deep-AR
gap? ViennaPS iterates etchant coverage to convergence each step; petch warm-starts with n_fp=1.
Test n_fp = 1, 3, 6 (warm-start on) on the trench ARDE vs ViennaPS [1, 0.861, 0.731]. cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE, SUB = 0.04, 0.5, 1.5, 0.3, 6.0
DURS = [0.7, 1.1, 1.5, 1.9, 2.3, 2.7]
VPS_AR = np.array([3.7, 6.1, 8.6]); VPS_NR = np.array([1.0, 0.861, 0.731])


def depth(dur, n_fp):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['n_fp'] = n_fp
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc", ion_reflection=True)
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 22)), par=p, flags=fl,
                       n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False, **GEO)
    return t3.center_depth_3d(g)


def arde(deps):
    deps = np.asarray(deps, float)
    ar = 0.5 * (deps[1:] + deps[:-1]) / W
    nr = np.diff(deps) / np.diff(DURS); nr = nr / nr[0]
    return ar, nr


for n_fp in [1, 3, 6]:
    deps = [depth(d, n_fp) for d in DURS]
    ar, nr = arde(deps)
    at = np.interp(VPS_AR, ar, nr)
    rmse = float(np.sqrt(np.mean((at - VPS_NR) ** 2)))
    print(f"\n  n_fp={n_fp}  depths {np.round(deps,2)}", flush=True)
    print(f"    @ ViennaPS AR {VPS_AR}: petch {np.round(at,3)} vs {VPS_NR}  deep={at[-1]:.3f}  RMSE {rmse:.3f}", flush=True)
print("\ndone", flush=True)
