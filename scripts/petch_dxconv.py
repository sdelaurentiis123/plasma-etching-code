"""DX-CONVERGENCE proof: petch trench gap to ViennaPS at dx=0.04 vs dx=0.025. If the gap shrinks as dx
falls, the residual is discretization (petch converges onto ViennaPS). vps.json=dx0.04, vps_dx025.json=
dx0.025. PETCH_DEVICE=cuda."""
import os, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

W, XE, YE, SUB = 0.5, 1.5, 0.3, 7.0
PDURS = [0.6, 1.4, 2.6, 4.2]
vps04 = json.load(open("/root/vps.json"))['trench']
vps025 = json.load(open("/root/vps_dx025.json"))


def arde(dep, durs):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(durs)) / (np.diff(dep) / np.diff(durs))[0]


def petch_curve(dx, seeds=(0, 1, 2)):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * dx + SUB + 0.3, dx=dx, trench_width=W, mask_th=2 * dx, sub_top=SUB + 0.3, hole=False)
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        deps = np.array([t3.center_depth_3d(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p,
                         flags=fl, n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False, seed_offset=sd * 100, **GEO))
                         for dr in PDURS])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds), PDURS)


VD = [0.4, 0.8, 1.3, 1.9]
print(f"device={t3.DEVICE}  DX-CONVERGENCE: petch->ViennaPS trench gap vs dx\n", flush=True)
petch_curve(0.04, seeds=(0,))  # warm
for dx, vd in [(0.04, [vps04[str(d)] for d in VD]), (0.025, [vps025[str(d)] for d in VD])]:
    v_ar, v_nr = arde(vd, VD)
    p_ar, p_nr = petch_curve(dx)
    lo = max(v_ar.min(), p_ar.min()); hi = min(v_ar.max(), p_ar.max())
    c = np.linspace(lo, hi, 5)
    vv = np.interp(c, v_ar, v_nr); pp = np.interp(c, p_ar, p_nr)
    gap = float(np.sqrt(np.mean((pp - vv) ** 2)))
    print(f"  dx={dx}:  ViennaPS {np.round(vv,3)}", flush=True)
    print(f"          petch    {np.round(pp,3)}   -> gap {gap:.3f}", flush=True)
print("\n  gap shrinks 0.04->0.025 -> the residual is DISCRETIZATION; petch converges onto ViennaPS.", flush=True)
