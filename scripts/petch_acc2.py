"""Fair AR-matched accuracy: petch (fudge-free, wrap-fixed) vs ViennaPS (/root/vps.json). Fudge-free petch
etches slower, so give it LONGER durations to span the same AR range, then compare nr-vs-AR on COMMON AR
points in the overlap (no extrapolation). PETCH_DEVICE=cuda."""
import os, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

vps = json.load(open("/root/vps.json"))
TR = dict(dx=0.04, W=0.5, XE=1.5, YE=0.3, vdurs=[0.4, 0.8, 1.3, 1.9], pdurs=[0.6, 1.4, 2.6, 4.2], SUB=7.0)
HO = dict(dx=0.05, W=0.5, XE=1.5, vdurs=[0.4, 1.0, 1.8, 2.8], pdurs=[0.6, 1.6, 3.0, 4.8], SUB=7.3)


def arde(dep, durs, W):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(durs)) / (np.diff(dep) / np.diff(durs))[0]


def petch_curve(hole, durs, seeds=(0, 1, 2)):
    g = HO if hole else TR
    GEO = dict(Lx=g['XE'], Ly=(g['XE'] if hole else g['YE']), Lz=2 * g['dx'] + g['SUB'] + 0.3, dx=g['dx'],
               trench_width=g['W'], mask_th=2 * g['dx'], sub_top=g['SUB'] + 0.3, hole=hole)
    accd = None
    for sd in seeds:
        p = dict(petch.PAR); p['rate_scale'] = 0.1
        if not hole:
            p['periodic_y'] = 1
        fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                         warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
        dfun = t3.max_depth_3d if hole else t3.center_depth_3d
        deps = np.array([dfun(t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl,
                         n_ion=60000, n_neu=60000, reinit_method="fsm", verbose=False, seed_offset=sd * 100, **GEO))
                         for dr in durs])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds), durs, g['W'])


print(f"device={t3.DEVICE}  AR-MATCHED accuracy: petch (fudge-free) vs ViennaPS, SAME box+geometry\n", flush=True)
petch_curve(False, TR['pdurs'], seeds=(0,))  # warm

for hole, g, lab in [(False, TR, "TRENCH"), (True, HO, "HOLE")]:
    vkey = 'hole' if hole else 'trench'
    vdep = [vps[vkey][str(dr)] for dr in g['vdurs']]
    v_ar, v_nr = arde(vdep, g['vdurs'], g['W'])
    p_ar, p_nr = petch_curve(hole, g['pdurs'])
    lo = max(v_ar.min(), p_ar.min()); hi = min(v_ar.max(), p_ar.max())
    common = np.linspace(lo, hi, 5)
    vv = np.interp(common, v_ar, v_nr); pp = np.interp(common, p_ar, p_nr)
    gap = float(np.sqrt(np.mean((pp - vv) ** 2)))
    print(f"  {lab} (overlap AR {lo:.1f}-{hi:.1f}):", flush=True)
    print(f"    AR       {np.round(common,1)}", flush=True)
    print(f"    ViennaPS {np.round(vv,3)}", flush=True)
    print(f"    petch    {np.round(pp,3)}   -> gap {gap:.3f}", flush=True)
