"""petch (fudge-free, wrap-fixed) head-to-head vs the ViennaPS reference in /root/vps.json, SAME box,
SAME geometry. Trench + hole ARDE + speed. PETCH_DEVICE=cuda."""
import os, time, json
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

vps = json.load(open("/root/vps.json"))
TR = dict(dx=0.04, W=0.5, XE=1.5, YE=0.3, durs=[0.4, 0.8, 1.3, 1.9], SUB=7.0)
HO = dict(dx=0.05, W=0.5, XE=1.5, durs=[0.4, 1.0, 1.8, 2.8], SUB=7.3)


def arde(dep, durs, W):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(durs)) / (np.diff(dep) / np.diff(durs))[0]


def petch_curve(hole, seeds=(0, 1, 2)):
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
                         for dr in g['durs']])
        accd = deps if accd is None else accd + deps
    return arde(accd / len(seeds), g['durs'], g['W'])


print(f"device={t3.DEVICE}  PETCH (fudge-free, wrap-fixed) vs ViennaPS -- SAME BOX, SAME GEOMETRY\n", flush=True)
petch_curve(False, seeds=(0,))  # warm

for hole, g, lab in [(False, TR, "TRENCH"), (True, HO, "HOLE")]:
    vkey = 'hole' if hole else 'trench'
    vdep = [vps[vkey][str(dr)] if str(dr) in vps[vkey] else vps[vkey][dr] for dr in g['durs']]
    v_ar, v_nr = arde(vdep, g['durs'], g['W'])
    ar, nr = petch_curve(hole)
    pp = np.interp(v_ar, ar, nr)
    gap = float(np.sqrt(np.mean((pp - v_nr) ** 2)))
    print(f"  {lab}:", flush=True)
    print(f"    ViennaPS nr {np.round(v_nr,3)} @ AR {np.round(v_ar,2)}", flush=True)
    print(f"    petch    nr {np.round(pp,3)} @ same AR   -> gap {gap:.3f}", flush=True)

# speed: petch wall time on the deep trench + hole, vs ViennaPS json timings
def ptime(hole):
    g = HO if hole else TR
    GEO = dict(Lx=g['XE'], Ly=(g['XE'] if hole else g['YE']), Lz=2 * g['dx'] + g['SUB'] + 0.3, dx=g['dx'],
               trench_width=g['W'], mask_th=2 * g['dx'], sub_top=g['SUB'] + 0.3, hole=hole)
    p = dict(petch.PAR); p['rate_scale'] = 0.1
    if not hole: p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
    dr = g['durs'][-1]
    t0 = time.time()
    t3.run_etch_3d(t_end=dr, n_steps=max(8, int(dr * 22)), par=p, flags=fl, n_ion=60000, n_neu=60000,
                   reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0


print("\n  SPEED (deepest point wall time, same box):", flush=True)
pt = ptime(False); print(f"    trench: petch {pt:.1f}s  vs  ViennaPS-CPU {vps['speed']['trench_dur1.9_s']:.1f}s  = {vps['speed']['trench_dur1.9_s']/pt:.0f}x", flush=True)
ph = ptime(True); print(f"    hole:   petch {ph:.1f}s  vs  ViennaPS-CPU {vps['speed']['hole_dur2.8_s']:.1f}s  = {vps['speed']['hole_dur2.8_s']/ph:.0f}x", flush=True)
