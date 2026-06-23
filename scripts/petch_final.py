"""Clean final petch-vs-ViennaPS trench ARDE: faithful ion reflection + converged coverage (n_fp=3),
AVERAGED over seeds to beat petch's deep-floor MC noise. Reports nr at ViennaPS AR [3.7,6.1,8.6]
(ground truth 1, 0.861, 0.731 -- deterministic) + per-point std + wall-clock. cuda."""
import os, time
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

DX, W, XE, YE, SUB = 0.04, 0.5, 1.5, 0.3, 6.0
DURS = [0.7, 1.1, 1.5, 1.9, 2.3, 2.7]
SEEDS = [0, 100, 200]
VPS_AR = np.array([3.7, 6.1, 8.6]); VPS_NR = np.array([1.0, 0.861, 0.731])


def depth(dur, seed):
    GEO = dict(Lx=XE, Ly=YE, Lz=2 * DX + SUB + 0.3, dx=DX, trench_width=W,
               mask_th=2 * DX, sub_top=SUB + 0.3, hole=False)
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1; p['n_fp'] = 3
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc", ion_reflection=True)
    t0 = time.time()
    g = t3.run_etch_3d(t_end=dur, n_steps=max(8, int(dur * 22)), par=p, flags=fl, n_ion=60000, n_neu=60000,
                       reinit_method="fsm", verbose=False, seed_offset=seed, **GEO)
    return t3.center_depth_3d(g), time.time() - t0


def nr_curve(seed):
    res = [depth(d, seed) for d in DURS]
    deps = np.array([r[0] for r in res]); wall = sum(r[1] for r in res)
    ar = 0.5 * (deps[1:] + deps[:-1]) / W
    nr = np.diff(deps) / np.diff(DURS); nr = nr / nr[0]
    return np.interp(VPS_AR, ar, nr), wall


curves = []; walls = []
for s in SEEDS:
    at, wall = nr_curve(s)
    curves.append(at); walls.append(wall)
    print(f"  seed {s}: nr {np.round(at,3)}  (sweep wall {wall:.1f}s)", flush=True)
curves = np.array(curves)
mean = curves.mean(0); std = curves.std(0)
rmse = float(np.sqrt(np.mean((mean - VPS_NR) ** 2)))
print(f"\n  ViennaPS (truth): {VPS_NR}", flush=True)
print(f"  petch mean:       {np.round(mean,3)}  +/- {np.round(std,3)}", flush=True)
print(f"  deep(AR8.6): petch {mean[-1]:.3f} +/- {std[-1]:.3f}  vs VPS 0.731", flush=True)
print(f"  RMSE {rmse:.3f}   (was 0.152 before faithful ion reflection)", flush=True)
print(f"  avg sweep wall {np.mean(walls):.1f}s over {len(DURS)} depths", flush=True)
import json
json.dump(dict(vps=VPS_NR.tolist(), petch_mean=mean.tolist(), petch_std=std.tolist(),
               rmse=rmse, ar=VPS_AR.tolist()), open("/root/petch_final.json", "w"))
print("wrote /root/petch_final.json", flush=True)
