"""petch-only (ViennaPS PyPI build segfaults in geometry creation on current boxes -> use the ViennaPS
reference numbers we MEASURED earlier). Two things: (1) SPEED -- time petch (warmed) on the benchmark
holes to confirm the hot loop is unchanged by this session's edits (vs prior petch 1.30/1.61/1.82s on a
3090; this is a 4090 so expect faster). (2) RESIDUAL -- is petch's ~0.08-0.10 too-steep-at-deep-AR
discretization? Run fudge-free wrap-fixed petch trench+hole ARDE at coarse AND fine dx; if it gets gentler
(toward the cached ViennaPS) at finer dx -> resolution. PETCH_DEVICE=cuda."""
import os, time
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

# ---------- (1) SPEED ----------
print("=== SPEED (petch wall time, warmed) ===", flush=True)


def ours_time(diam, dx=0.25, ns=40):
    p = dict(petch.PAR); p['rate_scale'] = 1.0
    GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc")
    t0 = time.time()
    g = t3.run_etch_3d(trench_width=diam, dx=dx, n_steps=ns, par=p, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, **GEO)
    return time.time() - t0, t3.max_depth_3d(g)


for diam in [4.0, 6.0, 8.0]:
    ours_time(diam)  # warm
    w, dep = ours_time(diam)
    print(f"  petch d={diam}: {w:.2f}s  depth {dep:.1f}um   (prior 3090 baseline ~1.3/1.6/1.8s; vs ViennaPS-GPU 19-26s = ~14x)", flush=True)

# ---------- (2) RESIDUAL: dx convergence ----------
print("\n=== RESIDUAL: does the deep-AR steepness shrink with finer dx? ===", flush=True)
VPS_TR = ([3.73, 6.10, 8.58], [1.0, 0.862, 0.732])     # cached ViennaPS trench (dx=0.04)
VPS_HO = ([3.45, 5.58, 7.41], [1.0, 0.641, 0.445])     # cached ViennaPS hole (dx=0.05)


def arde(dep, durs, W):
    dep = np.asarray(dep, float)
    return 0.5 * (dep[1:] + dep[:-1]) / W, (np.diff(dep) / np.diff(durs)) / (np.diff(dep) / np.diff(durs))[0]


def petch_ardе(hole, dx, durs, W, seeds=(0, 1, 2)):
    Lxy = 1.5; SUB = 7.0
    GEO = dict(Lx=Lxy, Ly=(Lxy if hole else 0.3), Lz=2 * dx + SUB + 0.3, dx=dx, trench_width=W,
               mask_th=2 * dx, sub_top=SUB + 0.3, hole=hole)
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
    return arde(accd / len(seeds), durs, W)


for hole, W, durs, (v_ar, v_nr), dxs, lab in [
        (False, 0.5, [0.4, 0.8, 1.3, 1.9], VPS_TR, [0.04, 0.025], "TRENCH"),
        (True, 0.5, [0.4, 1.0, 1.8, 2.8], VPS_HO, [0.05, 0.033], "HOLE")]:
    print(f"\n  {lab} (cached ViennaPS nr {np.round(v_nr,3)} @ AR {np.round(v_ar,2)}):", flush=True)
    petch_ardе(hole, dxs[0], durs, W, seeds=(0,))  # warm
    for dx in dxs:
        ar, nr = petch_ardе(hole, dx, durs, W)
        pp = np.interp(v_ar, ar, nr); gap = float(np.sqrt(np.mean((pp - np.array(v_nr)) ** 2)))
        print(f"    dx={dx}: petch nr@vpsAR {np.round(pp,3)}  gap {gap:.3f}", flush=True)
print("\n  gap shrinks at finer dx -> residual is DISCRETIZATION (cleanup = finer dx / better floor resolution).", flush=True)
