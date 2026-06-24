#!/usr/bin/env python3
"""Does ViennaPS-style ray scaling (rays_per_point) fix petch's grid-sensitivity? Etch the SAME trench
at dx=0.25/0.15/0.10 with (a) the old FIXED ray budget and (b) rays_per_point (rays scale with #facets),
report deep nr + max AR reached + wall time. Clean fix => nr converges monotonically (no stall) AND stays
fast. PETCH_DEVICE=cuda."""
import os, time
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

W = 1.0
GEO = dict(Lx=6.0, Ly=1.5, sub_top=14.0, hole=False)   # Ly=1.5 so even dx=0.25 has 6 y-cells (MC needs >=2)


def run(dx, rpp=None, nfix=40000):
    Lz = 2*dx + 14.0 + 0.3
    g = dict(GEO); g['Lz'] = Lz
    p = dict(petch.PAR); p['rate_scale'] = 0.1; p['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", neutral_transport="mc", ion_reflection=True)
    NS = 70
    t0 = time.time()
    geo = t3.run_etch_3d(trench_width=W, dx=dx, mask_th=2*dx, t_end=3.0, n_steps=NS, par=p, flags=fl,
                         n_ion=nfix, n_neu=nfix, rays_per_point=rpp, reinit_method="fsm", verbose=False,
                         record_depth_every=3, **g)
    wall = time.time() - t0
    h = geo['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    tm = st / NS * 3.0; r = np.gradient(dd, tm); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max(); nr = np.clip(r / max(r0, 1e-9), 0, 1.5)
    nr8 = float(np.interp(8, ar, nr)) if ar.max() >= 8 else float('nan')
    return nr8, float(ar.max()), wall


print(f"device={t3.DEVICE}\n", flush=True)
print("dx     | FIXED budget (40k)          | rays_per_point=80", flush=True)
print("       | nr@AR8  maxAR  wall          | nr@AR8  maxAR  wall", flush=True)
for dx in [0.25, 0.15, 0.10]:
    a8, am, aw = run(dx, rpp=None, nfix=40000)
    b8, bm, bw = run(dx, rpp=80)
    print(f"{dx:.2f}   | {a8:.3f}   {am:4.1f}  {aw:5.1f}s   | {b8:.3f}   {bm:4.1f}  {bw:5.1f}s", flush=True)
print("\n  CLEAN FIX => rays_per_point nr@AR8 stable across dx + maxAR doesn't collapse at fine dx.", flush=True)
print("  (fixed budget: nr drifts up + fine dx stalls/under-samples)", flush=True)
