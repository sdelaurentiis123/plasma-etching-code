#!/usr/bin/env python3
"""WAVE 2 GATE (G2): is the GPU Jacobi Godunov-Eikonal reinit (reinit_fsm) CORRECT?

Two checks, both runnable on the Warp CPU backend (no GPU needed for correctness):
  1. |grad phi| ~ 1 inside the band (the property the old PDE reinit_gpu got WRONG: |grad|=1.32 on
     masked fronts). Compare FSM vs skfmm-narrow on a deep-hole phi with a re-pinned mask.
  2. The phi=0 contour does NOT move (sub-cell distance preserved) -> a full etch lands at the same
     depth as skfmm-narrow within MC noise.

Run locally: PETCH_DEVICE=cpu python scripts/reinit_correctness_3d.py
"""
import os
os.environ.setdefault("PETCH_DEVICE", "cpu")
import numpy as np
import petch
from petch import threed as t3

dx = 0.25
# build a deep-hole geometry, etch a few steps to get a realistic warped front + re-pinned mask
GEO = dict(Lx=8, Ly=8, Lz=16, mask_th=2, sub_top=12, hole=True)
geo = t3.make_trench_3d(GEO['Lx'], GEO['Ly'], GEO['Lz'], dx, 4.0, GEO['mask_th'], GEO['sub_top'], hole=True)
phi0 = geo['phi'].copy()
# carve a non-trivial front: drop the central column region (gas) to mimic a partially etched hole
xs, ys, zs = geo['xs'], geo['ys'], geo['zs']
X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
r = np.sqrt((X - GEO['Lx']/2)**2 + (Y - GEO['Ly']/2)**2)
etched = (r < 2.0) & (Z < GEO['sub_top']) & (Z > 4.0)
phi0[etched] = np.abs(phi0[etched]) * -1.0 - 0.3      # make it gas, slightly negative
phi0 = t3.skfmm.distance(phi0, dx=dx)                  # clean SDF as the common starting point

band = 4 * dx + 2.0 * dx


def grad_in_band(phi):
    gx, gy, gz = np.gradient(phi, dx)
    g = np.sqrt(gx**2 + gy**2 + gz**2)
    m = np.abs(phi) < (band - dx)        # interior of the band (avoid the clamp edge)
    return g[m]


ref = t3.reinit_narrow(phi0, dx, band)
fsm = t3.reinit_fsm(phi0, dx, band)

gr = grad_in_band(ref); gf = grad_in_band(fsm)
# agreement of the signed distance inside the band
m = np.abs(ref) < (band - dx)
dd = np.abs(fsm[m] - ref[m])
# contour preservation: zero-crossing locations should match (sign agreement in the band)
sign_match = np.mean(np.sign(fsm[m]) == np.sign(ref[m]))

print(f"device={t3.DEVICE}  dx={dx}  band={band:.2f}  cells_in_band={m.sum()}")
print(f"  |grad| in band:   skfmm  mean={gr.mean():.3f} p95={np.percentile(gr,95):.3f}")
print(f"                    FSM    mean={gf.mean():.3f} p95={np.percentile(gf,95):.3f}  (target ~1.00)")
print(f"  FSM vs skfmm distance in band:  max|dphi|={dd.max():.3f}  mean|dphi|={dd.mean():.3f}  ({dd.max()/dx:.2f} cells)")
print(f"  sign agreement in band: {100*sign_match:.2f}%  (contour preserved)")

ok_grad = abs(gf.mean() - 1.0) < 0.05 and np.percentile(gf, 95) < 1.15
ok_dist = dd.max() < 1.0 * dx          # within ~1 cell of skfmm everywhere in band
ok_sign = sign_match > 0.999
print(f"\n  GATE: |grad|~1 {'PASS' if ok_grad else 'FAIL'} | dist<1cell {'PASS' if ok_dist else 'FAIL'} | "
      f"contour {'PASS' if ok_sign else 'FAIL'}")

# 2. full-etch depth parity (short run, CPU) -- FSM vs skfmm-narrow
print("\n  depth parity (8-step etch, CPU):", flush=True)
p = dict(petch.PAR)
res = {}
for meth in ["skfmm", "fsm"]:
    g = t3.run_etch_3d(trench_width=4.0, dx=dx, n_steps=8, par=p,
                       flags=petch.Flags(coverage_sticking=True, sampling="sobol"),
                       n_ion=6000, n_neu=6000, reinit_method=meth, verbose=False,
                       Lx=8, Ly=8, Lz=16, mask_th=2, sub_top=12, t_end=1.0, hole=True)
    res[meth] = t3.max_depth_3d(g)
    print(f"    {meth:6s} depth={res[meth]:.3f}um", flush=True)
ddep = abs(res['fsm'] - res['skfmm'])
print(f"    depth delta = {ddep:.3f}um  ({'PASS' if ddep < 0.5 else 'CHECK'})")
print(f"\n  OVERALL: {'PASS - FSM ready for the box speed run' if (ok_grad and ok_dist and ok_sign and ddep<0.5) else 'NOT READY - keep skfmm-narrow'}")
