#!/usr/bin/env python3
"""Isolate reinit FIDELITY (no feedback compounding): take one developed front, apply skfmm and
gpu reinit to the SAME phi, compare the resulting SDFs in the near band + the zero-contour drift.

Also a ground-truth test: a known analytic SDF (sphere) perturbed off |grad|=1, reinit, compare to
the exact known distance."""
import numpy as np, petch
from petch import threed as t3

# ---------- (A) ground-truth: perturbed sphere, exact distance known ----------
N = 64; dx = 0.25
xs = (np.arange(N) - N/2) * dx
X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
R = 4.0
exact = np.sqrt(X**2 + Y**2 + Z**2) - R            # true SDF of a sphere
perturbed = exact * (1.0 + 0.4*np.sin(X)*np.cos(Y))  # break |grad|=1, keep zero-contour
band = np.abs(exact) < 4*dx

import skfmm
sk = skfmm.distance(perturbed, dx=dx)
gp = t3.reinit_gpu(perturbed, dx)
print("(A) perturbed-sphere reinit vs EXACT distance (near band):")
print(f"    skfmm   max|err|={np.abs(sk-exact)[band].max():.3f}  mean|err|={np.abs(sk-exact)[band].mean():.3f}")
print(f"    gpu     max|err|={np.abs(gp-exact)[band].max():.3f}  mean|err|={np.abs(gp-exact)[band].mean():.3f}")
print(f"    gpu-nan={np.isnan(gp).any()}")

# ---------- (B) real developed etch front: skfmm vs gpu on the SAME phi ----------
GEO = dict(Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0)
p = dict(petch.PAR); p['rate_scale'] = 0.15
g = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=12, par=p,
                   flags=petch.Flags(coverage_sticking=True),
                   n_ion=20000, n_neu=20000, reinit_method="skfmm", verbose=False, **GEO)
phi_dev = g['phi'].copy()
# de-reinitialize slightly so both methods have real work to do
phi_in = phi_dev * 1.0
sk2 = skfmm.distance(phi_in, dx=0.25)
gp2 = t3.reinit_gpu(phi_in, 0.25)
bnd = np.abs(sk2) < 4*0.25
print("\n(B) developed etch front, skfmm vs gpu on identical phi (near band):")
print(f"    max|skfmm-gpu|={np.abs(sk2-gp2)[bnd].max():.3f}  mean={np.abs(sk2-gp2)[bnd].mean():.3f}  gpu-nan={np.isnan(gp2).any()}")
# zero-contour drift along the center column (hole axis)
ic = N//2 if False else g['phi'].shape[0]//2
jc = g['phi'].shape[1]//2
def zero_cross_z(phi, zs):
    col = phi[ic, jc, :]
    s = np.where(np.diff(np.sign(col)) != 0)[0]
    return zs[s[0]] if len(s) else np.nan
zs = g['zs']
print(f"    center-column lowest zero-cross z: skfmm={zero_cross_z(sk2,zs):.3f}  gpu={zero_cross_z(gp2,zs):.3f}")
