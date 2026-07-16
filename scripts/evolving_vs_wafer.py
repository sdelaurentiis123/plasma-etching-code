"""Legacy replay against the fitted Blauw/Clausing AR curve; not experimental validation.

The historical 1/.43/.29/.20 sequence was mislabeled as direct de Boer wafer data.  It is a
calculated model curve and every point is exposed.  This script remains only to reproduce the old
Knudsen calibration behavior; use the checksummed Figure-9 development runners for direct pixels.
"""
import os; os.environ.setdefault("PETCH_DEVICE", "cpu")
import time, numpy as np
from scipy.ndimage import uniform_filter1d
import petch
from petch import threed as t3
W, MASK = 2.0, 0.5
FIELD = 57.6 * 0.0226
MODEL_CURVE = {10.0: 0.43, 20.0: 0.29, 40.0: 0.20}
for seed in (0, 101):
    par = dict(petch.PAR); par['periodic_y'] = 1; par['rate_scale'] = 0.0226
    par['knudsen_wall_loss_scale'] = 2.9
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     sampling="sobol", ion_reflection=True, neutral_transport='knudsen',
                     warm_start_coverage=True)
    t0 = time.time()
    g = t3.run_etch_3d(Lx=10.0, Ly=5.0, Lz=90.0, dx=0.25, trench_width=W, mask_th=MASK, sub_top=86.0,
                       t_end=380.0, n_steps=420, hole=False, par=par, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, record_depth_every=1, seed_offset=seed)
    h = g['depth_history']; dd = np.maximum.accumulate(np.array([x[1] for x in h]))
    tm = np.array([x[0] for x in h]) / 420 * 380.0
    dd_s = uniform_filter1d(dd, 21, mode="nearest")
    rate = uniform_filter1d(np.gradient(dd_s, tm), 21, mode="nearest")
    ar_eff = (dd_s + MASK) / W
    vals = {A: (float(np.interp(A, ar_eff, rate)) / FIELD if A <= ar_eff.max() - 0.5 else np.nan) for A in MODEL_CURVE}
    tag = " | ".join(
        f"AR{int(A)}={vals[A]:.3f}(model {MODEL_CURVE[A]})" for A in MODEL_CURVE)
    print(f"seed{seed}: {tag}  maxAR={ar_eff.max():.1f}  ({time.time()-t0:.0f}s)", flush=True)
print("DONE", flush=True)
