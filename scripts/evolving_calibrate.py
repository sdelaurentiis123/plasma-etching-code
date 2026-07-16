"""Legacy calibration to the fitted Blauw/Clausing model curve, not de Boer pixels.

Retained for reproducibility only.  The 1/.43/.29 sequence is exposed calculated-model output and
cannot support a wafer-validation or held-out claim.
"""
import os; os.environ.setdefault("PETCH_DEVICE", "cpu")
import numpy as np
from scipy.ndimage import uniform_filter1d
import petch
from petch import threed as t3
W, MASK = 2.0, 0.5
FIELD = 57.6 * 0.0226
MODEL_CURVE = {10.0: 0.43, 20.0: 0.29}
for wls in [1.9, 2.4, 2.9]:
    par = dict(petch.PAR); par['periodic_y'] = 1; par['rate_scale'] = 0.0226
    par['knudsen_wall_loss_scale'] = wls
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     sampling="sobol", ion_reflection=True, neutral_transport='knudsen',
                     warm_start_coverage=True)
    g = t3.run_etch_3d(Lx=10.0, Ly=5.0, Lz=52.0, dx=0.25, trench_width=W, mask_th=MASK, sub_top=49.0,
                       t_end=85.0, n_steps=300, hole=False, par=par, flags=fl, n_ion=40000, n_neu=40000,
                       reinit_method="fsm", verbose=False, record_depth_every=1)
    h = g['depth_history']; dd = np.maximum.accumulate(np.array([x[1] for x in h]))
    tm = np.array([x[0] for x in h]) / 300 * 85.0
    dd_s = uniform_filter1d(dd, 15, mode="nearest")
    rate = uniform_filter1d(np.gradient(dd_s, tm), 15, mode="nearest")
    ar_eff = (dd_s + MASK) / W
    out = {}
    for A, wnr in MODEL_CURVE.items():
        out[A] = float(np.interp(A, ar_eff, rate)) / FIELD if A <= ar_eff.max() - 0.5 else np.nan
    errs = [out[a] - MODEL_CURVE[a] for a in MODEL_CURVE if np.isfinite(out[a])]
    rmse = float(np.sqrt(np.mean(np.array(errs) ** 2))) if errs else np.nan
    gate = "PASS" if rmse <= 0.05 else "fail"
    print(f"wls={wls}: nr@10={out[10.0]:.3f} (model .43) nr@20={out[20.0]:.3f} (model .29)  legacy-RMSE={rmse:.3f} [{gate}]  (max AR {ar_eff.max():.1f})", flush=True)
print("DONE", flush=True)
