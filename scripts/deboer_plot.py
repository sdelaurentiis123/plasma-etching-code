#!/usr/bin/env python3
"""Plot the de Boer trench ARDE: experiment vs petch (periodic-y on/off, MC + radiosity). Shows the
periodic-y fix (un-stall MC / kill radiosity cliff) and the residual transport gap. PETCH_DEVICE=cuda."""
import os
os.environ.setdefault("PETCH_DEVICE", "cuda")
import numpy as np
import petch
from petch import threed as t3

EXP_AR = np.array([0.0, 10.0, 20.0, 40.0]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])
W = 2.0; GEO = dict(Lx=10, Ly=5, Lz=72, mask_th=2, sub_top=66, hole=False, t_end=10.0)
DX, NS = 0.25, 140


def curve(nt, periodic, beta=0.47):
    p = dict(petch.PAR); p['rate_scale'] = 0.30; p['betaE'] = beta; p['periodic_y'] = periodic
    fl = petch.Flags(coverage_sticking=True, neutral_transport=nt, sampling="sobol", warm_start_coverage=True)
    g = t3.run_etch_3d(trench_width=W, dx=DX, n_steps=NS, par=p, flags=fl, n_ion=30000, n_neu=30000,
                       reinit_method="fsm", verbose=False, record_depth_every=5, **GEO)
    h = g['depth_history']; st = np.array([x[0] for x in h]); dd = np.array([x[1] for x in h])
    t = st / NS * GEO['t_end']; r = np.gradient(dd, t); ar = dd / W
    r0 = r[ar < 2].max() if (ar < 2).any() else r.max()
    return ar, np.clip(r / max(r0, 1e-9), 0, 1.6)


configs = [("mc", 0, "petch MC (no periodic)", "#9aa0a6", "--"),
           ("mc", 1, "petch MC + periodic-y", "#1a73e8", "-"),
           ("radiosity", 1, "petch radiosity + periodic-y", "#7a3b2e", "-")]
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.5, 5))
arr = np.linspace(0, 40, 100)
ax.plot(arr, np.interp(arr, EXP_AR, EXP_R), 'k-', lw=2.5, label='de Boer/Blauw experiment')
ax.plot(EXP_AR, EXP_R, 'ko', ms=8)
for nt, pk, lab, col, ls in configs:
    ar, nr = curve(nt, pk)
    o = np.argsort(ar)
    ax.plot(ar[o], nr[o], ls, color=col, lw=1.8, label=lab, alpha=0.9)
    print(f"{lab}: AR_max {ar.max():.1f}", flush=True)
ax.set_xlim(0, 40); ax.set_ylim(0, 1.05); ax.set_xlabel("aspect ratio"); ax.set_ylabel("normalized bottom etch rate")
ax.set_title("Trench ARDE: petch vs de Boer experiment (periodic-y fix)"); ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("deboer_plot.png", dpi=140)
print("wrote deboer_plot.png", flush=True)
