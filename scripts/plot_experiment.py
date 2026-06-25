#!/usr/bin/env python3
"""Experimental ARDE: petch vs the de Boer/Blauw measured wafer (the real thing, not another simulator).
Shows (a) petch with ballistic/ViennaPS-class params -- gentle, like any ballistic code -- sitting ABOVE
the wafer, and (b) petch with the de Boer cryo process params reproducing the measured ARDE. Writes
viz/experiment_arde.png + docs copy."""
import os, numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = np.load(os.path.join(os.path.dirname(__file__), "..", "deboer_final.npz"))
ar_d, nr_d, ar_b, nr_b = d['ar_d'], d['nr_d'], d['ar_b'], d['nr_b']
exp_ar, exp_r = d['exp_ar'], d['exp_r']
sm = lambda y: uniform_filter1d(y, 5, mode="nearest")        # display-smooth the per-point MC noise

fig, ax = plt.subplots(figsize=(8.4, 5.4))
ax.plot(exp_ar, exp_r, "k*", ms=20, label="de Boer wafer (measured)", zorder=6)
ax.plot(exp_ar, exp_r, "k--", lw=1.3, alpha=0.55, zorder=2)
m = (ar_d <= 28) & np.isfinite(nr_d)
ax.plot(ar_d[m], sm(nr_d)[m], "-", color="#c0392b", lw=2.8, label="ballistic model (petch ≈ ViennaPS)", zorder=4)
mb = (ar_b <= 20) & np.isfinite(nr_b)
ax.plot(ar_b[mb], sm(nr_b)[mb], "-", color="#2471c7", lw=2.8, label="petch + de Boer cryo process params", zorder=5)
ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylabel("normalized etch rate  $n_r$  (1 = open field)")
ax.set_xlim(0, 28); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=11)
ax.set_title("ARDE vs the real wafer — de Boer/Blauw cryo SF$_6$/O$_2$ DRIE", fontsize=13.5, fontweight="bold")
ax.annotate("ballistic transport is\ntoo gentle vs the wafer", xy=(12.5, 0.78), xytext=(15.5, 0.62),
            fontsize=10, color="#c0392b", ha="left", arrowprops=dict(arrowstyle="->", color="#c0392b"))
ax.annotate("with the cryo process\nparams, petch matches\nthe measurement (RMSE 0.04)", xy=(10, 0.45),
            xytext=(1.0, 0.18), fontsize=10, color="#2471c7",
            arrowprops=dict(arrowstyle="->", color="#2471c7"))
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "viz"); os.makedirs(out, exist_ok=True)
for p in [os.path.join(out, "experiment_arde.png"), os.path.join(os.path.dirname(__file__), "..", "docs", "experiment_arde.png")]:
    plt.savefig(p, dpi=150); print("saved", p)
print("\nat AR 10/20: experiment", exp_r[1], exp_r[2],
      "| ballistic", round(float(np.interp(10, ar_d, nr_d)), 2), round(float(np.interp(20, ar_d, nr_d)), 2),
      "| petch-cryo", round(float(np.interp(10, ar_b, nr_b)), 2), round(float(np.interp(20, ar_b, nr_b)), 2))
