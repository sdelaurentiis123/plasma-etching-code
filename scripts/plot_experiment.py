#!/usr/bin/env python3
"""Experimental ARDE: the measured de Boer/Blauw wafer vs PURE ViennaPS (run on the same W=2µm trench
geometry) and petch with the cryo process params. Shows the empirical data points, the open-source SOTA
sitting too gentle above the wafer, and petch reproducing the measurement. Writes viz/ + docs/."""
import os, numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
db = np.load(os.path.join(here, "..", "deboer_final.npz"))
ar_b, nr_b = db['ar_b'], db['nr_b']                       # petch + de Boer cryo process params
exp_ar, exp_r = db['exp_ar'], db['exp_r']                 # measured wafer
vd = np.load(os.path.join(here, "..", "vps_deboer.npz"))  # PURE ViennaPS on the W=2µm de Boer trench
W = float(vd['W']); dur, dep = vd['dur'], vd['dep']
vps_ar = 0.5*(dep[1:]+dep[:-1])/W; vps_nr = np.diff(dep)/np.diff(dur); vps_nr = vps_nr/vps_nr[0]
sm = lambda y: uniform_filter1d(y, 5, mode="nearest")

fig, ax = plt.subplots(figsize=(8.4, 5.4))
ax.plot(exp_ar, exp_r, "k*", ms=20, label="de Boer wafer (measured)", zorder=6)
ax.plot(exp_ar, exp_r, "k--", lw=1.3, alpha=0.55, zorder=2)
ax.plot(vps_ar, vps_nr, "s-", color="#c0392b", lw=2.6, ms=8, label="ViennaPS (GPU, run pure)", zorder=4)
mb = (ar_b <= 20) & np.isfinite(nr_b)
ax.plot(ar_b[mb], sm(nr_b)[mb], "-", color="#2471c7", lw=2.8, label="petch + de Boer cryo process params", zorder=5)
ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylabel("normalized etch rate  $n_r$  (1 = open field)")
ax.set_xlim(0, 28); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=11)
ax.set_title("ARDE vs the real wafer — de Boer/Blauw cryo SF$_6$/O$_2$ DRIE (W = 2 µm trench)",
             fontsize=13, fontweight="bold")
ax.annotate("ViennaPS (and any ballistic model)\nis too gentle vs the real wafer", xy=(10, 0.6), xytext=(13, 0.74),
            fontsize=10, color="#c0392b", ha="left", arrowprops=dict(arrowstyle="->", color="#c0392b"))
ax.annotate("with the cryo process params\npetch matches the measurement", xy=(10, 0.45), xytext=(0.8, 0.16),
            fontsize=10, color="#2471c7", arrowprops=dict(arrowstyle="->", color="#2471c7"))
plt.tight_layout()
for p in [os.path.join(here, "..", "viz", "experiment_arde.png"), os.path.join(here, "..", "docs", "experiment_arde.png")]:
    plt.savefig(p, dpi=150); print("saved", p)
print(f"\nAR10: wafer {exp_r[1]} | ViennaPS {np.interp(10,vps_ar,vps_nr):.2f} | petch-cryo {np.interp(10,ar_b,nr_b):.2f}")
