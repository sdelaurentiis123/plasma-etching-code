#!/usr/bin/env python3
"""Experimental ARDE: the measured de Boer/Blauw wafer vs PURE ViennaPS (run on the same W=2µm trench
geometry) and petch's deterministic discrete-ordinates (Sn) transport ("dda") on the same width.
All three are reproducible with current code: ViennaPS from vps_deboer.npz (real GPU run),
petch-DDA from dda_vs_mc_arde.npz (static nr(AR), W=2µm, neutral_transport='dda', ion_reflection).

The honest story this figure tells: both ballistic engines agree with each other and track the wafer
through the knee (petch-DDA 0.58 / ViennaPS 0.60 vs measured 0.43 at AR10), but the wafer's flat
high-AR tail is a structural gap (gas conductance / charging) that no ballistic knob closes.
(An earlier version showed a 'petch-cryo' curve built on the removed cal_F knob — retired.)
Writes viz/ + docs/."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
db = np.load(os.path.join(here, "..", "deboer_final.npz"))
exp_ar, exp_r = db['exp_ar'], db['exp_r']                 # measured wafer (de Boer/Blauw)
vd = np.load(os.path.join(here, "..", "vps_deboer.npz"))  # PURE ViennaPS on the W=2µm de Boer trench
W = float(vd['W']); dur, dep = vd['dur'], vd['dep']
vps_ar = 0.5*(dep[1:]+dep[:-1])/W; vps_nr = np.diff(dep)/np.diff(dur); vps_nr = vps_nr/vps_nr[0]
dd = np.load(os.path.join(here, "..", "dda_vs_mc_arde.npz"))  # petch Sn ("dda"), static, W=2µm
dda_ar, dda_nr = dd['ar2'], dd['petch_dda_w2']

fig, ax = plt.subplots(figsize=(8.4, 5.4))
ax.plot(exp_ar, exp_r, "k*", ms=20, label="de Boer wafer (measured)", zorder=6)
ax.plot(exp_ar, exp_r, "k--", lw=1.3, alpha=0.55, zorder=2)
ax.plot(vps_ar, vps_nr, "s-", color="#c0392b", lw=2.6, ms=8, label="ViennaPS (GPU, run pure)", zorder=4)
ax.plot(dda_ar, dda_nr, "o-", color="#2471c7", lw=2.8, ms=8,
        label="petch S$_n$ transport (deterministic)", zorder=5)
ax.set_xlabel("aspect ratio  (depth / width)"); ax.set_ylabel("normalized etch rate  $n_r$  (1 = open field)")
ax.set_xlim(0, 42); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=11)
ax.set_title("ARDE vs the real wafer — de Boer/Blauw cryo SF$_6$/O$_2$ DRIE (W = 2 µm trench)",
             fontsize=13, fontweight="bold")
ax.annotate("both ballistic engines agree\n(petch-S$_n$ 0.58, ViennaPS 0.60 @ AR 10)\nand track the wafer knee",
            xy=(10, 0.58), xytext=(14, 0.72), fontsize=10, color="#2471c7",
            arrowprops=dict(arrowstyle="->", color="#2471c7"))
ax.annotate("the flat high-AR tail is the structural\nfrontier (conductance / charging) —\nno ballistic knob closes it",
            xy=(30, 0.245), xytext=(17, 0.08), fontsize=10, color="0.25",
            arrowprops=dict(arrowstyle="->", color="0.4"))
plt.tight_layout()
for p in [os.path.join(here, "..", "viz", "experiment_arde.png"), os.path.join(here, "..", "docs", "experiment_arde.png")]:
    plt.savefig(p, dpi=150); print("saved", p)
print(f"\nAR10: wafer {exp_r[1]} | ViennaPS {np.interp(10,vps_ar,vps_nr):.2f} | petch-Sn {np.interp(10,dda_ar,dda_nr):.2f}")
