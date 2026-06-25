#!/usr/bin/env python3
"""Clean ARDE accuracy figure: petch vs ViennaPS, trench (left) + hole (right), nr vs AR, with the
de Boer real-wafer reference for context. Filters the shallow-depth metric transient. Reads
accuracy_sweep.npz (argv[1]); writes viz/accuracy_sweep.png (+ docs copy)."""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W = 0.5
AR_REF = 3.0
SUB = 6.3
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "accuracy_sweep.npz")
EXP_AR = np.array([0., 10., 20., 40.]); EXP_R = np.array([1.0, 0.43, 0.29, 0.20])  # de Boer wafer


def nr_curve(dur, dep):
    keep = dep < SUB - 0.5
    dur, dep = dur[keep], dep[keep]
    ar = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dur)
    ref = np.interp(AR_REF, ar, rate) if ar.min() <= AR_REF <= ar.max() else rate[0]
    nr = rate / ref
    good = (ar >= AR_REF - 0.6) & (ar <= 7.5) & (nr > 0) & (nr < 1.25)   # clean common range; drop noisy tail
    return ar[good], nr[good]


plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.25})
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
for ax, feat, title in [(axes[0], "trench", "TRENCH  (line / slot)"), (axes[1], "hole", "HOLE  (round via)")]:
    # real-wafer reference (context: both simulators sit above reality)
    ax.plot(EXP_AR, EXP_R, "k*--", ms=12, lw=1.2, alpha=0.55, label="real wafer (de Boer)", zorder=2)
    for eng, color, mk in [("vps", "#c0392b", "s"), ("petch", "#2471c7", "o")]:
        kd, kp = f"{eng}_{feat}_dur", f"{eng}_{feat}_dep"
        if kd in d and kp in d:
            ar, nr = nr_curve(d[kd], d[kp])
            ax.plot(ar, nr, mk + "-", color=color, lw=2.6, ms=8, label=("ViennaPS" if eng == "vps" else "petch"), zorder=5)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("aspect ratio  (depth / width)")
    ax.set_xlim(0, 9); ax.set_ylim(0, 1.08); ax.legend(loc="upper right", framealpha=0.95)
axes[0].set_ylabel("normalized etch rate  $n_r$")
fig.suptitle("Aspect-ratio-dependent etch rate: petch vs ViennaPS  —  SF$_6$/O$_2$, W = 0.5 µm",
             fontsize=14, fontweight="bold")
axes[0].text(0.04, 0.06, "petch ≈ ViennaPS on holes;\nslightly gentler on trenches.\nBoth above the real wafer.",
             transform=axes[0].transAxes, fontsize=9.5, va="bottom",
             bbox=dict(boxstyle="round,pad=0.4", fc="#f3f1ea", ec="#b9b29c", alpha=0.9))
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "viz"); os.makedirs(out, exist_ok=True)
for p in [os.path.join(out, "accuracy_sweep.png"),
          os.path.join(os.path.dirname(__file__), "..", "docs", "accuracy_sweep.png")]:
    plt.savefig(p, dpi=150); print("saved", p)

print("\nnr at AR = 4, 6:")
for feat in ["hole", "trench"]:
    for eng in ["vps", "petch"]:
        kd = f"{eng}_{feat}_dur"
        if kd in d:
            ar, nr = nr_curve(d[kd], d[f"{eng}_{feat}_dep"])
            v = [f"{np.interp(a, ar, nr):.2f}" if ar.max() >= a else "  - " for a in (4, 6)]
            print(f"  {feat:6s} {eng:5s}: AR4={v[0]} AR6={v[1]}")
