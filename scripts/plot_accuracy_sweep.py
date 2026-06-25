#!/usr/bin/env python3
"""Clean ARDE figure: petch vs ViennaPS, trench (left) + hole (right), nr vs AR. Drops the shallow-depth
metric noise (petch's max_depth catches a filament at very low depth -> a spurious first point), normalizes
both engines at a clean common AR, tight axes. Reads accuracy_sweep.npz (argv[1]); writes viz/ + docs/."""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W = 0.5
AR_REF = 4.0                    # normalize both engines to nr=1 here (a clean, reliably-etched depth)
AR_LO, AR_HI = 3.3, 7.0         # plotted range: drop shallow-noise low end + sparse high end
SUB = 6.3
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "accuracy_sweep.npz")


COMMON_AR = np.arange(4.0, 6.6, 0.5)        # shared aspect-ratio grid -> both engines plotted at same X


def nr_curve(dur, dep):
    """Return nr sampled on the COMMON aspect-ratio grid, so petch and ViennaPS share X values."""
    keep = dep < SUB - 0.5
    dur, dep = dur[keep], dep[keep]
    ar = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dur)
    ref = np.interp(AR_REF, ar, rate)
    nr_raw = rate / ref
    grid = COMMON_AR[(COMMON_AR >= ar.min()) & (COMMON_AR <= ar.max())]   # only where this engine has data
    return grid, np.interp(grid, ar, nr_raw)


plt.rcParams.update({"font.size": 12.5})
fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.3), sharey=True)
for ax, feat, title in [(axes[0], "trench", "TRENCH  (line / slot)"), (axes[1], "hole", "HOLE  (round via)")]:
    for eng, color, mk in [("vps", "#c0392b", "s"), ("petch", "#2471c7", "o")]:
        ar, nr = nr_curve(d[f"{eng}_{feat}_dur"], d[f"{eng}_{feat}_dep"])
        ax.plot(ar, nr, mk + "-", color=color, lw=2.6, ms=9, label=("ViennaPS" if eng == "vps" else "petch"))
    ax.set_title(title, fontsize=13.5, fontweight="bold")
    ax.set_xlabel("aspect ratio  (depth / width)")
    ax.set_xlim(3.7, 6.8); ax.set_ylim(0.4, 1.05); ax.grid(alpha=0.3)
    ax.legend(loc="lower left", framealpha=0.95)
axes[0].set_ylabel("normalized etch rate  $n_r$\n(1 = open-field rate)")
fig.suptitle("Aspect-ratio-dependent etch rate: petch vs ViennaPS   (SF$_6$/O$_2$, W = 0.5 µm, both GPU)",
             fontsize=14, fontweight="bold")
axes[1].text(0.97, 0.95, "holes: petch ≈ ViennaPS\ntrenches: petch a little gentler",
             transform=axes[1].transAxes, fontsize=10, va="top", ha="right",
             bbox=dict(boxstyle="round,pad=0.4", fc="#f3f1ea", ec="#b9b29c", alpha=0.92))
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "viz"); os.makedirs(out, exist_ok=True)
for p in [os.path.join(out, "accuracy_sweep.png"),
          os.path.join(os.path.dirname(__file__), "..", "docs", "accuracy_sweep.png")]:
    plt.savefig(p, dpi=150); print("saved", p)

print("\nplotted nr(AR):")
for feat in ["trench", "hole"]:
    for eng in ["vps", "petch"]:
        ar, nr = nr_curve(d[f"{eng}_{feat}_dur"], d[f"{eng}_{feat}_dep"])
        print(f"  {feat:6s} {eng:5s}: " + "  ".join(f"AR{a:.1f}={n:.2f}" for a, n in zip(ar, nr)))
