#!/usr/bin/env python3
"""Plot the full accuracy sweep: petch vs ViennaPS, trench (left) + hole (right), nr vs AR.
Reads accuracy_sweep.npz (pass path as argv[1]), writes accuracy_sweep.png."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W = 0.5
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "accuracy_sweep.npz")


AR_REF = 3.0          # normalize both engines to nr=1 at this aspect ratio (apples-to-apples)
SUB = 6.3             # domain top; drop samples that floored on the substrate bottom


def nr_curve(dur, dep):
    keep = dep < SUB - 0.5                         # trim domain-floored points (rate->0 artifact)
    dur, dep = dur[keep], dep[keep]
    ar = 0.5 * (dep[1:] + dep[:-1]) / W
    rate = np.diff(dep) / np.diff(dur)
    ref = np.interp(AR_REF, ar, rate) if ar.min() <= AR_REF <= ar.max() else rate[0]
    return ar, rate / ref


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, feat, title in [(axes[0], "trench", "TRENCH (periodic)"), (axes[1], "hole", "HOLE (3D)")]:
    for eng, color, mk in [("vps", "tab:red", "s"), ("petch", "tab:blue", "o")]:
        kd, kp = f"{eng}_{feat}_dur", f"{eng}_{feat}_dep"
        if kd in d and kp in d:
            ar, nr = nr_curve(d[kd], d[kp])
            label = "ViennaPS" if eng == "vps" else "petch"
            ax.plot(ar, nr, mk + "-", color=color, lw=2.2, ms=7, label=label)
    ax.set_title(title); ax.set_xlabel("aspect ratio (depth / width)")
    ax.set_xlim(0, None); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend()
axes[0].set_ylabel("normalized etch rate  nr")
fig.suptitle("ARDE accuracy: petch vs ViennaPS  —  W=0.5 µm  (SF6/O2, dx=0.04)", fontsize=13)
plt.tight_layout()
plt.savefig("accuracy_sweep.png", dpi=140)
print("saved accuracy_sweep.png")

# print the numbers at common AR points
print("\nnr at AR = 3, 5, 7:")
for feat in ["trench", "hole"]:
    for eng in ["vps", "petch"]:
        kd, kp = f"{eng}_{feat}_dur", f"{eng}_{feat}_dep"
        if kd in d:
            ar, nr = nr_curve(d[kd], d[kp])
            vals = [f"{np.interp(a, ar, nr):.2f}" if ar.max() >= a else "  - " for a in (3, 5, 7)]
            print(f"  {feat:6s} {eng:5s}: AR3={vals[0]} AR5={vals[1]} AR7={vals[2]}  (max AR {ar.max():.1f})")
