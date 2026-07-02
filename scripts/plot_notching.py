#!/usr/bin/env python3
"""Figure: the notching mechanism — deflected ions at the sidewall foot.
Left: foot-ion mean impact energy vs AR, model vs the HG tabulated points (the mechanism gate).
Right: foot etch-rate enhancement from the wired surface_charging="hg" path (the notch driver
live in petch's flux pipeline), floor throttle alongside.
Reads notching_gate_result.npz (scripts/notching_gate.py)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = np.load(os.path.join(HERE, "notching_gate_result.npz"))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.2))

axA.plot(g["ar"], g["foot_E"], "o-", color="#2471c7", lw=2.6, ms=8, label="petch 2-D charging solver")
axA.plot(g["ar"], g["hg_foot_E"], "k*", ms=16, label="Hwang–Giapis 1997 (tabulated)", zorder=6)
axA.plot(g["ar"], g["hg_foot_E"], "k--", lw=1.2, alpha=0.5)
ax2 = axA.twinx()
ax2.plot(g["ar"], g["foot_flux"], "s--", color="#e67e22", lw=1.8, ms=6, alpha=0.8)
ax2.set_ylabel("deflected-ion foot flux (normalized)", color="#e67e22")
ax2.tick_params(axis="y", labelcolor="#e67e22"); ax2.set_ylim(0, max(g["foot_flux"]) * 2.2)
axA.set_xlabel("aspect ratio"); axA.set_ylabel("foot-ion mean impact energy (eV)")
axA.set_ylim(0, 32); axA.grid(alpha=0.3); axA.legend(loc="lower right", fontsize=9)
okA = bool(g["okA"]); okB = bool(g["okB"])
axA.set_title(f"Deflected-ion energy at the foot — gate {'PASS' if okA else 'FAIL'}; "
              f"flux ~AR-independent {'PASS' if okB else 'FAIL'}", fontsize=11)

# right panel: the wired flux-path effect (from the smoke measurement, static W=0.5 trench)
AR_W = [2.0, 4.0]
floor_off = [43.40, 31.52]; floor_on = [24.75, 11.06]
foot_off = [0.347, 0.436]; foot_on = [15.276, 18.317]
x = np.arange(2); w = 0.35
axB.bar(x - w / 2, [f / o for f, o in zip(floor_on, floor_off)], w, color="#c0392b",
        label="floor rate (throttled)")
axB.bar(x + w / 2, [f / max(o, 1e-9) / 50 for f, o in zip(foot_on, foot_off)], w, color="#1e8449",
        label="foot rate enhancement (÷50 to fit)")
for i, (f, o) in enumerate(zip(foot_on, foot_off)):
    axB.text(i + w / 2, f / o / 50 + 0.02, f"×{f/o:.0f}", ha="center", fontsize=10, color="#1e8449")
for i, (f, o) in enumerate(zip(floor_on, floor_off)):
    axB.text(i - w / 2, f / o + 0.02, f"×{f/o:.2f}", ha="center", fontsize=10, color="#c0392b")
axB.set_xticks(x); axB.set_xticklabels(["AR 2", "AR 4"])
axB.set_ylabel("rate ratio (charging on / off)"); axB.set_ylim(0, 1.15)
axB.grid(alpha=0.25, axis="y"); axB.legend(loc="upper right", fontsize=9)
axB.set_title('surface_charging="hg" live in the flux path: floor throttled, foot lit', fontsize=11)

fig.suptitle("Notching mechanism — deflected-ion redistribution, gated vs Hwang–Giapis",
             fontweight="bold", fontsize=13)
plt.tight_layout()
p = os.path.join(HERE, "viz", "notching.png")
plt.savefig(p, dpi=150)
print("saved", p)
