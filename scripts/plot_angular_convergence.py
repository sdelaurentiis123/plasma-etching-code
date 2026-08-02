"""Plots for the P0 angular-convergence harness."""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/curated/angular_convergence_p0"
payload = json.load(open(os.path.join(OUT, "angular_convergence.json")))

# --- EXP B: wall-flux profile vs beam model ------------------------------
b = payload["experiment_b"]
labels = [("narrow_collisionless", "collisionless core only, sigma=0.148 deg", "C3"),
          ("measured_two_component", "measured core+collision tail (Kim 2025)", "C0"),
          ("krueger_digitized_iead", "Krueger IEAD, production lift (ml13 path)", "C2"),
          ("krueger_iead_sqrt2_corrected", "Krueger IEAD, closure-corrected", "C1")]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
ax = axes[0]
for key, label, color in labels:
    depth = np.asarray(b[key]["profile_depth"])
    profile = np.asarray(b[key]["wall_profile"])
    ax.semilogx(np.maximum(profile, 1e-12), depth.max() - depth, label=label, color=color)
ax.set_xlabel("wall flux per band / mouth flux")
ax.set_ylabel("depth below mask top (trench widths)")
ax.invert_yaxis()
ax.set_title("AR 9 sidewall delivery vs beam model")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

ax = axes[1]
walls = [b[key]["wall_fraction"] for key, _, _ in labels]
colors = [c for _, _, c in labels]
ax.bar(range(len(walls)), walls, color=colors)
ax.set_xticks(range(len(walls)))
ax.set_xticklabels(["collisionless\ncore only", "measured\ncore+tail",
                    "Krueger IEAD\n(production)", "Krueger IEAD\n(closure-fixed)"], fontsize=7)
ax.set_ylabel("total wall flux / mouth flux")
ax.set_title("Wall-flux budget (AR 9)")
for i, value in enumerate(walls):
    ax.text(i, value, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
ax.grid(alpha=0.3, axis="y")
ax = axes[2]
c = payload.get("experiment_c", [])
if c:
    ar = [r["aspect"] for r in c]
    ax.semilogx(ar, [r["wall_ratio"] for r in c], "o-", color="C1",
                label="wall flux ratio (closure-fixed / production)")
    ax.semilogx(ar, [r["floor_ratio"] for r in c], "s-", color="C4",
                label="floor flux ratio")
    ax.axhline(2 ** 0.5, ls="--", color="k", lw=0.8, label="sqrt(2)")
    ax.axhline(1.0, ls=":", color="k", lw=0.8)
    ax.set_xlabel("aspect ratio")
    ax.set_ylabel("ratio")
    ax.set_title("Azimuthal-closure error vs AR")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "expB_mouth_wall_flux.png"), dpi=150)

# --- EXP A: convergence curves ------------------------------------------
rows = payload.get("experiment_a")
if rows:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, experiment, xlabel in zip(
            axes, ("A1_gauss_hermite", "A2_azimuth", "A3_polar_bin"),
            ("Gauss-Hermite transverse order", "azimuthal quadrature order",
             "polar bin treatment")):
        for aspect, color in zip((30.0, 100.0, 200.0), ("C0", "C1", "C2")):
            subset = [r for r in rows if r["experiment"] == experiment
                      and r["aspect"] == aspect]
            if not subset:
                continue
            x = np.arange(len(subset))
            ax.plot(x, [r["wall_fraction"] for r in subset], "o-",
                    color=color, label=f"AR {aspect:.0f}")
            ax.set_xticks(x)
            ax.set_xticklabels([r["control"].split("=")[-1] for r in subset],
                               rotation=30, fontsize=7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("wall flux / mouth flux")
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("EXP A — what angular quadrature buys (static geometry, no chemistry)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "expA_quadrature_ablation.png"), dpi=150)
print("plots written to", OUT)
