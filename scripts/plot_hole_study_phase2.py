"""Figures for hole-study phase 2 (coupled rate-ARDE and the straight-wall envelope)."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / "results/curated/hole_study/phase2.json").read_text())
out = ROOT / "results/curated/hole_study"

rows = data["series_a_rate_arde"]
tails = sorted({r["tail_fraction"] for r in rows})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for tail in tails:
    sub = sorted((r for r in rows if r["tail_fraction"] == tail),
                 key=lambda r: r["aspect_ratio"])
    ar = [r["aspect_ratio"] for r in sub]
    rate = [r["floor_etch_nm_s"] for r in sub]
    axes[0].plot(ar, rate, marker="o", ms=3, label=f"tail {tail:.2f}")
    axes[1].plot(ar, [v / rate[0] for v in rate], marker="o", ms=3,
                 label=f"tail {tail:.2f}")
    if tail == max(tails):
        axes[1].plot(ar, [r["energetic_delivery"] for r in sub], ls="--", color="k",
                     lw=1, label="energetic delivery (tail 0.65)")
axes[0].set_xscale("log"); axes[0].set_xlabel("aspect ratio")
axes[0].set_ylabel("floor etch rate (nm/s)")
axes[0].set_title("Coupled rate vs aspect ratio")
axes[1].set_xscale("log"); axes[1].set_xlabel("aspect ratio")
axes[1].set_ylabel("rate / rate(AR 1)")
axes[1].set_title("Rate ARDE vs transport-only delivery")
for ax in axes:
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(out / "phase2_rate_arde.png", dpi=150)

fig, ax = plt.subplots(figsize=(5.5, 4.0))
for row in data["series_c_envelope"]:
    trace = row["trace"]
    ax.plot([t["time_s"] for t in trace], [t["straightness_deviation"] for t in trace],
            marker="o", ms=3, label=f"AR {row['aspect_ratio_initial']:.0f}")
ax.axhline(data["series_c_envelope"][0]["straightness_tolerance"], color="r", ls="--",
           lw=1, label="declared envelope")
ax.set_xlabel("process time (s)"); ax.set_ylabel("straightness deviation")
ax.set_title("Straight-wall validity horizon"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(out / "phase2_straight_wall_envelope.png", dpi=150)
print("wrote figures")
