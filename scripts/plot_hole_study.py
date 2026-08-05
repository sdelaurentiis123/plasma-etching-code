"""Figures for the HAR hole study phase-1 results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def plot_arde(payload, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, sweep in sorted(payload["series4_arde"].items()):
        fraction = float(key.split("_")[1])
        aspect = [row["aspect_ratio"] for row in sweep["rows"]]
        total = [row["total_bottom"] for row in sweep["rows"]]
        cascaded = [row["cascaded_share_of_bottom"] for row in sweep["rows"]]
        label = f"tail {fraction:.2f}"
        axes[0].plot(aspect, total, marker="o", ms=3, label=label)
        axes[1].plot(aspect, cascaded, marker="o", ms=3, label=label)
    axes[0].set(xlabel="aspect ratio", ylabel="energetic delivery to etch front",
                title="ARDE of total energetic delivery", xscale="log")
    axes[0].set_ylim(0.0, 1.05)
    axes[1].set(xlabel="aspect ratio",
                ylabel="cascaded share of delivered flux",
                title="hot-neutral share at the etch front", xscale="log")
    axes[1].set_ylim(0.0, 1.05)
    for axis in axes:
        axis.grid(alpha=0.3)
        axis.legend(fontsize=8, title="tail fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_wall_profile(payload, path):
    records = [r for r in payload["series3_cascade"]
               if r["aspect_ratio"] == max(x["aspect_ratio"]
                                           for x in payload["series3_cascade"])
               and r["max_bounces"] == 8]
    fig, axis = plt.subplots(figsize=(6.4, 4.2))
    for record in sorted(records, key=lambda r: r["tail_fraction"]):
        centres = np.asarray(record["wall_height_bin_centres"])
        profile = np.asarray(record["wall_rate_profile"])
        depth_from_top = record["aspect_ratio"] - centres
        width = centres[1] - centres[0] if centres.size > 1 else 1.0
        axis.plot(depth_from_top, profile / width, lw=1.4,
                  label=f"tail {record['tail_fraction']:.2f}")
    axis.set(xlabel="depth below entrance (hole diameters)",
             ylabel="wall energy deposition per unit depth",
             title=f"sidewall deposition, AR {records[0]['aspect_ratio']:.0f}",
             yscale="log")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8, title="tail fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_bottom_split(payload, path):
    deep = [r for r in payload["series3_cascade"] if r["max_bounces"] == 8]
    aspects = sorted({r["aspect_ratio"] for r in deep})
    fractions = sorted({r["tail_fraction"] for r in deep})
    fig, axis = plt.subplots(figsize=(6.8, 4.2))
    width = 0.8 / len(fractions)
    index = np.arange(len(aspects))
    for offset, fraction in enumerate(fractions):
        direct, cascaded = [], []
        for aspect in aspects:
            row = next(r for r in deep if r["aspect_ratio"] == aspect
                       and r["tail_fraction"] == fraction)
            direct.append(row["direct_bottom"])
            cascaded.append(row["cascaded_bottom"])
        position = index + offset * width - 0.4 + width / 2
        axis.bar(position, direct, width * 0.92, color=f"C{offset}", alpha=0.95,
                 label=f"tail {fraction:.2f} direct")
        axis.bar(position, cascaded, width * 0.92, bottom=direct,
                 color=f"C{offset}", alpha=0.45,
                 label=f"tail {fraction:.2f} cascaded")
    axis.set_xticks(index)
    axis.set_xticklabels([f"{a:.0f}" for a in aspects])
    axis.set(xlabel="aspect ratio", ylabel="delivered fraction of entering flux",
             title="etch-front delivery: direct ions vs cascaded hot particles")
    axis.grid(alpha=0.3, axis="y")
    axis.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/curated/hole_study/phase1.json")
    args = parser.parse_args()
    payload = json.loads(Path(args.input).read_text())
    out = Path(args.input).parent
    plot_arde(payload, out / "arde_delivery.png")
    plot_wall_profile(payload, out / "wall_deposition.png")
    plot_bottom_split(payload, out / "bottom_delivery_split.png")
    print(f"wrote figures to {out}")


if __name__ == "__main__":
    main()
