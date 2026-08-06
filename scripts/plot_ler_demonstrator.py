"""Figures for the LER demonstrator: Gate-2 transfer sweep and measured |T|^2."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from petch.ler_metrology import (  # noqa: E402
    EdgeStatistics, synthesize_edge_nm)
from petch.ler_transfer import estimate_transfer  # noqa: E402
from ler_gate2_shadowing import shadowed_substrate_edge_nm  # noqa: E402

DESTINATION = Path("results/curated/ler_demonstrator")


def figure_gate2():
    report = json.load(open(DESTINATION / "gate2_shadowing.json"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    for entry in report["beam_sweep"]:
        sweep = entry["sweep"]
        ax.plot([row["sigma_in_nm"] for row in sweep],
                [row["sigma_out_nm"] for row in sweep], "o-",
                label=f"beam {entry['sigma_theta_deg']:.0f}$^\\circ$")
    limit = max(row["sigma_in_nm"] for row in report["beam_sweep"][0]["sweep"])
    grid = np.linspace(0, limit, 5)
    ax.plot(grid, grid, "k--", lw=1, label="rigid transfer (slope 1)")
    ax.plot(grid, 0.5 * grid, "r:", lw=1.5, label="published slope 0.5")
    ax.set_xlabel(r"mask $\sigma_{in}$ (nm)")
    ax.set_ylabel(r"etched $\sigma_{out}$ (nm)")
    ax.set_title("Gate 2: static shadowing transfer")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for key, label, xkey in (("xi_scaling", r"$\xi$ (nm)", "correlation_length_nm"),
                             ("taper_scaling", r"$\theta_R$ (deg)",
                              "sidewall_angle_deg")):
        values = [entry[xkey] for entry in report[key]]
        ratios = [entry["fit"]["min_ratio"] for entry in report[key]]
        ax.plot(range(len(values)), ratios, "s-", label=label)
        for index, value in enumerate(values):
            ax.annotate(f"{value:g}", (index, ratios[index]),
                        textcoords="offset points", xytext=(0, 6), fontsize=7)
    ax.axhline(0.5, color="r", ls=":", lw=1.5, label="published (0.5)")
    ax.set_ylabel(r"min $\sigma_{out}/\sigma_{in}$")
    ax.set_xticks([])
    ax.set_title("structural scalings (correct sign, ~1/10 magnitude)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(DESTINATION / "gate2_transfer.png", dpi=150)
    plt.close(fig)


def figure_transfer_function():
    common = dict(roughness_exponent=0.6, sidewall_angle_deg=86.2,
                  mask_height_nm=150.0, n_points=2048, spacing_nm=1.0)
    seeds = tuple(range(101, 117))
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for xi, beam, style in ((15.0, 25.0, "-"), (30.0, 25.0, "--"),
                            (15.0, 2.0, ":")):
        statistics = EdgeStatistics(3.0, xi, common["roughness_exponent"])
        inputs, outputs = [], []
        for seed in seeds:
            edge = synthesize_edge_nm(statistics, n_points=common["n_points"],
                                      spacing_nm=common["spacing_nm"],
                                      seed=int(seed))
            inputs.append(edge)
            outputs.append(shadowed_substrate_edge_nm(
                edge, spacing_nm=common["spacing_nm"],
                mask_height_nm=common["mask_height_nm"],
                sidewall_angle_deg=common["sidewall_angle_deg"],
                sigma_theta_deg=beam))
        estimate = estimate_transfer(np.array(inputs), np.array(outputs),
                                     spacing_nm=common["spacing_nm"])
        frequency = estimate.frequencies_per_nm
        positive = frequency > 0
        ax.semilogx(frequency[positive], estimate.magnitude_squared[positive],
                    style, lw=1.4,
                    label=fr"$\xi$={xi:.0f} nm, beam {beam:.0f}$^\circ$")
    ax.axhline(1.0, color="k", lw=0.8, alpha=0.5)
    ax.set_xlabel("spatial frequency (1/nm)")
    ax.set_ylabel(r"$|T(f)|^2$")
    ax.set_ylim(0.9, 1.1)
    ax.set_title("measured shadowing transfer function (16 realizations)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(DESTINATION / "measured_transfer.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    DESTINATION.mkdir(parents=True, exist_ok=True)
    figure_gate2()
    figure_transfer_function()
    print("wrote gate2_transfer.png, measured_transfer.png")
