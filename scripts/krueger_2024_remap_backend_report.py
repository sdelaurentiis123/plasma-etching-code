#!/usr/bin/env python3
"""Render the bounded Krueger remap-backend audit as a static diagnostic plot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDIT = ROOT / "results" / "krueger_2024_remap_backend_audit" / "audit.json"
DEFAULT_OUTPUT = (
    ROOT / "results" / "krueger_2024_remap_backend_audit" / "comparison.png")


def build_summary(audit, reference_backend=None, candidate_backend="common_refinement"):
    cases = dict(audit["cases"])
    if reference_backend is None:
        reference_backend = "legacy_knn" if "legacy_knn" in cases else "indexed_knn"
    reference = cases[reference_backend]["steps"][-1]["metrics"]
    candidate = cases[candidate_backend]["steps"][-1]["metrics"]
    metrics = (
        ("Depth", "etch_depth_nm"),
        ("Opening", "mask_opening_nm"),
        ("Top width", "top_feature_width_nm"),
        ("Mask thickness", "remaining_mask_thickness_nm"),
    )
    relative_ppm = []
    absolute_nm = []
    for label, name in metrics:
        reference_value = float(reference[name])
        delta = float(candidate[name]) - reference_value
        relative_ppm.append((label, 1.0e6 * delta / reference_value))
        absolute_nm.append((label, delta))
    residuals = []
    for backend in ("legacy_knn", "indexed_knn", "common_refinement"):
        if backend not in cases or cases[backend].get("status") != "complete":
            continue
        residuals.append((
            backend,
            max(float(step["maximum_remap_relative_conservation_residual"])
                for step in cases[backend]["steps"]),
        ))
    return {
        "relative_change_ppm": relative_ppm,
        "absolute_change_nm": absolute_nm,
        "maximum_conservation_residual": residuals,
        "reference_backend": reference_backend,
        "candidate_backend": candidate_backend,
        "configuration": dict(audit.get("configuration", {})),
        "partitioned_overlap_status": cases.get(
            "partitioned_overlap", {}).get("status"),
        "partitioned_overlap_message": cases.get(
            "partitioned_overlap", {}).get("exception", {}).get("message"),
    }


def render(summary, output):
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    })
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    configuration = summary["configuration"]
    figure.suptitle(
        "Krueger remap gate — "
        f"{int(configuration.get('steps', 2))} paired "
        f"{configuration.get('dx_nm', 10):g} nm steps "
        f"({configuration.get('total_physical_time_s', 0.05):g} s)",
        fontsize=13, fontweight="bold")

    labels = [item[0] for item in summary["relative_change_ppm"]]
    values = np.asarray([item[1] for item in summary["relative_change_ppm"]])
    colors = ["#2878B5" if value >= 0.0 else "#D95F02" for value in values]
    bars = axes[0].barh(labels, values, color=colors, alpha=0.88)
    axes[0].axvline(0.0, color="#333333", linewidth=0.8)
    axes[0].set_xlabel(
        f"{summary['candidate_backend']} − {summary['reference_backend']} "
        f"(ppm of {summary['reference_backend']})")
    axes[0].set_title("Second-step profile response")
    axes[0].grid(axis="x", alpha=0.25)
    for bar, value in zip(bars, values):
        if value < -200.0:
            x, horizontal = value + 75.0, "left"
        elif value >= 0.0:
            x, horizontal = value + 35.0, "left"
        else:
            x, horizontal = value - 35.0, "right"
        axes[0].text(
            x, bar.get_y() + bar.get_height() / 2.0,
            f"{value:+.3g} ppm", va="center",
            ha=horizontal)
    axes[0].invert_yaxis()

    names = [item[0].replace("_", "\n")
             for item in summary["maximum_conservation_residual"]]
    residual = [item[1] for item in summary["maximum_conservation_residual"]]
    axes[1].bar(names, residual, color="#2A9D8F", alpha=0.88)
    axes[1].set_yscale("log")
    axes[1].axhline(1.0e-12, color="#B23A48", linestyle="--", linewidth=1.2,
                    label="gate 1×10⁻¹²")
    axes[1].set_ylabel("Maximum relative remap-ledger residual")
    axes[1].set_title("All completing backends conserve")
    axes[1].set_ylim(1.0e-17, 3.0e-12)
    axes[1].grid(axis="y", which="both", alpha=0.25)
    axes[1].legend(frameon=False, loc="upper left")
    for index, value in enumerate(residual):
        axes[1].text(index, value * 1.6, f"{value:.2e}", ha="center")
    if summary["partitioned_overlap_status"] == "refused":
        axes[1].text(
            0.5, -0.23,
            "Partitioned planar overlap: certified refusal on nonparallel faces",
            transform=axes[1].transAxes, ha="center", va="top", color="#B23A48")

    relative = dict(summary["relative_change_ppm"])
    figure.text(
        0.5, 0.01,
        "Common refinement changes global depth/opening by "
        f"{relative['Depth']:+.3g}/{relative['Opening']:+.3g} ppm; "
        "the local top-width diagnostic changes by "
        f"{relative['Top width'] / 1.0e4:+.3g}%.",
        ha="center", fontsize=9)
    figure.tight_layout(rect=(0.02, 0.08, 0.98, 0.93))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference-backend")
    parser.add_argument("--candidate-backend", default="common_refinement")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    audit = json.loads(args.audit.read_text(encoding="utf-8"))
    render(build_summary(
        audit, reference_backend=args.reference_backend,
        candidate_backend=args.candidate_backend), args.output)


if __name__ == "__main__":
    main()
