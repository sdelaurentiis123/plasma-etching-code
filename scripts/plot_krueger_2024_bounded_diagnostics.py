#!/usr/bin/env python3
"""Plot the bounded Krueger base-case endpoint diagnostics.

This plot is deliberately restricted to three development artifacts produced
with the Krueger base boundary.  It does not load calibration targets or any
held-out oxygen/power observations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/petch-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_CAMPAIGN = Path("results/krueger_2024_r19_response_check")
SOURCE_ARTIFACTS = {
    "frozen_checkpoint_2x2": Path("frozen_checkpoint_2x2/audit.json"),
    "frozen_surface_chemistry": Path("frozen_surface_chemistry/audit.json"),
    "frozen_radiosity_chemistry": Path("frozen_radiosity_chemistry/audit.json"),
}
OUTPUT_NAME = "krueger_bounded_endpoint_diagnostics.png"
SIDECAR_NAME = "krueger_bounded_endpoint_diagnostics_provenance.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    firewall = payload.get("data_firewall", {})
    if firewall.get("boundary_case") != "base":
        raise ValueError(f"refusing non-base boundary artifact: {path}")
    if firewall.get("held_out_observations_loaded") is not False:
        raise ValueError(f"held-out firewall is not closed: {path}")
    return payload


def _mask_rate_effects(checkpoint: dict[str, Any]) -> tuple[list[str], np.ndarray]:
    metric = "velocity.amorphous_carbon_mask.net_signed_volume_rate_m3_s"
    contrasts = checkpoint["contrasts"]
    entries = [
        (
            "R17→R19 parameters\n@ R17 checkpoint",
            contrasts["instantaneous_parameter_effect_at_r17_checkpoint"][metric],
        ),
        (
            "R17→R19 parameters\n@ R19 checkpoint",
            contrasts["instantaneous_parameter_effect_at_r19_checkpoint"][metric],
        ),
        (
            "R17→R19 checkpoint\n@ R17 parameters",
            contrasts["accumulated_checkpoint_effect_at_r17_parameters"][metric],
        ),
        (
            "R17→R19 checkpoint\n@ R19 parameters",
            contrasts["accumulated_checkpoint_effect_at_r19_parameters"][metric],
        ),
    ]
    labels = [label for label, _ in entries]
    reductions = np.asarray(
        [abs(float(record["fraction_of_right"])) * 100.0 for _, record in entries],
        dtype=float,
    )
    return labels, reductions


def _passing_horizons(
    coupled: dict[str, Any],
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    passing = [item for item in coupled["horizons"] if item.get("common_pass")]
    if len(passing) != 2:
        raise ValueError(f"expected two common passing horizons, found {len(passing)}")
    horizon_ms = np.asarray([float(item["horizon_s"]) * 1e3 for item in passing])
    oxide = {
        label: np.asarray(
            [
                float(
                    item["parameter_results"][label]["nominal"]["integrated_exchange"]
                    ["removed"]["SiO2_formula_unit"]
                )
                for item in passing
            ]
        )
        for label in ("r17", "r19")
    }

    limits = coupled["gates"]
    max_by_horizon = lambda key: np.asarray(  # noqa: E731 - compact audited extraction
        [
            max(float(result[key]) for result in item["parameter_results"].values())
            for item in passing
        ]
    )
    embedded = np.asarray(
        [
            max(
                float(result["nominal"]["accepted_embedded_error_maxima"]["maximum_relative_error"])
                for result in item["parameter_results"].values()
            )
            for item in passing
        ]
    )
    gate_use = {
        "Tolerance-halving\noxide disagreement":
            100.0
            * max_by_horizon("tolerance_halving_relative_oxide_error")
            / float(limits["maximum_tolerance_halving_relative_oxide_error"]),
        "Embedded chemistry\nstep error":
            100.0 * embedded / float(limits["maximum_embedded_relative_error"]),
        "Radiosity\nbalance error":
            100.0
            * max_by_horizon("maximum_radiosity_relative_balance_error")
            / float(limits["maximum_radiosity_relative_balance_error"]),
        "Frozen-geometry\ngross displacement":
            100.0
            * max_by_horizon("maximum_gross_displacement_dx")
            / float(limits["maximum_gross_displacement_dx"]),
    }
    return horizon_ms, oxide, gate_use


def _uncoupled_probability_gate_use(uncoupled: dict[str, Any]) -> float:
    first = uncoupled["horizons"][0]
    result = first["parameter_results"]["r17"]
    drift = float(result["fine"]["neutral_reaction_probability_drift"]["maximum_absolute_probability_drift"])
    limit = float(uncoupled["gates"]["maximum_neutral_reaction_probability_absolute_drift"])
    return 100.0 * drift / limit


def _timeout_record(coupled: dict[str, Any]) -> dict[str, Any]:
    timed_out = [item for item in coupled["horizons"] if not item.get("common_pass")]
    if len(timed_out) != 1:
        raise ValueError(f"expected one unresolved horizon, found {len(timed_out)}")
    failure = timed_out[0].get("first_failure", {})
    if failure.get("classification") != "implementation_or_controller_evidence_only":
        raise ValueError("timeout is not classified as implementation/controller evidence only")
    if failure.get("physics_conclusion_permitted") is not False:
        raise ValueError("timeout artifact unexpectedly permits a physics conclusion")
    return timed_out[0]


def _annotate_bar(ax: plt.Axes, bar: Any, suffix: str = "%") -> None:
    value = float(bar.get_width())
    ax.text(
        value + 1.0,
        bar.get_y() + bar.get_height() / 2.0,
        f"{value:.2f}{suffix}",
        va="center",
        ha="left",
        fontsize=9,
        color="#172033",
    )


def build_figure(campaign: Path) -> tuple[Path, Path]:
    source_paths = {name: campaign / relative for name, relative in SOURCE_ARTIFACTS.items()}
    payloads = {name: _load_json(path) for name, path in source_paths.items()}
    source_hashes = {name: _sha256(path) for name, path in source_paths.items()}

    checkpoint = payloads["frozen_checkpoint_2x2"]
    uncoupled = payloads["frozen_surface_chemistry"]
    coupled = payloads["frozen_radiosity_chemistry"]
    effect_labels, effect_values = _mask_rate_effects(checkpoint)
    horizon_ms, oxide, gate_use = _passing_horizons(coupled)
    uncoupled_gate_use = _uncoupled_probability_gate_use(uncoupled)
    timeout = _timeout_record(coupled)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    figure = plt.figure(figsize=(14.2, 11.0), constrained_layout=True, facecolor="#f8fafc")
    grid = figure.add_gridspec(3, 2, height_ratios=(1.0, 1.12, 0.055))
    ax_effect = figure.add_subplot(grid[0, 0])
    ax_oxide = figure.add_subplot(grid[0, 1])
    ax_gate = figure.add_subplot(grid[1, :])
    ax_footer = figure.add_subplot(grid[2, :])
    figure.suptitle(
        "Krüger base-case endpoint diagnostics: bounded evidence before another long run",
        fontsize=16,
        fontweight="bold",
        color="#172033",
    )

    # Panel A: the accumulated checkpoint state matters much more than the small
    # R17→R19 parameter update for the instantaneous mask net rate.
    colors = ["#2563eb", "#2563eb", "#f59e0b", "#f59e0b"]
    y = np.arange(len(effect_values))
    bars = ax_effect.barh(y, effect_values, color=colors, height=0.62)
    ax_effect.set_yticks(y, effect_labels)
    ax_effect.invert_yaxis()
    ax_effect.set_xlim(0.0, max(effect_values) * 1.22)
    ax_effect.set_xlabel("Reduction in |mask net volume rate| (%)")
    ax_effect.set_title("A. Geometry + accumulated surface state dominate the endpoint rate", loc="left")
    ax_effect.grid(axis="x", color="#d8dee9", linewidth=0.8, alpha=0.8)
    ax_effect.set_axisbelow(True)
    for bar in bars:
        _annotate_bar(ax_effect, bar)
    ax_effect.text(
        0.98,
        0.50,
        "Checkpoint effect is 3.8–7.4× larger\nthan the parameter update.",
        transform=ax_effect.transAxes,
        ha="right",
        va="center",
        fontsize=9,
        color="#172033",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#fff7ed", "edgecolor": "#fdba74"},
    )
    ax_effect.spines[["top", "right", "left"]].set_visible(False)

    # Panel B: paired oxide-removal direction under the coupled cached-radiosity
    # chemistry integrator.
    x = np.arange(len(horizon_ms))
    width = 0.34
    r17_bars = ax_oxide.bar(x - width / 2, oxide["r17"], width, label="R17 parameters", color="#2563eb")
    r19_bars = ax_oxide.bar(x + width / 2, oxide["r19"], width, label="R19 parameters", color="#14b8a6")
    ax_oxide.set_xticks(x, [f"{value:.2f} ms" for value in horizon_ms])
    ax_oxide.set_ylabel("Integrated SiO₂ removal (formula units)")
    ax_oxide.set_title("B. R19 gives slightly less oxide removal on the same R19 checkpoint", loc="left")
    ax_oxide.grid(axis="y", color="#d8dee9", linewidth=0.8, alpha=0.8)
    ax_oxide.set_axisbelow(True)
    ax_oxide.legend(frameon=False, loc="upper left")
    for index, (left_bar, right_bar) in enumerate(zip(r17_bars, r19_bars)):
        ratio = oxide["r19"][index] / oxide["r17"][index]
        ymax = max(left_bar.get_height(), right_bar.get_height())
        ax_oxide.text(
            index,
            ymax * 1.035,
            f"R19/R17 = {ratio * 100:.2f}%\n(R19 lower by {(1.0 - ratio) * 100:.2f}%)",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#172033",
        )
    ax_oxide.set_ylim(0.0, max(oxide["r17"].max(), oxide["r19"].max()) * 1.23)
    ax_oxide.spines[["top", "right"]].set_visible(False)

    # Panel C: show every diagnostic as utilization of its declared ceiling.
    # A log axis makes both the sub-percent refinement disagreement and the
    # intentionally failed one-shot probability drift readable.
    gate_labels = list(gate_use)
    y_gate = np.arange(len(gate_labels))
    for horizon_index, (marker, color, offset) in enumerate(
        [("o", "#2563eb", -0.10), ("s", "#14b8a6", 0.10)]
    ):
        values = np.asarray([gate_use[label][horizon_index] for label in gate_labels])
        ax_gate.scatter(
            values,
            y_gate + offset,
            s=70,
            marker=marker,
            color=color,
            edgecolor="#ffffff",
            linewidth=0.8,
            label=f"Coupled subcycling: {horizon_ms[horizon_index]:.2f} ms",
            zorder=3,
        )
        for value, ypos in zip(values, y_gate + offset):
            ax_gate.text(value * 1.14, ypos, f"{value:.2f}%", va="center", fontsize=8, color="#172033")

    failure_y = len(gate_labels)
    ax_gate.scatter(
        [uncoupled_gate_use],
        [failure_y],
        s=95,
        marker="X",
        color="#dc2626",
        label="One-shot frozen flux: 7.80 ms",
        zorder=4,
    )
    ax_gate.text(
        uncoupled_gate_use / 1.12,
        failure_y,
        f"{uncoupled_gate_use:.1f}%  FAIL",
        va="center",
        ha="right",
        fontsize=9,
        color="#991b1b",
        fontweight="bold",
    )
    all_labels = gate_labels + ["One-shot frozen-flux\nlocal reaction-p drift"]
    ax_gate.set_yticks(np.arange(len(all_labels)), all_labels)
    ax_gate.invert_yaxis()
    ax_gate.set_xscale("log")
    ax_gate.set_xlim(0.35, 900.0)
    ax_gate.axvline(100.0, color="#dc2626", linewidth=1.4, linestyle="--")
    ax_gate.text(100.0, -0.55, "declared limit", ha="center", va="bottom", fontsize=9, color="#991b1b")
    ax_gate.set_xlabel("Share of declared numerical/scope limit (%) — left of 100% passes")
    ax_gate.set_title(
        "C. Cached-radiosity chemistry subcycling passes through 15.59 ms",
        loc="left",
    )
    ax_gate.grid(axis="x", which="both", color="#d8dee9", linewidth=0.8, alpha=0.8)
    ax_gate.set_axisbelow(True)
    ax_gate.legend(frameon=False, loc="upper right", ncol=3, fontsize=9)
    ax_gate.spines[["top", "right", "left"]].set_visible(False)

    timeout_ms = float(timeout["horizon_s"]) * 1e3
    timeout_reason = timeout["first_failure"]["reason"]
    ax_gate.text(
        0.995,
        0.24,
        f"{timeout_ms:.2f} ms: bounded timeout\nIMPLEMENTATION / CONTROLLER EVIDENCE ONLY\nNo physics failure; no physics conclusion.\n{timeout_reason}",
        transform=ax_gate.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#334155",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f1f5f9", "edgecolor": "#94a3b8"},
    )

    ax_footer.axis("off")
    ax_footer.text(
        0.5,
        0.5,
        "Development evidence only • Krüger base boundary • held-out observations were not loaded • no profile evolution was run",
        transform=ax_footer.transAxes,
        ha="center",
        va="center",
        fontsize=9,
        color="#475569",
    )

    output = campaign / OUTPUT_NAME
    sidecar = campaign / SIDECAR_NAME
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "Title": "Krueger bounded endpoint diagnostics",
        "Description": "Base-case development evidence; held-out observations not loaded.",
        "Source artifact SHA-256": json.dumps(source_hashes, sort_keys=True),
    }
    figure.savefig(output, dpi=180, bbox_inches="tight", facecolor=figure.get_facecolor(), metadata=metadata)
    plt.close(figure)

    derived = {
        "mask_rate_reduction_percent": dict(zip(effect_labels, effect_values.tolist())),
        "oxide_removal_formula_units": {
            f"{horizon_ms[index]:.8f}_ms": {
                "r17": float(oxide["r17"][index]),
                "r19": float(oxide["r19"][index]),
                "r19_to_r17_ratio": float(oxide["r19"][index] / oxide["r17"][index]),
            }
            for index in range(len(horizon_ms))
        },
        "gate_use_percent_of_limit": {
            label: {
                f"{horizon_ms[index]:.8f}_ms": float(values[index])
                for index in range(len(horizon_ms))
            }
            for label, values in gate_use.items()
        },
        "one_shot_frozen_flux_probability_drift_percent_of_limit": uncoupled_gate_use,
        "unresolved_horizon": {
            "horizon_ms": timeout_ms,
            "classification": timeout["first_failure"]["classification"],
            "physics_conclusion_permitted": timeout["first_failure"]["physics_conclusion_permitted"],
            "reason": timeout_reason,
        },
    }
    sidecar_payload = {
        "schema": "petch.krueger_2024.bounded_endpoint_diagnostics_plot.v1",
        "claim_scope": "base-case development evidence; not held-out validation",
        "held_out_observations_loaded": False,
        "output_png": output.name,
        "source_artifacts": {
            name: {"path": str(path), "sha256": source_hashes[name]}
            for name, path in source_paths.items()
        },
        "derived_values": derived,
    }
    sidecar.write_text(json.dumps(sidecar_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output, sidecar


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campaign",
        type=Path,
        default=DEFAULT_CAMPAIGN,
        help=f"campaign directory (default: {DEFAULT_CAMPAIGN})",
    )
    args = parser.parse_args()
    output, sidecar = build_figure(args.campaign)
    print(output)
    print(sidecar)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
