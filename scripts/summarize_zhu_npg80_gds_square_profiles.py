#!/usr/bin/env python3
"""Reduce the exact Oxford GDS trajectory board to partner-readable envelopes."""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


INPUT = (
    ROOT / "results" / "curated" / "zhu_npg80_gds_square_profiles_v1"
    / "audit.json"
)
SUMMARY = INPUT.parent / "profile_envelopes.json"
TABLE = INPUT.parent / "profile_envelopes.csv"
FIGURE = INPUT.parent / "profile_envelopes.png"


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _range(profiles, section, field):
    values = [float(profile[section][field]) for profile in profiles]
    return [float(min(values)), float(max(values))]


def build(audit):
    if audit.get("schema") != "petch.zhu-npg80-exact-gds-square-profile-board.v1":
        raise ValueError("not an exact Oxford GDS square profile board")
    if audit.get("smoke_only"):
        raise ValueError("refusing to summarize a smoke-only profile board")
    if audit.get("target_sem_used") or audit.get("target_depth_used"):
        raise ValueError("blind-result reducer requires a target-free board")
    grouped = {}
    for profile in audit["profiles"]:
        grouped.setdefault(float(profile["width_nm"]), []).append(profile)
    expected = set(map(float, audit["geometry"]["square_width_nm"]))
    if set(grouped) != expected:
        raise ValueError("profile widths do not match the checksum-bound GDS board")
    rows = []
    for width in sorted(grouped):
        profiles = grouped[width]
        rows.append({
            "width_nm": width,
            "pitch_nm": float(audit["geometry"]["pitch_nm"]),
            "profile_count": len(profiles),
            "etch_depth_nm": _range(profiles, "profile", "etched_depth_nm"),
            "top_cd_nm": _range(profiles, "profile", "top_cd_nm"),
            "middle_cd_nm": _range(profiles, "profile", "middle_cd_nm"),
            "bottom_cd_nm": _range(profiles, "profile", "bottom_cd_nm"),
            "sidewall_angle_from_wafer_deg": _range(
                profiles, "profile", "sidewall_angle_from_wafer_deg"),
            "bow_nm": _range(profiles, "profile", "bow_nm"),
            "cr_center_remaining_nm": _range(
                profiles, "cr_mask", "center_remaining_thickness_nm"),
            "cr_center_exhausted_fraction": float(np.mean([
                bool(profile["cr_mask"]["mask_exhausted_at_center"])
                for profile in profiles
            ])),
            "cr_layer_retired_fraction": float(np.mean([
                bool(profile["cr_mask"]["material_layer_retired"])
                for profile in profiles
            ])),
            "terminal_reasons": sorted(set(
                str(profile["terminal_reason"]) for profile in profiles)),
            "maximum_transport_relative_particle_balance_error": max(
                float(profile[
                    "maximum_transport_relative_particle_balance_error"])
                for profile in profiles),
            "maximum_state_remap_relative_conservation_residual": max(
                float(profile[
                    "maximum_state_remap_relative_conservation_residual"])
                for profile in profiles),
        })
    return {
        "schema": "petch.zhu-npg80-gds-square-profile-envelopes.v1",
        "condition_id": audit["condition_id"],
        "source_profile_board_sha256": None,
        "gds_sha256": audit["geometry"]["gds_sha256"],
        "target_sem_used": False,
        "target_depth_used": False,
        "mask_polarity_assumption": audit["geometry"][
            "mask_polarity_assumption"],
        "mask_polarity_confirmed_by_operator": audit["geometry"][
            "mask_polarity_confirmed_by_operator"],
        "claim_boundary": audit["claim_boundary"],
        "rows": rows,
    }


def _json_text(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _csv_text(summary):
    fields = [
        "width_nm", "pitch_nm", "profile_count",
        "etch_depth_min_nm", "etch_depth_max_nm",
        "top_cd_min_nm", "top_cd_max_nm",
        "middle_cd_min_nm", "middle_cd_max_nm",
        "bottom_cd_min_nm", "bottom_cd_max_nm",
        "sidewall_angle_min_deg", "sidewall_angle_max_deg",
        "bow_min_nm", "bow_max_nm",
        "cr_center_remaining_min_nm", "cr_center_remaining_max_nm",
        "cr_center_exhausted_fraction", "cr_layer_retired_fraction",
        "terminal_reasons",
    ]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in summary["rows"]:
        writer.writerow({
            "width_nm": row["width_nm"],
            "pitch_nm": row["pitch_nm"],
            "profile_count": row["profile_count"],
            "etch_depth_min_nm": row["etch_depth_nm"][0],
            "etch_depth_max_nm": row["etch_depth_nm"][1],
            "top_cd_min_nm": row["top_cd_nm"][0],
            "top_cd_max_nm": row["top_cd_nm"][1],
            "middle_cd_min_nm": row["middle_cd_nm"][0],
            "middle_cd_max_nm": row["middle_cd_nm"][1],
            "bottom_cd_min_nm": row["bottom_cd_nm"][0],
            "bottom_cd_max_nm": row["bottom_cd_nm"][1],
            "sidewall_angle_min_deg": row[
                "sidewall_angle_from_wafer_deg"][0],
            "sidewall_angle_max_deg": row[
                "sidewall_angle_from_wafer_deg"][1],
            "bow_min_nm": row["bow_nm"][0],
            "bow_max_nm": row["bow_nm"][1],
            "cr_center_remaining_min_nm": row[
                "cr_center_remaining_nm"][0],
            "cr_center_remaining_max_nm": row[
                "cr_center_remaining_nm"][1],
            "cr_center_exhausted_fraction": row[
                "cr_center_exhausted_fraction"],
            "cr_layer_retired_fraction": row["cr_layer_retired_fraction"],
            "terminal_reasons": ";".join(row["terminal_reasons"]),
        })
    return stream.getvalue()


def _plot(summary, path):
    rows = summary["rows"]
    width = np.asarray([row["width_nm"] for row in rows])
    depth = np.asarray([row["etch_depth_nm"] for row in rows])
    bottom = np.asarray([row["bottom_cd_nm"] for row in rows])
    exhausted = np.asarray([
        row["cr_center_exhausted_fraction"] for row in rows])
    figure, axes = plt.subplots(3, 1, figsize=(7.2, 8.6), sharex=True)
    axes[0].fill_between(width, depth[:, 0], depth[:, 1], alpha=0.28)
    axes[0].plot(width, depth.mean(axis=1), marker="o", linewidth=1.2)
    axes[0].set_ylabel("Etch depth (nm)")
    axes[0].grid(alpha=0.25)
    axes[1].fill_between(width, bottom[:, 0], bottom[:, 1], alpha=0.28)
    axes[1].plot(width, bottom.mean(axis=1), marker="o", linewidth=1.2)
    axes[1].plot(width, width, linestyle="--", color="black", linewidth=0.9,
                 label="nominal GDS CD")
    axes[1].set_ylabel("Bottom CD (nm)")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)
    axes[2].plot(width, exhausted, marker="o", linewidth=1.2)
    axes[2].set_ylim(-0.03, 1.03)
    axes[2].set_ylabel("Cr exhausted fraction")
    axes[2].set_xlabel("Exact GDS square width (nm)")
    axes[2].grid(alpha=0.25)
    figure.suptitle("Oxford NPG80 TiO2 exact-GDS conditional profile board")
    figure.tight_layout()
    figure.savefig(path, dpi=180, metadata={"Software": "petch"})
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    audit = json.loads(INPUT.read_text(encoding="utf-8"))
    summary = build(audit)
    summary["source_profile_board_sha256"] = _hash(INPUT)
    json_text = _json_text(summary)
    csv_text = _csv_text(summary)
    if args.write:
        SUMMARY.write_text(json_text, encoding="utf-8")
        TABLE.write_text(csv_text, encoding="utf-8")
        _plot(summary, FIGURE)
        print(SUMMARY.relative_to(ROOT))
        print(TABLE.relative_to(ROOT))
        print(FIGURE.relative_to(ROOT))
        return
    if not SUMMARY.exists() or SUMMARY.read_text(encoding="utf-8") != json_text:
        raise SystemExit("exact-GDS profile envelope JSON is stale")
    if not TABLE.exists() or TABLE.read_text(encoding="utf-8") != csv_text:
        raise SystemExit("exact-GDS profile envelope CSV is stale")
    print(f"PASS {SUMMARY.relative_to(ROOT)}")
    print(f"PASS {TABLE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
