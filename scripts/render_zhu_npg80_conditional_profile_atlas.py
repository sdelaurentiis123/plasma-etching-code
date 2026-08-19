#!/usr/bin/env python3
"""Render every frozen Oxford conditional profile as a deterministic SVG atlas."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "results" / "curated" / "zhu_npg80_conditional_profiles_v1"
    / "audit.json"
)
OUTPUT_DIRECTORY = SOURCE.parent
SVG_PATH = OUTPUT_DIRECTORY / "profile_atlas.svg"
MANIFEST_PATH = OUTPUT_DIRECTORY / "profile_atlas_manifest.json"

WIDTHS = (80.0, 120.0, 160.0, 200.0, 240.0, 280.0, 320.0)
SCENARIOS = (
    "ion_low_tail_0p0",
    "ion_low_tail_0p65",
    "ion_high_tail_0p0",
    "ion_high_tail_0p65",
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_source() -> dict[str, object]:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = payload["profiles"]
    if (
        payload["target_sem_used"] is not False
        or payload["target_depth_used"] is not False
        or len(profiles) != 56
    ):
        raise RuntimeError("conditional profile source is not the frozen blind board")
    return payload


def _endpoint_relief_nm(profile: dict[str, object]) -> float:
    cross_section = profile["profile"]["cross_section"]
    fraction = np.asarray(
        [row["relief_fraction_from_top"] for row in cross_section], dtype=float
    )
    height_nm = 1.0e3 * np.asarray(
        [row["height_um"] for row in cross_section], dtype=float
    )
    slope, _ = np.polyfit(fraction, height_nm, 1)
    relief = -float(slope)
    if not np.isfinite(relief) or relief <= 0.0:
        raise RuntimeError("invalid cross-section relief")
    return relief


def _ordered_profiles(payload: dict[str, object]) -> list[dict[str, object]]:
    rates = tuple(sorted({
        float(row["blanket_rate_nm_min"]) for row in payload["profiles"]
    }))
    if len(rates) != 2:
        raise RuntimeError("expected two conditional rate endpoints")
    keyed = {
        (
            float(row["width_nm"]),
            str(row["transport_scenario"]["name"]),
            float(row["blanket_rate_nm_min"]),
        ): row
        for row in payload["profiles"]
    }
    keys = [
        (width, scenario, rate)
        for scenario in SCENARIOS
        for rate in rates
        for width in WIDTHS
    ]
    if set(keyed) != set(keys):
        raise RuntimeError("conditional board axes changed")
    return [keyed[key] for key in keys]


def _svg(payload: dict[str, object]) -> str:
    profiles = _ordered_profiles(payload)
    left = 205.0
    top = 122.0
    cell_w = 132.0
    cell_h = 126.0
    plot_w = 110.0
    plot_h = 94.0
    page_w = left + len(WIDTHS) * cell_w + 18.0
    page_h = top + 8 * cell_h + 58.0
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{page_w:.0f}" '
            f'height="{page_h:.0f}" viewBox="0 0 {page_w:.0f} {page_h:.0f}">'
        ),
        '<rect width="100%" height="100%" fill="#fbfbfc"/>',
        '<style>text{font-family:Inter,Arial,sans-serif;fill:#17202a}'
        '.small{font-size:9px}.axis{font-size:10px}.label{font-size:11px}'
        '.title{font-size:19px;font-weight:700}.subtitle{font-size:11px;fill:#4d5b6a}'
        '</style>',
        '<text x="18" y="28" class="title">Oxford NPG80 blind conditional profile atlas</text>',
        '<text x="18" y="48" class="subtitle">All 56 frozen endpoints · same physical scale in every panel · no target SEM or target depth used</text>',
        '<text x="18" y="67" class="subtitle">Filled blue: 1200 s endpoint. Orange outline: last pre-clear geometry; reported depth is film-capped at 700 nm.</text>',
        '<rect x="18" y="79" width="12" height="8" fill="#5b8fd1" fill-opacity="0.55" stroke="#2f5f98"/>',
        '<text x="36" y="87" class="small">not cleared</text>',
        '<rect x="106" y="79" width="12" height="8" fill="#f2a65a" fill-opacity="0.35" stroke="#bd6512"/>',
        '<text x="124" y="87" class="small">cleared / pre-clear shape</text>',
    ]
    for column, width in enumerate(WIDTHS):
        x = left + column * cell_w + plot_w / 2.0
        lines.append(
            f'<text x="{x:.2f}" y="108" text-anchor="middle" class="label">{width:.0f} nm CD</text>'
        )
    for index, row in enumerate(profiles):
        row_index = index // len(WIDTHS)
        column = index % len(WIDTHS)
        x0 = left + column * cell_w
        y0 = top + row_index * cell_h
        metrics = row["profile"]
        cross_section = metrics["cross_section"]
        relief_nm = _endpoint_relief_nm(row)
        # All panels use the same 400 nm lateral and 700 nm vertical scale.
        x_scale = plot_w / 400.0
        y_scale = plot_h / 700.0
        samples = [
            (0.0, float(metrics["top_cd_nm"])),
            *[
                (
                    float(item["relief_fraction_from_top"]),
                    float(item["mean_width_nm"]),
                )
                for item in cross_section
            ],
            (1.0, float(metrics["bottom_cd_nm"])),
        ]
        left_edge = [
            (
                x0 + plot_w / 2.0 - 0.5 * width * x_scale,
                y0 + fraction * relief_nm * y_scale,
            )
            for fraction, width in samples
        ]
        right_edge = [
            (
                x0 + plot_w / 2.0 + 0.5 * width * x_scale,
                y0 + fraction * relief_nm * y_scale,
            )
            for fraction, width in reversed(samples)
        ]
        points = " ".join(
            f"{x:.2f},{y:.2f}" for x, y in (*left_edge, *right_edge)
        )
        cleared = bool(row["tio2_clearance_detected"])
        fill = "#f2a65a" if cleared else "#5b8fd1"
        stroke = "#bd6512" if cleared else "#2f5f98"
        lines.extend([
            f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{plot_w:.2f}" height="{plot_h:.2f}" fill="#ffffff" stroke="#d5dbe3"/>',
            f'<line x1="{x0:.2f}" y1="{y0 + plot_h:.2f}" x2="{x0 + plot_w:.2f}" y2="{y0 + plot_h:.2f}" stroke="#9aa7b4" stroke-dasharray="2 2"/>',
            f'<polygon points="{points}" fill="{fill}" fill-opacity="0.48" stroke="{stroke}" stroke-width="1.2"/>',
            (
                f'<text x="{x0 + plot_w / 2.0:.2f}" y="{y0 + plot_h + 12:.2f}" '
                f'text-anchor="middle" class="small">{relief_nm:.0f} nm'
                f'{" pre-clear" if cleared else ""}</text>'
            ),
        ])
        if column == 0:
            scenario = str(row["transport_scenario"]["name"])
            energy = "low E" if "ion_low" in scenario else "high E"
            tail = "core" if scenario.endswith("0p0") else "65% tail"
            rate = float(row["blanket_rate_nm_min"])
            lines.extend([
                f'<text x="18" y="{y0 + 37:.2f}" class="label">{energy} · {tail}</text>',
                f'<text x="18" y="{y0 + 53:.2f}" class="axis">{rate:.3f} nm/min analog</text>',
            ])
    lines.extend([
        f'<text x="18" y="{page_h - 25:.2f}" class="subtitle">Geometry is conditional on a cross-machine blanket-rate interval and inferred IADF; Cr is pinned and passivation/charging/redeposition are omitted.</text>',
        '</svg>',
    ])
    return "\n".join(lines) + "\n"


def _manifest(payload: dict[str, object], svg_text: str) -> dict[str, object]:
    profiles = _ordered_profiles(payload)
    relief = [_endpoint_relief_nm(row) for row in profiles]
    return {
        "schema": "petch.zhu-npg80-conditional-profile-atlas.v1",
        "condition_id": payload["condition_id"],
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _hash(SOURCE),
        },
        "target_sem_used": False,
        "target_depth_used": False,
        "profile_count": len(profiles),
        "cleared_profile_count": sum(
            bool(row["tio2_clearance_detected"]) for row in profiles
        ),
        "not_cleared_profile_count": sum(
            not bool(row["tio2_clearance_detected"]) for row in profiles
        ),
        "endpoint_geometry_relief_range_nm": [min(relief), max(relief)],
        "axes": {
            "width_nm": list(WIDTHS),
            "transport_scenario": list(SCENARIOS),
            "blanket_rate_nm_min": sorted({
                float(row["blanket_rate_nm_min"]) for row in profiles
            }),
        },
        "render": {
            "path": str(SVG_PATH.relative_to(ROOT)),
            "sha256": sha256(svg_text.encode("utf-8")).hexdigest(),
            "same_lateral_scale_nm": 400.0,
            "same_vertical_scale_nm": 700.0,
        },
        "claim_boundary": (
            "complete visualization of the frozen conditional board; not a "
            "validated Oxford surface law or target SEM prediction"
        ),
    }


def _render_manifest(payload: dict[str, object], svg_text: str) -> str:
    return json.dumps(_manifest(payload, svg_text), indent=2) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    payload = _load_source()
    svg_text = _svg(payload)
    manifest_text = _render_manifest(payload, svg_text)
    if args.write:
        SVG_PATH.write_text(svg_text, encoding="utf-8")
        MANIFEST_PATH.write_text(manifest_text, encoding="utf-8")
        print(SVG_PATH.relative_to(ROOT))
        return
    if (
        not SVG_PATH.exists()
        or SVG_PATH.read_text(encoding="utf-8") != svg_text
        or not MANIFEST_PATH.exists()
        or MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text
    ):
        raise SystemExit("conditional profile atlas is stale")
    print(f"PASS {SVG_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
