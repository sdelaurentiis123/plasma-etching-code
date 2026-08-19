#!/usr/bin/env python3
"""Reproduce the Hong et al. 2023 TiO2 feature-response board.

The source pixels are intentionally not redistributed.  Figure 2 marker
centres were inspected on a 400-dpi Poppler render at original resolution;
Figure 3 values are author-printed annotations rather than scale-bar
measurements.  This board constrains feature transport and mask survival.  It
is not an absolute-chemistry transfer to the Zhu Oxford NPG80 process.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "hong_2023_tio2"
FIGURE2_CSV_PATH = OUTPUT_DIRECTORY / "figure2_feature_response.csv"
FIGURE3_CSV_PATH = OUTPUT_DIRECTORY / "figure3_reported_profile_labels.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "165aa9bb4bccb1212bba1f772bf273a7d13f4300f5e5e78a36fce9c16f73bc45"
)
FIGURE2_PAGE_SHA256 = (
    "ae7f30fe47476ca852c851f6632715e73ea70db4c2d9f12ab0deb261826ed1ef"
)
FIGURE2_CROP_SHA256 = (
    "08b48ce2f9c834ae2608c1931bd5d1bf3f0a97eebbd6d73ae7e1f4dd04519b97"
)
FIGURE3_PAGE_SHA256 = (
    "557e4896c18594b6f728a4431b7f105a5a0b1ea574ab235559ac60b2ab73d23a"
)
FIGURE3_CROP_SHA256 = (
    "95697f27f7ab6f395133069cac438156019c94b5d81a835d5e9dd569fd2e5d95"
)

FIGURE2_CROP_ORIGIN_FULL_PAGE_PX = (150.0, 200.0)
FIGURE2_CROP_SIZE_PX = (3000, 2300)
RENDER_SIZE_PX = (3308, 4410)
PLASMA_MODES = ("CW", "S30_B70", "S70_B30")


@dataclass(frozen=True)
class Axis:
    panel: str
    quantity: str
    units: str
    y_high_px: float
    value_high: float
    y_low_px: float
    value_low: float

    def value_at(self, y_px: float) -> float:
        fraction = (float(y_px) - self.y_high_px) / (
            self.y_low_px - self.y_high_px
        )
        return self.value_high + fraction * (
            self.value_low - self.value_high
        )

    @property
    def value_per_pixel(self) -> float:
        return (self.value_low - self.value_high) / (
            self.y_low_px - self.y_high_px
        )


AXES = {
    "a": Axis("a", "tio2_etch_rate", "nm/min", 124.5, 120.0, 859.5, 0.0),
    "b": Axis("b", "acl_etch_rate", "nm/min", 124.5, 80.0, 859.5, 0.0),
    "c": Axis("c", "tio2_acl_selectivity", "dimensionless", 1195.5, 15.0, 1932.0, 0.0),
    "d": Axis("d", "arde_p1_depth_over_p3_depth", "dimensionless", 1195.5, 2.5, 1932.0, 1.0),
}


@dataclass(frozen=True)
class PixelSeries:
    panel: str
    chemistry: str
    specimen: str
    marker: str
    x_px: tuple[float, float, float]
    y_px: tuple[float, float, float]
    error_high_y_px: tuple[float, float, float] | None = None
    error_low_y_px: tuple[float, float, float] | None = None


# Coordinates are relative to the declared Figure 2 crop.  The marker centres
# are the original-resolution visual picks; the ±5 pixel bound includes both
# centre placement and thick-axis placement.  Figure 2(d) error-bar cap centres
# are retained separately from this digitization allowance.
PIXEL_SERIES = (
    PixelSeries("a", "C4F8/SF6/Ar", "blank", "filled_square", (653, 960, 1266), (240, 557, 600)),
    PixelSeries("a", "C4F8/SF6/Ar", "P1", "filled_circle", (653, 960, 1266), (271, 586, 632)),
    PixelSeries("a", "BCl3/CF4/Ar", "blank", "filled_up_triangle", (653, 960, 1266), (490, 719, 746)),
    PixelSeries("a", "BCl3/CF4/Ar", "P1", "filled_down_triangle", (653, 960, 1266), (518, 760, 782)),
    PixelSeries("b", "C4F8/SF6/Ar", "blank", "filled_square", (1957, 2262, 2568), (185, 722, 762)),
    PixelSeries("b", "C4F8/SF6/Ar", "P1", "filled_circle", (1957, 2262, 2568), (236, 757, 789)),
    PixelSeries("b", "BCl3/CF4/Ar", "blank", "filled_up_triangle", (1957, 2262, 2568), (689, 796, 823)),
    PixelSeries("b", "BCl3/CF4/Ar", "P1", "filled_down_triangle", (1957, 2262, 2568), (720, 809, 818)),
    PixelSeries("c", "C4F8/SF6/Ar", "blank", "filled_square", (653, 960, 1266), (1847, 1682, 1594)),
    PixelSeries("c", "C4F8/SF6/Ar", "P1", "filled_circle", (653, 960, 1266), (1857, 1737, 1646)),
    PixelSeries("c", "BCl3/CF4/Ar", "blank", "filled_up_triangle", (653, 960, 1266), (1681, 1564, 1400)),
    PixelSeries("c", "BCl3/CF4/Ar", "P1", "filled_down_triangle", (653, 960, 1266), (1682, 1584, 1433)),
    PixelSeries(
        "d", "C4F8/SF6/Ar", "P1_over_P3", "filled_square",
        (1957, 2262, 2568), (1694, 1876, 1888),
        (1662, 1854, 1863.5), (1711.5, 1903, 1912.5),
    ),
    PixelSeries(
        "d", "BCl3/CF4/Ar", "P1_over_P3", "filled_circle",
        (1957, 2262, 2568), (1493, 1563, 1612),
        (1466, 1539.5, 1588.5), (1515, 1588.5, 1637.5),
    ),
)

FIGURE2_FIELDS = (
    "panel", "quantity", "chemistry", "specimen", "plasma_mode",
    "source_power_on_percent", "bias_power_on_percent", "value", "units",
    "experimental_error_low", "experimental_error_high",
    "marker_center_x_crop_px", "marker_center_y_crop_px",
    "marker_center_x_full_page_px", "marker_center_y_full_page_px", "marker",
    "digitization_absolute_bound",
)


def _duty(mode: str) -> tuple[int, int]:
    return {"CW": (100, 100), "S30_B70": (30, 70), "S70_B30": (70, 30)}[mode]


def figure2_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    x_offset, y_offset = FIGURE2_CROP_ORIGIN_FULL_PAGE_PX
    for series in PIXEL_SERIES:
        axis = AXES[series.panel]
        for index, mode in enumerate(PLASMA_MODES):
            x_px = series.x_px[index]
            y_px = series.y_px[index]
            value = axis.value_at(y_px)
            source_on, bias_on = _duty(mode)
            if series.error_high_y_px is None:
                error_low = error_high = ""
            else:
                error_high = f"{axis.value_at(series.error_high_y_px[index]):.9g}"
                error_low = f"{axis.value_at(series.error_low_y_px[index]):.9g}"
            rows.append(
                {
                    "panel": series.panel,
                    "quantity": axis.quantity,
                    "chemistry": series.chemistry,
                    "specimen": series.specimen,
                    "plasma_mode": mode,
                    "source_power_on_percent": str(source_on),
                    "bias_power_on_percent": str(bias_on),
                    "value": f"{value:.9g}",
                    "units": axis.units,
                    "experimental_error_low": error_low,
                    "experimental_error_high": error_high,
                    "marker_center_x_crop_px": f"{x_px:.1f}",
                    "marker_center_y_crop_px": f"{y_px:.1f}",
                    "marker_center_x_full_page_px": f"{x_px + x_offset:.1f}",
                    "marker_center_y_full_page_px": f"{y_px + y_offset:.1f}",
                    "marker": series.marker,
                    "digitization_absolute_bound": f"{5.0 * abs(axis.value_per_pixel):.9g}",
                }
            )
    return rows


def _csv_text(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def figure2_csv_text() -> str:
    return _csv_text(FIGURE2_FIELDS, figure2_rows())


FIGURE3_FIELDS = (
    "chemistry", "plasma_mode", "source_power_on_percent",
    "bias_power_on_percent", "reported_remaining_acl_nm",
    "reported_p1_tio2_depth_nm", "evidence_type",
)
FIGURE3_VALUES = {
    "C4F8/SF6/Ar": ((116, 392), (395, 450), (424, 454)),
    "BCl3/CF4/Ar": ((426, 370), (443, 340), (465, 374)),
}


def figure3_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for chemistry, values in FIGURE3_VALUES.items():
        for mode, (acl_nm, depth_nm) in zip(PLASMA_MODES, values):
            source_on, bias_on = _duty(mode)
            rows.append(
                {
                    "chemistry": chemistry,
                    "plasma_mode": mode,
                    "source_power_on_percent": str(source_on),
                    "bias_power_on_percent": str(bias_on),
                    "reported_remaining_acl_nm": str(acl_nm),
                    "reported_p1_tio2_depth_nm": str(depth_nm),
                    "evidence_type": "author_annotated_sem",
                }
            )
    return rows


def figure3_csv_text() -> str:
    return _csv_text(FIGURE3_FIELDS, figure3_rows())


def manifest(figure2_digest: str, figure3_digest: str) -> dict[str, object]:
    arde = {
        (row["chemistry"], row["plasma_mode"]): float(row["value"])
        for row in figure2_rows()
        if row["panel"] == "d"
    }
    c4_selectivity = {
        row["plasma_mode"]: float(row["value"])
        for row in figure2_rows()
        if row["panel"] == "c"
        and row["chemistry"] == "C4F8/SF6/Ar"
        and row["specimen"] == "P1"
    }
    return {
        "manifest_id": "HONG-2023-TIO2-FEATURE-RESPONSE-R1",
        "source": {
            "citation": (
                "J. W. Hong et al., Materials Science in Semiconductor "
                "Processing 164, 107617 (2023)"
            ),
            "doi": "10.1016/j.mssp.2023.107617",
            "author_copy_url": "https://swb.skku.edu/_res/pnpl/etc/2023-12.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "render_dpi": 400,
            "render_size_px": list(RENDER_SIZE_PX),
            "figure2_pdf_page": 3,
            "figure2_page_sha256": FIGURE2_PAGE_SHA256,
            "figure2_crop_geometry": "3000x2300+150+200",
            "figure2_crop_sha256": FIGURE2_CROP_SHA256,
            "figure3_pdf_page": 4,
            "figure3_page_sha256": FIGURE3_PAGE_SHA256,
            "figure3_crop_geometry": "1900x1050+550+150",
            "figure3_crop_sha256": FIGURE3_CROP_SHA256,
            "redistribution": "source PDF and pixels are not committed",
        },
        "experiment": {
            "reactor": "300 mm inductively coupled plasma etcher",
            "icp_source_power_W": 1000,
            "bias_power_W": 150,
            "pressure_mTorr": 2,
            "substrate_temperature_C": 70,
            "pulse_frequency_Hz": 1000,
            "c4f8_sf6_ar_sccm": [130, 20, 15],
            "bcl3_cf4_ar_sccm": [100, 25, 60],
            "measured_dc_bias_V": {
                "C4F8/SF6/Ar": -135,
                "BCl3/CF4/Ar": -190,
            },
            "nominal_tio2_thickness_nm": 500,
            "nominal_acl_thickness_nm": 500,
            "patterns": {
                "P1": {"diameter_nm": 80, "space_nm": 245},
                "P2": {"diameter_nm": 155, "space_nm": 160},
                "P3": {"diameter_nm": 220, "space_nm": 80},
            },
        },
        "pixel_calibration": {
            panel: {
                "quantity": axis.quantity,
                "units": axis.units,
                "crop_y_high_px": axis.y_high_px,
                "value_high": axis.value_high,
                "crop_y_low_px": axis.y_low_px,
                "value_low": axis.value_low,
                "value_per_vertical_pixel": axis.value_per_pixel,
            }
            for panel, axis in AXES.items()
        },
        "digitization": {
            "method": (
                "400-dpi Poppler render; original-resolution visual marker "
                "inspection; color-channel row-density cross-check"
            ),
            "marker_and_axis_allowance_px": 5,
            "figure2d_error_bars": (
                "cap centres digitized separately; these are experimental "
                "error bars, not the digitization allowance"
            ),
            "figure3_values": "author annotations transcribed exactly",
            "visual_audit_status": "passed_original_resolution",
        },
        "internal_checks": {
            "paper_text_c4f8_p1_selectivity_s30_b70": 3.9,
            "digitized_c4f8_p1_selectivity_s30_b70": c4_selectivity["S30_B70"],
            "paper_text_c4f8_p1_selectivity_s70_b30": 5.9,
            "digitized_c4f8_p1_selectivity_s70_b30": c4_selectivity["S70_B30"],
            "c4f8_arde_cw": arde[("C4F8/SF6/Ar", "CW")],
            "c4f8_arde_s70_b30": arde[("C4F8/SF6/Ar", "S70_B30")],
            "pulsing_moves_c4f8_arde_toward_unity": abs(
                arde[("C4F8/SF6/Ar", "S70_B30")] - 1.0
            ) < abs(arde[("C4F8/SF6/Ar", "CW")] - 1.0),
        },
        "claim_boundary": {
            "valid": [
                "external TiO2 feature-transport and ARDE response target",
                "external mask-survival/profile-form target",
                "evidence that layout-dependent feature response can occur under one reactor condition",
            ],
            "not_valid": [
                "an absolute surface-law calibration for the Zhu CHF3/SF6/O2 Oxford process",
                "a radial Oxford NPG80 flux map",
                "a chromium-mask erosion law",
                "evidence that clustered pillar collapse is caused by reactor nonuniformity",
            ],
            "transfer_differences_from_zhu": [
                "ICP versus capacitively coupled RIE",
                "C4F8/SF6/Ar versus CHF3/SF6/O2",
                "2 mTorr versus 30 mTorr",
                "ACL versus chromium mask",
                "circular holes/pillars versus square/rectangular metasurface pillars",
            ],
        },
        "outputs": {
            "figure2_csv": {
                "path": "data/experimental/hong_2023_tio2/figure2_feature_response.csv",
                "sha256": figure2_digest,
            },
            "figure3_csv": {
                "path": "data/experimental/hong_2023_tio2/figure3_reported_profile_labels.csv",
                "sha256": figure3_digest,
            },
        },
    }


def manifest_text(figure2_payload: str, figure3_payload: str) -> str:
    return json.dumps(
        manifest(
            sha256(figure2_payload.encode()).hexdigest(),
            sha256(figure3_payload.encode()).hexdigest(),
        ),
        indent=2,
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    figure2_payload = figure2_csv_text()
    figure3_payload = figure3_csv_text()
    payloads = {
        FIGURE2_CSV_PATH: figure2_payload,
        FIGURE3_CSV_PATH: figure3_payload,
        MANIFEST_PATH: manifest_text(figure2_payload, figure3_payload),
    }
    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        for path, payload in payloads.items():
            path.write_text(payload, encoding="utf-8")
    else:
        for path, payload in payloads.items():
            if not path.exists() or path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"stale or missing generated artifact: {path}")


if __name__ == "__main__":
    main()
