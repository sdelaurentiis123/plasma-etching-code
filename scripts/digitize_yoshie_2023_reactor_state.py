#!/usr/bin/env python3
"""Replay Yoshie et al. 2023 reactor-state diagnostics from Figures 12 and 14.

The two outputs are deliberately diagnostic constraints, not wafer fluxes:

* Figure 12 supplies electron density sampled at the start, midpoint, and end
  of each 0.25 s bias window.  Simpson averages preserve the fast SF6-driven
  transition in 8 s timing II.
* Figure 14 supplies phase-resolved CF, CF2, and F optical-emission ratios to
  Ar.  Those ratios constrain timing and relative shape only; no actinometry
  or excited-state kinetics was published that would turn them into absolute
  ground-state particle fluxes.

Publisher rasters are not redistributed.  This script checksum-verifies the
official images, reproduces the committed CSVs and manifest exactly, and can
draw full-resolution QA overlays.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "experimental" / "yoshie_2023"
ELECTRON_CSV = OUTPUT_DIR / "figure12_bias_window_electron_density.csv"
OES_CSV = OUTPUT_DIR / "figure14_phase_resolved_oes.csv"
MANIFEST_PATH = OUTPUT_DIR / "reactor_state_digitization_manifest.json"
DEFAULT_SOURCE_DIR = ROOT / "tmp" / "sources" / "yoshie_2023"

SOURCE = {
    12: {
        "name": "yoshie_fig12.jpg",
        "url": (
            "https://ars.els-cdn.com/content/image/"
            "1-s2.0-S0169433223016604-gr12_lrg.jpg"
        ),
        "sha256": (
            "9ee2992ef60c22debb6619e79ea97446ac0698c258e746322bd8858e70ee2109"
        ),
        "size": (2713, 1476),
    },
    14: {
        "name": "yoshie_fig14.jpg",
        "url": (
            "https://ars.els-cdn.com/content/image/"
            "1-s2.0-S0169433223016604-gr14_lrg.jpg"
        ),
        "sha256": (
            "e7a121f1b9be7543c217626978d0ae243f4e22cd87a00765877cbe65b7c547b0"
        ),
        "size": (1338, 4337),
    },
}


@dataclass(frozen=True)
class Axis:
    x_left_px: float
    x_right_px: float
    x_left_value: float
    x_right_value: float
    y_zero_px: float
    y_reference_px: float
    y_zero_value: float
    y_reference_value: float

    def x_px(self, value):
        fraction = (
            (float(value) - self.x_left_value)
            / (self.x_right_value - self.x_left_value)
        )
        return self.x_left_px + fraction * (
            self.x_right_px - self.x_left_px
        )

    def y_value(self, y_px):
        fraction = (
            (self.y_zero_px - float(y_px))
            / (self.y_zero_px - self.y_reference_px)
        )
        return self.y_zero_value + fraction * (
            self.y_reference_value - self.y_zero_value
        )


AXES = {
    "12a": Axis(225.0, 883.0, -2.0, 2.0, 1302.0, 496.0, 0.0, 1.8),
    "12b": Axis(1298.0, 2617.0, -4.0, 4.0, 1302.0, 496.0, 0.0, 2.5),
    "14b": Axis(216.0, 1259.0, -2.0, 2.0, 1790.0, 1228.0, -3.0, 17.0),
    "14e": Axis(216.0, 1259.0, -2.0, 2.0, 4197.0, 3558.0, -0.4, 1.6),
}


@dataclass(frozen=True)
class BiasWindow:
    panel: str
    cycle_duration_s: int
    timing: str
    bias_start_s: float
    sample_y_px: tuple[float, float, float]
    window_uncertainty_1e10_cm3: float

    @property
    def bias_end_s(self):
        return self.bias_start_s + 0.25

    @property
    def sample_times_s(self):
        midpoint = 0.5 * (self.bias_start_s + self.bias_end_s)
        return (self.bias_start_s, midpoint, self.bias_end_s)

    @property
    def sample_x_px(self):
        return tuple(AXES[self.panel].x_px(value)
                     for value in self.sample_times_s)

    @property
    def sample_density_1e10_cm3(self):
        return tuple(AXES[self.panel].y_value(value)
                     for value in self.sample_y_px)

    @property
    def simpson_average_1e10_cm3(self):
        start, midpoint, end = self.sample_density_1e10_cm3
        return (start + 4.0 * midpoint + end) / 6.0


# Full-resolution centerline transcriptions from the measured-marker/dashed
# traces.  The steep 8 s timing-II transition receives the largest allowance;
# the allowance covers trace interpolation and pixel placement, not the
# unreported uncertainty of the plasma diagnostic itself.
BIAS_WINDOWS = (
    BiasWindow("12a", 4, "I", -0.15, (702.0, 635.0, 870.0), 0.06),
    BiasWindow("12a", 4, "II", 0.35, (1158.0, 1179.0, 1178.0), 0.04),
    BiasWindow("12a", 4, "III", 1.35, (1077.0, 1064.0, 1044.0), 0.04),
    BiasWindow("12b", 8, "I", -0.15, (708.0, 707.0, 925.0), 0.08),
    BiasWindow("12b", 8, "II", 0.10, (780.0, 945.0, 1175.0), 0.15),
    BiasWindow("12b", 8, "III", 1.10, (1245.0, 1250.0, 1232.0), 0.04),
    BiasWindow("12b", 8, "IV", 2.10, (1120.0, 1093.0, 1065.0), 0.05),
)


@dataclass(frozen=True)
class OESPoint:
    panel: str
    signal: str
    time_s: float
    y_px: float

    @property
    def x_px(self):
        return AXES[self.panel].x_px(self.time_s)

    @property
    def ratio_to_ar(self):
        return AXES[self.panel].y_value(self.y_px)


OES_TIMES = tuple(-1.9 + 0.2 * index for index in range(20))


def _oes_series(panel, signal, y_pixels):
    return tuple(
        OESPoint(panel, signal, time, y_px)
        for time, y_px in zip(OES_TIMES, y_pixels)
    )


# Centers were isolated by full-resolution HSV color-component inspection and
# then reconciled visually against the error-bar markers.  The two overlapped
# CF/CF2 centers at -0.9 s share a plotted center.  Error-bar semantics are not
# stated in the article and are therefore not silently promoted to sigma.
OES_POINTS = (
    *_oes_series(
        "14b",
        "CF",
        (
            1609.0, 1458.2, 1298.4, 1530.3, 1521.7,
            1580.1, 1561.0, 1611.6, 1605.7, 1619.2,
            1627.6, 1632.0, 1676.6, 1687.8, 1689.4,
            1670.0, 1683.3, 1676.4, 1688.3, 1682.9,
        ),
    ),
    *_oes_series(
        "14b",
        "CF2",
        (
            1634.9, 1481.7, 1326.4, 1507.0, 1490.5,
            1580.1, 1525.3, 1559.6, 1532.4, 1597.7,
            1590.5, 1595.2, 1655.1, 1671.4, 1716.1,
            1690.9, 1704.7, 1723.9, 1757.0, 1742.5,
        ),
    ),
    *_oes_series(
        "14e",
        "F",
        (
            3992.4, 4105.6, 3964.0, 3960.0, 4005.2,
            3995.1, 4007.3, 3972.6, 3992.4, 3950.5,
            3963.8, 3980.6, 3708.0, 3660.4, 3810.9,
            3809.5, 3801.2, 3817.7, 3908.1, 3926.2,
        ),
    ),
)


ELECTRON_FIELDS = (
    "source_figure", "cycle_duration_s", "timing", "bias_start_s",
    "bias_end_s", "start_density_1e10_cm3", "mid_density_1e10_cm3",
    "end_density_1e10_cm3", "simpson_average_density_1e10_cm3",
    "digitization_window_uncertainty_1e10_cm3", "sample_x_px",
    "sample_y_px", "quantity_semantics", "supports_absolute_ion_flux",
    "boundary_evidence_tier", "observation_id", "source_image_sha256",
)
OES_FIELDS = (
    "source_figure", "panel", "signal", "time_s",
    "emission_intensity_ratio_to_Ar", "marker_center_x_px",
    "marker_center_y_px", "digitization_uncertainty_au",
    "quantity_semantics", "supports_absolute_ground_state_flux",
    "boundary_evidence_tier", "observation_id", "source_image_sha256",
)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def electron_rows():
    rows = []
    for window in BIAS_WINDOWS:
        densities = window.sample_density_1e10_cm3
        rows.append({
            "source_figure": f"Figure {window.panel}",
            "cycle_duration_s": str(window.cycle_duration_s),
            "timing": window.timing,
            "bias_start_s": f"{window.bias_start_s:.2f}",
            "bias_end_s": f"{window.bias_end_s:.2f}",
            "start_density_1e10_cm3": f"{densities[0]:.4f}",
            "mid_density_1e10_cm3": f"{densities[1]:.4f}",
            "end_density_1e10_cm3": f"{densities[2]:.4f}",
            "simpson_average_density_1e10_cm3": (
                f"{window.simpson_average_1e10_cm3:.4f}"
            ),
            "digitization_window_uncertainty_1e10_cm3": (
                f"{window.window_uncertainty_1e10_cm3:.2f}"
            ),
            "sample_x_px": ";".join(
                f"{value:.2f}" for value in window.sample_x_px
            ),
            "sample_y_px": ";".join(
                f"{value:.1f}" for value in window.sample_y_px
            ),
            "quantity_semantics": (
                "bulk_electron_density_window_not_positive_ion_flux"
            ),
            "supports_absolute_ion_flux": "false",
            "boundary_evidence_tier": "B_facility_diagnostic",
            "observation_id": (
                f"yoshie_ne_{window.cycle_duration_s}s_{window.timing}"
            ),
            "source_image_sha256": SOURCE[12]["sha256"],
        })
    return rows


def oes_rows():
    rows = []
    for point in OES_POINTS:
        rows.append({
            "source_figure": "Figure 14",
            "panel": point.panel,
            "signal": point.signal,
            "time_s": f"{point.time_s:.1f}",
            "emission_intensity_ratio_to_Ar": f"{point.ratio_to_ar:.4f}",
            "marker_center_x_px": f"{point.x_px:.2f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "digitization_uncertainty_au": (
                "0.08" if point.signal == "F" else "0.20"
            ),
            "quantity_semantics": (
                "optical_emission_ratio_not_ground_state_density_or_flux"
            ),
            "supports_absolute_ground_state_flux": "false",
            "boundary_evidence_tier": "B_facility_diagnostic",
            "observation_id": (
                f"yoshie_oes_{point.signal}_"
                f"{point.time_s:+.1f}".replace("+", "p").replace("-", "m")
            ),
            "source_image_sha256": SOURCE[14]["sha256"],
        })
    return rows


def _csv_text(fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def electron_csv_text():
    return _csv_text(ELECTRON_FIELDS, electron_rows())


def oes_csv_text():
    return _csv_text(OES_FIELDS, oes_rows())


def manifest(electron_payload, oes_payload):
    averages = {
        f"{row['cycle_duration_s']}s_{row['timing']}": float(
            row["simpson_average_density_1e10_cm3"]
        )
        for row in electron_rows()
    }
    return {
        "manifest_id": "YOSHIE-2023-FIGURES12-14-REACTOR-STATE-R1",
        "source": {
            "citation": (
                "T. Yoshie et al., Applied Surface Science 638 (2023) "
                "157981, DOI 10.1016/j.apsusc.2023.157981"
            ),
            "license": "CC BY 4.0",
            "figures": {
                str(number): {
                    "url": details["url"],
                    "sha256": details["sha256"],
                    "size_px": list(details["size"]),
                }
                for number, details in SOURCE.items()
            },
        },
        "pixel_calibration": {
            name: {
                "x_left_px": axis.x_left_px,
                "x_right_px": axis.x_right_px,
                "x_left_value": axis.x_left_value,
                "x_right_value": axis.x_right_value,
                "y_zero_px": axis.y_zero_px,
                "y_reference_px": axis.y_reference_px,
                "y_zero_value": axis.y_zero_value,
                "y_reference_value": axis.y_reference_value,
            }
            for name, axis in AXES.items()
        },
        "digitization": {
            "method": (
                "official publisher JPEG; PIL/NumPy checksum and dark-axis "
                "checks; full-resolution centerline/color-component "
                "transcription; visual overlay reconciliation"
            ),
            "figure12_sampling": (
                "start/mid/end of each published 0.25 s bias window; "
                "Simpson average"
            ),
            "figure14_sampling": (
                "full marker board at 0.2 s spacing for CF, CF2, and F"
            ),
            "figure14_error_bar_semantics": "not_reported_by_source",
            "source_figures_visually_inspected_at_full_resolution": True,
        },
        "derived_checks": {
            "electron_bias_windows": len(BIAS_WINDOWS),
            "oes_markers": len(OES_POINTS),
            "bias_window_average_density_1e10_cm3": averages,
            "eight_second_timing_II_crosses_fast_density_collapse": True,
        },
        "claim_boundary": {
            "valid": [
                (
                    "relative phase and timing constraints on the measured "
                    "bulk electron-density trace"
                ),
                (
                    "relative phase-resolved CF, CF2, and F optical-emission "
                    "constraints"
                ),
                (
                    "tests of a reactor model against independently measured "
                    "diagnostic shapes before feature-depth evaluation"
                ),
            ],
            "not_valid": [
                "positive-ion density or ion flux inferred from electron density alone",
                "species-resolved positive-ion flux",
                "ground-state CF, CF2, or F density or wafer flux inferred from OES alone",
                "an ion energy-angle distribution",
                "a depth calibration or permission to fit Figure 5 or Figure 6",
            ],
        },
        "outputs": {
            "electron_density_csv": {
                "path": (
                    "data/experimental/yoshie_2023/"
                    "figure12_bias_window_electron_density.csv"
                ),
                "sha256": hashlib.sha256(
                    electron_payload.encode("utf-8")
                ).hexdigest(),
            },
            "oes_csv": {
                "path": (
                    "data/experimental/yoshie_2023/"
                    "figure14_phase_resolved_oes.csv"
                ),
                "sha256": hashlib.sha256(
                    oes_payload.encode("utf-8")
                ).hexdigest(),
            },
        },
    }


def manifest_text(electron_payload, oes_payload):
    return json.dumps(
        manifest(electron_payload, oes_payload),
        indent=2,
        sort_keys=True,
    ) + "\n"


def expected_files():
    electron = electron_csv_text()
    oes = oes_csv_text()
    return {
        ELECTRON_CSV: electron,
        OES_CSV: oes,
        MANIFEST_PATH: manifest_text(electron, oes),
    }


def verify_sources(source_dir):
    images = {}
    for number, details in SOURCE.items():
        path = Path(source_dir) / details["name"]
        if _sha256(path) != details["sha256"]:
            raise RuntimeError(f"Figure {number} source checksum mismatch")
        image = Image.open(path).convert("RGB")
        if image.size != details["size"]:
            raise RuntimeError(
                f"Figure {number} size {image.size} != {details['size']}"
            )
        images[number] = image

    gray12 = np.asarray(images[12].convert("L"))
    gray14 = np.asarray(images[14].convert("L"))
    checks = {
        "12a left axis": np.mean(
            np.min(gray12[496:1303, 224:230], axis=1) < 96),
        "12a bottom axis": np.mean(
            np.min(gray12[1298:1305, 225:884], axis=0) < 96),
        "12b left axis": np.mean(
            np.min(gray12[496:1303, 1297:1303], axis=1) < 96),
        "12b bottom axis": np.mean(
            np.min(gray12[1298:1305, 1298:2618], axis=0) < 96),
        "14b left axis": np.mean(
            np.min(gray14[1151:1791, 214:221], axis=1) < 96),
        "14b bottom axis": np.mean(
            np.min(gray14[1786:1793, 216:1260], axis=0) < 96),
        "14e left axis": np.mean(
            np.min(gray14[3558:4198, 214:221], axis=1) < 96),
        "14e bottom axis": np.mean(
            np.min(gray14[4193:4200, 216:1260], axis=0) < 96),
    }
    failed = {
        name: float(value) for name, value in checks.items() if value < 0.45
    }
    if failed:
        raise RuntimeError(f"axis dark-pixel verification failed: {failed}")
    return images


def verify_committed_files():
    for path, expected in expected_files().items():
        if path.read_text(encoding="utf-8") != expected:
            raise RuntimeError(
                f"{path.relative_to(ROOT)} is not reproduced by this script"
            )


def write_overlay(images, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    overlay12 = images[12].copy()
    draw12 = ImageDraw.Draw(overlay12)
    timing_colors = {
        (4, "I"): "#ff00ff",
        (4, "II"): "#ff7f00",
        (4, "III"): "#377eb8",
        (8, "I"): "#ff00ff",
        (8, "II"): "#e41a1c",
        (8, "III"): "#377eb8",
        (8, "IV"): "#00aa76",
    }
    for window in BIAS_WINDOWS:
        color = timing_colors[(window.cycle_duration_s, window.timing)]
        for x_px, y_px in zip(window.sample_x_px, window.sample_y_px):
            radius = 11
            draw12.ellipse(
                (
                    x_px - radius, y_px - radius,
                    x_px + radius, y_px + radius,
                ),
                outline=color,
                width=3,
            )
            draw12.line(
                (x_px - 8, y_px, x_px + 8, y_px),
                fill=color,
                width=3,
            )
            draw12.line(
                (x_px, y_px - 8, x_px, y_px + 8),
                fill=color,
                width=3,
            )
    overlay12.crop((180, 470, 2660, 1325)).save(
        output_dir / "yoshie_figure12_bias_window_overlay.png"
    )

    overlay14 = images[14].copy()
    draw14 = ImageDraw.Draw(overlay14)
    signal_colors = {"CF": "#ff00ff", "CF2": "#0066ff", "F": "#ff00ff"}
    for point in OES_POINTS:
        radius = 9
        draw14.ellipse(
            (
                point.x_px - radius, point.y_px - radius,
                point.x_px + radius, point.y_px + radius,
            ),
            outline=signal_colors[point.signal],
            width=3,
        )
    overlay14.crop((180, 1130, 1280, 4210)).save(
        output_dir / "yoshie_figure14_oes_overlay.png"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--overlay-dir", type=Path)
    parser.add_argument(
        "--print-files",
        action="store_true",
        help="print generated file payloads instead of verifying committed files",
    )
    args = parser.parse_args()
    if args.print_files:
        for path, payload in expected_files().items():
            print(f"===== {path.relative_to(ROOT)} =====")
            print(payload, end="")
    else:
        verify_committed_files()
    if args.source_dir is not None:
        images = verify_sources(args.source_dir)
        if args.overlay_dir is not None:
            write_overlay(images, args.overlay_dir)
    elif args.overlay_dir is not None:
        raise SystemExit("--overlay-dir requires --source-dir")
    print(json.dumps({
        "status": "verified",
        "electron_bias_windows": len(BIAS_WINDOWS),
        "oes_markers": len(OES_POINTS),
        "claim": "diagnostic constraints, not absolute wafer fluxes",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
