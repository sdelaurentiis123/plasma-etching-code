#!/usr/bin/env python3
"""Replay Yoshie et al. 2023 Figures 4--6 with checksum and pixel QA.

Figure 4 supplies same-reactor blanket poly-Si rates for the four selected
bias timings.  Figures 5 and 6 supply the complete width-by-timing feature
board.  The paper reports etch rate per cumulative bias-on time; fixed cycle
counts and a 0.25 s bias pulse make the corresponding absolute depths
deterministic.

The publisher rasters are not redistributed.  Download the three official
figure assets into one directory and pass it with ``--source-dir``.  This
script verifies their SHA-256 values and dimensions before inspecting axes or
drawing an overlay.
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
BLANKET_CSV = OUTPUT_DIR / "figure4_blanket_poly_si_rates.csv"
FEATURE_CSV = OUTPUT_DIR / "figures5_6_feature_depths.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"
DEFAULT_SOURCE_DIR = ROOT / "tmp" / "sources" / "yoshie_2023"

SOURCE = {
    4: {
        "name": "yoshie_fig4.jpg",
        "url": (
            "https://ars.els-cdn.com/content/image/"
            "1-s2.0-S0169433223016604-gr4_lrg.jpg"
        ),
        "sha256": (
            "297013dbb00ef11e7152355e8ec91fdb7e4bc751f6c3d39d56c8c00dde29b3df"
        ),
        "size": (3150, 1778),
    },
    5: {
        "name": "yoshie_fig5.jpg",
        "url": (
            "https://ars.els-cdn.com/content/image/"
            "1-s2.0-S0169433223016604-gr5_lrg.jpg"
        ),
        "sha256": (
            "de722d00eadba115a690562b519387ea79e3a634f19a12f6252df168c6056aef"
        ),
        "size": (2389, 3750),
    },
    6: {
        "name": "yoshie_fig6.jpg",
        "url": (
            "https://ars.els-cdn.com/content/image/"
            "1-s2.0-S0169433223016604-gr6_lrg.jpg"
        ),
        "sha256": (
            "202e04674548182e40d2be4c7677cb1b1124cd9056f87baa13a1525495a9f7fe"
        ),
        "size": (2206, 4337),
    },
}


@dataclass(frozen=True)
class Axis:
    y_zero_px: float
    y_reference_px: float
    reference_rate: float

    def rate(self, y_px):
        return self.reference_rate * (
            self.y_zero_px - y_px
        ) / (
            self.y_zero_px - self.y_reference_px
        )


AXES = {
    "4a": Axis(1581.0, 732.0, 500.0),
    "4b": Axis(1581.0, 720.0, 600.0),
    "5e": Axis(778.0, 83.5, 300.0),
    "6f": Axis(718.5, 77.5, 700.0),
}


@dataclass(frozen=True)
class BlanketPoint:
    panel: str
    cycle_duration_s: int
    timing: str
    bias_start_s: float
    x_px: float
    y_px: float


BLANKET_POINTS = (
    BlanketPoint("4a", 4, "I", -0.15, 623.0, 1373.5),
    BlanketPoint("4a", 4, "II", 0.35, 718.0, 791.5),
    BlanketPoint("4a", 4, "III", 1.35, 908.0, 1184.5),
    BlanketPoint("4b", 8, "I", -0.15, 2245.0, 1368.0),
    BlanketPoint("4b", 8, "II", 0.10, 2292.5, 730.5),
    BlanketPoint("4b", 8, "III", 1.10, 2482.5, 1112.0),
    BlanketPoint("4b", 8, "IV", 2.10, 2674.0, 1336.0),
)


@dataclass(frozen=True)
class FeaturePoint:
    panel: str
    cycle_duration_s: int
    timing: str
    initial_width_nm: int
    x_px: float
    y_px: float
    upper_cap_y_px: float
    lower_cap_y_px: float


WIDTHS = (80, 90, 100, 110, 120, 130, 140)
X5 = (1513.0, 1655.0, 1795.0, 1936.0, 2073.0, 2216.0, 2347.0)
X6 = (1397.0, 1521.0, 1656.0, 1787.0, 1919.0, 2042.0, 2171.0)


def _feature_series(panel, cycle, timing, xs, centers, uppers, lowers):
    return tuple(
        FeaturePoint(panel, cycle, timing, width, x, center, upper, lower)
        for width, x, center, upper, lower in zip(
            WIDTHS, xs, centers, uppers, lowers
        )
    )


FEATURE_POINTS = (
    *_feature_series(
        "5e", 4, "I", X5,
        (435.5, 444.0, 451.0, 451.0, 450.5, 450.5, 451.0),
        (412.0, 420.0, 427.0, 427.0, 426.0, 426.0, 427.0),
        (461.0, 468.0, 475.0, 475.0, 475.0, 475.0, 475.0),
    ),
    *_feature_series(
        "5e", 4, "II", X5,
        (343.0, 272.5, 239.0, 214.5, 199.5, 184.5, 166.0),
        (306.0, 246.0, 214.0, 189.0, 174.0, 158.0, 140.0),
        (364.0, 299.0, 264.0, 240.0, 225.0, 211.0, 192.0),
    ),
    *_feature_series(
        "5e", 4, "III", X5,
        (395.0, 378.0, 369.0, 368.5, 369.5, 368.5, 370.5),
        (370.0, 356.0, 345.0, 338.0, 344.0, 344.0, 344.0),
        (419.0, 404.0, 395.0, 393.0, 395.0, 393.0, 396.0),
    ),
    *_feature_series(
        "6f", 8, "I", X6,
        (363.0, 348.5, 343.0, 338.5, 333.5, 333.0, 327.0),
        (351.0, 337.0, 332.0, 326.0, 322.0, 322.0, 316.0),
        (375.0, 360.0, 354.0, 350.0, 344.0, 344.0, 339.0),
    ),
    *_feature_series(
        "6f", 8, "II", X6,
        (329.0, 279.0, 255.0, 241.0, 227.5, 217.0, 205.5),
        (301.0, 265.0, 241.0, 227.0, 215.0, 203.0, 192.0),
        (343.0, 293.0, 269.0, 254.0, 240.0, 231.0, 219.0),
    ),
    *_feature_series(
        "6f", 8, "III", X6,
        (424.5, 395.5, 383.0, 368.5, 358.0, 353.5, 349.0),
        (411.0, 384.0, 369.0, 356.0, 347.0, 343.0, 338.0),
        (438.0, 407.0, 395.0, 380.0, 369.0, 364.0, 359.0),
    ),
    *_feature_series(
        "6f", 8, "IV", X6,
        (519.5, 514.5, 514.5, 509.5, 509.5, 509.5, 510.0),
        (505.0, 502.0, 502.0, 496.0, 496.0, 496.0, 496.0),
        (534.0, 527.0, 527.0, 523.0, 523.0, 523.0, 524.0),
    ),
)

BLANKET_FIELDS = (
    "source_figure", "cycle_duration_s", "timing", "bias_start_s",
    "blanket_etch_rate_nm_per_bias_min", "marker_center_x_px",
    "marker_center_y_px", "digitization_rate_uncertainty_nm_per_bias_min",
    "split", "boundary_evidence_tier", "observation_id",
    "source_image_sha256",
)
FEATURE_FIELDS = (
    "source_figure", "cycle_duration_s", "timing", "initial_width_nm",
    "etch_rate_nm_per_bias_min", "cumulative_bias_time_s",
    "derived_etch_depth_nm", "plotted_lower_rate_nm_per_bias_min",
    "plotted_upper_rate_nm_per_bias_min", "marker_center_x_px",
    "marker_center_y_px", "upper_cap_y_px", "lower_cap_y_px",
    "digitization_rate_uncertainty_nm_per_bias_min",
    "digitization_depth_uncertainty_nm", "measurement_uncertainty_semantics",
    "split", "boundary_evidence_tier", "observation_id",
    "source_image_sha256",
)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def blanket_rows():
    result = []
    for point in BLANKET_POINTS:
        figure = int(point.panel[0])
        rate = AXES[point.panel].rate(point.y_px)
        result.append({
            "source_figure": f"Figure {point.panel}",
            "cycle_duration_s": str(point.cycle_duration_s),
            "timing": point.timing,
            "bias_start_s": f"{point.bias_start_s:.2f}",
            "blanket_etch_rate_nm_per_bias_min": f"{rate:.3f}",
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "digitization_rate_uncertainty_nm_per_bias_min": "3.0",
            "split": "calibration",
            "boundary_evidence_tier": "B_facility_conditioned",
            "observation_id": (
                f"yoshie_blanket_{point.cycle_duration_s}s_{point.timing}"
            ),
            "source_image_sha256": SOURCE[figure]["sha256"],
        })
    return result


def feature_rows():
    result = []
    for point in FEATURE_POINTS:
        figure = int(point.panel[0])
        cumulative_bias_s = 168.75 if figure == 5 else 112.5
        axis = AXES[point.panel]
        rate = axis.rate(point.y_px)
        upper = axis.rate(point.upper_cap_y_px)
        lower = axis.rate(point.lower_cap_y_px)
        rate_uncertainty = 2.0 if figure == 5 else 4.0
        depth = rate * cumulative_bias_s / 60.0
        depth_uncertainty = rate_uncertainty * cumulative_bias_s / 60.0
        result.append({
            "source_figure": f"Figure {point.panel}",
            "cycle_duration_s": str(point.cycle_duration_s),
            "timing": point.timing,
            "initial_width_nm": str(point.initial_width_nm),
            "etch_rate_nm_per_bias_min": f"{rate:.3f}",
            "cumulative_bias_time_s": f"{cumulative_bias_s:g}",
            "derived_etch_depth_nm": f"{depth:.3f}",
            "plotted_lower_rate_nm_per_bias_min": f"{lower:.3f}",
            "plotted_upper_rate_nm_per_bias_min": f"{upper:.3f}",
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "upper_cap_y_px": f"{point.upper_cap_y_px:.1f}",
            "lower_cap_y_px": f"{point.lower_cap_y_px:.1f}",
            "digitization_rate_uncertainty_nm_per_bias_min": (
                f"{rate_uncertainty:.1f}"
            ),
            "digitization_depth_uncertainty_nm": f"{depth_uncertainty:.3f}",
            "measurement_uncertainty_semantics": "not_reported",
            "split": "held_out_transfer",
            "boundary_evidence_tier": "B_facility_conditioned",
            "observation_id": (
                f"yoshie_fig{figure}_{point.timing}_w"
                f"{point.initial_width_nm}"
            ),
            "source_image_sha256": SOURCE[figure]["sha256"],
        })
    return result


def _csv_text(fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def blanket_csv_text():
    return _csv_text(BLANKET_FIELDS, blanket_rows())


def feature_csv_text():
    return _csv_text(FEATURE_FIELDS, feature_rows())


def manifest(blanket_payload, feature_payload):
    feature = feature_rows()
    maxima = {
        f"{row['cycle_duration_s']}s_{row['timing']}": max(
            float(item["etch_rate_nm_per_bias_min"])
            for item in feature
            if item["cycle_duration_s"] == row["cycle_duration_s"]
            and item["timing"] == row["timing"]
        )
        for row in feature
    }
    return {
        "manifest_id": "YOSHIE-2023-FIGURES4-6-DEPTH-R1",
        "preregistration_commit_sha256": (
            "dec49c464af246b04f0168588dc436a76442e64c0294248b324ceb96cf9d0c62"
        ),
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
        "experiment": {
            "reactor": "NLD ICP, 150 mm wafer",
            "pressure_Pa": 4.0,
            "continuous_Ar_sccm": 80.0,
            "C4F8_pulse_sccm_s": [15.0, 1.0],
            "SF6_pulse_sccm_s": [15.0, 1.0],
            "source_power_W": 400.0,
            "substrate_bias_power_W": 50.0,
            "bias_pulse_duration_s": 0.25,
            "initial_trench_widths_nm": list(WIDTHS),
            "figure5": {
                "cycle_duration_s": 4,
                "cycles": 675,
                "wall_time_min": 45,
                "cumulative_bias_time_s": 168.75,
            },
            "figure6": {
                "cycle_duration_s": 8,
                "cycles": 450,
                "wall_time_min": 60,
                "cumulative_bias_time_s": 112.5,
            },
        },
        "pixel_calibration": {
            name: {
                "y_zero_px": axis.y_zero_px,
                "y_reference_px": axis.y_reference_px,
                "reference_rate_nm_per_bias_min": axis.reference_rate,
            }
            for name, axis in AXES.items()
        },
        "digitization": {
            "method": (
                "official publisher JPEG; PIL/NumPy dark-axis localization; "
                "full-resolution color-component marker and error-cap "
                "transcription; visual overlay reconciliation"
            ),
            "figure4_rate_uncertainty_nm_per_bias_min": 3.0,
            "figure5_rate_uncertainty_nm_per_bias_min": 2.0,
            "figure6_rate_uncertainty_nm_per_bias_min": 4.0,
            "plotted_error_bar_semantics": "not reported by the source",
            "machine_digitization_created_before_preregistration": False,
            "source_figure_visually_inspected_before_preregistration": True,
        },
        "text_cross_checks": {
            "figure4a_reported_maximum_nm_per_bias_min": 466.0,
            "figure4b_reported_maximum_nm_per_bias_min": 591.0,
            "digitized_figure4a_timing_II": float(
                next(
                    row["blanket_etch_rate_nm_per_bias_min"]
                    for row in blanket_rows()
                    if row["cycle_duration_s"] == "4"
                    and row["timing"] == "II"
                )
            ),
            "digitized_figure4b_timing_II": float(
                next(
                    row["blanket_etch_rate_nm_per_bias_min"]
                    for row in blanket_rows()
                    if row["cycle_duration_s"] == "8"
                    and row["timing"] == "II"
                )
            ),
        },
        "derived_checks": {
            "blanket_point_count": len(BLANKET_POINTS),
            "feature_point_count": len(FEATURE_POINTS),
            "feature_series_maxima_nm_per_bias_min": maxima,
        },
        "claim_boundary": {
            "valid": (
                "Figure 4 may condition a same-reactor boundary observable "
                "before fixed-duration feature prediction only inside an "
                "independently validated material- and cycle-history-resolved "
                "mechanism"
            ),
            "not_valid": [
                "a multiplicative blanket-to-feature depth scale",
                "a measured species-resolved wafer flux",
                "a measured ion energy-angle distribution",
                "a first-principles knobs-to-flux validation",
                "permission to fit any Figure 5 or Figure 6 depth",
            ],
        },
        "outputs": {
            "blanket_csv": {
                "path": (
                    "data/experimental/yoshie_2023/"
                    "figure4_blanket_poly_si_rates.csv"
                ),
                "sha256": hashlib.sha256(
                    blanket_payload.encode("utf-8")
                ).hexdigest(),
            },
            "feature_csv": {
                "path": (
                    "data/experimental/yoshie_2023/"
                    "figures5_6_feature_depths.csv"
                ),
                "sha256": hashlib.sha256(
                    feature_payload.encode("utf-8")
                ).hexdigest(),
            },
        },
    }


def manifest_text(blanket_payload, feature_payload):
    return json.dumps(
        manifest(blanket_payload, feature_payload),
        indent=2,
        sort_keys=True,
    ) + "\n"


def expected_files():
    blanket = blanket_csv_text()
    feature = feature_csv_text()
    return {
        BLANKET_CSV: blanket,
        FEATURE_CSV: feature,
        MANIFEST_PATH: manifest_text(blanket, feature),
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

    gray5 = np.asarray(images[5].convert("L"))
    gray6 = np.asarray(images[6].convert("L"))
    checks = {
        "5e left axis": np.mean(
            np.min(gray5[83:780, 1509:1514], axis=1) < 64),
        "5e right axis": np.mean(
            np.min(gray5[83:780, 2348:2353], axis=1) < 64),
        "5e top axis": np.mean(
            np.min(gray5[82:86, 1510:2352], axis=0) < 96),
        "5e bottom axis": np.mean(
            np.min(gray5[777:781, 1510:2352], axis=0) < 96),
        "6f left axis": np.mean(
            np.min(gray6[76:721, 1394:1399], axis=1) < 96),
        "6f right axis": np.mean(
            np.min(gray6[76:721, 2168:2173], axis=1) < 96),
        "6f top axis": np.mean(
            np.min(gray6[75:81, 1395:2172], axis=0) < 96),
        "6f bottom axis": np.mean(
            np.min(gray6[717:722, 1395:2172], axis=0) < 96),
    }
    failed = {
        name: float(value) for name, value in checks.items() if value < 0.65
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
    colors = {
        (5, "I"): "#ff7f00",
        (5, "II"): "#e41a1c",
        (5, "III"): "#377eb8",
        (6, "I"): "#ff7f00",
        (6, "II"): "#e41a1c",
        (6, "III"): "#377eb8",
        (6, "IV"): "#00aa76",
    }
    for number, crop_box in {
        4: (200, 620, 3080, 1600),
        5: (1450, 60, 2380, 800),
        6: (1330, 50, 2200, 740),
    }.items():
        overlay = images[number].copy()
        draw = ImageDraw.Draw(overlay)
        if number == 4:
            points = BLANKET_POINTS
        else:
            points = [
                point for point in FEATURE_POINTS
                if int(point.panel[0]) == number
            ]
        for point in points:
            color = "#ff00ff" if number == 4 else colors[
                (number, point.timing)
            ]
            radius = 12
            draw.ellipse(
                (
                    point.x_px - radius,
                    point.y_px - radius,
                    point.x_px + radius,
                    point.y_px + radius,
                ),
                outline=color,
                width=3,
            )
            draw.line(
                (
                    point.x_px - 10, point.y_px,
                    point.x_px + 10, point.y_px,
                ),
                fill=color,
                width=3,
            )
            draw.line(
                (
                    point.x_px, point.y_px - 10,
                    point.x_px, point.y_px + 10,
                ),
                fill=color,
                width=3,
            )
            if isinstance(point, FeaturePoint):
                draw.line(
                    (
                        point.x_px + 15, point.upper_cap_y_px,
                        point.x_px + 15, point.lower_cap_y_px,
                    ),
                    fill=color,
                    width=2,
                )
        overlay.crop(crop_box).save(
            output_dir / f"yoshie_figure{number}_overlay.png"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--overlay-dir", type=Path)
    parser.add_argument(
        "--print-files", action="store_true",
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
        "preregistration_commit_sha256": (
            "dec49c464af246b04f0168588dc436a76442e64c0294248b324ceb96cf9d0c62"
        ),
        "blanket_points": len(BLANKET_POINTS),
        "held_out_feature_points": len(FEATURE_POINTS),
        "figure4a_maximum_cross_check": AXES["4a"].rate(791.5),
        "figure4b_maximum_cross_check": AXES["4b"].rate(730.5),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
