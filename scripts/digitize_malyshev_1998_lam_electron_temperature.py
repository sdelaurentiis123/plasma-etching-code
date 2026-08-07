#!/usr/bin/env python3
"""Replay measured electron-temperature markers in Malyshev Figure 3.

The publisher PDF embeds Figure 3 as a native one-bit raster. Only marker
centers reconciled against the original pixels are transcribed; the smooth
guide curves are excluded. The paper says the 11 cm values were reported
elsewhere and the 6.5 cm values were previously unpublished.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "malyshev_1998_lam"
CSV_PATH = OUTPUT_DIRECTORY / "figure3_electron_temperature.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "electron_temperature_manifest.json"
PACKAGE_CSV_PATH = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "malyshev_1998_lam_electron_temperature.csv"
)

SOURCE_PDF_SHA256 = (
    "569ab180bb8cab71dda0860e0350f03d24b391e084e19116fefb74f1c719789c")
SOURCE_DOI = "10.1063/1.368010"
SOURCE_PDF_PAGE = 4
NATIVE_IMAGE_SIZE = (859, 1057)
NATIVE_PIXEL_SHA256 = (
    "84d7d77205e51cf61a33391aad96bdbe1a8c7827fe8cb4961f41e4bc455053ad")


@dataclass(frozen=True)
class PanelAxis:
    gap_cm: float
    top_px: float
    bottom_px: float
    minimum_temperature_eV: float
    maximum_temperature_eV: float
    left_px: float = 118.0
    right_px: float = 856.0
    maximum_power_W: float = 950.0

    def power_at_pixel(self, x_px: float) -> float:
        return (
            (float(x_px) - self.left_px)
            / (self.right_px - self.left_px)
            * self.maximum_power_W
        )

    def temperature_at_pixel(self, y_px: float) -> float:
        fraction = (
            (self.bottom_px - float(y_px))
            / (self.bottom_px - self.top_px)
        )
        return (
            self.minimum_temperature_eV
            + fraction
            * (
                self.maximum_temperature_eV
                - self.minimum_temperature_eV
            )
        )


PANELS = {
    11.0: PanelAxis(
        gap_cm=11.0,
        top_px=0.0,
        bottom_px=475.0,
        minimum_temperature_eV=1.1,
        maximum_temperature_eV=3.7,
    ),
    6.5: PanelAxis(
        gap_cm=6.5,
        top_px=485.0,
        bottom_px=960.0,
        minimum_temperature_eV=1.4,
        maximum_temperature_eV=4.0,
    ),
}


@dataclass(frozen=True)
class PixelPoint:
    gap_cm: float
    pressure_mTorr: float
    x_px: float
    y_px: float
    marker: str


PIXEL_POINTS = (
    # 11 cm: open diamonds, 0.5 mTorr.
    PixelPoint(11.0, 0.5, 143.3, 251.8, "open_diamond"),
    PixelPoint(11.0, 0.5, 350.6, 160.7, "open_diamond"),
    PixelPoint(11.0, 0.5, 777.6, 65.1, "open_diamond"),
    # 11 cm: filled down-triangles, 1 mTorr.
    PixelPoint(11.0, 1.0, 144.3, 265.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 222.7, 343.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 273.2, 300.9, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 350.7, 236.8, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 505.9, 211.2, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 661.2, 154.9, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 785.5, 163.0, "filled_down_triangle"),
    # 11 cm: open circles, 2 mTorr.
    PixelPoint(11.0, 2.0, 144.9, 302.9, "open_circle"),
    PixelPoint(11.0, 2.0, 194.9, 337.1, "open_circle"),
    PixelPoint(11.0, 2.0, 233.7, 337.1, "open_circle"),
    PixelPoint(11.0, 2.0, 272.6, 334.6, "open_circle"),
    PixelPoint(11.0, 2.0, 349.8, 298.5, "open_circle"),
    PixelPoint(11.0, 2.0, 505.3, 254.6, "open_circle"),
    PixelPoint(11.0, 2.0, 660.6, 249.9, "open_circle"),
    PixelPoint(11.0, 2.0, 792.7, 230.3, "open_circle"),
    # 11 cm: filled squares, 10 mTorr.
    PixelPoint(11.0, 10.0, 148.2, 435.6, "filled_square"),
    PixelPoint(11.0, 10.0, 195.0, 425.0, "filled_square"),
    PixelPoint(11.0, 10.0, 234.5, 407.5, "filled_square"),
    PixelPoint(11.0, 10.0, 273.5, 403.6, "filled_square"),
    PixelPoint(11.0, 10.0, 351.0, 391.0, "filled_square"),
    PixelPoint(11.0, 10.0, 505.9, 372.5, "filled_square"),
    PixelPoint(11.0, 10.0, 583.9, 354.2, "filled_square"),
    PixelPoint(11.0, 10.0, 661.0, 356.5, "filled_square"),
    PixelPoint(11.0, 10.0, 738.5, 364.5, "filled_square"),
    PixelPoint(11.0, 10.0, 809.5, 322.0, "filled_square"),
    # 11 cm: open up-triangles, 20 mTorr. Low-power overlaps omitted.
    PixelPoint(11.0, 20.0, 195.7, 442.1, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 234.4, 436.9, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 273.3, 439.2, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 350.8, 421.7, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 505.9, 406.1, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 661.4, 387.2, "open_up_triangle"),
    PixelPoint(11.0, 20.0, 816.6, 376.2, "open_up_triangle"),
    # 6.5 cm: open diamonds, 0.5 mTorr. The overlapping first pair is omitted.
    PixelPoint(6.5, 0.5, 157.5, 674.5, "open_diamond"),
    PixelPoint(6.5, 0.5, 662.0, 559.0, "open_diamond"),
    PixelPoint(6.5, 0.5, 817.0, 525.5, "open_diamond"),
    # 6.5 cm: filled down-triangles, 1 mTorr.
    PixelPoint(6.5, 1.0, 157.5, 760.2, "filled_down_triangle"),
    PixelPoint(6.5, 1.0, 196.5, 760.0, "filled_down_triangle"),
    PixelPoint(6.5, 1.0, 662.0, 629.0, "filled_down_triangle"),
    PixelPoint(6.5, 1.0, 817.0, 550.2, "filled_down_triangle"),
    PixelPoint(6.5, 1.0, 816.9, 621.9, "filled_down_triangle"),
    # 6.5 cm: open circles, 2 mTorr.
    PixelPoint(6.5, 2.0, 137.2, 835.5, "open_circle"),
    PixelPoint(6.5, 2.0, 147.8, 789.4, "open_circle"),
    PixelPoint(6.5, 2.0, 176.2, 809.0, "open_circle"),
    PixelPoint(6.5, 2.0, 254.0, 798.0, "open_circle"),
    PixelPoint(6.5, 2.0, 303.5, 821.9, "open_circle"),
    PixelPoint(6.5, 2.0, 389.8, 774.2, "open_circle"),
    PixelPoint(6.5, 2.0, 506.0, 756.8, "open_circle"),
    PixelPoint(6.5, 2.0, 661.4, 702.9, "open_circle"),
    PixelPoint(6.5, 2.0, 816.6, 736.4, "open_circle"),
    # 6.5 cm: filled squares, 10 mTorr.
    PixelPoint(6.5, 10.0, 137.5, 894.5, "filled_square"),
    PixelPoint(6.5, 10.0, 150.4, 910.6, "filled_square"),
    PixelPoint(6.5, 10.0, 192.3, 926.6, "filled_square"),
    PixelPoint(6.5, 10.0, 235.0, 925.9, "filled_square"),
    PixelPoint(6.5, 10.0, 273.9, 911.8, "filled_square"),
    PixelPoint(6.5, 10.0, 351.5, 894.5, "filled_square"),
    PixelPoint(6.5, 10.0, 506.6, 876.9, "filled_square"),
    PixelPoint(6.5, 10.0, 662.0, 825.6, "filled_square"),
    PixelPoint(6.5, 10.0, 759.0, 812.0, "filled_square"),
    PixelPoint(6.5, 10.0, 800.0, 754.0, "filled_square"),
)


FIELDNAMES = (
    "source_figure",
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "electron_temperature_eV",
    "marker",
    "marker_center_x_px",
    "marker_center_y_px",
    "digitization_power_uncertainty_W",
    "digitization_temperature_uncertainty_eV",
    "measurement_method",
    "reported_measurement_uncertainty",
    "tcp_power_semantics",
    "supports_absorbed_power",
    "supports_wafer_flux",
    "validation_role",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _pixel_sha256(path: Path) -> str:
    image = Image.open(path).convert("1")
    return hashlib.sha256(image.tobytes()).hexdigest()


def rows() -> list[dict[str, str]]:
    output = []
    for point in PIXEL_POINTS:
        axis = PANELS[point.gap_cm]
        power_uncertainty = (
            2.0 * axis.maximum_power_W
            / (axis.right_px - axis.left_px)
        )
        temperature_uncertainty = (
            3.0
            * (
                axis.maximum_temperature_eV
                - axis.minimum_temperature_eV
            )
            / (axis.bottom_px - axis.top_px)
        )
        output.append({
            "source_figure": "Figure 3",
            "window_to_wafer_gap_cm": f"{point.gap_cm:g}",
            "pressure_mTorr": f"{point.pressure_mTorr:g}",
            "tcp_source_power_W": f"{axis.power_at_pixel(point.x_px):.3f}",
            "electron_temperature_eV": (
                f"{axis.temperature_at_pixel(point.y_px):.5f}"
            ),
            "marker": point.marker,
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "digitization_power_uncertainty_W": f"{power_uncertainty:.4f}",
            "digitization_temperature_uncertainty_eV": (
                f"{temperature_uncertainty:.5f}"
            ),
            "measurement_method": "optical_emission_spectroscopy",
            "reported_measurement_uncertainty": "not_reported_in_article",
            "tcp_power_semantics": (
                "power_into_matching_network_not_absorbed_power"
            ),
            "supports_absorbed_power": "false",
            "supports_wafer_flux": "false",
            "validation_role": "measured_electron_state_conditioning_input",
        })
    return output


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def manifest(csv_sha256: str) -> dict[str, object]:
    decoded = rows()
    return {
        "manifest_id": "MALYSHEV-1998-LAM-ELECTRON-TEMPERATURE-R1",
        "source": {
            "citation": (
                "M. V. Malyshev, V. M. Donnelly, A. Kornblit, and "
                "N. A. Ciampa, Journal of Applied Physics 84, 137-146 (1998)"
            ),
            "doi": SOURCE_DOI,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": SOURCE_PDF_PAGE,
            "figure": 3,
            "native_image": {
                "pdfimages_index_on_page": 1,
                "size_px": list(NATIVE_IMAGE_SIZE),
                "one_bit_pixel_sha256": NATIVE_PIXEL_SHA256,
            },
        },
        "measurement": {
            "observable": "electron_temperature_eV",
            "method": "optical emission spectroscopy",
            "reported_uncertainty": "not reported in this article",
            "eleven_cm_provenance": "reported elsewhere in source ref. 3",
            "six_point_five_cm_provenance": "previously unpublished",
        },
        "digitization": {
            "marker_count": len(PIXEL_POINTS),
            "method": (
                "native one-bit figure raster; original-resolution "
                "human-vision marker reconciliation; PIL axis and point replay"
            ),
            "horizontal_pixel_allowance": 2.0,
            "vertical_pixel_allowance": 3.0,
            "power_range_W": [
                min(float(row["tcp_source_power_W"]) for row in decoded),
                max(float(row["tcp_source_power_W"]) for row in decoded),
            ],
            "temperature_range_eV": [
                min(float(row["electron_temperature_eV"]) for row in decoded),
                max(float(row["electron_temperature_eV"]) for row in decoded),
            ],
            "csv_sha256": csv_sha256,
        },
        "exclusions": {
            "smooth_curves": "guide curves; never digitized as measurements",
            "low_power_overlap": (
                "markers whose identity could not be reconciled independently "
                "at original pixels were omitted"
            ),
            "20_mTorr_6_point_5_cm": "not present in source panel b",
        },
        "use_boundary": {
            "supports_measured_Te_conditioning": True,
            "supports_absorbed_power_validation": False,
            "supports_wafer_flux_validation": False,
            "supports_feature_depth_validation": False,
            "interpolation_rule": (
                "interpolate only within a fixed gap and pressure series; "
                "report interpolation separately from measurement"
            ),
        },
    }


def _extract_native_image(source_pdf: Path, directory: Path) -> Path:
    if _sha256(source_pdf) != SOURCE_PDF_SHA256:
        raise ValueError("source PDF checksum does not match Malyshev 1998")
    prefix = directory / "malyshev_figure3_native"
    subprocess.run(
        [
            "pdfimages", "-f", str(SOURCE_PDF_PAGE), "-l",
            str(SOURCE_PDF_PAGE), "-png", str(source_pdf), str(prefix),
        ],
        check=True,
    )
    matches = [
        path for path in sorted(directory.glob("malyshev_figure3_native-*.png"))
        if Image.open(path).size == NATIVE_IMAGE_SIZE
    ]
    if len(matches) != 1:
        raise ValueError("native Figure 3 raster was not recovered uniquely")
    if _pixel_sha256(matches[0]) != NATIVE_PIXEL_SHA256:
        raise ValueError("native Figure 3 pixels do not match")
    return matches[0]


def _write_overlay(source: Path, overlay_path: Path) -> None:
    image = Image.open(source).convert("RGB")
    draw = ImageDraw.Draw(image)
    colors = {0.5: "magenta", 1.0: "green", 2.0: "blue", 10.0: "red", 20.0: "orange"}
    for point in PIXEL_POINTS:
        radius = 6
        draw.ellipse(
            (
                point.x_px - radius,
                point.y_px - radius,
                point.x_px + radius,
                point.y_px + radius,
            ),
            outline=colors[point.pressure_mTorr],
            width=2,
        )
    overlay_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(overlay_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--overlay", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    csv_payload = csv_text()
    csv_digest = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    manifest_payload = manifest(csv_digest)

    if args.source_pdf is not None:
        with tempfile.TemporaryDirectory() as temporary:
            image = _extract_native_image(
                args.source_pdf.resolve(), Path(temporary))
            if args.overlay is not None:
                _write_overlay(image, args.overlay.resolve())
    elif args.overlay is not None:
        parser.error("--overlay requires --source-pdf")

    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(csv_payload, encoding="utf-8")
        PACKAGE_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
        PACKAGE_CSV_PATH.write_text(csv_payload, encoding="utf-8")
        MANIFEST_PATH.write_text(
            json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    else:
        print(csv_payload, end="")


if __name__ == "__main__":
    main()
