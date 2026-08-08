#!/usr/bin/env python3
"""Replay volume-average electron-density markers in Malyshev Figure 11.

The publisher PDF embeds Figure 11 as a native one-bit raster. Only marker
centers reconciled against the original pixels are transcribed. The linear
fits are excluded, as is the inseparable low-power overlap in panel b.
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
CSV_PATH = OUTPUT_DIRECTORY / "figure11_electron_density.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "electron_density_manifest.json"

SOURCE_PDF_SHA256 = (
    "569ab180bb8cab71dda0860e0350f03d24b391e084e19116fefb74f1c719789c")
SOURCE_DOI = "10.1063/1.368010"
SOURCE_PDF_PAGE = 8
NATIVE_IMAGE_SIZE = (872, 1097)
NATIVE_PIXEL_SHA256 = (
    "b941d2ab8d77dfc8034541234cfde780fc4f958e2f83ccf7a802eb789d645cbb")


@dataclass(frozen=True)
class PanelAxis:
    gap_cm: float
    top_px: float
    bottom_px: float
    maximum_electron_density_cm3: float
    left_px: float = 182.0
    right_px: float = 870.0
    maximum_power_W: float = 950.0

    def power_at_pixel(self, x_px: float) -> float:
        return (
            (float(x_px) - self.left_px)
            / (self.right_px - self.left_px)
            * self.maximum_power_W
        )

    def density_at_pixel(self, y_px: float) -> float:
        return (
            (self.bottom_px - float(y_px))
            / (self.bottom_px - self.top_px)
            * self.maximum_electron_density_cm3
        )


PANELS = {
    11.0: PanelAxis(
        gap_cm=11.0,
        top_px=13.0,
        bottom_px=487.0,
        maximum_electron_density_cm3=2.5e11,
    ),
    6.5: PanelAxis(
        gap_cm=6.5,
        top_px=523.0,
        bottom_px=997.0,
        maximum_electron_density_cm3=1.1e11,
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
    # 11 cm, 0.5 mTorr: open diamonds.
    PixelPoint(11.0, 0.5, 325.0, 381.0, "open_diamond"),
    PixelPoint(11.0, 0.5, 396.0, 323.0, "open_diamond"),
    PixelPoint(11.0, 0.5, 679.0, 154.0, "open_diamond"),
    PixelPoint(11.0, 0.5, 796.0, 95.0, "open_diamond"),
    # 11 cm, 1 mTorr: filled down-triangles.
    PixelPoint(11.0, 1.0, 326.0, 394.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 397.0, 345.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 543.0, 237.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 687.0, 155.0, "filled_down_triangle"),
    PixelPoint(11.0, 1.0, 824.0, 88.0, "filled_down_triangle"),
    # 11 cm, 2 mTorr: open circles.
    PixelPoint(11.0, 2.0, 289.0, 435.0, "open_circle"),
    PixelPoint(11.0, 2.0, 326.0, 411.0, "open_circle"),
    PixelPoint(11.0, 2.0, 399.0, 378.0, "open_circle"),
    PixelPoint(11.0, 2.0, 543.0, 290.0, "open_circle"),
    PixelPoint(11.0, 2.0, 687.0, 196.0, "open_circle"),
    PixelPoint(11.0, 2.0, 833.0, 109.0, "open_circle"),
    # 11 cm, 10 mTorr: filled squares.
    PixelPoint(11.0, 10.0, 290.0, 454.0, "filled_square"),
    PixelPoint(11.0, 10.0, 399.0, 365.0, "filled_square"),
    PixelPoint(11.0, 10.0, 543.0, 357.0, "filled_square"),
    PixelPoint(11.0, 10.0, 823.0, 244.0, "filled_square"),
    # 6.5 cm, 2 mTorr: open circles. The merged first pair is omitted.
    PixelPoint(6.5, 2.0, 399.0, 872.0, "open_circle"),
    PixelPoint(6.5, 2.0, 543.0, 772.0, "open_circle"),
    PixelPoint(6.5, 2.0, 689.0, 668.0, "open_circle"),
    PixelPoint(6.5, 2.0, 806.0, 611.0, "open_circle"),
    PixelPoint(6.5, 2.0, 833.0, 572.0, "open_circle"),
    # 6.5 cm, 10 mTorr: filled squares. The merged first pair is omitted.
    PixelPoint(6.5, 10.0, 399.0, 914.0, "filled_square"),
    PixelPoint(6.5, 10.0, 543.0, 877.0, "filled_square"),
    PixelPoint(6.5, 10.0, 812.0, 811.0, "filled_square"),
)


FIELDNAMES = (
    "source_figure",
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "volume_average_electron_density_cm3",
    "marker",
    "marker_center_x_px",
    "marker_center_y_px",
    "digitization_power_uncertainty_W",
    "digitization_electron_density_uncertainty_cm3",
    "measurement_method",
    "volume_average_conversion",
    "reported_measurement_uncertainty",
    "tcp_power_semantics",
    "supports_local_wafer_electron_density",
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
        density_uncertainty = (
            3.0 * axis.maximum_electron_density_cm3
            / (axis.bottom_px - axis.top_px)
        )
        output.append({
            "source_figure": "Figure 11",
            "window_to_wafer_gap_cm": f"{point.gap_cm:g}",
            "pressure_mTorr": f"{point.pressure_mTorr:g}",
            "tcp_source_power_W": f"{axis.power_at_pixel(point.x_px):.3f}",
            "volume_average_electron_density_cm3": (
                f"{axis.density_at_pixel(point.y_px):.8e}"
            ),
            "marker": point.marker,
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "digitization_power_uncertainty_W": f"{power_uncertainty:.4f}",
            "digitization_electron_density_uncertainty_cm3": (
                f"{density_uncertainty:.8e}"
            ),
            "measurement_method": "Langmuir_probe_analysis_reported_elsewhere",
            "volume_average_conversion": (
                "radial_symmetry_and_axial_sin_pi_h_over_gap"
            ),
            "reported_measurement_uncertainty": "not_reported_in_article",
            "tcp_power_semantics": (
                "power_into_matching_network_not_absorbed_power"
            ),
            "supports_local_wafer_electron_density": "false",
            "supports_wafer_flux": "false",
            "validation_role": (
                "measured_volume_average_electron_state_conditioning_input"
            ),
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
        "manifest_id": "MALYSHEV-1998-LAM-ELECTRON-DENSITY-R1",
        "source": {
            "citation": (
                "M. V. Malyshev, V. M. Donnelly, A. Kornblit, and "
                "N. A. Ciampa, Journal of Applied Physics 84, 137-146 (1998)"
            ),
            "doi": SOURCE_DOI,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": SOURCE_PDF_PAGE,
            "figure": 11,
            "native_image": {
                "pdfimages_index_on_page": 0,
                "size_px": list(NATIVE_IMAGE_SIZE),
                "one_bit_pixel_sha256": NATIVE_PIXEL_SHA256,
            },
        },
        "measurement": {
            "observable": "volume_average_electron_density_cm3",
            "method": "Langmuir probe analysis reported elsewhere",
            "reported_uncertainty": "not reported in this article",
            "conversion": (
                "measured along a line 1.35 cm above wafer/chuck; converted "
                "assuming radial symmetry and axial sin(pi*h/gap)"
            ),
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
            "electron_density_range_cm3": [
                min(
                    float(row["volume_average_electron_density_cm3"])
                    for row in decoded
                ),
                max(
                    float(row["volume_average_electron_density_cm3"])
                    for row in decoded
                ),
            ],
            "csv_sha256": csv_sha256,
        },
        "exclusions": {
            "linear_lines": "fits/guides; never digitized as measurements",
            "low_power_panel_b_overlap": (
                "open-circle and filled-square identities merge in the "
                "one-bit raster and were omitted"
            ),
            "unreported_pressure_series_panel_b": (
                "0.5 and 1 mTorr are absent from panel b"
            ),
        },
        "use_boundary": {
            "supports_measured_volume_average_ne_conditioning": True,
            "supports_local_wafer_electron_density": False,
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
    prefix = directory / "malyshev_figure11_native"
    subprocess.run(
        [
            "pdfimages", "-f", str(SOURCE_PDF_PAGE), "-l",
            str(SOURCE_PDF_PAGE), "-png", str(source_pdf), str(prefix),
        ],
        check=True,
    )
    matches = [
        path for path in sorted(directory.glob("malyshev_figure11_native-*.png"))
        if Image.open(path).size == NATIVE_IMAGE_SIZE
    ]
    if len(matches) != 1:
        raise ValueError("native Figure 11 raster was not recovered uniquely")
    if _pixel_sha256(matches[0]) != NATIVE_PIXEL_SHA256:
        raise ValueError("native Figure 11 pixels do not match")
    return matches[0]


def _write_overlay(source: Path, overlay_path: Path) -> None:
    image = Image.open(source).convert("RGB")
    draw = ImageDraw.Draw(image)
    colors = {0.5: "magenta", 1.0: "green", 2.0: "blue", 10.0: "red"}
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
        MANIFEST_PATH.write_text(
            json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    else:
        print(csv_payload, end="")


if __name__ == "__main__":
    main()
