#!/usr/bin/env python3
"""Replay Malyshev et al. (1998) Lam chlorine-dissociation data.

Figures 7 and 8 report measured Cl2 density relative to the plasma-off
density in a Lam Research Alliance metal etcher.  The publisher PDF embeds
each figure as a native one-bit raster.  This replay extracts those images
directly with ``pdfimages`` so no PDF-page resampling enters the data path.

Only visually unambiguous marker centers are transcribed.  The smooth curves
are the paper's model and are deliberately excluded.  The paper's bars span
the independently reduced Ar and Xe actinometry estimates; they are not
silently treated as statistical standard deviations.
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
CSV_PATH = OUTPUT_DIRECTORY / "figures7_8_chlorine_dissociation.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_PDF_SHA256 = (
    "569ab180bb8cab71dda0860e0350f03d24b391e084e19116fefb74f1c719789c")
SOURCE_DOI = "10.1063/1.368010"
SOURCE_PDF_PAGE = 6


@dataclass(frozen=True)
class Axis:
    figure: int
    native_image_index: int
    image_size: tuple[int, int]
    pixel_sha256: str
    gap_cm: float
    left_px: float
    right_px: float
    top_px: float
    bottom_px: float
    power_max_W: float
    relative_cl2_max_percent: float

    def power_at_pixel(self, x_px: float) -> float:
        return (
            (float(x_px) - self.left_px)
            / (self.right_px - self.left_px)
            * self.power_max_W
        )

    def relative_cl2_at_pixel(self, y_px: float) -> float:
        return (
            (self.bottom_px - float(y_px))
            / (self.bottom_px - self.top_px)
            * self.relative_cl2_max_percent
        )


AXES = {
    7: Axis(
        figure=7,
        native_image_index=1,
        image_size=(892, 692),
        pixel_sha256=(
            "b8deff3a79739243d4eeaabad6f9f7a294a44ffc3b3fa20600cd4fe882171ce3"),
        gap_cm=11.0,
        left_px=109.0,
        right_px=890.0,
        top_px=1.0,
        bottom_px=604.0,
        power_max_W=950.0,
        relative_cl2_max_percent=140.0,
    ),
    8: Axis(
        figure=8,
        native_image_index=2,
        image_size=(893, 649),
        pixel_sha256=(
            "d58263f2dd8c924996cdc298a40f0eff548bb3a8e81689dca199f39eb8fef34c"),
        gap_cm=6.5,
        left_px=109.0,
        right_px=891.0,
        top_px=1.0,
        bottom_px=560.0,
        power_max_W=950.0,
        relative_cl2_max_percent=150.0,
    ),
}


@dataclass(frozen=True)
class PixelPoint:
    figure: int
    pressure_mTorr: float
    x_px: float
    y_px: float
    marker: str
    flow_condition: str = "normal"


PIXEL_POINTS = (
    # Figure 7, 11 cm gap. Filled squares: 10 mTorr.
    PixelPoint(7, 10.0, 191.0, 139.0, "filled_square"),
    PixelPoint(7, 10.0, 232.0, 161.0, "filled_square"),
    PixelPoint(7, 10.0, 274.0, 174.0, "filled_square"),
    PixelPoint(7, 10.0, 356.0, 159.0, "filled_square"),
    PixelPoint(7, 10.0, 520.0, 208.0, "filled_square"),
    PixelPoint(7, 10.0, 603.0, 160.0, "filled_square"),
    PixelPoint(7, 10.0, 685.0, 247.0, "filled_square"),
    PixelPoint(7, 10.0, 767.0, 367.0, "filled_square"),
    PixelPoint(7, 10.0, 849.0, 255.0, "filled_square"),
    # Open circles: 2 mTorr.
    PixelPoint(7, 2.0, 232.0, 208.0, "open_circle"),
    PixelPoint(7, 2.0, 356.0, 307.0, "open_circle"),
    PixelPoint(7, 2.0, 520.0, 359.0, "open_circle"),
    PixelPoint(7, 2.0, 685.0, 423.0, "open_circle"),
    PixelPoint(7, 2.0, 822.0, 480.0, "open_circle"),
    # Filled down-triangles: 1 mTorr.
    PixelPoint(7, 1.0, 232.0, 281.0, "filled_down_triangle"),
    PixelPoint(7, 1.0, 356.0, 324.0, "filled_down_triangle"),
    PixelPoint(7, 1.0, 520.0, 445.0, "filled_down_triangle"),
    PixelPoint(7, 1.0, 685.0, 475.0, "filled_down_triangle"),
    PixelPoint(7, 1.0, 808.0, 496.0, "filled_down_triangle"),
    # Open diamonds: 0.5 mTorr.
    PixelPoint(7, 0.5, 356.0, 410.0, "open_diamond"),
    PixelPoint(7, 0.5, 808.0, 518.0, "open_diamond"),
    # Figure 8, 6.5 cm gap. Filled squares: normal-flow 10 mTorr.
    PixelPoint(8, 10.0, 233.0, 220.0, "filled_square"),
    PixelPoint(8, 10.0, 274.0, 215.0, "filled_square"),
    PixelPoint(8, 10.0, 356.0, 214.0, "filled_square"),
    PixelPoint(8, 10.0, 520.0, 222.0, "filled_square"),
    PixelPoint(8, 10.0, 685.0, 275.0, "filled_square"),
    PixelPoint(8, 10.0, 788.0, 358.0, "filled_square"),
    PixelPoint(8, 10.0, 833.0, 322.0, "filled_square"),
    # Open square: explicitly one-third of the normal flow, diagnostic only.
    PixelPoint(8, 10.0, 520.0, 243.0, "open_square", "one_third_flow"),
    # Open circles: 2 mTorr.
    PixelPoint(8, 2.0, 397.0, 251.0, "open_circle"),
    PixelPoint(8, 2.0, 520.0, 288.0, "open_circle"),
    PixelPoint(8, 2.0, 685.0, 314.0, "open_circle"),
    PixelPoint(8, 2.0, 850.0, 378.0, "open_circle"),
    # Filled down-triangles: 1 mTorr.
    PixelPoint(8, 1.0, 233.0, 254.0, "filled_down_triangle"),
    PixelPoint(8, 1.0, 685.0, 396.0, "filled_down_triangle"),
    PixelPoint(8, 1.0, 809.0, 411.0, "filled_down_triangle"),
    # Open diamonds: 0.5 mTorr.
    PixelPoint(8, 0.5, 685.0, 433.0, "open_diamond"),
    PixelPoint(8, 0.5, 801.0, 441.0, "open_diamond"),
)


FIELDNAMES = (
    "source_figure",
    "window_to_wafer_gap_cm",
    "pressure_mTorr",
    "tcp_source_power_W",
    "relative_cl2_density_percent",
    "cl2_dissociation_percent",
    "cl2_flow_sccm",
    "rare_gas_flow_sccm",
    "flow_condition",
    "marker",
    "marker_center_x_px",
    "marker_center_y_px",
    "digitization_power_uncertainty_W",
    "digitization_relative_cl2_uncertainty_percentage_point",
    "reported_absolute_density_relative_uncertainty_percent",
    "error_bar_semantics",
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


def _flows(point: PixelPoint) -> tuple[float, float]:
    if point.flow_condition == "one_third_flow":
        return 38.0, 2.0
    return {
        10.0: (114.0, 6.0),
        2.0: (114.0, 6.0),
        1.0: (57.0, 3.0),
        0.5: (38.0, 2.0),
    }[point.pressure_mTorr]


def rows() -> list[dict[str, str]]:
    output = []
    for point in PIXEL_POINTS:
        axis = AXES[point.figure]
        power = axis.power_at_pixel(point.x_px)
        relative_cl2 = axis.relative_cl2_at_pixel(point.y_px)
        cl2_flow, rare_flow = _flows(point)
        power_uncertainty = (
            2.0 * axis.power_max_W
            / (axis.right_px - axis.left_px)
        )
        relative_cl2_uncertainty = (
            3.0 * axis.relative_cl2_max_percent
            / (axis.bottom_px - axis.top_px)
        )
        output.append({
            "source_figure": f"Figure {point.figure}",
            "window_to_wafer_gap_cm": f"{axis.gap_cm:g}",
            "pressure_mTorr": f"{point.pressure_mTorr:g}",
            "tcp_source_power_W": f"{power:.3f}",
            "relative_cl2_density_percent": f"{relative_cl2:.4f}",
            "cl2_dissociation_percent": f"{100.0 - relative_cl2:.4f}",
            "cl2_flow_sccm": f"{cl2_flow:g}",
            "rare_gas_flow_sccm": f"{rare_flow:g}",
            "flow_condition": point.flow_condition,
            "marker": point.marker,
            "marker_center_x_px": f"{point.x_px:.1f}",
            "marker_center_y_px": f"{point.y_px:.1f}",
            "digitization_power_uncertainty_W": f"{power_uncertainty:.4f}",
            "digitization_relative_cl2_uncertainty_percentage_point": (
                f"{relative_cl2_uncertainty:.4f}"
            ),
            "reported_absolute_density_relative_uncertainty_percent": "25",
            "error_bar_semantics": "range_between_Ar_and_Xe_reductions_not_sigma",
            "tcp_power_semantics": "power_into_matching_network_not_absorbed_power",
            "supports_absorbed_power": "false",
            "supports_wafer_flux": "false",
            "validation_role": (
                "diagnostic_flow_check"
                if point.flow_condition == "one_third_flow"
                else "reactor_dissociation_validation_candidate"
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
    powers = [float(row["tcp_source_power_W"]) for row in decoded]
    dissociations = [
        float(row["cl2_dissociation_percent"]) for row in decoded]
    return {
        "manifest_id": "MALYSHEV-1998-LAM-CL2-DISSOCIATION-R1",
        "source": {
            "citation": (
                "M. V. Malyshev, V. M. Donnelly, A. Kornblit, and "
                "N. A. Ciampa, Journal of Applied Physics 84, 137-146 "
                "(1998)"
            ),
            "doi": SOURCE_DOI,
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": SOURCE_PDF_PAGE,
            "figures": [7, 8],
            "native_images": {
                str(figure): {
                    "pdfimages_index_on_page": axis.native_image_index,
                    "size_px": list(axis.image_size),
                    "one_bit_pixel_sha256": axis.pixel_sha256,
                }
                for figure, axis in AXES.items()
            },
        },
        "experiment": {
            "tool": "Lam Research Alliance metal etcher, high-flow chamber",
            "plasma": "Cl2 plus 5 percent equal rare-gas mixture",
            "wafer": "SiO2-coated Si wafer",
            "substrate_rf_bias": "off",
            "wall": "anodized aluminum; chamber wall temperature 333 K",
            "source_power_measurement": "into matching network with inline meters",
            "gaps_cm": [6.5, 11.0],
            "pressures_mTorr": [0.5, 1.0, 2.0, 10.0],
        },
        "measurement": {
            "observable": "100 * n_Cl2 / n_Cl2_plasma_off",
            "derived_observable": "100 minus relative Cl2 percent",
            "method": "Cl2 305 nm emission normalized independently to Ar and Xe",
            "reported_absolute_density_relative_uncertainty_percent": 25.0,
            "error_bars": "ends are Ar-derived and Xe-derived values; not sigma",
        },
        "digitization": {
            "marker_count": len(PIXEL_POINTS),
            "method": (
                "native one-bit figure rasters extracted from publisher PDF; "
                "original-resolution human-vision marker reconciliation; "
                "PIL axis and point replay"
            ),
            "power_range_W": [min(powers), max(powers)],
            "dissociation_range_percent": [
                min(dissociations), max(dissociations)],
            "horizontal_pixel_allowance": 2.0,
            "vertical_pixel_allowance": 3.0,
            "csv_sha256": csv_sha256,
        },
        "exclusions": {
            "low_power_dense_overlap": (
                "markers not visually separable from overlapping traces were omitted"
            ),
            "documented_10_mTorr_approximately_100_W_anomaly": (
                "excluded; paper attributes apparent n_Cl2 above 100 percent "
                "to enhanced emission during a discharge-mode transition"
            ),
            "smooth_curves": "paper model outputs; never digitized as measurements",
        },
        "use_boundary": {
            "supports_reactor_dissociation_validation": True,
            "supports_absorbed_power_validation": False,
            "supports_wafer_flux_validation": False,
            "supports_feature_depth_validation": False,
            "primary_grade_observable": "relative_cl2_density_percent",
            "unphysical_derived_values": (
                "retain negative 100-minus-relative values; never clip the "
                "primary measured observable"
            ),
            "calibration_rule": (
                "feature depth and these measured markers may not tune the same "
                "reactor parameter set"
            ),
        },
    }


def _extract_native_images(source_pdf: Path, directory: Path) -> dict[int, Path]:
    if _sha256(source_pdf) != SOURCE_PDF_SHA256:
        raise ValueError("source PDF checksum does not match Malyshev 1998")
    prefix = directory / "malyshev_native"
    subprocess.run(
        [
            "pdfimages", "-f", str(SOURCE_PDF_PAGE), "-l",
            str(SOURCE_PDF_PAGE), "-png", str(source_pdf), str(prefix),
        ],
        check=True,
    )
    by_size = {
        Image.open(path).size: path
        for path in sorted(directory.glob("malyshev_native-*.png"))
    }
    images = {}
    for figure, axis in AXES.items():
        path = by_size.get(axis.image_size)
        if path is None:
            raise ValueError(f"native Figure {figure} raster was not recovered")
        if _pixel_sha256(path) != axis.pixel_sha256:
            raise ValueError(f"native Figure {figure} pixels do not match")
        images[figure] = path
    return images


def _write_overlays(images: dict[int, Path], overlay_directory: Path) -> None:
    overlay_directory.mkdir(parents=True, exist_ok=True)
    colors = {10.0: "red", 2.0: "blue", 1.0: "green", 0.5: "magenta"}
    for figure, source in images.items():
        image = Image.open(source).convert("RGB")
        draw = ImageDraw.Draw(image)
        for point in PIXEL_POINTS:
            if point.figure != figure:
                continue
            radius = 6
            draw.ellipse(
                (
                    point.x_px - radius, point.y_px - radius,
                    point.x_px + radius, point.y_px + radius,
                ),
                outline=colors[point.pressure_mTorr],
                width=2,
            )
        image.save(overlay_directory / f"malyshev_1998_figure{figure}_qa.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--overlay-directory", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    csv_payload = csv_text()
    csv_digest = hashlib.sha256(csv_payload.encode("utf-8")).hexdigest()
    manifest_payload = manifest(csv_digest)

    if args.source_pdf is not None:
        with tempfile.TemporaryDirectory() as temporary:
            images = _extract_native_images(
                args.source_pdf.resolve(), Path(temporary))
            if args.overlay_directory is not None:
                _write_overlays(images, args.overlay_directory.resolve())
    elif args.overlay_directory is not None:
        parser.error("--overlay-directory requires --source-pdf")

    if args.write:
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        CSV_PATH.write_text(csv_payload, encoding="utf-8")
        MANIFEST_PATH.write_text(
            json.dumps(manifest_payload, indent=2) + "\n", encoding="utf-8")
    else:
        print(csv_payload, end="")


if __name__ == "__main__":
    main()
