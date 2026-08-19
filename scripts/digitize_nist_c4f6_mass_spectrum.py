#!/usr/bin/env python3
"""Vision-audit the NIST 70 eV electron-ionization spectrum of 1,3-C4F6.

NIST WebBook SRD 69 exposes a printable 800 x 600 PNG but explicitly forbids
redistributing the spectrum.  The source pixels therefore remain local.  This
script checksum-locks that render, verifies every red stick against the source
image, and reproduces only the numerical digitization and its audit manifest.

The spectrum is a molecular-beam fragmentation prior, not a plasma ion
mixture.  In particular, it cannot replace Benck's absolute reactor currents
or identify Krueger's unpublished species-resolved boundary.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRECTORY = ROOT / "data" / "experimental" / "nist_c4f6_mass_spectrum"
CSV_PATH = OUTPUT_DIRECTORY / "electron_ionization_sticks.csv"
MANIFEST_PATH = OUTPUT_DIRECTORY / "digitization_manifest.json"

SOURCE_URL = (
    "https://webbook.nist.gov/cgi/cbook.cgi?"
    "Index=0&Large=on&Spec=C685632&Type=Mass"
)
SOURCE_IMAGE_SHA256 = (
    "3efe05fb3a4c27cbc74e46a6b342cce9010a17ff4a9f645b635c4ea04c96e011"
)
SOURCE_IMAGE_SIZE = (800, 600)
NIST_MS_NUMBER = 5987

# Full-image plot calibration.  The red base peak at nominal m/z 93 supplies
# the 100-intensity top and the printed plot frame supplies 0/200 and 0/100.
PLOT_LEFT_X_PX = 102.0
PLOT_RIGHT_X_PX = 743.0
PLOT_TOP_Y_PX = 66.0
PLOT_BASELINE_Y_PX = 506.0


@dataclass(frozen=True)
class Stick:
    nominal_m_over_z: int
    x_px: float
    top_y_px: float
    assignment: str
    assignment_class: str


STICKS = (
    Stick(31, 202.0, 439.0, "CF+", "monoisotopic_fragment"),
    Stick(32, 205.0, 504.0, "CF+ isotope", "isotope_satellite"),
    Stick(36, 218.0, 504.0, "C3+", "monoisotopic_fragment"),
    Stick(43, 240.5, 503.0, "C2F+", "monoisotopic_fragment"),
    Stick(50, 262.5, 499.0, "CF2+", "monoisotopic_fragment"),
    Stick(55, 279.0, 494.0, "C3F+", "monoisotopic_fragment"),
    Stick(62, 301.5, 494.0, "C2F2+", "monoisotopic_fragment"),
    Stick(69, 323.5, 485.0, "CF3+", "monoisotopic_fragment"),
    Stick(74, 339.5, 474.0, "C3F2+", "monoisotopic_fragment"),
    Stick(75, 343.0, 498.0, "C3F2+ isotope", "isotope_satellite"),
    Stick(81, 362.0, 497.0, "C2F3+", "monoisotopic_fragment"),
    Stick(93, 400.5, 66.0, "C3F3+", "monoisotopic_fragment"),
    Stick(94, 403.5, 490.0, "C3F3+ isotope", "isotope_satellite"),
    Stick(112, 461.5, 454.0, "C3F4+", "monoisotopic_fragment"),
    Stick(113, 465.0, 504.0, "C3F4+ isotope", "isotope_satellite"),
    Stick(124, 500.0, 488.0, "C4F4+", "monoisotopic_fragment"),
    Stick(125, 503.0, 504.0, "C4F4+ isotope 1", "isotope_satellite"),
    Stick(126, 506.0, 504.0, "C4F4+ isotope 2", "isotope_satellite"),
    Stick(131, 522.5, 498.0, "C3F5+", "monoisotopic_fragment"),
    Stick(143, 560.5, 482.0, "C4F5+", "monoisotopic_fragment"),
    Stick(144, 564.0, 502.0, "C4F5+ isotope", "isotope_satellite"),
    Stick(162, 621.5, 313.0, "C4F6+", "monoisotopic_parent"),
    Stick(163, 624.5, 497.0, "C4F6+ isotope", "isotope_satellite"),
)

FIELDNAMES = (
    "nominal_m_over_z",
    "assignment",
    "assignment_class",
    "relative_intensity_percent",
    "stick_center_x_full_image_px",
    "stick_top_y_full_image_px",
    "m_over_z_from_axis",
    "digitization_absolute_intensity_bound_percent",
)


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def m_over_z_at_pixel(x_px: float) -> float:
    return (
        (float(x_px) - PLOT_LEFT_X_PX)
        * 200.0
        / (PLOT_RIGHT_X_PX - PLOT_LEFT_X_PX)
    )


def intensity_at_pixel(y_px: float) -> float:
    return (
        (PLOT_BASELINE_Y_PX - float(y_px))
        * 100.0
        / (PLOT_BASELINE_Y_PX - PLOT_TOP_Y_PX)
    )


def rows() -> list[dict[str, str]]:
    return [
        {
            "nominal_m_over_z": str(stick.nominal_m_over_z),
            "assignment": stick.assignment,
            "assignment_class": stick.assignment_class,
            "relative_intensity_percent": f"{intensity_at_pixel(stick.top_y_px):.6f}",
            "stick_center_x_full_image_px": f"{stick.x_px:.1f}",
            "stick_top_y_full_image_px": f"{stick.top_y_px:.1f}",
            "m_over_z_from_axis": f"{m_over_z_at_pixel(stick.x_px):.6f}",
            # One vertical pixel is 0.2273 intensity point.  Two pixels cover
            # antialiasing and manual top-row selection conservatively.
            "digitization_absolute_intensity_bound_percent": "0.455",
        }
        for stick in STICKS
    ]


def csv_text() -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows())
    return stream.getvalue()


def _channel_intensity(assignment: str) -> float:
    return next(
        intensity_at_pixel(stick.top_y_px)
        for stick in STICKS
        if stick.assignment == assignment
    )


def manifest(csv_digest: str) -> dict[str, object]:
    cf = _channel_intensity("CF+")
    cf2 = _channel_intensity("CF2+")
    cf3 = _channel_intensity("CF3+")
    monoisotopic_total = sum(
        intensity_at_pixel(stick.top_y_px)
        for stick in STICKS
        if stick.assignment_class.startswith("monoisotopic")
    )
    return {
        "manifest_id": "NIST-WEBBOOK-C4F6-EI-MASS-SPECTRUM-R1",
        "source": {
            "database": "NIST Chemistry WebBook SRD 69",
            "species": "1,3-butadiene, 1,1,2,3,4,4-hexafluoro-",
            "formula": "C4F6",
            "cas_registry_number": "685-63-2",
            "nist_ms_number": NIST_MS_NUMBER,
            "spectrum_type": "electron ionization mass spectrum",
            "url": SOURCE_URL,
            "source_image_sha256": SOURCE_IMAGE_SHA256,
            "source_image_size_px": list(SOURCE_IMAGE_SIZE),
            "redistribution": (
                "NIST states that this spectrum cannot be downloaded; source "
                "pixels and overlays are not committed"
            ),
        },
        "pixel_calibration": {
            "plot_left_x_px_at_m_over_z_0": PLOT_LEFT_X_PX,
            "plot_right_x_px_at_m_over_z_200": PLOT_RIGHT_X_PX,
            "plot_top_y_px_at_relative_intensity_100": PLOT_TOP_Y_PX,
            "plot_baseline_y_px_at_relative_intensity_0": PLOT_BASELINE_Y_PX,
            "maximum_nominal_mass_axis_error": max(
                abs(m_over_z_at_pixel(stick.x_px) - stick.nominal_m_over_z)
                for stick in STICKS
            ),
            "absolute_intensity_bound_percent": 0.455,
        },
        "digitization": {
            "method": (
                "PIL full-image red-stick segmentation followed by original-"
                "resolution visual assignment and two-axis affine replay"
            ),
            "point_count": len(STICKS),
            "visual_audit_status": "passed_original_resolution",
            "raw_source_committed": False,
        },
        "derived_checks": {
            "base_peak": "C3F3+",
            "parent_relative_intensity_percent": _channel_intensity("C4F6+"),
            "cfx_direct_ei_relative_intensity": {
                "CF+": cf,
                "CF2+": cf2,
                "CF3+": cf3,
            },
            "cfx_direct_ei_ratios_to_cf": {
                "CF2+/CF+": cf2 / cf,
                "CF3+/CF+": cf3 / cf,
            },
            "heavy_parent_or_c3f3_fraction_of_monoisotopic_intensity": (
                (_channel_intensity("C3F3+") + _channel_intensity("C4F6+"))
                / monoisotopic_total
            ),
        },
        "model_consequence": {
            "required": [
                "retain C3F3+ and C4F6+ direct-parent ionization products",
                "resolve secondary fragment ionization and ion-neutral chemistry",
                "grade reactor composition against Benck mixture and pressure boards",
            ],
            "forbidden": [
                "map aggregate C4F6 ionization directly only to CF+/CF2+/CF3+",
                "use EI relative intensity as an absolute plasma ion flux",
                "transplant the 70 eV spectrum into Krueger as a wafer mixture",
                "fit a C4F6 branching ratio to Krueger's 825 nm endpoint",
            ],
        },
        "output": {
            "path": (
                "data/experimental/nist_c4f6_mass_spectrum/"
                "electron_ionization_sticks.csv"
            ),
            "sha256": csv_digest,
        },
    }


def manifest_text(csv_payload: str) -> str:
    digest = sha256(csv_payload.encode("utf-8")).hexdigest()
    return json.dumps(manifest(digest), indent=2) + "\n"


def verify_committed_files() -> None:
    payload = csv_text()
    if CSV_PATH.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed NIST C4F6 mass-spectrum CSV is stale")
    if MANIFEST_PATH.read_text(encoding="utf-8") != manifest_text(payload):
        raise RuntimeError("committed NIST C4F6 mass-spectrum manifest is stale")


def verify_source_image(path: Path) -> Image.Image:
    if _sha256(path) != SOURCE_IMAGE_SHA256:
        raise RuntimeError("NIST C4F6 spectrum image checksum changed")
    image = Image.open(path).convert("RGB")
    if image.size != SOURCE_IMAGE_SIZE:
        raise RuntimeError(f"unexpected NIST spectrum size: {image.size}")
    pixels = np.asarray(image)
    red = (
        (pixels[:, :, 0] > 220)
        & (pixels[:, :, 1] < 180)
        & (pixels[:, :, 2] < 180)
    )
    for stick in STICKS:
        x = int(round(stick.x_px))
        top = int(round(stick.top_y_px))
        baseline = int(round(PLOT_BASELINE_Y_PX))
        support = int(np.sum(red[max(top - 1, 0):baseline + 1, x - 1:x + 2]))
        expected = max(baseline - top, 1)
        if support < 0.75 * expected:
            raise RuntimeError(
                f"insufficient red-stick support at m/z {stick.nominal_m_over_z}: "
                f"{support}"
            )
    return image


def draw_overlay(image: Image.Image, output: Path) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for stick in STICKS:
        x = stick.x_px
        y = stick.top_y_px
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), outline="#0066cc", width=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-image", type=Path)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()

    verify_committed_files()
    image = (
        verify_source_image(args.source_image)
        if args.source_image is not None else None
    )
    if args.overlay is not None:
        if image is None:
            raise ValueError("--overlay requires --source-image")
        draw_overlay(image, args.overlay)
    print(json.dumps({
        "status": "verified",
        "point_count": len(STICKS),
        "source_image_verified": image is not None,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
