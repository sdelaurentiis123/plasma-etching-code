#!/usr/bin/env python3
"""Reproduce the quantitative Woo 2024 C4F6 reactor/profile board.

Woo's public Korea University thesis reports patterned-SiO2 and ACL etch
rates, plasma diagnostics, OES trends, and SEM profiles in one
CF4/C4F6/He ICP.  Figure 4.1 is digitized from an archived 600 dpi render by
fitting the black line segments between the five large markers.  Fitting the
segments avoids choosing a subjective center inside the square/triangle
markers and provides an independent check against the two endpoint values
printed in the body text.

This is a reactor/surface validation board, not a Krueger boundary transplant.
The thesis does not publish species-resolved ion fluxes, an IEAD, or absolute
neutral fluxes.  Its SEM exposure times were deliberately adjusted to obtain
similar depths, so the SEM panels are not blind absolute-depth targets.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PDF = ROOT / "research_sources" / "woo_2024_c4f6_thesis.pdf"
RENDER = (
    ROOT / "data" / "experimental" / "woo_2024_c4f6"
    / "figure4_1_page70_600dpi.png"
)
OUTPUT_DIR = ROOT / "data" / "experimental" / "woo_2024_c4f6"
CSV_PATH = OUTPUT_DIR / "figure4_1_patterned_etch_rates.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"
RESULT_DIR = ROOT / "results" / "curated" / "woo_2024_c4f6_board"
AUDIT_PATH = RESULT_DIR / "audit.json"
OVERLAY_PATH = RESULT_DIR / "figure4_1_overlay.png"

SOURCE_PDF_SHA256 = (
    "16bdf4843d0218fe6801ed418fcdd342ef2f825d29a7074d0a1e334b77229523"
)
RENDER_SHA256 = (
    "7b35d03cdd9edd15eacbce1bdac74b469c9756ce35c18fbc5c956493ef5f09a7"
)
RENDER_SIZE = (4300, 6072)

# Full-page 600 dpi axis centers.  Repeated horizontal/vertical black-pixel
# runs locate the labeled major ticks.  Two ticks on each axis are sufficient
# because both axes are linear.
X_AT_35_PERCENT = 1501.0
X_AT_40_PERCENT = 1751.0
Y_AT_60_NM_MIN = 2813.0
Y_AT_10_NM_MIN = 4179.0

GAS_FRACTIONS_PERCENT = np.asarray(
    [37.5, 43.75, 50.0, 56.25, 62.5], dtype=float
)
SERIES_BANDS = {
    "SiO2": (3050, 3650),
    "ACL": (3650, 4150),
}
TEXT_ENDPOINTS_NM_MIN = {
    "SiO2": (48.155, 31.922),
    "ACL": (25.12, 12.902),
}
FIELDNAMES = (
    "c4f6_fraction_of_cf4_plus_c4f6_percent",
    "material",
    "patterned_etch_rate_nm_min",
    "digitization_uncertainty_nm_min",
    "marker_center_x_px",
    "line_fit_center_y_px",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _x_at_fraction(fraction_percent: float) -> float:
    return X_AT_35_PERCENT + (
        (float(fraction_percent) - 35.0)
        * (X_AT_40_PERCENT - X_AT_35_PERCENT)
        / 5.0
    )


def _rate_at_y(y_px: float) -> float:
    return 60.0 + (
        (float(y_px) - Y_AT_60_NM_MIN)
        * (10.0 - 60.0)
        / (Y_AT_10_NM_MIN - Y_AT_60_NM_MIN)
    )


def _verify_sources() -> Image.Image:
    if _sha256(SOURCE_PDF) != SOURCE_PDF_SHA256:
        raise RuntimeError("Woo thesis PDF checksum does not match")
    if _sha256(RENDER) != RENDER_SHA256:
        raise RuntimeError("Woo Figure-4.1 render checksum does not match")
    image = Image.open(RENDER).convert("RGB")
    if image.size != RENDER_SIZE:
        raise RuntimeError(
            f"unexpected Woo render size {image.size}; expected {RENDER_SIZE}"
        )
    gray = np.asarray(image.convert("L"))
    axis_checks = {
        "left": np.mean(np.min(gray[2800:4240, 1297:1313], axis=1) < 64),
        "right": np.mean(np.min(gray[2800:4240, 3230:3246], axis=1) < 64),
        "top": np.mean(np.min(gray[2805:2821, 1297:3246], axis=0) < 64),
        "bottom": np.mean(np.min(gray[4221:4237, 1297:3246], axis=0) < 64),
    }
    failed = {
        name: float(value)
        for name, value in axis_checks.items()
        if value < 0.90
    }
    if failed:
        raise RuntimeError(f"Woo Figure-4.1 axis integrity failed: {failed}")
    return image


def _segment_fits(
    black: np.ndarray, x_points: np.ndarray, y_band: tuple[int, int]
) -> tuple[list[np.ndarray], list[dict[str, float]]]:
    """Fit the four connecting lines, excluding 55 px around each marker."""
    fits = []
    diagnostics = []
    y_low, y_high = y_band
    for left, right in zip(x_points[:-1], x_points[1:]):
        x_samples = []
        y_samples = []
        for x_px in range(int(left + 55), int(right - 55) + 1):
            candidates = np.flatnonzero(black[y_low:y_high, x_px]) + y_low
            if candidates.size:
                x_samples.append(float(x_px))
                y_samples.append(float(np.median(candidates)))
        if len(x_samples) < 180:
            raise RuntimeError("too few Woo line pixels for a stable segment fit")
        x_array = np.asarray(x_samples)
        y_array = np.asarray(y_samples)
        initial = np.polyfit(x_array, y_array, 1)
        residual = y_array - np.polyval(initial, x_array)
        centered = residual - np.median(residual)
        keep = np.abs(centered) < 6.0
        if np.count_nonzero(keep) < 180:
            raise RuntimeError("Woo robust line fit rejected too many columns")
        fit = np.polyfit(x_array[keep], y_array[keep], 1)
        mad = float(np.median(np.abs(centered[keep])))
        if mad > 1.5:
            raise RuntimeError(f"Woo line-fit residual is too large: {mad}")
        fits.append(fit)
        diagnostics.append(
            {
                "left_x_px": float(left),
                "right_x_px": float(right),
                "retained_column_count": int(np.count_nonzero(keep)),
                "median_absolute_residual_px": mad,
                "slope_px_per_px": float(fit[0]),
                "intercept_px": float(fit[1]),
            }
        )
    return fits, diagnostics


def _point_centers(fits: list[np.ndarray], x_points: np.ndarray) -> np.ndarray:
    centers = []
    for index, x_px in enumerate(x_points):
        adjacent = []
        if index:
            adjacent.append(float(np.polyval(fits[index - 1], x_px)))
        if index < len(fits):
            adjacent.append(float(np.polyval(fits[index], x_px)))
        centers.append(float(np.mean(adjacent)))
    return np.asarray(centers)


def digitization() -> tuple[list[dict[str, str]], dict[str, object], Image.Image]:
    image = _verify_sources()
    rgb = np.asarray(image)
    black = np.all(rgb < 100, axis=2)
    x_points = np.asarray(
        [_x_at_fraction(value) for value in GAS_FRACTIONS_PERCENT]
    )

    rows = []
    fit_diagnostics = {}
    digitized_endpoints = {}
    point_centers = {}
    for material, band in SERIES_BANDS.items():
        fits, diagnostics = _segment_fits(black, x_points, band)
        y_points = _point_centers(fits, x_points)
        rates = np.asarray([_rate_at_y(value) for value in y_points])
        fit_diagnostics[material] = diagnostics
        point_centers[material] = y_points.tolist()
        digitized_endpoints[material] = [
            float(rates[0]),
            float(rates[-1]),
        ]
        for fraction, x_px, y_px, rate in zip(
            GAS_FRACTIONS_PERCENT, x_points, y_points, rates
        ):
            rows.append(
                {
                    "c4f6_fraction_of_cf4_plus_c4f6_percent": f"{fraction:g}",
                    "material": material,
                    "patterned_etch_rate_nm_min": f"{rate:.4f}",
                    # 0.2 nm/min is wider than six vertical pixels and the
                    # largest body-text/line-fit endpoint discrepancy.
                    "digitization_uncertainty_nm_min": "0.2",
                    "marker_center_x_px": f"{x_px:.1f}",
                    "line_fit_center_y_px": f"{y_px:.3f}",
                }
            )

    endpoint_differences = {
        material: [
            digitized_endpoints[material][index] - reported[index]
            for index in range(2)
        ]
        for material, reported in TEXT_ENDPOINTS_NM_MIN.items()
    }
    maximum_endpoint_difference = max(
        abs(value)
        for differences in endpoint_differences.values()
        for value in differences
    )
    if maximum_endpoint_difference > 0.2:
        raise RuntimeError(
            "Woo pixel digitization does not reconcile to printed endpoints"
        )

    diagnostics = {
        "x_points_px": x_points.tolist(),
        "point_centers_y_px": point_centers,
        "segment_fits": fit_diagnostics,
        "digitized_endpoint_minus_text_nm_min": endpoint_differences,
        "maximum_abs_endpoint_difference_nm_min": maximum_endpoint_difference,
    }
    return rows, diagnostics, image


def _csv_text(rows: list[dict[str, str]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _manifest(
    rows: list[dict[str, str]],
    diagnostics: dict[str, object],
    csv_payload: str,
) -> dict[str, object]:
    return {
        "manifest_id": "WOO-2024-C4F6-FIG4.1-R1",
        "source": {
            "citation": (
                "Byungjun Woo, A Study on the Etching Characteristics of "
                "High Aspect Ratio Oxide Etching using C4F6 Plasma in "
                "Inductively Coupled Plasma with Low Frequency Bias Power, "
                "M.S. thesis, Korea University (2024)"
            ),
            "doi": "10.23186/korea.000000288984.11009.0001569",
            "repository_record": (
                "https://dcollection.korea.ac.kr/srch/srchDetail/000000288984"
            ),
            "local_pdf": "research_sources/woo_2024_c4f6_thesis.pdf",
            "pdf_sha256": SOURCE_PDF_SHA256,
            "pdf_page": 70,
            "print_page": 53,
            "figure": "Figure 4.1",
            "render_dpi": 600,
            "rendered_page_sha256": RENDER_SHA256,
            "rendered_page_size_px": list(RENDER_SIZE),
        },
        "experiment": {
            "reactor": "planar ICP",
            "pressure_mTorr": 10.0,
            "total_flow_sccm": 140.0,
            "he_flow_sccm": 100.0,
            "total_cf4_plus_c4f6_flow_sccm": 40.0,
            "source_power_W": 150.0,
            "bias_power_W": 400.0,
            "lower_electrode_temperature_C": 17.0,
            "sample": (
                "2 x 2 cm center coupon; 2400 nm SiO2, 1400 nm ACL, "
                "50 nm SiON; nominal 100 nm line / 500 nm pitch"
            ),
        },
        "digitization": {
            "method": (
                "600-dpi Poppler render; PIL/NumPy axis integrity; robust "
                "least-squares fits to black connecting segments with 55 px "
                "excluded around each large marker; adjacent segment "
                "intersection averaged at interior markers; visual overlay"
            ),
            "axis_calibration": {
                "x_at_35_percent_px": X_AT_35_PERCENT,
                "x_at_40_percent_px": X_AT_40_PERCENT,
                "y_at_60_nm_min_px": Y_AT_60_NM_MIN,
                "y_at_10_nm_min_px": Y_AT_10_NM_MIN,
            },
            "uncertainty_nm_min": 0.2,
            "diagnostics": diagnostics,
        },
        "body_text_cross_checks": {
            "figure4_1_endpoint_rates_nm_min": {
                material: list(values)
                for material, values in TEXT_ENDPOINTS_NM_MIN.items()
            },
            "figure4_3_ion_current_density_mA_cm2": [0.0168, 0.05255],
            "figure4_3_printed_percent_increase": 21.0,
            "figure4_3_arithmetic_percent_increase": (
                (0.05255 / 0.0168 - 1.0) * 100.0
            ),
            "power_sweep_nominal_fraction_percent": 56.25,
            "power_sweep_printed_flows_sccm": {"CF4": 15.0, "C4F6": 25.0},
            "power_sweep_fraction_from_printed_flows_percent": 62.5,
        },
        "output": {
            "path": (
                "data/experimental/woo_2024_c4f6/"
                "figure4_1_patterned_etch_rates.csv"
            ),
            "sha256": hashlib.sha256(csv_payload.encode("utf-8")).hexdigest(),
            "row_count": len(rows),
        },
        "claim_boundary": {
            "supports": [
                "absolute patterned SiO2 and ACL rate trends in one C4F6 ICP",
                "same-reactor ion-current, Te, self-bias, OES, XPS, and "
                "profile-shape conditioning for a future reactor model",
            ],
            "does_not_support": [
                "species-resolved ion flux or IEAD",
                "absolute neutral flux",
                "a boundary transplant to Krueger's CCP",
                "a blind feature-depth test; exposure time was adjusted to "
                "obtain 2100-2200 nm in the SEM comparison",
            ],
        },
    }


def _audit(manifest: dict[str, object]) -> dict[str, object]:
    cross = manifest["body_text_cross_checks"]
    return {
        "schema": "petch.woo-2024-c4f6-board.v1",
        "audit_id": "WOO-2024-C4F6-REACTOR-PROFILE-2026-08-06-R1",
        "source_pdf_sha256": SOURCE_PDF_SHA256,
        "figure4_1_render_sha256": RENDER_SHA256,
        "quantitative_rate_points": manifest["output"]["row_count"],
        "pixel_text_reconciliation": {
            "limit_nm_min": 0.2,
            "maximum_abs_difference_nm_min": manifest["digitization"][
                "diagnostics"
            ]["maximum_abs_endpoint_difference_nm_min"],
            "passed": True,
        },
        "source_internal_consistency": {
            "ion_current_printed_percent_increase": cross[
                "figure4_3_printed_percent_increase"
            ],
            "ion_current_arithmetic_percent_increase": cross[
                "figure4_3_arithmetic_percent_increase"
            ],
            "ion_current_percentage_consistent": False,
            "power_sweep_nominal_fraction_percent": cross[
                "power_sweep_nominal_fraction_percent"
            ],
            "power_sweep_fraction_from_printed_flows_percent": cross[
                "power_sweep_fraction_from_printed_flows_percent"
            ],
            "power_sweep_fraction_consistent": False,
        },
        "feature_depth_classification": {
            "reported_depth_range_nm": [2100.0, 2200.0],
            "exposure_time_adjusted_to_equalize_depth": True,
            "value_blind_held_out_depth": False,
            "may_calibrate_simulation_time_from_depth": False,
            "formal_absolute_feature_depth_pass": False,
        },
        "boundary_identifiability": {
            "aggregate_ion_current_density_measured": True,
            "electron_temperature_measured": True,
            "self_bias_measured": True,
            "species_resolved_ion_flux_measured": False,
            "iead_measured": False,
            "absolute_neutral_flux_measured": False,
            "oes_is_absolute_flux": False,
            "knobs_to_feature_depth_identified": False,
        },
        "verdict": (
            "Strong same-reactor C4F6 multi-observable board and an "
            "original-pixel absolute patterned-rate dataset. It is not a "
            "blind depth validation: the SEM times were selected to equalize "
            "depth, the kinetic boundary remains aggregate, and two printed "
            "condition/arithmetic inconsistencies must be resolved before a "
            "reactor model can receive a formal pass."
        ),
    }


def _json_text(value: dict[str, object]) -> str:
    return json.dumps(value, indent=2) + "\n"


def _draw_overlay(
    image: Image.Image,
    rows: list[dict[str, str]],
    output: Path,
) -> None:
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    colors = {"SiO2": "#00a651", "ACL": "#0066ff"}
    for row in rows:
        x_px = float(row["marker_center_x_px"])
        y_px = float(row["line_fit_center_y_px"])
        color = colors[row["material"]]
        radius = 22
        draw.ellipse(
            (x_px - radius, y_px - radius, x_px + radius, y_px + radius),
            outline=color,
            width=5,
        )
        draw.line(
            (x_px - 14, y_px, x_px + 14, y_px), fill=color, width=4
        )
        draw.line(
            (x_px, y_px - 14, x_px, y_px + 14), fill=color, width=4
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    overlay.crop((1180, 2700, 3370, 4450)).save(output)


def _payloads() -> tuple[str, str, str, list[dict[str, str]], Image.Image]:
    rows, diagnostics, image = digitization()
    csv_payload = _csv_text(rows)
    manifest = _manifest(rows, diagnostics, csv_payload)
    audit = _audit(manifest)
    return (
        csv_payload,
        _json_text(manifest),
        _json_text(audit),
        rows,
        image,
    )


def write_outputs() -> None:
    csv_payload, manifest_payload, audit_payload, rows, image = _payloads()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    CSV_PATH.write_text(csv_payload, encoding="utf-8")
    MANIFEST_PATH.write_text(manifest_payload, encoding="utf-8")
    AUDIT_PATH.write_text(audit_payload, encoding="utf-8")
    _draw_overlay(image, rows, OVERLAY_PATH)


def verify_outputs() -> None:
    csv_payload, manifest_payload, audit_payload, _, _ = _payloads()
    expected = {
        CSV_PATH: csv_payload,
        MANIFEST_PATH: manifest_payload,
        AUDIT_PATH: audit_payload,
    }
    for path, payload in expected.items():
        if path.read_text(encoding="utf-8") != payload:
            raise RuntimeError(f"committed Woo output is stale: {path}")
    if not OVERLAY_PATH.is_file():
        raise RuntimeError("committed Woo visual overlay is missing")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        write_outputs()
    else:
        verify_outputs()
    print(AUDIT_PATH.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
