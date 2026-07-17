#!/usr/bin/env python3
"""Digitize the held-out reactor inputs in Krüger et al. 2024 Figure 16.

Figure 16(a) reports HPEM wafer fluxes for the oxygen-ratio transfer cases. Figure 16(b)
reports normalized joint ion energy-angle distributions for the low-frequency-power transfer
cases.  These are reactor-model outputs, not experimental measurements.  They are nevertheless
the boundary inputs used by the paper's held-out MCFPM comparison and must be reconstructed before
petch can test the same transfer without inventing gas-flow or power scaling laws.

The default calibration is tied to the 1870x2475 PNG rendering of manuscript page 15 listed in
the output metadata.  Only derived numerical tables and provenance are written to the repository.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from PIL import Image


EXPECTED_PDF_SHA256 = (
    "65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51")
EXPECTED_PAGE_SHA256 = (
    "f944841595bbdcdfe1851753fc4028830470181afcbcc6ab888f139558ff55cd")
EXPECTED_PAGE_SHAPE = (2475, 1870, 3)
RATIOS = (0.5, 1.0, 1.5, 2.5)
POWERS_KW = (0.0, 4.0, 6.0, 8.0)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _quantile(value, weight, probability):
    order = np.argsort(value)
    value = value[order]
    cumulative = np.cumsum(weight[order])
    cumulative /= cumulative[-1]
    return float(np.interp(float(probability), cumulative, value))


def _nearest_legend_index(rgb, legend):
    result = []
    for start in range(0, len(rgb), 4096):
        selected = rgb[start:start + 4096]
        distance = np.sum(
            (selected[:, None, :] - legend[None, :, :]) ** 2, axis=2)
        result.append(np.argmin(distance, axis=1))
    return np.concatenate(result)


def _flux_color_masks(image):
    red, green, blue = np.moveaxis(image.astype(int), -1, 0)
    nearly_gray = np.max(image, axis=2) - np.min(image, axis=2) < 20
    return {
        "CF2": (green > red + 35) & (green > blue + 35) & (blue > 35),
        "C3F4": (blue > red + 70) & (red > 55) & (green < 80),
        "C2F3": (red > green + 20) & (red > blue + 20) & (red < 190),
        "CF": (red < 80) & (green < 80) & (blue < 80) & nearly_gray,
        "O": (np.abs(red - green) < 25) & (blue > red + 40) & (blue > 150),
        "ions": (red > green + 70) & (red > blue + 70) & (red > 150),
        "CF3": (blue > red + 45) & (blue > green + 45) & (blue < 200),
    }


def _extract_transfer_fluxes(image, pixel_half_width=10):
    # Full-page coordinates.  The three log-decade tick anchors are independently readable.
    x_left, x_right = 280.0, 882.0
    y_log_anchor = np.array([[407.0, 17.0], [606.0, 16.0], [804.0, 15.0]])
    log_slope, log_intercept = np.polyfit(
        y_log_anchor[:, 0], y_log_anchor[:, 1], 1)
    y_ranges = {
        "CF2": (350, 405),
        "C3F4": (385, 435),
        "C2F3": (405, 460),
        "CF": (455, 505),
        "O": (350, 500),
        "ions": (575, 625),
        "CF3": (580, 635),
    }
    masks = _flux_color_masks(image)
    rows = []
    extraction = []
    for ratio in RATIOS:
        nominal_x = x_left + (ratio - RATIOS[0]) / (RATIOS[-1] - RATIOS[0]) * (
            x_right - x_left)
        center_x = int(round(np.clip(nominal_x, x_left + 6, x_right - 6)))
        for species, mask in masks.items():
            low, high = y_ranges[species]
            selected_y, _ = np.where(mask[
                low:high + 1,
                center_x - int(pixel_half_width):center_x + int(pixel_half_width) + 1,
            ])
            selected_y = selected_y + low
            if selected_y.size < 4:
                raise RuntimeError(
                    f"Figure-16 flux trace for {species} at ratio {ratio} is unresolved")
            y = float(np.median(selected_y))
            log10_flux = float(log_slope * y + log_intercept)
            digitized = 10.0 ** log10_flux
            text_reported = {
                (0.5, "O"): 4.1e16,
                (2.5, "O"): 1.5e17,
            }.get((ratio, species))
            selected_flux = digitized if text_reported is None else text_reported
            # Two rendered pixels cover line thickness and calibration selection.  This is a
            # digitization interval, not uncertainty in the underlying HPEM calculation.
            lower = 10.0 ** (log_slope * (y + 2.0) + log_intercept)
            upper = 10.0 ** (log_slope * (y - 2.0) + log_intercept)
            if text_reported is not None:
                lower = min(lower, text_reported)
                upper = max(upper, text_reported)
            rows.append((
                ratio, species, selected_flux, digitized,
                "article_text" if text_reported is not None else "figure_16a_trace",
                lower, upper))
            extraction.append({
                "oxygen_to_fluorocarbon_ratio": ratio,
                "species": species,
                "nominal_x_pixel": nominal_x,
                "sample_center_x_pixel": center_x,
                "median_y_pixel": y,
                "selected_colored_pixel_count": int(selected_y.size),
            })
    return rows, extraction, {
        "x_ratio_0p5_pixel": x_left,
        "x_ratio_2p5_pixel": x_right,
        "log10_flux_tick_anchors": y_log_anchor.tolist(),
        "trace_half_width_pixels": int(pixel_half_width),
        "line_and_calibration_half_interval_pixels": 2.0,
    }


def _extract_iead_panel(image, panel_index, colorfulness_threshold=15.0):
    panel_x = ((279, 409), (440, 571), (600, 731), (761, 892))[panel_index]
    angle_zero_x = (344.0, 505.0, 665.0, 826.0)[panel_index]
    x_grid, y_grid = np.meshgrid(
        np.arange(panel_x[0], panel_x[1] + 1), np.arange(944, 1577))
    plot_rgb = image[944:1577, panel_x[0]:panel_x[1] + 1].astype(float)
    colorfulness = np.max(plot_rgb, axis=2) - np.min(plot_rgb, axis=2)
    selected = colorfulness > float(colorfulness_threshold)
    rgb = plot_rgb[selected]
    if rgb.shape[0] < 500:
        raise RuntimeError("Figure-16 IEAD color mask retained too few pixels")
    legend = image[1660, 486:723].astype(float)
    legend_index = _nearest_legend_index(rgb, legend)
    log10_density = -2.0 + 2.0 * legend_index / (legend.shape[0] - 1)
    weight = np.power(10.0, log10_density)
    energy = (1576.0 - y_grid[selected]) * 5000.0 / (1576.0 - 967.0)
    angle = (x_grid[selected] - angle_zero_x) * 5.0 / 55.0
    physical = (
        (energy >= 0.0) & (energy <= 5200.0)
        & (angle >= -6.0) & (angle <= 6.0)
        & np.isfinite(weight) & (weight > 0.0))
    energy, angle, weight = energy[physical], angle[physical], weight[physical]
    weight /= weight.sum()
    return energy, angle, weight


def _bin_iead(power, energy, angle, weight, energy_bin_eV, angle_bin_deg):
    energy_key = np.floor(energy / float(energy_bin_eV)).astype(int)
    angle_key = np.floor((angle + 10.0) / float(angle_bin_deg)).astype(int)
    key = energy_key * 10000 + angle_key
    rows = []
    for value in np.unique(key):
        selected = key == value
        total = float(weight[selected].sum())
        rows.append((
            float(power),
            float(np.sum(energy[selected] * weight[selected]) / total),
            float(np.sum(angle[selected] * weight[selected]) / total),
            total,
        ))
    return rows


def _iead_summary(energy, angle, weight):
    mean_angle = float(np.dot(weight, angle))
    return {
        "source_pixel_count": int(weight.size),
        "mean_energy_eV": float(np.dot(weight, energy)),
        "energy_quantile_eV": {
            name: _quantile(energy, weight, probability)
            for name, probability in (("q05", 0.05), ("q50", 0.50), ("q95", 0.95))
        },
        "mean_signed_angle_deg": mean_angle,
        "angle_standard_deviation_deg": float(np.sqrt(
            np.dot(weight, (angle - mean_angle) ** 2))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--page-image", type=Path, required=True)
    parser.add_argument("--source-pdf", type=Path, required=True)
    parser.add_argument("--output-flux-csv", type=Path, required=True)
    parser.add_argument("--output-iead-csv", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--energy-bin-eV", type=float, default=250.0)
    parser.add_argument("--angle-bin-deg", type=float, default=0.25)
    args = parser.parse_args()
    if _sha(args.source_pdf) != EXPECTED_PDF_SHA256:
        raise ValueError("Krüger source PDF checksum mismatch")
    if _sha(args.page_image) != EXPECTED_PAGE_SHA256:
        raise ValueError("Krüger rendered transfer page checksum mismatch")
    image = np.asarray(Image.open(args.page_image).convert("RGB"))
    if image.shape != EXPECTED_PAGE_SHAPE:
        raise ValueError("Krüger rendered transfer page shape mismatch")
    if (not np.isfinite(args.energy_bin_eV) or args.energy_bin_eV <= 0.0
            or not np.isfinite(args.angle_bin_deg) or args.angle_bin_deg <= 0.0):
        raise ValueError("IEAD bins must be finite and positive")

    flux_rows, flux_extraction, flux_calibration = _extract_transfer_fluxes(image)
    iead_rows = []
    iead_summaries = {}
    for index, power in enumerate(POWERS_KW):
        energy, angle, weight = _extract_iead_panel(image, index)
        iead_summaries[str(int(power))] = _iead_summary(energy, angle, weight)
        iead_rows.extend(_bin_iead(
            power, energy, angle, weight,
            args.energy_bin_eV, args.angle_bin_deg))

    args.output_flux_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_flux_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "oxygen_to_fluorocarbon_ratio", "species", "flux_cm2_s",
            "figure_digitized_flux_cm2_s", "selected_value_source",
            "digitization_lower_cm2_s", "digitization_upper_cm2_s"))
        writer.writerows(flux_rows)
    args.output_iead_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_iead_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "low_frequency_power_kw", "energy_eV", "signed_angle_deg",
            "probability_weight"))
        writer.writerows(iead_rows)

    metadata = {
        "source": (
            "Krüger et al., JVST A 42, 043008 (2024), Figure 16, "
            "https://doi.org/10.1116/6.0003554"),
        "evidence_kind": "published_HPEM_reactor_model_output",
        "source_pdf_sha256": EXPECTED_PDF_SHA256,
        "rendered_page_sha256": EXPECTED_PAGE_SHA256,
        "rendered_page_shape": list(image.shape),
        "flux_digitization": {
            "pixel_calibration": flux_calibration,
            "trace_records": flux_extraction,
            "units": "cm^-2 s^-1",
            "article_text_constraints": {
                "O_ratio_0p5_cm2_s": 4.1e16,
                "O_ratio_2p5_cm2_s": 1.5e17,
                "source_location": "Section IX.A immediately following Figure 16",
            },
        },
        "iead_digitization": {
            "panel_powers_kw": list(POWERS_KW),
            "plot_y_energy_zero_pixel": 1576.0,
            "plot_y_energy_5000eV_pixel": 967.0,
            "panel_angle_zero_pixels": [344.0, 505.0, 665.0, 826.0],
            "angle_five_degree_pixel_offset": 55.0,
            "legend_x_pixels": [486, 722],
            "legend_y_pixel": 1660,
            "legend_log10_density": [-2.0, 0.0],
            "colorfulness_threshold": 15.0,
            "joint_bins": {
                "energy_eV": float(args.energy_bin_eV),
                "angle_deg": float(args.angle_bin_deg),
            },
            "source_summaries": iead_summaries,
        },
        "limitations": [
            "all Figure-16 boundary values are HPEM predictions, not measured wafer fluxes",
            "normalized IEAD density below 1e-2 of panel maximum is clipped by the figure",
            "combined positive-ion composition is unresolved",
            "the plotted signed-angle plane does not determine a 3-D azimuthal law",
            "rendered line/color thickness defines digitization uncertainty but not HPEM model error",
        ],
    }
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.write_text(
        json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
