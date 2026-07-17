#!/usr/bin/env python3
"""Digitize Krüger et al. 2024 Figure 4(a) into a weighted IEAD sample.

The source figure reports a normalized joint energy/angle density on a logarithmic color scale
from 1e-2 to 1. This script maps plot pixels to the colorbar, then performs deterministic
systematic resampling. It writes only derived numerical facts; the copyrighted source raster is
not redistributed.

The default pixel calibration is for the 1000x1253 RGB image embedded on PDF page 5 of the
author manuscript (PDF sha256 65b7750b...). The calibration and extraction sensitivities are
reported in the metadata and must travel with any boundary built from the table.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from PIL import Image

from petch.surface_kinetics import ReducedSiO2FluorocarbonParameters


EXPECTED_PDF_SHA256 = (
    "65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51")
EXPECTED_IMAGE_SHAPE = (1253, 1000, 3)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _quantile(value, weight, probability):
    order = np.argsort(value)
    selected = value[order]
    cumulative = np.cumsum(weight[order])
    cumulative /= cumulative[-1]
    return float(np.interp(float(probability), cumulative, selected))


def _nearest_legend_index(rgb, legend):
    output = []
    for start in range(0, len(rgb), 4096):
        selected = rgb[start:start + 4096]
        distance = np.sum(
            (selected[:, None, :] - legend[None, :, :]) ** 2, axis=2)
        output.append(np.argmin(distance, axis=1))
    return np.concatenate(output)


def _extract(
        image, *, colorfulness_threshold=15.0,
        energy_zero_y=1072.0, energy_5000_y=46.0,
        angle_zero_x=294.0, angle_plus5_x=384.0,
        legend_left_x=122, legend_right_x=444, legend_y=1205):
    y_grid, x_grid = np.mgrid[46:1073, 203:385]
    plot_rgb = image[46:1073, 203:385].astype(float)
    colorfulness = np.max(plot_rgb, axis=2) - np.min(plot_rgb, axis=2)
    selected = colorfulness > float(colorfulness_threshold)
    rgb = plot_rgb[selected]
    if len(rgb) < 1000:
        raise RuntimeError("IEAD color mask retained too few source pixels")
    legend = image[
        int(legend_y), int(legend_left_x):int(legend_right_x) + 1].astype(float)
    legend_index = _nearest_legend_index(rgb, legend)
    log10_density = -2.0 + 2.0 * legend_index / (len(legend) - 1)
    weight = np.power(10.0, log10_density)
    energy = (
        (float(energy_zero_y) - y_grid[selected])
        * 5000.0 / (float(energy_zero_y) - float(energy_5000_y)))
    angle = (
        (x_grid[selected] - float(angle_zero_x))
        * 5.0 / (float(angle_plus5_x) - float(angle_zero_x)))
    physical = (
        (energy >= 0.0) & (energy <= 5200.0)
        & (angle >= -6.5) & (angle <= 6.5)
        & np.isfinite(weight) & (weight > 0.0))
    energy = energy[physical]
    angle = angle[physical]
    weight = weight[physical]
    weight /= np.sum(weight)
    return energy, angle, weight


def _summary(energy, angle, weight):
    mean_energy = float(np.sum(weight * energy))
    mean_angle = float(np.sum(weight * angle))
    cosine = np.cos(np.deg2rad(angle))
    parameters = ReducedSiO2FluorocarbonParameters.krueger_2024_reduced_projection()
    return {
        "source_pixel_count": int(len(weight)),
        "mean_energy_eV": mean_energy,
        "energy_quantile_eV": {
            f"q{int(probability * 100):02d}": _quantile(
                energy, weight, probability)
            for probability in (0.01, 0.05, 0.50, 0.95, 0.99)
        },
        "mean_signed_angle_deg": mean_angle,
        "angle_standard_deviation_deg": float(np.sqrt(
            np.sum(weight * (angle - mean_angle) ** 2))),
        "mean_cosine_incidence": float(np.sum(weight * cosine)),
        "normal_and_oblique_mean_yield": {
            "bare_sio2": float(np.sum(
                weight * parameters.bare_sio2_yield.evaluate(energy, cosine))),
            "complex_sio2": float(np.sum(
                weight * parameters.complex_sio2_yield.evaluate(energy, cosine))),
            "polymer_sputter": float(np.sum(
                weight * parameters.polymer_sputter_yield.evaluate(energy, cosine))),
        },
    }


def _bin_joint(energy, angle, weight, energy_bin_eV, angle_bin_deg):
    """Compress pixels to weighted centroids while preserving energy-angle correlation."""
    energy_bin = np.floor(energy / float(energy_bin_eV)).astype(int)
    angle_bin = np.floor((angle + 10.0) / float(angle_bin_deg)).astype(int)
    key = energy_bin * 10000 + angle_bin
    rows = []
    for value in np.unique(key):
        selected = key == value
        total = float(np.sum(weight[selected]))
        rows.append((
            float(np.sum(energy[selected] * weight[selected]) / total),
            float(np.sum(angle[selected] * weight[selected]) / total),
            total,
        ))
    return tuple(sorted(rows))


def _uncertainty_summaries(image):
    summaries = []
    for threshold in (10.0, 15.0, 20.0):
        for pixel_shift in (-2.0, 0.0, 2.0):
            energy, angle, weight = _extract(
                image,
                colorfulness_threshold=threshold,
                energy_zero_y=1072.0 + pixel_shift,
                energy_5000_y=46.0 - pixel_shift,
                angle_zero_x=294.0 + 0.5 * pixel_shift,
                angle_plus5_x=384.0 - 0.5 * pixel_shift,
                legend_left_x=122 + int(pixel_shift),
                legend_right_x=444 - int(pixel_shift))
            summaries.append(_summary(energy, angle, weight))
    keys = (
        ("mean_energy_eV", lambda item: item["mean_energy_eV"]),
        ("q05_energy_eV", lambda item: item["energy_quantile_eV"]["q05"]),
        ("q50_energy_eV", lambda item: item["energy_quantile_eV"]["q50"]),
        ("q95_energy_eV", lambda item: item["energy_quantile_eV"]["q95"]),
        ("angle_standard_deviation_deg",
         lambda item: item["angle_standard_deviation_deg"]),
        ("bare_sio2_mean_yield",
         lambda item: item["normal_and_oblique_mean_yield"]["bare_sio2"]),
        ("complex_sio2_mean_yield",
         lambda item: item["normal_and_oblique_mean_yield"]["complex_sio2"]),
        ("polymer_sputter_mean_yield",
         lambda item: item["normal_and_oblique_mean_yield"]["polymer_sputter"]),
    )
    return {
        name: {
            "minimum": float(min(extractor(item) for item in summaries)),
            "maximum": float(max(extractor(item) for item in summaries)),
        }
        for name, extractor in keys
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--source-pdf", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--energy-bin-eV", type=float, default=100.0)
    parser.add_argument("--angle-bin-deg", type=float, default=0.25)
    args = parser.parse_args()
    if (not np.isfinite(args.energy_bin_eV) or args.energy_bin_eV <= 0.0
            or not np.isfinite(args.angle_bin_deg) or args.angle_bin_deg <= 0.0):
        raise ValueError("IEAD bin sizes must be positive and finite")
    source_sha = _sha(args.source_pdf)
    if source_sha != EXPECTED_PDF_SHA256:
        raise ValueError("Krüger source PDF checksum does not match the audited manuscript")
    image = np.asarray(Image.open(args.image).convert("RGB"))
    if image.shape != EXPECTED_IMAGE_SHAPE:
        raise ValueError(
            f"embedded Figure-4 image shape is {image.shape}, expected "
            f"{EXPECTED_IMAGE_SHAPE}")

    energy, angle, weight = _extract(image)
    records = _bin_joint(
        energy, angle, weight, args.energy_bin_eV, args.angle_bin_deg)
    sampled_energy = np.asarray([item[0] for item in records])
    sampled_angle = np.asarray([item[1] for item in records])
    sampled_weight = np.asarray([item[2] for item in records])
    source_summary = _summary(energy, angle, weight)
    sample_summary = _summary(sampled_energy, sampled_angle, sampled_weight)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("energy_eV", "signed_angle_deg", "probability_weight"))
        writer.writerows(records)
    metadata = {
        "source": (
            "Krüger et al., JVST A 42, 043008 (2024), Figure 4(a), "
            "https://doi.org/10.1116/6.0003554"),
        "source_pdf_sha256": source_sha,
        "embedded_image_sha256": _sha(args.image),
        "embedded_image_shape": list(image.shape),
        "digitization_method": (
            "nearest RGB point on the published log10 colorbar, followed by deterministic "
            "systematic weighted resampling"),
        "central_pixel_calibration": {
            "plot_region_x": [203, 384],
            "plot_region_y": [46, 1072],
            "energy_zero_y": 1072.0,
            "energy_5000_y": 46.0,
            "angle_zero_x": 294.0,
            "angle_plus5_x": 384.0,
            "colorfulness_threshold": 15.0,
            "legend_x": [122, 444],
            "legend_y": 1205,
            "legend_log10_density": [-2.0, 0.0],
        },
        "source_density_summary": source_summary,
        "resampled_node_count": len(records),
        "joint_binning": {
            "energy_bin_eV": float(args.energy_bin_eV),
            "angle_bin_deg": float(args.angle_bin_deg),
            "representative": "probability-weighted energy/angle centroid",
        },
        "resampled_summary": sample_summary,
        "resampling_error": {
            "mean_energy_relative": abs(
                sample_summary["mean_energy_eV"] - source_summary["mean_energy_eV"])
            / source_summary["mean_energy_eV"],
            "angle_standard_deviation_relative": abs(
                sample_summary["angle_standard_deviation_deg"]
                - source_summary["angle_standard_deviation_deg"])
            / source_summary["angle_standard_deviation_deg"],
            "bare_sio2_mean_yield_relative": abs(
                sample_summary["normal_and_oblique_mean_yield"]["bare_sio2"]
                - source_summary["normal_and_oblique_mean_yield"]["bare_sio2"])
            / source_summary["normal_and_oblique_mean_yield"]["bare_sio2"],
        },
        "digitization_sensitivity": _uncertainty_summaries(image),
        "limitations": [
            "the published density is normalized and clipped below 1e-2 of its maximum",
            "the combined positive-ion distribution has no species-resolved composition",
            "the signed-angle plot is a two-dimensional boundary; a 3-D azimuthal law "
            "requires a separately declared closure",
            "the HPEM distribution is reactor-model output, not a wafer measurement",
        ],
    }
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.output_metadata.write_text(
        json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
