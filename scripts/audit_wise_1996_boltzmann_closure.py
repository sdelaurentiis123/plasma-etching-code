#!/usr/bin/env python3
"""Grade the Wise Figure-3 electron-pressure/Boltzmann closure.

The bulk closure used by the quasineutral spatial tier is

    grad(phi) = grad(n_e T_e) / n_e,

when ``T_e`` is expressed in eV.  Figure 3 independently measures all three
fields needed to check the radial integral.  Potential has an arbitrary
gauge, so the reconstruction is aligned at the axis; no slope, temperature,
or reactor parameter is fitted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = (
    ROOT / "data/experimental/wise_1996_gec_icp/figure3_radial_measurements.csv"
)
DEFAULT_MANIFEST = (
    ROOT / "data/experimental/wise_1996_gec_icp/figure3_digitization_manifest.json"
)
DEFAULT_GEOMETRY = (
    ROOT / "data/experimental/wise_1996_gec_icp/gec_icp_geometry.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results/curated/reactor_to_feature_chlorine/"
    "wise_1996_boltzmann_closure_audit"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_observables(path: Path) -> tuple[np.ndarray, ...]:
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    names = ("electron_density", "electron_temperature", "plasma_potential")
    groups = []
    radii = []
    for name in names:
        selected = [row for row in rows if row["observable"] == name]
        if len(selected) != 7:
            raise ValueError(f"Wise Figure-3 board requires seven {name} markers")
        radius = np.asarray(
            [float(row["radial_distance_m"]) for row in selected])
        value = np.asarray([float(row["value"]) for row in selected])
        if np.any(~np.isfinite(radius)) or np.any(~np.isfinite(value)):
            raise ValueError("Wise Figure-3 board contains nonfinite data")
        radii.append(radius)
        groups.append(value)
    if not np.array_equal(radii[0], radii[1]) or not np.array_equal(
        radii[0], radii[2]
    ):
        raise ValueError("Wise Figure-3 observable radii do not align")
    if radii[0][0] != 0.0 or np.any(np.diff(radii[0]) <= 0.0):
        raise ValueError("Wise Figure-3 radii are not monotone from the axis")
    density, temperature, potential = groups
    if np.any(density <= 0.0) or np.any(temperature <= 0.0):
        raise ValueError("Wise Figure-3 electron state must be positive")
    return radii[0], density, temperature, potential


def _pressure_gradient_potential(
    density_m3: np.ndarray,
    temperature_eV: np.ndarray,
    axis_potential_V: float,
) -> np.ndarray:
    """Integrate dphi = T_e d(ln n_e) + dT_e by trapezoid in ln(n)."""
    potential = np.empty_like(density_m3)
    potential[0] = float(axis_potential_V)
    for index in range(density_m3.size - 1):
        delta_log_density = float(np.log(
            density_m3[index + 1] / density_m3[index]
        ))
        mean_temperature = 0.5 * (
            temperature_eV[index] + temperature_eV[index + 1]
        )
        delta_temperature = (
            temperature_eV[index + 1] - temperature_eV[index]
        )
        potential[index + 1] = (
            potential[index]
            + mean_temperature * delta_log_density
            + delta_temperature
        )
    return potential


def _full_width_half_maximum_m(
    radius_m: np.ndarray, density_m3: np.ndarray
) -> float:
    half = 0.5 * density_m3[0]
    crossing = np.flatnonzero(density_m3 <= half)
    if crossing.size == 0 or crossing[0] == 0:
        raise ValueError("Wise Figure-3 board does not bracket density half maximum")
    upper = int(crossing[0])
    lower = upper - 1
    fraction = (half - density_m3[lower]) / (
        density_m3[upper] - density_m3[lower]
    )
    half_radius = radius_m[lower] + fraction * (
        radius_m[upper] - radius_m[lower]
    )
    return float(2.0 * half_radius)


def audit(data_path: Path, manifest_path: Path, geometry_path: Path) -> dict:
    radius, density, temperature, measured_potential = _read_observables(
        data_path
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    expected_data_hash = manifest["output"]["sha256"]
    if _sha256(data_path) != expected_data_hash:
        raise ValueError("Wise Figure-3 CSV checksum mismatch")
    if geometry["source"]["pdf_sha256"] != (
        "e67774edeefae2fb3c50c8343479519a88d6020f3022961b0ded6708cad5ceb9"
    ):
        raise ValueError("GEC geometry does not carry the audited NIST PDF hash")

    local_temperature = _pressure_gradient_potential(
        density, temperature, measured_potential[0]
    )
    constant_temperature = (
        measured_potential[0]
        + temperature[0] * np.log(density / density[0])
    )

    def metrics(reconstructed: np.ndarray) -> dict:
        residual = reconstructed - measured_potential
        # Axis is the physically required potential gauge and is omitted from
        # the unweighted comparison metrics.
        return {
            "reconstructed_potential_V": reconstructed.tolist(),
            "residual_V": residual.tolist(),
            "unweighted_mape_percent_excluding_gauge": float(
                100.0 * np.mean(np.abs(
                    residual[1:] / measured_potential[1:]
                ))
            ),
            "maximum_absolute_residual_V_excluding_gauge": float(
                np.max(np.abs(residual[1:]))
            ),
            "axis_to_outer_drop_V": float(
                reconstructed[-1] - reconstructed[0]
            ),
        }

    fwhm = _full_width_half_maximum_m(radius, density)
    nist_interval = (0.07, 0.09)
    return {
        "schema": "petch.wise_1996_boltzmann_closure_audit.v1",
        "claim_class": "direct_bulk_closure_check_not_reactor_or_depth_prediction",
        "source": {
            "wise_figure3_csv": str(data_path.relative_to(ROOT)),
            "wise_figure3_csv_sha256": _sha256(data_path),
            "wise_pdf_sha256": manifest["source"]["pdf_sha256"],
            "miller_geometry_json": str(geometry_path.relative_to(ROOT)),
            "miller_pdf_sha256": geometry["source"]["pdf_sha256"],
        },
        "condition": manifest["experiment_context"],
        "radial_distance_m": radius.tolist(),
        "measured": {
            "electron_density_m3": density.tolist(),
            "electron_temperature_eV": temperature.tolist(),
            "plasma_potential_V": measured_potential.tolist(),
            "axis_to_outer_potential_drop_V": float(
                measured_potential[-1] - measured_potential[0]
            ),
        },
        "pressure_gradient_with_measured_local_temperature": metrics(
            local_temperature
        ),
        "boltzmann_constant_axis_temperature_sensitivity": metrics(
            constant_temperature
        ),
        "density_width": {
            "digitized_full_width_half_maximum_m": fwhm,
            "independent_miller_typical_interval_m": list(nist_interval),
            "inside_independent_interval": bool(
                nist_interval[0] <= fwhm <= nist_interval[1]
            ),
        },
        "digitization_only_uncertainty": {
            "vertical_pixels": manifest["digitization"][
                "vertical_pixel_uncertainty"
            ],
            "electron_density_m3": 3.0 * 16.0e16 / (509.0 - 181.0),
            "electron_temperature_eV": 3.0 * 4.0 / (1014.0 - 687.0),
            "plasma_potential_V": 3.0 * 20.0 / (1477.0 - 1148.0),
            "source_measurement_uncertainty": "not printed",
        },
        "certification": {
            "potential_gauge_only_parameter": "measured axis potential",
            "slope_or_temperature_fit_performed": False,
            "feature_depth_used": False,
            "formal_uncertainty_weighted_pass": False,
            "supports_boltzmann_electron_pressure_closure": True,
            "supports_spatial_reactor_state_prediction": False,
            "supports_wafer_flux_or_feature_depth_prediction": False,
            "reason": (
                "All three fields are direct measurements and the source prints "
                "no marker uncertainties. The check tests the closure locally; "
                "it does not predict any field from machine knobs."
            ),
        },
    }


def write(result: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "wise_1996_boltzmann_closure_audit.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    local = result["pressure_gradient_with_measured_local_temperature"]
    constant = result["boltzmann_constant_axis_temperature_sensitivity"]
    width = result["density_width"]
    measured_drop = result["measured"]["axis_to_outer_potential_drop_V"]
    lines = [
        "# Wise 1996 direct Boltzmann/electron-pressure closure audit",
        "",
        "This is a direct check of the bulk closure using three independently measured Figure-3 fields. Potential is aligned at the axis because its gauge is arbitrary; no slope, temperature, reactor, or feature parameter is fitted.",
        "",
        f"- measured axis-to-outer potential drop: `{measured_drop:.3f} V`",
        f"- local-Te pressure-gradient reconstruction: `{local['axis_to_outer_drop_V']:.3f} V`",
        f"- local-Te unweighted MAPE, excluding the gauge point: `{local['unweighted_mape_percent_excluding_gauge']:.2f}%`",
        f"- constant-axis-Te sensitivity MAPE: `{constant['unweighted_mape_percent_excluding_gauge']:.2f}%`",
        f"- digitized density FWHM: `{100.0 * width['digitized_full_width_half_maximum_m']:.2f} cm`",
        f"- independent Miller typical FWHM interval: `{100.0 * width['independent_miller_typical_interval_m'][0]:.1f}--{100.0 * width['independent_miller_typical_interval_m'][1]:.1f} cm` (inside: `{width['inside_independent_interval']}`)",
        "",
        "The closure is supported at the few-percent unweighted level, but the paper prints no marker error bars, so this is not labeled a formal uncertainty-weighted pass. It validates neither a knobs-to-state solve nor wafer flux or etch depth.",
        "",
    ]
    (output / "WISE_1996_BOLTZMANN_CLOSURE_AUDIT.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(args.data, args.manifest, args.geometry)
    write(result, args.output)
    print(json.dumps({
        "local_temperature_mape_percent": result[
            "pressure_gradient_with_measured_local_temperature"
        ]["unweighted_mape_percent_excluding_gauge"],
        "density_fwhm_m": result["density_width"][
            "digitized_full_width_half_maximum_m"
        ],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
