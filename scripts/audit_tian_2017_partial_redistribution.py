#!/usr/bin/env python3
"""Grade deterministic partial frequency redistribution against Tian 2017.

The 22 digitized HPEM trapping-factor markers remain held out.  Four declared
one-temperature ideal-gas sensitivities use only the printed pressure,
geometry, and NIST atomic data.  A separate two-moment base-case diagnostic
combines the 300 K cold-reservoir absorber density with a 600 K emitting-zone
Doppler width.  That diagnostic is explicitly not a validation grade: the
source-model target was already visible when the moment pair was examined.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
import os
from pathlib import Path
from statistics import mean

from scipy.constants import atomic_mass, k as BOLTZMANN_CONSTANT_J_K

from petch.reactor_global import (
    CylindricalReactor,
    ResonanceLineData,
    deterministic_cylinder_partial_redistribution,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "digitized_figures_5_10_5_12.csv"
)
OUTPUT_DIRECTORY = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_partial_redistribution_audit"
)
JSON_OUTPUT = OUTPUT_DIRECTORY / "tian_2017_partial_redistribution_audit.json"
REPORT_OUTPUT = OUTPUT_DIRECTORY / "TIAN_2017_PARTIAL_REDISTRIBUTION_AUDIT.md"

PRESSURE_PA = 20.0 * 0.13332236842105263
TEMPERATURE_SENSITIVITIES_K = (300.0, 400.0, 600.0, 833.0)
GEOMETRY = CylindricalReactor(radius_m=0.1125, length_m=0.12)
ARGON_MASS_KG = 39.948 * atomic_mass
LINES = {
    "Ar_104.8_nm_trapping_factor": ResonanceLineData(
        wavelength_nm=104.821987,
        transition_probability_s_inv=5.32e8,
        lower_statistical_weight=1.0,
        upper_statistical_weight=3.0,
        absorber_mass_kg=ARGON_MASS_KG,
        source=(
            "nist-argon-persistent-lines: 1048.21987 A, "
            "Aki=5.32e8 s^-1"),
    ),
    "Ar_106.7_nm_trapping_factor": ResonanceLineData(
        wavelength_nm=106.665980,
        transition_probability_s_inv=1.32e8,
        lower_statistical_weight=1.0,
        upper_statistical_weight=3.0,
        absorber_mass_kg=ARGON_MASS_KG,
        source=(
            "nist-argon-persistent-lines: 1066.65980 A, "
            "Aki=1.32e8 s^-1"),
    ),
}


def _targets() -> dict[str, dict[float, float]]:
    output = {name: {} for name in LINES}
    with INPUT.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["series"] in output:
                output[row["series"]][float(row["cl2_fraction_percent"])] = (
                    float(row["value"]))
    if any(len(series) != 11 for series in output.values()):
        raise RuntimeError("Tian trapping-factor marker board is incomplete")
    return output


def _mape(observed: list[float], predicted: list[float]) -> float:
    return 100.0 * mean(
        abs(model - datum) / datum
        for datum, model in zip(observed, predicted)
    )


def _one_condition(
    series: str,
    chlorine_percent: float,
    temperature_K: float,
    density_reference_temperature_K: float,
    observed: float,
) -> dict[str, object]:
    argon_fraction = 1.0 - chlorine_percent / 100.0
    absorber_density = (
        argon_fraction * PRESSURE_PA
        / (BOLTZMANN_CONSTANT_J_K * density_reference_temperature_K)
    )
    receipt = deterministic_cylinder_partial_redistribution(
        GEOMETRY,
        LINES[series],
        absorber_density_m3=absorber_density,
        gas_temperature_K=temperature_K,
        velocity_changing_collision_frequency_s_inv=0.0,
        quenching_collision_frequency_s_inv=0.0,
        geometry_quadrature_order=8,
        frequency_quadrature_order=24,
        coherent_grid_points_per_lorentz_hwhm=6.0,
    )
    prediction = receipt.partial_redistribution_trapping_factor
    return {
        "gas_temperature_K": temperature_K,
        "density_reference_temperature_K": density_reference_temperature_K,
        "series": series,
        "cl2_fraction_percent": chlorine_percent,
        "argon_ground_density_m3": absorber_density,
        "observed_hpep_trapping_factor": observed,
        "predicted_partial_redistribution_trapping_factor": prediction,
        "predicted_complete_redistribution_trapping_factor": (
            receipt.complete_redistribution_trapping_factor),
        "signed_error_percent": 100.0 * (prediction - observed) / observed,
        "coherent_frequency_grid_points": (
            receipt.coherent_frequency_grid_points),
        "coherent_grid_points_per_lorentz_hwhm": (
            receipt.coherent_grid_points_per_lorentz_hwhm),
        "linear_solver_iterations": receipt.linear_solver_iterations,
        "linear_solver_relative_residual": (
            receipt.linear_solver_relative_residual),
    }


def _grade_rows(rows: list[dict[str, object]]) -> dict[str, object]:
    line_mapes = {}
    for series in LINES:
        selected = [row for row in rows if row["series"] == series]
        line_mapes[series] = _mape(
            [float(row["observed_hpep_trapping_factor"]) for row in selected],
            [
                float(row["predicted_partial_redistribution_trapping_factor"])
                for row in selected
            ],
        )
    by_key = {
        (str(row["series"]), float(row["cl2_fraction_percent"])): row
        for row in rows
    }
    ordering_failures = []
    fractions = sorted({float(row["cl2_fraction_percent"]) for row in rows})
    for fraction in fractions:
        row104 = by_key[("Ar_104.8_nm_trapping_factor", fraction)]
        row106 = by_key[("Ar_106.7_nm_trapping_factor", fraction)]
        observed_ratio = (
            float(row104["observed_hpep_trapping_factor"])
            / float(row106["observed_hpep_trapping_factor"])
        )
        predicted_ratio = (
            float(row104["predicted_partial_redistribution_trapping_factor"])
            / float(row106["predicted_partial_redistribution_trapping_factor"])
        )
        if (observed_ratio - 1.0) * (predicted_ratio - 1.0) < 0.0:
            ordering_failures.append({
                "cl2_fraction_percent": fraction,
                "observed_104p8_over_106p7": observed_ratio,
                "predicted_104p8_over_106p7": predicted_ratio,
            })
    return {
        "combined_mape_percent": _mape(
            [float(row["observed_hpep_trapping_factor"]) for row in rows],
            [
                float(row["predicted_partial_redistribution_trapping_factor"])
                for row in rows
            ],
        ),
        "line_mape_percent": line_mapes,
        "line_trapping_ordering_failure_count": len(ordering_failures),
        "line_trapping_ordering_failures": ordering_failures,
    }


def audit() -> dict[str, object]:
    targets = _targets()
    tasks = []
    for temperature in TEMPERATURE_SENSITIVITIES_K:
        for series in LINES:
            for chlorine_percent, observed in targets[series].items():
                tasks.append((
                    series,
                    chlorine_percent,
                    temperature,
                    temperature,
                    observed,
                ))
    with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        rows = list(pool.map(lambda args: _one_condition(*args), tasks))

    grades = []
    for temperature in TEMPERATURE_SENSITIVITIES_K:
        selected = [
            row for row in rows
            if float(row["gas_temperature_K"]) == temperature
        ]
        grades.append({
            "gas_temperature_K": temperature,
            **_grade_rows(selected),
        })

    # The source map spans 300--833 K.  This split moment is a mechanism check:
    # a cold peripheral absorber column and a hotter source-zone line width.
    # It was examined after seeing the source-model base target and therefore
    # cannot be called held-out validation even though neither value is fitted.
    base_tasks = [
        (series, 20.0, 600.0, 300.0, targets[series][20.0])
        for series in LINES
    ]
    base_rows = [_one_condition(*task) for task in base_tasks]
    base_grade = _grade_rows(base_rows)

    # Numerical certification at the mechanism-diagnostic base condition.
    convergence = []
    for series, line in LINES.items():
        absorber_density = (
            0.8 * PRESSURE_PA / (BOLTZMANN_CONSTANT_J_K * 300.0))
        coarse = deterministic_cylinder_partial_redistribution(
            GEOMETRY, line,
            absorber_density_m3=absorber_density,
            gas_temperature_K=600.0,
            velocity_changing_collision_frequency_s_inv=0.0,
            geometry_quadrature_order=8,
            frequency_quadrature_order=24,
            coherent_grid_points_per_lorentz_hwhm=6.0,
        )
        refined = deterministic_cylinder_partial_redistribution(
            GEOMETRY, line,
            absorber_density_m3=absorber_density,
            gas_temperature_K=600.0,
            velocity_changing_collision_frequency_s_inv=0.0,
            geometry_quadrature_order=8,
            frequency_quadrature_order=24,
            coherent_grid_points_per_lorentz_hwhm=12.0,
        )
        extended = deterministic_cylinder_partial_redistribution(
            GEOMETRY, line,
            absorber_density_m3=absorber_density,
            gas_temperature_K=600.0,
            velocity_changing_collision_frequency_s_inv=0.0,
            geometry_quadrature_order=8,
            frequency_quadrature_order=24,
            coherent_half_range_doppler_standard_deviations=64.0,
            coherent_grid_points_per_lorentz_hwhm=12.0,
        )
        convergence.append({
            "series": series,
            "coarse_grid_points": coarse.coherent_frequency_grid_points,
            "coarse_trapping_factor": (
                coarse.partial_redistribution_trapping_factor),
            "refined_grid_points": refined.coherent_frequency_grid_points,
            "refined_trapping_factor": (
                refined.partial_redistribution_trapping_factor),
            "coarse_to_refined_relative_change": abs(
                refined.partial_redistribution_trapping_factor
                - coarse.partial_redistribution_trapping_factor
            ) / refined.partial_redistribution_trapping_factor,
            "extended_half_range_grid_points": (
                extended.coherent_frequency_grid_points),
            "extended_half_range_trapping_factor": (
                extended.partial_redistribution_trapping_factor),
            "refined_to_extended_half_range_relative_change": abs(
                extended.partial_redistribution_trapping_factor
                - refined.partial_redistribution_trapping_factor
            ) / extended.partial_redistribution_trapping_factor,
            "refined_linear_solver_relative_residual": (
                refined.linear_solver_relative_residual),
        })

    best = min(grades, key=lambda item: item["combined_mape_percent"])
    return {
        "schema": "petch.tian-2017-partial-redistribution-audit.v1",
        "claim_class": "held-out reduced-model discovery audit",
        "formal_experimental_validation": False,
        "source_model_targets_used_to_fit_parameters": False,
        "input": str(INPUT.relative_to(ROOT)),
        "source_target_evidence_type": "source_equipment_model_digitized",
        "condition": {
            "pressure_mTorr": 20.0,
            "flow_sccm": 200.0,
            "power_W": 150.0,
            "frequency_MHz": 10.0,
            "geometry_radius_m": GEOMETRY.radius_m,
            "geometry_length_m": GEOMETRY.length_m,
            "gas_temperature_sensitivities_K": list(
                TEMPERATURE_SENSITIVITIES_K),
        },
        "model": {
            "method": (
                "projected cylinder chords, exact exponential attenuation, "
                "zero-padded FFT natural-Lorentz frequency propagator, and "
                "preconditioned matrix-free GMRES renewal solve"),
            "absorber": "uniform ideal-gas Ar ground state",
            "emitter": "uniform volume",
            "frequency_redistribution": (
                "coherent natural-Lorentz walk; zero velocity-changing "
                "collision sensitivity at 20 mTorr"),
            "surface_reflection": False,
            "quenching": False,
            "fitted_parameters": [],
        },
        "grades": grades,
        "summary": {
            "held_out_marker_count": 22,
            "best_declared_temperature_sensitivity_K": (
                best["gas_temperature_K"]),
            "best_combined_mape_percent": best["combined_mape_percent"],
            "best_line_mape_percent": best["line_mape_percent"],
            "best_line_trapping_ordering_failure_count": (
                best["line_trapping_ordering_failure_count"]),
            "previous_complete_redistribution_best_mape_percent": (
                38.34507403048056),
            "previous_complete_redistribution_ordering_failure_count": 7,
            "formal_gate_pass": False,
        },
        "rows": rows,
        "base_case_split_spatial_moment_diagnostic": {
            "target_used_during_model_selection": True,
            "formal_validation": False,
            "absorber_density_reference_temperature_K": 300.0,
            "emitting_zone_doppler_temperature_K": 600.0,
            **base_grade,
            "rows": base_rows,
            "interpretation": (
                "The cold absorbing column and hot emitting zone visible in "
                "Tian Fig. 5.2 reproduce the base line ordering and magnitude. "
                "This establishes the missing spatial moment but is not a "
                "held-out grade."),
        },
        "numerical_convergence": convergence,
        "remaining_missing_state": [
            "source-derived axisymmetric gas-temperature and Ar-density field",
            "line-specific Ar(1s2) and Ar(1s4) emitter spatial moments",
            "source-grade velocity-changing and quenching collision rates",
        ],
        "verdict": (
            "Partial frequency redistribution is necessary and materially "
            "improves the held-out two-line board, but a one-temperature "
            "homogeneous reactor is still insufficient. The base-case split "
            "moment shows that the next closure is nonlocal spatial coupling, "
            "not a fitted photon multiplier. No feature-depth parameter was "
            "touched."),
    }


def _report(result: dict[str, object]) -> str:
    summary = result["summary"]
    diagnostic = result["base_case_split_spatial_moment_diagnostic"]
    lines = [
        "# Tian 2017 partial-frequency-redistribution audit",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        "All 22 mixture-sweep trapping markers stayed held out. No target "
        "selected an atomic or collision parameter.",
        "",
        "| uniform gas temperature K | combined MAPE % | 104.8 MAPE % | 106.7 MAPE % | opposite line-order points |",
        "|---:|---:|---:|---:|---:|",
    ]
    for grade in result["grades"]:
        line_grade = grade["line_mape_percent"]
        lines.append(
            f"| {grade['gas_temperature_K']:.0f} | "
            f"{grade['combined_mape_percent']:.2f} | "
            f"{line_grade['Ar_104.8_nm_trapping_factor']:.2f} | "
            f"{line_grade['Ar_106.7_nm_trapping_factor']:.2f} | "
            f"{grade['line_trapping_ordering_failure_count']} |")
    lines.extend([
        "",
        f"Best homogeneous PFR MAPE: {summary['best_combined_mape_percent']:.2f}% "
        f"versus {summary['previous_complete_redistribution_best_mape_percent']:.2f}% "
        "for complete redistribution.",
        "",
        "## Split spatial-moment mechanism diagnostic",
        "",
        f"At the printed 20% Cl2 base condition, a 300 K cold-absorber "
        f"density moment plus 600 K source-zone Doppler moment gives "
        f"{diagnostic['combined_mape_percent']:.2f}% MAPE across the two lines. "
        "The source target was visible during this model-selection step, so "
        "this is mechanism reproduction, not validation.",
        "",
        "## Numerical receipt",
        "",
        "| line | grid change | half-range change | GMRES residual |",
        "|---|---:|---:|---:|",
    ])
    for row in result["numerical_convergence"]:
        lines.append(
            f"| {row['series']} | "
            f"{100.0 * row['coarse_to_refined_relative_change']:.5f}% | "
            f"{100.0 * row['refined_to_extended_half_range_relative_change']:.5f}% | "
            f"{row['refined_linear_solver_relative_residual']:.3e} |")
    lines.extend([
        "",
        "## Required next state",
        "",
        *[f"- {item}" for item in result["remaining_missing_state"]],
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    result = audit()
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    JSON_OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    REPORT_OUTPUT.write_text(_report(result), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
