#!/usr/bin/env python3
"""Grade a deterministic homogeneous resonance-escape closure against Tian.

This is a no-fit model-discovery audit.  NIST line data, the printed chamber
geometry, ideal-gas absorber density, and four declared gas-temperature
sensitivities are passed through the fixed-quadrature Voigt/cylinder solver.
All 22 Ar trapping-factor markers are held out.  The purpose is to determine
whether a homogeneous complete-frequency-redistribution closure is sufficient,
not to tune it into agreement with HPEM.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean

from scipy.constants import atomic_mass, k as BOLTZMANN_CONSTANT_J_K

from petch.reactor_global import (
    CylindricalReactor,
    ResonanceLineData,
    deterministic_cylinder_resonance_escape,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "digitized_figures_5_10_5_12.csv"
)
OUTPUT_DIRECTORY = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_resonance_escape_audit"
)
JSON_OUTPUT = OUTPUT_DIRECTORY / "tian_2017_resonance_escape_audit.json"
REPORT_OUTPUT = OUTPUT_DIRECTORY / "TIAN_2017_RESONANCE_ESCAPE_AUDIT.md"

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


def audit() -> dict[str, object]:
    targets = _targets()
    rows = []
    grades = []
    for temperature in TEMPERATURE_SENSITIVITIES_K:
        observed_all = []
        predicted_all = []
        line_values: dict[str, dict[float, float]] = {
            name: {} for name in LINES}
        for series, line in LINES.items():
            for chlorine_percent, observed in targets[series].items():
                argon_fraction = 1.0 - chlorine_percent / 100.0
                absorber_density = (
                    argon_fraction * PRESSURE_PA
                    / (BOLTZMANN_CONSTANT_J_K * temperature)
                )
                receipt = deterministic_cylinder_resonance_escape(
                    GEOMETRY,
                    line,
                    absorber_density_m3=absorber_density,
                    gas_temperature_K=temperature,
                    geometry_quadrature_order=12,
                    frequency_quadrature_order=40,
                )
                prediction = receipt.trapping_factor
                line_values[series][chlorine_percent] = prediction
                observed_all.append(observed)
                predicted_all.append(prediction)
                rows.append({
                    "gas_temperature_K": temperature,
                    "series": series,
                    "cl2_fraction_percent": chlorine_percent,
                    "argon_ground_density_m3": absorber_density,
                    "observed_hpep_trapping_factor": observed,
                    "predicted_homogeneous_trapping_factor": prediction,
                    "signed_error_percent": 100.0 * (
                        prediction - observed) / observed,
                    "escape_probability": receipt.escape_probability,
                    "line_center_mean_path_optical_depth": (
                        receipt.line_center_mean_path_optical_depth),
                    "source_weighted_mean_exit_path_m": (
                        receipt.source_weighted_mean_exit_path_m),
                    "source_axial_scale_length_m": None,
                    "target_used_to_select_any_parameter": False,
                })
        line_grades = {}
        for series in LINES:
            fractions = sorted(targets[series])
            line_grades[series] = _mape(
                [targets[series][fraction] for fraction in fractions],
                [line_values[series][fraction] for fraction in fractions],
            )
        ratio_ordering_failures = []
        for fraction in sorted(targets["Ar_104.8_nm_trapping_factor"]):
            observed_ratio = (
                targets["Ar_104.8_nm_trapping_factor"][fraction]
                / targets["Ar_106.7_nm_trapping_factor"][fraction]
            )
            model_ratio = (
                line_values["Ar_104.8_nm_trapping_factor"][fraction]
                / line_values["Ar_106.7_nm_trapping_factor"][fraction]
            )
            if (observed_ratio - 1.0) * (model_ratio - 1.0) < 0.0:
                ratio_ordering_failures.append({
                    "cl2_fraction_percent": fraction,
                    "observed_104p8_over_106p7": observed_ratio,
                    "model_104p8_over_106p7": model_ratio,
                })
        grades.append({
            "gas_temperature_K": temperature,
            "combined_mape_percent": _mape(observed_all, predicted_all),
            "line_mape_percent": line_grades,
            "line_trapping_ordering_failure_count": len(
                ratio_ordering_failures),
            "line_trapping_ordering_failures": ratio_ordering_failures,
        })

    best = min(grades, key=lambda item: item["combined_mape_percent"])
    return {
        "schema": "petch.tian-2017-deterministic-resonance-escape-audit.v1",
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
        "atomic_data": {
            series: {
                "wavelength_nm": line.wavelength_nm,
                "transition_probability_s_inv": (
                    line.transition_probability_s_inv),
                "absorption_oscillator_strength": (
                    line.absorption_oscillator_strength),
                "source": line.source,
            }
            for series, line in LINES.items()
        },
        "model": {
            "method": (
                "fixed-quadrature 3-D cylinder paths plus normalized Voigt "
                "frequency integral and analytic far-wing tail"),
            "absorber": "uniform ideal-gas Ar ground state",
            "emitter": "uniform volume, line independent",
            "frequency_redistribution": "complete",
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
            "formal_gate_pass": False,
        },
        "rows": rows,
        "discovered_missing_state": [
            "spatial gas-temperature and ground-state-density field",
            "line-specific Ar(1s2) and Ar(1s4) emitter spatial moments",
            "partial frequency redistribution after resonant absorption",
            "state-specific collisional broadening and quenching",
        ],
        "verdict": (
            "A homogeneous complete-redistribution cylinder is not the next "
            "predictive rung. It misses the two-line trapping board and, over "
            "part of the mixture sweep, predicts the opposite ordering of "
            "104.8- and 106.7-nm trapping. Tian's own fields explain why: the "
            "gas spans 300--833 K and the two radiating states have different, "
            "mixture-dependent spatial support. The deterministic transport "
            "kernel is viable, but it must consume spatial moments from a "
            "two-zone or axisymmetric reactor rather than one homogeneous 0-D "
            "state. No depth parameter was touched."),
    }


def _report(result: dict[str, object]) -> str:
    summary = result["summary"]
    lines = [
        "# Tian 2017 deterministic resonance-escape audit",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        "All 22 digitized trapping markers were held out; no target selected "
        "a parameter. The best of the four declared homogeneous-temperature "
        f"sensitivities was {summary['best_declared_temperature_sensitivity_K']:.0f} K "
        f"at {summary['best_combined_mape_percent']:.2f}% combined MAPE. This "
        "is a failed reduced-model gate, not experimental validation.",
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
        "## Required next state",
        "",
        *[f"- {item}" for item in result["discovered_missing_state"]],
        "",
        "The appropriate deterministic replacement for Tian's photon Monte "
        "Carlo is fixed-quadrature ray/frequency transport over those spatial "
        "moments. A scalar escape factor or fitted broadband photon yield is "
        "not supported.",
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
