#!/usr/bin/env python3
"""Audit coupled zonal/partial-frequency transport at Tian's base case."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.constants import atomic_mass, k as BOLTZMANN_CONSTANT_J_K

from petch.reactor_global import (
    AxisymmetricRadiationZoneField,
    ResonanceLineData,
    deterministic_zonal_partial_redistribution,
)


ROOT = Path(__file__).resolve().parents[1]
FIELD_INPUT = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "base_case_spatial_moments.json"
)
TARGET_INPUT = (
    ROOT / "data" / "experimental" / "tian_2017_vuv_figures"
    / "digitized_figures_5_10_5_12.csv"
)
OUTPUT_DIRECTORY = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "tian_2017_zonal_partial_redistribution_audit"
)
JSON_OUTPUT = (
    OUTPUT_DIRECTORY / "tian_2017_zonal_partial_redistribution_audit.json")
REPORT_OUTPUT = (
    OUTPUT_DIRECTORY / "TIAN_2017_ZONAL_PARTIAL_REDISTRIBUTION_AUDIT.md")
ARGON_MASS_KG = 39.948 * atomic_mass
LINES = {
    "Ar_104.8_nm_trapping_factor": ResonanceLineData(
        wavelength_nm=104.821987,
        transition_probability_s_inv=5.32e8,
        lower_statistical_weight=1.0,
        upper_statistical_weight=3.0,
        absorber_mass_kg=ARGON_MASS_KG,
        source="nist-argon-persistent-lines 1048.21987 A",
    ),
    "Ar_106.7_nm_trapping_factor": ResonanceLineData(
        wavelength_nm=106.665980,
        transition_probability_s_inv=1.32e8,
        lower_statistical_weight=1.0,
        upper_statistical_weight=3.0,
        absorber_mass_kg=ARGON_MASS_KG,
        source="nist-argon-persistent-lines 1066.65980 A",
    ),
}


def _field() -> tuple[AxisymmetricRadiationZoneField, dict[str, object]]:
    source = json.loads(FIELD_INPUT.read_text(encoding="utf-8"))
    moment = source["axisymmetric_zone_field"]
    condition = source["condition"]
    temperatures = np.asarray([
        zone["gas_temperature_K"] for zone in moment["zones"]
    ])
    pressure_pa = condition["pressure_mTorr"] * 0.13332236842105263
    absorber = (
        condition["Ar_fraction"] * pressure_pa
        / (BOLTZMANN_CONSTANT_J_K * temperatures))
    emitter = np.asarray([
        zone["Ar_1s4_emitter_density_cm3"] for zone in moment["zones"]
    ]) * 1.0e6
    return AxisymmetricRadiationZoneField(
        radial_edges_m=np.asarray(moment["radial_edges_cm"]) * 1.0e-2,
        axial_edges_m=np.asarray(moment["axial_edges_cm"]) * 1.0e-2,
        cell_zone_index=np.asarray(
            moment["cell_zone_index_radial_by_axial"]),
        gas_temperature_K=temperatures,
        absorber_density_m3=absorber,
        emitter_density_m3=emitter,
        source=(
            "Tian 2017 Figs. 5.2(b), 5.3(b) contour-bracketed "
            "three-zone moment field"),
    ), source


def _targets() -> dict[str, float]:
    target = {}
    with TARGET_INPUT.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if (
                row["series"] in LINES
                and float(row["cl2_fraction_percent"]) == 20.0
            ):
                target[row["series"]] = float(row["value"])
    if set(target) != set(LINES):
        raise RuntimeError("Tian base trapping targets are incomplete")
    return target


def _run(field, line, order, points_per_hwhm):
    return deterministic_zonal_partial_redistribution(
        field,
        line,
        velocity_changing_collision_frequency_s_inv=0.0,
        quenching_collision_frequency_s_inv=0.0,
        surface_quadrature_order=order,
        direction_quadrature_order=order,
        frequency_quadrature_order=24,
        coherent_grid_points_per_lorentz_hwhm=points_per_hwhm,
    )


def audit() -> dict[str, object]:
    field, source = _field()
    targets = _targets()
    rows = []
    convergence = []
    for series, line in LINES.items():
        by_order = {
            order: _run(field, line, order, 6.0)
            for order in (8, 10, 12)
        }
        frequency_refined = _run(field, line, 12, 12.0)
        primary = by_order[12]
        observed = targets[series]
        rows.append({
            "series": series,
            "observed_hpep_trapping_factor": observed,
            "predicted_zonal_partial_redistribution_trapping_factor": (
                primary.trapping_factor),
            "predicted_zonal_complete_redistribution_trapping_factor": (
                primary.complete_frequency_redistribution_trapping_factor),
            "signed_error_percent": 100.0 * (
                primary.trapping_factor - observed) / observed,
            "initial_emission_zone_probability": (
                primary.initial_emission_zone_probability.tolist()),
            "escape_boundary_labels": list(primary.escape_boundary_labels),
            "partial_redistribution_escape_boundary_probability": (
                primary.partial_redistribution_escape_boundary_probability
                .tolist()),
            "complete_redistribution_escape_boundary_probability": (
                primary.complete_redistribution_escape_boundary_probability
                .tolist()),
            "partial_redistribution_wafer_escape_probability": (
                primary.partial_redistribution_wafer_escape_probability),
            "partial_redistribution_quench_probability": (
                primary.partial_redistribution_quench_probability),
            "terminal_probability_conservation_error_maximum": (
                primary.terminal_probability_conservation_error_maximum),
            "transition_probability_conservation_error_maximum": (
                primary.transition_probability_conservation_error_maximum),
            "zone_source_measure_relative_volume_error_maximum": (
                primary.zone_source_measure_relative_volume_error_maximum),
            "linear_solver_iterations": primary.linear_solver_iterations,
            "linear_solver_relative_residual": (
                primary.linear_solver_relative_residual),
            "coherent_frequency_grid_points": (
                primary.coherent_frequency_grid_points),
        })
        convergence.append({
            "series": series,
            "order_8_trapping_factor": by_order[8].trapping_factor,
            "order_10_trapping_factor": by_order[10].trapping_factor,
            "order_12_trapping_factor": by_order[12].trapping_factor,
            "order_10_to_12_relative_change": abs(
                by_order[12].trapping_factor - by_order[10].trapping_factor
            ) / by_order[12].trapping_factor,
            "order_10_to_12_wafer_escape_relative_change": abs(
                by_order[12].partial_redistribution_wafer_escape_probability
                - by_order[10].partial_redistribution_wafer_escape_probability
            ) / by_order[12].partial_redistribution_wafer_escape_probability,
            "frequency_coarse_grid_points": (
                by_order[12].coherent_frequency_grid_points),
            "frequency_refined_grid_points": (
                frequency_refined.coherent_frequency_grid_points),
            "frequency_refinement_relative_change": abs(
                frequency_refined.trapping_factor
                - by_order[12].trapping_factor
            ) / frequency_refined.trapping_factor,
            "transition_probability_conservation_error_maximum": (
                frequency_refined
                .transition_probability_conservation_error_maximum),
            "terminal_probability_conservation_error_maximum": (
                frequency_refined
                .terminal_probability_conservation_error_maximum),
            "refined_linear_solver_relative_residual": (
                frequency_refined.linear_solver_relative_residual),
        })
    mape = 100.0 * mean(
        abs(
            row["predicted_zonal_partial_redistribution_trapping_factor"]
            - row["observed_hpep_trapping_factor"]
        ) / row["observed_hpep_trapping_factor"]
        for row in rows
    )
    return {
        "schema": "petch.tian-2017-zonal-partial-redistribution-audit.v1",
        "claim_class": "source-model mechanism reproduction",
        "formal_experimental_validation": False,
        "formal_gate_pass": False,
        "source_model_target_visible_during_model_selection": True,
        "atomic_or_depth_parameter_fitted": False,
        "spatial_field_input": str(FIELD_INPUT.relative_to(ROOT)),
        "target_input": str(TARGET_INPUT.relative_to(ROOT)),
        "condition": source["condition"],
        "model": {
            "spatial_method": (
                "conservative projected-chord finite-volume zone-to-zone "
                "absorption propagator"),
            "frequency_method": (
                "zero-padded FFT natural-Lorentz propagator plus matrix-free "
                "preconditioned GMRES renewal solve"),
            "boundary_method": (
                "conservative lower-endcap/upper-endcap/sidewall terminal "
                "ledger including quenching and analytic far-wing escape"),
            "absorber_density": (
                "local ideal-gas Ar density from printed pressure, mixture, "
                "and contour-bracketed zone temperature"),
            "emitter_density": (
                "contour-bracketed Ar(1s4) field; reused for the unprinted "
                "Ar(1s2) spatial shape as an explicit sensitivity"),
            "velocity_changing_collisions_s_inv": 0.0,
            "quenching_collisions_s_inv": 0.0,
            "fitted_parameters": [],
        },
        "summary": {
            "base_case_line_count": 2,
            "combined_mape_percent": mape,
            "maximum_absolute_line_error_percent": max(
                abs(row["signed_error_percent"]) for row in rows),
            "line_order_reproduced": (
                rows[0]["predicted_zonal_partial_redistribution_trapping_factor"]
                < rows[1]["predicted_zonal_partial_redistribution_trapping_factor"]
            ),
            "formal_gate_pass": False,
        },
        "rows": rows,
        "numerical_convergence": convergence,
        "remaining_blocker": (
            "Tian does not publish the Ar(1s2) emitter field or its numerical "
            "mesh. A raw source-state field export is required for a blind "
            "base-case grade and mixture-dependent spatial fields are required "
            "for the 22-point sweep."),
        "verdict": (
            "Coupled spatial and partial-frequency transport closes the base "
            "mechanism without fitting atomic or depth parameters. Because the "
            "source target was visible during model selection and one line's "
            "emitter field is unpublished, this is a high-accuracy mechanism "
            "reproduction, not experimental validation. The boundary ledger "
            "is numerically usable, but absolute wafer photon flux additionally "
            "requires a line-specific upper-state population field."),
    }


def _report(result: dict[str, object]) -> str:
    lines = [
        "# Tian 2017 zonal partial-redistribution audit",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        "| line | observed HPEM | deterministic zonal PFR | error % | complete-redistribution result |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in result["rows"]:
        lines.append(
            f"| {row['series']} | {row['observed_hpep_trapping_factor']:.3f} | "
            f"{row['predicted_zonal_partial_redistribution_trapping_factor']:.3f} | "
            f"{row['signed_error_percent']:.2f} | "
            f"{row['predicted_zonal_complete_redistribution_trapping_factor']:.3f} |")
    lines.extend([
        "",
        f"Combined two-line MAPE: {result['summary']['combined_mape_percent']:.2f}%.",
        "",
        "## Numerical gate",
        "",
        "| line | order 10→12 trapping change | order 10→12 wafer change | frequency refinement change | transition conservation error | terminal conservation error | GMRES residual |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in result["numerical_convergence"]:
        lines.append(
            f"| {row['series']} | "
            f"{100.0 * row['order_10_to_12_relative_change']:.3f}% | "
            f"{100.0 * row['order_10_to_12_wafer_escape_relative_change']:.3f}% | "
            f"{100.0 * row['frequency_refinement_relative_change']:.5f}% | "
            f"{row['transition_probability_conservation_error_maximum']:.3e} | "
            f"{row['terminal_probability_conservation_error_maximum']:.3e} | "
            f"{row['refined_linear_solver_relative_residual']:.3e} |")
    lines.extend([
        "",
        "## Exact blocker",
        "",
        result["remaining_blocker"],
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
