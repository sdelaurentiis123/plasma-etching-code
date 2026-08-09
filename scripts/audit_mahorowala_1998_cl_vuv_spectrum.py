#!/usr/bin/env python3
"""Resolve the missing Mahorowala depth against a line-resolved Cl I source.

This is a fail-closed source and measurement-design audit.  It uses no feature
depth to alter the reactor, atomic, radiation, or surface parameters.  Only
the 104.82--106.67 nm source is passed through Du's measured photoetch yield.
For the calculated 106.67--112 and 112--120 nm bands, the script reports the
surface yield that an independent experiment would need to measure; it never
borrows the 105-nm yield across wavelength.

OPEN-ADAS source files are caller supplied and license restricted.  They are
read only after the explicit CLI acknowledgement and are never copied into the
result or repository.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import mean

from petch import Du2022ShortwavePhotoEtchYield
from petch.reactor_global import (
    OPEN_ADAS_CL0_COLLISION_RECORDS_SHA256,
    OPEN_ADAS_CL0_NIST_RECORDS_SHA256,
    chlorine_direct_coronal_spectrum,
    deterministic_uniform_chlorine_line_wafer_boundary,
    load_open_adas_cl0_personal_research,
    malyshev_1998_lam_geometry,
    uniform_isotropic_cylinder_direct_lamellar_floor_transfer,
    uniform_isotropic_cylinder_to_disk_transfer,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "coupled_center_anchored_table4_feedback"
    / "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
)
FEATURE_INPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_center_anchored_feature_depth_dx40nm_no_reflection.json"
)
OUTPUT_DIR = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "photoassisted_spectral_audit"
)
JSON_OUTPUT = OUTPUT_DIR / "mahorowala_1998_cl_vuv_spectrum.json"
REPORT_OUTPUT = OUTPUT_DIR / "MAHOROWALA_1998_CL_VUV_SPECTRUM.md"

ETCH_TIME_S = 75.0
WAFER_RADIUS_M = 0.100
OPENING_WIDTH_M = 310.0e-9
OXIDE_MASK_THICKNESS_M = 200.0e-9
DU_SUPPORTED_BAND_NM = (104.82, 106.67)
UNMEASURED_SHORT_BANDS_NM = ((106.67, 112.0), (112.0, 120.0))


def _mape(observed: list[float], predicted: list[float]) -> float:
    return float(100.0 * mean(
        abs(model - datum) / datum
        for datum, model in zip(observed, predicted)
    ))


def _finite_range(values: list[float | None]) -> list[float] | None:
    finite = [float(value) for value in values if value is not None]
    return [min(finite), max(finite)] if finite else None


def _positive_range(values: list[float | None]) -> list[float] | None:
    finite = [
        float(value) for value in values
        if value is not None and float(value) > 0.0
    ]
    return [min(finite), max(finite)] if finite else None


def audit(collision_path: Path, nist_level_path: Path) -> dict[str, object]:
    collision, observed_levels = load_open_adas_cl0_personal_research(
        collision_path,
        nist_level_path,
        accept_restricted_personal_use=True,
    )
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    feature = json.loads(FEATURE_INPUT.read_text(encoding="utf-8"))
    feature_by_run = {int(row["run"]): row for row in feature["rows"]}
    geometry = malyshev_1998_lam_geometry(6.5).active_geometry
    unit_transfer = uniform_isotropic_cylinder_to_disk_transfer(
        geometry,
        wafer_radius_m=WAFER_RADIUS_M,
        volume_emissivity_m3_s=1.0,
        quadrature_order=24,
    )
    du_lower = Du2022ShortwavePhotoEtchYield.measured_lower_yield()
    du_upper = Du2022ShortwavePhotoEtchYield.measured_upper_yield()
    silicon_density = du_lower.silicon_atom_density_m3

    rows: list[dict[str, object]] = []
    observed_depths: list[float] = []
    baseline_depths: list[float] = []
    supported_low_depths: list[float] = []
    supported_high_depths: list[float] = []
    resonance_quadrature_receipt = None
    for source_row in source["rows"]:
        if source_row["quantitative_status"] != "usable":
            continue
        feature_row = feature_by_run[int(source_row["run"])]
        observed = float(source_row["observed_feature_depth_nm"])
        if not math.isclose(
            observed,
            float(feature_row["observed_feature_depth_nm"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError("reactor and feature boards disagree on target depth")
        baseline = float(feature_row["predicted_center_feature_depth_nm"])
        missing = max(0.0, observed - baseline)
        missing_atom_rate_m2_s = missing * 1.0e-9 / ETCH_TIME_S * silicon_density
        electron_density = float(source_row[
            "coupled_table4_product_reflective_wall_reactor_electron_density_m3"
        ])
        chlorine_density = float(source_row[
            "coupled_table4_product_reflective_wall_reactor_atomic_chlorine_density_m3"
        ])
        temperature_eV = 2.0 / 3.0 * float(source_row[
            "coupled_table4_product_reflective_wall_reactor_mean_electron_energy_eV"
        ])
        spectrum = chlorine_direct_coronal_spectrum(
            collision,
            observed_levels,
            electron_temperature_eV=temperature_eV,
        )
        target_floor_transfer = (
            uniform_isotropic_cylinder_direct_lamellar_floor_transfer(
                geometry,
                wafer_radius_m=WAFER_RADIUS_M,
                opening_width_m=OPENING_WIDTH_M,
                optical_path_depth_m=(
                    OXIDE_MASK_THICKNESS_M + observed * 1.0e-9
                ),
                volume_emissivity_m3_s=1.0,
                radial_quadrature_order=32,
                axial_quadrature_order=32,
            )
        )
        direct_floor_fraction = target_floor_transfer.direct_floor_fraction

        def band_values(band: tuple[float, float]) -> dict[str, float | None]:
            coefficient_cm3_s = spectrum.band_rate_coefficient_cm3_s(*band)
            emissivity_m3_s = coefficient_cm3_s * 1.0e-6 * (
                electron_density * chlorine_density
            )
            clear_flux_m2_s = (
                emissivity_m3_s * unit_transfer.geometry_flux_length_m
            )
            required_yield = (
                missing_atom_rate_m2_s / clear_flux_m2_s
                if clear_flux_m2_s > 0.0 else None
            )
            return {
                "minimum_wavelength_nm": band[0],
                "maximum_wavelength_nm": band[1],
                "direct_coronal_rate_coefficient_cm3_s": coefficient_cm3_s,
                "uniform_clear_volume_emissivity_m3_s": emissivity_m3_s,
                "uniform_clear_wafer_flux_cm2_s": clear_flux_m2_s / 1.0e4,
                "residual_closing_yield_si_per_photon_if_unit_feature_delivery": (
                    required_yield
                ),
                "residual_rate_divided_by_target_depth_direct_floor_flux_si_per_photon": (
                    None if required_yield is None
                    else required_yield / direct_floor_fraction
                ),
            }

        supported = band_values(DU_SUPPORTED_BAND_NM)
        unmeasured = [band_values(band) for band in UNMEASURED_SHORT_BANDS_NM]
        supported_flux_m2_s = float(supported["uniform_clear_wafer_flux_cm2_s"]) * 1.0e4
        supported_low_addition_nm = (
            du_lower.etch_velocity_m_s(
                supported_flux_m2_s, photon_wavelength_nm=106.0
            ) * ETCH_TIME_S * 1.0e9
        )
        supported_high_addition_nm = (
            du_upper.etch_velocity_m_s(
                supported_flux_m2_s, photon_wavelength_nm=106.0
            ) * ETCH_TIME_S * 1.0e9
        )
        sub120_lines = sorted(
            (
                line for line in spectrum.lines
                if 104.82 <= line.wavelength_nm < 120.0
            ),
            key=lambda line: line.photon_rate_coefficient_cm3_s,
            reverse=True,
        )
        dominant_118p88 = max(
            (
                line for line in sub120_lines
                if 118.86 <= line.wavelength_nm <= 118.90
                and line.transition_probability_s_inv >= 1.0e8
            ),
            key=lambda line: line.photon_rate_coefficient_cm3_s,
        )
        resonance_boundary = (
            deterministic_uniform_chlorine_line_wafer_boundary(
                geometry,
                dominant_118p88,
                observed_levels,
                wafer_radius_m=WAFER_RADIUS_M,
                electron_density_m3=electron_density,
                chlorine_atom_density_m3=chlorine_density,
                gas_temperature_K=500.0,
                velocity_changing_collision_frequency_s_inv=0.0,
                nonradiative_quenching_frequency_s_inv=0.0,
                surface_quadrature_order=12,
                direction_quadrature_order=12,
                frequency_quadrature_order=24,
                coherent_grid_points_per_lorentz_hwhm=8.0,
            )
        )
        if int(source_row["run"]) == 1:
            receipt_levels = []
            for label, surface_order, direction_order in (
                ("coarse", 6, 6),
                ("production", 12, 12),
                ("fine", 16, 16),
            ):
                receipt = (
                    resonance_boundary
                    if label == "production"
                    else deterministic_uniform_chlorine_line_wafer_boundary(
                        geometry,
                        dominant_118p88,
                        observed_levels,
                        wafer_radius_m=WAFER_RADIUS_M,
                        electron_density_m3=electron_density,
                        chlorine_atom_density_m3=chlorine_density,
                        gas_temperature_K=500.0,
                        surface_quadrature_order=surface_order,
                        direction_quadrature_order=direction_order,
                        frequency_quadrature_order=24,
                        coherent_grid_points_per_lorentz_hwhm=8.0,
                    )
                )
                receipt_levels.append({
                    "label": label,
                    "surface_quadrature_order": surface_order,
                    "direction_quadrature_order": direction_order,
                    "frequency_quadrature_order": 24,
                    "coherent_grid_points_per_lorentz_hwhm": 8.0,
                    "wafer_escape_probability": (
                        receipt.radiation
                        .partial_redistribution_wafer_escape_probability),
                    "wafer_photon_flux_m2_s": receipt.wafer_photon_flux_m2_s,
                })
            production_to_fine = abs(
                receipt_levels[1]["wafer_photon_flux_m2_s"]
                / receipt_levels[2]["wafer_photon_flux_m2_s"] - 1.0
            )
            resonance_quadrature_receipt = {
                "condition_run": 1,
                "levels": receipt_levels,
                "production_to_fine_relative_change": production_to_fine,
                "passed_one_percent_gate": production_to_fine < 0.01,
                "partial_wafer_radial_discontinuity_split_exactly": True,
            }
        resonance_required_yield = (
            missing_atom_rate_m2_s / resonance_boundary.wafer_photon_flux_m2_s
            if resonance_boundary.wafer_photon_flux_m2_s > 0.0 else None
        )
        rows.append({
            "run": int(source_row["run"]),
            "inductive_power_W": float(source_row["inductive_power_W"]),
            "rf_bias_power_W": float(source_row["rf_bias_power_W"]),
            "cl2_flow_sccm": float(source_row["cl2_flow_sccm"]),
            "observed_depth_nm": observed,
            "conserved_nonphoton_feature_depth_nm": baseline,
            "corresponding_surface_plane_depth_nm": float(
                feature_row["surface_plane_depth_nm"]),
            "unexplained_depth_nm": missing,
            "target_optical_path_depth_nm": (
                OXIDE_MASK_THICKNESS_M * 1.0e9 + observed
            ),
            "direct_geometrical_floor_fraction_at_target_depth": (
                direct_floor_fraction
            ),
            "electron_temperature_proxy_eV": temperature_eV,
            "electron_density_m3": electron_density,
            "atomic_chlorine_density_m3": chlorine_density,
            "du_supported_band": supported,
            "du_supported_band_depth_addition_at_90_yield_nm": (
                supported_low_addition_nm
            ),
            "du_supported_band_depth_addition_at_244_yield_nm": (
                supported_high_addition_nm
            ),
            "unmeasured_shortwave_bands": unmeasured,
            "dominant_118p88_resonance_boundary": {
                "wavelength_nm": resonance_boundary.wavelength_nm,
                "primary_line_emissivity_m3_s": (
                    resonance_boundary.primary_line_emissivity_m3_s),
                "primary_line_emission_rate_s": (
                    resonance_boundary.primary_line_emission_rate_s),
                "lower_ground_population_fraction": (
                    resonance_boundary.lower_ground_population_fraction),
                "lower_state_absorber_density_m3": (
                    resonance_boundary.lower_state_absorber_density_m3),
                "alternate_branch_loss_frequency_s_inv": (
                    resonance_boundary.alternate_branch_loss_frequency_s_inv),
                "partial_redistribution_wafer_escape_probability": (
                    resonance_boundary.radiation
                    .partial_redistribution_wafer_escape_probability),
                "partial_redistribution_quench_probability": (
                    resonance_boundary.radiation
                    .partial_redistribution_quench_probability),
                "terminal_probability_conservation_error_maximum": (
                    resonance_boundary.radiation
                    .terminal_probability_conservation_error_maximum),
                "transition_probability_conservation_error_maximum": (
                    resonance_boundary.radiation
                    .transition_probability_conservation_error_maximum),
                "linear_solver_relative_residual": (
                    resonance_boundary.radiation
                    .linear_solver_relative_residual),
                "linear_solver_iterations": (
                    resonance_boundary.radiation.linear_solver_iterations),
                "wafer_photon_flux_cm2_s": (
                    resonance_boundary.wafer_photon_flux_m2_s / 1.0e4),
                "residual_closing_yield_at_wafer_si_per_photon": (
                    resonance_required_yield),
                "residual_closing_yield_at_direct_target_floor_si_per_photon": (
                    None if resonance_required_yield is None
                    else resonance_required_yield / direct_floor_fraction
                ),
                "velocity_changing_collision_frequency_s_inv": 0.0,
                "nonradiative_quenching_frequency_s_inv": 0.0,
                "prediction_supported": False,
                "known_limitations": list(
                    resonance_boundary.known_limitations),
            },
            "dominant_sub120nm_lines": [
                {
                    "wavelength_nm": line.wavelength_nm,
                    "upper_observed_index": line.upper_observed_index,
                    "lower_observed_index": line.lower_observed_index,
                    "transition_probability_s_inv": (
                        line.transition_probability_s_inv
                    ),
                    "photon_rate_coefficient_cm3_s": (
                        line.photon_rate_coefficient_cm3_s
                    ),
                }
                for line in sub120_lines[:8]
            ],
            "feature_depth_used_to_select_any_model_parameter": False,
            "feature_scale_photon_delivery_assumed": (
                "unity only to report an optimistic lower bound on required yield"
            ),
        })
        observed_depths.append(observed)
        baseline_depths.append(baseline)
        supported_low_depths.append(baseline + supported_low_addition_nm)
        supported_high_depths.append(baseline + supported_high_addition_nm)

    supported_coefficients = [
        float(row["du_supported_band"]["direct_coronal_rate_coefficient_cm3_s"])
        for row in rows
    ]
    first_unmeasured_coefficients = [
        float(row["unmeasured_shortwave_bands"][0][
            "direct_coronal_rate_coefficient_cm3_s"
        ]) for row in rows
    ]
    second_unmeasured_coefficients = [
        float(row["unmeasured_shortwave_bands"][1][
            "direct_coronal_rate_coefficient_cm3_s"
        ]) for row in rows
    ]
    first_required_yields = [
        row["unmeasured_shortwave_bands"][0][
            "residual_closing_yield_si_per_photon_if_unit_feature_delivery"
        ] for row in rows
    ]
    second_required_yields = [
        row["unmeasured_shortwave_bands"][1][
            "residual_closing_yield_si_per_photon_if_unit_feature_delivery"
        ] for row in rows
    ]
    first_target_floor_yields = [
        row["unmeasured_shortwave_bands"][0][
            "residual_rate_divided_by_target_depth_direct_floor_flux_si_per_photon"
        ] for row in rows
    ]
    second_target_floor_yields = [
        row["unmeasured_shortwave_bands"][1][
            "residual_rate_divided_by_target_depth_direct_floor_flux_si_per_photon"
        ] for row in rows
    ]
    target_floor_fractions = [
        row["direct_geometrical_floor_fraction_at_target_depth"] for row in rows
    ]
    resonance_wafer_fluxes = [
        row["dominant_118p88_resonance_boundary"]["wafer_photon_flux_cm2_s"]
        for row in rows
    ]
    resonance_wafer_yields = [
        row["dominant_118p88_resonance_boundary"][
            "residual_closing_yield_at_wafer_si_per_photon"
        ] for row in rows
    ]
    resonance_floor_yields = [
        row["dominant_118p88_resonance_boundary"][
            "residual_closing_yield_at_direct_target_floor_si_per_photon"
        ] for row in rows
    ]
    return {
        "schema": "petch.mahorowala-1998-cl-vuv-spectrum-audit.v2",
        "claim_class": "atomic-source sensitivity and exact measurement target",
        "formal_depth_prediction_or_validation": False,
        "feature_depth_fit_used": False,
        "input_boards": {
            "reactor_boundary": str(INPUT.relative_to(ROOT)),
            "deterministic_feature": str(FEATURE_INPUT.relative_to(ROOT)),
        },
        "atomic_source": {
            "method": (
                "observed NIST level separations plus OPEN-ADAS AUTOSTRUCTURE "
                "distorted-wave effective collision strengths and branching"
            ),
            "collision_physical_records_sha256": (
                OPEN_ADAS_CL0_COLLISION_RECORDS_SHA256
            ),
            "nist_level_physical_records_sha256": (
                OPEN_ADAS_CL0_NIST_RECORDS_SHA256
            ),
            "raw_files_redistributed": False,
            "open_adas_personal_use_restriction_acknowledged": True,
            "matched_observed_level_count": (
                chlorine_direct_coronal_spectrum(
                    collision,
                    observed_levels,
                    electron_temperature_eV=2.5,
                ).matched_observed_level_count
            ),
            "calculated_level_count": len(collision.levels),
            "resolved_observed_level_count": len(observed_levels),
            "included_physics": [
                "ground-fine-structure Boltzmann populations",
                "direct electron excitation",
                "radiative branching over all calculated lower levels",
            ],
            "missing_physics": [
                "excited-state cascades",
                "radiation trapping and escape outside the separately solved dominant 118.88 nm line",
                "collisional and wall quenching",
                "measured transition probabilities for every line",
            ],
            "prediction_supported": False,
        },
        "radiation_transfer": {
            "method": (
                "transparent one-pass band integral plus conservative zonal "
                "partial-frequency redistribution for dominant 118.88 nm"
            ),
            "geometry_flux_length_m": unit_transfer.geometry_flux_length_m,
            "wafer_intercept_probability": unit_transfer.wafer_intercept_probability,
            "uniform_emissivity_assumed": True,
            "resonance_transfer_included_for_dominant_118p88_nm": True,
            "feature_scale_delivery_included": False,
            "direct_absorbing_lamellar_floor_sensitivity_included": True,
            "line_opening_width_nm": OPENING_WIDTH_M * 1.0e9,
            "oxide_mask_thickness_nm": OXIDE_MASK_THICKNESS_M * 1.0e9,
            "wave_electromagnetics_included": False,
            "dominant_118p88_quadrature_receipt": (
                resonance_quadrature_receipt),
        },
        "surface_response": {
            "du_measured_band_nm": list(DU_SUPPORTED_BAND_NM),
            "du_measured_yield_si_per_photon": [
                du_lower.silicon_atoms_per_photon,
                du_upper.silicon_atoms_per_photon,
            ],
            "du_yield_applied_outside_measured_band": False,
            "rf_ion_pae_antisynergy_included": False,
        },
        "summary": {
            "usable_run_count": len(rows),
            "conserved_nonphoton_feature_mape_percent": _mape(
                observed_depths, baseline_depths
            ),
            "du_supported_atomic_band_mape_range_percent": [
                _mape(observed_depths, supported_low_depths),
                _mape(observed_depths, supported_high_depths),
            ],
            "du_supported_104p82_106p67_rate_coefficient_range_cm3_s": (
                _finite_range(supported_coefficients)
            ),
            "unmeasured_106p67_112_rate_coefficient_range_cm3_s": (
                _finite_range(first_unmeasured_coefficients)
            ),
            "unmeasured_112_120_rate_coefficient_range_cm3_s": (
                _finite_range(second_unmeasured_coefficients)
            ),
            "required_106p67_112_yield_if_unit_feature_delivery_range_si_per_photon": (
                _positive_range(first_required_yields)
            ),
            "required_112_120_yield_if_unit_feature_delivery_range_si_per_photon": (
                _positive_range(second_required_yields)
            ),
            "target_depth_direct_geometrical_floor_fraction_range": (
                _finite_range(target_floor_fractions)
            ),
            "dominant_118p88_partial_redistribution_wafer_flux_range_cm2_s": (
                _finite_range(resonance_wafer_fluxes)
            ),
            "dominant_118p88_residual_closing_yield_at_wafer_range_si_per_photon": (
                _positive_range(resonance_wafer_yields)
            ),
            "dominant_118p88_residual_closing_yield_at_direct_target_floor_range_si_per_photon": (
                _positive_range(resonance_floor_yields)
            ),
            "residual_rate_over_target_depth_direct_floor_flux_106p67_112_range_si_per_photon": (
                _positive_range(first_target_floor_yields)
            ),
            "residual_rate_over_target_depth_direct_floor_flux_112_120_range_si_per_photon": (
                _positive_range(second_target_floor_yields)
            ),
            "underpredicted_run_count_with_positive_residual": sum(
                float(row["unexplained_depth_nm"]) > 0.0 for row in rows
            ),
            "nonunderpredicted_runs_not_repairable_by_additive_photons": [
                int(row["run"]) for row in rows
                if float(row["unexplained_depth_nm"]) == 0.0
            ],
            "formal_depth_closure_granted": False,
        },
        "rows": rows,
        "exact_experimental_blockers": [
            "absolute line-resolved 109.2--110.9 and 118.88 nm wafer flux under each target condition",
            "chlorinated-poly-Si photoetch yield at 109.7, 110.8, and 118.9 nm at 60 C",
            "the same yields under the target 13.56 MHz RF ion waveform to resolve PAE/IAE anti-synergy",
            "electromagnetic validation of the direct geometrical floor sensitivity for the 310 nm opening (2.6--2.8 wavelengths wide)",
        ],
        "verdict": (
            "The atomic-Cl source does not populate the band where the large "
            "90--244 Si/photon yield was measured, so that measurement cannot "
            "close Mahorowala depth. The source instead predicts strong 109--111 "
            "and 118.88 nm lines, but their wavelength-resolved chlorinated-Si "
            "yield is unpublished. Depth is therefore experimentally blocked at "
            "a much narrower interface: those line yields, their RF anti-synergy, "
            "and feature-floor photon delivery."
        ),
    }


def _format_range(values: list[float] | None) -> str:
    if values is None:
        return "unavailable"
    return f"{values[0]:.3g}--{values[1]:.3g}"


def _report(result: dict[str, object]) -> str:
    summary = result["summary"]
    quadrature = result["radiation_transfer"][
        "dominant_118p88_quadrature_receipt"
    ]
    lines = [
        "# Mahorowala 1998 line-resolved chlorine VUV audit",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        f"The conserved non-photon full-feature board remains at "
        f"{summary['conserved_nonphoton_feature_mape_percent']:.3f}% MAPE. "
        "This is the deterministic 40 nm feature calculation, not the planar "
        "surface-rate projection. Passing only "
        "the calculated 104.82--106.67 nm Cl I source through Du's measured "
        "yield gives a MAPE interval of "
        f"{summary['du_supported_atomic_band_mape_range_percent'][0]:.3f}--"
        f"{summary['du_supported_atomic_band_mape_range_percent'][1]:.3f}%. "
        "It does not close the board.",
        "",
        "| band nm | calculated direct-coronal rate coefficient cm^3/s | yield at unit delivery | residual rate / target-depth direct-floor flux | status |",
        "|---|---:|---:|---:|---|",
        f"| 104.82--106.67 | {_format_range(summary['du_supported_104p82_106p67_rate_coefficient_range_cm3_s'])} | measured 90--244 | not a reconstruction | supported but negligible source |",
        f"| 106.67--112 | {_format_range(summary['unmeasured_106p67_112_rate_coefficient_range_cm3_s'])} | {_format_range(summary['required_106p67_112_yield_if_unit_feature_delivery_range_si_per_photon'])} | {_format_range(summary['residual_rate_over_target_depth_direct_floor_flux_106p67_112_range_si_per_photon'])} Si/photon | source sensitivity; surface response unmeasured |",
        f"| 112--120 | {_format_range(summary['unmeasured_112_120_rate_coefficient_range_cm3_s'])} | {_format_range(summary['required_112_120_yield_if_unit_feature_delivery_range_si_per_photon'])} | {_format_range(summary['residual_rate_over_target_depth_direct_floor_flux_112_120_range_si_per_photon'])} Si/photon | source sensitivity; surface response unmeasured |",
        "",
        "The dominant calculated shortwave lines are 118.88, 109.74, and "
        "110.75 nm. OPEN-ADAS raw data are not redistributed; the exact "
        "physical-record hashes and license boundary are retained in the JSON.",
        "",
        "The dominant 118.88-nm line is now propagated with conservative "
        "partial-frequency redistribution, ground-fine-structure absorption, "
        "alternate radiative branches, and a finite 200-mm wafer. Its wafer "
        "flux spans "
        f"{_format_range(summary['dominant_118p88_partial_redistribution_wafer_flux_range_cm2_s'])} cm^-2 s^-1. "
        "Closing only the positive depth residual would require "
        f"{_format_range(summary['dominant_118p88_residual_closing_yield_at_wafer_range_si_per_photon'])} Si/photon at the wafer, or "
        f"{_format_range(summary['dominant_118p88_residual_closing_yield_at_direct_target_floor_range_si_per_photon'])} Si/photon after the target-depth direct-floor sensitivity. "
        "Those are measurement targets, not fitted yields; velocity-changing "
        "and nonradiative collision frequencies are still unset.",
        "",
        "The released 118.88-nm transport is quadrature checked: the 12x12 "
        "production surface/direction rule changes by only "
        f"{100.0 * quadrature['production_to_fine_relative_change']:.3f}% "
        "against the 16x16 rule (1% gate: "
        f"`{quadrature['passed_one_percent_gate']}`). The finite-wafer radial "
        "discontinuity is split exactly rather than crossed by an unsplit "
        "Gaussian panel.",
        "",
        "The exact absorbing-ray cylinder-to-line integral delivers only "
        f"{summary['target_depth_direct_geometrical_floor_fraction_range'][0]:.3f}--"
        f"{summary['target_depth_direct_geometrical_floor_fraction_range'][1]:.3f} "
        "of wafer-plane photons to the target-depth floors. This is a useful "
        "deterministic transport sensitivity, not the final photon boundary: "
        "the 310 nm opening is only 2.6--2.8 wavelengths across at the dominant "
        "lines, so wave-optical validation remains mandatory.",
        "",
        "A strictly additive photon channel also cannot be the complete board "
        "repair: runs "
        + ", ".join(str(run) for run in summary[
            "nonunderpredicted_runs_not_repairable_by_additive_photons"])
        + " are already at or above their measured depths. Any successful "
        "spectrum-resolved mechanism must reproduce the condition dependence "
        "and RF ion/photon anti-synergy, not add one global depth offset.",
        "",
        "## Experiment that decides depth",
        "",
    ]
    lines.extend(f"- {item}" for item in result["exact_experimental_blockers"])
    lines.extend([
        "",
        "No target depth was used to select a reactor, atomic, radiative, or "
        "surface parameter. The reported unmeasured-band yields are measurement "
        "targets, not fitted constants. Both a unity-delivery lower bound and a "
        "target-depth absorbing-ray sensitivity are shown; neither substitutes "
        "for wavelength-resolved electromagnetic feature validation.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adas-collision", type=Path, required=True)
    parser.add_argument("--nist-levels", type=Path, required=True)
    parser.add_argument(
        "--accept-open-adas-personal-use",
        action="store_true",
        help="acknowledge that the caller-provided files are restricted to personal use",
    )
    args = parser.parse_args()
    if not args.accept_open_adas_personal_use:
        parser.error("--accept-open-adas-personal-use is required")
    result = audit(args.adas_collision, args.nist_levels)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    JSON_OUTPUT.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    REPORT_OUTPUT.write_text(_report(result), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
