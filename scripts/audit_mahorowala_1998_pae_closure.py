#!/usr/bin/env python3
"""Compute the 105-nm-equivalent boundary implied by the Cl2 depth gap.

This is an identifiability audit, not a depth fit.  It asks what independently
measured 105-nm response would be required to close each residual on
the best conserved reactor/surface board.  It then compares that requirement
with broadband VUV measurements and with a separately calculated 139-nm
atomic-Cl excitation sensitivity.  The latter wavelength is never fed through
the 105-nm surface yield.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import mean

from petch import Du2022ShortwavePhotoEtchYield
from petch.reactor_global import (
    KemaneciCl139nmEmissionSensitivity,
    malyshev_1998_lam_geometry,
    uniform_isotropic_cylinder_to_disk_transfer,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "coupled_center_anchored_table4_feedback"
    / "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
)
OUTPUT_DIR = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "photoassisted_depth_identifiability"
)
JSON_OUTPUT = OUTPUT_DIR / "mahorowala_1998_pae_closure.json"
REPORT_OUTPUT = OUTPUT_DIR / "MAHOROWALA_1998_PAE_CLOSURE.md"

ETCH_TIME_S = 75.0
DU_WAVELENGTH_NM = 106.0
WAFER_RADIUS_M = 0.100

# Independent amplitude comparators.  They are not treated as interchangeable
# spectra or as values measured in Mahorowala's reactor.
DU_2022_EXTERNAL_VUV_FLUX_CM2_S = (7.6e12, 3.15e13)
DU_2022_YIELD_SI_PER_PHOTON = (90.0, 244.0)
EDELBERG_1999_REACTIVE_WAFER_VUV_FLUX_CM2_S = 0.7e14
EDELBERG_1999_NONREACTIVE_WAFER_VUV_FLUX_CM2_S = 4.0e14
TIAN_2017_95PCT_CL2_TOTAL_VUV_FLUX_CM2_S = 2.3e13
TIAN_2017_95PCT_CL2_VUV_TO_ION_FLUX_RATIO = 0.02
ZHU_2014_LONGWAVE_ANALOG_UPPER_YIELD_SI_PER_PHOTON = 10.0


def _mape(observed: list[float], predicted: list[float]) -> float:
    return float(100.0 * mean(
        abs(model - datum) / datum
        for datum, model in zip(observed, predicted)
    ))


def audit() -> dict[str, object]:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    lower_yield = Du2022ShortwavePhotoEtchYield.measured_lower_yield()
    upper_yield = Du2022ShortwavePhotoEtchYield.measured_upper_yield()
    geometry = malyshev_1998_lam_geometry(6.5).active_geometry
    transparent_unit_transfer = uniform_isotropic_cylinder_to_disk_transfer(
        geometry,
        wafer_radius_m=WAFER_RADIUS_M,
        volume_emissivity_m3_s=1.0,
        quadrature_order=24,
    )
    atomic_emitter = KemaneciCl139nmEmissionSensitivity(
        radiative_survival_fraction=1.0)

    rows = []
    observed_values = []
    baseline_values = []
    for source_row in source["rows"]:
        if source_row["quantitative_status"] != "usable":
            continue
        observed = float(source_row["observed_feature_depth_nm"])
        baseline = float(source_row[
            "coupled_table4_product_reflective_wall_surface_plane_depth_nm_75s"
        ])
        missing_depth = max(0.0, observed - baseline)
        missing_velocity = missing_depth * 1.0e-9 / ETCH_TIME_S
        # Lower yield requires the larger photon flux; upper yield requires the
        # smaller one.
        required_high = lower_yield.required_photon_flux_m2_s(
            missing_velocity, photon_wavelength_nm=DU_WAVELENGTH_NM)
        required_low = upper_yield.required_photon_flux_m2_s(
            missing_velocity, photon_wavelength_nm=DU_WAVELENGTH_NM)
        required_139nm_at_analog_upper_yield = (
            missing_velocity * lower_yield.silicon_atom_density_m3
            / ZHU_2014_LONGWAVE_ANALOG_UPPER_YIELD_SI_PER_PHOTON
        )

        electron_density = float(source_row[
            "coupled_table4_product_reflective_wall_reactor_electron_density_m3"
        ])
        chlorine_density = float(source_row[
            "coupled_table4_product_reflective_wall_reactor_atomic_chlorine_density_m3"
        ])
        temperature = float(source_row[
            "coupled_table4_product_reflective_wall_reactor_mean_electron_energy_eV"
        ]) * 2.0 / 3.0
        excitation_rate = atomic_emitter.primary_excitation_rate_m3_s(
            electron_density_m3=electron_density,
            chlorine_atom_density_m3=chlorine_density,
            electron_temperature_eV=temperature,
        )
        cl139_clear_flux = excitation_rate * transparent_unit_transfer.geometry_flux_length_m
        ion_flux = float(source_row[
            "coupled_table4_product_reflective_wall_total_ion_flux_m2_s"
        ])

        rows.append({
            "run": int(source_row["run"]),
            "inductive_power_W": float(source_row["inductive_power_W"]),
            "rf_bias_power_W": float(source_row["rf_bias_power_W"]),
            "cl2_flow_sccm": float(source_row["cl2_flow_sccm"]),
            "observed_depth_nm": observed,
            "conserved_baseline_depth_nm": baseline,
            "unexplained_depth_nm": missing_depth,
            "unexplained_etch_velocity_m_s": missing_velocity,
            "required_106nm_flux_at_244_yield_cm2_s": required_low / 1.0e4,
            "required_106nm_flux_at_90_yield_cm2_s": required_high / 1.0e4,
            "required_106nm_to_total_ion_flux_ratio_at_244_yield": (
                required_low / ion_flux),
            "required_106nm_to_total_ion_flux_ratio_at_90_yield": (
                required_high / ion_flux),
            "required_139nm_flux_at_10_yield_analog_cm2_s": (
                required_139nm_at_analog_upper_yield / 1.0e4),
            "kemaneci_reaction18_temperature_proxy_eV": temperature,
            "kemaneci_reaction18_excitation_rate_m3_s": excitation_rate,
            "clear_139nm_flux_if_every_excitation_escaped_cm2_s": (
                cl139_clear_flux / 1.0e4),
            "required_effective_shortwave_fraction_of_clear_139nm_proxy_at_244_yield": (
                required_low / cl139_clear_flux),
            "required_effective_shortwave_fraction_of_clear_139nm_proxy_at_90_yield": (
                required_high / cl139_clear_flux),
            "required_139nm_fraction_of_clear_139nm_proxy_at_10_yield_analog": (
                required_139nm_at_analog_upper_yield / cl139_clear_flux),
            "wavelength_transfer_performed": False,
        })
        observed_values.append(observed)
        baseline_values.append(baseline)

    minimum_required = min(
        row["required_106nm_flux_at_244_yield_cm2_s"] for row in rows)
    maximum_required = max(
        row["required_106nm_flux_at_90_yield_cm2_s"] for row in rows)
    maximum_required_beta = max(
        row["required_106nm_to_total_ion_flux_ratio_at_90_yield"]
        for row in rows
    )
    minimum_effective_fraction = min(
        row[
            "required_effective_shortwave_fraction_of_clear_139nm_proxy_at_244_yield"
        ] for row in rows
    )
    maximum_effective_fraction = max(
        row[
            "required_effective_shortwave_fraction_of_clear_139nm_proxy_at_90_yield"
        ] for row in rows
    )
    maximum_longwave_proxy_fraction = max(
        row["required_139nm_fraction_of_clear_139nm_proxy_at_10_yield_analog"]
        for row in rows
    )
    return {
        "schema": "petch.mahorowala-1998-pae-identifiability.v1",
        "claim_class": (
            "source-bounded omitted-mechanism amplitude and identifiability audit"
        ),
        "input_board": str(INPUT.relative_to(ROOT)),
        "feature_depth_used_to_select_photon_boundary": True,
        "formal_prediction_or_validation": False,
        "etch_time_s": ETCH_TIME_S,
        "radiation_transfer": {
            "method": (
                "deterministic disk-overlap convolution with analytic axial "
                "cosine-law integration; no photon Monte Carlo"
            ),
            "geometry_radius_m": geometry.radius_m,
            "geometry_length_m": geometry.length_m,
            "wafer_radius_m": WAFER_RADIUS_M,
            "clear_geometry_flux_length_m": (
                transparent_unit_transfer.geometry_flux_length_m),
            "clear_wafer_intercept_probability": (
                transparent_unit_transfer.wafer_intercept_probability),
            "wall_reflection_included": False,
            "resonance_trapping_or_quenching_included": False,
        },
        "independent_evidence": {
            "du_2022_surface_wavelength_nm": [104.82, 106.67],
            "du_2022_measured_yield_si_per_photon": list(
                DU_2022_YIELD_SI_PER_PHOTON),
            "du_2022_external_vuv_flux_cm2_s": list(
                DU_2022_EXTERNAL_VUV_FLUX_CM2_S),
            "edelberg_1999_commercial_cl2_bcl3_reactive_wafer_total_vuv_flux_cm2_s": (
                EDELBERG_1999_REACTIVE_WAFER_VUV_FLUX_CM2_S),
            "edelberg_1999_commercial_cl2_bcl3_nonreactive_wafer_total_vuv_flux_cm2_s": (
                EDELBERG_1999_NONREACTIVE_WAFER_VUV_FLUX_CM2_S),
            "tian_2017_95pct_cl2_modeled_total_vuv_flux_cm2_s": (
                TIAN_2017_95PCT_CL2_TOTAL_VUV_FLUX_CM2_S),
            "tian_2017_95pct_cl2_modeled_vuv_to_ion_flux_ratio": (
                TIAN_2017_95PCT_CL2_VUV_TO_ION_FLUX_RATIO),
            "hirsch_2020_rf_or_pulsed_ion_anti_synergy": (
                "energetic-ion exposure suppresses PAE; target RF waveform "
                "has no published photon-resolved suppression law"
            ),
            "zhu_2014_longwave_response": (
                "photons above 120 nm contributed about 10% of PAE in the "
                "Cl2/Kr window experiment; cited analogous 130--135 nm "
                "photoetch yields were about 0--10 atoms/photon"
            ),
        },
        "summary": {
            "usable_run_count": len(rows),
            "conserved_baseline_mape_percent": _mape(
                observed_values, baseline_values),
            "required_effective_106nm_flux_range_cm2_s": [
                minimum_required, maximum_required],
            "maximum_required_effective_106nm_to_ion_flux_ratio": (
                maximum_required_beta),
            "required_effective_shortwave_fraction_of_clear_139nm_excitation_proxy_range": [
                minimum_effective_fraction, maximum_effective_fraction],
            "maximum_required_139nm_fraction_of_clear_139nm_proxy_at_10_yield_analog": (
                maximum_longwave_proxy_fraction),
            "clear_139nm_proxy_closes_all_runs_at_10_yield_analog": (
                maximum_longwave_proxy_fraction <= 1.0),
            "broadband_amplitude_overlaps_independent_vuv_evidence": (
                maximum_required
                <= EDELBERG_1999_NONREACTIVE_WAFER_VUV_FLUX_CM2_S
            ),
            "spectrally_supported_depth_closure": False,
            "formal_depth_closure_granted": False,
        },
        "rows": rows,
        "missing_for_formal_depth_closure": [
            "absolute line-resolved 109--120-nm photon flux at Mahorowala's wafer for every run",
            "chlorinated-Si photoetch yields at the 109.7, 110.8, and 118.9 nm Cl I lines",
            "RF-waveform-resolved PAE/IAE anti-synergy on chlorinated poly-Si",
            "surface photo-etch yield at the measured reactor spectrum and 60 C",
            "feature-scale VUV view-factor and photo-product escape validation",
        ],
        "verdict": (
            "The residual can be written as a 105-nm-equivalent photon target, "
            "but broadband amplitude overlap is not spectral closure. A subsequent "
            "line-resolved Cl I audit finds negligible atomic emission in Du's "
            "104.82--106.67 nm response band and strong unmeasured 109--120 nm "
            "lines. This scalar audit therefore grants no photo-assisted depth "
            "closure."
        ),
    }


def _report(result: dict[str, object]) -> str:
    summary = result["summary"]
    rows = result["rows"]
    lines = [
        "# Mahorowala 1998 photo-assisted depth identifiability",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        f"The conserved no-PAE board has {summary['conserved_baseline_mape_percent']:.2f}% "
        "MAPE. Expressing its positive residuals in units of Du's independently "
        "measured 90--244 Si/photon response gives a 106-nm-equivalent flux spanning "
        f"{summary['required_effective_106nm_flux_range_cm2_s'][0]:.3g}--"
        f"{summary['required_effective_106nm_flux_range_cm2_s'][1]:.3g} cm^-2 s^-1. "
        "That is a target-variable identity, not a fitted validation score or a "
        "spectrally supported source.",
        "",
        "| run | source W | bias W | flow sccm | observed nm | baseline nm | missing nm | required 106 nm flux, 244 yield | required 106 nm flux, 90 yield | max beta |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['run']} | {row['inductive_power_W']:.0f} | "
            f"{row['rf_bias_power_W']:.0f} | {row['cl2_flow_sccm']:.0f} | "
            f"{row['observed_depth_nm']:.1f} | "
            f"{row['conserved_baseline_depth_nm']:.1f} | "
            f"{row['unexplained_depth_nm']:.1f} | "
            f"{row['required_106nm_flux_at_244_yield_cm2_s']:.3g} | "
            f"{row['required_106nm_flux_at_90_yield_cm2_s']:.3g} | "
            f"{row['required_106nm_to_total_ion_flux_ratio_at_90_yield']:.4f} |"
        )
    lines.extend([
        "",
        "## Guardrail",
        "",
        "Kemaneci reaction 18 produces a 139-nm excitation proxy. Du's absolute "
        "surface yield is restricted to 104.82--106.67 nm, so the code does not "
        "apply that yield at 139 nm. The reported proxy fractions combine an "
        "unknown shortwave branch, radiative escape, and spectral surface response; "
        "they are target requirements, not calibrated reactor constants. The "
        "line-resolved follow-up also forbids treating all sub-120-nm photons as "
        "though they carried Du's 105-nm response.",
        "",
        "The independent wavelength result is also carried as a negative bound: "
        "using 10 Si/photon as a nonpredictive 130--139 nm analog requires as much "
        f"as {summary['maximum_required_139nm_fraction_of_clear_139nm_proxy_at_10_yield_analog']:.2f}x "
        "the transparent one-pass Kemaneci reaction-18 proxy. Thus 139-nm light "
        "alone is not an amplitude closure for all runs under that response.",
        "",
        "## Exact experiment that closes the board",
        "",
        "Measure the absolute sub-120-nm spectrum at the wafer synchronously with "
        "the species-resolved ion flux/IED under each Table-2.2 condition, then "
        "measure the RF-waveform-resolved Si photo-etch response with ions optically "
        "blocked/unblocked. Those independent boundaries make the depth board "
        "out-of-sample. No feature depth may be used to choose their normalization.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    result = audit()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    JSON_OUTPUT.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    REPORT_OUTPUT.write_text(_report(result), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
