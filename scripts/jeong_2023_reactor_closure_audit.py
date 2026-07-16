#!/usr/bin/env python3
"""Diagnose the missing Jeong-2023 reactor-to-feature closure without profile fitting.

This audit is intentionally zero-dimensional and consumes the 200 nm flux sweep only
as *development diagnostics*.  It cannot convert those points into held-out validation.
Its job is to quantify which missing boundary response would be required, combine the
already completed source-backed kill tests, and decide whether another moving-profile
run is scientifically earned.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq

from petch.experimental_boundary import build_jeong_2023_boundary_state
from petch.experimental_data import (
    load_jeong_2023_etch_depths,
    load_jeong_2023_radical_densities,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
    SurfaceFluxes,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "jeong_2023"
DEFAULT_OUTPUT = (
    ROOT / "results" / "jeong_2023_reactor_closure_audit" / "audit.json")


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git_revision():
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True).stdout.strip()


def _atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _write_plot(path, payload):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    bounded = payload["bounded_source_backed_tests"][
        "low_energy_fc_ion_direct_polymer_deposition"]
    inverse = payload["inverse_missing_response_diagnostic"]["rows"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)

    labels = ("Experiment", "Current 0-D", "20% low-E CF+")
    gains = (
        bounded["experimental_endpoint_gain_nm"],
        bounded["baseline_endpoint_gain_nm"],
        bounded["strongest_endpoint_gain_nm"],
    )
    colors = ("#1f9d8a", "#d95f59", "#e5a93d")
    axes[0].bar(labels, gains, color=colors)
    axes[0].set_ylabel("200 nm density-sweep endpoint gain (nm)")
    axes[0].set_title("Source-backed additions do not close the slope")
    axes[0].tick_params(axis="x", rotation=14)
    for index, value in enumerate(gains):
        axes[0].text(index, value + 18.0, f"{value:.0f}", ha="center", fontsize=9)
    axes[0].set_ylim(0.0, 1.18 * max(gains))

    density = np.asarray([item["electron_density_m3"] for item in inverse]) / 1.0e15
    response = np.asarray([
        item["required_ion_response_multiplier_relative_to_declared_bohm_flux"]
        for item in inverse])
    axes[1].plot(density, response, "o-", color="#386cb0", linewidth=2.0)
    axes[1].axhline(1.0, color="#777777", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel(r"Measured electron density ($10^{15}$ m$^{-3}$)")
    axes[1].set_ylabel("Required effective ion-response multiplier")
    axes[1].set_title(
        "Data inversion diagnoses missing boundary response\n"
        r"$\Gamma_\mathrm{effective}\propto n_e^{0.463}$ (not a prediction)")
    axes[1].grid(alpha=0.25)
    for x_value, y_value in zip(density, response):
        axes[1].annotate(
            f"{y_value:.2f}", (x_value, y_value),
            xytext=(0, 8), textcoords="offset points", ha="center", fontsize=9)
    figure.suptitle("Jeong 2023 reactor-to-feature closure audit", fontsize=14)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _condition_boundary(target, radicals):
    support_density = (
        1.9e15 if target.control_mode == "ion_energy"
        else target.electron_density_m3)
    selected = tuple(
        item for item in radicals
        if np.isclose(item.electron_density_m3, support_density))
    return build_jeong_2023_boundary_state(
        target, selected, reference_plane_m=3.0e-6,
        radical_channel_mode="huang_2019_reduced")


def _zero_d_depth_nm(
        target, radicals, energetic_response_scale, ion_flux_multiplier=1.0,
        max_step_s=1.0):
    boundary = _condition_boundary(target, radicals)
    neutral = {
        item.name: item.flux_m2_s for item in boundary.species
        if item.charge_number == 0}
    positive = tuple(
        item for item in boundary.species if item.charge_number > 0)
    energetic = tuple(EnergeticFlux(
        item.name,
        item.flux_m2_s * ion_flux_multiplier,
        [target.self_bias_magnitude_v],
        [1.0],
        [1.0],
    ) for item in positive)
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
            energetic_response_scale=energetic_response_scale))
    result = mechanism.advance(
        SiO2SurfaceState.bare(),
        SurfaceFluxes(neutral, energetic),
        target.etch_duration_s,
        max_step_s=max_step_s)
    return float(
        result.state.removed_formula_units_m2
        / mechanism.parameters.bulk_formula_density_m3 * 1.0e9)


def _calibrate_current_anchor(depths, radicals):
    anchor = next(item for item in depths if item.split == "calibration")
    scale = brentq(
        lambda value: (
            _zero_d_depth_nm(anchor, radicals, value) - anchor.etch_depth_nm),
        0.2, 4.0, xtol=2e-7)
    return anchor, float(scale), _zero_d_depth_nm(anchor, radicals, scale)


def _inverse_flux_response(depths, radicals, energetic_response_scale):
    targets = sorted(
        (item for item in depths
         if item.control_mode == "ion_flux"
         and np.isclose(item.trench_width_nm, 200.0)),
        key=lambda item: item.electron_density_m3)
    rows = []
    for target in targets:
        multiplier = brentq(
            lambda value: (
                _zero_d_depth_nm(
                    target, radicals, energetic_response_scale,
                    ion_flux_multiplier=value)
                - target.etch_depth_nm),
            0.05, 4.0, xtol=2e-7)
        rows.append({
            "electron_density_m3": float(target.electron_density_m3),
            "experimental_depth_nm": float(target.etch_depth_nm),
            "required_ion_response_multiplier_relative_to_declared_bohm_flux": (
                float(multiplier)),
            "replayed_depth_nm": _zero_d_depth_nm(
                target, radicals, energetic_response_scale,
                ion_flux_multiplier=multiplier),
        })
    density = np.asarray([item["electron_density_m3"] for item in rows])
    multiplier = np.asarray([
        item["required_ion_response_multiplier_relative_to_declared_bohm_flux"]
        for item in rows])
    effective_flux = density * multiplier
    beta, log_prefactor = np.polyfit(np.log(density), np.log(effective_flux), 1)
    low_multiplier = multiplier[0]
    for row, value in zip(rows, multiplier):
        row["shape_only_equivalent_mass_amu_if_low_density_is_40amu"] = float(
            40.0 * (low_multiplier / value) ** 2)
    return {
        "rows": rows,
        "effective_flux_power_law_exponent_beta": float(beta),
        "effective_flux_power_law_prefactor": float(np.exp(log_prefactor)),
        "development_data_consumed": True,
        "claim_limitation": (
            "This inversion uses all three 200 nm experimental flux-sweep depths. "
            "It diagnoses the missing response and is not a predictive closure."),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()

    depth_path = DATA / "digitized_figure7_depths.csv"
    radical_path = DATA / "digitized_figure6_radicals.csv"
    historical_path = (
        ROOT / "results" / "jeong_2023_predictive_validation" / "summary.json")
    feedback_path = (
        ROOT / "results" / "jeong_2023_surface_feedback_audit" / "audit.json")
    sheath_path = (
        ROOT / "results" / "jeong_2023_virtual_sheath_audit" / "audit.json")
    depths = load_jeong_2023_etch_depths(depth_path)
    radicals = load_jeong_2023_radical_densities(radical_path)
    historical = json.loads(historical_path.read_text(encoding="utf-8"))
    feedback = json.loads(feedback_path.read_text(encoding="utf-8"))
    sheath = json.loads(sheath_path.read_text(encoding="utf-8"))

    anchor, current_scale, anchor_depth = _calibrate_current_anchor(depths, radicals)
    inverse = _inverse_flux_response(depths, radicals, current_scale)
    direct_rows = [
        item for item in feedback["endpoint_gains"]
        if item["model"] == "published_fc_ion_direct_deposition"]
    baseline_gain = next(
        item["predicted_zero_d_endpoint_gain_nm"] for item in direct_rows
        if item["synthetic_low_energy_fraction"] == 0.0)
    strongest_direct = min(
        direct_rows, key=lambda item: item["predicted_zero_d_endpoint_gain_nm"])
    experimental_gain = strongest_direct["experimental_feature_endpoint_gain_nm"]
    required_reduction = baseline_gain - experimental_gain
    supplied_reduction = (
        baseline_gain - strongest_direct["predicted_zero_d_endpoint_gain_nm"])

    payload = {
        "campaign": "jeong_2023_reactor_to_feature_closure_audit",
        "status": "complete",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "scope": (
            "bounded zero-D mechanism and boundary diagnosis; no moving-profile "
            "claim and no held-out validation"),
        "inputs": {
            "checksums_sha256": {
                str(path.relative_to(ROOT)): _sha256(path)
                for path in (
                    depth_path, radical_path, historical_path,
                    feedback_path, sheath_path)},
            "jeong_gas_flow_sccm": {"C4F8": 80.0, "Ar": 20.0},
            "pressure_mTorr": 20.0,
            "etch_duration_s": 1200.0,
        },
        "historical_profile_campaign": {
            "status": historical["status"],
            "verdict": historical["verdict"],
            "current_operator": historical["operator_replay"]["current_operator"],
            "historical_energy_endpoint_gain_nm": historical["trend_checks"][
                "ion_energy"]["predicted_endpoint_gain_nm"],
            "experimental_energy_endpoint_gain_nm": historical["trend_checks"][
                "ion_energy"]["experimental_endpoint_gain_nm"],
            "historical_flux_endpoint_gain_nm": historical["trend_checks"][
                "ion_flux"]["predicted_endpoint_gain_nm"],
            "experimental_flux_endpoint_gain_nm": historical["trend_checks"][
                "ion_flux"]["experimental_endpoint_gain_nm"],
        },
        "current_operator_anchor_calibration": {
            "calibration_point": {
                "control_mode": anchor.control_mode,
                "trench_width_nm": anchor.trench_width_nm,
                "self_bias_magnitude_v": anchor.self_bias_magnitude_v,
                "experimental_depth_nm": anchor.etch_depth_nm,
            },
            "zero_d_energetic_response_scale": current_scale,
            "zero_d_replayed_depth_nm": anchor_depth,
            "profile_scale_is_not_yet_frozen": True,
            "reason": (
                "The current source-faithful operator removes the historical implicit "
                "substrate polymer-nucleation projection; a moving-profile anchor has "
                "not passed timestep/remap refinement."),
        },
        "bounded_source_backed_tests": {
            "collisionless_virtual_sheath": {
                "maximum_density_induced_yield_flattening_fraction": max(
                    item["density_induced_yield_flattening_fraction"]
                    for item in sheath["amplitude_summary"].values()),
                "maximum_5_to_70_eV_fraction": sheath["decision"][
                    "maximum_5_to_70_eV_fraction"],
                "sufficient": sheath["decision"][
                    "collisionless_iedf_is_sufficient_flux_slope_closure"],
            },
            "low_energy_fc_ion_direct_polymer_deposition": {
                "synthetic_fraction_maximum": strongest_direct[
                    "synthetic_low_energy_fraction"],
                "baseline_endpoint_gain_nm": baseline_gain,
                "strongest_endpoint_gain_nm": strongest_direct[
                    "predicted_zero_d_endpoint_gain_nm"],
                "experimental_endpoint_gain_nm": experimental_gain,
                "fraction_of_required_gain_reduction_supplied": float(
                    supplied_reduction / required_reduction),
                "sufficient": feedback["decision"][
                    "fc_ion_direct_deposition_closes_feature_slope"],
            },
        },
        "inverse_missing_response_diagnostic": inverse,
        "scientific_findings": [
            (
                "The Jeong paper reports electron density and self-bias, not a measured "
                "species-resolved ion flux, IEAD, or complete radical wall-flux vector."),
            (
                "The historical all-Ar adapter is not a reactor model. The common engine "
                "now accepts explicit multi-ion density mixtures or reactor-supplied "
                "species fluxes without changing transport."),
            (
                "Neither collisionless sheath broadening nor the omitted published "
                "low-energy fluorocarbon-ion deposition channel has enough bounded "
                "leverage to close the 200 nm density-sweep slope."),
            (
                "The experimental depths can be reproduced only by a density-dependent "
                "effective response inferred from those same depths; that is diagnostic "
                "calibration, not independent prediction."),
        ],
        "decision": {
            "moving_profile_matrix_earned": False,
            "strict_jeong_flux_validation_closed": False,
            "blocking_boundary_evidence": [
                "species-resolved positive-ion flux versus discharge condition",
                "species-resolved ion energy-angle distributions or validated reactor/sheath output",
                "complete radical wall-flux vector rather than selected volume densities",
                "hot-neutral production and energy-angle distributions after surface neutralization",
            ],
            "next_engine_action": (
                "Retain the new species-resolved boundary contract and direct FC-ion "
                "deposition reaction. Do not fit the held-out Jeong flux sweep or spend "
                "moving-profile compute until an evidenced reactor boundary is supplied. "
                "Proceed to benchmarks whose source boundary is sufficiently specified."),
        },
        "primary_sources": {
            "jeong_2023": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10222222/",
            "huang_kushner_2019": (
                "https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf"),
            "li_oehrlein_kushner_2004": (
                "https://cpseg.eecs.umich.edu/pub/articles/jvsta_22_500_2004.pdf"),
            "kim_2021_ion_mass_energy": "https://doi.org/10.3390/coatings11080993",
        },
    }
    _atomic_json(args.output, payload)
    plot_path = args.plot or args.output.with_name("closure_audit.png")
    _write_plot(plot_path, payload)
    print(json.dumps({
        "output": str(args.output),
        "plot": str(plot_path),
        "current_zero_d_anchor_scale": current_scale,
        "required_effective_flux_beta": inverse[
            "effective_flux_power_law_exponent_beta"],
        "moving_profile_matrix_earned": payload["decision"][
            "moving_profile_matrix_earned"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
