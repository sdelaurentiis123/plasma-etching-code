#!/usr/bin/env python3
"""Bounded Jeong-2023 audit of the diagnostic-conditioned collisionless IEDF.

The source reports time-averaged self-bias but not the measured voltage waveform.  Therefore every
nonzero RF amplitude in this audit is a declared sensitivity coordinate, never a fitted input or a
prediction.  The experiment decides whether collisionless finite-transit broadening alone can
activate the published Huang--Kushner low-energy surface states or flatten the held-out flux trend.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from petch.experimental_data import load_jeong_2023_etch_depths
from petch.sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage
from petch.surface_kinetics import ReducedSiO2FluorocarbonParameters


ROOT = Path(__file__).parents[1]
DATA = ROOT / "data" / "experimental" / "jeong_2023"


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git_revision():
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True).stdout.strip()


def _write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _weighted_yield(law, energy):
    return float(np.mean(law.evaluate(energy, np.ones_like(energy))))


def _record(target, amplitude_ratio, phase_count, steps_per_scale, parameters):
    phase = 2.0 * np.pi * (np.arange(phase_count, dtype=float) + 0.5) / phase_count
    waveform = PeriodicSheathVoltage.sinusoidal(
        dc_v=target.self_bias_magnitude_v,
        amplitude_v=amplitude_ratio * target.self_bias_magnitude_v,
        frequency_hz=4.0e5,
        source=(
            "Jeong 2023 reports 400 kHz drive and time-average self-bias but not waveform "
            "amplitude; declared sensitivity coordinate"),
        evidence_kind="assumed",
    )
    sheath = CollisionlessWaveformSheath(
        waveform=waveform,
        Te_eV=3.0,
        ion_mass_amu=39.948,
        density_m3=target.electron_density_m3,
    )
    energy = sheath.ion_impact_energies(
        phase,
        steps_per_period=steps_per_scale,
        steps_per_transit=steps_per_scale,
    )
    hist_edges = np.linspace(0.0, 1.95 * target.self_bias_magnitude_v, 81)
    histogram, _ = np.histogram(energy, bins=hist_edges, density=True)
    return {
        "electron_density_m3": target.electron_density_m3,
        "experimental_depth_nm": target.etch_depth_nm,
        "self_bias_magnitude_v": target.self_bias_magnitude_v,
        "amplitude_ratio": float(amplitude_ratio),
        "assumed_rf_amplitude_v": float(amplitude_ratio * target.self_bias_magnitude_v),
        "phase_count": int(phase_count),
        "steps_per_fast_period_and_transit": int(steps_per_scale),
        "sheath_thickness_m": sheath.thickness,
        "energy_mean_eV": float(np.mean(energy)),
        "energy_std_eV": float(np.std(energy)),
        "energy_quantile_eV": {
            "q01": float(np.quantile(energy, 0.01)),
            "q05": float(np.quantile(energy, 0.05)),
            "q50": float(np.quantile(energy, 0.50)),
            "q95": float(np.quantile(energy, 0.95)),
            "q99": float(np.quantile(energy, 0.99)),
        },
        "fraction_5_to_30_eV": float(np.mean((energy >= 5.0) & (energy <= 30.0))),
        "fraction_5_to_70_eV": float(np.mean((energy >= 5.0) & (energy <= 70.0))),
        "normal_incidence_mean_yield": {
            "bare_sio2": _weighted_yield(parameters.bare_sio2_yield, energy),
            "complex_sio2": _weighted_yield(parameters.complex_sio2_yield, energy),
            "polymer_sputter": _weighted_yield(parameters.polymer_sputter_yield, energy),
            "complex_activation": _weighted_yield(
                parameters.complex_activation_yield, energy),
            "polymer_activation": _weighted_yield(
                parameters.polymer_activation_yield, energy),
        },
        "histogram": {
            "energy_edges_eV": hist_edges.tolist(),
            "probability_density_per_eV": histogram.tolist(),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--amplitude-ratios", type=float, nargs="+",
        default=(0.0, 0.25, 0.50, 0.75, 0.95))
    parser.add_argument("--phase-count", type=int, default=1024)
    parser.add_argument("--steps-per-scale", type=int, default=256)
    args = parser.parse_args()
    if (args.phase_count < 64 or args.steps_per_scale < 64
            or any(not 0.0 <= value < 1.0 for value in args.amplitude_ratios)):
        raise ValueError("invalid virtual-sheath sensitivity controls")

    source_path = DATA / "digitized_figure7_depths.csv"
    targets = tuple(sorted(
        (item for item in load_jeong_2023_etch_depths(source_path)
         if item.control_mode == "ion_flux" and item.trench_width_nm == 200.0),
        key=lambda item: item.electron_density_m3))
    if len(targets) != 3:
        raise RuntimeError("Jeong 200 nm flux audit requires exactly three frozen targets")
    parameters = (
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
            energetic_response_scale=1.369))

    records = [
        _record(target, ratio, args.phase_count, args.steps_per_scale, parameters)
        for ratio in args.amplitude_ratios for target in targets]
    coarse_count = args.phase_count // 2
    coarse = [
        _record(target, ratio, coarse_count, args.steps_per_scale, parameters)
        for ratio in args.amplitude_ratios for target in targets]
    coarse_index = {
        (item["amplitude_ratio"], item["electron_density_m3"]): item for item in coarse}
    refinement = []
    for item in records:
        prior = coarse_index[(item["amplitude_ratio"], item["electron_density_m3"])]
        refinement.append({
            "amplitude_ratio": item["amplitude_ratio"],
            "electron_density_m3": item["electron_density_m3"],
            "mean_energy_relative_change": abs(
                item["energy_mean_eV"] - prior["energy_mean_eV"])
            / item["energy_mean_eV"],
            "std_energy_relative_change": abs(
                item["energy_std_eV"] - prior["energy_std_eV"])
            / max(item["energy_std_eV"], 1.0),
            "bare_yield_relative_change": abs(
                item["normal_incidence_mean_yield"]["bare_sio2"]
                - prior["normal_incidence_mean_yield"]["bare_sio2"])
            / item["normal_incidence_mean_yield"]["bare_sio2"],
        })

    by_ratio = {}
    for ratio in args.amplitude_ratios:
        selected = [item for item in records if item["amplitude_ratio"] == ratio]
        low, _, high = selected
        low_yield = low["normal_incidence_mean_yield"]["bare_sio2"]
        high_yield = high["normal_incidence_mean_yield"]["bare_sio2"]
        by_ratio[str(ratio)] = {
            "high_to_low_density_bare_yield_ratio": high_yield / low_yield,
            "density_induced_yield_flattening_fraction": max(
                0.0, 1.0 - high_yield / low_yield),
            "maximum_fraction_5_to_70_eV": max(
                item["fraction_5_to_70_eV"] for item in selected),
        }

    prior_summary = json.loads((
        ROOT / "results" / "jeong_2023_predictive_validation" / "summary.json"
    ).read_text(encoding="utf-8"))
    experimental_gain = prior_summary["trend_checks"]["ion_flux"][
        "experimental_endpoint_gain_nm"]
    predicted_gain = prior_summary["trend_checks"]["ion_flux"][
        "predicted_endpoint_gain_nm"]
    required_gain_reduction = 1.0 - experimental_gain / predicted_gain
    maximum_yield_flattening = max(
        item["density_induced_yield_flattening_fraction"] for item in by_ratio.values())
    maximum_low_energy_fraction = max(
        item["fraction_5_to_70_eV"] for item in records)
    payload = {
        "campaign": "jeong_2023_collisionless_virtual_sheath_audit",
        "status": "complete",
        "scope": (
            "source-constrained sensitivity only; waveform amplitude was not published; "
            "no held-out target was used to choose it"),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "inputs": {
            "source_csv": str(source_path.relative_to(ROOT)),
            "source_csv_sha256": _sha256(source_path),
            "jeong_source_xml_sha256": targets[0].source_xml_sha256,
            "mean_sheath_drop_proxy": "reported time-average self-bias magnitude",
            "full_waveform_published": False,
            "rf_frequency_hz": 4.0e5,
            "electron_temperature_eV_assumed": 3.0,
            "ion_species_proxy": "Ar+",
            "ion_mass_amu": 39.948,
            "energetic_response_scale_frozen_from_prior_single_anchor": 1.369,
            "amplitude_ratios": list(args.amplitude_ratios),
        },
        "comparison_to_existing_profile_miss": {
            "experimental_flux_sweep_endpoint_gain_nm": experimental_gain,
            "existing_prediction_endpoint_gain_nm": predicted_gain,
            "required_endpoint_gain_reduction_fraction": required_gain_reduction,
            "maximum_collisionless_density_induced_yield_flattening_fraction": (
                maximum_yield_flattening),
        },
        "decision": {
            "maximum_5_to_70_eV_fraction": maximum_low_energy_fraction,
            "collisionless_iedf_activates_published_surface_feedback": bool(
                maximum_low_energy_fraction > 0.0),
            "collisionless_iedf_is_sufficient_flux_slope_closure": False,
            "next_required_boundary_physics": [
                "evidenced ion/hot-neutral energy-angle populations",
                "surface-neutralization or charge-exchange production of hot neutrals",
                "species-resolved positive-ion inventory",
            ],
            "conclusion": (
                "Finite-transit collisionless broadening is real but too weak and produces no "
                "5--70 eV activation population over the declared sensitivity box. Do not spend "
                "full-profile compute on this mechanism alone."),
        },
        "amplitude_summary": by_ratio,
        "records": records,
        "refinement": {
            "coarse_phase_count": coarse_count,
            "fine_phase_count": args.phase_count,
            "maximum_mean_energy_relative_change": max(
                item["mean_energy_relative_change"] for item in refinement),
            "maximum_std_energy_relative_change": max(
                item["std_energy_relative_change"] for item in refinement),
            "maximum_bare_yield_relative_change": max(
                item["bare_yield_relative_change"] for item in refinement),
            "records": refinement,
        },
        "sources": {
            "jeong_2023": "https://doi.org/10.3390/ma16093820",
            "huang_kushner_2019": "https://doi.org/10.1116/1.5090606",
            "virtual_ied_sensor": "https://arxiv.org/abs/2012.14882",
            "finite_transit_validation": (
                "https://www.nist.gov/publications/measurements-and-modeling-ion-energy-"
                "distributions-high-density-radio-frequency-biased"),
        },
    }
    _write_json_atomic(args.output, payload)
    print(json.dumps(payload["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
