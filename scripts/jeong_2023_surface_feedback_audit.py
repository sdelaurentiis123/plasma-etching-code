#!/usr/bin/env python3
"""Bounded zero-D audit of the Huang--Kushner activated-site feedback.

This script does not fit Jeong's held-out flux sweep.  It reuses the already frozen 200 nm
boundary artifacts, preserves total ion flux, and asks two narrower questions:

1. Does the published activated-site mechanism change anything under the declared 740 eV
   monoenergy Jeong closure?
2. If a declared fraction of that same ion flux is moved to a 15 eV diagnostic population,
   does the source-backed feedback act in the expected direction and with numerically refined
   behavior?

The low-energy fractions are sensitivity coordinates, not inferred IEAD components.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from dataclasses import replace
from pathlib import Path

import numpy as np

from petch.surface_kinetics import (
    EnergeticFlux,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
    SurfaceFluxes,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "jeong_2023_predictive_validation"
DEFAULT_OUTPUT = ROOT / "results" / "jeong_2023_surface_feedback_audit" / "audit.json"
INPUT_FILES = (
    "flux_1p1e15_width200_medium.json",
    "flux_1p9e15_width200_medium.json",
    "flux_3p1e15_width200_medium.json",
)
LOW_ENERGY_FRACTIONS = (0.0, 0.01, 0.05, 0.10, 0.20)
MAX_STEP_LEVELS_S = (0.5, 0.25)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_revision() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        capture_output=True, text=True).stdout.strip()


def _legacy_projection(parameters):
    """Reconstruct the prior implicit-activation projection for a paired comparison."""
    return replace(
        parameters,
        polymer_deposition_probability_on_substrate=dict(
            parameters.activated_polymer_deposition_probability_on_substrate),
        activated_polymer_deposition_probability_on_substrate={},
        activated_polymer_deposition_probability_on_polymer={},
        complex_activation_yield=None,
        polymer_activation_yield=None,
        activation_energetic_species=(),
    )


def _surface_fluxes(boundary, low_energy_fraction, low_energy_species):
    ion_flux = float(boundary["ion_flux_m2_s"])
    fraction = float(low_energy_fraction)
    energetic = []
    if fraction < 1.0:
        energetic.append(EnergeticFlux(
            "Ar+", ion_flux * (1.0 - fraction), [740.0], [1.0], [1.0]))
    if fraction > 0.0:
        energetic.append(EnergeticFlux(
            low_energy_species, ion_flux * fraction, [15.0], [1.0], [1.0]))
    return SurfaceFluxes(
        boundary["radical_channel_flux_m2_s"], tuple(energetic))


def _run(
        mechanism, boundary, duration_s, max_step_s, low_energy_fraction,
        low_energy_species="Ar+"):
    result = mechanism.advance(
        SiO2SurfaceState.bare(),
        _surface_fluxes(boundary, low_energy_fraction, low_energy_species),
        float(duration_s), max_step_s=float(max_step_s))
    state = result.state
    return {
        "removed_depth_nm": float(
            state.removed_formula_units_m2
            / mechanism.parameters.bulk_formula_density_m3 * 1.0e9),
        "polymer_units_m2": float(state.polymer_units_m2),
        "polymer_monolayers": float(
            state.polymer_units_m2
            / mechanism.parameters.polymer_monolayer_density_m2),
        "complex_fraction": float(state.complex_fraction),
        "activated_complex_fraction": float(state.activated_complex_fraction),
        "activated_polymer_fraction": float(state.activated_polymer_fraction),
        "material_ledger_residual_units_m2": float(
            result.material_exchange.residual_units_m2("SiO2_formula_unit")),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--energetic-response-scale", type=float, default=1.369)
    args = parser.parse_args()

    parameters = ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
        energetic_response_scale=args.energetic_response_scale)
    activated = ReducedSiO2FluorocarbonMechanism(parameters)
    legacy = ReducedSiO2FluorocarbonMechanism(_legacy_projection(parameters))
    source_path = ROOT / "src" / "petch" / "surface_kinetics.py"

    inputs = []
    for filename in INPUT_FILES:
        path = INPUT / filename
        artifact = json.loads(path.read_text())
        run = artifact["runs"][0]
        inputs.append({
            "filename": filename,
            "sha256": _sha256(path),
            "electron_density_m3": float(run["target"]["electron_density_m3"]),
            "experimental_depth_nm": float(run["target"]["etch_depth_nm"]),
            "etch_duration_s": float(run["target"]["etch_duration_s"]),
            "boundary": run["boundary"],
        })

    records = []
    for max_step_s in MAX_STEP_LEVELS_S:
        for item in inputs:
            common = {
                "electron_density_m3": item["electron_density_m3"],
                "experimental_depth_nm": item["experimental_depth_nm"],
                "max_step_s": max_step_s,
            }
            records.append({
                **common,
                "model": "legacy_implicit_activated_complex_projection",
                "synthetic_low_energy_fraction": None,
                "synthetic_low_energy_species": None,
                **_run(
                    legacy, item["boundary"], item["etch_duration_s"],
                    max_step_s, 0.0),
            })
            for fraction in LOW_ENERGY_FRACTIONS:
                for model, species in (
                        ("published_explicit_activation", "Ar+"),
                        ("published_fc_ion_direct_deposition", "CF+")):
                    records.append({
                        **common,
                        "model": model,
                        "synthetic_low_energy_fraction": fraction,
                        "synthetic_low_energy_species": species,
                        **_run(
                            activated, item["boundary"], item["etch_duration_s"],
                            max_step_s, fraction, species),
                    })

    fine = [record for record in records if record["max_step_s"] == min(MAX_STEP_LEVELS_S)]
    coarse = [record for record in records if record["max_step_s"] == max(MAX_STEP_LEVELS_S)]
    coarse_by_key = {
        (record["electron_density_m3"], record["model"],
         record["synthetic_low_energy_fraction"],
         record["synthetic_low_energy_species"]): record
        for record in coarse
    }
    refinement = []
    for record in fine:
        key = (record["electron_density_m3"], record["model"],
               record["synthetic_low_energy_fraction"],
               record["synthetic_low_energy_species"])
        earlier = coarse_by_key[key]
        refinement.append({
            "electron_density_m3": key[0],
            "model": key[1],
            "synthetic_low_energy_fraction": key[2],
            "synthetic_low_energy_species": key[3],
            "depth_relative_change": abs(
                record["removed_depth_nm"] - earlier["removed_depth_nm"])
            / max(abs(record["removed_depth_nm"]), np.finfo(float).tiny),
            "polymer_monolayer_absolute_change": abs(
                record["polymer_monolayers"] - earlier["polymer_monolayers"]),
        })

    experiment_gain = inputs[-1]["experimental_depth_nm"] - inputs[0][
        "experimental_depth_nm"]
    endpoint_gains = []
    for model in (
            "published_explicit_activation",
            "published_fc_ion_direct_deposition"):
        for fraction in LOW_ENERGY_FRACTIONS:
            selected = sorted(
                (record for record in fine
                 if record["model"] == model
                 and record["synthetic_low_energy_fraction"] == fraction),
                key=lambda record: record["electron_density_m3"])
            endpoint_gains.append({
                "model": model,
                "synthetic_low_energy_fraction": fraction,
                "synthetic_low_energy_species": selected[0][
                    "synthetic_low_energy_species"],
                "predicted_zero_d_endpoint_gain_nm": (
                    selected[-1]["removed_depth_nm"] - selected[0]["removed_depth_nm"]),
                "experimental_feature_endpoint_gain_nm": experiment_gain,
            })

    output = {
        "campaign": "jeong_2023_surface_feedback_zero_d_audit",
        "status": "complete",
        "scope": (
            "zero-D source-mechanism response; sensitivity only; not a profile prediction"),
        "git_revision": _git_revision(),
        "surface_kinetics_sha256": _sha256(source_path),
        "inputs": [{key: value for key, value in item.items() if key != "boundary"}
                   for item in inputs],
        "configuration": {
            "energetic_response_scale_frozen": args.energetic_response_scale,
            "high_energy_eV": 740.0,
            "diagnostic_low_energy_eV": 15.0,
            "synthetic_low_energy_fractions": list(LOW_ENERGY_FRACTIONS),
            "synthetic_low_energy_species": ["Ar+", "CF+"],
            "max_step_levels_s": list(MAX_STEP_LEVELS_S),
            "initial_surface": "bare",
            "low_energy_fraction_semantics": (
                "unfitted sensitivity coordinate preserving total ion flux; not an IEAD claim"),
        },
        "source_contract": {
            "complex_activation_window_eV": [5.0, 70.0],
            "polymer_activation_window_eV": [5.0, 30.0],
            "primary_source": "https://doi.org/10.1116/1.5090606",
            "jeong_boundary_missing": [
                "ion_energy_distribution", "hot_neutral_energy_angle_distribution",
                "species_resolved_fluorocarbon_ion_fluxes", "atomic_F"],
        },
        "records": records,
        "refinement": {
            "records": refinement,
            "maximum_depth_relative_change": max(
                item["depth_relative_change"] for item in refinement),
            "maximum_polymer_monolayer_absolute_change": max(
                item["polymer_monolayer_absolute_change"] for item in refinement),
        },
        "endpoint_gains": endpoint_gains,
        "decision": {
            "declared_740ev_closure_activates_feedback": False,
            "low_energy_sensitivity_has_expected_sign": {
                model: all(
                    selected[index + 1]["predicted_zero_d_endpoint_gain_nm"]
                    < selected[index]["predicted_zero_d_endpoint_gain_nm"]
                    for index in range(len(selected) - 1))
                for model in (
                    "published_explicit_activation",
                    "published_fc_ion_direct_deposition")
                for selected in [[
                    item for item in endpoint_gains if item["model"] == model]]
            },
            "fc_ion_direct_deposition_closes_feature_slope": min(
                item["predicted_zero_d_endpoint_gain_nm"]
                for item in endpoint_gains
                if item["model"] == "published_fc_ion_direct_deposition")
                <= experiment_gain,
            "moving_profile_rerun_earned": False,
            "reason": (
                "The published activation and direct fluorocarbon-ion deposition channels are "
                "real, but the declared Jeong monoenergy boundary cannot populate them. The "
                "bounded synthetic fractions determine whether even a generous low-energy "
                "population has enough leverage before a profile rerun. A reactor-resolved "
                "species IEAD/hot-neutral boundary is still required; do not fit the synthetic "
                "fraction to held-out depths."),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "maximum_depth_relative_change": output["refinement"][
            "maximum_depth_relative_change"],
        "endpoint_gains": endpoint_gains,
        "decision": output["decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
