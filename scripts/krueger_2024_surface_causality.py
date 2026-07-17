#!/usr/bin/env python3
"""Bounded flat-surface causality audit for the reduced Krüger 2024 replay.

This is not profile validation. It checks whether the provenance-locked published flux/IEAD
boundary and reduced oxide/mask mechanisms have the source-required directions before a trench
run is allowed: oxygen removes passivation, energetic bombardment removes film and oxide, the
mask remains protected by its film, and all material ledgers close.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np

from petch.amorphous_carbon_mask import (
    AmorphousCarbonMaskMechanism,
    AmorphousCarbonMaskParameters,
)
from petch.reactor_boundary import (
    load_krueger_2024_digitized_iead,
    load_krueger_2024_reactor_flux_deck,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SurfaceFluxes,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "krueger_2024"


def _git_revision():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def _write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _fluxes(neutral, iead, ion_flux, *, oxygen_multiplier=1.0, energy_scale=1.0):
    selected = dict(neutral)
    selected["O"] *= float(oxygen_multiplier)
    ion = EnergeticFlux(
        "ions", ion_flux, iead.energy_eV * float(energy_scale),
        np.cos(np.deg2rad(iead.signed_angle_deg)), iead.probability_weight)
    return SurfaceFluxes(selected, (ion,))


def _oxide_record(fluxes, duration_s, max_step_s):
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.krueger_2024_reduced_projection())
    result = mechanism.advance(
        mechanism.initial_state(), fluxes, duration_s, max_step_s=max_step_s)
    removed = np.asarray(
        result.material_exchange.removed_units_m2["SiO2_formula_unit"])
    film_removed = np.asarray(
        result.material_exchange.removed_units_m2["fluorocarbon_film_unit"])
    return {
        "etch_depth_nm": float(
            result.state.removed_formula_units_m2
            / mechanism.parameters.bulk_formula_density_m3 * 1e9),
        "polymer_inventory_monolayers": float(
            result.state.polymer_units_m2
            / mechanism.parameters.polymer_monolayer_density_m2),
        "complex_fraction": float(result.state.complex_fraction),
        "formed_complex_units_m2": float(result.formed_complex_units_m2),
        "removed_complex_units_m2": float(result.removed_complex_units_m2),
        "removed_bare_formula_units_m2": float(
            result.removed_bare_formula_units_m2),
        "maximum_ledger_residual_units_m2": float(max(
            np.max(np.abs(
                result.material_exchange.residual_units_m2("SiO2_formula_unit"))),
            np.max(np.abs(
                result.material_exchange.residual_units_m2(
                    "fluorocarbon_film_unit"))),
        )),
        "removed_sio2_formula_units_m2": float(removed),
        "removed_film_units_m2": float(film_removed),
    }


def _mask_record(fluxes, duration_s):
    mechanism = AmorphousCarbonMaskMechanism(
        AmorphousCarbonMaskParameters.krueger_2024_reduced_projection())
    result = mechanism.advance(mechanism.initial_state(), fluxes, duration_s)
    net_velocity = result.etch_velocity_m_s - result.normal_growth_velocity_m_s
    return {
        "polymer_inventory_monolayers": float(
            result.state.polymer_units_m2
            / mechanism.parameters.polymer_monolayer_density_m2),
        "carbon_recession_nm": float(
            result.removed_carbon_atoms_m2
            / mechanism.parameters.bulk_carbon_atom_density_m3 * 1e9),
        "gross_film_growth_nm": float(
            result.deposited_polymer_units_m2
            / mechanism.parameters.polymer_unit_density_m3 * 1e9),
        "gross_film_removal_nm": float(
            result.removed_polymer_units_m2
            / mechanism.parameters.polymer_unit_density_m3 * 1e9),
        "net_geometry_recession_nm": float(net_velocity * duration_s * 1e9),
        "maximum_ledger_residual_units_m2": float(max(
            np.max(np.abs(
                result.material_exchange.residual_units_m2(
                    "amorphous_carbon_atom"))),
            np.max(np.abs(
                result.material_exchange.residual_units_m2(
                    "fluorocarbon_film_unit"))),
        )),
    }


def _condition(neutral, iead, ion_flux, duration_s, oxygen, energy, max_step):
    fluxes = _fluxes(
        neutral, iead, ion_flux,
        oxygen_multiplier=oxygen, energy_scale=energy)
    return {
        "oxygen_flux_multiplier": float(oxygen),
        "ion_energy_multiplier": float(energy),
        "oxide": _oxide_record(fluxes, duration_s, max_step),
        "mask": _mask_record(fluxes, duration_s),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, default=60.0)
    parser.add_argument("--oxide-max-step-s", type=float, default=0.05)
    args = parser.parse_args()
    if (not np.isfinite(args.duration_s) or args.duration_s <= 0.0
            or not np.isfinite(args.oxide_max_step_s)
            or args.oxide_max_step_s <= 0.0):
        raise ValueError("invalid Krüger flat-audit integration controls")

    deck = load_krueger_2024_reactor_flux_deck(DATA)
    iead = load_krueger_2024_digitized_iead(DATA)
    neutral = {
        item.name: item.flux_m2_s
        for item in deck.species_fluxes if item.role == "neutral"}
    ion_flux = deck.get("ions").flux_m2_s
    base = _condition(
        neutral, iead, ion_flux, args.duration_s,
        1.0, 1.0, args.oxide_max_step_s)
    oxygen_sweep = [
        _condition(
            neutral, iead, ion_flux, args.duration_s,
            multiplier, 1.0, args.oxide_max_step_s)
        for multiplier in (0.5, 1.0, 2.0)
    ]
    energy_sweep = [
        _condition(
            neutral, iead, ion_flux, args.duration_s,
            1.0, multiplier, args.oxide_max_step_s)
        for multiplier in (0.5, 0.75, 1.0, 1.25)
    ]
    refinement = [
        {
            "oxide_max_step_s": step,
            "oxide": _oxide_record(
                _fluxes(neutral, iead, ion_flux),
                args.duration_s, step),
        }
        for step in (0.2, 0.1, 0.05)
    ]
    coarse, medium, fine = refinement
    depth_refinement_relative = abs(
        medium["oxide"]["etch_depth_nm"] - fine["oxide"]["etch_depth_nm"]
    ) / fine["oxide"]["etch_depth_nm"]
    oxygen_polymer = [
        item["oxide"]["polymer_inventory_monolayers"] for item in oxygen_sweep]
    oxygen_depth = [item["oxide"]["etch_depth_nm"] for item in oxygen_sweep]
    energy_depth = [item["oxide"]["etch_depth_nm"] for item in energy_sweep]
    energy_mask_film = [
        item["mask"]["polymer_inventory_monolayers"] for item in energy_sweep]
    maximum_ledger_residual = max(
        item[material]["maximum_ledger_residual_units_m2"]
        for item in oxygen_sweep + energy_sweep
        for material in ("oxide", "mask"))
    experimental_depth_nm = 825.0
    payload = {
        "campaign": "krueger_2024_reduced_surface_causality",
        "status": "complete",
        "claim": "development_replay_only",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "inputs": {
            "duration_s": args.duration_s,
            "published_flux_deck_sha256": deck.source_sha256,
            "digitized_iead_table_sha256": iead.table_sha256,
            "digitization_metadata_sha256": iead.metadata_sha256,
            "digitized_mean_energy_eV": iead.mean_energy_eV,
            "aggregate_ion_flux_m2_s": ion_flux,
            "neutral_flux_m2_s": neutral,
            "oxide_max_step_s": args.oxide_max_step_s,
        },
        "base_flat_surface": base,
        "experimental_context": {
            "base_trench_depth_nm": experimental_depth_nm,
            "flat_depth_is_not_profile_prediction": True,
            "trench_to_flat_depth_ratio_required": (
                experimental_depth_nm / base["oxide"]["etch_depth_nm"]),
            "initial_mask_height_nm": 850.0,
        },
        "oxygen_sweep": oxygen_sweep,
        "ion_energy_sensitivity": energy_sweep,
        "timestep_refinement": {
            "records": refinement,
            "medium_to_fine_depth_relative_change": depth_refinement_relative,
        },
        "gates": {
            "oxygen_reduces_oxide_polymer": bool(
                oxygen_polymer[0] > oxygen_polymer[1] > oxygen_polymer[2]),
            "oxygen_increases_oxide_removal": bool(
                oxygen_depth[0] < oxygen_depth[1] < oxygen_depth[2]),
            "ion_energy_increases_oxide_removal": bool(
                all(left < right for left, right in zip(
                    energy_depth[:-1], energy_depth[1:]))),
            "ion_energy_reduces_mask_film": bool(
                all(left > right for left, right in zip(
                    energy_mask_film[:-1], energy_mask_film[1:]))),
            "oxide_timestep_refined_below_1e-5": bool(
                depth_refinement_relative < 1e-5),
            "material_ledgers_roundoff_closed": bool(
                maximum_ledger_residual == 0.0),
        },
        "maximum_material_ledger_residual_units_m2": maximum_ledger_residual,
        "decision": {
            "short_trench_pilot_earned": bool(
                oxygen_polymer[0] > oxygen_polymer[1] > oxygen_polymer[2]
                and oxygen_depth[0] < oxygen_depth[1] < oxygen_depth[2]
                and all(left < right for left, right in zip(
                    energy_depth[:-1], energy_depth[1:]))
                and depth_refinement_relative < 1e-5
                and maximum_ledger_residual == 0.0),
            "predictive_validation_earned": False,
            "next_action": (
                "Run one coarse, wall-clock-bounded 90 nm trench pilot with the exact "
                "digitized IEAD; compare depth, opening, mask height, and profile causality. "
                "Do not tune held-out oxygen or power sweeps."),
        },
    }
    _write_json_atomic(args.output, payload)


if __name__ == "__main__":
    main()
