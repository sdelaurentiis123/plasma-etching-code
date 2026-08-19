#!/usr/bin/env python3
"""Materialize the conserved Oxford reactor state as a feature boundary.

This audit closes a software interface, not an experimental evidence gap.  It
preserves every positive-ion and thermal-neutral identity carried by the
reactor, the central-optic ion flux partition, ion charge and mass, and an
explicit deterministic energy-angle measure.  The target tool's IEAD, neutral
radial transfer, and TiO2/Cr surface probabilities remain evidence gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.iadf_two_component import kim_2025_reference_iadf
from petch.reactor_global import zhu_reactor_species
from petch.species_resolved_feature_boundary import (
    build_species_resolved_feature_boundary,
)


COUPLED_STATE = (
    ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
    / "central_276V.json"
)
AXISYMMETRIC_STATE = (
    ROOT / "results" / "curated" / "zhu_npg80_sheath_coupled_v1"
    / "axisymmetric_audit.json"
)
OUTPUT_DIR = (
    ROOT / "results" / "curated"
    / "zhu_npg80_species_resolved_feature_boundary_v1"
)
OUTPUT = OUTPUT_DIR / "audit.json"
CONDITION_ID = "zhu-2026-npg80-tio2-chf3-sf6-o2-20min"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _boundary_hash(boundary) -> str:
    digest = hashlib.sha256()
    digest.update(np.float64(boundary.reference_plane_m).tobytes())
    for item in boundary.species:
        digest.update(item.name.encode("utf-8"))
        digest.update(np.int64(item.charge_number).tobytes())
        digest.update(np.float64(item.mass_amu).tobytes())
        digest.update(np.float64(item.flux_m2_s).tobytes())
        digest.update(np.ascontiguousarray(item.velocity_sqrt_eV).tobytes())
        digest.update(np.ascontiguousarray(item.weight).tobytes())
    return digest.hexdigest()


def _ranked_fraction(flux: dict[str, float], limit: int = 8) -> list[dict]:
    total = float(sum(flux.values()))
    return [
        {
            "species": name,
            "flux_m2_s": float(value),
            "particle_flux_fraction": float(value / total),
        }
        for name, value in sorted(
            flux.items(), key=lambda item: (-item[1], item[0])
        )[:limit]
    ]


def build_receipt() -> dict:
    coupled = _load(COUPLED_STATE)
    axisymmetric = _load(AXISYMMETRIC_STATE)
    if (
        coupled["condition_id"] != CONDITION_ID
        or axisymmetric["condition_id"] != CONDITION_ID
        or coupled["input"]["feature_or_sem_target_used"]
        or axisymmetric["sem_target_used"]
        or axisymmetric["measured_depth_target_used"]
    ):
        raise ValueError("target outcome or mismatched condition entered boundary")

    species = {item.name: item for item in zhu_reactor_species()}
    state = coupled["state"]
    radial = axisymmetric["central_48x16_result"]
    total_local_ion_flux = float(radial[
        "central_3mm_optic_average_flux_m2_s"
    ])
    fractions = {
        str(name): float(value)
        for name, value in radial["species_flux_fraction"].items()
    }
    ion_flux = {
        name: total_local_ion_flux * fraction
        for name, fraction in fractions.items()
    }
    neutral_flux = {
        str(name): float(value)
        for name, value in state["neutral_thermal_flux_m2_s"].items()
    }
    if not np.isclose(sum(fractions.values()), 1.0, rtol=0.0, atol=2.0e-12):
        raise ValueError("axisymmetric ion fractions do not close")
    if set(ion_flux) != set(state["axial_positive_ion_flux_m2_s"]):
        raise ValueError("axisymmetric ion inventory does not match reactor")
    if set(neutral_flux) - set(species):
        raise ValueError("reactor neutral inventory lacks species metadata")

    ion_mass = {name: species[name].mass_amu for name in ion_flux}
    ion_charge = {name: species[name].charge_number for name in ion_flux}
    neutral_mass = {name: species[name].mass_amu for name in neutral_flux}
    electron_temperature_eV = (
        2.0 / 3.0 * float(state["mean_electron_energy_eV"])
    )
    sheath_drop_V = float(
        coupled["input"]["powered_electrode_sheath_drop_V"]
    )
    ion_energy = {
        name: charge * sheath_drop_V + 0.5 * electron_temperature_eV
        for name, charge in ion_charge.items()
    }

    cases = []
    for tail_fraction in (0.0, 0.65):
        boundary = build_species_resolved_feature_boundary(
            ion_flux_m2_s=ion_flux,
            ion_mass_amu=ion_mass,
            ion_charge_number=ion_charge,
            ion_energy_eV=ion_energy,
            ion_iadf=kim_2025_reference_iadf(
                tail_fraction=tail_fraction
            ),
            neutral_flux_m2_s=neutral_flux,
            neutral_mass_amu=neutral_mass,
            neutral_temperature_K=float(coupled["input"][
                "gas_temperature_K"
            ]),
            reference_plane_m=1.0e-6,
            ion_polar_order=12,
            ion_azimuthal_order=8,
            neutral_transverse_order=3,
            neutral_normal_order=4,
            provenance={
                "condition_id": CONDITION_ID,
                "tail_fraction_is_target_measurement": False,
                "target_sem_or_depth_used": False,
            },
        )
        recovered_ion_flux = sum(
            item.flux_m2_s
            for item in boundary.species if item.charge_number > 0
        )
        recovered_neutral_flux = sum(
            item.flux_m2_s
            for item in boundary.species if item.charge_number == 0
        )
        weighted_energy = sum(
            ion_flux[name] * boundary.get(name).mean_energy_eV
            for name in ion_flux
        ) / sum(ion_flux.values())
        cases.append({
            "name": f"kim_widths_tail_fraction_{tail_fraction:g}",
            "tail_fraction": tail_fraction,
            "tail_fraction_is_target_measurement": False,
            "ion_species_count": len(ion_flux),
            "neutral_species_count": len(neutral_flux),
            "quadrature_node_count": int(sum(
                item.weight.size for item in boundary.species
            )),
            "positive_ion_flux_m2_s": recovered_ion_flux,
            "neutral_thermal_flux_m2_s": recovered_neutral_flux,
            "ion_flux_relative_conservation_residual": (
                recovered_ion_flux / total_local_ion_flux - 1.0
            ),
            "neutral_flux_relative_conservation_residual": (
                recovered_neutral_flux / sum(neutral_flux.values()) - 1.0
            ),
            "current_density_A_m2": boundary.current_density_A_m2,
            "particle_flux_weighted_mean_ion_energy_eV": weighted_energy,
            "boundary_sha256": _boundary_hash(boundary),
            "deterministic_quadrature": True,
            "monte_carlo": False,
        })

    double_ions = {
        name: {
            "charge_number": ion_charge[name],
            "impact_energy_eV": ion_energy[name],
            "flux_m2_s": ion_flux[name],
        }
        for name in sorted(ion_flux) if ion_charge[name] > 1
    }
    return {
        "schema": "petch.zhu-npg80-species-resolved-feature-boundary.v1",
        "condition_id": CONDITION_ID,
        "inputs": {
            "coupled_reactor_state": str(COUPLED_STATE.relative_to(ROOT)),
            "coupled_reactor_state_sha256": _sha256(COUPLED_STATE),
            "axisymmetric_ion_state": str(
                AXISYMMETRIC_STATE.relative_to(ROOT)
            ),
            "axisymmetric_ion_state_sha256": _sha256(AXISYMMETRIC_STATE),
            "powered_electrode_sheath_drop_V": sheath_drop_V,
            "electron_temperature_eV": electron_temperature_eV,
            "gas_temperature_K": coupled["input"]["gas_temperature_K"],
            "impact_energy_rule": "Z_i*V_s + Te/2",
            "impact_energy_rule_is_target_iead_measurement": False,
        },
        "inventory": {
            "positive_ion_species_count": len(ion_flux),
            "thermal_neutral_species_count": len(neutral_flux),
            "dominant_positive_ions": _ranked_fraction(ion_flux),
            "dominant_thermal_neutrals": _ranked_fraction(neutral_flux),
            "multiply_charged_ions": double_ions,
        },
        "boundary_cases": cases,
        "certification": {
            "zero_fit_to_sem_or_depth": True,
            "reactor_species_identity_preserved": True,
            "ion_mass_and_charge_preserved": True,
            "absolute_flux_preserved": True,
            "deterministic_non_monte_carlo_boundary": True,
            "species_specific_iadf_interface_supported": True,
            "conditional_feature_transport_boundary_executable": True,
            "target_machine_self_bias_measured": False,
            "target_machine_iead_measured": False,
            "target_neutral_radial_transfer_solved": False,
            "target_tio2_cr_surface_probabilities_validated": False,
            "supports_unique_absolute_profile_prediction": False,
            "result": (
                "The reactor-to-feature software join is closed without "
                "species aggregation. Experimental boundary and surface-law "
                "gates remain open and are not converted into fitted claims."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    encoded = json.dumps(build_receipt(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != encoded:
            raise SystemExit("committed species-boundary audit is stale")
        print(encoded, end="")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded, encoding="utf-8")
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
