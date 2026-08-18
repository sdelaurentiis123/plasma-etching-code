#!/usr/bin/env python3
"""Build/check the target-free measured CF3+ reactive-collision scale."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from petch.reactor_global.chf3_ion_collisions import (
    ANGSTROM2_TO_M2,
    load_peko_2002_cf3_chf3_reactive_collision_model,
)
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "recipe_manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "cf3_reactive_collision_scale.json"
)
BOLTZMANN_J_K = 1.380649e-23
PASCAL_PER_MTORR = 0.133322368


def load_pre_sem_manifest(path: Path = DEFAULT_MANIFEST) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if (
        payload.get("schema") != "petch.experimental-recipe.v1"
        or payload.get("measurement_state") != "pre_sem_specific_condition"
        or payload["outcomes"]["specific_condition_sem_received"]
        or payload["outcomes"]["post_etch_tio2_depth_nm"] is not None
    ):
        raise ValueError("collision scale requires the unrevealed pre-SEM case")
    return payload


def build_receipt(*, gas_temperature_K: float = 350.0) -> dict:
    manifest = load_pre_sem_manifest(DEFAULT_MANIFEST)
    process = manifest["process"]
    model = load_peko_2002_cf3_chf3_reactive_collision_model()
    flows = process["gases_sccm"]
    total_flow = float(sum(flows.values()))
    chf3_feed_fraction = float(flows["CHF3"] / total_flow)
    pressure_mTorr = float(process["pressure_Torr"] * 1000.0)
    total_density = (
        pressure_mTorr * PASCAL_PER_MTORR
        / (BOLTZMANN_J_K * float(gas_temperature_K))
    )
    chf3_density_proxy = total_density * chf3_feed_fraction
    path_length_m = 1.0e-3
    samples = []
    for laboratory_energy_eV in (50.0, 100.0, 200.0, 300.0, 380.0):
        item = model.slab_sensitivity(
            laboratory_energy_eV=laboratory_energy_eV,
            chf3_number_density_m3=chf3_density_proxy,
            path_length_m=path_length_m,
            feed_fraction_used_as_density_proxy=True,
        )
        samples.append({
            "laboratory_energy_eV": item.laboratory_energy_eV,
            "relative_collision_energy_eV": (
                item.relative_collision_energy_eV),
            "reactive_cross_section_A2": {
                "lower": item.cross_section.lower_m2 / ANGSTROM2_TO_M2,
                "central": item.cross_section.central_m2 / ANGSTROM2_TO_M2,
                "upper": item.cross_section.upper_m2 / ANGSTROM2_TO_M2,
            },
            "reactive_optical_depth_per_1mm": {
                "lower": item.optical_depth_lower,
                "central": item.optical_depth_central,
                "upper": item.optical_depth_upper,
            },
            "reactive_destruction_probability_per_1mm": {
                "lower": item.destruction_probability_lower,
                "central": item.destruction_probability_central,
                "upper": item.destruction_probability_upper,
            },
        })
    return {
        "schema": "petch.cf3-reactive-collision-scale.v1",
        "condition_id": manifest["condition_id"],
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "source": "peko-2002-chf3-ion-molecule",
        "model": "Peko2002CF3CHF3ReactiveCollisionModel",
        "energy_convention": {
            "source_axis": "relative collision energy",
            "projectile": "CF3+",
            "stationary_target": "CHF3",
            "projectile_mass_amu": model.projectile_mass_amu,
            "target_mass_amu": model.target_mass_amu,
            "laboratory_to_relative_factor": (
                model.laboratory_to_relative_energy_factor),
            "supported_laboratory_energy_eV": list(
                model.laboratory_energy_support_eV),
        },
        "declared_scale_assumptions": {
            "pressure_mTorr": pressure_mTorr,
            "gas_temperature_K": gas_temperature_K,
            "ideal_total_neutral_density_m3": total_density,
            "chf3_feed_fraction": chf3_feed_fraction,
            "chf3_density_proxy_m3": chf3_density_proxy,
            "feed_fraction_used_as_sheath_edge_density_proxy": True,
            "comparison_path_length_m": path_length_m,
            "path_is_not_an_inferred_sheath_thickness": True,
        },
        "samples": samples,
        "interpretation": (
            "The measured CF3+ destruction scale is non-negligible at the "
            "target pressure. This rejects a collisionless final boundary, "
            "but the feed-fraction density proxy and fixed 1 mm slab are "
            "sensitivity assumptions rather than a target-tool transport "
            "solution."
        ),
        "identifiability_gates": {
            "measured_cf3_chf3_reactive_kernel_available": True,
            "elastic_or_momentum_transfer_kernel_available": False,
            "angular_scattering_kernel_available": False,
            "all_target_ion_neutral_pairs_available": False,
            "target_sheath_neutral_composition_validated": False,
            "supports_complete_molecular_transport": False,
            "supports_target_iead": False,
            "supports_absolute_depth_prediction": False,
        },
    }


def _canonical_json(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--gas-temperature-K", type=float, default=350.0)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _canonical_json(build_receipt(
        gas_temperature_K=args.gas_temperature_K))
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing committed receipt: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed CF3 reactive-collision receipt is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
