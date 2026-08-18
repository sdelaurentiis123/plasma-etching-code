#!/usr/bin/env python3
"""Build/check the target-free measured CHF2+ mobility scale receipt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from petch.reactor_global.chf3_ion_mobility import (
    load_basurto_2002_chf2_chf3_mobility_model,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "chf2_mobility_scale.json"
)
MANIFEST = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "recipe_manifest.json"
)
BOLTZMANN_J_K = 1.380649e-23
TORR_TO_PA = 133.322368
TABLE_TEMPERATURE_PROXY_K = 293.15


def load_pre_sem_manifest() -> dict:
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if (
        payload.get("schema") != "petch.experimental-recipe.v1"
        or payload.get("measurement_state") != "pre_sem_specific_condition"
        or payload["outcomes"]["specific_condition_sem_received"]
        or payload["outcomes"]["post_etch_tio2_depth_nm"] is not None
    ):
        raise ValueError("mobility scale requires the unrevealed pre-SEM case")
    return payload


def build_receipt() -> dict:
    manifest = load_pre_sem_manifest()
    pressure_Pa = manifest["process"]["pressure_Torr"] * TORR_TO_PA
    density = pressure_Pa / (BOLTZMANN_J_K * TABLE_TEMPERATURE_PROXY_K)
    model = load_basurto_2002_chf2_chf3_mobility_model()
    samples = []
    for reduced_field in (50.0, 100.0, 200.0, 400.0):
        state = model.evaluate(
            reduced_field_Td=reduced_field,
            total_neutral_density_m3=density,
        )
        samples.append({
            "reduced_field_Td": state.reduced_field_Td,
            "reduced_mobility_cm2_V_s": (
                state.reduced_mobility_cm2_V_s),
            "actual_mobility_m2_V_s": state.actual_mobility_m2_V_s,
            "drift_speed_m_s": state.drift_speed_m_s,
            "effective_momentum_relaxation_frequency_s_inv": (
                state.effective_momentum_relaxation_frequency_s_inv),
            "drift_relaxation_length_m": (
                state.drift_relaxation_length_m),
        })
    return {
        "schema": "petch.target-free-chf2-mobility-scale.v1",
        "condition_id": manifest["condition_id"],
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "coefficient_selected_from_target": None,
        "source": "basurto-2002-chf3-ion-mobility",
        "species_pair": "CHF2+ in CHF3",
        "density_proxy": {
            "pressure_Torr": manifest["process"]["pressure_Torr"],
            "temperature_K": TABLE_TEMPERATURE_PROXY_K,
            "temperature_basis": (
                "20 C electrode setpoint used only as an ideal-gas scale; "
                "the plasma gas temperature is not measured"
            ),
            "total_number_density_m3": density,
        },
        "measured_support": {
            "digitized_reduced_field_Td": list(
                model.reduced_field_support_Td),
            "source_relative_uncertainty": model.source_relative_uncertainty,
            "digitization_reduced_field_relative_bound": (
                model.digitization_reduced_field_relative_bound),
            "digitization_reduced_mobility_relative_bound": (
                model.digitization_reduced_mobility_relative_bound),
        },
        "samples": samples,
        "identifiability_gates": {
            "measured_chf2_chf3_swarm_mobility_available": True,
            "elastic_differential_cross_section_available": False,
            "angular_scattering_kernel_available": False,
            "target_ion_fraction_known": False,
            "supports_target_iead": False,
            "supports_absolute_depth_prediction": False,
        },
        "interpretation": (
            "The measured mobility closes a target-pressure momentum-"
            "relaxation scale for CHF2+ in CHF3. It is not a unique elastic "
            "cross section and cannot be substituted for the missing mixed-"
            "ion angular kernels or target species fractions."
        ),
    }


def _render(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = _render(build_receipt())
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed CHF2+ mobility receipt is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
