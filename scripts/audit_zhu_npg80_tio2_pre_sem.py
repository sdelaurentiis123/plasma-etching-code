#!/usr/bin/env python3
"""Build/check the target-free Zhu NPG80 TiO2/Cr pre-SEM receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
DEFAULT_MANIFEST = DATA_DIR / "recipe_manifest.json"
DEFAULT_OUTPUT = DATA_DIR / "pre_sem_receipt.json"
SELF_BIAS_TRANSFER = DATA_DIR / "oxford80_self_bias_transfer.json"
CF3_COLLISION_SCALE = DATA_DIR / "cf3_reactive_collision_scale.json"
CHF2_MOBILITY_SCALE = DATA_DIR / "chf2_mobility_scale.json"
TIO2_ANALOG_BOARD = DATA_DIR / "janissen_tio2_analog_board.json"
TIO2_DEPTH_GATE = DATA_DIR / "pre_sem_depth_gate.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != "petch.experimental-recipe.v1":
        raise ValueError("unexpected recipe-manifest schema")
    if payload.get("measurement_state") != "pre_sem_specific_condition":
        raise ValueError("condition is not marked as a pre-SEM commitment")
    if payload["outcomes"]["specific_condition_sem_received"]:
        raise ValueError("pre-SEM receipt cannot consume a revealed SEM")
    if payload["outcomes"]["post_etch_tio2_depth_nm"] is not None:
        raise ValueError("pre-SEM receipt cannot consume a measured depth")
    image = ROOT / payload["source_image"]["path"]
    if _sha256(image) != payload["source_image"]["sha256"]:
        raise ValueError("recipe screenshot checksum mismatch")
    return payload


def build_receipt(manifest: dict) -> dict:
    process = manifest["process"]
    stack = manifest["stack"]
    gases = process["gases_sccm"]

    total_flow = sum(gases.values())
    if total_flow <= 0.0:
        raise ValueError("total process-gas flow must be positive")
    gas_fractions = {
        species: flow / total_flow for species, flow in gases.items()
    }

    etch_time_min = process["etch_time_s"] / 60.0
    film_nm = stack["film_initial_thickness_nm"]
    mask_nm = stack["mask_initial_thickness_nm"]
    required_rate = film_nm / etch_time_min
    required_selectivity = film_nm / mask_nm

    # Literature comparison only: Janissen et al. report approximately
    # 40 nm/min and 14:1 for an adjacent CHF3/O2 TiO2/Cr RIE process. These
    # values are never installed as coefficients in the reactor or surface
    # model and are explicitly not the supplied SF6-containing condition.
    adjacent_rate = 40.0
    adjacent_selectivity = 14.0
    adjacent_supported_relief = mask_nm * adjacent_selectivity
    adjacent_mask_required = film_nm / adjacent_selectivity
    adjacent_clear_time = film_nm / adjacent_rate
    adjacent_mask_exhaustion_time = adjacent_supported_relief / adjacent_rate

    self_bias_transfer = None
    if SELF_BIAS_TRANSFER.exists():
        self_bias_transfer = json.loads(
            SELF_BIAS_TRANSFER.read_text(encoding="utf-8"))
        if self_bias_transfer.get("condition_id") != manifest["condition_id"]:
            raise ValueError("self-bias transfer belongs to another condition")

    cf3_collision_scale = None
    if CF3_COLLISION_SCALE.exists():
        cf3_collision_scale = json.loads(
            CF3_COLLISION_SCALE.read_text(encoding="utf-8"))
        if cf3_collision_scale.get("condition_id") != manifest["condition_id"]:
            raise ValueError("CF3 collision receipt belongs to another condition")

    chf2_mobility_scale = None
    if CHF2_MOBILITY_SCALE.exists():
        chf2_mobility_scale = json.loads(
            CHF2_MOBILITY_SCALE.read_text(encoding="utf-8"))
        if chf2_mobility_scale.get("condition_id") != manifest["condition_id"]:
            raise ValueError("CHF2 mobility receipt belongs to another condition")

    tio2_analog_board = None
    if TIO2_ANALOG_BOARD.exists():
        tio2_analog_board = json.loads(
            TIO2_ANALOG_BOARD.read_text(encoding="utf-8"))
        if tio2_analog_board.get("condition_id") != manifest["condition_id"]:
            raise ValueError("TiO2 analog board belongs to another condition")

    tio2_depth_gate = None
    if TIO2_DEPTH_GATE.exists():
        tio2_depth_gate = json.loads(
            TIO2_DEPTH_GATE.read_text(encoding="utf-8"))
        if tio2_depth_gate.get("condition_id") != manifest["condition_id"]:
            raise ValueError("TiO2 depth gate belongs to another condition")
        if (
            tio2_depth_gate.get("sem_target_used") is not False
            or tio2_depth_gate.get("measured_depth_target_used") is not False
        ):
            raise ValueError("TiO2 depth gate is not target-free")

    return {
        "schema": "petch.pre-sem-depth-receipt.v1",
        "condition_id": manifest["condition_id"],
        "frozen_on": manifest["received_on"],
        "claim_class": (
            "target-free recipe arithmetic plus adjacent-experiment "
            "comparison; not an absolute-depth prediction"
        ),
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "coefficient_selected_from_target": None,
        "authoritative_recipe": {
            "etch_time_min": etch_time_min,
            "pressure_Torr": process["pressure_Torr"],
            "table_rf_forward_power_setpoint_W": (
                process["table_rf_forward_power_setpoint_W"]
            ),
            "table_temperature_C": process["table_temperature_C"],
            "total_process_flow_sccm": total_flow,
            "gas_flow_fractions": gas_fractions,
            "film_initial_thickness_nm": film_nm,
            "film_deposition_method": stack["film_deposition_method"],
            "mask_initial_thickness_nm": mask_nm,
            "substrate_material": stack["substrate_material"],
            "device_family": manifest["experimental_program"]["device_family"],
        },
        "necessary_conditions_for_full_clear": {
            "minimum_effective_tio2_rate_nm_min": required_rate,
            "minimum_zero_margin_tio2_to_cr_selectivity": required_selectivity,
            "warning": (
                "Both gates are necessary. Practical mask survival needs "
                "selectivity above the zero-margin threshold."
            ),
        },
        "adjacent_literature_comparison": {
            "source": "janissen-2016-tio2-rie",
            "chemistry_boundary": (
                "CHF3/O2 on TiO2 with Cr mask; not the supplied "
                "CHF3/SF6/O2 condition and not machine-transferable"
            ),
            "approximate_tio2_rate_nm_min": adjacent_rate,
            "approximate_tio2_to_cr_selectivity": adjacent_selectivity,
            "rate_implied_clear_time_min": adjacent_clear_time,
            "mask_required_for_700nm_clear_nm": adjacent_mask_required,
            "mask_supported_tio2_relief_nm": adjacent_supported_relief,
            "mask_exhaustion_time_at_comparison_rate_min": (
                adjacent_mask_exhaustion_time
            ),
            "mask_shortfall_for_full_clear_nm": (
                adjacent_mask_required - mask_nm
            ),
            "time_after_comparison_mask_exhaustion_min": (
                etch_time_min - adjacent_mask_exhaustion_time
            ),
            "interpretation": (
                "The comparison rate clears the film in time, but the "
                "comparison selectivity exhausts 45 nm Cr before full clear."
            ),
        },
        "adjacent_same_group_device_evidence": {
            **manifest["experimental_program"]["adjacent_publication"],
            "used_as_condition_matched_target": False,
        },
        "machine_family_self_bias_evidence": {
            "available": self_bias_transfer is not None,
            "receipt": (
                str(SELF_BIAS_TRANSFER.relative_to(ROOT))
                if self_bias_transfer is not None else None
            ),
            "matched_chemistry_reduced_drive_anchor_V": (
                self_bias_transfer["mechanical_anchor_selection"]["anchor_V"]
                if self_bias_transfer is not None else None
            ),
            "exact_tool_conditioning_drift_is_censored": (
                True if self_bias_transfer is not None else None
            ),
            "interpretation": (
                "Enables a deterministic machine-family sensitivity ensemble, "
                "not a unique target-condition voltage or depth prediction."
                if self_bias_transfer is not None else
                "No machine-family voltage transfer is committed."
            ),
        },
        "measured_molecular_collision_evidence": {
            "available": cf3_collision_scale is not None,
            "receipt": (
                str(CF3_COLLISION_SCALE.relative_to(ROOT))
                if cf3_collision_scale is not None else None
            ),
            "covered_pair_and_channels": (
                "CF3+ + CHF3 summed CID and summed DCT destruction"
                if cf3_collision_scale is not None else None
            ),
            "complete_molecular_transport": False,
            "interpretation": (
                "Measured reactive destruction proves collisions are "
                "non-negligible; elastic/angular kernels and the remaining "
                "ion-neutral pairs are still unresolved."
                if cf3_collision_scale is not None else
                "No target-relevant molecular collision kernel committed."
            ),
        },
        "measured_molecular_mobility_evidence": {
            "available": chf2_mobility_scale is not None,
            "receipt": (
                str(CHF2_MOBILITY_SCALE.relative_to(ROOT))
                if chf2_mobility_scale is not None else None
            ),
            "covered_pair_and_observable": (
                "CHF2+ in CHF3 reduced mobility"
                if chf2_mobility_scale is not None else None
            ),
            "elastic_differential_cross_section": False,
            "interpretation": (
                "Measured mass-resolved mobility closes a bulk/presheath "
                "momentum-relaxation scale, not an angular collision kernel "
                "or a target sheath IEAD."
                if chf2_mobility_scale is not None else
                "No target-relevant molecular-ion mobility is committed."
            ),
        },
        "audited_tio2_process_analogs": {
            "available": tio2_analog_board is not None,
            "receipt": (
                str(TIO2_ANALOG_BOARD.relative_to(ROOT))
                if tio2_analog_board is not None else None
            ),
            "closest_source_dc_bias_V_signed": (
                tio2_analog_board["closest_stack_witness"]
                ["dc_bias_V_signed"]
                if tio2_analog_board is not None else None
            ),
            "source_feature_depth_board_available": (
                True if tio2_analog_board is not None else False
            ),
            "transferred_as_target_coefficient": False,
        },
        "blind_tio2_clearance_forecast": {
            "available": tio2_depth_gate is not None,
            "receipt": (
                str(TIO2_DEPTH_GATE.relative_to(ROOT))
                if tio2_depth_gate is not None else None
            ),
            "forecast_type": (
                tio2_depth_gate["blind_forecast"]["forecast_type"]
                if tio2_depth_gate is not None else None
            ),
            "primary_outcome": (
                tio2_depth_gate["blind_forecast"]["primary_outcome"]
                if tio2_depth_gate is not None else None
            ),
            "predicted_film_capped_tio2_depth_nm": (
                tio2_depth_gate["blind_forecast"]["predicted_tio2_depth_nm"]
                if tio2_depth_gate is not None else None
            ),
            "promoted_to_absolute_feature_profile_prediction": False,
            "interpretation": (
                "A preregistered binary film-clearance call is available, "
                "while the continuous feature profile remains unidentified."
                if tio2_depth_gate is not None else
                "No target-free clearance call is committed."
            ),
        },
        "identifiability_gates": {
            "recipe_and_stack_frozen": True,
            "specific_condition_sem_withheld": True,
            "achieved_dc_self_bias_measured": (
                process["measured_dc_self_bias_V"] is not None
            ),
            "machine_family_self_bias_sensitivity_available": (
                self_bias_transfer is not None
            ),
            "absorbed_plasma_power_measured": False,
            "measured_cf3_chf3_reactive_collision_kernel_available": (
                cf3_collision_scale is not None
            ),
            "measured_chf2_chf3_swarm_mobility_available": (
                chf2_mobility_scale is not None
            ),
            "complete_molecular_ion_transport_available": False,
            "species_resolved_wafer_fluxes_measured_or_validated": False,
            "tio2_surface_law_measured_or_validated_for_condition": False,
            "adjacent_tio2_process_response_board_available": (
                tio2_analog_board is not None
            ),
            "target_free_binary_clearance_forecast_frozen": (
                tio2_depth_gate is not None
            ),
            "cr_surface_law_measured_or_validated_for_condition": False,
            "feature_geometry_and_loading_known": False,
            "supports_absolute_depth_prediction": False,
        },
        "next_discriminating_observations": [
            "achieved DC self-bias and, if available, RF voltage/current",
            "blanket TiO2 removal and Cr loss or residual Cr",
            "TiO2 phase/density and actual wafer-surface temperature",
            "mask geometry, pitch, loading, chip/carrier size and position",
            "scale-bearing SEM for the exact condition",
        ],
    }


def _canonical_json(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed receipt differs from a clean rebuild",
    )
    args = parser.parse_args()

    receipt = build_receipt(load_manifest(args.manifest))
    rendered = _canonical_json(receipt)
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing committed receipt: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed pre-SEM receipt is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
