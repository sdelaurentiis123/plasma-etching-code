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
        "identifiability_gates": {
            "recipe_and_stack_frozen": True,
            "specific_condition_sem_withheld": True,
            "achieved_dc_self_bias_measured": (
                process["measured_dc_self_bias_V"] is not None
            ),
            "absorbed_plasma_power_measured": False,
            "species_resolved_wafer_fluxes_measured_or_validated": False,
            "tio2_surface_law_measured_or_validated_for_condition": False,
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
