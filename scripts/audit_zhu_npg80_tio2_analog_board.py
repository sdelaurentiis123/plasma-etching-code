#!/usr/bin/env python3
"""Build/check the target-free Janissen TiO2 process-analogy board."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[1]
TARGET_DIR = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
SOURCE_DIR = ROOT / "data" / "experimental" / "janissen_2016_tio2"
TARGET_MANIFEST = TARGET_DIR / "recipe_manifest.json"
TABLE_S31 = SOURCE_DIR / "table_s3_1_chf3_rie_optimization.csv"
TABLE_S32 = SOURCE_DIR / "table_s3_2_tio2_nanofabrication.csv"
TABLE_S33 = SOURCE_DIR / "table_s3_3_feature_dimensions.csv"
OUTPUT = TARGET_DIR / "janissen_tio2_analog_board.json"
MILLITORR_PER_MICROBAR = 0.750061683


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_receipt() -> dict:
    manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    if (
        manifest.get("measurement_state") != "pre_sem_specific_condition"
        or manifest["outcomes"]["specific_condition_sem_received"]
        or manifest["outcomes"]["post_etch_tio2_depth_nm"] is not None
    ):
        raise ValueError("analogy board requires the unrevealed pre-SEM case")
    process = manifest["process"]
    stack = manifest["stack"]
    target_drive = (
        process["table_rf_forward_power_setpoint_W"]
        / (1000.0 * process["pressure_Torr"])
    )

    s31 = _rows(TABLE_S31)
    ra = sorted(
        (row for row in s31 if row["batch"] == "Ra"
         and float(row["pressure_ubar"]) == 50.0),
        key=lambda row: float(row["rf_power_W"]),
    )
    power = np.asarray([float(row["rf_power_W"]) for row in ra])
    tio2_rate = np.asarray([
        float(row["tio2_etch_rate_nm_min"]) for row in ra])
    cr_rate = np.asarray([float(row["cr_etch_rate_nm_min"]) for row in ra])
    target_power = float(process["table_rf_forward_power_setpoint_W"])
    if not power[0] <= target_power <= power[-1]:
        raise ValueError("target power is outside the source power sweep")
    source_tio2_at_target_power = float(PchipInterpolator(
        power, tio2_rate, extrapolate=False)(target_power))
    source_cr_at_target_power = float(PchipInterpolator(
        power, cr_rate, extrapolate=False)(target_power))
    source_selectivity = source_tio2_at_target_power / source_cr_at_target_power
    film_nm = float(stack["film_initial_thickness_nm"])
    mask_nm = float(stack["mask_initial_thickness_nm"])

    s32 = {row["figure"]: row for row in _rows(TABLE_S32)}
    closest = s32["3.2a"]
    source_pressure_mTorr = (
        float(closest["pressure_ubar"]) * MILLITORR_PER_MICROBAR)
    source_drive = float(closest["rf_power_W"]) / source_pressure_mTorr
    s33 = {row["figure"]: row for row in _rows(TABLE_S33)}
    feature_rates = {
        figure: float(row["average_height_nm"]) / float(row["etch_time_min"])
        for figure, row in s33.items()
    }

    clear_time = film_nm / source_tio2_at_target_power
    cr_for_clear = film_nm / source_selectivity
    return {
        "schema": "petch.target-free-tio2-analog-board.v1",
        "condition_id": manifest["condition_id"],
        "sem_target_used": False,
        "measured_depth_target_used": False,
        "coefficient_selected_from_target": None,
        "source": "janissen-2016-tio2-rie",
        "closest_stack_witness": {
            "figure": "3.2a",
            "system_model": closest["system_model"],
            "source_material": "single-crystal rutile TiO2",
            "target_material": "ALD TiO2; phase and density unreported",
            "mask_height_nm": float(closest["mask_height_nm"]),
            "mask_diameter_nm": float(closest["mask_diameter_nm"]),
            "gas_flows_sccm": {
                "CHF3": float(closest["CHF3_sccm"]),
                "O2": float(closest["O2_sccm"]),
            },
            "rf_power_W": float(closest["rf_power_W"]),
            "pressure_mTorr": source_pressure_mTorr,
            "forward_power_per_pressure_W_mTorr": source_drive,
            "dc_bias_V_signed": float(closest["dc_bias_V_signed"]),
            "etch_time_min": float(closest["etch_time_min"]),
            "reported_feature_height_nm": 430.0,
            "implied_feature_rate_nm_min": 430.0 / 11.0,
            "reported_approximate_selectivity": 14.0,
        },
        "target_similarity": {
            "target_forward_power_per_pressure_W_mTorr": target_drive,
            "source_to_target_reduced_drive_ratio": source_drive / target_drive,
            "same_cr_mask_thickness": (
                float(closest["mask_height_nm"]) == mask_nm),
            "same_machine": False,
            "same_tio2_material_state": False,
            "same_active_gas_set": False,
            "target_added_SF6_sccm": process["gases_sccm"]["SF6"],
        },
        "source_power_sweep_interpolation": {
            "scope": (
                "PCHIP interpolation inside the Janissen Ra batch only; "
                "not transferred to the target"
            ),
            "power_support_W": power.tolist(),
            "target_power_W": target_power,
            "source_system_tio2_rate_nm_min": source_tio2_at_target_power,
            "source_system_cr_rate_nm_min": source_cr_at_target_power,
            "source_system_selectivity": source_selectivity,
            "source_system_clear_time_for_700nm_min": clear_time,
            "source_system_cr_consumption_for_700nm_nm": cr_for_clear,
            "source_system_cr_residual_after_700nm_nm": mask_nm - cr_for_clear,
            "supports_target_prediction": False,
        },
        "source_feature_depth_board": {
            "same_F2_nominal_recipe_different_mask_and_pattern": [
                {
                    "figure": figure,
                    "etch_time_min": float(s33[figure]["etch_time_min"]),
                    "height_nm": float(s33[figure]["average_height_nm"]),
                    "implied_rate_nm_min": feature_rates[figure],
                    "height_global_rsd_percent": float(
                        s33[figure]["height_global_rsd_percent"]),
                }
                for figure in ("S3.4", "3.3")
            ],
            "minimum_implied_rate_nm_min": min(feature_rates.values()),
            "maximum_implied_rate_nm_min": max(feature_rates.values()),
            "target_required_clear_rate_nm_min": (
                film_nm / (process["etch_time_s"] / 60.0)),
            "warning": (
                "The non-proportional source heights differ in mask and "
                "pattern as well as time; they are validation outcomes, not "
                "a temporal-rate fit."
            ),
        },
        "identifiability_gates": {
            "exact_source_tables_visually_audited": True,
            "source_power_response_available": True,
            "source_feature_depths_available": True,
            "target_machine_boundary_validated": False,
            "target_material_surface_law_validated": False,
            "supports_absolute_target_depth_prediction": False,
        },
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
            raise SystemExit("committed TiO2 analog board is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
