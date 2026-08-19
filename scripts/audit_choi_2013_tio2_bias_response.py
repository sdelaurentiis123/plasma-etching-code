#!/usr/bin/env python3
"""Freeze the source-reported TiO2 DC-bias response and its claim boundary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "data" / "experimental" / "choi_2013_tio2_cf4"
    / "bias_response_audit.json"
)


def build() -> dict[str, object]:
    voltage = np.asarray([50.0, 250.0])
    rate = np.asarray([130.9, 197.2])
    selectivity = np.asarray([0.65, 0.56])
    root_voltage = np.sqrt(voltage)
    root_slope = float((rate[1] - rate[0]) / np.diff(root_voltage)[0])
    root_intercept = float(rate[0] - root_slope * root_voltage[0])
    return {
        "schema": "petch.choi-2013-tio2-bias-response.v1",
        "source": {
            "citation": (
                "K.-R. Choi et al., Dry etching properties of TiO2 thin films "
                "in O2/CF4/Ar plasma, Vacuum 92, 85-89 (2013)"
            ),
            "doi": "10.1016/j.vacuum.2012.11.009",
            "figure": "Figure 5 and source-reported endpoint values",
            "full_text_read_online": True,
            "local_source_pdf_archived": False,
            "source_measurement_uncertainty_reported": False,
        },
        "condition": {
            "reactor": "planar inductively coupled plasma",
            "film": "approximately 200 nm electron-beam-evaporated TiO2",
            "feed_sccm": {"O2": 3.0, "CF4": 16.0, "Ar": 4.0},
            "source_power_W": 700.0,
            "pressure_Pa": 2.0,
        },
        "source_reported_endpoints": {
            "dc_bias_magnitude_V": voltage.tolist(),
            "tio2_etch_rate_nm_min": rate.tolist(),
            "tio2_to_sio2_selectivity": selectivity.tolist(),
        },
        "direct_response": {
            "bias_magnitude_ratio": float(voltage[1] / voltage[0]),
            "etch_rate_ratio": float(rate[1] / rate[0]),
            "selectivity_ratio": float(selectivity[1] / selectivity[0]),
            "etch_rate_monotone_in_bias_over_reported_interval": True,
            "selectivity_decreases_over_reported_interval": True,
        },
        "two_point_sqrt_bias_decomposition": {
            "equation": "rate_nm_min = intercept + slope*sqrt(abs(dc_bias_V))",
            "intercept_nm_min": root_intercept,
            "slope_nm_min_sqrtV": root_slope,
            "endpoint_replay_max_abs_error_nm_min": float(np.max(np.abs(
                root_intercept + root_slope * root_voltage - rate
            ))),
            "identified_as_surface_law": False,
            "reason": (
                "two endpoints cannot distinguish a thresholded sputter law, "
                "ion-assisted chemical desorption, or changing radical flux"
            ),
        },
        "surface_evidence": {
            "xps_ti_f_bond_formation_reported": True,
            "ion_bombardment_enhanced_bond_breaking_reported": True,
            "ion_bombardment_enhanced_product_desorption_reported": True,
            "oxygen_nonmonotonicity_reported": True,
            "minimum_model_topology": [
                "fluorinated TiOxFy or TiFx surface inventory",
                "neutral fluorine supply",
                "ion-energy-dependent bond breaking and product desorption",
                "competitive oxygen or fluorocarbon blocking/passivation state",
            ],
        },
        "model_discrimination": {
            "energy_independent_rate_normalized_removal_sufficient": False,
            "pure_physical_sputtering_only_supported": False,
            "ion_assisted_surface_chemistry_required": True,
        },
        "freddie_boundary": {
            "same_target_film_state": False,
            "same_feed_chemistry": False,
            "same_reactor_topology": False,
            "coefficient_transfer_allowed": False,
            "valid_use": (
                "response-sign and surface-state topology gate for a future "
                "TiO2 mechanism"
            ),
            "changes_absolute_oxford_depth_forecast": False,
        },
    }


def _render(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    rendered = _render(build())
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
        return
    if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("Choi TiO2 bias-response audit is stale")
    print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
