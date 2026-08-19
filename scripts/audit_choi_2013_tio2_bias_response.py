#!/usr/bin/env python3
"""Freeze Choi's source-reported TiO2 response axes and claim boundary.

The filename is retained because the original receipt was a bias-only board.
Schema v2 adds the source-reported oxygen, source-power, pressure, and AFM
endpoints.  These axes discriminate model topology; none identifies a
transferable Oxford coefficient because wafer fluxes and IEADs were not
reported.
"""
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
    oxygen_flow = np.asarray([0.0, 3.0, 9.0])
    oxygen_rate = np.asarray([154.1, 179.4, 137.5])
    source_power = np.asarray([600.0, 800.0])
    power_rate = np.asarray([136.0, 208.3])
    power_selectivity = np.asarray([0.66, 0.83])
    pressure = np.asarray([1.2, 2.8])
    pressure_rate = np.asarray([187.7, 138.7])
    pressure_selectivity = np.asarray([0.60, 0.50])
    root_voltage = np.sqrt(voltage)
    root_slope = float((rate[1] - rate[0]) / np.diff(root_voltage)[0])
    root_intercept = float(rate[0] - root_slope * root_voltage[0])
    return {
        "schema": "petch.choi-2013-tio2-multiaxis-response.v2",
        "source": {
            "citation": (
                "K.-R. Choi et al., Dry etching properties of TiO2 thin films "
                "in O2/CF4/Ar plasma, Vacuum 92, 85-89 (2013)"
            ),
            "doi": "10.1016/j.vacuum.2012.11.009",
            "figures": "Figures 2-6 and source-reported endpoint values",
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
            "oxygen_sweep": {
                "fixed_feed_sccm": {"CF4": 16.0, "Ar": 4.0},
                "source_power_W": 700.0,
                "dc_bias_magnitude_V": 150.0,
                "pressure_Pa": 2.0,
                "o2_flow_sccm": oxygen_flow.tolist(),
                "tio2_etch_rate_nm_min": oxygen_rate.tolist(),
                "tio2_to_sio2_selectivity_over_0_to_9_sccm": [0.57, 0.81],
            },
            "source_power_sweep": {
                "feed_sccm": {"O2": 3.0, "CF4": 16.0, "Ar": 4.0},
                "dc_bias_magnitude_V": 150.0,
                "pressure_Pa": 2.0,
                "source_power_W": source_power.tolist(),
                "tio2_etch_rate_nm_min": power_rate.tolist(),
                "tio2_to_sio2_selectivity": power_selectivity.tolist(),
            },
            "dc_bias_sweep": {
                "feed_sccm": {"O2": 3.0, "CF4": 16.0, "Ar": 4.0},
                "source_power_W": 700.0,
                "pressure_Pa": 2.0,
                "dc_bias_magnitude_V": voltage.tolist(),
                "tio2_etch_rate_nm_min": rate.tolist(),
                "tio2_to_sio2_selectivity": selectivity.tolist(),
            },
            "pressure_sweep": {
                "feed_sccm": {"O2": 3.0, "CF4": 16.0, "Ar": 4.0},
                "source_power_W": 700.0,
                "dc_bias_magnitude_V": 150.0,
                "pressure_Pa": pressure.tolist(),
                "tio2_etch_rate_nm_min": pressure_rate.tolist(),
                "tio2_to_sio2_selectivity": pressure_selectivity.tolist(),
            },
            "afm_rms_roughness_angstrom": {
                "as_deposited": 36.5,
                "cf4_ar": 59.8,
                "o2_cf4_ar": 29.8,
            },
        },
        "bias_direct_response": {
            "bias_magnitude_ratio": float(voltage[1] / voltage[0]),
            "etch_rate_ratio": float(rate[1] / rate[0]),
            "selectivity_ratio": float(selectivity[1] / selectivity[0]),
            "etch_rate_monotone_in_bias_over_reported_interval": True,
            "selectivity_decreases_over_reported_interval": True,
        },
        "multiaxis_direct_response": {
            "oxygen_rate_is_nonmonotonic": bool(
                oxygen_rate[1] > oxygen_rate[0]
                and oxygen_rate[1] > oxygen_rate[2]
            ),
            "oxygen_peak_flow_sccm": float(oxygen_flow[np.argmax(oxygen_rate)]),
            "oxygen_peak_rate_nm_min": float(np.max(oxygen_rate)),
            "source_power_rate_ratio": float(power_rate[1] / power_rate[0]),
            "source_power_selectivity_ratio": float(
                power_selectivity[1] / power_selectivity[0]
            ),
            "pressure_rate_ratio": float(pressure_rate[1] / pressure_rate[0]),
            "pressure_selectivity_ratio": float(
                pressure_selectivity[1] / pressure_selectivity[0]
            ),
            "oxygen_addition_reverses_cf4_ar_roughening": True,
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
            "cf4_ar_roughening_and_oxygen_assisted_smoothing_reported": True,
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
            "neutral_supply_and_competitive_oxygen_state_required": True,
            "collisional_sheath_or_ion_delivery_pressure_response_required": True,
            "chemistry_dependent_surface_morphology_response_required": True,
            "axis_interpretation": {
                "dc_bias": (
                    "closest source axis to an ion-energy perturbation, but the "
                    "unreported ion flux and IEAD prevent microscopic-yield inversion"
                ),
                "source_power": (
                    "entangles ion density, ion energy, electron kinetics, and neutral-F supply"
                ),
                "oxygen_flow": (
                    "requires competing F production/availability and oxygen blocking or cleanup"
                ),
                "pressure": (
                    "entangles neutral residence/supply with sheath collisions and ion delivery"
                ),
                "roughness": (
                    "requires a chemistry-dependent morphology observable, not only net depth; "
                    "the three AFM endpoints do not uniquely identify its microscopic cause"
                ),
            },
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
