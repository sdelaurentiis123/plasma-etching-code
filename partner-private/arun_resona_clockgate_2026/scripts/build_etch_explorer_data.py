#!/usr/bin/env python3
"""Build the receipt-backed Arun clock-gate process-explorer payload.

This is a deterministic, geometry-native entrance-transport calculation.  It
does not replace the common moving-interface engine and it deliberately does
not infer wafer fluxes from ICP/RF powers.  It answers the narrower question
that is already identifiable from the supplied mask: which incident
directions can reach the silicon through the 30 um printed entrance mask?

The resulting ion and direct-neutral transmission fields are then exposed to
the browser together with the provenance-bearing Belen SF6/O2 silicon law.
The browser may evaluate conditional surface-law transfers, normalized to the
independent Miao two-minute depth interval, but the payload refuses to label
those transfers as a target-tool prediction until the mask-wall return and
wafer boundary are measured.
"""
from __future__ import annotations

import argparse
import base64
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MAIN_REPO = ROOT.parents[1]
sys.path.insert(0, str(MAIN_REPO / "src"))

from petch.extruded_mask_transport import (  # noqa: E402
    cosine_flux_hemisphere_ordinates,
    direct_extruded_mask_transmission,
    gaussian_transverse_angle_ordinates,
)
from petch.mask_footprints import polygon_union_footprint_levelset  # noqa: E402


POLYGONS = ROOT / "derived" / "clkgate_x1_m1_polygons_um.json"
FEATURE_AUDIT = ROOT / "results" / "feature_geometry_audit.json"
RECIPE_BOARD = ROOT / "results" / "sf6_o2_recipe_board.json"
OUTPUT = ROOT / "results" / "etch_explorer_data.json"

MASK_HEIGHT_UM = 30.0
DISPLAY_DX_UM = 0.4
ION_SIGMAS_DEG = (0.5, 1.0, 1.5, 2.0, 3.0)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_geometry(dx_um: float):
    source = json.loads(POLYGONS.read_text(encoding="utf-8"))
    field = polygon_union_footprint_levelset(
        cell_width=source["cell_width_um"],
        cell_length=source["cell_length_um"],
        dx=dx_um,
        polygons=source["polygons"],
    )
    # The level-set helper returns duplicate periodic endpoints.  Transport
    # operates on the unique periodic cell so scipy's wrap mode is exact.
    gas = np.asarray(field[:-1, :-1] < 0.0, dtype=bool)
    return source, gas


def _line_of_sight_map(gas, *, dx_um, directions):
    return direct_extruded_mask_transmission(
        gas,
        mask_height=MASK_HEIGHT_UM,
        grid_spacing=dx_um,
        ordinates=directions,
        periodic_lateral=True,
        subdivisions_per_crossed_cell=2.0,
    ).transmission


def _ion_directions(sigma_deg: float, order: int):
    return gaussian_transverse_angle_ordinates(
        np.deg2rad(float(sigma_deg)), order_per_component=order)


def _cosine_neutral_directions(mu_order: int, azimuth_count: int):
    return cosine_flux_hemisphere_ordinates(
        polar_cosine_order=mu_order, azimuth_count=azimuth_count)


def _statistics(field, gas):
    values = np.asarray(field, dtype=float)[gas]
    quantile = np.quantile(values, [0.05, 0.5, 0.95])
    return {
        "mean": float(np.mean(values)),
        "minimum": float(np.min(values)),
        "p05": float(quantile[0]),
        "median": float(quantile[1]),
        "p95": float(quantile[2]),
        "maximum": float(np.max(values)),
    }


def _encode_fraction(field):
    quantized = np.rint(np.clip(field, 0.0, 1.0) * 255.0).astype(np.uint8)
    return base64.b64encode(np.ascontiguousarray(quantized).tobytes()).decode("ascii")


def _belen_rate_ratio(ion_transmission, neutral_direct, *, f_recovery, o_recovery,
                      ion_energy_eV=100.0):
    """Conditional Belen-law spatial rate normalized to its open field."""
    ion_flux = 12.0e19 * np.asarray(ion_transmission, dtype=float)
    f_flux = 1800.0e19 * (
        neutral_direct + float(f_recovery) * (1.0 - neutral_direct))
    o_flux = 100.0e19 * (
        neutral_direct + float(o_recovery) * (1.0 - neutral_direct))
    energy = float(ion_energy_eV)
    physical = ion_flux * 0.0337 * max(
        np.sqrt(energy) - np.sqrt(20.0), 0.0)
    ion_enhanced = ion_flux * 7.0 * max(
        np.sqrt(energy) - np.sqrt(15.0), 0.0)
    oxygen_sputter = ion_flux * 3.0 * max(
        np.sqrt(energy) - np.sqrt(10.0), 0.0)
    f_ratio = 0.5 * f_flux / (3.0e21 + 2.0 * ion_enhanced)
    o_ratio = o_flux / (4.0e17 + oxygen_sputter)
    theta_f = f_ratio / (1.0 + f_ratio + o_ratio)
    rate = 3.0e21 * theta_f / 4.0 + physical + theta_f * ion_enhanced

    open_ion = 12.0e19
    open_physical = open_ion * 0.0337 * max(
        np.sqrt(energy) - np.sqrt(20.0), 0.0)
    open_enhanced = open_ion * 7.0 * max(
        np.sqrt(energy) - np.sqrt(15.0), 0.0)
    open_o_sputter = open_ion * 3.0 * max(
        np.sqrt(energy) - np.sqrt(10.0), 0.0)
    open_f_ratio = 0.5 * 1800.0e19 / (3.0e21 + 2.0 * open_enhanced)
    open_o_ratio = 100.0e19 / (4.0e17 + open_o_sputter)
    open_theta_f = open_f_ratio / (1.0 + open_f_ratio + open_o_ratio)
    open_rate = (
        3.0e21 * open_theta_f / 4.0
        + open_physical
        + open_theta_f * open_enhanced
    )
    return rate / open_rate


def build():
    feature = json.loads(FEATURE_AUDIT.read_text(encoding="utf-8"))
    recipes = json.loads(RECIPE_BOARD.read_text(encoding="utf-8"))
    polygons, gas = _load_geometry(DISPLAY_DX_UM)

    ion_maps = {}
    ion_stats = {}
    for sigma in ION_SIGMAS_DEG:
        field = _line_of_sight_map(
            gas,
            dx_um=DISPLAY_DX_UM,
            directions=_ion_directions(sigma, 7),
        )
        key = f"{sigma:g}"
        ion_maps[key] = _encode_fraction(field)
        ion_stats[key] = _statistics(field, gas)

    neutral = _line_of_sight_map(
        gas,
        dx_um=DISPLAY_DX_UM,
        directions=_cosine_neutral_directions(8, 16),
    )
    neutral_stats = _statistics(neutral, gas)

    # Two independent numerical checks on the central ion-distribution map.
    ion_coarse_q5 = _line_of_sight_map(
        gas, dx_um=DISPLAY_DX_UM, directions=_ion_directions(1.5, 5))
    _, production_gas = _load_geometry(0.2)
    ion_production_q5 = _line_of_sight_map(
        production_gas, dx_um=0.2, directions=_ion_directions(1.5, 5))
    central_q7 = np.frombuffer(
        base64.b64decode(ion_maps["1.5"]), dtype=np.uint8
    ).reshape(gas.shape).astype(float) / 255.0

    central_rate = _belen_rate_ratio(
        central_q7, neutral, f_recovery=0.5, o_recovery=0.1)
    central_rate[~gas] = 0.0
    central_depth = 3.7 * central_rate

    nonpredictive = [
        "fluorine_sticking_probability",
        "oxygen_sticking_probability",
        "spontaneous_fluorine_removal_rate_m2_s",
        "oxygen_desorption_rate_m2_s",
        "physical_sputter_yield",
        "ion_enhanced_yield",
        "oxygen_sputter_yield",
        "target-tool wafer boundary",
        "printed-polymer F/O wall return",
    ]
    return {
        "schema": "petch.partner.clockgate-etch-explorer.v1",
        "claim_status": (
            "deterministic exact-topology entrance transport plus a conditional "
            "surface-law transfer; not a target-tool absolute-profile prediction"
        ),
        "target_sem_used": False,
        "inputs": {
            "polygon_receipt": {
                "path": str(POLYGONS.relative_to(ROOT)),
                "sha256": _sha(POLYGONS),
                "source_gds_sha256": polygons["source_gds_sha256"],
            },
            "feature_audit_sha256": _sha(FEATURE_AUDIT),
            "recipe_board_sha256": _sha(RECIPE_BOARD),
            "footprint_um": [
                polygons["cell_width_um"], polygons["cell_length_um"]],
            "mask_height_um": MASK_HEIGHT_UM,
            "minimum_opening_um": recipes["geometry_contract"]["minimum_opening_um"],
            "mask_area_fraction": feature["physical_geometry"][
                "exact_mask_area_fraction"],
            "working_polarity": feature["physical_geometry"]["working_polarity"],
        },
        "geometry": {
            "polygons_um": polygons["polygons"],
            "grid": {
                "dx_um": DISPLAY_DX_UM,
                "shape_xy": list(gas.shape),
                "fraction_encoding": "uint8_base64_row_major_x_then_y; value/255",
                "gas_opening": _encode_fraction(gas.astype(float)),
            },
        },
        "transport": {
            "operator": (
                "deterministic Gauss quadrature of complete straight characteristics "
                "through the exact periodic extruded mask"
            ),
            "ion": {
                "distribution": (
                    "independent zero-mean Gaussian transverse angular components"
                ),
                "component_sigma_deg": list(ION_SIGMAS_DEG),
                "gauss_hermite_order_per_component": 7,
                "transmission_maps": ion_maps,
                "statistics_over_opening": ion_stats,
            },
            "neutral_direct": {
                "distribution": "cosine incident flux hemisphere",
                "mu_gauss_legendre_order": 8,
                "azimuth_count": 16,
                "transmission_map": _encode_fraction(neutral),
                "statistics_over_opening": neutral_stats,
                "scope": (
                    "unreacted line-of-sight component only; diffuse wall return is "
                    "an exposed conditional input, not silently guessed"
                ),
            },
            "numerical_checks": {
                "ion_sigma_1p5_deg": {
                    "display_dx_0p4_q5_mean": _statistics(
                        ion_coarse_q5, gas)["mean"],
                    "display_dx_0p4_q7_mean": _statistics(
                        central_q7, gas)["mean"],
                    "production_dx_0p2_q5_mean": _statistics(
                        ion_production_q5, production_gas)["mean"],
                    "absolute_q5_to_q7_mean_change": abs(
                        _statistics(ion_coarse_q5, gas)["mean"]
                        - _statistics(central_q7, gas)["mean"]),
                    "absolute_dx_0p4_to_0p2_mean_change_q5": abs(
                        _statistics(ion_coarse_q5, gas)["mean"]
                        - _statistics(ion_production_q5, production_gas)["mean"]),
                },
            },
        },
        "surface_transfer": {
            "model": "BelenSiliconSF6O2Mechanism reference equations",
            "status": "literature-calibrated development law; not independently predictive",
            "open_boundary": {
                "ion_flux_m2_s": 12.0e19,
                "f_flux_m2_s": 1800.0e19,
                "o_flux_m2_s": 100.0e19,
                "mean_ion_energy_eV": 100.0,
            },
            "parameters": {
                "fluorine_sticking_probability": 0.5,
                "oxygen_sticking_probability": 1.0,
                "spontaneous_fluorine_removal_rate_m2_s": 3.0e21,
                "oxygen_desorption_rate_m2_s": 4.0e17,
                "physical_sputter_prefactor_per_sqrt_eV": 0.0337,
                "ion_enhanced_prefactor_per_sqrt_eV": 7.0,
                "oxygen_sputter_prefactor_per_sqrt_eV": 3.0,
            },
            "normalization": {
                "source": "miao-2016-cryo-grating",
                "condition": recipes["selected_baseline"],
                "two_minute_open_feature_depth_um": [3.5, 3.9],
                "use": (
                    "independent thin-mask process anchor; not Arun target fitting"
                ),
            },
            "nonpredictive_inputs": nonpredictive,
            "central_scenario_for_ui_only": {
                "ion_component_sigma_deg": 1.5,
                "f_wall_recovery_fraction": 0.5,
                "o_wall_recovery_fraction": 0.1,
                "open_feature_depth_um": 3.7,
                "depth_statistics_over_opening_um": _statistics(
                    central_depth, gas),
                "warning": (
                    "the recovery fractions are an interactive scenario, not measured values"
                ),
            },
        },
        "recipe_board": recipes,
        "refusal": {
            "unique_absolute_profile_ready": False,
            "missing_measurements": [
                "exact etcher model and achieved DC self-bias or ion-energy distribution",
                "wafer F, O, and positive-ion fluxes or a same-run blanket Si depth",
                "printed-polymer F/O reaction or return probability and mask loss",
                "pre/post cross-section for a blind held-out profile score",
            ],
            "next_physics_step": (
                "solve diffuse reactive wall return and moving sidewalls through the common "
                "3-D engine over the declared wafer-boundary and polymer-response ensemble"
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))
    elif not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit("etch-explorer data are stale")
    else:
        print(f"PASS {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
