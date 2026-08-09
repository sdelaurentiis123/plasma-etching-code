#!/usr/bin/env python3
"""Propagate the 2-D chlorine ion tier through sheath and clean surface.

The existing Mahorowala board is normalized on the thesis's independent
400 W/100 sccm center-current estimate.  This audit replaces the global
Lee--Lieberman edge partition with each declared axisymmetric source moment,
renormalizes *once* at that same center-current anchor, and recomputes the
species-resolved RF sheath and clean surface response.  No depth selects a
source moment or normalization.

The emitted reactor receipts can be passed directly to the deterministic
feature script with ``--product-wall none``.  Product-return fields are left in
the copied rows for provenance only and are explicitly invalidated because
their nonlinear fixed point has not been recomputed with the spatial boundary.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from forecast_mahorowala_1998_cl2_depth import (
    ELECTRODE_AREA_M2,
    PLASMA_POTENTIAL_EV,
    _combined_surface_projection,
)
from petch.chlorine_species_resolved_si import (
    SpeciesResolvedChlorineSiMechanism,
)
from petch.reactor_global import DiagnosticConditionedRFSheathTransfer


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REACTOR = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "mahorowala_1998_diagnostic_conditioned_depth_projection.json"
)
DEFAULT_SPATIAL = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "axisymmetric_ion_transport_audit"
    / "mahorowala_1998_axisymmetric_ion_transport.json"
)
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "axisymmetric_depth_reachability"
)
POSITIVE_IONS = ("Cl+", "Cl2+")


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _receipt_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def audit(reactor_path: Path, spatial_path: Path):
    base = json.loads(reactor_path.read_text(encoding="utf-8"))
    spatial = json.loads(spatial_path.read_text(encoding="utf-8"))
    spatial_rows = {
        (float(row["inductive_power_W"]), float(row["cl2_flow_sccm"])): row
        for row in spatial["rows"]
    }
    center_key = (400.0, 100.0)
    center = spatial_rows[center_key]
    source_modes = tuple(center["spatial_source_modes"])
    reference_total_flux = float(
        base["conditioning"]["reference_wafer_total_ion_flux_m2_s"])
    mechanism = SpeciesResolvedChlorineSiMechanism()
    sheath_transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=ELECTRODE_AREA_M2,
        plasma_potential_eV=PLASMA_POTENTIAL_EV,
        frequency_hz=13.56e6,
        collapse_fraction=1.0,
        phase_count=96,
        steps_per_period=128,
        steps_per_transit=128,
        source=(
            "Mahorowala axisymmetric charged transport propagated through "
            "the deterministic RF sheath"
        ),
    )
    receipts = {}
    summaries = {}
    for mode_name in source_modes:
        center_raw = center["spatial_source_modes"][mode_name][
            "wafer_flux_m2_s"]
        center_raw_total = float(sum(
            float(center_raw[name]) for name in POSITIVE_IONS))
        normalization = reference_total_flux / center_raw_total
        receipt = deepcopy(base)
        receipt["schema"] = (
            "petch.mahorowala_1998_axisymmetric_diagnostic_depth_projection.v1"
        )
        receipt["claim_class"] = (
            "diagnostic-center-conditioned axisymmetric reactor-to-surface "
            "projection; source-field sensitivity, not formal feature depth"
        )
        receipt["axisymmetric_transport"] = {
            "source_moment": mode_name,
            "source_moment_provenance": center[
                "spatial_source_modes"][mode_name]["source_provenance"],
            "source_moment_conditioned_or_measured": center[
                "spatial_source_modes"][mode_name][
                    "source_conditioned_or_measured"],
            "single_center_current_normalization": normalization,
            "center_condition": "400 W, 100 sccm",
            "center_reference_total_positive_ion_flux_m2_s": (
                reference_total_flux),
            "feature_depth_used": False,
            "supports_implicit_differentiation": True,
            "supports_absolute_wafer_flux_prediction": False,
            "product_return_fixed_point_recomputed": False,
        }
        errors = []
        positive_residual = 0
        negative_residual = 0
        for row in receipt["rows"]:
            key = (
                float(row["inductive_power_W"]),
                float(row["cl2_flow_sccm"]),
            )
            spatial_mode = spatial_rows[key]["spatial_source_modes"][mode_name]
            positive_flux = {
                name: normalization * float(
                    spatial_mode["wafer_flux_m2_s"][name])
                for name in POSITIVE_IONS
            }
            total_flux = float(sum(positive_flux.values()))
            boundary = SimpleNamespace(
                atomic_chlorine_flux_m2_s=float(
                    row["wafer_atomic_chlorine_flux_m2_s"]),
                positive_ion_flux_m2_s=positive_flux,
            )
            sheath = sheath_transfer.predict(
                positive_ion_flux_m2_s=positive_flux,
                electron_temperature_eV=(
                    (2.0 / 3.0)
                    * float(row["reactor_mean_electron_energy_eV"])
                ),
                electron_density_m3=float(
                    row["reactor_electron_density_m3"]),
                delivered_bias_power_W=float(row["rf_bias_power_W"]),
            )
            surface = _combined_surface_projection(
                mechanism,
                boundary,
                sheath.distributions["Cl+"],
                sheath.distributions["Cl2+"],
            )
            depth = float(surface.etch_velocity_m_s * 75.0 * 1.0e9)
            row.update({
                "wafer_total_positive_ion_flux_m2_s": total_flux,
                "wafer_clplus_flux_m2_s": positive_flux["Cl+"],
                "wafer_cl2plus_flux_m2_s": positive_flux["Cl2+"],
                "wafer_clplus_fraction": positive_flux["Cl+"] / total_flux,
                "wafer_neutral_to_total_ion_ratio": (
                    float(row["wafer_atomic_chlorine_flux_m2_s"])
                    / total_flux
                ),
                "rf_sheath_bias_dc_component_v": sheath.bias_dc_component_v,
                "rf_sheath_dc_v": sheath.sheath_dc_v,
                "rf_sheath_amplitude_v": sheath.sheath_rf_amplitude_v,
                "rf_sheath_clplus_mean_energy_eV": (
                    sheath.distributions["Cl+"].mean_energy_eV),
                "rf_sheath_cl2plus_mean_energy_eV": (
                    sheath.distributions["Cl2+"].mean_energy_eV),
                "rf_sheath_power_closure_relative_residual": (
                    sheath.power_closure_relative_residual),
                "joined_species_resolved_surface_plane_depth_nm_75s": depth,
                "axisymmetric_source_moment": mode_name,
                "axisymmetric_raw_total_wafer_ion_flux_m2_s": float(sum(
                    spatial_mode["wafer_flux_m2_s"].values())),
                "axisymmetric_center_normalized_total_ion_flux_scale": (
                    spatial_mode[
                        "diagnostic_center_renormalized_total_ion_flux_scale"
                    ]),
                "axisymmetric_product_return_fields_valid": False,
                "formal_feature_depth_pass": False,
            })
            observed = row["observed_feature_depth_nm"]
            error = None if observed is None else 100.0 * (
                depth / float(observed) - 1.0)
            row["joined_species_resolved_signed_error_percent"] = error
            if error is not None:
                errors.append(error)
                positive_residual += error > 0.0
                negative_residual += error < 0.0
        summary = {
            "usable_depth_count": len(errors),
            "clean_surface_plane_mape_percent": float(
                np.mean(np.abs(errors))),
            "clean_surface_plane_signed_error_range_percent": [
                float(min(errors)), float(max(errors))],
            "overpredicted_condition_count": int(positive_residual),
            "underpredicted_condition_count": int(negative_residual),
            "formal_feature_depth_passes": 0,
            "feature_depth_used_to_select_source_moment": False,
        }
        receipt["axisymmetric_transport"]["summary"] = summary
        receipts[mode_name] = receipt
        summaries[mode_name] = summary
    base_errors = [
        float(row["joined_species_resolved_signed_error_percent"])
        for row in base["rows"]
        if row["joined_species_resolved_signed_error_percent"] is not None
    ]
    return {
        "schema": "petch.mahorowala_1998_axisymmetric_depth_reachability.v1",
        "reactor_receipt": _receipt_path(reactor_path),
        "reactor_receipt_sha256": _hash(reactor_path),
        "spatial_receipt": _receipt_path(spatial_path),
        "spatial_receipt_sha256": _hash(spatial_path),
        "base_global_edge_surface_plane_mape_percent": float(
            np.mean(np.abs(base_errors))),
        "source_mode_summary": summaries,
        "source_modes_selected_from_depth": [],
        "formal_feature_depth_passes": 0,
        "conclusion": (
            "The center-conditioned spatial tier can change sweep trends and "
            "ion composition, but no unmeasured source moment is an authorized "
            "absolute-depth calibration. Mixed-sign residuals prevent a single "
            "positive rate channel from closing the board globally."
        ),
        "next_measurement": (
            "radially and species-resolved Cl+/Cl2+ current/flux at the wafer "
            "plane across at least two source-power/flow conditions"
        ),
    }, receipts


def _attach_evolving_feature_receipts(result, output: Path):
    candidates = {
        "global_edge": (
            output / "feature_global_edge_clean_dx40nm_no_reflection.json"),
        "top_annular_wise_class": (
            output
            / "feature_top_annular_wise_class_clean_dx40nm_no_reflection.json"
        ),
        "top_center_compact": (
            output
            / "feature_top_center_compact_clean_dx40nm_no_reflection.json"
        ),
        "inductive_em_annular": (
            output
            / "feature_inductive_em_annular_clean_dx40nm_no_reflection.json"
        ),
    }
    board = {}
    for name, path in candidates.items():
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        errors = [
            float(row["signed_error_percent"])
            for row in payload["rows"]
            if row["signed_error_percent"] is not None
        ]
        board[name] = {
            "receipt": _receipt_path(path),
            "receipt_sha256": _hash(path),
            "run_count": len(errors),
            "mape_percent": float(payload["mape_percent"]),
            "signed_error_range_percent": [
                float(min(errors)), float(max(errors))],
            "overpredicted_condition_count": int(sum(
                error > 0.0 for error in errors)),
            "underpredicted_condition_count": int(sum(
                error < 0.0 for error in errors)),
            "dx_um": float(payload["dx_um"]),
            "product_wall_limit": payload["product_wall_limit"],
            "chlorine_specular_reflection": payload[
                "chlorine_specular_reflection"],
            "formal_feature_depth_pass": False,
        }
    if board:
        result["evolving_feature_board"] = board
        compared = [
            value["mape_percent"] for key, value in board.items()
            if key != "global_edge"
        ]
        if "global_edge" in board and compared:
            result["evolving_feature_board"][
                "maximum_spatial_mape_change_from_global_edge_points"
            ] = float(max(
                abs(value - board["global_edge"]["mape_percent"])
                for value in compared
            ))


def write(result, receipts, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    (output / "mahorowala_1998_axisymmetric_depth_reachability.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for mode_name, receipt in receipts.items():
        (output / f"reactor_receipt_{mode_name}.json").write_text(
            json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Mahorowala axisymmetric depth reachability",
        "",
        "The same independent 400 W/100 sccm current estimate normalizes every source moment once. No feature depth selects a source field. Product-return fields are invalidated; the emitted receipts are valid for the clean-surface/`--product-wall none` feature board.",
        "",
        f"- global-edge clean-surface MAPE: `{result['base_global_edge_surface_plane_mape_percent']:.2f}%`",
        "- formal feature-depth passes: `0`",
        "",
        "| source moment | clean-surface MAPE | signed-error range | over / under |",
        "|---|---:|---:|---:|",
    ]
    for mode_name, summary in result["source_mode_summary"].items():
        interval = summary[
            "clean_surface_plane_signed_error_range_percent"]
        lines.append(
            f"| {mode_name} | {summary['clean_surface_plane_mape_percent']:.2f}% | "
            f"{interval[0]:+.1f}% to {interval[1]:+.1f}% | "
            f"{summary['overpredicted_condition_count']} / "
            f"{summary['underpredicted_condition_count']} |"
        )
    feature_board = result.get("evolving_feature_board", {})
    feature_rows = {
        key: value for key, value in feature_board.items()
        if isinstance(value, dict)
    }
    if feature_rows:
        lines.extend((
            "",
            "## Evolving-feature check",
            "",
            "All cases use the same 40 nm grid, clean-surface boundary, no product return, and no ion reflection so only the reactor boundary changes.",
            "",
            "| boundary | feature MAPE | signed-error range | over / under |",
            "|---|---:|---:|---:|",
        ))
        for name, summary in feature_rows.items():
            interval = summary["signed_error_range_percent"]
            lines.append(
                f"| {name} | {summary['mape_percent']:.2f}% | "
                f"{interval[0]:+.1f}% to {interval[1]:+.1f}% | "
                f"{summary['overpredicted_condition_count']} / "
                f"{summary['underpredicted_condition_count']} |"
            )
        lines.extend((
            "",
            "The spatial boundary changes feature MAPE by at most "
            f"`{feature_board['maximum_spatial_mape_change_from_global_edge_points']:.2f}` percentage points on this controlled board; it does not remove the mixed-sign residual.",
        ))
    lines.extend((
        "",
        result["conclusion"],
        "",
        "Required discriminating measurement: " + result["next_measurement"] + ".",
        "",
    ))
    (output / "MAHOROWALA_1998_AXISYMMETRIC_DEPTH_REACHABILITY.md").write_text(
        "\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reactor", type=Path, default=DEFAULT_REACTOR)
    parser.add_argument("--spatial", type=Path, default=DEFAULT_SPATIAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    result, receipts = audit(arguments.reactor, arguments.spatial)
    _attach_evolving_feature_receipts(result, arguments.output)
    write(result, receipts, arguments.output)
    print(json.dumps(result["source_mode_summary"], sort_keys=True))


if __name__ == "__main__":
    main()
