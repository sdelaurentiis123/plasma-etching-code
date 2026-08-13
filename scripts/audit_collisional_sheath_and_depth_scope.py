#!/usr/bin/env python3
"""Freeze the deterministic collisional-sheath and Krueger depth-scope audit."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    DeterministicCollisionalRFSheath,
)
from petch.sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "results" / "curated" / "collisional_sheath_depth_scope"
    / "audit.json")
KRUEGER_DATA = ROOT / "data" / "experimental" / "krueger_2024"
KRUEGER_META = KRUEGER_DATA / "digitized_figure4_iead_metadata.json"
KRUEGER_IEAD = KRUEGER_DATA / "digitized_figure4_iead.csv"
DEPTH_AUDIT = (
    ROOT / "results" / "curated" / "depth_identifiability" / "audit.json")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _argon_sensitivity() -> dict[str, object]:
    waveform = PeriodicSheathVoltage.sinusoidal(
        dc_v=1000.0,
        amplitude_v=0.0,
        frequency_hz=2.0e6,
        source="declared 1 keV static Ar collisional sensitivity",
    )
    sheath = CollisionlessWaveformSheath(
        waveform=waveform,
        Te_eV=4.0,
        ion_mass_amu=39.948,
        thickness_m=0.01,
    )
    rows = []
    for pressure_mTorr in (0.0, 1.0, 5.0, 10.0):
        density = (
            pressure_mTorr * 0.13332236842105263
            / (1.380649e-23 * 500.0)
        )
        result = DeterministicCollisionalRFSheath(
            sheath=sheath,
            collision_model=ArgonBornMayerPhelpsCollisionModel(),
            gas_number_density_m3=density,
            neutral_gas_temperature_K=500.0,
            source_ion_flux_m2_s=1.0e19,
            phase_count=4,
            initial_thermal_radial_order=1,
            initial_thermal_azimuth_order=2,
            position_quadrature_order=4,
            hazard_quadrature_order=5,
            impact_quadrature_order=2,
            collision_azimuth_order=2,
            maximum_collision_order=2,
            steps_per_period=64,
            steps_per_transit=128,
        ).solve()
        rows.append({
            "pressure_mTorr": pressure_mTorr,
            "mean_optical_depth": result.mean_total_optical_depth,
            "resolved_ion_arrival_probability": (
                result.ion_arrival_probability),
            "unresolved_ion_probability": result.unresolved_probability,
            "resolved_ion_mean_energy_eV": (
                result.distribution.mean_energy_eV),
            "resolved_ion_rms_angle_deg": float(np.rad2deg(np.sqrt(
                result.distribution.mean_squared_polar_angle_rad2))),
            "resolved_fast_neutral_arrivals_per_source_ion": (
                result.resolved_fast_neutral_arrivals_per_source_ion),
            "unresolved_fast_neutral_collisions_per_source_ion": (
                result.unresolved_fast_neutral_collisions_per_source_ion),
            "fast_neutral_lineage_ledger_relative_residual": (
                result.fast_neutral_lineage_ledger_relative_residual),
            "maximum_resolved_energy_ledger_relative_residual": (
                result.maximum_resolved_energy_ledger_relative_residual),
            "probability_ledger_relative_residual": (
                result.probability_ledger_relative_residual),
        })
    convergence = {}
    for pressure_mTorr in (1.0, 10.0):
        density = (
            pressure_mTorr * 0.13332236842105263
            / (1.380649e-23 * 500.0)
        )
        order_rows = []
        for collision_order in (1, 2, 3):
            result = DeterministicCollisionalRFSheath(
                sheath=sheath,
                collision_model=ArgonBornMayerPhelpsCollisionModel(),
                gas_number_density_m3=density,
                neutral_gas_temperature_K=500.0,
                source_ion_flux_m2_s=1.0e19,
                phase_count=4,
                initial_thermal_radial_order=1,
                initial_thermal_azimuth_order=2,
                position_quadrature_order=3,
                hazard_quadrature_order=4,
                impact_quadrature_order=2,
                collision_azimuth_order=2,
                maximum_collision_order=collision_order,
                steps_per_period=64,
                steps_per_transit=128,
            ).solve()
            order_rows.append({
                "collision_order": collision_order,
                "unresolved_ion_probability": result.unresolved_probability,
                "resolved_ion_mean_energy_eV": (
                    result.distribution.mean_energy_eV),
                "resolved_fast_neutral_arrivals_per_source_ion": (
                    result.resolved_fast_neutral_arrivals_per_source_ion),
            })
        relative_change = abs(
            order_rows[-1]["resolved_ion_mean_energy_eV"]
            / order_rows[-2]["resolved_ion_mean_energy_eV"]
            - 1.0
        )
        converged = (
            order_rows[-1]["unresolved_ion_probability"] < 1.0e-3
            and relative_change < 5.0e-3
        )
        convergence[f"{pressure_mTorr:g}_mTorr"] = {
            "rows": order_rows,
            "order2_to_order3_mean_energy_relative_change": relative_change,
            "gate": {
                "unresolved_ion_probability_limit": 1.0e-3,
                "mean_energy_relative_change_limit": 5.0e-3,
                "passed": converged,
            },
            "status": "converged" if converged else "not_converged",
        }
    return {
        "condition": {
            "gas": "pure Ar",
            "sheath_voltage_V": 1000.0,
            "sheath_thickness_m": 0.01,
            "gas_temperature_K": 500.0,
            "ion_energy_angle_target_used": None,
            "feature_depth_target_used": None,
        },
        "rows": rows,
        "collision_order_convergence_probe": {
            "conditions": convergence,
            "quadrature": {
                "phase_count": 4,
                "initial_thermal_radial_order": 1,
                "initial_thermal_azimuth_order": 2,
                "position_quadrature_order": 3,
                "hazard_quadrature_order": 4,
                "impact_quadrature_order": 2,
                "collision_azimuth_order": 2,
            },
        },
    }


def build_audit() -> dict[str, object]:
    metadata = json.loads(KRUEGER_META.read_text())
    depth = json.loads(DEPTH_AUDIT.read_text())
    exact_pdf_hash = "65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51"
    exact_image_hash = "99ac8dd2f916071b51de67c729e131a6e20d3b92688d92693d1a467cd0782f9c"
    if (
        metadata["source_pdf_sha256"] != exact_pdf_hash
        or metadata["embedded_image_sha256"] != exact_image_hash
    ):
        raise ValueError("Krueger IEAD source-pixel receipt changed")
    return {
        "schema": "petch.collisional-sheath-depth-scope.v1",
        "audit_id": "COLLISIONAL-SHEATH-DEPTH-SCOPE-2026-08-12",
        "source_receipts": {
            "krueger_iead_table_sha256": _sha(KRUEGER_IEAD),
            "krueger_iead_metadata_sha256": _sha(KRUEGER_META),
            "krueger_author_pdf_sha256": exact_pdf_hash,
            "krueger_embedded_figure_png_sha256": exact_image_hash,
            "depth_identifiability_audit_sha256": _sha(DEPTH_AUDIT),
        },
        "digitization_visual_gate": {
            "author_pdf_re_rendered_at_dpi": 400,
            "visual_result": (
                "Figure 4a axes, 0--5000 eV support, -5--5 degree angle "
                "support, narrow central plume, low-energy lobe, and "
                "high-energy maximum agree with the archived calibration"),
            "exact_source_pdf_hash_match": True,
            "exact_embedded_image_file_hash_match": True,
            "mean_energy_eV": metadata["resampled_summary"]["mean_energy_eV"],
            "angle_standard_deviation_deg": metadata[
                "resampled_summary"]["angle_standard_deviation_deg"],
            "claim": (
                "the checksum-bound table is a faithful deterministic "
                "digitization of the archived HPEM image; it remains model "
                "output, not a wafer measurement"),
        },
        "argon_collision_operator": _argon_sensitivity(),
        "krueger_depth_impact": {
            "published_feature_boundary_kind": (
                "HPEM wafer-plane combined-positive-ion IEAD and aggregate flux"),
            "sheath_collisions_already_upstream_of_boundary": True,
            "postprocessing_with_new_sheath_would_double_count": True,
            "argon_cross_sections_transferable_to_CFx_positive_ions": False,
            "existing_simulated_depth_nm": depth["inputs"]["simulated_depth_nm"],
            "experimental_depth_nm": depth["inputs"]["target_depth_nm"],
            "depth_changed_by_this_audit": False,
            "exact_depth_prediction_authorized": False,
            "remaining_boundary_requirements": [
                "species-resolved positive-ion fluxes and IEADs",
                "stable C4F6 wafer flux",
                "species-resolved ion-neutral elastic/CX cross sections",
                "measured or validated reactor voltage/current waveform",
            ],
        },
        "verdict": (
            "The deterministic collisional sheath is now available for "
            "source-supported pure-Ar knobs-to-wafer prediction and exposes "
            "an unscattered fast-neutral lower-bound boundary. It does not "
            "repair Krueger depth: the paper already supplies a downstream "
            "HPEM wafer IEAD, and the missing molecular-ion composition and "
            "stable-parent boundary remain depth-identifying data."),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    audit = build_audit()
    payload = json.dumps(audit, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != payload:
            raise SystemExit("collisional-sheath depth-scope audit is stale")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(payload)


if __name__ == "__main__":
    main()
