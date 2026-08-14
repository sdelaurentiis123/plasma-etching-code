#!/usr/bin/env python3
"""Freeze a deterministic global-Ar -> moving-RF-sheath wafer audit."""
from __future__ import annotations

import argparse
from datetime import date
import json
from pathlib import Path

import numpy as np

from petch.reactor_global import (
    ArgonBornMayerPhelpsCollisionModel,
    ArgonCurrentDrivenReactorToWaferCondition,
    ArgonGlobalCondition,
    ArgonTransportState,
    CylindricalReactor,
    DeterministicCurrentDrivenArgonReactorToWaferModel,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    PASCAL_PER_MTORR,
    PeriodicCurrentDensity,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "current_driven_argon_reactor_stack"
    / "audit.json"
)


def _model(
    *,
    phase_node_count: int = 24,
    steps_per_period: int = 64,
    steps_per_transit: int = 64,
) -> DeterministicCurrentDrivenArgonReactorToWaferModel:
    transport = FixedArgonTransportProvider(ArgonTransportState(
        ion_mean_free_path_m=0.004,
        ambipolar_diffusion_m2_s=50.0,
        metastable_effective_diffusion_m2_s=1.0,
        source="manufactured Ar transport for numerical closure audit",
        evidence_kind="assumed",
    ))
    return DeterministicCurrentDrivenArgonReactorToWaferModel(
        global_model=LeeLiebermanArgonGlobalModel(transport),
        collision_model=ArgonBornMayerPhelpsCollisionModel(),
        sheath_phase_quadrature_count=2048,
        phase_node_count=phase_node_count,
        position_node_count=4,
        total_energy_node_count=4,
        transverse_fraction_node_count=5,
        initial_thermal_radial_order=1,
        output_azimuth_order=4,
        impact_quadrature_order=2,
        collision_azimuth_order=4,
        collision_event_quadrature_order=3,
        steps_per_period=steps_per_period,
        steps_per_transit=steps_per_transit,
        maximum_transit_periods=16.0,
    )


def _condition() -> ArgonCurrentDrivenReactorToWaferCondition:
    return ArgonCurrentDrivenReactorToWaferCondition(
        global_condition=ArgonGlobalCondition(
            condition_id="manufactured-current-driven-Ar-audit",
            absorbed_power_W=500.0,
            pressure_Pa=2.0 * PASCAL_PER_MTORR,
            gas_temperature_K=500.0,
            geometry=CylindricalReactor(radius_m=0.15, length_m=0.075),
            ion_wall_energy_factor_Te=5.0,
            ion_wall_energy_source=(
                "Lee-Lieberman published-range member; numerical audit"),
            ion_wall_energy_evidence="published_range_member",
            absorbed_power_source="manufactured numerical-audit scalar",
            absorbed_power_evidence="assumed",
            absorbed_power_boundary_kind="manufactured_test",
        ),
        sheath_current_density=PeriodicCurrentDensity(
            fundamental_frequency_hz=2.0e6,
            harmonic_number=np.array([1, 2]),
            sine_A_m2=np.array([-40.0, 4.5]),
            cosine_A_m2=np.array([0.0, -2.5]),
            source="manufactured two-harmonic sheath current",
            evidence_kind="assumed",
        ),
    )


def build_receipt() -> dict[str, object]:
    solution = _model().solve(_condition())
    coarse = _model(
        steps_per_period=48, steps_per_transit=48).solve(_condition())
    fine = _model(
        steps_per_period=96, steps_per_transit=96).solve(_condition())
    distribution = solution.wafer.distribution
    tangent = solution.sheath.current_scale_jvp(1.0)

    def relative(left: float, right: float) -> float:
        return abs(right - left) / max(abs(left), abs(right), 1.0e-300)

    coarse_rms = float(np.rad2deg(np.sqrt(
        coarse.wafer.distribution.mean_squared_polar_angle_rad2)))
    fine_rms = float(np.rad2deg(np.sqrt(
        fine.wafer.distribution.mean_squared_polar_angle_rad2)))
    temporal_changes = {
        "mean_energy_relative_change": relative(
            coarse.wafer.distribution.mean_energy_eV,
            fine.wafer.distribution.mean_energy_eV),
        "rms_angle_relative_change": relative(coarse_rms, fine_rms),
        "collision_count_relative_change": relative(
            coarse.wafer.expected_collision_count_lower_bound,
            fine.wafer.expected_collision_count_lower_bound),
        "arrival_probability_relative_change": relative(
            coarse.wafer.ion_arrival_probability,
            fine.wafer.ion_arrival_probability),
    }
    temporal_limits = {
        "mean_energy_relative_change": 2.0e-3,
        "rms_angle_relative_change": 2.0e-3,
        "collision_count_relative_change": 2.0e-3,
        "arrival_probability_relative_change": 2.0e-3,
    }
    temporal_passed = all(
        temporal_changes[key] <= temporal_limits[key]
        for key in temporal_limits
    )
    phase_coarse = _model(phase_node_count=12).solve(_condition())
    # The reported solution already is the 24-phase-node member.
    phase_fine = solution
    phase_coarse_rms = float(np.rad2deg(np.sqrt(
        phase_coarse.wafer.distribution.mean_squared_polar_angle_rad2)))
    phase_fine_rms = float(np.rad2deg(np.sqrt(
        phase_fine.wafer.distribution.mean_squared_polar_angle_rad2)))
    phase_changes = {
        "mean_energy_relative_change": relative(
            phase_coarse.wafer.distribution.mean_energy_eV,
            phase_fine.wafer.distribution.mean_energy_eV),
        "rms_angle_relative_change": relative(
            phase_coarse_rms, phase_fine_rms),
        "collision_count_relative_change": relative(
            phase_coarse.wafer.expected_collision_count_lower_bound,
            phase_fine.wafer.expected_collision_count_lower_bound),
        "arrival_probability_relative_change": relative(
            phase_coarse.wafer.ion_arrival_probability,
            phase_fine.wafer.ion_arrival_probability),
    }
    phase_limits = {
        "mean_energy_relative_change": 1.0e-2,
        "rms_angle_relative_change": 2.0e-2,
        "collision_count_relative_change": 1.0e-2,
        "arrival_probability_relative_change": 1.0e-2,
    }
    phase_passed = all(
        phase_changes[key] <= phase_limits[key] for key in phase_limits)
    receipt = {
        "schema": "petch.current_driven_argon_reactor_stack.v1",
        "generated_on": str(date.today()),
        "purpose": (
            "numerical and conservation audit; no equipment or depth target "
            "was fitted or used"
        ),
        "inputs": {
            "absorbed_bulk_power_W": 500.0,
            "pressure_mTorr": 2.0,
            "gas_temperature_K": 500.0,
            "reactor_radius_m": 0.15,
            "reactor_length_m": 0.075,
            "sheath_current_frequency_hz": 2.0e6,
            "sheath_current_harmonics_A_m2": {
                "sine": [-40.0, 4.5],
                "cosine": [0.0, -2.5],
            },
            "feature_depth_target_used": None,
            "ion_energy_angle_target_used": None,
        },
        "global_plasma": {
            "electron_temperature_eV": (
                solution.global_plasma.electron_temperature_eV),
            "electron_density_m3": (
                solution.global_plasma.electron_density_m3),
            "axial_ion_flux_m2_s": (
                solution.global_plasma.axial_ion_flux_m2_s),
            "maximum_normalized_balance_residual": (
                solution.global_plasma.maximum_normalized_residual),
        },
        "bohm_flux_seam": {
            "sheath_edge_density_m3": solution.sheath_edge_density_m3,
            "relative_residual": solution.bohm_flux_seam_relative_residual,
        },
        "moving_rf_sheath": {
            "xi": solution.sheath.xi,
            "maximum_width_m": solution.sheath.maximum_width_m,
            "maximum_voltage_v": solution.sheath.maximum_voltage_v,
            "mean_voltage_v": solution.sheath.mean_voltage_v,
            "child_current_relative_residual": (
                solution.sheath.child_current_relative_residual),
            "charge_voltage_relative_residual": (
                solution.sheath.charge_voltage_relative_residual),
            "exact_uniform_current_scale_jvp": {
                "d_width_d_log_current_m": (
                    tangent.maximum_width_tangent_m),
                "d_voltage_d_log_current_v": (
                    tangent.maximum_voltage_tangent_v),
            },
        },
        "collisional_wafer_boundary": {
            "source_ion_flux_m2_s": solution.wafer.source_ion_flux_m2_s,
            "arriving_ion_flux_m2_s": solution.wafer.arriving_ion_flux_m2_s,
            "ion_arrival_probability": (
                solution.wafer.ion_arrival_probability),
            "mean_impact_energy_eV": distribution.mean_energy_eV,
            "rms_impact_angle_deg": float(np.rad2deg(np.sqrt(
                distribution.mean_squared_polar_angle_rad2))),
            "expected_collision_count": (
                solution.wafer.expected_collision_count_lower_bound),
            "resolved_fast_neutral_flux_lower_bound_m2_s": (
                solution.wafer.resolved_fast_neutral_flux_m2_s),
            "unresolved_fast_neutral_collisions_per_source_ion": (
                solution.wafer
                .unresolved_fast_neutral_collisions_per_source_ion),
            "low_energy_phelps_collision_count_per_source_ion": (
                solution.wafer
                .below_born_mayer_support_collision_probability_lower_bound),
            "probability_ledger_relative_residual": (
                solution.wafer.probability_ledger_relative_residual),
            "maximum_collision_energy_ledger_relative_residual": (
                solution.wafer
                .maximum_resolved_energy_ledger_relative_residual),
            "linear_solve_relative_residual": (
                solution.wafer.provenance[
                    "linear_solve_relative_residual"]),
        },
        "temporal_characteristic_convergence": {
            "coarse_steps_per_period_and_transit": 48,
            "reported_steps_per_period_and_transit": 64,
            "fine_steps_per_period_and_transit": 96,
            "coarse_to_fine_relative_changes": temporal_changes,
            "limits": temporal_limits,
            "passed": temporal_passed,
            "note": (
                "This certifies the time characteristic at a fixed kinetic "
                "ordinate grid. A production prediction must additionally "
                "carry a condition-specific ordinate-grid receipt."
            ),
        },
        "rf_phase_ordinate_convergence": {
            "coarse_phase_node_count": 12,
            "reported_phase_node_count": 24,
            "coarse_to_reported_relative_changes": phase_changes,
            "limits": phase_limits,
            "passed": phase_passed,
            "angle_is_numerically_gradeable": True,
            "angle_is_equipment_validated": False,
            "angle_gate_blocker": (
                "none for the declared Phelps/LXCat low-energy model; an "
                "independent measured IEAD remains required for equipment "
                "validation"
            ),
        },
        "support_gates": {
            "bulk_condition_is_equipment_evidenced": (
                solution.global_plasma.supports_prediction),
            "sheath_current_is_measured_or_circuit_validated": (
                solution.sheath.current.supports_predictive_boundary),
            "resolved_ion_boundary_predictive": (
                solution.supports_resolved_ion_boundary_prediction),
            "moving_electron_front_resolved": True,
            "time_dependent_poisson_field_resolved": True,
            "RF_phase_is_kinetic_state": True,
            "all_ion_collision_orders_closed": True,
            "low_energy_ion_angular_scattering_closed": True,
            "low_energy_ion_angular_model": (
                "Phelps LXCat isotropic/backscatter decomposition"
            ),
            "repeated_fast_neutral_transport_closed": False,
            "generator_matching_network_inversion_closed": False,
            "molecular_chemistry_closed": False,
            "feature_surface_law_closed": False,
            "supports_equipment_prediction": False,
            "supports_feature_depth": False,
        },
        "maximum_conservation_residual": (
            solution.maximum_conservation_residual),
        "source": {
            "moving_sheath": (
                "Turner & Chabert, APL 104, 164102 (2014), "
                "DOI 10.1063/1.4872172, equations 1-19"
            ),
            "argon_scattering": solution.wafer.model_source,
        },
    }
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
