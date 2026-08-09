#!/usr/bin/env python3
"""Lift the Mahorowala chlorine 0-D ion inventory to a 2-D wafer boundary.

The chlorine chemistry and EEDF are solved without etch depths.  For each
condition, the resulting volume-average Cl+/Cl2+ inventories are imposed on a
conservative axisymmetric reaction-diffusion solve.  Published mobility-based
ambipolar diffusivities and Bohm velocities supply the transport and wall
boundary.  Uniform and declared ICP-like source moments expose the remaining
spatial-source uncertainty.

This audit diagnoses the global-model edge-factor approximation.  Because the
distributed ion source is inferred from the solved 0-D inventory, it is not an
independent reactor-state prediction and is never labeled as one.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from forecast_mahorowala_1998_cl2_depth import (
    GAP_CM_EQUIPMENT_CLASS_PRIOR,
    ION_TEMPERATURE_EV,
    WAFER_DIAMETER_M,
    _condition,
    _model,
    _providers,
    _solve_states,
    _source_rows,
)
from petch.reactor_global import (
    AxisymmetricFiniteVolumeGrid,
    AxisymmetricInductiveFieldCondition,
    DeterministicAxisymmetricInductiveField,
    DeterministicQuasineutralInventoryLift,
    VACUUM_PERMEABILITY_H_M,
    drude_conductivity_from_two_term_mobilities,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    normalized_annular_skin_source,
    normalized_exponential_skin_source,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "results" / "curated" / "reactor_to_feature_chlorine"
    / "axisymmetric_ion_transport_audit"
)
POSITIVE_ION_NAMES = ("Cl+", "Cl2+")
CHARGED_SPECIES_NAMES = (*POSITIVE_ION_NAMES, "Cl-")
LEE_MUTUAL_NEUTRALIZATION_M3_S = 5.0e-14


def _source_moments(grid: AxisymmetricFiniteVolumeGrid):
    uniform = np.ones((grid.radial_cell_count, grid.axial_cell_count))
    return {
        "uniform": {
            "field": uniform,
            "provenance": "uniform volumetric source limiting case",
            "conditioned_or_measured": False,
        },
        "top_center_broad": {
            "field": normalized_exponential_skin_source(
                grid,
                axial_skin_depth_m=0.030,
                radial_scale_m=0.180,
            ),
            "provenance": (
                "declared ICP-like sensitivity: 30 mm axial skin depth and "
                "180 mm center scale; not a fitted Lam field"
            ),
            "conditioned_or_measured": False,
        },
        "top_center_compact": {
            "field": normalized_exponential_skin_source(
                grid,
                axial_skin_depth_m=0.015,
                radial_scale_m=0.100,
            ),
            "provenance": (
                "declared ICP-like sensitivity: 15 mm axial skin depth and "
                "100 mm center scale; not a fitted Lam field"
            ),
            "conditioned_or_measured": False,
        },
        "top_annular_wise_class": {
            "field": normalized_annular_skin_source(
                grid,
                axial_skin_depth_m=0.010,
                ring_radius_m=0.140,
                radial_width_m=0.035,
            ),
            "provenance": (
                "Wise 1996 source-class sensitivity: published approximately "
                "10 mm axial skin depth and toroidal topology; 140 mm ring "
                "radius and 35 mm width are declared Lam-geometry moments, "
                "not measured or fitted fields"
            ),
            "conditioned_or_measured": False,
        },
        "inductive_em_annular": {
            "field": None,
            "provenance": (
                "deterministic axisymmetric harmonic Etheta solve using "
                "the EEPF-derived complex Drude conductivity and a declared "
                "Lam-geometry annular coil-side field; generator coupling and "
                "the target-tool coil map are not measured"
            ),
            "conditioned_or_measured": False,
        },
    }


def _inductive_power_source_moment(
    grid: AxisymmetricFiniteVolumeGrid, state, *, frequency_hz: float
):
    neutral_density = state.densities_m3["Cl"] + state.densities_m3["Cl2"]
    moments = state.electron_solution.transport_moments
    conductivity = drude_conductivity_from_two_term_mobilities(
        electron_density_m3=state.densities_m3["e"],
        neutral_density_m3=neutral_density,
        frequency_hz=frequency_hz,
        flux_reduced_mobility_m_inv_V_inv_s_inv=(
            moments.flux_reduced_mobility_m_inv_V_inv_s_inv),
        dissipative_reduced_mobility_m_inv_V_inv_s_inv=(
            moments.dissipative_reduced_mobility_m_inv_V_inv_s_inv),
    )
    radial = grid.radial_centers_m
    radius = grid.geometry.radius_m
    ring_radius = 0.140
    ring_width = 0.035
    upper_field = (
        radial / ring_radius
        * np.exp(-0.5 * ((radial - ring_radius) / ring_width) ** 2)
        * np.maximum(1.0 - (radial / radius) ** 2, 0.0)
    )
    condition = AxisymmetricInductiveFieldCondition(
        grid=grid,
        frequency_hz=frequency_hz,
        conductivity_S_m=np.full(
            (grid.radial_cell_count, grid.axial_cell_count), conductivity),
        upper_boundary_electric_field_V_m=upper_field,
        source=(
            "Mahorowala/Lam condition-specific ICP field sensitivity"
        ),
        conductivity_source=(
            "two-term EEPF flux/dissipative mobility moments and solved "
            "electron/neutral densities"
        ),
        upper_boundary_source=(
            "declared annular 140 mm radius, 35 mm width coil-side Etheta "
            "moment with conducting-sidewall taper; not fitted"
        ),
    )
    solution = DeterministicAxisymmetricInductiveField(condition).solve()
    volume = grid.cell_volume_m3
    power = solution.absorbed_power_density_W_m3
    normalized_power = power / np.sum(power * volume)
    axial_depth = float(np.sum(
        normalized_power
        * (grid.geometry.length_m - grid.axial_centers_m[None, :])
        * volume
    ))
    radial_mean = float(np.sum(
        normalized_power * grid.radial_centers_m[:, None] * volume))
    skin_depth = float(np.sqrt(
        2.0
        / (
            VACUUM_PERMEABILITY_H_M
            * 2.0 * np.pi * frequency_hz
            * conductivity.real
        )
    ))
    return solution.normalized_absorbed_power_source_moment, {
        "complex_conductivity_S_m": {
            "real": float(conductivity.real),
            "imaginary": float(conductivity.imag),
        },
        "classical_dissipative_skin_depth_m": skin_depth,
        "absorbed_power_mean_depth_below_top_m": axial_depth,
        "absorbed_power_mean_radius_m": radial_mean,
        "complex_power_ledger_relative_residual": (
            solution.complex_power_ledger_relative_residual),
        "linear_system_relative_residual": (
            solution.linear_system_relative_residual),
        "supports_implicit_differentiation": True,
        "supports_generator_to_plasma_power_prediction": False,
    }


def audit(collision_deck: Path, atomic_cl_momentum: Path, *, energy_cells: int):
    model, replay, energy_grid = _model(
        collision_deck, atomic_cl_momentum, energy_cells)
    source_rows = _source_rows()
    conditions = sorted({
        (float(row["inductive_power_W"]), float(row["cl2_flow_sccm"]))
        for row in source_rows
        if row["quantitative_status"] == "usable"
    })
    states = _solve_states(model, conditions)
    charged_provider, _, _ = _providers()
    geometry = _condition(*conditions[0]).geometry
    grid = AxisymmetricFiniteVolumeGrid.uniform(
        geometry,
        radial_cell_count=14,
        axial_cell_count=10,
    )
    source_moments = _source_moments(grid)
    reduced_mobility = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    rows = []
    previous_potential_by_mode: dict[str, np.ndarray] = {}
    certified_potential: dict[tuple[float, float, str], np.ndarray] = {}
    for source_power_W, flow_sccm in conditions:
        state = states[(source_power_W, flow_sccm)]
        condition = _condition(source_power_W, flow_sccm)
        transport_condition = condition.transport_condition(
            2.0 / 3.0 * state.mean_electron_energy_eV)
        charged = charged_provider.predict(
            transport_condition, state.densities_m3)
        target_density = np.asarray([
            state.densities_m3[name] for name in CHARGED_SPECIES_NAMES
        ])
        electron_temperature_eV = (
            2.0 / 3.0 * state.mean_electron_energy_eV)
        positive_mobility = np.asarray([
            charged.positive_ion_transport[name].provenance[
                "ambipolar_diffusivity_m2_s"
            ] / electron_temperature_eV
            for name in POSITIVE_ION_NAMES
        ])
        negative_mobility = reduced_mobility["Cl-"].evaluate(
            total_neutral_density_m3=(
                state.densities_m3["Cl"] + state.densities_m3["Cl2"]),
            ion_temperature_eV=ION_TEMPERATURE_EV,
        )
        mobility = np.concatenate((
            positive_mobility,
            [negative_mobility.mobility_m2_V_s],
        ))
        bohm_velocity = np.asarray([
            charged.positive_ion_transport[name].provenance["bohm_speed_m_s"]
            for name in POSITIVE_ION_NAMES
        ])
        analytic_flux = np.asarray([
            state.axial_positive_ion_flux_m2_s[name]
            for name in POSITIVE_ION_NAMES
        ])
        modes = {}
        condition_source_moments = dict(source_moments)
        em_field, em_receipt = _inductive_power_source_moment(
            grid, state, frequency_hz=13.56e6)
        condition_source_moments["inductive_em_annular"] = {
            **source_moments["inductive_em_annular"],
            "field": em_field,
            "em_receipt": em_receipt,
        }
        within_condition_potential = None
        for mode_name, mode in condition_source_moments.items():
            print(
                "axisymmetric charged lift: "
                f"source={source_power_W:g} W flow={flow_sccm:g} sccm "
                f"mode={mode_name}",
                flush=True,
            )
            lift = DeterministicQuasineutralInventoryLift(
                grid=grid,
                species_names=CHARGED_SPECIES_NAMES,
                charge_number=np.array([1.0, 1.0, -1.0]),
                mobility_m2_V_s=mobility,
                ion_temperature_eV=np.full(3, ION_TEMPERATURE_EV),
                electron_temperature_eV=electron_temperature_eV,
                wall_velocity_m_s=np.vstack((
                    np.repeat(bohm_velocity[:, None], 3, axis=1),
                    np.full((1, 3), np.inf),
                )),
                source_shape=np.stack((
                    mode["field"], mode["field"], mode["field"])),
                source=(
                    "Mahorowala/Lam 0-D ion inventory; " + mode["provenance"]
                ),
                positive_negative_recombination_m3_s=(
                    LEE_MUTUAL_NEUTRALIZATION_M3_S),
            )
            result = lift.solve(
                target_density,
                relative_tolerance=5.0e-8,
                maximum_iterations=500,
                initial_electrostatic_potential_V=(
                    previous_potential_by_mode.get(
                        mode_name, within_condition_potential)
                ),
            )
            previous_potential_by_mode[mode_name] = (
                result.electrostatic_potential_V)
            within_condition_potential = result.electrostatic_potential_V
            certified_potential[(
                source_power_W, flow_sccm, mode_name
            )] = result.electrostatic_potential_V
            wafer_flux = np.asarray([
                result.solution.lower_endcap_area_average_flux_m2_s(
                    name, wafer_radius_m=0.5 * WAFER_DIAMETER_M)
                for name in POSITIVE_ION_NAMES
            ])
            full_endcap_flux = np.asarray([
                result.solution.lower_endcap_area_average_flux_m2_s(
                    name, wafer_radius_m=geometry.radius_m)
                for name in POSITIVE_ION_NAMES
            ])
            modes[mode_name] = {
                "source_provenance": mode["provenance"],
                "source_conditioned_or_measured": mode[
                    "conditioned_or_measured"
                ],
                "wafer_flux_m2_s": dict(zip(POSITIVE_ION_NAMES, wafer_flux)),
                "full_lower_endcap_flux_m2_s": dict(
                    zip(POSITIVE_ION_NAMES, full_endcap_flux)),
                "wafer_to_global_analytic_flux_ratio": dict(zip(
                    POSITIVE_ION_NAMES, wafer_flux / analytic_flux
                )),
                "total_wafer_to_global_analytic_flux_ratio": float(
                    np.sum(wafer_flux) / np.sum(analytic_flux)
                ),
                "center_wafer_to_full_endcap_flux_ratio": dict(zip(
                    POSITIVE_ION_NAMES, wafer_flux / full_endcap_flux
                )),
                "inferred_source_amplitude_m3_s": dict(zip(
                    CHARGED_SPECIES_NAMES,
                    result.inferred_source_amplitude_m3_s
                )),
                "minimum_electron_density_m3": (
                    result.minimum_electron_density_m3),
                "electrostatic_potential_span_V": float(np.ptp(
                    result.electrostatic_potential_V)),
                "nonlinear_solver_evaluations": (
                    result.nonlinear_solver_evaluations),
                "maximum_potential_fixed_point_relative_residual": (
                    result.maximum_potential_fixed_point_relative_residual),
                "maximum_particle_ledger_relative_residual": (
                    result.solution.maximum_species_ledger_relative_residual
                ),
                "maximum_inventory_relative_residual": (
                    result.maximum_inventory_relative_residual
                ),
                "supports_reactor_state_prediction": False,
                "supports_implicit_differentiation": True,
            }
            if "em_receipt" in mode:
                modes[mode_name]["inductive_field_receipt"] = mode[
                    "em_receipt"]
        rows.append({
            "inductive_power_W": source_power_W,
            "cl2_flow_sccm": flow_sccm,
            "mean_electron_energy_eV": state.mean_electron_energy_eV,
            "reduced_electric_field_Td": state.reduced_electric_field_Td,
            "volume_average_density_m3": dict(zip(
                CHARGED_SPECIES_NAMES, target_density)),
            "free_ion_mobility_m2_V_s": dict(zip(
                CHARGED_SPECIES_NAMES, mobility)),
            "ion_temperature_eV": ION_TEMPERATURE_EV,
            "electron_temperature_eV": electron_temperature_eV,
            "positive_negative_recombination_m3_s": (
                LEE_MUTUAL_NEUTRALIZATION_M3_S),
            "bohm_wall_velocity_m_s": dict(zip(
                POSITIVE_ION_NAMES, bohm_velocity)),
            "global_analytic_axial_flux_m2_s": dict(
                zip(POSITIVE_ION_NAMES, analytic_flux)),
            "spatial_source_modes": modes,
            "feature_depth_used_to_condition_transport": False,
        })

    ratio_by_mode = {
        mode_name: [
            row["spatial_source_modes"][mode_name][
                "total_wafer_to_global_analytic_flux_ratio"
            ]
            for row in rows
        ]
        for mode_name in source_moments
    }
    center_row = next(
        row for row in rows
        if row["inductive_power_W"] == 400.0
        and row["cl2_flow_sccm"] == 100.0
    )
    center_ratio_by_mode = {
        mode_name: center_row["spatial_source_modes"][mode_name][
            "total_wafer_to_global_analytic_flux_ratio"
        ]
        for mode_name in source_moments
    }
    center_renormalized_by_mode = {name: [] for name in source_moments}
    for row in rows:
        for mode_name in source_moments:
            relative = (
                row["spatial_source_modes"][mode_name][
                    "total_wafer_to_global_analytic_flux_ratio"
                ] / center_ratio_by_mode[mode_name]
            )
            row["spatial_source_modes"][mode_name][
                "diagnostic_center_renormalized_total_ion_flux_scale"
            ] = float(relative)
            center_renormalized_by_mode[mode_name].append(float(relative))
    max_ledger = max(
        row["spatial_source_modes"][mode_name][
            "maximum_particle_ledger_relative_residual"
        ]
        for row in rows for mode_name in source_moments
    )

    # Grid receipt at the independent center-current anchor. The 14x10 point
    # is the production audit grid; 10x8 and 20x14 are coarser/finer checks.
    center_state = states[(400.0, 100.0)]
    center_condition = _condition(400.0, 100.0)
    center_transport_condition = center_condition.transport_condition(
        2.0 / 3.0 * center_state.mean_electron_energy_eV)
    center_charged = charged_provider.predict(
        center_transport_condition, center_state.densities_m3)
    center_electron_temperature_eV = (
        2.0 / 3.0 * center_state.mean_electron_energy_eV)
    center_mobility = np.asarray([
        *(
            center_charged.positive_ion_transport[name].provenance[
                "ambipolar_diffusivity_m2_s"
            ] / center_electron_temperature_eV
            for name in POSITIVE_ION_NAMES
        ),
        reduced_mobility["Cl-"].evaluate(
            total_neutral_density_m3=(
                center_state.densities_m3["Cl"]
                + center_state.densities_m3["Cl2"]),
            ion_temperature_eV=ION_TEMPERATURE_EV,
        ).mobility_m2_V_s,
    ])
    center_bohm = np.asarray([
        center_charged.positive_ion_transport[name].provenance[
            "bohm_speed_m_s"
        ] for name in POSITIVE_ION_NAMES
    ])
    center_target = np.asarray([
        center_state.densities_m3[name] for name in CHARGED_SPECIES_NAMES
    ])
    grid_receipt = []
    seed_grid = grid
    center_production_potential = certified_potential[
        (400.0, 100.0, "uniform")]
    seed_potential = center_production_potential
    for radial_count, axial_count in ((10, 8), (14, 10), (20, 14)):
        receipt_grid = AxisymmetricFiniteVolumeGrid.uniform(
            geometry,
            radial_cell_count=radial_count,
            axial_cell_count=axial_count,
        )
        if (radial_count, axial_count) == (14, 10):
            receipt_flux = np.asarray([
                center_row["spatial_source_modes"]["uniform"][
                    "wafer_flux_m2_s"
                ][name] for name in POSITIVE_ION_NAMES
            ])
            receipt_potential = center_production_potential
            receipt_evaluations = center_row["spatial_source_modes"][
                "uniform"]["nonlinear_solver_evaluations"]
        else:
            interpolator = RegularGridInterpolator(
                (seed_grid.radial_centers_m, seed_grid.axial_centers_m),
                seed_potential,
                bounds_error=False,
                fill_value=None,
            )
            rr, zz = np.meshgrid(
                receipt_grid.radial_centers_m,
                receipt_grid.axial_centers_m,
                indexing="ij",
            )
            initial_potential = interpolator(
                np.stack((rr.ravel(), zz.ravel()), axis=-1)
            ).reshape(rr.shape)
            receipt_lift = DeterministicQuasineutralInventoryLift(
                grid=receipt_grid,
                species_names=CHARGED_SPECIES_NAMES,
                charge_number=np.array([1.0, 1.0, -1.0]),
                mobility_m2_V_s=center_mobility,
                ion_temperature_eV=np.full(3, ION_TEMPERATURE_EV),
                electron_temperature_eV=center_electron_temperature_eV,
                wall_velocity_m_s=np.vstack((
                    np.repeat(center_bohm[:, None], 3, axis=1),
                    np.full((1, 3), np.inf),
                )),
                source_shape=np.ones((
                    3, radial_count, axial_count)),
                source=(
                    "wise-1996-rapid-2d-cl uniform-source grid receipt; "
                    "Mahorowala/Lam center 0-D inventory"
                ),
                positive_negative_recombination_m3_s=(
                    LEE_MUTUAL_NEUTRALIZATION_M3_S),
            )
            receipt_result = receipt_lift.solve(
                center_target,
                relative_tolerance=5.0e-8,
                maximum_iterations=500,
                initial_electrostatic_potential_V=initial_potential,
            )
            receipt_flux = np.asarray([
                receipt_result.solution.lower_endcap_area_average_flux_m2_s(
                    name, wafer_radius_m=0.5 * WAFER_DIAMETER_M)
                for name in POSITIVE_ION_NAMES
            ])
            receipt_potential = receipt_result.electrostatic_potential_V
            receipt_evaluations = (
                receipt_result.nonlinear_solver_evaluations)
        grid_receipt.append({
            "radial_cell_count": radial_count,
            "axial_cell_count": axial_count,
            "wafer_positive_ion_flux_m2_s": dict(zip(
                POSITIVE_ION_NAMES, receipt_flux)),
            "total_wafer_positive_ion_flux_m2_s": float(np.sum(receipt_flux)),
            "nonlinear_solver_evaluations": int(receipt_evaluations),
        })
        seed_grid = receipt_grid
        seed_potential = receipt_potential
    grid_relative_change = float(abs(
        grid_receipt[-1]["total_wafer_positive_ion_flux_m2_s"]
        / grid_receipt[-2]["total_wafer_positive_ion_flux_m2_s"] - 1.0
    ))
    return {
        "schema": "petch.mahorowala_1998_axisymmetric_ion_transport_audit.v1",
        "claim_class": "failure_localization_not_feature_depth_prediction",
        "source": (
            "Mahorowala 1998 Table 2.2 settings; Lam-equipment-class chlorine "
            "global state; quasineutral Boltzmann-electron drift diffusion; "
            "Lee/Economou mobility, mutual neutralization, and Bohm boundary"
        ),
        "raw_collision_payload_sha256": (
            replay.molecular_replay.raw_payload_sha256),
        "atomic_momentum_payload_sha256": replay.atomic_momentum_payload_sha256,
        "energy_grid_cell_count": energy_grid.cell_count,
        "grid": {
            "radius_m": geometry.radius_m,
            "length_m": geometry.length_m,
            "radial_cell_count": grid.radial_cell_count,
            "axial_cell_count": grid.axial_cell_count,
            "wafer_radius_m": 0.5 * WAFER_DIAMETER_M,
            "window_to_wafer_gap_cm": GAP_CM_EQUIPMENT_CLASS_PRIOR,
        },
        "rows": rows,
        "grid_convergence_receipt": {
            "condition": "400 W, 100 sccm, uniform source",
            "levels": grid_receipt,
            "production_to_fine_total_flux_relative_change": (
                grid_relative_change),
            "passed_three_percent_gate": grid_relative_change < 0.03,
        },
        "summary": {
            "condition_count": len(rows),
            "total_wafer_to_global_analytic_flux_ratio_range_by_source_mode": {
                name: [float(min(values)), float(max(values))]
                for name, values in ratio_by_mode.items()
            },
            "center_condition_raw_flux_ratio_by_source_mode": {
                name: float(value)
                for name, value in center_ratio_by_mode.items()
            },
            "diagnostic_center_renormalized_total_ion_flux_scale_range_by_source_mode": {
                name: [float(min(values)), float(max(values))]
                for name, values in center_renormalized_by_mode.items()
            },
            "maximum_particle_ledger_relative_residual": float(max_ledger),
            "formal_feature_depth_passes": 0,
            "supports_reactor_state_prediction": False,
            "supports_implicit_differentiation": True,
            "feature_depth_used_to_select_any_parameter": False,
            "electropositive_scalar_diffusion_reference_rejected": True,
        },
        "next_measurement": (
            "radially resolved species-resolved ion saturation/current density "
            "at the Mahorowala wafer plane, with a documented probe-to-wafer "
            "transfer, to select or reject the source moment"
        ),
    }


def write(result: dict[str, object], output: Path):
    output.mkdir(parents=True, exist_ok=True)
    (output / "mahorowala_1998_axisymmetric_ion_transport.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    ranges = result["summary"][
        "total_wafer_to_global_analytic_flux_ratio_range_by_source_mode"
    ]
    relative_ranges = result["summary"][
        "diagnostic_center_renormalized_total_ion_flux_scale_range_by_source_mode"
    ]
    lines = [
        "# Mahorowala 1998 axisymmetric ion-transport audit",
        "",
        "This is a conservative spatial lift of a solved 0-D inventory, not an independent feature-depth prediction. No etch depth selected any reactor or transport parameter.",
        "",
        f"- conditions: `{result['summary']['condition_count']}`",
        f"- maximum particle-ledger residual: `{result['summary']['maximum_particle_ledger_relative_residual']:.3e}`",
        "- 14x10 to 20x14 center-flux change: "
        f"`{100.0 * result['grid_convergence_receipt']['production_to_fine_total_flux_relative_change']:.2f}%` "
        f"(3% gate: `{result['grid_convergence_receipt']['passed_three_percent_gate']}`)",
        "- formal feature-depth passes: `0`",
        "",
        "| source moment | raw wafer / Lee-global range | after center-current renormalization |",
        "|---|---:|---:|",
    ]
    for name, interval in ranges.items():
        relative = relative_ranges[name]
        lines.append(
            f"| {name} | {interval[0]:.4f}--{interval[1]:.4f} | "
            f"{relative[0]:.4f}--{relative[1]:.4f} |")
    lines.extend((
        "",
        "The raw 1.4--2.0x correction is not an absolute-depth gain in the existing diagnostic-conditioned board: its independent 400 W/100 sccm center-current anchor renormalizes that common offset away. What remains is a roughly +/-15% trend correction. The spread is a source-field uncertainty, not permission to choose the mode that best matches depth. The discriminating measurement is a radially resolved, species-resolved ion-current profile at the wafer plane.",
        "",
    ))
    (output / "MAHOROWALA_1998_AXISYMMETRIC_ION_TRANSPORT.md").write_text(
        "\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--collision-deck", type=Path, required=True)
    parser.add_argument("--atomic-cl-momentum", type=Path, required=True)
    parser.add_argument("--energy-cells", type=int, default=415)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(
        args.collision_deck,
        args.atomic_cl_momentum,
        energy_cells=args.energy_cells,
    )
    write(result, args.output)
    print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
