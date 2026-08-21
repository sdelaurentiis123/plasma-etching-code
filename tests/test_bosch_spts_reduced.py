from pathlib import Path

import numpy as np
import pytest

from petch.bosch_process_data import (
    BoschProcessTrace, C4F8_FLOW_CHANNEL, SF6_FLOW_CHANNEL,
    load_bosch_process_traces,
)
from petch.reactor_global.bosch_spts_reduced import (
    BoschSPTSDynamicWallLaw, BoschSPTSDynamicWallState,
    BoschSPTSReducedParameters, BoschSPTSWallConditioningLaw,
    DeterministicBoschSPTSReactorToWafer, _bosch_dynamic_wall_interval,
    advance_bosch_spts_dynamic_wall, conditioned_bosch_spts_parameters,
    solve_bosch_spts_reduced_reactor,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"


@pytest.fixture(scope="module")
def trace():
    return load_bosch_process_traces(
        DATA / "Process_data.nc", DATA / "Dictionary_process.nc")[0]


def test_measured_waveform_resolves_bosch_phases_and_closes_inventory(trace):
    solution = solve_bosch_spts_reduced_reactor(
        trace, BoschSPTSReducedParameters())
    sf6 = trace.channels[SF6_FLOW_CHANNEL][:-1] > 300.0
    c4f8 = trace.channels[C4F8_FLOW_CHANNEL][:-1] > 150.0

    assert solution.volume_average_density_m3.shape == (
        trace.elapsed_s.size - 1, 3)
    assert np.mean(solution.source_rate_m3_s[sf6, 0] > 0.0) > 0.95
    assert np.mean(solution.source_rate_m3_s[sf6, 0]) > 20.0 * np.mean(
        solution.source_rate_m3_s[c4f8, 0])
    assert np.mean(solution.source_rate_m3_s[c4f8, 1]) > 20.0 * np.mean(
        solution.source_rate_m3_s[sf6, 1])
    assert np.mean(solution.source_rate_m3_s[c4f8, 1] > 0.0) > 0.90
    assert np.mean(solution.source_rate_m3_s[sf6 | c4f8, 2] > 0.0) > 0.90
    assert solution.maximum_inventory_ledger_relative_residual < 1.0e-12
    residual = np.abs(
        solution.integrated_production_m3
        - solution.integrated_loss_m3 - solution.final_density_m3)
    scale = np.maximum(solution.integrated_production_m3, 1.0)
    assert np.max(residual / scale) < 1.0e-12


def test_zero_source_trace_has_zero_plasma_inventory(trace):
    channels = {name: np.array(value, copy=True) for name, value in trace.channels.items()}
    for name in (
            SF6_FLOW_CHANNEL, C4F8_FLOW_CHANNEL,
            "Stat3_Etch_MV_SourceRFLoadPower",
            "Stat3_Etch_MV_SourceRFReflectedPower"):
        channels[name][:] = 0.0
    dark = BoschProcessTrace(
        experiment_key="manufactured_dark_trace", source_group=trace.source_group,
        process_date=trace.process_date, wafer_number=trace.wafer_number,
        elapsed_s=trace.elapsed_s, channels=channels)

    solution = solve_bosch_spts_reduced_reactor(
        dark, BoschSPTSReducedParameters())

    assert np.all(solution.source_rate_m3_s == 0.0)
    assert np.all(solution.volume_average_density_m3 == 0.0)
    assert np.all(solution.final_density_m3 == 0.0)
    assert np.all(solution.integrated_loss_m3 == 0.0)


def test_measured_vpp_conditions_ion_energy_without_using_platen_power(trace):
    parameters = BoschSPTSReducedParameters()
    solution = solve_bosch_spts_reduced_reactor(trace, parameters)
    vpp = np.abs(trace.channels["Stat3_Etch_MV_PlatenRFPeakToPeak"][:-1])
    expected = (
        parameters.plasma_potential_per_electron_temperature
        * parameters.electron_temperature_eV
        + parameters.collisional_ion_energy_transmission
        * parameters.sheath_bias_fraction_of_vpp * vpp)

    assert np.array_equal(solution.ion_energy_eV, expected)
    assert np.ptp(solution.ion_energy_eV) > 100.0


def test_axisymmetric_lift_predicts_radial_flux_and_exact_linear_jvp(trace):
    model = DeterministicBoschSPTSReactorToWafer(BoschSPTSReducedParameters())
    solution = model.solve(trace)

    assert solution.maximum_axisymmetric_species_ledger_relative_residual < 1.0e-12
    assert solution.inventory_lift_condition_number < 10.0
    assert np.all(solution.wafer_area_average_flux_m2_s >= 0.0)
    assert np.any(np.ptp(solution.radial_flux_m2_s[:, 0, :], axis=1) > 0.0)
    assert model.density_to_radial_flux_jvp(
        solution.reactor.volume_average_density_m3
    ) == pytest.approx(solution.radial_flux_m2_s, rel=2.0e-15)


def test_absorbed_power_coupling_changes_flux_without_changing_surface_law(trace):
    low = DeterministicBoschSPTSReactorToWafer(
        BoschSPTSReducedParameters(absorbed_source_power_fraction=0.60)).solve(trace)
    high = DeterministicBoschSPTSReactorToWafer(
        BoschSPTSReducedParameters(absorbed_source_power_fraction=0.80)).solve(trace)
    sf6 = trace.channels[SF6_FLOW_CHANNEL][:-1] > 300.0

    assert np.mean(high.wafer_area_average_flux_m2_s[sf6, 0]) > np.mean(
        low.wafer_area_average_flux_m2_s[sf6, 0])
    assert np.mean(high.wafer_area_average_flux_m2_s[sf6, 2]) > np.mean(
        low.wafer_area_average_flux_m2_s[sf6, 2])


def test_species_resolved_dual_zone_source_changes_radial_shapes_conservatively(trace):
    annular = DeterministicBoschSPTSReactorToWafer(BoschSPTSReducedParameters())
    split = DeterministicBoschSPTSReactorToWafer(BoschSPTSReducedParameters(
        source_central_fraction=(0.0, 0.5, 1.0)))
    annular_solution = annular.solve(trace)
    split_solution = split.solve(trace)

    assert split_solution.maximum_axisymmetric_species_ledger_relative_residual < 1e-12
    assert not np.allclose(
        split_solution.radial_flux_m2_s[:, 1],
        annular_solution.radial_flux_m2_s[:, 1])
    assert not np.allclose(
        split_solution.radial_flux_m2_s[:, 2],
        annular_solution.radial_flux_m2_s[:, 2])
    assert np.array_equal(
        split_solution.reactor.volume_average_density_m3,
        annular_solution.reactor.volume_average_density_m3)


def test_shared_wall_conditioning_law_uses_only_declared_lot_features():
    law = BoschSPTSWallConditioningLaw(
        log_carbon_cycle_coefficient=0.4,
        silicon_precondition_coefficient=-0.2,
        silicon_oxide_precondition_coefficient=0.3,
    )

    assert law.multiplier("3C") == 1.0
    assert law.multiplier("9C") == pytest.approx(np.exp(0.4 * np.log(3.0)))
    assert law.multiplier("3C-Si") == pytest.approx(np.exp(-0.2))
    assert law.multiplier("3C-SiO2") == pytest.approx(np.exp(0.3))
    with pytest.raises(ValueError, match="undeclared"):
        law.multiplier("custom-lot")
    with pytest.raises(ValueError, match="invalid"):
        BoschSPTSWallConditioningLaw(log_carbon_cycle_coefficient=1.6)


def test_conditioning_changes_neutral_reactor_state_but_not_ions(trace):
    base_parameters = BoschSPTSReducedParameters()
    law = BoschSPTSWallConditioningLaw(
        silicon_precondition_coefficient=np.log(2.0))
    conditioned_parameters = conditioned_bosch_spts_parameters(
        base_parameters, law, "3C-Si")
    base = solve_bosch_spts_reduced_reactor(trace, base_parameters)
    conditioned = solve_bosch_spts_reduced_reactor(trace, conditioned_parameters)

    assert conditioned_parameters.neutral_wall_loss_multiplier == pytest.approx(2.0)
    assert np.all(
        np.mean(conditioned.volume_average_density_m3[:, :2], axis=0)
        < np.mean(base.volume_average_density_m3[:, :2], axis=0)
    )
    assert np.array_equal(
        conditioned.volume_average_density_m3[:, 2],
        base.volume_average_density_m3[:, 2],
    )
    assert np.array_equal(conditioned.ion_energy_eV, base.ion_energy_eV)


def test_dynamic_wall_exact_interval_preserves_bounds_and_one_sided_limits():
    mean, end = _bosch_dynamic_wall_interval(0.2, 0.7, 0.0)
    assert 0.2 < mean < end < 1.0
    expected_end = 1.0 + (0.2 - 1.0) * np.exp(-0.7)
    expected_mean = 1.0 + (0.2 - 1.0) * (-np.expm1(-0.7) / 0.7)
    assert end == pytest.approx(expected_end)
    assert mean == pytest.approx(expected_mean)

    mean, end = _bosch_dynamic_wall_interval(0.8, 0.0, 0.5)
    assert 0.0 < end < mean < 0.8
    assert end == pytest.approx(0.8 * np.exp(-0.5))
    assert mean == pytest.approx(0.8 * (-np.expm1(-0.5) / 0.5))

    assert _bosch_dynamic_wall_interval(0.37, 0.0, 0.0) == (0.37, 0.37)


def test_dynamic_wall_state_is_dose_driven_and_carried_between_wafers(trace):
    law = BoschSPTSDynamicWallLaw(
        deposition_rate_per_reference_wafer=0.20,
        cleaning_rate_per_reference_wafer=0.05,
        log_wall_loss_response=np.log(2.0),
    )
    state = BoschSPTSDynamicWallState()
    steps = []
    for _ in range(10):
        step = advance_bosch_spts_dynamic_wall(trace, law, state)
        steps.append(step)
        state = step.end_state

    assert all(step.normalized_c4f8_dose > 0.0 for step in steps)
    assert all(step.normalized_sf6_dose > 0.0 for step in steps)
    assert all(
        right.end_state.occupancy > left.end_state.occupancy
        for left, right in zip(steps, steps[1:])
    )
    assert all(
        right.combined_wall_loss_multiplier
        > left.combined_wall_loss_multiplier
        for left, right in zip(steps, steps[1:])
    )
    assert steps[1].start_state == steps[0].end_state
    assert steps[-1].end_state.occupancy < 1.0
    manifest = law.manifest()
    assert manifest["target_depth_used"] is False
    assert manifest["wafer_number_used"] is False
    assert manifest["per_lot_initial_state_fitted"] is False


def test_dynamic_wall_multiplier_changes_neutrals_but_not_ion_channel(trace):
    law = BoschSPTSDynamicWallLaw(
        deposition_rate_per_reference_wafer=0.8,
        cleaning_rate_per_reference_wafer=0.01,
        log_wall_loss_response=np.log(4.0),
    )
    step = advance_bosch_spts_dynamic_wall(
        trace, law, BoschSPTSDynamicWallState())
    assert 1.0 < step.combined_wall_loss_multiplier <= 4.0
    base_parameters = BoschSPTSReducedParameters()
    dynamic_parameters = BoschSPTSReducedParameters(
        neutral_wall_loss_multiplier=step.combined_wall_loss_multiplier)
    base = solve_bosch_spts_reduced_reactor(trace, base_parameters)
    dynamic = solve_bosch_spts_reduced_reactor(trace, dynamic_parameters)

    assert np.all(
        np.mean(dynamic.volume_average_density_m3[:, :2], axis=0)
        < np.mean(base.volume_average_density_m3[:, :2], axis=0)
    )
    assert np.array_equal(
        dynamic.volume_average_density_m3[:, 2],
        base.volume_average_density_m3[:, 2],
    )
    assert np.array_equal(dynamic.ion_energy_eV, base.ion_energy_eV)


def test_dynamic_wall_law_rejects_nonphysical_rates_and_states():
    with pytest.raises(ValueError, match="invalid Bosch dynamic wall law"):
        BoschSPTSDynamicWallLaw(
            deposition_rate_per_reference_wafer=0.0,
            cleaning_rate_per_reference_wafer=0.1,
            log_wall_loss_response=0.0,
        )
    with pytest.raises(ValueError, match="invalid Bosch dynamic wall law"):
        BoschSPTSDynamicWallLaw(
            deposition_rate_per_reference_wafer=0.1,
            cleaning_rate_per_reference_wafer=0.1,
            log_wall_loss_response=2.0,
        )
    with pytest.raises(ValueError, match="invalid Bosch dynamic wall state"):
        BoschSPTSDynamicWallState(1.01)
