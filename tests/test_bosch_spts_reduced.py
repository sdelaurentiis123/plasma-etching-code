from pathlib import Path

import numpy as np
import pytest

from petch.bosch_process_data import (
    BoschProcessTrace, C4F8_FLOW_CHANNEL, SF6_FLOW_CHANNEL,
    load_bosch_process_traces,
)
from petch.reactor_global.bosch_spts_reduced import (
    BoschSPTSReducedParameters, BoschSPTSWallConditioningLaw,
    DeterministicBoschSPTSReactorToWafer, conditioned_bosch_spts_parameters,
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
