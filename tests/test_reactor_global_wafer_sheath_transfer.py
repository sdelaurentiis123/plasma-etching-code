import numpy as np
import pytest

from petch.reactor_global import DiagnosticConditionedRFSheathTransfer
from petch.reactor_global.network import E_CHARGE_C


def _power(area, fluxes, gain_eV):
    return area * E_CHARGE_C * sum(fluxes.values()) * gain_eV


def test_static_waveform_recovers_exact_power_closed_monoenergetic_limit():
    area = 0.04
    fluxes = {"Cl+": 3.0e19, "Cl2+": 7.0e19}
    transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=area,
        plasma_potential_eV=20.0,
        frequency_hz=13.56e6,
        collapse_fraction=0.0,
        phase_count=32,
        steps_per_period=64,
        steps_per_transit=128,
        source="manufactured static sheath",
    )

    result = transfer.predict(
        positive_ion_flux_m2_s=fluxes,
        electron_temperature_eV=4.0,
        electron_density_m3=1.0e17,
        delivered_bias_power_W=_power(area, fluxes, 50.0),
    )

    assert np.isclose(result.bias_dc_component_v, 50.0, atol=0.03)
    assert np.isclose(result.realized_mean_bias_energy_gain_eV, 50.0, atol=1e-7)
    assert abs(result.power_closure_relative_residual) < 2e-9
    for distribution in result.distributions.values():
        assert np.isclose(distribution.mean_energy_eV, 72.0, atol=0.04)
        assert distribution.standard_deviation_eV < 0.02


def test_fully_modulated_sheath_is_deterministic_mass_resolved_and_power_closed():
    area = 0.04
    fluxes = {"Cl+": 2.5e19, "Cl2+": 7.5e19}
    transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=area,
        plasma_potential_eV=20.0,
        frequency_hz=13.56e6,
        collapse_fraction=1.0,
        phase_count=48,
        steps_per_period=96,
        steps_per_transit=128,
        source="manufactured fully modulated sheath",
    )
    arguments = dict(
        positive_ion_flux_m2_s=fluxes,
        electron_temperature_eV=3.6,
        electron_density_m3=1.0e17,
        delivered_bias_power_W=_power(area, fluxes, 80.0),
    )

    first = transfer.predict(**arguments)
    second = transfer.predict(**arguments)

    assert abs(first.power_closure_relative_residual) < 2e-8
    assert np.isclose(first.realized_mean_bias_energy_gain_eV, 80.0, atol=2e-6)
    assert first.sheath_rf_amplitude_v == pytest.approx(
        first.bias_dc_component_v
    )
    assert first.distributions["Cl+"].standard_deviation_eV > 1.0
    assert first.distributions["Cl2+"].standard_deviation_eV > 1.0
    assert not np.array_equal(
        first.distributions["Cl+"].energy_eV,
        first.distributions["Cl2+"].energy_eV,
    )
    np.testing.assert_array_equal(
        first.distributions["Cl+"].energy_eV,
        second.distributions["Cl+"].energy_eV,
    )
    assert not first.evidence_supports_prediction
    assert not first.supports_feature_depth


def test_zero_delivered_bias_power_recovers_static_plasma_sheath():
    transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=0.04,
        plasma_potential_eV=20.0,
        frequency_hz=13.56e6,
        collapse_fraction=1.0,
        phase_count=32,
        steps_per_period=64,
        steps_per_transit=128,
        source="manufactured zero-bias limit",
    )
    result = transfer.predict(
        positive_ion_flux_m2_s={"Cl+": 3.0e19, "Cl2+": 7.0e19},
        electron_temperature_eV=4.0,
        electron_density_m3=1.0e17,
        delivered_bias_power_W=0.0,
    )
    # Finite trajectory steps shift the numerically inverted static voltage by
    # only a few millivolts; energy and power closure are exact below.
    assert result.bias_dc_component_v == pytest.approx(0.0, abs=0.01)
    assert result.sheath_rf_amplitude_v == pytest.approx(0.0, abs=0.01)
    assert result.realized_mean_bias_energy_gain_eV == pytest.approx(
        0.0, abs=1.0e-7
    )
    assert result.reconstructed_bias_power_W == pytest.approx(0.0, abs=1.0e-8)


def test_sheath_transfer_rejects_an_unmapped_ion_species():
    transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45},
        electrode_area_m2=0.04,
        plasma_potential_eV=20.0,
        frequency_hz=13.56e6,
    )
    with pytest.raises(ValueError, match="invalid plasma state"):
        transfer.predict(
            positive_ion_flux_m2_s={"Cl+": 1.0e20, "Cl2+": 1.0e20},
            electron_temperature_eV=4.0,
            electron_density_m3=1.0e17,
            delivered_bias_power_W=80.0,
        )


def test_bias_voltage_projection_round_trips_power_closure():
    area = 0.04
    fluxes = {"Cl+": 2.5e19, "Cl2+": 7.5e19}
    transfer = DiagnosticConditionedRFSheathTransfer(
        ion_mass_amu={"Cl+": 35.45, "Cl2+": 70.90},
        electrode_area_m2=area,
        plasma_potential_eV=20.0,
        frequency_hz=13.56e6,
        collapse_fraction=1.0,
        phase_count=32,
        steps_per_period=64,
        steps_per_transit=128,
        source="manufactured voltage projection",
    )
    projected = transfer.project_from_bias_dc_component(
        positive_ion_flux_m2_s=fluxes,
        electron_temperature_eV=4.0,
        electron_density_m3=1.0e17,
        bias_dc_component_v=75.0,
    )
    replayed = transfer.predict(
        positive_ion_flux_m2_s=fluxes,
        electron_temperature_eV=4.0,
        electron_density_m3=1.0e17,
        delivered_bias_power_W=projected.delivered_bias_power_W,
    )
    assert replayed.bias_dc_component_v == pytest.approx(75.0, abs=2.0e-7)
    assert replayed.realized_mean_bias_energy_gain_eV == pytest.approx(
        projected.realized_mean_bias_energy_gain_eV, rel=2.0e-9
    )
    assert projected.power_closure_relative_residual == 0.0
