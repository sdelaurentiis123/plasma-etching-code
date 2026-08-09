import numpy as np
import pytest

from petch.reactor_global import (
    ChlorineWaferTransferEvidence,
    DiagnosticConditionedChlorineWaferTransfer,
    isotropic_thermal_particle_flux_m2_s,
)
from petch.reactor_global.network import E_CHARGE_C


def _evidence(*, measured=False):
    return ChlorineWaferTransferEvidence(
        reference_total_ion_flux="independent wafer-current diagnostic",
        reference_neutral_to_ion_ratio="independent radical diagnostic",
        electrode_area="measured powered electrode area",
        bias_power_coupling="delivered bias-power diagnostic",
        plasma_potential="independent plasma-potential diagnostic",
        equipment_transfer="same conditioned equipment",
        reference_facts_measured=measured,
        bias_coupling_measured=measured,
        plasma_potential_measured=measured,
    )


def _transfer(*, measured=False):
    return DiagnosticConditionedChlorineWaferTransfer(
        reference_model_positive_ion_flux_m2_s={
            "Cl+": 3.0e19,
            "Cl2+": 2.0e19,
        },
        reference_model_atomic_chlorine_density_m3=2.0e19,
        reference_gas_temperature_K=500.0,
        reference_wafer_total_ion_flux_m2_s=1.25e20,
        reference_wafer_neutral_to_ion_flux_ratio=100.0,
        electrode_area_m2=0.04,
        bias_power_to_ion_fraction=1.0,
        plasma_potential_eV=20.0,
        evidence=_evidence(measured=measured),
    )


def test_thermal_flux_is_exact_n_vbar_over_four():
    density = 2.0e19
    temperature = 500.0
    result = isotropic_thermal_particle_flux_m2_s(density, temperature)
    mass = 35.45 * 1.66053906660e-27
    expected = 0.25 * density * np.sqrt(
        8.0 * 1.380649e-23 * temperature / (np.pi * mass))
    assert result == pytest.approx(expected, rel=1e-15)


def test_reference_condition_closes_diagnostic_fluxes_and_ion_power():
    transfer = _transfer()
    result = transfer.predict(
        model_positive_ion_flux_m2_s={"Cl+": 3.0e19, "Cl2+": 2.0e19},
        model_atomic_chlorine_density_m3=2.0e19,
        gas_temperature_K=500.0,
        applied_bias_power_W=80.0,
    )

    assert result.total_positive_ion_flux_m2_s == pytest.approx(1.25e20)
    assert result.neutral_to_total_ion_flux_ratio == pytest.approx(100.0)
    expected_gain = 80.0 / (0.04 * E_CHARGE_C * 1.25e20)
    assert result.mean_sheath_energy_gain_eV == pytest.approx(expected_gain)
    assert result.mean_impact_energy_eV == pytest.approx(expected_gain + 20.0)
    assert result.reconstructed_ion_power_W == pytest.approx(80.0)
    assert abs(result.power_closure_relative_residual) < 1.0e-15
    assert not result.evidence_supports_prediction
    assert not result.supports_feature_depth


def test_frozen_diagnostic_transfer_preserves_species_and_state_trends():
    transfer = _transfer(measured=True)
    result = transfer.predict(
        model_positive_ion_flux_m2_s={"Cl+": 6.0e19, "Cl2+": 2.0e19},
        model_atomic_chlorine_density_m3=1.0e19,
        gas_temperature_K=500.0,
        applied_bias_power_W=80.0,
    )

    assert result.positive_ion_flux_m2_s["Cl+"] == pytest.approx(1.5e20)
    assert result.positive_ion_flux_m2_s["Cl2+"] == pytest.approx(5.0e19)
    assert result.atomic_chlorine_flux_m2_s == pytest.approx(6.25e21)
    assert result.neutral_to_total_ion_flux_ratio == pytest.approx(31.25)
    assert result.evidence_supports_prediction
    assert not result.supports_feature_depth


def test_transfer_rejects_species_drift_and_zero_reference_radicals():
    transfer = _transfer()
    with pytest.raises(ValueError, match="invalid state"):
        transfer.predict(
            model_positive_ion_flux_m2_s={"Cl+": 1.0e20},
            model_atomic_chlorine_density_m3=1.0e19,
            gas_temperature_K=500.0,
            applied_bias_power_W=80.0,
        )
    with pytest.raises(ValueError, match="reference model atomic"):
        DiagnosticConditionedChlorineWaferTransfer(
            reference_model_positive_ion_flux_m2_s={"Cl+": 1.0e20},
            reference_model_atomic_chlorine_density_m3=0.0,
            reference_gas_temperature_K=500.0,
            reference_wafer_total_ion_flux_m2_s=1.0e20,
            reference_wafer_neutral_to_ion_flux_ratio=100.0,
            electrode_area_m2=0.04,
            bias_power_to_ion_fraction=1.0,
            plasma_potential_eV=20.0,
            evidence=_evidence(),
        )
