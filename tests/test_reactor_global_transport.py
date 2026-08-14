import numpy as np
import pytest

from petch.reactor_global import (
    ArgonGlobalCondition,
    CylindricalReactor,
    LeeLiebermanArgonGlobalModel,
    LeeLiebermanArgonTransportProvider,
    PASCAL_PER_MTORR,
    STANDARD_PRESSURE_PA,
    argon_relative_temperature_eV,
    lee_lieberman_argon_ion_temperature_eV,
    nist_argon_self_diffusion_m2_s,
    phelps_argon_momentum_transfer_cross_section_m2,
    phelps_argon_momentum_transfer_rate_m3_s,
)


def _condition(*, pressure_mTorr=10.0, gas_temperature_K=600.0):
    return ArgonGlobalCondition(
        condition_id=f"argon-{pressure_mTorr:g}-mTorr",
        absorbed_power_W=1000.0,
        pressure_Pa=pressure_mTorr * PASCAL_PER_MTORR,
        gas_temperature_K=gas_temperature_K,
        geometry=CylindricalReactor(radius_m=0.1525, length_m=0.075),
        ion_wall_energy_factor_Te=5.0,
        ion_wall_energy_source="lee-lieberman-1994-global range 5--8 Te",
        ion_wall_energy_evidence="published_range_member",
    )


def test_phelps_momentum_transfer_law_uses_projectile_lab_energy():
    energies = np.array([0.01, 1.0, 100.0])
    expected = (
        1.15e-18
        * energies ** -0.1
        * (1.0 + 0.015 / energies) ** 0.6
    )
    assert np.allclose(
        phelps_argon_momentum_transfer_cross_section_m2(energies),
        expected,
        rtol=1.0e-15,
        atol=0.0,
    )
    assert np.isclose(
        phelps_argon_momentum_transfer_cross_section_m2(1.0),
        expected[1],
        rtol=1.0e-15,
        atol=0.0,
    )
    with pytest.raises(ValueError, match="laboratory"):
        phelps_argon_momentum_transfer_cross_section_m2(0.0)


def test_relative_temperature_and_momentum_rate_are_finite():
    relative = argon_relative_temperature_eV(0.5, 600.0)
    assert 0.25 < relative < 0.30
    rate = phelps_argon_momentum_transfer_rate_m3_s(0.5, 600.0)
    assert 1.0e-16 < rate < 1.0e-14
    assert rate == phelps_argon_momentum_transfer_rate_m3_s(0.5, 600.0)


def test_lee_lieberman_ion_temperature_transition_is_continuous():
    threshold = PASCAL_PER_MTORR
    assert lee_lieberman_argon_ion_temperature_eV(
        0.1 * threshold, 600.0) == 0.5
    assert lee_lieberman_argon_ion_temperature_eV(
        threshold, 600.0) == 0.5
    assert np.isclose(
        lee_lieberman_argon_ion_temperature_eV(
            threshold * (1.0 + 1.0e-10), 600.0),
        0.5,
        rtol=1.0e-10,
        atol=0.0,
    )
    assert (
        lee_lieberman_argon_ion_temperature_eV(100.0 * threshold, 600.0)
        < lee_lieberman_argon_ion_temperature_eV(10.0 * threshold, 600.0)
        < 0.5
    )


def test_nist_argon_self_diffusion_reproduces_reference_value_and_pressure_law():
    at_standard_pressure = nist_argon_self_diffusion_m2_s(
        298.15, STANDARD_PRESSURE_PA)
    assert np.isclose(at_standard_pressure, 0.182e-4, rtol=4.0e-3)
    at_half_pressure = nist_argon_self_diffusion_m2_s(
        298.15, 0.5 * STANDARD_PRESSURE_PA)
    assert np.isclose(
        at_half_pressure, 2.0 * at_standard_pressure,
        rtol=1.0e-15, atol=0.0)


def test_source_backed_provider_exposes_every_intermediate_and_no_fit_target():
    state = LeeLiebermanArgonTransportProvider().predict(
        _condition(), electron_temperature_eV=3.0)
    assert state.evidence_kind == "published_model"
    assert not state.supports_prediction
    assert 1.0e-4 < state.ion_mean_free_path_m < 0.1
    assert state.ambipolar_diffusion_m2_s > 0.0
    assert state.metastable_effective_diffusion_m2_s > 0.0
    assert state.provenance["nist_temperature_extrapolated"] is True
    assert state.provenance["coefficient_selection_target"] is None
    assert (
        state.metastable_effective_diffusion_m2_s
        < state.provenance["bulk_metastable_diffusion_m2_s"]
    )
    assert (
        state.metastable_effective_diffusion_m2_s
        < state.provenance["knudsen_metastable_diffusion_m2_s"]
    )


def test_transport_scales_with_pressure_and_closes_global_balances():
    provider = LeeLiebermanArgonTransportProvider()
    low = provider.predict(_condition(pressure_mTorr=1.0), 3.0)
    high = provider.predict(_condition(pressure_mTorr=100.0), 3.0)
    assert high.ion_mean_free_path_m < low.ion_mean_free_path_m
    assert high.ambipolar_diffusion_m2_s < low.ambipolar_diffusion_m2_s
    assert (
        high.metastable_effective_diffusion_m2_s
        < low.metastable_effective_diffusion_m2_s
    )

    solution = LeeLiebermanArgonGlobalModel(provider).solve(_condition())
    assert solution.maximum_normalized_residual <= 1.0e-8
    assert solution.transport.evidence_kind == "published_model"
    assert not solution.supports_prediction
