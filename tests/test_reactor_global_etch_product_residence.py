import pytest

from petch.reactor_global import (
    DiagnosticConditionedEtchProductResidenceTransfer,
)


def _transfer(wall_reactivity):
    return DiagnosticConditionedEtchProductResidenceTransfer(
        reference_total_ion_flux_m2_s=1.0e20,
        reference_sicl2_to_total_ion_ratio=0.3,
        reference_gross_si_source_rate_m2_s=2.0e20,
        reference_exhaust_loss_frequency_s_inv=10.0,
        reactor_volume_m3=0.04,
        reactor_physical_area_m2=1.0,
        gas_temperature_K=500.0,
        wall_reactivity=wall_reactivity,
    )


def test_reference_condition_reproduces_independent_sicl2_ratio_exactly():
    transfer = _transfer(0.0)
    result = transfer.predict(
        gross_si_source_rate_m2_s=2.0e20,
        exhaust_loss_frequency_s_inv=10.0,
    )
    assert result.sicl2_flux_m2_s == pytest.approx(3.0e19)
    assert result.residence_time_s == pytest.approx(0.1)
    assert result.gross_si_source_scale == 1.0
    assert result.residence_time_scale == 1.0


def test_reflective_wall_has_exact_inverse_exhaust_scaling():
    transfer = _transfer(0.0)
    low_flow = transfer.predict(
        gross_si_source_rate_m2_s=2.0e20,
        exhaust_loss_frequency_s_inv=5.0,
    )
    high_flow = transfer.predict(
        gross_si_source_rate_m2_s=2.0e20,
        exhaust_loss_frequency_s_inv=20.0,
    )
    assert low_flow.sicl2_flux_m2_s / high_flow.sicl2_flux_m2_s == pytest.approx(4.0)


def test_reactive_wall_bounds_flow_response_without_erasing_source_scaling():
    reflective = _transfer(0.0)
    reactive = _transfer(1.0)
    assert reactive.wall_loss_frequency_s_inv > 0.0
    reflective_ratio = (
        reflective.predict(
            gross_si_source_rate_m2_s=2.0e20,
            exhaust_loss_frequency_s_inv=5.0,
        ).sicl2_flux_m2_s
        / reflective.predict(
            gross_si_source_rate_m2_s=2.0e20,
            exhaust_loss_frequency_s_inv=20.0,
        ).sicl2_flux_m2_s
    )
    reactive_ratio = (
        reactive.predict(
            gross_si_source_rate_m2_s=2.0e20,
            exhaust_loss_frequency_s_inv=5.0,
        ).sicl2_flux_m2_s
        / reactive.predict(
            gross_si_source_rate_m2_s=2.0e20,
            exhaust_loss_frequency_s_inv=20.0,
        ).sicl2_flux_m2_s
    )
    assert 1.0 < reactive_ratio < reflective_ratio
    doubled = reactive.predict(
        gross_si_source_rate_m2_s=4.0e20,
        exhaust_loss_frequency_s_inv=10.0,
    )
    assert doubled.sicl2_flux_m2_s == pytest.approx(6.0e19)
