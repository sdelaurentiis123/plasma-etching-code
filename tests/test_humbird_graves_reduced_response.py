import numpy as np
import pytest

from petch.humbird_graves_reduced_response import (
    HumbirdGravesReducedResponse,
    HumbirdGravesResponseParameters,
)


def test_source_domain_and_zero_fluence_are_exact():
    response = HumbirdGravesReducedResponse()
    for method in (
            response.carbon_inventory_ml,
            response.fluorine_inventory_ml,
            response.si_yield_per_ion,
            response.cumulative_si_etch_ml):
        value = method(
            0.0,
            energy_eV=20.0,
            cf2_per_ion=9.0,
            atomic_f_per_ion=0.0,
        )
        assert value == 0.0
    with pytest.raises(ValueError, match="20 or 200 eV"):
        response.carbon_inventory_ml(
            100.0,
            energy_eV=100.0,
            cf2_per_ion=9.0,
            atomic_f_per_ion=0.0,
        )


def test_nrt_opportunities_and_layer_capacity_are_source_derived():
    response = HumbirdGravesReducedResponse()
    assert response.nrt_displacements_per_ion(20.0) == pytest.approx(
        0.18645280090855146, rel=2e-6)
    assert response.nrt_displacements_per_ion(200.0) == pytest.approx(
        3.563976623139938, rel=2e-6)
    assert response.carbon_capacity_ml(200.0) == pytest.approx(
        2.0 * response.carbon_capacity_ml(20.0))


def test_atomic_f_response_is_site_capped_not_linear_without_bound():
    response = HumbirdGravesReducedResponse()
    fluence = np.array([200.0, 800.0, 1400.0])
    zero = response.si_yield_per_ion(
        fluence, energy_eV=200.0,
        cf2_per_ion=9.0, atomic_f_per_ion=0.0)
    one = response.si_yield_per_ion(
        fluence, energy_eV=200.0,
        cf2_per_ion=8.0, atomic_f_per_ion=1.0)
    two = response.si_yield_per_ion(
        fluence, energy_eV=200.0,
        cf2_per_ion=7.0, atomic_f_per_ion=2.0)
    assert np.all(one > zero)
    assert np.all(two > one)
    assert np.all(np.diff(two) < 0.0)


def test_cumulative_removal_integrates_instantaneous_yield_and_stays_positive():
    response = HumbirdGravesReducedResponse()
    fluence = np.array([200.0, 400.0, 800.0, 1400.0])
    cumulative = response.cumulative_si_etch_ml(
        fluence,
        energy_eV=200.0,
        cf2_per_ion=7.0,
        atomic_f_per_ion=2.0,
    )
    assert np.all(np.diff(cumulative) > 0.0)
    assert cumulative[-1] == pytest.approx(79.15, rel=2e-3)


def test_all_kinetic_coefficients_are_labeled_as_regressed_or_derived():
    evidence = HumbirdGravesResponseParameters().evidence
    assert "regressed" in evidence["carbon_kinetics"]
    assert "regressed" in evidence["fluorine_kinetics"]
    assert "regressed" in evidence["silicon_yield"]
    assert "feature depth" in evidence["silicon_yield"]
