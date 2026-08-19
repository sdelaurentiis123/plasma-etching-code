import numpy as np
import pytest

from petch.reactor_global.c4f6_ion_sources import (
    C4F6PositiveIonSourceModel,
    load_nist_c4f6_direct_ion_branches,
)


def test_nist_parent_partition_retains_heavy_products_and_closes_rate():
    branches = load_nist_c4f6_direct_ion_branches()
    rate = 3.7e20
    partition = branches.partition_event_rate_m3_s(rate)

    assert len(partition) == 15
    assert "C3F3+" in partition
    assert "C4F6+" in partition
    assert np.isclose(sum(partition.values()), rate, rtol=2e-15)
    assert not branches.supports_energy_resolved_branching
    assert not branches.supports_absolute_reactor_flux
    assert not branches.supports_krueger_boundary


def test_nist_direct_light_ion_ratios_replay_digitized_spectrum():
    branches = load_nist_c4f6_direct_ion_branches()
    partition = branches.partition_event_rate_m3_s(1.0)

    assert partition["CF2+"] / partition["CF+"] == pytest.approx(
        0.1044776119402985
    )
    assert partition["CF3+"] / partition["CF+"] == pytest.approx(
        0.31343283582089554
    )


def test_direct_plus_secondary_source_ledger_is_additive_and_fail_closed():
    model = C4F6PositiveIonSourceModel()
    result = model.evaluate(
        aggregate_parent_ionization_rate_m3_s=2.5e20,
        electron_density_m3=8.0e16,
        neutral_cfx_densities_m3={
            "CF": 1.0e18,
            "CF2": 4.0e18,
            "CF3": 2.0e18,
        },
        electron_temperature_eV=3.0,
        gas_temperature_K=350.0,
    )

    assert all(value > 0.0 for value in result.secondary_cfx_sources_m3_s.values())
    for name, secondary in result.secondary_cfx_sources_m3_s.items():
        assert result.combined_sources_m3_s[name] == pytest.approx(
            result.direct_parent_sources_m3_s[name] + secondary
        )
    assert not result.supports_steady_reactor_composition
    assert not result.supports_wafer_flux
    assert not result.supports_krueger_boundary
    assert "ion-neutral conversion and charge exchange" in (
        result.known_missing_operators
    )


def test_secondary_source_model_requires_complete_nonnegative_cfx_state():
    model = C4F6PositiveIonSourceModel()

    with pytest.raises(ValueError, match="invalid CFx"):
        model.evaluate(
            aggregate_parent_ionization_rate_m3_s=1.0,
            electron_density_m3=1.0e16,
            neutral_cfx_densities_m3={"CF": 1.0, "CF2": 1.0},
            electron_temperature_eV=3.0,
        )
