import math

import pytest

from petch.depth_identifiability import (
    effective_yield_from_depth,
    ideal_gas_impingement_flux,
    krueger_2024_depth_identifiability,
    parent_fraction_for_flux_ratio,
)


def test_depth_normalization_distinguishes_wafer_and_floor_ions():
    target = effective_yield_from_depth(825.0, 60.0, 2.2e28, 1.2e20)
    assert target == pytest.approx(2.5208333333333335)
    assert target / 0.70 == pytest.approx(3.6011904761904763)


def test_parent_impingement_sensitivity_is_dimensionally_consistent():
    flux = ideal_gas_impingement_flux(0.541, 162.03, 300.0)
    assert flux == pytest.approx(6.46e21, rel=0.01)
    fraction = parent_fraction_for_flux_ratio(
        0.25, flux, 1.2e20, parent_floor_delivery=0.10,
        ion_floor_delivery=0.70)
    assert fraction == pytest.approx(0.0325, rel=0.02)


@pytest.mark.parametrize("value", [0.0, -1.0, math.inf, math.nan])
def test_parent_impingement_refuses_nonphysical_mass(value):
    with pytest.raises(ValueError):
        ideal_gas_impingement_flux(1.0, value, 300.0)


def test_audit_retracts_ceiling_without_promoting_analog_to_prediction():
    audit = krueger_2024_depth_identifiability(
        simulated_depth_nm=346.83264081620524,
        takada_400eV_peak_yield=1.1969,
        karahashi_pure_ion_peak_yield=1.8736)
    verdict = audit["verdict"]
    assert verdict["former_universal_1p5_ceiling_proof_valid"] is False
    assert verdict["published_inputs_identify_absolute_depth"] is False
    assert verdict["exact_825nm_prediction_authorized"] is False
    assert verdict["measured_missing_mechanism_is_plausible"] is True
    surface = audit["direct_surface_evidence"]
    assert surface[
        "takada_900eV_to_target_run_average_ratio_range"] == pytest.approx([
            2.4 / 2.5208333333333335,
            2.5 / 2.5208333333333335,
        ])
    assert "cannot be transplanted" in surface["interpretation"]
    c4f6 = audit["external_c4f6_boundary_evidence"]
    assert c4f6["neutral_mass_spectrum_contains_c4f6_parent_signal"] is True
    assert c4f6["reported_positive_ion_density_order"] == [
        "CF+", "C3F3+", "CF3+", "CF2+", "C3F5+",
    ]
    assert c4f6["absolute_wafer_flux_calibrated"] is False
    assert c4f6["transferable_to_krueger_boundary"] is False
