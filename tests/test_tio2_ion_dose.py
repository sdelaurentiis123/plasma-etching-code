import math

import pytest

from petch.tio2_ion_dose import (
    build_clearance_gate,
    depth_nm_from_positive_ion_dose,
    minimum_feature_transmission_for_depth,
    required_formula_units_per_incident_ion,
    tio2_formula_unit_density_m3,
)


def test_tio2_formula_unit_density_is_atom_counted_from_mass_density():
    assert math.isclose(
        tio2_formula_unit_density_m3(3250.0),
        2.450599437808329e28,
        rel_tol=0.0,
        abs_tol=1.0e13,
    )


def test_required_yield_and_forward_depth_are_exact_inverses():
    flux = 2.2891234955809833e19
    required = required_formula_units_per_incident_ion(
        700.0, 3250.0, flux, 1200.0)
    depth = depth_nm_from_positive_ion_dose(
        required, 3250.0, flux, 1200.0)

    assert math.isclose(required, 0.6244819650320838)
    assert math.isclose(depth, 700.0)


def test_feature_transmission_is_not_silently_folded_into_surface_yield():
    flux = 2.2891234955809833e19
    blanket = required_formula_units_per_incident_ion(
        700.0, 4150.0, flux, 1200.0)
    feature = required_formula_units_per_incident_ion(
        700.0, 4150.0, flux, 1200.0, feature_transmission=0.5)

    assert math.isclose(feature, 2.0 * blanket)
    assert math.isclose(
        minimum_feature_transmission_for_depth(
            700.0, 4150.0, flux, 1200.0, 2.0),
        blanket / 2.0,
    )


def test_clearance_gate_closes_formula_unit_ledger():
    gate = build_clearance_gate(
        700.0,
        3940.0,
        2.2891234955809833e19,
        1200.0,
        feature_transmission=0.75,
    )
    removed_formula_units_m2 = (
        gate.required_formula_units_per_incident_ion
        * gate.delivered_positive_ion_dose_m2
    )
    required_formula_units_m2 = (
        gate.depth_nm * 1.0e-9 * gate.formula_unit_density_m3
    )

    assert math.isclose(
        removed_formula_units_m2,
        required_formula_units_m2,
        rel_tol=2.0e-16,
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"depth_nm": 0.0},
        {"mass_density_kg_m3": -1.0},
        {"positive_ion_flux_m2_s": float("nan")},
        {"duration_s": 0.0},
        {"feature_transmission": 1.01},
    ],
)
def test_clearance_gate_rejects_nonphysical_inputs(kwargs):
    values = {
        "depth_nm": 700.0,
        "mass_density_kg_m3": 3940.0,
        "positive_ion_flux_m2_s": 2.0e19,
        "duration_s": 1200.0,
        "feature_transmission": 1.0,
    }
    values.update(kwargs)
    with pytest.raises(ValueError):
        build_clearance_gate(**values)
