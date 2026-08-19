import numpy as np

from petch.iadf_two_component import kim_2025_reference_iadf
from petch.species_resolved_feature_boundary import (
    build_species_resolved_feature_boundary,
)


def _boundary():
    return build_species_resolved_feature_boundary(
        ion_flux_m2_s={"A+": 2.0e18, "B++": 3.0e18},
        ion_mass_amu={"A+": 20.0, "B++": 40.0},
        ion_charge_number={"A+": 1, "B++": 2},
        ion_energy_eV={
            "A+": np.asarray([100.0, 120.0]),
            "B++": np.asarray([200.0, 240.0]),
        },
        ion_energy_weight={
            "A+": np.asarray([0.25, 0.75]),
            "B++": np.asarray([0.25, 0.75]),
        },
        ion_iadf=kim_2025_reference_iadf(tail_fraction=0.0),
        neutral_flux_m2_s={"F": 4.0e20, "CF2": 5.0e20},
        neutral_mass_amu={"F": 18.9984, "CF2": 50.008},
        neutral_temperature_K=293.15,
        reference_plane_m=1.0e-6,
        ion_polar_order=8,
        ion_azimuthal_order=4,
    )


def test_species_identity_flux_charge_and_energy_are_preserved():
    boundary = _boundary()
    assert [item.name for item in boundary.species] == [
        "A+", "B++", "CF2", "F"
    ]
    assert boundary.get("A+").charge_number == 1
    assert boundary.get("B++").charge_number == 2
    assert boundary.get("B++").flux_m2_s == 3.0e18
    np.testing.assert_allclose(
        boundary.get("A+").mean_energy_eV, 115.0, rtol=2.0e-12
    )


def test_absolute_current_and_deterministic_replay_close():
    first = _boundary()
    second = _boundary()
    expected = 1.602176634e-19 * (2.0e18 + 2.0 * 3.0e18)
    assert np.isclose(first.current_density_A_m2, expected)
    for left, right in zip(first.species, second.species):
        np.testing.assert_array_equal(left.velocity_sqrt_eV, right.velocity_sqrt_eV)
        np.testing.assert_array_equal(left.weight, right.weight)
    assert first.provenance["monte_carlo"] is False
