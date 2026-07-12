from pathlib import Path

import numpy as np
import pytest

from petch.experimental_boundary import (
    Jeon2022BoundaryClosure, build_jeon_2022_boundary_state,
)
from petch.experimental_data import (
    load_jeon_2022_electron_bias_controls, load_jeon_2022_plasma_controls,
)


DATA = Path(__file__).parents[1] / "data" / "experimental" / "jeon_2022"


def _condition():
    plasma = load_jeon_2022_plasma_controls(DATA / "digitized_plasma_controls.csv")
    electron = load_jeon_2022_electron_bias_controls(
        DATA / "digitized_electron_bias_controls.csv")
    return (next(item for item in plasma if item.condition_family == "gas_fraction_cw"
                 and item.c4f8_fraction == 0.2),
            next(item for item in electron if item.condition_family == "gas_fraction_cw"
                 and item.c4f8_fraction == 0.2))


def _closure():
    return Jeon2022BoundaryClosure(
        ion_name="Ar+", ion_mass_amu=39.948,
        ion_normal_energy_eV=[800.0, 900.0], ion_normal_energy_weight=[0.25, 0.75],
        ion_tangential_temperature_eV=0.026,
        neutral_flux_fraction={"F": 0.3, "CF2": 0.7},
        neutral_mass_amu={"F": 18.998, "CF2": 50.005}, neutral_temperature_K=300.0,
        provenance={"source": "manufactured adapter gate"})


def test_jeon_boundary_preserves_measured_integrated_ratio_and_explicit_iedf():
    plasma, electron = _condition()
    boundary = build_jeon_2022_boundary_state(
        plasma, electron, _closure(), reference_plane_m=2e-6)

    ion = boundary.get("Ar+")
    neutral_flux = boundary.get("F").flux_m2_s + boundary.get("CF2").flux_m2_s
    assert np.isclose(neutral_flux / ion.flux_m2_s, plasma.neutral_to_ion_flux_ratio)
    assert np.isclose(ion.mean_energy_eV, 875.0 + 0.026)
    assert boundary.provenance["self_bias_is_not_iedf"] is True
    assert boundary.provenance["closure_supports_prediction"] is False


def test_jeon_boundary_refuses_mismatched_diagnostic_conditions():
    plasma, _ = _condition()
    electron = load_jeon_2022_electron_bias_controls(
        DATA / "digitized_electron_bias_controls.csv")[-1]
    with pytest.raises(ValueError, match="different experimental conditions"):
        build_jeon_2022_boundary_state(
            plasma, electron, _closure(), reference_plane_m=2e-6)


def test_jeon_boundary_closure_requires_complete_normalized_species_input():
    with pytest.raises(ValueError, match="invalid Jeon boundary closure"):
        Jeon2022BoundaryClosure(
            ion_name="Ar+", ion_mass_amu=39.948,
            ion_normal_energy_eV=[900.0], ion_normal_energy_weight=[1.0],
            ion_tangential_temperature_eV=0.026,
            neutral_flux_fraction={"CF2": 0.8}, neutral_mass_amu={"CF2": 50.005},
            neutral_temperature_K=300.0, provenance={})
