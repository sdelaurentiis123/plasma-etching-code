from pathlib import Path

import pytest

from petch.reactor_global.chlorine_vuv_spectrum import (
    Adf04Level,
    chlorine_direct_coronal_spectrum,
    load_open_adas_cl0_personal_research,
    map_observed_to_collision_levels,
    parse_adf04_levels_bytes,
    parse_adf04_type3_bytes,
)


def _collision_fixture() -> bytes:
    return b"""Cl+ 0        17         1      104591.00
    1  522563524555        ( 2)1(  1.5)            0.00
    2  522563524555        ( 2)1(  0.5)          632.31
    3  522563524545517     ( 4)1(  2.5)       114967.74
    -1
 1.00    3         1.00+02 1.00+04 1.00+06
    2    1 4.55-03 2.62-08 3.20-01 1.93+00 1.65+00
    3    1 2.05+05 5.44-09 1.42-01 4.09-02-5.33-04
    3    2 8.68-07 1.01-09 2.66-02 5.74-03 2.97-07
   -1   -1
"""


def _nist_fixture() -> bytes:
    return b"""Cl+ 0        17         1      104591.00
    1  521522563524555      ( 2)1(  1.5)            0.00
    2  521522563524555      ( 2)1(  0.5)          882.35
    3  521522563524545517   ( 4)1(  2.5)        71958.36
    -1
   -1   -1
"""


def test_fixed_width_parser_preserves_concatenated_signed_fields():
    dataset = parse_adf04_type3_bytes(_collision_fixture())
    transition = next(
        item for item in dataset.transitions
        if (item.upper_index, item.lower_index) == (3, 1)
    )
    assert transition.transition_probability_s_inv == 2.05e5
    assert transition.effective_collision_strengths == pytest.approx(
        (5.44e-9, 1.42e-1, 4.09e-2)
    )
    assert transition.high_energy_parameter == pytest.approx(-5.33e-4)


def test_observed_states_map_after_nist_core_prefix_is_removed():
    calculated = parse_adf04_type3_bytes(_collision_fixture())
    observed = parse_adf04_levels_bytes(_nist_fixture())
    assert map_observed_to_collision_levels(calculated.levels, observed) == {
        1: 1, 2: 2, 3: 3,
    }


def test_direct_coronal_rate_uses_observed_energy_and_is_fail_closed():
    calculated = parse_adf04_type3_bytes(_collision_fixture())
    observed = parse_adf04_levels_bytes(_nist_fixture())
    spectrum = chlorine_direct_coronal_spectrum(
        calculated, observed, electron_temperature_eV=2.5
    )
    line = next(item for item in spectrum.lines if item.lower_observed_index == 1)
    assert line.wavelength_nm == pytest.approx(138.9692, rel=1.0e-5)
    assert line.photon_rate_coefficient_cm3_s > 0.0
    assert not spectrum.prediction_supported
    quenched = chlorine_direct_coronal_spectrum(
        calculated,
        observed,
        electron_temperature_eV=2.5,
        nonradiative_loss_s_inv=2.05e5,
    )
    quenched_line = next(
        item for item in quenched.lines if item.lower_observed_index == 1
    )
    assert quenched_line.photon_rate_coefficient_cm3_s < line.photon_rate_coefficient_cm3_s


def test_loader_requires_license_acknowledgement_before_reading(tmp_path: Path):
    collision = tmp_path / "collision.dat"
    nist = tmp_path / "nist.dat"
    collision.write_bytes(_collision_fixture())
    nist.write_bytes(_nist_fixture())
    with pytest.raises(PermissionError, match="personal use"):
        load_open_adas_cl0_personal_research(collision, nist)
    with pytest.raises(ValueError, match="hash mismatch"):
        load_open_adas_cl0_personal_research(
            collision, nist, accept_restricted_personal_use=True
        )


def test_mapping_pairs_repeated_states_by_energy_rank():
    calculated = [
        Adf04Level(10, "abc", 2, "1", 0.5, 200.0),
        Adf04Level(11, "abc", 2, "1", 0.5, 100.0),
    ]
    observed = [
        Adf04Level(20, "521abc", 2, "1", 0.5, 20.0),
        Adf04Level(21, "521abc", 2, "1", 0.5, 10.0),
    ]
    assert map_observed_to_collision_levels(calculated, observed) == {
        21: 11, 20: 10,
    }
