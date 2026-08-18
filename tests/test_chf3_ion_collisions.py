import math

import numpy as np
import pytest

from petch.reactor_global.chf3_ion_collisions import (
    ANGSTROM2_TO_M2,
    load_peko_2002_cf3_chf3_reactive_collision_model,
)
from scripts.digitize_peko_2002_chf3_dct import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    manifest_text,
)


def test_committed_digitization_replays_from_pixel_coordinates():
    payload = csv_text()
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)


def test_measured_dct_nodes_replay_and_interpolant_is_local_positive_c1():
    model = load_peko_2002_cf3_chf3_reactive_collision_model()
    replay = model.dct_sum_cross_section_m2(model.relative_energy_eV)
    np.testing.assert_allclose(
        replay / ANGSTROM2_TO_M2,
        model.dct_sum_cross_section_A2,
        rtol=0.0,
        atol=2.0e-15,
    )
    probe = np.geomspace(*model.relative_energy_support_eV, 1001)
    assert np.all(model.dct_sum_cross_section_m2(probe) > 0.0)
    derivative = model.dct_sum_cross_section_derivative_m2_per_eV(probe)
    assert np.all(np.isfinite(derivative))


def test_lab_to_relative_energy_conversion_and_support_fail_closed():
    model = load_peko_2002_cf3_chf3_reactive_collision_model()
    assert model.laboratory_to_relative_energy_factor == 70.0 / 139.0
    assert math.isclose(
        model.relative_energy_from_laboratory_eV(200.0),
        200.0 * 70.0 / 139.0,
    )
    lower_lab, upper_lab = model.laboratory_energy_support_eV
    assert 40.0 < lower_lab < 41.0
    assert 388.0 < upper_lab < 390.0
    with pytest.raises(ValueError, match="outside measured"):
        model.reactive_destruction_cross_section(20.0)
    with pytest.raises(ValueError, match="outside measured"):
        model.slab_sensitivity(
            laboratory_energy_eV=400.0,
            chf3_number_density_m3=1.0e20,
            path_length_m=1.0e-3,
            feed_fraction_used_as_density_proxy=False,
        )


def test_target_recipe_reactive_destruction_is_non_negligible_per_mm():
    model = load_peko_2002_cf3_chf3_reactive_collision_model()
    total_density = (
        30.0 * 0.133322368 / (1.380649e-23 * 350.0)
    )
    chf3_density_proxy = total_density * 55.0 / 61.0
    result = model.slab_sensitivity(
        laboratory_energy_eV=200.0,
        chf3_number_density_m3=chf3_density_proxy,
        path_length_m=1.0e-3,
        feed_fraction_used_as_density_proxy=True,
    )

    assert 0.18 < result.optical_depth_central < 0.20
    assert 0.16 < result.destruction_probability_central < 0.18
    assert (
        result.destruction_probability_lower
        < result.destruction_probability_central
        < result.destruction_probability_upper
    )
    assert result.feed_fraction_used_as_density_proxy is True
    assert result.supports_complete_molecular_transport is False
    assert result.supports_target_iead is False
    assert result.supports_absolute_depth_prediction is False


def test_cid_and_dct_uncertainties_are_not_collapsed_into_a_fake_fit():
    model = load_peko_2002_cf3_chf3_reactive_collision_model()
    energy = 100.0
    dct_A2 = model.dct_sum_cross_section_m2(energy) / ANGSTROM2_TO_M2
    band = model.reactive_destruction_cross_section(energy)

    assert math.isclose(
        band.central_m2 / ANGSTROM2_TO_M2,
        18.0 + dct_A2,
    )
    assert band.lower_m2 < band.central_m2 < band.upper_m2
    assert model.provenance["coefficient_selected_from_depth_target"] is None
