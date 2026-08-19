import numpy as np
import pytest

from petch.tio2_square_pillar import (
    integrate_square_pillar_depth,
    integrate_square_pillar_depth_from_blanket_rate,
)


def _run(**overrides):
    arguments = dict(
        depth_nodes_nm=[0.0, 50.0, 90.0],
        floor_transmission=[1.0, 1.0, 1.0],
        mask_transmission=[1.0, 1.0, 1.0],
        film_thickness_nm=100.0,
        mask_thickness_nm=20.0,
        positive_ion_flux_m2_s=1.0e19,
        duration_s=100.0,
        mass_density_kg_m3=3250.0,
        formula_units_per_incident_ion=1.0,
        tio2_to_cr_selectivity=20.0,
        integration_step_s=1.0,
    )
    arguments.update(overrides)
    return integrate_square_pillar_depth(**arguments)


def test_constant_transmission_matches_atom_counted_blanket_rate():
    result = _run(mask_thickness_nm=100.0)
    expected = min(100.0, result.blanket_rate_nm_s * 100.0)
    assert np.isclose(result.mask_pinned_depth_nm, expected)
    assert np.isclose(result.controlled_depth_nm, expected)
    assert result.mask_survives_duration


def test_floor_attenuation_reduces_depth_without_changing_blanket_rate():
    open_result = _run(mask_thickness_nm=100.0)
    attenuated = _run(
        mask_thickness_nm=100.0,
        floor_transmission=[0.5, 0.5, 0.5],
    )
    assert attenuated.blanket_rate_nm_s == open_result.blanket_rate_nm_s
    assert np.isclose(
        attenuated.mask_pinned_depth_nm,
        0.5 * open_result.mask_pinned_depth_nm,
    )


def test_direct_blanket_rate_entry_matches_atom_counted_entry():
    atom_counted = _run()
    direct = integrate_square_pillar_depth_from_blanket_rate(
        depth_nodes_nm=[0.0, 50.0, 90.0],
        floor_transmission=[1.0, 1.0, 1.0],
        mask_transmission=[1.0, 1.0, 1.0],
        film_thickness_nm=100.0,
        mask_thickness_nm=20.0,
        blanket_tio2_rate_nm_s=atom_counted.blanket_rate_nm_s,
        duration_s=100.0,
        tio2_to_cr_selectivity=20.0,
    )
    assert direct == atom_counted


def test_mask_exhaustion_stops_controlled_profile_but_not_pinned_diagnostic():
    result = _run(mask_thickness_nm=0.1, tio2_to_cr_selectivity=1.0)
    assert not result.mask_survives_duration
    assert result.mask_exhaustion_time_s is not None
    assert result.controlled_depth_nm == result.depth_at_mask_exhaustion_nm
    assert result.controlled_depth_nm < result.mask_pinned_depth_nm
    assert result.residual_mask_nm == 0.0


def test_transmission_curve_refuses_unphysical_or_post_clear_nodes():
    with pytest.raises(ValueError, match="transmission curve"):
        _run(depth_nodes_nm=[0.0, 100.0], floor_transmission=[1.0, 1.0])
    with pytest.raises(ValueError, match="transmission curve"):
        _run(floor_transmission=[1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="transmission curve"):
        _run(mask_transmission=[1.0, 2.1, 1.0])
