from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest


sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
from deboer_2002_direct_validation import _integrated_curve, _key, predict_depths  # noqa: E402
from deboer_feature3d import floor_rate, ion_species  # noqa: E402


def _row(*, opening, depth, split):
    return SimpleNamespace(
        series_time_min=10.0, mask_opening_um=opening, etch_depth_um=depth,
        split=split,
        role=("open_rate_anchor" if split == "boundary_input" else "held_out_prediction"),
        digitization_uncertainty_y_um=0.5, measurement_uncertainty_um=None)


def test_constant_rate_profile_integral_reproduces_boundary_rate_at_every_width():
    curve = {
        "aspect_ratio": [0.0, 2.0, 4.0, 8.0, 12.0, 20.0],
        "normalized_rate": [1.0] * 6,
    }
    predictions = predict_depths([
        _row(opening=100.0, depth=20.0, split="boundary_input"),
        _row(opening=10.0, depth=20.0, split="held_out_transfer"),
    ], curve)

    assert np.allclose([item["predicted_depth_um"] for item in predictions], 20.0)


def test_profile_transit_integral_is_positive_and_strictly_increasing():
    _, rate, integral = _integrated_curve({
        "aspect_ratio": [0.0, 1.0, 2.0, 4.0],
        "normalized_rate": [1.0, 1.1, 0.8, 0.4],
    })
    assert np.all(rate > 0.0)
    assert np.all(np.diff(integral) > 0.0)


def test_rate_cache_separates_legacy_face_and_physical_area_observables():
    common = dict(s_f=0.08, dx_um=0.01, seed=0, aspect_ratio=2.0)
    assert _key(**common, floor_average="face") != _key(
        **common, floor_average="area")


def test_floor_rate_refuses_undeclared_surface_average_before_engine_execution():
    with pytest.raises(ValueError, match="floor_average"):
        floor_rate(0.0, None, floor_average="vertex")


def test_deboer_ion_population_declares_energy_and_unambiguous_angular_width():
    ion = ion_species(
        2.0e19, 1.0e-6, energy_eV=40.0, iad_sigma_deg=3.0, log2=4, seed=7)

    assert ion.provenance["representative_normal_energy_eV"] == pytest.approx(40.0)
    assert ion.provenance["iad_component_sigma_deg"] == pytest.approx(3.0)
    assert ion.provenance["approximate_polar_iad_fwhm_deg"] == pytest.approx(4.8075)
    assert ion.density_model.tangential_temperature_eV == pytest.approx(
        2.0 * 40.0 * np.deg2rad(3.0) ** 2)


@pytest.mark.parametrize("energy, sigma", [(1.0, 3.0), (40.0, 0.0)])
def test_deboer_ion_population_refuses_invalid_declared_inputs(energy, sigma):
    with pytest.raises(ValueError, match="invalid ion energy"):
        ion_species(2.0e19, 1.0e-6, energy_eV=energy, iad_sigma_deg=sigma, log2=4)
