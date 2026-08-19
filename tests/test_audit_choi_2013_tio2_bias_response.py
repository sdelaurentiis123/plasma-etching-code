import pytest

from scripts.audit_choi_2013_tio2_bias_response import build


def test_choi_bias_board_preserves_source_reported_response():
    audit = build()

    assert audit["source_reported_endpoints"]["dc_bias_magnitude_V"] == [50.0, 250.0]
    assert audit["source_reported_endpoints"]["tio2_etch_rate_nm_min"] == [130.9, 197.2]
    assert audit["direct_response"]["etch_rate_ratio"] == pytest.approx(197.2 / 130.9)
    assert audit["direct_response"]["selectivity_decreases_over_reported_interval"] is True


def test_choi_board_requires_ion_assisted_chemistry_without_transferring_coefficients():
    audit = build()

    assert audit["model_discrimination"][
        "energy_independent_rate_normalized_removal_sufficient"
    ] is False
    assert audit["model_discrimination"]["pure_physical_sputtering_only_supported"] is False
    assert audit["model_discrimination"]["ion_assisted_surface_chemistry_required"] is True
    assert audit["two_point_sqrt_bias_decomposition"]["identified_as_surface_law"] is False
    assert audit["two_point_sqrt_bias_decomposition"][
        "endpoint_replay_max_abs_error_nm_min"
    ] < 1.0e-12
    assert audit["freddie_boundary"]["coefficient_transfer_allowed"] is False
    assert audit["freddie_boundary"]["changes_absolute_oxford_depth_forecast"] is False
