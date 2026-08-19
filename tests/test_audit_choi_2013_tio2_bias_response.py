import pytest

from scripts.audit_choi_2013_tio2_bias_response import build


def test_choi_bias_board_preserves_source_reported_response():
    audit = build()
    bias = audit["source_reported_endpoints"]["dc_bias_sweep"]

    assert bias["dc_bias_magnitude_V"] == [50.0, 250.0]
    assert bias["tio2_etch_rate_nm_min"] == [130.9, 197.2]
    assert audit["bias_direct_response"]["etch_rate_ratio"] == pytest.approx(
        197.2 / 130.9
    )
    assert audit["bias_direct_response"][
        "selectivity_decreases_over_reported_interval"
    ] is True


def test_choi_board_preserves_orthogonal_source_reported_axes():
    audit = build()
    endpoint = audit["source_reported_endpoints"]

    assert endpoint["oxygen_sweep"]["o2_flow_sccm"] == [0.0, 3.0, 9.0]
    assert endpoint["oxygen_sweep"]["tio2_etch_rate_nm_min"] == [154.1, 179.4, 137.5]
    assert endpoint["source_power_sweep"]["tio2_etch_rate_nm_min"] == [136.0, 208.3]
    assert endpoint["pressure_sweep"]["tio2_etch_rate_nm_min"] == [187.7, 138.7]
    assert endpoint["afm_rms_roughness_angstrom"] == {
        "as_deposited": 36.5,
        "cf4_ar": 59.8,
        "o2_cf4_ar": 29.8,
    }
    response = audit["multiaxis_direct_response"]
    assert response["oxygen_rate_is_nonmonotonic"] is True
    assert response["oxygen_peak_flow_sccm"] == 3.0
    assert response["source_power_rate_ratio"] == pytest.approx(208.3 / 136.0)
    assert response["pressure_rate_ratio"] == pytest.approx(138.7 / 187.7)


def test_choi_board_requires_ion_assisted_chemistry_without_transferring_coefficients():
    audit = build()

    assert audit["model_discrimination"][
        "energy_independent_rate_normalized_removal_sufficient"
    ] is False
    assert audit["model_discrimination"]["pure_physical_sputtering_only_supported"] is False
    assert audit["model_discrimination"]["ion_assisted_surface_chemistry_required"] is True
    assert audit["model_discrimination"][
        "neutral_supply_and_competitive_oxygen_state_required"
    ] is True
    assert audit["model_discrimination"][
        "collisional_sheath_or_ion_delivery_pressure_response_required"
    ] is True
    assert audit["two_point_sqrt_bias_decomposition"]["identified_as_surface_law"] is False
    assert audit["two_point_sqrt_bias_decomposition"][
        "endpoint_replay_max_abs_error_nm_min"
    ] < 1.0e-12
    assert audit["freddie_boundary"]["coefficient_transfer_allowed"] is False
    assert audit["freddie_boundary"]["changes_absolute_oxford_depth_forecast"] is False
