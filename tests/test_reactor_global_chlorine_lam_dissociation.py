import hashlib
from pathlib import Path

import pytest

from petch.reactor_global import (
    MALYSHEV_1998_CHLORINE_DISSOCIATION_CSV_SHA256,
    MalyshevMeasuredChlorineDissociationProvider,
    RateContext,
    build_hamilton_dissociation_chlorine_particle_network,
    malyshev_1998_eq7_transport_diagnostic,
    malyshev_1998_eq7_wall_return_inversion,
    malyshev_1998_eq11_relative_cl2_density_percent,
)


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTAL_CSV = (
    ROOT / "data" / "experimental" / "malyshev_1998_lam"
    / "figures7_8_chlorine_dissociation.csv"
)
PACKAGE_CSV = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "malyshev_1998_lam_chlorine_dissociation.csv"
)


def _provider():
    return MalyshevMeasuredChlorineDissociationProvider.from_package_data()


def test_eq11_observable_is_not_the_neutral_particle_fraction():
    molecule = 6.0e19
    atom = 8.0e19
    assert malyshev_1998_eq11_relative_cl2_density_percent(
        molecule, atom) == pytest.approx(60.0)
    assert 100.0 * molecule / (molecule + atom) == pytest.approx(
        42.857142857142854)

    with pytest.raises(ValueError, match="Eq. 11"):
        malyshev_1998_eq11_relative_cl2_density_percent(0.0, 0.0)


def _marker(*, gap_cm, pressure_mTorr, power_W):
    return next(
        marker for marker in _provider().markers
        if marker.window_to_wafer_gap_cm == gap_cm
        and marker.pressure_mTorr == pressure_mTorr
        and marker.tcp_source_power_W == power_W
    )


def test_packaged_lam_dissociation_board_is_the_audited_data_exactly():
    assert PACKAGE_CSV.read_bytes() == EXPERIMENTAL_CSV.read_bytes()
    assert hashlib.sha256(PACKAGE_CSV.read_bytes()).hexdigest() == (
        MALYSHEV_1998_CHLORINE_DISSOCIATION_CSV_SHA256
    )
    provider = _provider()
    assert len(provider.markers) == 38
    assert sum(
        marker.validation_role
        == "reactor_dissociation_validation_candidate"
        for marker in provider.markers
    ) == 37
    assert sum(
        marker.supports_eq7_inversion for marker in provider.markers
    ) == 33


def test_eq7_inversion_reproduces_marker_and_rate_ledger_exactly():
    marker = _marker(gap_cm=11.0, pressure_mTorr=2.0, power_W=499.936)
    state = malyshev_1998_eq7_wall_return_inversion(marker)
    context = RateContext(
        state.electron_temperature_state.electron_temperature.value)
    network = build_hamilton_dissociation_chlorine_particle_network()
    expected_hamilton = sum(
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
        if reaction.name.startswith("e_Cl2_dissociation_")
    )
    expected_attachment = next(
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
        if reaction.name == "e_Cl2_dissociative_attachment"
    )
    expected_frequency = (
        (expected_hamilton + expected_attachment)
        * state.electron_density_state.volume_average_electron_density.value
    )
    relative = marker.relative_cl2_density_percent / 100.0

    assert state.hamilton_neutral_dissociation_rate_m3_s == pytest.approx(
        expected_hamilton, rel=4.0e-15, abs=0.0)
    assert state.lee_dissociative_attachment_rate_m3_s == pytest.approx(
        expected_attachment, rel=4.0e-15, abs=0.0)
    assert state.electron_driven_cl2_destruction_frequency_s_inv == (
        pytest.approx(expected_frequency))
    assert state.required_wall_return_frequency_s_inv == pytest.approx(
        expected_frequency * relative / (2.0 * (1.0 - relative)))
    assert state.cl_to_cl2_number_density_ratio == pytest.approx(
        2.0 * (1.0 - relative) / relative)
    assert state.reproduced_relative_cl2_density_percent == pytest.approx(
        marker.relative_cl2_density_percent)


def test_reported_cl2_envelope_is_propagated_without_inventing_sigma():
    finite = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=11.0, pressure_mTorr=2.0, power_W=499.936))
    assert finite.reported_cl2_uncertainty_lower_frequency_s_inv < (
        finite.required_wall_return_frequency_s_inv)
    assert finite.reported_cl2_uncertainty_upper_frequency_s_inv > (
        finite.required_wall_return_frequency_s_inv)

    near_zero = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=11.0, pressure_mTorr=10.0, power_W=200.704))
    assert near_zero.reported_cl2_uncertainty_upper_frequency_s_inv is None
    assert near_zero.required_wall_return_frequency_s_inv > 1.0e4


def test_inversion_refuses_unphysical_marker_and_missing_electron_board():
    negative = _marker(
        gap_cm=11.0, pressure_mTorr=10.0, power_W=99.744)
    assert not negative.supports_eq7_inversion
    with pytest.raises(ValueError, match="cannot support"):
        malyshev_1998_eq7_wall_return_inversion(negative)

    unsupported_density = _marker(
        gap_cm=6.5, pressure_mTorr=1.0, power_W=699.744)
    with pytest.raises(ValueError, match="gap/pressure"):
        malyshev_1998_eq7_wall_return_inversion(unsupported_density)


def test_eq7_diagnostic_cannot_masquerade_as_wall_flux_or_depth_prediction():
    state = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=11.0, pressure_mTorr=1.0, power_W=499.936))

    assert not state.supports_prediction
    assert not state.supports_wall_probability_inference
    assert not state.supports_wafer_flux
    assert not state.supports_feature_depth
    assert not state.electron_temperature_state.supports_prediction
    assert not state.electron_density_state.supports_prediction


def test_eq7_transport_diagnostic_closes_exact_cylinder_at_declared_T():
    inversion = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=11.0, pressure_mTorr=2.0, power_W=499.936))
    state = malyshev_1998_eq7_transport_diagnostic(
        inversion,
        gas_temperature_K=333.0,
        gas_temperature_basis=(
            "source-reported initial 60 C wall/gas condition; powered gas "
            "temperature unmeasured"
        ),
    )
    relative_cl2 = (
        inversion.dissociation_marker.relative_cl2_density_percent / 100.0)
    assert state.initial_cl2_particle_fraction == pytest.approx(0.95)
    assert state.particle_pressure_multiplier == pytest.approx(
        1.0 + 0.95 * (1.0 - relative_cl2))
    assert state.diffusivity_state.provenance[
        "source_correction_factor"] == 1.25
    assert state.diffusivity_state.provenance[
        "feature_depth_target"] is None
    assert state.target_is_transport_attainable
    assert state.effective_wall_recombination_probability == pytest.approx(
        0.06544819172183121)
    assert state.matched_wall_state.exact_loss_frequency_s_inv == (
        pytest.approx(inversion.required_wall_return_frequency_s_inv))
    assert state.matched_wall_state.numerical_closure_passes
    assert not state.supports_prediction
    assert not state.supports_local_wall_probability_prediction
    assert not state.supports_wafer_flux
    assert not state.supports_feature_depth


def test_eq7_transport_diagnostic_fails_closed_above_absorbing_limit():
    inversion = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=11.0, pressure_mTorr=10.0, power_W=200.704))
    state = malyshev_1998_eq7_transport_diagnostic(
        inversion,
        gas_temperature_K=333.0,
        gas_temperature_basis="source initial temperature sensitivity",
    )
    assert not state.target_is_transport_attainable
    assert state.effective_wall_recombination_probability is None
    assert state.matched_wall_state is None
    assert state.status == "target_exceeds_absorbing_wall_limit"
    assert state.required_wall_return_frequency_s_inv > (
        state.absorbing_wall_state.exact_loss_frequency_s_inv)


def test_eq7_transport_diagnostic_requires_temperature_evidence_label():
    inversion = malyshev_1998_eq7_wall_return_inversion(
        _marker(gap_cm=6.5, pressure_mTorr=2.0, power_W=499.297))
    with pytest.raises(ValueError, match="evidence basis"):
        malyshev_1998_eq7_transport_diagnostic(
            inversion,
            gas_temperature_K=333.0,
            gas_temperature_basis="",
        )
