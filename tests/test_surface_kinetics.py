import numpy as np
import pytest

from petch.surface_kinetics import (
    EnergeticFlux,
    EnergeticYield,
    FaceResolvedEnergeticFlux,
    ParameterEvidence,
    ReducedSiO2FluorocarbonMechanism,
    ReducedSiO2FluorocarbonParameters,
    SiO2SurfaceState,
    SurfaceFluxes,
)


def _evidence():
    names = {
        "site_density_m2", "bulk_formula_density_m3", "polymer_monolayer_density_m2",
        "complex_formation_probability", "polymer_deposition_probability_on_substrate",
        "polymer_deposition_probability_on_polymer", "oxygen_polymer_etch_probability",
        "bare_sio2_yield", "complex_sio2_yield", "polymer_sputter_yield",
    }
    return {name: ParameterEvidence("synthetic analytic test", "manufactured") for name in names}


def _mechanism(**overrides):
    values = dict(
        site_density_m2=5.0e18,
        bulk_formula_density_m3=2.2e28,
        polymer_monolayer_density_m2=4.0e18,
        complex_formation_probability={"CF2": 0.2},
        polymer_deposition_probability_on_substrate={"CF2": 0.1},
        polymer_deposition_probability_on_polymer={"CF2": 0.05},
        oxygen_species="O",
        oxygen_polymer_etch_probability=0.1,
        bare_sio2_yield=EnergeticYield(0.05, 20.0, 100.0),
        complex_sio2_yield=EnergeticYield(0.2, 20.0, 100.0),
        polymer_sputter_yield=EnergeticYield(0.1, 10.0, 100.0),
        evidence=_evidence(),
    )
    values.update(overrides)
    return ReducedSiO2FluorocarbonMechanism(ReducedSiO2FluorocarbonParameters(**values))


def _ions(flux=1.0e19, energy=100.0, cosine=1.0):
    return EnergeticFlux("Ar+", flux, [energy], [cosine], [1.0])


def test_energetic_yield_reproduces_threshold_reference_and_angular_limits():
    law = EnergeticYield(
        0.2, 20.0, 100.0, energy_exponent=1.0, angular_model="chang_sawin_1997")
    assert law.evaluate(19.0, 1.0) == 0.0
    assert np.isclose(law.evaluate(100.0, 1.0), 0.2)
    assert np.isclose(law.evaluate(100.0, 0.5), 0.2)
    assert np.isclose(law.evaluate(100.0, 0.0), 0.0, atol=1e-15)

    sputter = EnergeticYield(
        0.2, 20.0, 100.0, angular_model="kress_1999", angular_parameter=9.3)
    assert np.isclose(sputter.evaluate(100.0, 1.0), 0.2)
    assert sputter.evaluate(100.0, 0.5) > 0.2
    assert sputter.evaluate(100.0, 0.0) == 0.0


def test_face_resolved_events_preserve_nonlinear_energy_angle_yield_without_averaging():
    events = FaceResolvedEnergeticFlux(
        "Ar+", 2, event_face=[0, 0, 1], event_flux_m2_s=[1e18, 2e18, 4e18],
        event_energy_eV=[20.0, 100.0, 60.0], event_cosine_incidence=[1.0, 0.5, 1.0])
    law = EnergeticYield(
        0.2, 20.0, 100.0, energy_exponent=2.0, angular_model="chang_sawin_1997")

    expected_event_yield = law.evaluate(
        events.event_energy_eV, events.event_cosine_incidence)
    expected = np.bincount(
        events.event_face, weights=events.event_flux_m2_s * expected_event_yield, minlength=2)
    assert np.array_equal(events.flux_m2_s, [3e18, 4e18])
    assert np.allclose(events.yield_rate_m2_s(law), expected)


def test_no_flux_is_an_exact_identity_and_zero_velocity():
    state = SiO2SurfaceState([0.2, 0.8], [1e18, 2e18], [3e18, 4e18])
    result = _mechanism().advance(state, SurfaceFluxes({}), 10.0)
    assert np.array_equal(result.state.complex_fraction, state.complex_fraction)
    assert np.array_equal(result.state.polymer_units_m2, state.polymer_units_m2)
    assert np.array_equal(result.state.removed_formula_units_m2, state.removed_formula_units_m2)
    assert np.array_equal(result.etch_velocity_m_s, np.zeros(2))


def test_neutral_transport_loss_is_sum_of_the_same_state_dependent_reaction_channels():
    mechanism = _mechanism(
        polymer_deposition_probability_on_substrate={"CF2": 0.3},
        polymer_deposition_probability_on_polymer={"CF2": 0.1},
        oxygen_polymer_etch_probability=0.4)
    state = SiO2SurfaceState(
        complex_fraction=[0.0, 0.5],
        polymer_units_m2=[0.0, np.log(2.0) * 4e18])
    probability = mechanism.neutral_reaction_probability(state)

    # Bare: CF2 complex formation (0.2) competes with substrate deposition (0.3).
    # Half-accessible: 0.2*.5*.5 complex + .3*.5+.1*.5 deposition.
    assert np.allclose(probability["CF2"], [0.5, 0.25])
    assert np.allclose(probability["O"], [0.0, 0.2])


def test_neutral_transport_refuses_competing_channel_probability_above_one():
    mechanism = _mechanism(
        complex_formation_probability={"CF2": 0.8},
        polymer_deposition_probability_on_substrate={"CF2": 0.4})
    with pytest.raises(ValueError, match="exceed one"):
        mechanism.neutral_reaction_probability(SiO2SurfaceState.bare())


def test_complex_formation_is_bounded_and_conserves_surface_sites():
    mechanism = _mechanism(
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    state = SiO2SurfaceState.bare((2,))
    fluxes = SurfaceFluxes({"CF2": np.array([1e18, 2e18])})
    result = mechanism.advance(state, fluxes, 3.0)

    expected = 1.0 - np.exp(-np.array([1e18, 2e18]) * 0.2 * 3.0 / 5e18)
    assert np.allclose(result.state.complex_fraction, expected)
    assert np.allclose(
        result.formed_complex_units_m2,
        result.state.complex_fraction * mechanism.parameters.site_density_m2)
    assert np.all(result.state.complex_fraction < 1.0)


def test_complex_mixed_flux_preserves_inactive_face_as_bitwise_identity():
    mechanism = _mechanism(
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState([0.7, 0.0], [0.0, 0.0])
    result = mechanism.advance(
        initial, SurfaceFluxes({"CF2": np.array([0.0, 1e18])}), 10.0)

    assert result.state.complex_fraction[0] == initial.complex_fraction[0]
    assert result.formed_complex_units_m2[0] == 0.0
    assert result.removed_complex_units_m2[0] == 0.0


def test_sub_ulp_complex_events_obey_float64_site_resolution_without_false_failure():
    mechanism = _mechanism(
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.695, 0.0)
    result = mechanism.advance(initial, SurfaceFluxes({"CF2": 1e3}), 10.0)

    assert abs(result.state.complex_fraction - initial.complex_fraction) <= np.spacing(0.695)
    assert result.formed_complex_units_m2 <= (
        8 * np.finfo(float).eps * mechanism.parameters.site_density_m2)


def test_complex_sputtering_conserves_sites_and_removed_formula_units():
    mechanism = _mechanism(
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    state = SiO2SurfaceState(np.ones(2), np.zeros(2))
    fluxes = SurfaceFluxes({}, (_ions(flux=np.array([1e18, 2e18])),))
    result = mechanism.advance(state, fluxes, 2.0)

    expected_fraction = np.exp(-np.array([1e18, 2e18]) * 0.2 * 2.0 / 5e18)
    assert np.allclose(result.state.complex_fraction, expected_fraction)
    assert np.allclose(
        result.removed_complex_units_m2,
        (1.0 - expected_fraction) * mechanism.parameters.site_density_m2)
    assert np.allclose(
        result.state.removed_formula_units_m2,
        result.removed_complex_units_m2 + result.removed_bare_formula_units_m2)
    assert np.allclose(
        result.etch_velocity_m_s,
        (result.removed_complex_units_m2 + result.removed_bare_formula_units_m2)
        / mechanism.parameters.bulk_formula_density_m3 / 2.0)


def test_polymer_inventory_has_exact_deposition_removal_balance_and_never_goes_negative():
    mechanism = _mechanism(complex_formation_probability={})
    initial = SiO2SurfaceState(0.0, 2.0e18)
    fluxes = SurfaceFluxes({"CF2": 3e18, "O": 4e18}, (_ions(2e18),))
    result = mechanism.advance(initial, fluxes, 4.0)

    assert result.state.polymer_units_m2 >= 0.0
    assert np.allclose(
        result.state.polymer_units_m2 - initial.polymer_units_m2,
        result.deposited_polymer_units_m2 - result.removed_polymer_units_m2,
        rtol=2e-11, atol=1e-8)


def test_polymer_removal_limit_reaches_exact_nonnegative_inventory_floor():
    mechanism = _mechanism(
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.0, 3.4e21)
    result = mechanism.advance(initial, SurfaceFluxes({}, (_ions(2e21),)), 1000.0)

    assert result.state.polymer_units_m2 == 0.0
    assert result.deposited_polymer_units_m2 == 0.0
    assert result.removed_polymer_units_m2 == initial.polymer_units_m2


def test_polymer_mixed_flux_preserves_inactive_face_as_bitwise_identity():
    mechanism = _mechanism(complex_formation_probability={})
    initial = SiO2SurfaceState([0.0, 0.0], [4.920299e18, 0.0])
    fluxes = SurfaceFluxes({"CF2": np.array([0.0, 1e18])})
    result = mechanism.advance(initial, fluxes, 10.0)

    assert result.state.polymer_units_m2[0] == initial.polymer_units_m2[0]
    assert result.deposited_polymer_units_m2[0] == 0.0
    assert result.removed_polymer_units_m2[0] == 0.0


def test_strang_coupling_converges_when_polymer_shields_substrate():
    mechanism = _mechanism()
    fluxes = SurfaceFluxes({"CF2": 5e18, "O": 1e18}, (_ions(2e18),))
    coarse = mechanism.advance(SiO2SurfaceState.bare(), fluxes, 10.0, max_step_s=1.0)
    medium = mechanism.advance(SiO2SurfaceState.bare(), fluxes, 10.0, max_step_s=0.5)
    fine = mechanism.advance(SiO2SurfaceState.bare(), fluxes, 10.0, max_step_s=0.25)
    coarse_error = abs(coarse.state.complex_fraction - fine.state.complex_fraction)
    medium_error = abs(medium.state.complex_fraction - fine.state.complex_fraction)
    assert medium_error < coarse_error
    assert medium_error < 2e-4


def test_unsupported_species_and_missing_evidence_refuse_strict_execution():
    mechanism = _mechanism(evidence={})
    fluxes = SurfaceFluxes({"C3F4": 1e18})
    validity = mechanism.validity(fluxes)
    assert not validity.within_declared_scope
    assert validity.unsupported_neutral_species == ("C3F4",)
    assert "missing parameter evidence" in validity.reasons[1]
    with pytest.raises(ValueError, match="outside declared scope"):
        mechanism.advance(SiO2SurfaceState.bare(), fluxes, 1.0)


def test_reduced_mechanism_reports_known_model_form_omissions():
    result = _mechanism().advance(SiO2SurfaceState.bare(), SurfaceFluxes({}), 0.0)
    assert result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert set(result.validity.nonpredictive_parameters) == set(_evidence())
    assert "polymer_crosslinking" in result.validity.known_model_form_omissions
