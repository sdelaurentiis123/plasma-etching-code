import numpy as np
import pytest

from petch.surface_kinetics import (
    EnergeticFlux,
    EnergeticYield,
    FaceResolvedEnergeticFlux,
    LowEnergyActivationYield,
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
        "complex_removal_reaction_order",
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


def test_reduced_mechanism_serializes_its_complete_physical_input_provenance():
    mechanism = _mechanism()
    provenance = dict(mechanism.provenance)

    assert provenance["model"] == "reduced-sio2-fluorocarbon-common-engine-v1"
    assert provenance["parameters"]["complex_formation_probability"] == {"CF2": 0.2}
    assert provenance["parameters"]["complex_sio2_yield"]["reference_yield"] == 0.2
    assert set(provenance["sources"]) == set(_evidence())
    assert provenance["known_omissions"] == list(mechanism.parameters.known_omissions)


def test_huang_2019_reduced_projection_preserves_published_reaction_table_and_one_scale():
    parameters = ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
        energetic_response_scale=1.25)
    mechanism = ReducedSiO2FluorocarbonMechanism(parameters)

    assert parameters.complex_formation_probability == {
        "CF": 0.4, "CF2": 0.3, "FC_complex_02": 0.2}
    assert parameters.polymer_deposition_probability_on_polymer["CF2"] == 0.0015
    assert parameters.polymer_deposition_probability_on_substrate == {}
    assert parameters.activated_polymer_deposition_probability_on_substrate[
        "CF2"] == 0.0015
    assert parameters.activated_polymer_deposition_probability_on_polymer[
        "CF2"] == 0.015
    assert parameters.complex_activation_yield.minimum_energy_eV == 5.0
    assert parameters.complex_activation_yield.maximum_energy_eV == 70.0
    assert parameters.polymer_activation_yield.maximum_energy_eV == 30.0
    assert parameters.energetic_polymer_deposition_yield.zero_energy_yield == 0.1
    assert parameters.energetic_polymer_deposition_yield.minimum_energy_eV == 5.0
    assert parameters.energetic_polymer_deposition_yield.maximum_energy_eV == 70.0
    assert "CF+" in parameters.energetic_polymer_deposition_species
    assert "Ar+" not in parameters.energetic_polymer_deposition_species
    assert parameters.bare_sio2_yield.reference_yield == 1.125
    assert parameters.complex_sio2_yield.reference_yield == 0.9375
    assert parameters.polymer_sputter_yield.reference_yield == 0.375
    assert parameters.bare_sio2_yield.threshold_energy_eV == 70.0
    assert parameters.complex_sio2_yield.reference_energy_eV == 140.0
    assert "common scale=1.25" in mechanism.provenance["sources"][
        "complex_sio2_yield"]["note"]
    assert "Jeong Figure 6 omits the atomic-F boundary flux" in parameters.known_omissions
    assert mechanism.provenance["sources"]["site_density_m2"]["source"].endswith(
        "aa6f40")

    with pytest.raises(ValueError, match="energetic_response_scale"):
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection(
            energetic_response_scale=0.0)


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


def test_low_energy_activation_yield_obeys_published_energy_window():
    law = LowEnergyActivationYield(0.3, 5.0, 30.0)

    assert law.evaluate(4.999, 1.0) == 0.0
    assert np.isclose(law.evaluate(5.0, 1.0), 0.3 * (1.0 - 5.0 / 30.0))
    assert np.isclose(law.evaluate(15.0, 1.0), 0.15)
    assert law.evaluate(30.0, 1.0) == 0.0
    assert law.evaluate(740.0, 1.0) == 0.0


def test_huang_low_energy_impacts_activate_complex_and_polymer_but_740ev_does_not():
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection())
    initial = SiO2SurfaceState(
        1.0, mechanism.parameters.polymer_monolayer_density_m2,
        activated_complex_fraction=0.0, activated_polymer_fraction=0.0)
    low = mechanism.advance(
        initial, SurfaceFluxes({}, (_ions(2.0e19, energy=15.0),)), 1.0)
    high = mechanism.advance(
        initial, SurfaceFluxes({}, (_ions(2.0e19, energy=740.0),)), 1.0)

    assert low.state.activated_complex_fraction > 0.0
    assert low.state.activated_polymer_fraction > 0.0
    assert high.state.activated_complex_fraction == 0.0
    assert high.state.activated_polymer_fraction == 0.0


def test_huang_activated_polymer_has_exact_tenfold_cf2_sticking_probability():
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection())
    monolayer = mechanism.parameters.polymer_monolayer_density_m2
    inactive = SiO2SurfaceState(1.0, 100.0 * monolayer)
    active = SiO2SurfaceState(
        1.0, 100.0 * monolayer, activated_polymer_fraction=1.0)

    inactive_probability = mechanism.neutral_reaction_probability(inactive)["CF2"]
    active_probability = mechanism.neutral_reaction_probability(active)["CF2"]

    assert np.isclose(active_probability / inactive_probability, 10.0, rtol=2e-15)


def test_huang_activated_complex_is_required_for_substrate_polymer_nucleation():
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection())
    inactive = SiO2SurfaceState(1.0, 0.0)
    active = SiO2SurfaceState(1.0, 0.0, activated_complex_fraction=1.0)

    inactive_probability = mechanism.neutral_reaction_probability(inactive)["CF2"]
    active_probability = mechanism.neutral_reaction_probability(active)["CF2"]

    assert inactive_probability == 0.0
    assert active_probability == 0.0015


def test_huang_low_energy_fluorocarbon_ion_deposits_polymer_on_complex_only():
    mechanism = ReducedSiO2FluorocarbonMechanism(
        ReducedSiO2FluorocarbonParameters.huang_kushner_2019_reduced_projection())
    complex_surface = SiO2SurfaceState(1.0, 0.0)
    bare_surface = SiO2SurfaceState.bare()
    fc_low = EnergeticFlux("CF+", 2.0e19, [15.0], [1.0], [1.0])
    ar_low = EnergeticFlux("Ar+", 2.0e19, [15.0], [1.0], [1.0])
    fc_high = EnergeticFlux("CF+", 2.0e19, [740.0], [1.0], [1.0])

    deposited = mechanism.advance(
        complex_surface, SurfaceFluxes({}, (fc_low,)), 1.0)
    wrong_species = mechanism.advance(
        complex_surface, SurfaceFluxes({}, (ar_low,)), 1.0)
    wrong_energy = mechanism.advance(
        complex_surface, SurfaceFluxes({}, (fc_high,)), 1.0)
    no_complex = mechanism.advance(
        bare_surface, SurfaceFluxes({}, (fc_low,)), 1.0)

    assert deposited.state.polymer_units_m2 > 0.0
    assert deposited.deposited_polymer_units_m2 > 0.0
    assert np.isclose(
        deposited.state.polymer_units_m2,
        deposited.deposited_polymer_units_m2
        - deposited.removed_polymer_units_m2,
        rtol=5e-13, atol=64.0)
    assert wrong_species.state.polymer_units_m2 == 0.0
    assert wrong_energy.state.polymer_units_m2 == 0.0
    assert no_complex.state.polymer_units_m2 == 0.0


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


def test_face_resolved_events_preserve_optional_impact_phase_space_immutably():
    position = np.array([[0.25, 0.5, 0.0], [0.75, 0.5, 0.0]])
    direction = np.array([[0.0, 0.0, -1.0], [0.6, 0.0, -0.8]])
    events = FaceResolvedEnergeticFlux(
        "electron", 2, event_face=[0, 1], event_flux_m2_s=[1e18, 2e18],
        event_energy_eV=[4.0, 8.0], event_cosine_incidence=[1.0, 0.8],
        event_position=position, event_incident_direction=direction)

    position[:] = -1.0; direction[:] = 0.0
    assert np.array_equal(events.event_position, [[0.25, 0.5, 0.0], [0.75, 0.5, 0.0]])
    assert np.array_equal(
        events.event_incident_direction, [[0.0, 0.0, -1.0], [0.6, 0.0, -0.8]])
    assert not events.event_position.flags.writeable
    assert not events.event_incident_direction.flags.writeable
    with pytest.raises(ValueError, match="unit vectors"):
        FaceResolvedEnergeticFlux(
            "electron", 1, [0], [1.0], [4.0], [1.0],
            event_incident_direction=[[0.0, 0.0, -2.0]])


def test_no_flux_is_an_exact_identity_and_zero_velocity():
    state = SiO2SurfaceState(
        [0.2, 0.8], [1e18, 2e18], [3e18, 4e18],
        [0.1, 0.4], [0.3, 0.7])
    result = _mechanism().advance(state, SurfaceFluxes({}), 10.0)
    assert np.array_equal(result.state.complex_fraction, state.complex_fraction)
    assert np.array_equal(result.state.polymer_units_m2, state.polymer_units_m2)
    assert np.array_equal(result.state.removed_formula_units_m2, state.removed_formula_units_m2)
    assert np.array_equal(
        result.state.activated_complex_fraction, state.activated_complex_fraction)
    assert np.array_equal(
        result.state.activated_polymer_fraction, state.activated_polymer_fraction)
    assert np.array_equal(result.etch_velocity_m_s, np.zeros(2))


def test_sio2_activation_state_enforces_dependent_bounds():
    with pytest.raises(ValueError, match="invalid SiO2 surface state"):
        SiO2SurfaceState(0.2, 1.0, activated_complex_fraction=0.3)
    with pytest.raises(ValueError, match="invalid SiO2 surface state"):
        SiO2SurfaceState(0.2, 1.0, activated_polymer_fraction=1.1)


def test_high_energy_turnover_removes_activated_sites_with_material_ledgers_closed():
    activation = LowEnergyActivationYield(0.3, 5.0, 30.0)
    mechanism = _mechanism(
        complex_activation_yield=activation,
        polymer_activation_yield=activation,
        activation_energetic_species=("Ar+",),
        activated_polymer_deposition_probability_on_substrate={"CF2": 0.2},
        activated_polymer_deposition_probability_on_polymer={"CF2": 0.5},
        evidence={
            **_evidence(),
            "activated_polymer_deposition_probability_on_substrate": ParameterEvidence(
                "synthetic analytic test", "manufactured"),
            "activated_polymer_deposition_probability_on_polymer": ParameterEvidence(
                "synthetic analytic test", "manufactured"),
            "complex_activation_yield": ParameterEvidence(
                "synthetic analytic test", "manufactured"),
            "polymer_activation_yield": ParameterEvidence(
                "synthetic analytic test", "manufactured"),
            "activation_energetic_species": ParameterEvidence(
                "synthetic analytic test", "manufactured"),
        })
    initial = SiO2SurfaceState(
        1.0, 4.0e18, activated_complex_fraction=1.0,
        activated_polymer_fraction=1.0)
    result = mechanism.advance(
        initial, SurfaceFluxes({}, (_ions(2.0e19, energy=100.0),)), 2.0)

    assert result.state.activated_complex_fraction < initial.activated_complex_fraction
    assert result.state.activated_polymer_fraction < initial.activated_polymer_fraction
    assert np.allclose(
        result.state.removed_formula_units_m2,
        result.removed_complex_units_m2 + result.removed_bare_formula_units_m2)
    assert np.allclose(
        result.state.polymer_units_m2 - initial.polymer_units_m2,
        -result.removed_polymer_units_m2, rtol=5e-13, atol=128.0)


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


def test_near_saturated_complex_uses_stable_complementary_rate_integrals():
    mechanism = _mechanism(
        site_density_m2=1.0e19,
        complex_formation_probability={"CF2": 1.0},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(1.0, 0.0)
    # Reproduces the scale separation that previously made duration-integral(theta)
    # lose the bare-site exposure in a saturated Jeong anchor face.
    formation_rate = 1.547e19
    desired_removal_rate = 6.516e13
    ion_yield = mechanism.parameters.complex_sio2_yield.evaluate(100.0, 1.0)
    fluxes = SurfaceFluxes(
        {"CF2": formation_rate},
        (_ions(desired_removal_rate / ion_yield),))

    result = mechanism.advance(initial, fluxes, 30.0)

    stored_site_change = (
        result.state.complex_fraction - initial.complex_fraction
    ) * mechanism.parameters.site_density_m2
    assert result.state.complex_fraction < initial.complex_fraction
    assert np.isclose(
        result.formed_complex_units_m2 - result.removed_complex_units_m2,
        stored_site_change, rtol=0.0,
        atol=8 * np.finfo(float).eps * mechanism.parameters.site_density_m2)


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
    assert not result.material_exchange.product_routing_complete
    assert np.allclose(
        result.material_exchange.removed_units_m2["SiO2_formula_unit"],
        result.removed_complex_units_m2 + result.removed_bare_formula_units_m2)
    assert np.all(result.material_exchange.residual_units_m2("SiO2_formula_unit") == 0.0)
    assert np.allclose(
        result.etch_velocity_m_s,
        (result.removed_complex_units_m2 + result.removed_bare_formula_units_m2)
        / mechanism.parameters.bulk_formula_density_m3 / 2.0)


def test_quadratic_complex_activation_is_exact_bounded_and_conservative():
    mechanism = _mechanism(
        complex_removal_reaction_order=2,
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.35, 0.0)
    fluxes = SurfaceFluxes({"CF2": 1.7e18}, (_ions(2.3e18),))

    full = mechanism.advance(initial, fluxes, 3.0)
    first_half = mechanism.advance(initial, fluxes, 1.5)
    two_halves = mechanism.advance(first_half.state, fluxes, 1.5)

    assert 0.0 < full.state.complex_fraction < 1.0
    assert np.isclose(
        full.state.complex_fraction, two_halves.state.complex_fraction,
        rtol=2e-14, atol=2e-14)
    assert np.isclose(
        full.formed_complex_units_m2,
        first_half.formed_complex_units_m2 + two_halves.formed_complex_units_m2,
        rtol=2e-14, atol=32.0)
    assert np.isclose(
        full.removed_complex_units_m2,
        first_half.removed_complex_units_m2 + two_halves.removed_complex_units_m2,
        rtol=2e-14, atol=32.0)
    site_change = (
        (full.state.complex_fraction - initial.complex_fraction)
        * mechanism.parameters.site_density_m2)
    assert np.isclose(
        site_change,
        full.formed_complex_units_m2 - full.removed_complex_units_m2,
        rtol=0.0, atol=32.0)
    assert full.formed_complex_units_m2 >= 0.0
    assert full.removed_complex_units_m2 >= 0.0


def test_quadratic_complex_removal_only_matches_rational_decay():
    mechanism = _mechanism(
        complex_removal_reaction_order=2,
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.6, 0.0)
    fluxes = SurfaceFluxes({}, (_ions(2.0e18),))
    duration = 4.0
    removal_hazard = (
        2.0e18 * mechanism.parameters.complex_sio2_yield.evaluate(100.0, 1.0)
        / mechanism.parameters.site_density_m2)

    result = mechanism.advance(initial, fluxes, duration)

    expected = initial.complex_fraction / (
        1.0 + removal_hazard * initial.complex_fraction * duration)
    assert np.isclose(result.state.complex_fraction, expected, rtol=2e-14, atol=0.0)
    assert np.isclose(
        result.removed_complex_units_m2,
        (initial.complex_fraction - expected) * mechanism.parameters.site_density_m2,
        rtol=2e-14, atol=32.0)


def test_complex_removal_reaction_order_is_declared_and_bounded():
    with pytest.raises(ValueError, match="reaction order"):
        _mechanism(complex_removal_reaction_order=3)

    mechanism = _mechanism(complex_removal_reaction_order=2, evidence={})
    validity = mechanism.validity(SurfaceFluxes({}))
    assert "complex_removal_reaction_order" in validity.nonpredictive_parameters


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


def test_zero_polymer_without_substrate_nucleation_is_an_exact_invariant():
    mechanism = _mechanism(
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={"CF2": 0.9})
    initial = SiO2SurfaceState(0.0, 0.0)
    fluxes = SurfaceFluxes({"CF2": 1.0e20}, (_ions(2.0e20),))

    result = mechanism.advance(initial, fluxes, 1000.0)

    assert result.state.polymer_units_m2 == 0.0
    assert result.deposited_polymer_units_m2 == 0.0
    assert result.removed_polymer_units_m2 == 0.0


def test_polymer_mixed_flux_preserves_inactive_face_as_bitwise_identity():
    mechanism = _mechanism(complex_formation_probability={})
    initial = SiO2SurfaceState([0.0, 0.0], [4.920299e18, 0.0])
    fluxes = SurfaceFluxes({"CF2": np.array([0.0, 1e18])})
    result = mechanism.advance(initial, fluxes, 10.0)

    assert result.state.polymer_units_m2[0] == initial.polymer_units_m2[0]
    assert result.deposited_polymer_units_m2[0] == 0.0
    assert result.removed_polymer_units_m2[0] == 0.0


def test_polymer_subrepresentational_throughput_is_exactly_ledgered_as_zero():
    mechanism = _mechanism(
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.0, 5.483335e16)

    result = mechanism.advance(initial, SurfaceFluxes({}, (_ions(10.0),)), 30.0)

    assert result.state.polymer_units_m2 == initial.polymer_units_m2
    assert result.deposited_polymer_units_m2 == 0.0
    assert result.removed_polymer_units_m2 == 0.0


def test_thick_polymer_slow_removal_keeps_the_exact_decay_on_the_physical_side():
    mechanism = _mechanism(
        complex_formation_probability={},
        polymer_deposition_probability_on_substrate={},
        polymer_deposition_probability_on_polymer={})
    initial = SiO2SurfaceState(0.0, 2.800839e20)
    removal_rate_m2_s = 2.835616e13
    ion_yield = mechanism.parameters.polymer_sputter_yield.evaluate(100.0, 1.0)
    fluxes = SurfaceFluxes(
        {}, (_ions(removal_rate_m2_s / ion_yield, energy=100.0),))

    result = mechanism.advance(initial, fluxes, 5.0)

    assert result.state.polymer_units_m2 < initial.polymer_units_m2
    assert 0.0 < result.removed_polymer_units_m2 <= removal_rate_m2_s * 5.0
    stored_change = result.state.polymer_units_m2 - initial.polymer_units_m2
    assert np.isclose(
        stored_change, -result.removed_polymer_units_m2,
        rtol=5e-13, atol=64.0 * np.spacing(float(initial.polymer_units_m2)))


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
    assert set(result.validity.nonpredictive_parameters) == (
        set(_evidence()) - {"complex_removal_reaction_order"})
    assert "polymer_crosslinking" in result.validity.known_model_form_omissions
