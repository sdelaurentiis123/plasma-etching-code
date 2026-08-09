from pathlib import Path

import numpy as np
import pytest

from petch.chang_sawin_chlorine_si import (
    BaloochCl2IonSiMechanism,
    BaloochCl2IonSiParameters,
    ChangSawinArClSiMechanism,
    ChangSawinArClSiParameters,
    ChangSawinClIonSiMechanism,
    ChangSawinClIonSiCl2SuppressionMechanism,
    ChangSawinClIonSiParameters,
    ChangSawinSiCl2SuppressionParameters,
    _chang_sawin_chemical_angular_factor,
)
from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


DATA = (
    Path(__file__).parents[1]
    / "data" / "surface_interactions" / "kounis_melas_2024")


def _ions(flux, *, cosine=1.0, energy=100.0):
    return EnergeticFlux(
        "Ar+", flux, np.array([energy]), np.array([cosine]), np.array([1.0]))


def _published_yield(ratio, *, y0=0.07, s=0.07, beta=0.83):
    theta = s * ratio / (s * ratio + 4.0 * beta)
    return theta, y0 * (1.0 - theta) + beta * theta


def test_molecular_card_reproduces_chang_table_3_4_site_balance_exactly():
    mechanism = ChangSawinArClSiMechanism()
    ratios = np.array([0.0, 10.0, 50.0, 100.0, 200.0, 1e8])
    ion_flux = 2.0e20
    result = mechanism.advance(
        mechanism.initial_state(ratios.shape),
        SurfaceFluxes(
            {"Cl2": ion_flux * ratios},
            (_ions(np.full(ratios.shape, ion_flux)),),
        ),
        1.0,
    )
    theta, expected_yield = _published_yield(ratios)

    assert np.allclose(result.chlorination_fraction, theta)
    assert np.allclose(
        result.etch_velocity_m_s,
        ion_flux * expected_yield / 5.0e28,
    )
    assert result.chlorination_fraction[-1] > 0.99999
    assert expected_yield[-1] == pytest.approx(0.83, rel=5e-7)


def test_normal_incidence_closure_conserves_si_and_chlorine_atoms():
    mechanism = ChangSawinArClSiMechanism()
    ion_flux = 2e20
    neutral_flux = 2e22
    result = mechanism.advance(
        mechanism.initial_state((1,)),
        SurfaceFluxes(
            {"Cl2": np.array([neutral_flux])},
            (_ions(np.array([ion_flux])),),
        ),
        3.0,
    )

    removed = (
        result.physical_removed_si_atoms_m2
        + result.chemical_removed_si_atoms_m2)
    assert np.array_equal(
        result.material_exchange.removed_units_m2["Si_atom"], removed)
    assert np.array_equal(
        result.material_exchange.outgoing_units_m2["Si_atom"], removed)
    assert np.allclose(
        result.consumed_chlorine_atoms_m2,
        4.0 * result.chemical_removed_si_atoms_m2)
    assert np.allclose(
        result.consumed_neutral_particles_m2,
        result.consumed_chlorine_atoms_m2 / 2.0)
    assert np.allclose(result.steady_site_balance_residual_cl_atoms_m2_s, 0.0)
    assert result.material_exchange.product_routing_complete
    assert {
        population.name for population in result.product_populations
    } == {"Si_physical", "SiCl4_chemical"}
    assert all(
        not population.transport_ready
        for population in result.product_populations)


def test_face_events_preserve_angle_measure_and_report_source_boundary():
    mechanism = ChangSawinArClSiMechanism(
        ChangSawinArClSiParameters.atomic_chlorine_100eV())
    events = FaceResolvedEnergeticFlux(
        "Ar+", 2,
        event_face=np.array([0, 0, 1]),
        event_flux_m2_s=np.array([1e20, 2e20, 3e20]),
        event_energy_eV=np.full(3, 100.0),
        event_cosine_incidence=np.cos(np.deg2rad([0.0, 60.0, 70.0])),
    )
    result = mechanism.advance(
        mechanism.initial_state((2,)),
        SurfaceFluxes({"Cl": np.array([3e22, 3e22])}, (events,)),
        1.0,
        strict=False,
    )

    factors = _chang_sawin_chemical_angular_factor(
        events.event_cosine_incidence)
    assert factors[0] == 1.0
    assert factors[1] == pytest.approx(1.0 / np.sqrt(2.0))
    assert factors[2] == pytest.approx(
        np.cos(np.deg2rad(70.0)) / np.cos(np.deg2rad(45.0)))
    assert result.etch_velocity_m_s[0] > result.etch_velocity_m_s[1]
    assert not result.validity.within_declared_scope
    assert any("measured angular range" in reason
               for reason in result.validity.reasons)


def test_molecular_card_refuses_energy_angle_and_hidden_species_extrapolation():
    mechanism = ChangSawinArClSiMechanism()
    state = mechanism.initial_state()
    with pytest.raises(ValueError, match="100 eV"):
        mechanism.advance(
            state, SurfaceFluxes({"Cl2": 1e22}, (_ions(1e20, energy=90.0),)),
            1.0)
    with pytest.raises(ValueError, match="molecular-Cl2 angular"):
        mechanism.advance(
            state, SurfaceFluxes(
                {"Cl2": 1e22}, (_ions(1e20, cosine=0.9),)),
            1.0)
    with pytest.raises(ValueError, match="no declared"):
        mechanism.advance(
            state, SurfaceFluxes(
                {"Cl2": 1e22, "F": 1e20}, (_ions(1e20),)),
            1.0)


def test_zero_duration_preserves_cumulative_state_exactly():
    mechanism = ChangSawinArClSiMechanism()
    state = mechanism.initial_state((2,))
    result = mechanism.advance(
        state,
        SurfaceFluxes(
            {"Cl2": np.array([1e22, 2e22])},
            (_ions(np.array([1e20, 1e20])),),
        ),
        0.0,
    )
    assert np.array_equal(
        result.state.removed_si_atoms_m2, state.removed_si_atoms_m2)
    assert np.array_equal(
        result.state.consumed_chlorine_atoms_m2,
        state.consumed_chlorine_atoms_m2)


def test_chang_card_is_an_independent_no_fit_diagnostic_against_deepmd():
    mechanism = ChangSawinArClSiMechanism()
    table = load_kounis_melas_2024_tables(DATA).reactive_ion_etch
    ratios = table.axes[0].values
    theta, prediction = _published_yield(ratios)
    observed = table.evaluate({
        "cl2_to_ar_flux_ratio": ratios,
    }).values["reactive_etch_yield"]

    # This is a retrospective cross-source diagnostic, not a preregistered
    # pass gate.  It freezes the actual disagreement so it cannot be described
    # later as exact agreement: the maximum error is 23.3%, RMSE is 0.0748.
    relative_error = np.abs(prediction - observed) / observed
    rmse = float(np.sqrt(np.mean((prediction - observed) ** 2)))
    assert np.max(relative_error) == pytest.approx(0.23229851, abs=1e-8)
    assert rmse == pytest.approx(0.07478523, abs=1e-8)
    assert np.all(np.diff(theta) > 0.0)


def test_parameter_provenance_does_not_promote_regression_to_first_principles():
    parameters = ChangSawinArClSiParameters.molecular_chlorine_100eV()
    mechanism = ChangSawinArClSiMechanism(parameters)

    assert parameters.evidence[
        "surface_chlorination_coefficient"].evidence_type == (
            "beam-yield regression")
    assert not parameters.evidence[
        "product_chlorine_atoms_per_si"].supports_prediction_within_declared_domain
    assert "not a first-principles" in mechanism.provenance["claim"]
    with pytest.raises(TypeError):
        parameters.evidence["surface_chlorination_coefficient"] = None


def test_cl_ion_cards_reproduce_chang_table_5_2_without_a_fit():
    parameters = ChangSawinClIonSiParameters.chang_thesis_table_5_2()

    chlorination, enhanced = parameters.coefficients(
        np.array([35.0, 60.0, 100.0]))

    np.testing.assert_allclose(chlorination, [0.18, 0.32, 0.45])
    np.testing.assert_allclose(enhanced, [1.14, 2.42, 3.61])
    with pytest.raises(ValueError, match="35--100 eV"):
        parameters.coefficients(120.0)
    extrapolated_s, extrapolated_beta = parameters.coefficients(
        120.0, allow_extrapolation=True)
    assert extrapolated_s > 0.45
    assert extrapolated_beta > 3.61


def test_cl_ion_normal_incidence_reproduces_eqs_5_1_to_5_3():
    mechanism = ChangSawinClIonSiMechanism()
    ion_flux = 2.0e20
    ratio = 200.0
    fluxes = SurfaceFluxes(
        {"Cl": ratio * ion_flux},
        (EnergeticFlux("Cl+", ion_flux, [35.0], [1.0], [1.0]),),
    )

    result = mechanism.advance(mechanism.initial_state(), fluxes, 2.0)

    expected_coverage = (
        0.18 * ratio + 1.0
    ) / (0.18 * ratio + 1.0 + 4.0 * 1.14)
    expected_rate = ion_flux * 1.14 * expected_coverage
    assert np.isclose(result.chlorination_fraction, expected_coverage)
    assert np.isclose(result.etch_velocity_m_s, expected_rate / 5.0e28)
    assert np.isclose(result.removed_si_atoms_m2, 2.0 * expected_rate)
    assert result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert result.validity.nonpredictive_parameters == (
        "product_chlorine_atoms_per_si",)
    assert result.steady_site_balance_residual_cl_atoms_m2_s < 1.0e6
    assert np.isclose(
        result.chlorine_atoms_supplied_by_neutrals_m2
        + result.chlorine_atoms_supplied_by_ions_m2,
        4.0 * result.removed_si_atoms_m2,
    )
    assert result.material_exchange.product_routing_complete
    assert result.product_populations[0].name == "SiCl4_chemical"


def test_cl_ion_quadrature_integrates_drives_before_site_balance():
    mechanism = ChangSawinClIonSiMechanism()
    events = FaceResolvedEnergeticFlux(
        "Cl+",
        2,
        np.array([0, 0, 1]),
        np.array([1.0e20, 3.0e20, 2.0e20]),
        np.array([35.0, 100.0, 60.0]),
        np.cos(np.deg2rad([0.0, 70.0, 60.0])),
    )
    neutral = np.array([4.0e22, 2.0e22])

    result = mechanism.advance(
        mechanism.initial_state((2,)),
        SurfaceFluxes({"Cl": neutral}, (events,)),
        1.0,
    )

    angular = np.minimum(
        events.event_cosine_incidence / np.cos(np.deg2rad(45.0)), 1.0)
    s, beta = mechanism.parameters.coefficients(events.event_energy_eV)
    total_ion = np.bincount(
        events.event_face, weights=events.event_flux_m2_s, minlength=2)
    ionic = np.bincount(
        events.event_face,
        weights=events.event_flux_m2_s * angular,
        minlength=2,
    )
    effective_s = np.bincount(
        events.event_face,
        weights=events.event_flux_m2_s * s,
        minlength=2,
    ) / total_ion
    capacity = np.bincount(
        events.event_face,
        weights=events.event_flux_m2_s * angular * beta,
        minlength=2,
    )
    supply = effective_s * neutral + ionic
    coverage = supply / (supply + 4.0 * capacity)
    np.testing.assert_allclose(result.chlorination_fraction, coverage)
    np.testing.assert_allclose(result.removed_si_atoms_m2, capacity * coverage)


def test_cl_ion_mechanism_fails_closed_on_unmeasured_boundary_axes():
    mechanism = ChangSawinClIonSiMechanism()
    low_ratio = SurfaceFluxes(
        {"Cl": 4.0e20},
        (EnergeticFlux("Cl+", 1.0e20, [60.0], [1.0], [1.0]),),
    )
    with pytest.raises(ValueError, match="ratio falls below"):
        mechanism.advance(mechanism.initial_state(), low_ratio, 1.0)

    molecular_ion = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux("Cl2+", 1.0e20, [60.0], [1.0], [1.0]),),
    )
    with pytest.raises(ValueError, match="no declared"):
        mechanism.advance(mechanism.initial_state(), molecular_ion, 1.0)

    outside_angle = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux(
            "Cl+", 1.0e20, [60.0],
            [np.cos(np.deg2rad(75.0))], [1.0]),),
    )
    with pytest.raises(ValueError, match="70 degree"):
        mechanism.advance(mechanism.initial_state(), outside_angle, 1.0)

    extrapolated_energy = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux("Cl+", 1.0e20, [120.0], [1.0], [1.0]),),
    )
    sensitivity = mechanism.advance(
        mechanism.initial_state(), extrapolated_energy, 1.0, strict=False)
    assert not sensitivity.validity.within_declared_scope
    assert sensitivity.removed_si_atoms_m2 > 0.0


def test_sicl2_suppression_reproduces_printed_equation_5_6_exactly():
    mechanism = ChangSawinClIonSiCl2SuppressionMechanism()
    ion_flux = 2.0e20
    cl_ratio = 200.0
    sicl2_ratio = np.array([0.0, 5.0, 10.0, 20.0])
    result = mechanism.advance(
        mechanism.initial_state(sicl2_ratio.shape),
        SurfaceFluxes(
            {
                "Cl": np.full(sicl2_ratio.shape, cl_ratio * ion_flux),
                "SiCl2": sicl2_ratio * ion_flux,
            },
            (
                FaceResolvedEnergeticFlux(
                    "Cl+",
                    sicl2_ratio.size,
                    np.arange(sicl2_ratio.size),
                    np.full(sicl2_ratio.size, ion_flux),
                    np.full(sicl2_ratio.size, 35.0),
                    np.ones(sicl2_ratio.size),
                ),
            ),
        ),
        1.0,
    )

    numerator = 0.18 * cl_ratio + 1.0 + 0.3 * sicl2_ratio
    denominator = (
        numerator + 4.0 * 1.14 + 3.0 * 10.0 * 0.3 * sicl2_ratio
    )
    coverage = numerator / denominator
    np.testing.assert_allclose(result.chlorination_fraction, coverage)
    np.testing.assert_allclose(
        result.removed_si_atoms_m2, ion_flux * 1.14 * coverage
    )
    np.testing.assert_allclose(
        result.sicl2_to_clplus_flux_ratio, sicl2_ratio
    )
    assert np.max(result.site_balance_residual_sites_m2_s) < 2.0e6
    assert result.material_exchange.product_routing_complete
    assert result.product_populations[0].name == (
        "SiCl4_substrate_with_SiCl2_suppression"
    )
    assert np.all(np.diff(result.etch_velocity_m_s) < 0.0)


def test_sicl2_suppression_evidence_and_scope_fail_closed():
    parameters = ChangSawinSiCl2SuppressionParameters.chang_thesis_equation_5_6()
    mechanism = ChangSawinClIonSiCl2SuppressionMechanism(
        suppression_parameters=parameters
    )
    assert parameters.sicl2_sticking_coefficient == 0.3
    assert parameters.chlorinated_sicl2_reaction_coefficient == 10.0
    assert not parameters.evidence[
        "single_coverage_closure"
    ].supports_prediction_within_declared_domain
    outside = SurfaceFluxes(
        {"Cl": 2.0e22, "SiCl2": 1.0e20},
        (EnergeticFlux("Cl+", 1.0e20, [60.0], [1.0], [1.0]),),
    )
    with pytest.raises(ValueError, match="only at 35 eV"):
        mechanism.advance(mechanism.initial_state(), outside, 1.0)
    sensitivity = mechanism.advance(
        mechanism.initial_state(), outside, 1.0, strict=False
    )
    assert not sensitivity.validity.within_declared_scope
    assert "single_coverage_closure" in (
        sensitivity.validity.nonpredictive_parameters
    )


def test_balooch_cl2plus_yield_replays_printed_slope_and_pixel_intercept():
    parameters = BaloochCl2IonSiParameters.chang_figure5_7()

    energy = np.array([parameters.threshold_energy_eV, 100.0, 400.0])
    expected = 0.22 * np.maximum(
        np.sqrt(energy) - np.sqrt(25.998846756576185), 0.0)

    np.testing.assert_allclose(parameters.yield_si_per_ion(energy), expected)
    assert parameters.yield_si_per_ion(parameters.threshold_energy_eV) == 0.0
    with pytest.raises(ValueError, match="26--625 eV"):
        parameters.yield_si_per_ion(700.0)


def test_balooch_cl2plus_absolute_removal_and_si_routing_close():
    mechanism = BaloochCl2IonSiMechanism()
    ion_flux = 3.0e20
    energy = 100.0
    fluxes = SurfaceFluxes(
        {"Cl": 3.0e22},
        (EnergeticFlux("Cl2+", ion_flux, [energy], [1.0], [1.0]),),
    )

    result = mechanism.advance(
        mechanism.initial_state(),
        fluxes,
        2.0,
        surface_chlorination_fraction=0.95,
    )

    expected_yield = 0.22 * (
        np.sqrt(energy) - np.sqrt(25.998846756576185))
    expected_removed = 2.0 * ion_flux * expected_yield
    assert np.isclose(result.removed_si_atoms_m2, expected_removed)
    assert np.isclose(result.incident_cl2plus_ions_m2, 2.0 * ion_flux)
    assert np.isclose(
        result.etch_velocity_m_s, ion_flux * expected_yield / 5.0e28)
    assert result.material_exchange.product_routing_complete
    assert result.product_populations[0].name == (
        "SiClx_from_Cl2plus_unresolved")
    assert result.validity.within_declared_scope


def test_balooch_cl2plus_fails_closed_on_surface_angle_and_energy_scope():
    mechanism = BaloochCl2IonSiMechanism()
    off_normal = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux(
            "Cl2+", 1.0e20, [100.0],
            [np.cos(np.deg2rad(1.0))], [1.0]),),
    )
    with pytest.raises(ValueError, match="normal-incidence"):
        mechanism.advance(
            mechanism.initial_state(), off_normal, 1.0,
            surface_chlorination_fraction=0.95)

    normal = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux("Cl2+", 1.0e20, [100.0], [1.0], [1.0]),),
    )
    with pytest.raises(ValueError, match="high-chlorination"):
        mechanism.advance(
            mechanism.initial_state(), normal, 1.0,
            surface_chlorination_fraction=0.80)

    outside = SurfaceFluxes(
        {"Cl": 1.0e22},
        (EnergeticFlux("Cl2+", 1.0e20, [700.0], [1.0], [1.0]),),
    )
    with pytest.raises(ValueError, match="26--625 eV"):
        mechanism.advance(
            mechanism.initial_state(), outside, 1.0,
            surface_chlorination_fraction=0.95)
