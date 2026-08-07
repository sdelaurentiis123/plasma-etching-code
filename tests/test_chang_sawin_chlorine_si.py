from pathlib import Path

import numpy as np
import pytest

from petch.chang_sawin_chlorine_si import (
    ChangSawinArClSiMechanism,
    ChangSawinArClSiParameters,
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
