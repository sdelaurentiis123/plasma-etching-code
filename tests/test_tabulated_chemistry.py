from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_interaction_table import SurfaceInteractionDomainError
from petch.surface_kinetics import (
    EnergeticFlux, FaceResolvedEnergeticFlux, ParameterEvidence, SurfaceFluxes,
)
from petch.tabulated_chemistry import (
    TabulatedSiClArMechanism, TabulatedSiPhysicalSputterMechanism,
    TabulatedSiSurfaceState,
)


DATA = (
    Path(__file__).parents[1] / "data" / "surface_interactions" / "kounis_melas_2024")
SI_ATOM_DENSITY_M3 = 8.0 / (5.43e-10) ** 3


def _mechanism():
    table = load_kounis_melas_2024_tables(DATA).reactive_ion_etch
    return TabulatedSiClArMechanism(
        table, SI_ATOM_DENSITY_M3,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032 RIE in.lammps: diamond-Si lattice a=5.43 angstrom",
            "source_derived", supports_prediction_within_declared_domain=True))


def _sputter_mechanism():
    table = load_kounis_melas_2024_tables(DATA).sputtering
    return TabulatedSiPhysicalSputterMechanism(
        table, SI_ATOM_DENSITY_M3,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032 diamond-Si lattice a=5.43 angstrom",
            "source_derived", supports_prediction_within_declared_domain=True))


def _ions(energy=100.0, cosine=1.0, flux=2e21):
    return EnergeticFlux("Ar+", flux, np.array([energy]), np.array([cosine]), np.array([1.0]))


def test_tabulated_si_rie_replays_source_yield_and_propagates_md_uncertainty():
    mechanism = _mechanism(); state = mechanism.initial_state((4,))
    ratios = np.array([10.0, 50.0, 100.0, 200.0])
    fluxes = SurfaceFluxes({"Cl2": 2e21 * ratios}, (_ions(flux=np.full(4, 2e21)),))
    result = mechanism.advance(state, fluxes, 2.0)
    table = mechanism.table.evaluate({"cl2_to_ar_flux_ratio": ratios})

    expected_rate = 2e21 * table.values["reactive_etch_yield"]
    assert np.allclose(result.removed_atoms_m2, 2.0 * expected_rate)
    assert np.allclose(result.etch_velocity_m_s, expected_rate / SI_ATOM_DENSITY_M3)
    assert np.allclose(
        result.etch_velocity_standard_uncertainty_m_s,
        2e21 * table.standard_uncertainty["reactive_etch_yield"] / SI_ATOM_DENSITY_M3)
    assert result.table_fingerprint == mechanism.table.fingerprint
    assert result.validity.parameter_evidence_supports_prediction
    assert result.validity.nonpredictive_parameters == ()
    assert not result.material_exchange.product_routing_complete
    assert np.array_equal(
        result.material_exchange.unresolved_units_m2["Si_atom"], result.removed_atoms_m2)


def test_tabulated_si_rie_refuses_unreleased_energy_angle_ratio_and_species():
    mechanism = _mechanism(); state = TabulatedSiSurfaceState.bare()
    for ions, match in ((_ions(energy=90.0), "100 eV"), (_ions(cosine=0.9), "normal")):
        with pytest.raises(ValueError, match=match):
            mechanism.advance(state, SurfaceFluxes({"Cl2": 2e22}, (ions,)), 1.0)
    with pytest.raises(SurfaceInteractionDomainError, match="flux_ratio"):
        mechanism.advance(state, SurfaceFluxes({"Cl2": 2e21}, (_ions(),)), 1.0)
    with pytest.raises(ValueError, match="no Si-Cl2-Ar"):
        mechanism.advance(
            state, SurfaceFluxes({"F": 2e22, "Cl2": 2e22}, (_ions(),)), 1.0)


def test_tabulated_si_rie_zero_flux_is_exact_identity():
    mechanism = _mechanism(); state = TabulatedSiSurfaceState(3.0e18)
    result = mechanism.advance(state, SurfaceFluxes({}), 10.0)
    assert np.array_equal(result.state.removed_atoms_m2, state.removed_atoms_m2)
    assert np.array_equal(result.etch_velocity_m_s, 0.0)


def test_tabulated_si_rie_refuses_unparameterized_one_sided_exposures():
    mechanism = _mechanism(); state = TabulatedSiSurfaceState.bare()
    with pytest.raises(ValueError, match="simultaneous"):
        mechanism.advance(state, SurfaceFluxes({"Cl2": 2e22}), 1.0)
    with pytest.raises(ValueError, match="simultaneous"):
        mechanism.advance(state, SurfaceFluxes({}, (_ions(),)), 1.0)


def test_tabulated_si_physical_sputter_routes_every_removed_atom_to_named_product():
    mechanism = _sputter_mechanism(); state = mechanism.initial_state((3,))
    flux = np.array([1e20, 2e20, 3e20])
    result = mechanism.advance(
        state, SurfaceFluxes({}, (_ions(energy=100.0, flux=flux),)), 2.0)
    source = mechanism.table.evaluate({"ion_energy": 100.0})
    expected = 2.0 * flux * source.values["physical_sputter_yield"]

    assert np.allclose(result.removed_atoms_m2, expected)
    assert result.material_exchange.product_routing_complete
    assert np.array_equal(result.material_exchange.outgoing_units_m2["Si_atom"], expected)
    assert len(result.product_populations) == 1
    assert result.product_populations[0].name == "Si"
    assert np.array_equal(result.product_populations[0].integrated_particle_count_m2, expected)
    assert not result.product_populations[0].transport_ready
    assert "energy and angular" in " ".join(result.validity.known_model_form_omissions)


def test_tabulated_si_physical_sputter_integrates_face_resolved_energy_events():
    mechanism = _sputter_mechanism(); state = mechanism.initial_state((2,))
    events = FaceResolvedEnergeticFlux(
        "Ar+", 2, event_face=np.array([0, 0, 1]),
        event_flux_m2_s=np.array([1e20, 2e20, 4e20]),
        event_energy_eV=np.array([50.0, 200.0, 100.0]),
        event_cosine_incidence=np.ones(3))
    result = mechanism.advance(state, SurfaceFluxes({}, (events,)), 1.0)
    yields = mechanism.table.evaluate({"ion_energy": np.array([50.0, 200.0, 100.0])})
    expected = np.array([
        1e20 * yields.values["physical_sputter_yield"][0]
        + 2e20 * yields.values["physical_sputter_yield"][1],
        4e20 * yields.values["physical_sputter_yield"][2],
    ])

    assert np.allclose(result.removed_atoms_m2, expected)


def test_tabulated_si_physical_sputter_refuses_angle_and_energy_extrapolation():
    mechanism = _sputter_mechanism(); state = mechanism.initial_state()
    with pytest.raises(ValueError, match="normal-incidence"):
        mechanism.advance(state, SurfaceFluxes({}, (_ions(cosine=0.9),)), 1.0)
    with pytest.raises(SurfaceInteractionDomainError, match="ion_energy"):
        mechanism.advance(state, SurfaceFluxes({}, (_ions(energy=300.0),)), 1.0)
