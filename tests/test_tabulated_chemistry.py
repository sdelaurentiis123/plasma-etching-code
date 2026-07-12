from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_interaction_table import SurfaceInteractionDomainError
from petch.surface_kinetics import EnergeticFlux, ParameterEvidence, SurfaceFluxes
from petch.tabulated_chemistry import TabulatedSiClArMechanism, TabulatedSiSurfaceState


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
