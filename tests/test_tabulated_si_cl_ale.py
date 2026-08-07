from pathlib import Path

import numpy as np
import pytest

from petch.interaction_data import load_kounis_melas_2024_tables
from petch.surface_interaction_table import SurfaceInteractionDomainError
from petch.surface_kinetics import (
    EnergeticFlux,
    ParameterEvidence,
    SurfaceFluxes,
)
from petch.tabulated_si_cl_ale import (
    TabulatedSiClAleProductMechanism,
    TabulatedSiClAleState,
)


DATA = (
    Path(__file__).parents[1]
    / "data" / "surface_interactions" / "kounis_melas_2024")
SI_ATOM_DENSITY_M3 = 8.0 / (5.43e-10) ** 3


def _mechanism():
    table = load_kounis_melas_2024_tables(DATA).ale_products
    return TabulatedSiClAleProductMechanism(
        table,
        SI_ATOM_DENSITY_M3,
        ParameterEvidence(
            "Kounis-Melas OSTI 2589032: diamond-Si lattice a=5.43 angstrom",
            "source_derived",
            supports_prediction_within_declared_domain=True,
        ),
    )


def _ions(*, energy=215.0, cosine=1.0, flux=1.0e19):
    return EnergeticFlux(
        "Ar+",
        flux,
        np.array([energy]),
        np.array([cosine]),
        np.array([1.0]),
    )


def test_ale_sequence_integrates_released_dose_windows_without_interpolation():
    mechanism = _mechanism()
    table = mechanism.table
    state = mechanism.initial_state(1.0e19)
    total_dose = mechanism.dose_bin_edges_m2[-1]
    flux = 1.0e19
    result = mechanism.advance(
        state, SurfaceFluxes({}, (_ions(flux=flux),)),
        total_dose / flux)

    width = np.diff(mechanism.dose_bin_edges_m2)
    expected = {
        name: np.sum(width * table.outputs[name])
        for name in table.outputs
    }
    for name, count in expected.items():
        assert result.product_counts_m2[name] == pytest.approx(count)
    expected_si = (
        expected["si_yield"]
        + expected["sicl_yield"]
        + expected["sicl2_yield"]
    )
    expected_cl = (
        expected["sicl_yield"]
        + 2.0 * expected["sicl2_yield"]
        + expected["cl_yield"]
    )
    assert result.removed_si_atoms_m2 == pytest.approx(expected_si)
    assert result.emitted_chlorine_atoms_m2 == pytest.approx(expected_cl)
    assert result.state.retained_chlorine_atoms_m2 == pytest.approx(
        1.0e19 - expected_cl)
    assert result.state.loaded_chlorine_atoms_m2 == pytest.approx(
        result.state.retained_chlorine_atoms_m2
        + result.state.emitted_chlorine_atoms_m2)
    assert result.material_exchange.product_routing_complete
    assert set(result.material_exchange.outgoing_units_m2) == {
        "Si_atom", "Cl_atom"}
    assert {item.name for item in result.product_populations} == {
        "Si", "SiCl", "SiCl2", "Cl"}
    sicl2 = next(
        item for item in result.product_populations
        if item.name == "SiCl2")
    assert sicl2.additional_source_inventories_per_particle == {
        "Cl_atom": 2.0}
    assert not sicl2.transport_ready


def test_ale_one_shot_and_binwise_advances_are_identical():
    mechanism = _mechanism()
    flux = 2.0e19
    initial = mechanism.initial_state(1.0e19)
    total_duration = mechanism.dose_bin_edges_m2[-1] / flux
    one_shot = mechanism.advance(
        initial, SurfaceFluxes({}, (_ions(flux=flux),)), total_duration)

    stepped = initial
    accumulated = {name: 0.0 for name in mechanism.table.outputs}
    for width in np.diff(mechanism.dose_bin_edges_m2):
        part = mechanism.advance(
            stepped, SurfaceFluxes({}, (_ions(flux=flux),)), width / flux)
        stepped = part.state
        for name in accumulated:
            accumulated[name] += float(part.product_counts_m2[name])

    assert stepped.ar_ion_dosage_m2 == pytest.approx(
        one_shot.state.ar_ion_dosage_m2)
    assert stepped.removed_si_atoms_m2 == pytest.approx(
        one_shot.state.removed_si_atoms_m2)
    assert stepped.retained_chlorine_atoms_m2 == pytest.approx(
        one_shot.state.retained_chlorine_atoms_m2)
    for name, count in accumulated.items():
        assert count == pytest.approx(one_shot.product_counts_m2[name])


def test_ale_refuses_wrong_condition_sequence_extension_and_hidden_chlorine():
    mechanism = _mechanism()
    state = mechanism.initial_state(1.0e19)
    with pytest.raises(ValueError, match="215 eV"):
        mechanism.advance(
            state, SurfaceFluxes({}, (_ions(energy=80.0),)), 1.0)
    with pytest.raises(ValueError, match="normal-incidence"):
        mechanism.advance(
            state, SurfaceFluxes({}, (_ions(cosine=0.9),)), 1.0)
    with pytest.raises(ValueError, match="no released"):
        mechanism.advance(
            state,
            SurfaceFluxes({"Cl2": 1.0e19}, (_ions(),)),
            1.0,
        )
    maximum = mechanism.dose_bin_edges_m2[-1]
    with pytest.raises(SurfaceInteractionDomainError, match="dose exceeds"):
        mechanism.advance(
            state, SurfaceFluxes({}, (_ions(flux=maximum),)), 1.01)
    with pytest.raises(ValueError, match="more Cl"):
        mechanism.advance(
            mechanism.initial_state(1.0e16),
            SurfaceFluxes({}, (_ions(flux=maximum),)),
            1.0,
        )


def test_ale_state_rejects_element_creation_and_reports_remap_physics():
    with pytest.raises(ValueError, match="does not close"):
        TabulatedSiClAleState(
            0.0,
            loaded_chlorine_atoms_m2=10.0,
            retained_chlorine_atoms_m2=8.0,
            emitted_chlorine_atoms_m2=1.0,
            removed_si_atoms_m2=0.0,
        )
    state = TabulatedSiClAleState.chlorinated(1.0e19, shape=(2,))
    assert state.surface_field_remap_modes()["ar_ion_dosage_m2"] == "intensive"
    assert state.surface_field_remap_modes()[
        "retained_chlorine_atoms_m2"] == "conservative"


def test_ale_zero_duration_is_exact_state_identity():
    mechanism = _mechanism()
    state = mechanism.initial_state(1.0e19, shape=(2,))
    result = mechanism.advance(state, SurfaceFluxes({}), 0.0)

    for name, before in state.conservative_surface_fields().items():
        assert np.array_equal(
            result.state.conservative_surface_fields()[name], before)
    assert np.array_equal(result.etch_velocity_m_s, np.zeros(2))
    assert all(
        np.array_equal(value, np.zeros(2))
        for value in result.product_counts_m2.values())
