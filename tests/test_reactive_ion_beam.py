from pathlib import Path

import numpy as np
import pytest

from petch.experimental_data import (
    KARAHASHI_2007_FIGURE4_SHA256,
    load_karahashi_2007_reactive_ion_yields,
)
from petch.mixed_layer_mechanism import (
    build_krueger_2024_mixed_layer_mechanisms,
)
from petch.reactive_ion_beam import Karahashi2007ReactiveIonYieldTable
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


DATA = (
    Path(__file__).parents[1] / "data" / "experimental" / "karahashi_2007"
    / "figure4_reactive_ion_yields.csv")


@pytest.fixture(scope="module")
def table():
    return Karahashi2007ReactiveIonYieldTable.from_observations(
        load_karahashi_2007_reactive_ion_yields(DATA),
        source_table_sha256=KARAHASHI_2007_FIGURE4_SHA256)


def test_table_reproduces_all_digitized_markers(table):
    observations = load_karahashi_2007_reactive_ion_yields(DATA)
    for row in observations:
        value = table.evaluate(row.species, row.energy_eV, 1.0)
        assert value == pytest.approx(row.yield_sio2_per_ion, abs=5e-13)


def test_table_interpolates_only_inside_each_species_support(table):
    value = table.evaluate("CF3+", np.array([1125.0]), np.array([1.0]))
    assert value == pytest.approx([(1.4703 + 1.5847) / 2.0])
    with pytest.raises(ValueError, match="measured support"):
        table.evaluate("CF+", 500.0, 1.0)
    with pytest.raises(ValueError, match="measured support"):
        table.evaluate("CF3+", 2500.0, 1.0)


def test_table_refuses_unmeasured_angle_and_unknown_species(table):
    with pytest.raises(ValueError, match="normal-incidence"):
        table.evaluate("CF3+", 1000.0, np.cos(np.deg2rad(1.0)))
    with pytest.raises(ValueError, match="no Karahashi"):
        table.evaluate("Ar+", 1000.0, 1.0)


def test_table_provenance_fingerprint_is_stable_and_nonempty(table):
    assert table.source_table_sha256 == KARAHASHI_2007_FIGURE4_SHA256
    assert len(table.fingerprint) == 64
    assert set(table.supported_species) == {"F+", "CF+", "CF2+", "CF3+"}


def test_table_uncertainty_broadcasts_and_closure_fingerprint_carries_tolerance(
        table):
    uncertainty = table.evaluate_uncertainty(
        "CF3+", 1000.0, np.ones(3))
    assert uncertainty.shape == (3,)
    assert uncertainty == pytest.approx(np.full(3, 0.011))
    changed_tolerance = Karahashi2007ReactiveIonYieldTable.from_observations(
        load_karahashi_2007_reactive_ion_yields(DATA),
        source_table_sha256=KARAHASHI_2007_FIGURE4_SHA256,
        cosine_tolerance=2e-5)
    assert changed_tolerance.fingerprint != table.fingerprint


def _beam(species, energy_eV, flux_m2_s=1e18):
    return SurfaceFluxes(
        neutral_flux_m2_s={},
        energetic_fluxes=(EnergeticFlux(
            name=species,
            flux_m2_s=flux_m2_s,
            energy_eV=np.array([energy_eV]),
            cosine_incidence=np.array([1.0]),
            weight=np.array([1.0]),
        ),))


def test_mixed_layer_end_to_end_reproduces_species_resolved_1000eV_ladder(table):
    oxide, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    measured = {}
    for species in ("F+", "CF+", "CF2+", "CF3+"):
        flux = 1e18
        duration = 1e-4
        result = oxide.advance(
            oxide.initial_state(()), _beam(species, 1000.0, flux),
            duration)
        measured[species] = (
            float(result.state.removed_formula_units_m2) / (flux * duration))
    assert [measured[name] for name in ("F+", "CF+", "CF2+", "CF3+")] == (
        sorted(measured.values()))
    expected = {"F+": 0.3232, "CF+": 0.6751, "CF2+": 1.1957, "CF3+": 1.4703}
    for species, target in expected.items():
        assert measured[species] == pytest.approx(target, rel=2e-12)


def test_reactive_beam_removal_closes_unresolved_material_exchange(table):
    oxide, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    flux = 1e18
    duration = 1e-4
    result = oxide.advance(
        oxide.initial_state(()), _beam("CF3+", 1000.0, flux), duration)
    expected = flux * duration * 1.4703
    exchange = result.material_exchange
    assert float(exchange.removed_units_m2["sio2_formula"]) == pytest.approx(
        expected, rel=2e-12)
    assert float(exchange.unresolved_units_m2["sio2_formula"]) == pytest.approx(
        expected, rel=2e-12)
    assert float(exchange.residual_units_m2("sio2_formula")) == 0.0
    assert not exchange.product_routing_complete


def test_mixed_layer_reactive_beam_closure_is_explicit_and_does_not_move_ar(table):
    measured, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    legacy, _ = build_krueger_2024_mixed_layer_mechanisms()
    state = measured.initial_state(())
    ar_flux = _beam("Ar+", 1000.0)
    measured_ar = measured.advance(state, ar_flux, 1e-4)
    legacy_ar = legacy.advance(legacy.initial_state(()), ar_flux, 1e-4)
    assert float(measured_ar.state.removed_formula_units_m2) == pytest.approx(
        float(legacy_ar.state.removed_formula_units_m2), rel=0.0, abs=0.0)
    assert measured.provenance["reactive_ion_yield_table"]["evidence_role"].startswith(
        "data-anchored closure")


def test_aggregate_krueger_ion_boundary_cannot_activate_species_table(table):
    measured, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    legacy, _ = build_krueger_2024_mixed_layer_mechanisms()
    aggregate = _beam("ions", 1000.0)
    measured_result = measured.advance(
        measured.initial_state(()), aggregate, 1e-4)
    legacy_result = legacy.advance(legacy.initial_state(()), aggregate, 1e-4)
    assert float(measured_result.state.removed_formula_units_m2) == pytest.approx(
        float(legacy_result.state.removed_formula_units_m2), rel=0.0, abs=0.0)
    assert any(
        "outside the Karahashi species table" in omission
        and "ions" in omission
        for omission in measured_result.validity.known_model_form_omissions)


def test_legacy_species_agnostic_path_is_declared_nonvalidating():
    legacy, _ = build_krueger_2024_mixed_layer_mechanisms()
    values = {}
    for species in ("F+", "CF+", "CF2+", "CF3+"):
        flux = 1e18
        duration = 1e-4
        beam = _beam(species, 1000.0, flux)
        result = legacy.advance(legacy.initial_state(()), beam, duration)
        values[species] = (
            float(result.state.removed_formula_units_m2) / (flux * duration))
    assert len(set(values.values())) == 1
    validity = legacy.validity(_beam("CF3+", 1000.0))
    assert any(
        "energetic ion identity is ignored" in omission
        for omission in validity.known_model_form_omissions)
    assert any(
        "stable parent-molecule/ion co-incidence is not represented" in omission
        for omission in validity.known_model_form_omissions)
    assert validity.parameter_evidence_supports_prediction is False


def test_mixed_layer_reactive_beam_closure_refuses_unmeasured_support(table):
    oxide, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    with pytest.raises(ValueError, match="measured support"):
        oxide.advance(
            oxide.initial_state(()), _beam("CF3+", 3406.0), 1e-4)


def test_zero_flux_reactive_population_does_not_trigger_support_refusal(table):
    oxide, _ = build_krueger_2024_mixed_layer_mechanisms(
        reactive_ion_yield_table=table)
    result = oxide.advance(
        oxide.initial_state(()), _beam("CF3+", 3406.0, 0.0), 1e-4)
    assert result.validity.within_declared_scope
    assert float(result.state.removed_formula_units_m2) == 0.0
