import numpy as np
import pytest

from petch.chlorine_poly_si import (
    HwangGiapisClSiMechanism, HwangGiapisClSiParameters, HwangGiapisClSiYield,
)
from petch.surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes


def test_hwang_giapis_yield_reproduces_eq_4_1_threshold_and_angle_cap():
    law = HwangGiapisClSiYield()
    energy = np.array([9.0, 10.0, 177.0, 177.0, 177.0])
    cosine = np.cos(np.deg2rad([0.0, 0.0, 0.0, 45.0, 80.0]))

    result = law.evaluate(energy, cosine)

    assert np.allclose(result[:2], 0.0)
    assert np.isclose(result[2], 0.1 * (np.sqrt(177.0) - np.sqrt(10.0)))
    assert np.isclose(result[3], result[2])
    assert np.isclose(
        result[4], result[2] * np.cos(np.deg2rad(80.0))
        / np.cos(np.deg2rad(45.0)))
    assert result[4] < result[3]
    assert 1.0 < result[2] < 1.02


def test_hwang_giapis_mechanism_conserves_removed_si_as_unresolved_product():
    mechanism = HwangGiapisClSiMechanism()
    state = mechanism.initial_state((2,))
    events = FaceResolvedEnergeticFlux(
        "Cl+", 2, np.array([0, 1]), np.array([2e20, 3e20]),
        np.array([177.0, 10.0]), np.ones(2))

    result = mechanism.advance(state, SurfaceFluxes({}, (events,)), 2.0)

    expected_rate = np.array([
        2e20 * HwangGiapisClSiYield().evaluate(177.0, 1.0), 0.0])
    assert np.allclose(result.removed_si_atoms_m2, 2.0 * expected_rate)
    assert np.allclose(result.etch_velocity_m_s, expected_rate / 5.0e28)
    assert np.allclose(
        result.material_exchange.removed_units_m2["poly_si_atoms"],
        result.material_exchange.unresolved_units_m2["poly_si_atoms"])
    assert not result.material_exchange.product_routing_complete
    assert set(result.validity.nonpredictive_parameters) == {
        "critical_angle_deg", "yield_prefactor_per_sqrt_eV"}


def test_hwang_giapis_mechanism_accepts_declared_fast_neutral_but_refuses_hidden_flux():
    mechanism = HwangGiapisClSiMechanism()
    neutral_event = FaceResolvedEnergeticFlux(
        "Cl_fast_neutral", 1, [0], [1e20], [100.0], [1.0])
    result = mechanism.advance(
        mechanism.initial_state((1,)), SurfaceFluxes({}, (neutral_event,)), 0.0)
    assert np.allclose(result.removed_si_atoms_m2, 0.0)

    with pytest.raises(ValueError, match="outside declared scope"):
        mechanism.advance(
            mechanism.initial_state((1,)),
            SurfaceFluxes({"Cl": np.array([1e20])}, (neutral_event,)), 1.0)


def test_hwang_giapis_parameter_evidence_is_complete_and_immutable():
    parameters = HwangGiapisClSiParameters.hwang_giapis_1997()
    assert set(parameters.evidence) == {
        "bulk_si_atom_density_m3", "yield_prefactor_per_sqrt_eV",
        "threshold_energy_eV", "critical_angle_deg"}
    with pytest.raises(TypeError):
        parameters.evidence["threshold_energy_eV"] = None
