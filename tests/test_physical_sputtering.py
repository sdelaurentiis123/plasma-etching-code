import numpy as np
import pytest

from petch.physical_sputtering import (
    PhysicalSputterMechanism, PhysicalSputterParameters, PhysicalSputterState,
)
from petch.surface_kinetics import EnergeticFlux, EnergeticYield, ParameterEvidence, SurfaceFluxes


def _evidence(predictive=True):
    return {
        name: ParameterEvidence(
            "manufactured physical-sputter gate", "analytic",
            supports_prediction_within_declared_domain=predictive)
        for name in (
            "bulk_material_unit_density_m3", "sputter_yield",
            "emitted_product_mass_amu", "emission_angular_model", "emission_energy_model")}


def _mechanism(predictive=True):
    return PhysicalSputterMechanism(PhysicalSputterParameters(
        material_name="SiO2", material_inventory_name="SiO2_formula_unit",
        projectile_species=("Ar+",), bulk_material_unit_density_m3=2.2e28,
        sputter_yield=EnergeticYield(0.2, 20.0, 100.0),
        emitted_product_name="sputtered_SiO2_unit", emitted_product_mass_amu=60.084,
        emitted_material_units_per_particle=1.0, emission_angular_model="diffuse_cosine",
        emission_energy_model="thompson", emission_energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        evidence=_evidence(predictive)))


def _ions(flux):
    return EnergeticFlux(
        "Ar+", flux, np.array([100.0]), np.array([1.0]), np.array([1.0]))


def test_physical_sputter_routes_every_removed_target_unit_to_outgoing_population():
    mechanism = _mechanism(); state = PhysicalSputterState.bare((2,))
    result = mechanism.advance(
        state, SurfaceFluxes({}, (_ions(np.array([1e18, 2e18])),)), duration_s=2.0)

    expected = np.array([4e17, 8e17])
    assert np.allclose(result.removed_material_units_m2, expected)
    assert np.allclose(result.state.removed_material_units_m2, expected)
    assert np.allclose(result.etch_velocity_m_s, np.array([2e17, 4e17]) / 2.2e28)
    assert result.material_exchange.product_routing_complete
    assert np.array_equal(
        result.product_populations[0].integrated_material_units_m2,
        result.material_exchange.outgoing_units_m2["SiO2_formula_unit"])
    assert result.product_populations[0].transport_ready
    assert result.validity.parameter_evidence_supports_prediction


def test_physical_sputter_refuses_undeclared_incident_species():
    mechanism = _mechanism(); state = PhysicalSputterState.bare()
    with pytest.raises(ValueError, match="no physical-sputter channel"):
        mechanism.advance(state, SurfaceFluxes({"F": 1e18}, (_ions(1e18),)), 1.0)


def test_physical_sputter_keeps_parameter_evidence_honest():
    result = _mechanism(predictive=False).advance(
        PhysicalSputterState.bare(), SurfaceFluxes({}, (_ions(1e18),)), 1.0)

    assert not result.validity.parameter_evidence_supports_prediction
    assert set(result.validity.nonpredictive_parameters) == set(_evidence())
