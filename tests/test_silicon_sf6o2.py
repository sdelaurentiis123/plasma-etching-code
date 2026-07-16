import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import FeatureGeometry3D, advance_feature_step_3d
from petch.silicon_sf6o2 import (
    BelenSiliconParameters, BelenSiliconSF6O2Mechanism, BelenSiliconState,
)
from petch.surface_kinetics import (
    EnergeticFlux, ParameterEvidence, SteinbruchelYield, SurfaceFluxes,
)


_INPUTS = {
    "site_density_m2",
    "bulk_si_atom_density_m3",
    "fluorine_sticking_probability",
    "oxygen_sticking_probability",
    "spontaneous_fluorine_removal_rate_m2_s",
    "oxygen_desorption_rate_m2_s",
    "physical_sputter_yield",
    "ion_enhanced_yield",
    "oxygen_sputter_yield",
    "fluorine_atoms_per_removed_si",
    "ion_enhanced_fluorine_release_per_si",
}


def _mechanism(*, predictive=True):
    evidence = {
        name: ParameterEvidence(
            "manufactured Belen silicon gate", "analytic",
            supports_prediction_within_declared_domain=predictive)
        for name in _INPUTS}
    bounds = {name: (0.0, 1.0e30) for name in _INPUTS}
    bounds.update({
        "physical_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 1.0),
            "threshold_energy_eV": (0.0, 100.0)},
        "ion_enhanced_yield": {
            "prefactor_per_sqrt_eV": (0.0, 10.0),
            "threshold_energy_eV": (0.0, 100.0)},
        "oxygen_sputter_yield": {
            "prefactor_per_sqrt_eV": (0.0, 10.0),
            "threshold_energy_eV": (0.0, 100.0)},
    })
    return BelenSiliconSF6O2Mechanism(BelenSiliconParameters(
        material_name="Si", material_inventory_name="Si_atom",
        fluorine_species="F", oxygen_species="O", projectile_species=("ion",),
        site_density_m2=5e18, bulk_si_atom_density_m3=5e28,
        fluorine_sticking_probability=0.5, oxygen_sticking_probability=0.25,
        spontaneous_fluorine_removal_rate_m2_s=4e19,
        oxygen_desorption_rate_m2_s=2e19,
        physical_sputter_yield=SteinbruchelYield(0.1, 25.0),
        ion_enhanced_yield=SteinbruchelYield(0.2, 25.0),
        oxygen_sputter_yield=SteinbruchelYield(0.3, 25.0),
        fluorine_atoms_per_removed_si=4.0,
        ion_enhanced_fluorine_release_per_si=2.0,
        evidence=evidence, parameter_bounds=bounds))


def _ions(flux=1e18, energy=100.0):
    return EnergeticFlux("ion", flux, [energy], [1.0], [1.0])


def test_steinbruchel_yield_is_exact_square_root_threshold_law():
    law = SteinbruchelYield(
        0.2, 25.0, angular_model="chang_sawin_1997")

    assert law.evaluate(24.0, 1.0) == 0.0
    assert np.isclose(law.evaluate(100.0, 1.0), 1.0)
    assert np.isclose(law.evaluate(100.0, 0.5), 1.0)
    assert law.evaluate(100.0, 0.0) == 0.0


def test_belen_coverages_solve_the_coupled_site_balances_and_update_at_zero_duration():
    mechanism = _mechanism()
    state = BelenSiliconState.bare((2,))
    fluxes = SurfaceFluxes(
        {"F": np.array([2e20, 4e20]), "O": np.array([1e20, 2e20])},
        (_ions(np.array([1e18, 2e18])),))

    result = mechanism.advance(state, fluxes, 0.0)

    assert np.allclose(
        result.fluorine_coverage + result.oxygen_coverage
        + result.available_site_fraction, 1.0)
    assert np.max(np.abs(result.fluorine_site_balance_residual_m2_s)) / 4e20 < 2e-16
    assert np.max(np.abs(result.oxygen_site_balance_residual_m2_s)) / 2e20 < 2e-16
    assert np.array_equal(result.state.removed_si_atoms_m2, np.zeros(2))
    assert np.array_equal(result.state.available_site_fraction, result.available_site_fraction)
    assert np.all(result.etch_velocity_m_s > 0.0)


def test_neutral_transport_sink_uses_the_same_available_sites_as_the_chemistry():
    mechanism = _mechanism()
    state = BelenSiliconState([0.2, 0.8], [0.0, 0.0])
    probability = mechanism.neutral_reaction_probability(state)

    assert np.allclose(probability["F"], [0.1, 0.4])
    assert np.allclose(probability["O"], [0.05, 0.2])


def test_direct_fluorine_channel_etches_without_ions_and_closes_target_ledger():
    mechanism = _mechanism()
    result = mechanism.advance(
        BelenSiliconState.bare(), SurfaceFluxes({"F": 2e20, "O": 0.0}), 2.0)

    expected_theta_f = (0.5 * 2e20 / 4e19) / (1.0 + 0.5 * 2e20 / 4e19)
    expected_rate = 4e19 * expected_theta_f / 4.0
    assert np.isclose(result.fluorine_coverage, expected_theta_f)
    assert np.isclose(result.chemical_removal_rate_m2_s, expected_rate)
    assert result.physical_sputter_rate_m2_s == 0.0
    assert result.ion_enhanced_removal_rate_m2_s == 0.0
    assert np.isclose(result.removed_si_atoms_m2, 2.0 * expected_rate)
    assert np.isclose(result.etch_velocity_m_s, expected_rate / 5e28)
    assert not result.material_exchange.product_routing_complete
    assert result.material_exchange.residual_units_m2("Si_atom") == 0.0


def test_ion_only_limit_retains_physical_sputtering_but_not_fluorine_enhancement():
    result = _mechanism().advance(
        BelenSiliconState.bare(), SurfaceFluxes({}, (_ions(),)), 1.0)

    assert result.fluorine_coverage == 0.0
    assert result.chemical_removal_rate_m2_s == 0.0
    assert result.ion_enhanced_removal_rate_m2_s == 0.0
    assert result.physical_sputter_rate_m2_s > 0.0


def test_belen_mechanism_refuses_undeclared_flux_and_exposes_sources_and_bounds():
    mechanism = _mechanism(predictive=False)
    with pytest.raises(ValueError, match="no declared silicon reaction channel"):
        mechanism.advance(
            BelenSiliconState.bare(), SurfaceFluxes({"photon": 1e18}), 1.0)

    provenance = mechanism.provenance
    assert set(provenance["sources"]) == _INPUTS
    assert set(provenance["bounds"]) == _INPUTS
    assert provenance["parameters"]["physical_sputter_yield"]["type"] == (
        "SteinbruchelYield")
    assert "no photon-assisted reaction channel is declared" in provenance["known_omissions"]
    assert not mechanism.advance(
        BelenSiliconState.bare(), SurfaceFluxes({"F": 1e18}), 0.0
    ).validity.parameter_evidence_supports_prediction

    parameters = mechanism.parameters
    invalid_bounds = dict(parameters.parameter_bounds)
    invalid_bounds["fluorine_sticking_probability"] = (0.0, 0.1)
    with pytest.raises(ValueError, match="outside declared bounds"):
        BelenSiliconParameters(**{
            name: getattr(parameters, name)
            for name in parameters.__dataclass_fields__
            if name not in {"parameter_bounds"}}, parameter_bounds=invalid_bounds)


def test_common_engine_converges_declared_neutral_surface_fixed_point_before_motion():
    dx = 0.25
    z = np.arange(8) * dx
    phi = np.broadcast_to(0.95 - z, (4, 4, 8)).copy()
    geometry = FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), dx, 1e-6)
    boundary = PlasmaBoundaryState((
        SpeciesBoundaryState("F", 0, 19.0, 2e20, [[0.0, 0.0, 1.0]], [1.0]),
        SpeciesBoundaryState("O", 0, 16.0, 1e20, [[0.0, 0.0, 1.0]], [1.0]),
    ), reference_plane_m=1.75e-6)

    result = advance_feature_step_3d(
        geometry, boundary, {"F": "neutral_reactant", "O": "neutral_reactant"},
        _mechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 0.75, 0.0, 0.75), source_z=1.75,
        n_position=32, seed=3, reinitialize=False, transport_device="cpu",
        neutral_radiosity_options={"rays_per_face": 8, "seed": 5},
        neutral_surface_fixed_point_tolerance=1e-8,
        neutral_surface_fixed_point_max_iterations=8)

    assert result.diagnostics["neutral_surface_fixed_point_iterations"] >= 2
    assert result.diagnostics["neutral_surface_fixed_point_residual"] <= 1e-8
    assert np.max(np.abs(result.surface.transport_fixed_point_change)) <= 1e-8


def test_common_engine_refuses_fixed_point_request_for_non_quasisteady_mechanism():
    class NotQuasiSteady:
        pass

    with pytest.raises(ValueError, match="does not declare"):
        advance_feature_step_3d(
            object(), object(), {}, NotQuasiSteady(), etchable_material_ids=(1,),
            duration_s=0.0, source_bounds=(0.0, 1.0, 0.0, 1.0), source_z=1.0,
            neutral_radiosity_options={}, neutral_surface_fixed_point_tolerance=1e-4)
