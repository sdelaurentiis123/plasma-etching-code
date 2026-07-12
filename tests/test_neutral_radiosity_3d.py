import numpy as np
import pytest

from petch.neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    solve_diffuse_neutral_radiosity_3d,
    transport_diffuse_surface_emission_3d,
    transport_surface_product_population_3d,
)
from petch.surface_exchange import SurfaceProductPopulation


def test_open_plane_reacts_and_escapes_without_artificial_reemission():
    result = solve_diffuse_neutral_radiosity_3d(
        [4.0], [2.0], [], [], [], [1.0], [0.25])

    assert np.allclose(result.incident_flux_m2_s, [4.0])
    assert np.allclose(result.reacted_flux_m2_s, [1.0])
    assert np.isclose(result.source_rate_s, 8.0)
    assert np.isclose(result.reacted_rate_s, 2.0)
    assert np.isclose(result.escaped_rate_s, 6.0)
    assert result.relative_balance_error < 1e-14


def test_unequal_area_cavity_uses_reciprocity_factor_and_conserves_projectiles():
    # F[0->1]=0.4 and F[1->0]=0.2 satisfy A0*F01=A1*F10 for A=[1,2].
    result = solve_diffuse_neutral_radiosity_3d(
        direct_flux_m2_s=[0.6, 0.8], face_area_m2=[1.0, 2.0],
        source_face=[0, 1], target_face=[1, 0], transfer_fraction=[0.4, 0.2],
        escape_fraction=[0.6, 0.8], reaction_probability=[0.3, 0.7])

    # Direct substitution into H0=.6 + (A1/A0)*.2*(1-.7)H1,
    # H1=.8 + (A0/A1)*.4*(1-.3)H0.
    expected_h0 = 0.696 / 0.9832
    expected_h1 = 0.8 + 0.14 * expected_h0
    assert np.allclose(result.incident_flux_m2_s, [expected_h0, expected_h1])
    assert np.isclose(
        result.source_rate_s, result.reacted_rate_s + result.escaped_rate_s,
        rtol=2e-13)
    assert result.relative_linear_residual < 1e-12


def test_form_factor_rows_must_close_projectile_balance_before_solving():
    with pytest.raises(ValueError, match="sum to one"):
        solve_diffuse_neutral_radiosity_3d(
            [1.0], [1.0], [], [], [], [0.9], [0.5])


def test_perfectly_reflecting_closed_cavity_refuses_singular_steady_state():
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_diffuse_neutral_radiosity_3d(
            [1.0, 0.0], [1.0, 1.0], [0, 1], [1, 0], [1.0, 1.0],
            [0.0, 0.0], [0.0, 0.0], maximum_iterations=20)


def test_surface_emission_accounts_for_escape_before_first_impact():
    # A0=1, A1=2 and F[0->1]=0.4. One emitted unit/area/time on face 0 therefore gives
    # first incident density 1*0.4/2=0.2 on face 1; 0.6 escapes directly.
    factors = DiffuseFormFactors3D(
        face_count=2, source_face=np.array([0, 1]), target_face=np.array([1, 0]),
        transfer_fraction=np.array([0.4, 0.2]), escape_fraction=np.array([0.6, 0.8]),
        rays_per_face=10)
    result = transport_diffuse_surface_emission_3d(
        emitted_flux_m2_s=np.array([1.0, 0.0]), face_area_m2=np.array([1.0, 2.0]),
        form_factors=factors, reaction_probability=np.array([1.0, 1.0]))

    assert np.allclose(result.first_incident_flux_m2_s, [0.0, 0.2])
    assert np.isclose(result.emitted_rate_s, 1.0)
    assert np.isclose(result.escaped_without_impact_rate_s, 0.6)
    assert np.isclose(result.reacted_rate_s, 0.4)
    assert np.isclose(result.escaped_after_reflection_rate_s, 0.0)
    assert result.relative_balance_error < 1e-14


def test_surface_emission_multiple_impacts_conserve_material():
    factors = DiffuseFormFactors3D(
        face_count=2, source_face=np.array([0, 1]), target_face=np.array([1, 0]),
        transfer_fraction=np.array([0.4, 0.2]), escape_fraction=np.array([0.6, 0.8]),
        rays_per_face=10)
    result = transport_diffuse_surface_emission_3d(
        emitted_flux_m2_s=np.array([1.0, 0.0]), face_area_m2=np.array([1.0, 2.0]),
        form_factors=factors, reaction_probability=np.array([0.25, 0.5]))

    assert np.all(result.total_incident_flux_m2_s >= result.first_incident_flux_m2_s)
    assert np.isclose(
        result.emitted_rate_s,
        result.reacted_rate_s + result.escaped_without_impact_rate_s
        + result.escaped_after_reflection_rate_s,
        rtol=2e-13)
    assert result.relative_balance_error < 2e-13


def test_named_surface_product_population_reuses_conservative_emission_operator():
    factors = DiffuseFormFactors3D(
        face_count=2, source_face=np.array([0, 1]), target_face=np.array([1, 0]),
        transfer_fraction=np.array([0.4, 0.2]), escape_fraction=np.array([0.6, 0.8]),
        rays_per_face=10)
    product = SurfaceProductPopulation(
        "Si", "Si_atom", integrated_particle_count_m2=np.array([2.0, 0.0]),
        material_units_per_particle=1.0, mass_amu=28.085,
        angular_model="diffuse_cosine", energy_model="thompson",
        energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0},
        provenance={"source": "manufactured transport bridge"})
    result = transport_surface_product_population_3d(
        product, duration_s=2.0, face_area_m2=np.array([1.0, 2.0]),
        form_factors=factors, reaction_probability=np.ones(2))

    assert np.isclose(result.emitted_rate_s, 1.0)
    assert np.isclose(result.reacted_rate_s, 0.4)
    assert np.isclose(result.escaped_without_impact_rate_s, 0.6)


def test_surface_product_transport_refuses_missing_or_unsupported_launch_law():
    factors = DiffuseFormFactors3D(
        face_count=1, source_face=np.array([], dtype=int),
        target_face=np.array([], dtype=int), transfer_fraction=np.array([]),
        escape_fraction=np.array([1.0]), rays_per_face=1)
    unresolved = SurfaceProductPopulation(
        "Si", "Si_atom", [1.0], 1.0, 28.085,
        provenance={"source": "yield-only table"})
    with pytest.raises(ValueError, match="lacks"):
        transport_surface_product_population_3d(
            unresolved, 1.0, [1.0], factors, [1.0])
    unsupported = SurfaceProductPopulation(
        "Si", "Si_atom", [1.0], 1.0, 28.085,
        angular_model="cosine_power_3", energy_model="thompson",
        energy_parameters={
            "surface_binding_energy_eV": 4.7, "maximum_energy_eV": 100.0})
    with pytest.raises(ValueError, match="cannot consume"):
        transport_surface_product_population_3d(
            unsupported, 1.0, [1.0], factors, [1.0])
