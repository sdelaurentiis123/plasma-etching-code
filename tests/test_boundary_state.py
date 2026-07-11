import numpy as np
import pytest

from petch.boundary_state import (
    PlasmaBoundaryState,
    IonEnergyTransverseMaxwellianDensity,
    MaxwellianFluxVelocityDensity,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    instantaneous_sinusoidal_ion_boundary_state,
    maxwellian_electron_boundary_state,
)
from petch.sheath import CollisionlessRFSheath, ECHARGE


def test_species_boundary_normalizes_weights_and_energy():
    species = SpeciesBoundaryState(
        name="Ar+", charge_number=1, mass_amu=40.0, flux_m2_s=2e19,
        velocity_sqrt_eV=np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 4.0]]),
        weight=np.array([1.0, 3.0]),
    )
    assert np.isclose(species.weight.sum(), 1.0)
    assert np.isclose(species.mean_energy_eV, 0.25 * 14.0 + 0.75 * 16.0)
    with pytest.raises(ValueError, match="nonnegative"):
        SpeciesBoundaryState("bad", 1, 1.0, 1.0, [[0, 0, -1]], [1])


def test_plasma_boundary_current_uses_signed_species_fluxes():
    velocity = np.array([[0.0, 0.0, 1.0]])
    ion = SpeciesBoundaryState("ion", 1, 40.0, 2e19, velocity, [1.0])
    electron = SpeciesBoundaryState("electron", -1, 5.4858e-4, 1.5e19, velocity, [1.0])
    state = PlasmaBoundaryState((ion, electron), reference_plane_m=0.0)
    assert np.isclose(state.current_density_A_m2, ECHARGE * 0.5e19)


def test_both_sheath_models_produce_the_same_boundary_contract():
    finite = collisionless_sheath_boundary_state(
        CollisionlessRFSheath(80.0, 0.0, 1e6, 4.0, 40.0, thickness_m=1e-3),
        1e19, n_phase=32, ion_name="Ar+")
    instant = instantaneous_sinusoidal_ion_boundary_state(
        80.0, 0.0, 4.0, 40.0, 1e19, n_phase=32, ion_name="Ar+")
    assert isinstance(finite, PlasmaBoundaryState)
    assert isinstance(instant, PlasmaBoundaryState)
    assert finite.get("Ar+").velocity_sqrt_eV.shape == instant.get("Ar+").velocity_sqrt_eV.shape
    assert np.isclose(finite.get("Ar+").mean_energy_eV, 82.0, atol=0.08)
    assert np.isclose(instant.get("Ar+").mean_energy_eV, 82.0, atol=1e-12)


def test_histogram_density_is_normalized_joint_support_for_adjoint_scoring():
    edges = (np.array([-1.0, 0.0, 2.0]), np.array([-2.0, 1.0]), np.array([0.0, 1.0, 3.0]))
    mass = np.arange(1, 2 * 1 * 2 + 1, dtype=float).reshape(2, 1, 2)
    density = RectilinearVelocityHistogramDensity(edges, mass)
    integral = 0.0
    for i in range(2):
        for j in range(1):
            for k in range(2):
                midpoint = np.array([[(edges[0][i] + edges[0][i + 1]) / 2,
                                      (edges[1][j] + edges[1][j + 1]) / 2,
                                      (edges[2][k] + edges[2][k + 1]) / 2]])
                value = np.exp(density.log_flux_density(midpoint))[0]
                volume = np.diff(edges[0])[i] * np.diff(edges[1])[j] * np.diff(edges[2])[k]
                integral += value * volume
    assert np.isclose(integral, 1.0)
    assert np.isneginf(density.log_flux_density([[0.0, 0.0, -0.1]])[0])


def test_species_exposes_same_density_contract_to_adjoint_consumers():
    density = RectilinearVelocityHistogramDensity(
        (np.array([-1, 1]), np.array([-1, 1]), np.array([0, 2])), np.ones((1, 1, 1)))
    species = SpeciesBoundaryState("ion", 1, 40.0, 1e19, [[0, 0, 1]], [1], density_model=density)
    assert np.isfinite(species.log_flux_density([[0.0, 0.0, 1.0]])[0])


def test_finite_transit_sheath_builds_normalized_continuous_ion_density():
    sheath = CollisionlessRFSheath(80.0, 20.0, 4e5, 4.0, 40.0, thickness_m=1e-3)
    state = collisionless_sheath_boundary_state(
        sheath, 1e19, n_phase=64, tangential_temperature_eV=0.1,
        n_transverse=3, normal_energy_bins=16)
    ion = state.get("ion")
    assert isinstance(ion.density_model, IonEnergyTransverseMaxwellianDensity)
    assert ion.velocity_sqrt_eV.shape == (64 * 3 * 3, 3)
    assert np.isclose(ion.weight.sum(), 1.0)
    assert np.all(np.isfinite(ion.log_flux_density(ion.velocity_sqrt_eV)))


def test_electron_boundary_is_analytic_half_maxwellian_flux_quadrature():
    state = maxwellian_electron_boundary_state(
        4.0, 2e19, n_transverse=5, n_normal=8, electron_name="e-")
    electron = state.get("e-")
    assert isinstance(electron.density_model, MaxwellianFluxVelocityDensity)
    assert electron.velocity_sqrt_eV.shape == (5 * 5 * 8, 3)
    assert np.isclose(electron.weight.sum(), 1.0)
    # Flux-weighted half-Maxwellian: T/2 in each tangent plus T in the normal direction.
    assert np.isclose(electron.mean_energy_eV, 8.0, atol=1e-12)
    assert np.all(np.isfinite(electron.log_flux_density(electron.velocity_sqrt_eV)))
    RectilinearVelocityHistogramDensity,
