import numpy as np
import pytest
from scipy.stats import qmc

from petch.boundary_state import (
    FoldedNormalTangentialDensity,
    PlasmaBoundaryState,
    IonEnergyTransverseMaxwellianDensity,
    MaxwellianFluxVelocityDensity,
    MixtureBoundaryDensity,
    RectilinearVelocityHistogramDensity,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    folded_normal_tangential_proposal,
    instantaneous_sinusoidal_ion_boundary_state,
    maxwellian_electron_boundary_state,
    mixture_boundary_proposal,
    qmc_boundary_proposal,
    qmc_boundary_proposal_with_auxiliary,
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
    # A continuous phase-to-energy map has connected energy support; density quadrature must not
    # manufacture empty internal bins merely because the transport rule uses few phase nodes.
    assert np.all(ion.density_model.probability_mass > 0.0)


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


def test_mixture_proposal_weights_are_numerical_and_density_is_exactly_scored():
    cold = maxwellian_electron_boundary_state(4.0, 1e19, n_transverse=3, n_normal=4).get("electron")
    broad = maxwellian_electron_boundary_state(40.0, 1e19, n_transverse=3, n_normal=4).get("electron")
    proposal = mixture_boundary_proposal((cold, broad), (0.8, 0.2), name="electron-proposal")
    assert isinstance(proposal.density_model, MixtureBoundaryDensity)
    assert np.isclose(proposal.weight.sum(), 1.0)
    velocity = proposal.velocity_sqrt_eV
    expected = np.logaddexp(
        np.log(0.8) + cold.log_flux_density(velocity),
        np.log(0.2) + broad.log_flux_density(velocity))
    assert np.allclose(proposal.log_flux_density(velocity), expected)
    assert proposal.provenance["role"] == "numerical_multiple_importance_proposal"


def test_qmc_proposals_sample_supported_densities_reproducibly():
    electron = maxwellian_electron_boundary_state(4.0, 1e19).get("electron")
    first = qmc_boundary_proposal(electron, 8, seed=17)
    second = qmc_boundary_proposal(electron, 8, seed=17)
    assert np.array_equal(first.velocity_sqrt_eV, second.velocity_sqrt_eV)
    assert np.all(np.isfinite(first.log_flux_density(first.velocity_sqrt_eV)))
    assert np.isclose(first.mean_energy_eV, 8.0, rtol=0.03)

    mixture = mixture_boundary_proposal((
        electron, maxwellian_electron_boundary_state(40.0, 1.0).get("electron")), (0.9, 0.1))
    sampled = qmc_boundary_proposal(mixture, 6, seed=23)
    assert sampled.velocity_sqrt_eV.shape == (2 * 64, 3)
    assert np.isclose(sampled.weight[:64].sum(), 0.9)
    assert np.isclose(sampled.weight[64:].sum(), 0.1)


def test_joint_qmc_proposal_preserves_mixture_strata_and_auxiliary_alignment():
    cold = maxwellian_electron_boundary_state(4.0, 1.0).get("electron")
    hot = maxwellian_electron_boundary_state(40.0, 1.0).get("electron")
    mixture = mixture_boundary_proposal((cold, hot), (0.8, 0.2))
    first, first_auxiliary = qmc_boundary_proposal_with_auxiliary(
        mixture, 8, 2, seed=31)
    second, second_auxiliary = qmc_boundary_proposal_with_auxiliary(
        mixture, 8, 2, seed=31)

    assert first.velocity_sqrt_eV.shape == (2 * 256, 3)
    assert first_auxiliary.shape == (2 * 256, 2)
    assert np.array_equal(first.velocity_sqrt_eV, second.velocity_sqrt_eV)
    assert np.array_equal(first_auxiliary, second_auxiliary)
    assert np.all((first_auxiliary > 0.0) & (first_auxiliary < 1.0))
    assert np.isclose(first.weight[:256].sum(), 0.8)
    assert np.isclose(first.weight[256:].sum(), 0.2)


@pytest.mark.parametrize("tangent_sign", [-1, 1])
def test_folded_normal_tangential_proposal_is_normalized_source_pushforward(tangent_sign):
    temperature = 4.0
    source = maxwellian_electron_boundary_state(temperature, 1.0).get("electron")
    template = folded_normal_tangential_proposal(source, tangent_sign)
    assert isinstance(template.density_model, FoldedNormalTangentialDensity)
    sampled = qmc_boundary_proposal(template, 15, seed=29)
    velocity = sampled.velocity_sqrt_eV

    assert np.all(tangent_sign * velocity[:, 0] > 0.0)
    assert np.all(velocity[:, 2] > 0.0)
    assert np.all(np.isfinite(sampled.log_flux_density(velocity)))
    # The pushforward maps source normal -> signed local tangent and folds source x -> local normal.
    assert np.allclose(
        np.mean(velocity ** 2, axis=0),
        [temperature, temperature / 2.0, temperature / 2.0], rtol=3e-3)

    probe = np.array([[tangent_sign * 3.0, 0.4, 0.2]])
    source_positive = np.array([[0.2, 0.4, 3.0]])
    source_negative = np.array([[-0.2, 0.4, 3.0]])
    expected = np.logaddexp(
        source.log_flux_density(source_positive),
        source.log_flux_density(source_negative))
    assert np.allclose(template.log_flux_density(probe), expected)


def test_maxwellian_continuous_sampler_reproduces_energy_moments_and_tail():
    temperature = 4.0
    density = MaxwellianFluxVelocityDensity(temperature)
    unit = qmc.Sobol(
        density.sampling_dimension, scramble=True, seed=101).random_base2(14)
    velocity = density.sample_flux_velocity(unit)
    component_energy = np.mean(velocity ** 2, axis=0)
    tail = np.mean(velocity[:, 2] ** 2 > temperature * np.log(10.0))

    assert np.allclose(component_energy, [temperature / 2, temperature / 2, temperature], rtol=3e-3)
    assert np.isclose(tail, 0.1, atol=2e-4)


def test_ion_and_histogram_continuous_samplers_preserve_declared_probability_mass():
    ion_density = IonEnergyTransverseMaxwellianDensity(
        np.array([0.0, 1.0, 3.0]), np.array([0.25, 0.75]), 0.2)
    ion_unit = qmc.Sobol(
        ion_density.sampling_dimension, scramble=True, seed=103).random_base2(14)
    ion_velocity = ion_density.sample_flux_velocity(ion_unit)
    normal_energy = ion_velocity[:, 2] ** 2
    assert np.isclose(np.mean(normal_energy < 1.0), 0.25, atol=1e-4)
    assert np.isclose(np.mean(normal_energy), 1.625, rtol=2e-3)
    assert np.allclose(np.mean(ion_velocity[:, :2] ** 2, axis=0), 0.1, rtol=4e-3)

    histogram = RectilinearVelocityHistogramDensity(
        (np.array([-1.0, 0.0, 2.0]), np.array([-2.0, 1.0]), np.array([0.0, 1.0])),
        np.array([[[0.3]], [[0.7]]]))
    histogram_unit = qmc.Sobol(
        histogram.sampling_dimension, scramble=True, seed=107).random_base2(14)
    histogram_velocity = histogram.sample_flux_velocity(histogram_unit)
    assert np.isclose(np.mean(histogram_velocity[:, 0] < 0.0), 0.3, atol=1e-4)
    assert np.all((histogram_velocity[:, 2] >= 0.0) & (histogram_velocity[:, 2] <= 1.0))


def test_mixture_continuous_sampler_preserves_physical_component_weights():
    density = MixtureBoundaryDensity(
        (MaxwellianFluxVelocityDensity(1.0), MaxwellianFluxVelocityDensity(9.0)),
        np.array([0.75, 0.25]))
    unit = qmc.Sobol(
        density.sampling_dimension, scramble=True, seed=109).random_base2(15)
    velocity = density.sample_flux_velocity(unit)

    # Flux-weighted Maxwellian total kinetic energy is 2T.
    assert np.isclose(np.mean(np.sum(velocity ** 2, axis=1)), 2.0 * 3.0, rtol=4e-3)
