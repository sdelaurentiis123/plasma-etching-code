import numpy as np
import pytest

from petch.boundary_state import (
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    instantaneous_sinusoidal_ion_boundary_state,
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
