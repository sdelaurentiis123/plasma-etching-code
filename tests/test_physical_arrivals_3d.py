import numpy as np
import pytest

from petch.physical_arrivals_3d import sample_physical_poisson_arrivals_3d
from petch.sheath import ECHARGE
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _population():
    return FaceResolvedEnergeticFlux(
        "ion", 2,
        np.array([0, 0, 1]),
        np.array([2.0, 3.0, 5.0]),
        np.array([10.0, 20.0, 30.0]),
        np.array([1.0, 0.8, 0.6]),
        event_position=np.array([
            [0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        event_incident_direction=np.array([
            [0.0, 0.0, -1.0],
            [0.0, 0.2, -np.sqrt(0.96)],
            [0.0, 0.4, -np.sqrt(0.84)]]))


def test_physical_arrival_sample_closes_count_flux_and_charge_identities():
    area = np.array([0.5, 0.2])
    duration = 2.0
    result = sample_physical_poisson_arrivals_3d(
        _population(), area, duration, seed=71)

    reconstructed = np.bincount(
        result.population.event_face,
        weights=(result.population.event_flux_m2_s
                 * area[result.population.event_face] * duration),
        minlength=2)
    np.testing.assert_array_equal(reconstructed, result.realized_face_count)
    np.testing.assert_allclose(result.expected_face_count, [5.0, 2.0])
    assert result.landed_charge_c(1) == pytest.approx(
        ECHARGE * result.realized_particle_count)
    assert result.landed_charge_c(-1) == pytest.approx(
        -ECHARGE * result.realized_particle_count)
    assert np.array_equal(
        sample_physical_poisson_arrivals_3d(
            _population(), area, duration, seed=71).source_event_count,
        result.source_event_count)


def test_physical_arrival_ensemble_has_poisson_mean_and_variance():
    population = FaceResolvedEnergeticFlux(
        "ion", 1, np.array([0]), np.array([4.0]), np.array([10.0]), np.array([1.0]))
    counts = np.asarray([
        sample_physical_poisson_arrivals_3d(
            population, np.array([0.5]), 1.0, seed=seed).realized_particle_count
        for seed in range(3000)])

    assert counts.mean() == pytest.approx(2.0, abs=0.08)
    assert counts.var(ddof=1) == pytest.approx(2.0, abs=0.15)


def test_physical_arrival_zero_count_population_preserves_schema():
    population = FaceResolvedEnergeticFlux(
        "electron", 2, np.array([1]), np.array([0.0]), np.array([1.0]), np.array([1.0]))
    result = sample_physical_poisson_arrivals_3d(
        population, np.ones(2), 1.0, seed=4)

    assert result.realized_particle_count == 0
    assert result.population.face_count == 2
    assert result.population.event_face.size == 0
    np.testing.assert_array_equal(result.realized_face_count, [0, 0])
