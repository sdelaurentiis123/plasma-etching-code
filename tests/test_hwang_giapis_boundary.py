from pathlib import Path

import numpy as np
import pytest
from scipy.special import betainc
from scipy.stats import qmc

from petch.boundary_state import (
    EnergyCosineAngleDensity2D, IonEnergyTransverseDensity2D,
    IonEnergyTransverseMaxwellianDensity, MaxwellianFluxVelocityDensity,
)
from petch.experimental_boundary import build_hwang_giapis_1997_boundary_state


DATA = Path(__file__).resolve().parents[1] / "data" / "experimental" / "hwang_giapis_1997"
IEDF = DATA / "fig4a_ion_energy_distribution.csv"
EEDF = DATA / "fig4b_electron_energy_distribution.csv"


def test_hwang_giapis_boundary_uses_source_energy_angle_laws_and_equal_currents():
    boundary = build_hwang_giapis_1997_boundary_state(
        IEDF, EEDF, reference_plane_m=3.7e-6)
    ion = boundary.get("Cl+")
    electron = boundary.get("electron")

    assert isinstance(ion.density_model, IonEnergyTransverseMaxwellianDensity)
    assert isinstance(ion.density_model_2d, IonEnergyTransverseDensity2D)
    assert isinstance(electron.density_model, MaxwellianFluxVelocityDensity)
    assert isinstance(electron.density_model_2d, EnergyCosineAngleDensity2D)
    assert ion.flux_m2_s == electron.flux_m2_s
    assert abs(boundary.current_density_A_m2) < 1e-12
    assert np.isclose(ion.density_model.probability_mass.sum(), 1.0)
    edges = ion.density_model.normal_energy_edges_eV
    center = 0.5 * (edges[:-1] + edges[1:])
    mass = ion.density_model.probability_mass
    assert np.isclose(np.dot(center, mass), 32.90013896, rtol=2e-10)
    assert np.isclose(mass[center < 39.0].sum(), 0.5960852083, rtol=2e-10)
    assert ion.provenance["reported_iadf_hwhm_deg"] == 4.3
    assert electron.provenance["eadf_cosine_power"] == 0.6
    assert boundary.provenance["two_dimensional_source_is_source_faithful"] is True


def test_hwang_giapis_2d_electron_sampler_recovers_digitized_eedf_and_cosine_eadf():
    electron = build_hwang_giapis_1997_boundary_state(
        IEDF, EEDF, reference_plane_m=3.7e-6).get("electron")
    unit = qmc.Sobol(
        electron.flux_sampling_dimension_2d,
        scramble=True,
        seed=117).random_base2(16)
    velocity = electron.sample_flux_velocity_2d(unit)
    energy = np.einsum("rc,rc->r", velocity, velocity)
    angle = np.arctan2(velocity[:, 0], velocity[:, 1])

    assert np.isclose(np.mean(energy), 3.5999861565, rtol=3e-4)
    assert np.isclose(np.mean(energy < 5.0), 0.7453325639, atol=3e-4)
    theta = np.deg2rad(45.0)
    beta_shape = 0.8
    expected_center_fraction = (
        betainc(beta_shape, beta_shape, 0.5 * (1.0 + np.sin(theta)))
        - betainc(beta_shape, beta_shape, 0.5 * (1.0 - np.sin(theta))))
    assert np.isclose(
        np.mean(np.abs(angle) < theta),
        expected_center_fraction,
        atol=3e-4)
    assert abs(np.mean(angle)) < 2e-4


def test_hwang_giapis_boundary_without_fig4b_keeps_legacy_3d_electron_projection_declared():
    boundary = build_hwang_giapis_1997_boundary_state(
        IEDF, reference_plane_m=3.7e-6)
    electron = boundary.get("electron")
    assert electron.density_model_2d is None
    assert boundary.provenance["two_dimensional_source_is_source_faithful"] is False
    assert electron.provenance["supports_prediction_within_declared_benchmark"] is False


def test_hwang_giapis_boundary_checksum_gate_refuses_modified_iedf(tmp_path):
    modified = tmp_path / "fig4a.csv"
    modified.write_bytes(IEDF.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        build_hwang_giapis_1997_boundary_state(
            modified, reference_plane_m=3.7e-6)


def test_hwang_giapis_boundary_checksum_gate_refuses_modified_eedf(tmp_path):
    modified = tmp_path / "fig4b.csv"
    modified.write_bytes(EEDF.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="EEDF checksum mismatch"):
        build_hwang_giapis_1997_boundary_state(
            IEDF, modified, reference_plane_m=3.7e-6)
