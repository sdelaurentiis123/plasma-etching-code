import numpy as np
import pytest

from petch.boundary_transport_3d import trace_charged_surface_events_field_3d
from petch.charged_surface_response_3d import (
    ChargedSurfaceContext3D,
    OutgoingChargedParticleEvents3D,
    Sobolewski2021ArKineticSEE3D,
    account_charged_surface_transfer_3d,
    perfect_absorber_surface_transfer_3d,
    sobolewski_2021_ar_kinetic_see_yield,
)
from petch.sheath import ECHARGE
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _incident(name, face, flux):
    count = len(face)
    return FaceResolvedEnergeticFlux(
        name, 2, event_face=face, event_flux_m2_s=flux,
        event_energy_eV=np.full(count, 10.0),
        event_cosine_incidence=np.ones(count),
        event_position=np.column_stack((np.zeros(count), np.zeros(count), np.zeros(count))),
        event_incident_direction=np.tile([0.0, 0.0, -1.0], (count, 1)))


def _outgoing(name, charge, face, rate):
    count = len(face)
    return OutgoingChargedParticleEvents3D(
        name, charge, 2, source_face=face, event_rate_s=rate,
        event_position=np.column_stack((np.zeros(count), np.zeros(count), np.zeros(count))),
        event_velocity_sqrt_eV=np.tile([0.0, 0.0, 1.0], (count, 1)))


def test_perfect_absorber_exactly_reproduces_incident_signed_current_density():
    area = np.array([2e-12, 5e-12])
    ion = _incident("Ar+", [0, 1], [3e18, 4e18])
    electron = _incident("electron", [0, 1], [1e18, 6e18])
    result = perfect_absorber_surface_transfer_3d(
        (ion, electron), {"Ar+": 1, "electron": -1}, area)

    assert np.array_equal(
        result.positive_deposition_current_density_a_m2, ECHARGE * np.array([3e18, 4e18]))
    assert np.array_equal(
        result.negative_deposition_current_density_a_m2, ECHARGE * np.array([1e18, 6e18]))
    assert np.array_equal(
        result.face_current_density_a_m2, ECHARGE * np.array([2e18, -2e18]))
    assert result.outgoing == ()
    assert result.charge_balance_residual_c_s == 0.0
    assert result.relative_charge_balance_error == 0.0


def test_one_reflected_electron_deposits_zero_net_charge_without_losing_throughput():
    area = np.array([2e-12, 5e-12])
    flux = 3e18
    incident = _incident("electron", [0], [flux])
    reflected = _outgoing("electron_reflected", -1, [0], [flux * area[0]])
    result = account_charged_surface_transfer_3d(
        (incident,), {"electron": -1}, area, outgoing=(reflected,))

    assert abs(result.face_current_density_a_m2[0]) <= (
        2.0 * np.finfo(float).eps * ECHARGE * flux)
    assert np.isclose(
        result.positive_deposition_current_density_a_m2[0], ECHARGE * flux, rtol=2e-16)
    assert np.isclose(
        result.negative_deposition_current_density_a_m2[0], ECHARGE * flux, rtol=2e-16)
    assert np.isclose(
        result.incident_charge_rate_c_s, result.outgoing_charge_rate_c_s, rtol=2e-16)
    assert abs(result.deposited_charge_rate_c_s) <= (
        2.0 * np.finfo(float).eps * ECHARGE * flux * area[0])
    assert result.relative_charge_balance_error < 5e-16


def test_absorbed_electron_emitting_two_true_secondaries_deposits_positive_e():
    area = np.array([2e-12, 5e-12])
    flux = 3e18
    incident = _incident("electron", [0], [flux])
    secondaries = _outgoing("secondary_electron", -1, [0], [2.0 * flux * area[0]])
    result = account_charged_surface_transfer_3d(
        (incident,), {"electron": -1}, area, outgoing=(secondaries,))

    assert np.isclose(result.face_current_density_a_m2[0], ECHARGE * flux)
    assert np.isclose(result.deposited_charge_rate_c_s, ECHARGE * flux * area[0])
    assert result.relative_charge_balance_error < 5e-16


def test_neutralized_ion_emitting_one_electron_deposits_two_positive_charges():
    area = np.array([2e-12, 5e-12])
    flux = 3e18
    incident = _incident("Ar+", [1], [flux])
    secondary = _outgoing("secondary_electron", -1, [1], [flux * area[1]])
    result = account_charged_surface_transfer_3d(
        (incident,), {"Ar+": 1}, area, outgoing=(secondary,))

    assert np.isclose(result.face_current_density_a_m2[1], 2.0 * ECHARGE * flux)
    assert np.isclose(result.deposited_charge_rate_c_s, 2.0 * ECHARGE * flux * area[1])
    assert result.relative_charge_balance_error < 5e-16


def test_nearly_neutral_incident_ledger_scales_roundoff_by_absolute_throughput():
    area = np.array([2e-12, 5e-12])
    ion = _incident("Ar+", [0, 1], [3e18, 4e18])
    electron = _incident("electron", [0, 1], [3e18, 4e18 * (1.0 - 1e-14)])
    result = account_charged_surface_transfer_3d(
        (ion, electron), {"Ar+": 1, "electron": -1}, area)

    assert abs(result.incident_charge_rate_c_s) < 1e-13
    assert result.relative_charge_balance_error < 5e-16


def _parallel_triangle_transport_geometry():
    verts = np.array([
        [0.25, 0.0, 0.0], [0.25, 1.0, 0.0], [0.25, 0.0, 1.0],
        [0.75, 0.0, 0.0], [0.75, 1.0, 0.0], [0.75, 0.0, 1.0],
    ])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.5)
    normals = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    return verts, faces, areas, normals


def test_surface_emitted_charge_reimpacts_with_particle_rate_conserved():
    verts, faces, areas, normals = _parallel_triangle_transport_geometry()
    emitted = OutgoingChargedParticleEvents3D(
        "secondary_electron", -1, 2, source_face=[0], event_rate_s=[2.5e7],
        event_position=[[0.25, 1.0 / 3.0, 1.0 / 3.0]],
        event_velocity_sqrt_eV=[[1.0, 0.0, 0.0]])
    result, = trace_charged_surface_events_field_3d(
        (emitted,), verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6, launch_offset=1e-4,
        fixed_dt=0.05, max_steps=40, device="cpu")

    assert result.incident.event_face.tolist() == [1]
    assert np.isclose(
        result.incident.event_flux_m2_s[0], emitted.event_rate_s[0] / (0.5e-12))
    assert result.emitted_rate_s == result.landed_rate_s
    assert result.escaped_rate_s == 0.0
    assert result.truncated_rate_s == 0.0
    assert result.relative_particle_balance_error == 0.0


def test_surface_emitted_charge_truncation_is_explicit():
    verts, faces, areas, normals = _parallel_triangle_transport_geometry()
    emitted = OutgoingChargedParticleEvents3D(
        "secondary_electron", -1, 2, source_face=[1], event_rate_s=[2.5e7],
        event_position=[[0.75, 1.0 / 3.0, 1.0 / 3.0]],
        event_velocity_sqrt_eV=[[-1.0, 0.0, 0.0]])
    # Launching from face 1 toward face 0 lands; use a one-step horizon to prove that an unresolved
    # flight is rejected rather than disappearing from the charge balance.
    with pytest.raises(RuntimeError, match="exhausted max_steps"):
        trace_charged_surface_events_field_3d(
            (emitted,), verts, faces, areas, normals,
            nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
            potential_spacing=0.5, launch_offset=1e-4, fixed_dt=0.01, max_steps=1,
            device="cpu")

    diagnostic, = trace_charged_surface_events_field_3d(
        (emitted,), verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, launch_offset=1e-4, fixed_dt=0.01, max_steps=1,
        allow_truncation=True, device="cpu")
    assert diagnostic.truncated_rate_s == emitted.event_rate_s[0]
    assert diagnostic.relative_particle_balance_error == 0.0


def test_surface_emitted_charge_escape_is_explicit_and_conservative():
    verts, faces, areas, normals = _parallel_triangle_transport_geometry()
    # Remove the downstream target so the outward flight exits the finite field domain.
    verts = verts[:3]
    faces = faces[:1]
    areas = areas[:1]
    normals = normals[:1]
    emitted = OutgoingChargedParticleEvents3D(
        "secondary_electron", -1, 1, source_face=[0], event_rate_s=[2.5e7],
        event_position=[[0.25, 1.0 / 3.0, 1.0 / 3.0]],
        event_velocity_sqrt_eV=[[1.0, 0.0, 0.0]])
    result, = trace_charged_surface_events_field_3d(
        (emitted,), verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, launch_offset=1e-4, fixed_dt=0.05, max_steps=40,
        device="cpu")

    assert result.landed_rate_s == 0.0
    assert result.escaped_rate_s == emitted.event_rate_s[0]
    assert result.truncated_rate_s == 0.0
    assert result.relative_particle_balance_error == 0.0


def test_sobolewski_2021_kinetic_yield_replays_published_equation_8():
    energy = np.array([10.0, 100.0, 1000.0, 10000.0])
    expected = 0.030 * energy ** 2 / (200.0 + energy) ** 1.5
    assert np.array_equal(sobolewski_2021_ar_kinetic_see_yield(energy), expected)


def test_sobolewski_ar_see_is_material_tagged_lambertian_and_rate_conservative():
    area = np.array([2e-12, 5e-12])
    context = ChargedSurfaceContext3D(
        area, np.tile([0.0, 0.0, 1.0], (2, 1)), np.array(["SiO2", "mask"]))
    ion = _incident("Ar+", [0, 1], [3e18, 4e18])
    response = Sobolewski2021ArKineticSEE3D(
        "SiO2", emission_energy_eV=3.0,
        emission_energy_evidence="bounded 1--5 eV sensitivity around literature 'few eV'",
        angular_log2_samples=8, angular_seed=17)
    result = response.evaluate((ion,), {"Ar+": 1}, context)

    emitted, = result.outgoing
    expected_rate = (
        ion.event_flux_m2_s[0] * area[0]
        * sobolewski_2021_ar_kinetic_see_yield(ion.event_energy_eV[0]))
    assert np.isclose(np.sum(emitted.event_rate_s), expected_rate, rtol=3e-16)
    assert np.all(emitted.source_face == 0)
    cosine = emitted.event_velocity_sqrt_eV[:, 2] / np.sqrt(3.0)
    assert np.all(cosine > 0.0)
    assert np.isclose(np.mean(cosine), 2.0 / 3.0, rtol=2e-3)
    assert response.provenance["yield_scope"].startswith("kinetic Ar+")
    assert result.relative_charge_balance_error < 5e-16


def test_sobolewski_response_is_exact_perfect_absorber_when_material_is_absent():
    area = np.array([2e-12, 5e-12])
    context = ChargedSurfaceContext3D(
        area, np.tile([0.0, 0.0, 1.0], (2, 1)), np.array(["mask", "mask"]))
    ion = _incident("Ar+", [0, 1], [3e18, 4e18])
    response = Sobolewski2021ArKineticSEE3D(
        "SiO2", emission_energy_eV=3.0,
        emission_energy_evidence="bounded literature sensitivity")
    result = response.evaluate((ion,), {"Ar+": 1}, context)
    absorber = perfect_absorber_surface_transfer_3d((ion,), {"Ar+": 1}, area)

    assert np.array_equal(result.face_current_density_a_m2,
                          absorber.face_current_density_a_m2)
    assert result.outgoing == ()
    assert result.charge_balance_residual_c_s == 0.0
