import numpy as np
import pytest

from petch.boundary_transport_3d import (
    BoundaryTransport3DResult,
    trace_charged_surface_events_field_3d,
)
from petch.charged_surface_cascade_3d import (
    augment_transport_with_charged_reimpacts_3d,
    solve_charged_surface_cascade_3d,
)
from petch.charged_surface_response_3d import (
    ChargedSurfaceContext3D,
    GrazingSpecularIonReflection3D,
    OutgoingChargedParticleEvents3D,
    PerfectAbsorberChargedSurfaceResponse3D,
    account_charged_surface_transfer_3d,
)
from petch.surface_kinetics import ParameterEvidence, SurfaceFluxes
from petch.charging_poisson_3d import lump_triangle_sheet_charge_3d
from petch.sheath import ECHARGE
from petch.surface_kinetics import FaceResolvedEnergeticFlux


def _parallel_triangle_geometry():
    verts = np.array([
        [0.25, 0.0, 0.0], [0.25, 1.0, 0.0], [0.25, 0.0, 1.0],
        [0.75, 0.0, 0.0], [0.75, 1.0, 0.0], [0.75, 0.0, 1.0],
    ])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.5)
    normals = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    context = ChargedSurfaceContext3D(
        areas * 1e-12, normals, np.array(["SiO2", "SiO2"]))
    return verts, faces, areas, context


def _impact(name, face, rate_s, context, direction, energy_eV=1.0):
    position = np.array([[0.25 if face == 0 else 0.75, 1.0 / 3.0, 1.0 / 3.0]])
    return FaceResolvedEnergeticFlux(
        name, 2, event_face=[face],
        event_flux_m2_s=[rate_s / context.face_area_m2[face]],
        event_energy_eV=[energy_eV], event_cosine_incidence=[1.0],
        event_position=position, event_incident_direction=[direction])


class _PerfectSpecularResponse:
    def evaluate(self, incident_populations, charge_number_by_species, context):
        outgoing = []
        for population in incident_populations:
            normal = context.face_gas_normal[population.event_face]
            direction = population.event_incident_direction
            reflected = direction - 2.0 * np.sum(direction * normal, axis=1)[:, None] * normal
            outgoing.append(OutgoingChargedParticleEvents3D(
                population.name, charge_number_by_species[population.name],
                population.face_count, population.event_face,
                population.event_flux_m2_s * context.face_area_m2[population.event_face],
                population.event_position,
                np.sqrt(population.event_energy_eV)[:, None] * reflected))
        return account_charged_surface_transfer_3d(
            incident_populations, charge_number_by_species, context.face_area_m2,
            outgoing=tuple(outgoing))


class _LambertianElectronPerIonResponse:
    """Manufactured one-electron yield with deterministic cosine-weighted directions."""

    def __init__(self, directions_per_impact=64, energy_eV=2.0):
        self.directions_per_impact = int(directions_per_impact)
        self.energy_eV = float(energy_eV)
        self.last_outgoing = None

    def evaluate(self, incident_populations, charge_number_by_species, context):
        source_face = []
        event_rate = []
        event_position = []
        event_velocity = []
        count = self.directions_per_impact
        for population in incident_populations:
            if charge_number_by_species[population.name] <= 0:
                continue
            for event in range(len(population.event_face)):
                face = population.event_face[event]
                normal = context.face_gas_normal[face]
                reference = np.array([0.0, 0.0, 1.0]) if abs(normal[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
                tangent_a = np.cross(reference, normal)
                tangent_a /= np.linalg.norm(tangent_a)
                tangent_b = np.cross(normal, tangent_a)
                u = (np.arange(count) + 0.5) / count
                azimuth = 2.0 * np.pi * np.mod(np.arange(count) * 0.6180339887498949, 1.0)
                cosine = np.sqrt(u)
                sine = np.sqrt(1.0 - u)
                direction = (
                    cosine[:, None] * normal
                    + sine[:, None] * (
                        np.cos(azimuth)[:, None] * tangent_a
                        + np.sin(azimuth)[:, None] * tangent_b))
                rate = (
                    population.event_flux_m2_s[event] * context.face_area_m2[face] / count)
                source_face.extend([face] * count)
                event_rate.extend([rate] * count)
                event_position.extend([population.event_position[event]] * count)
                event_velocity.extend(np.sqrt(self.energy_eV) * direction)
        outgoing = ()
        if event_rate:
            self.last_outgoing = OutgoingChargedParticleEvents3D(
                "secondary_electron", -1, len(context.face_area_m2), source_face,
                event_rate, event_position, event_velocity)
            outgoing = (self.last_outgoing,)
        return account_charged_surface_transfer_3d(
            incident_populations, charge_number_by_species, context.face_area_m2,
            outgoing=outgoing)


def _solve(incident, charge, response, *, max_bounces=16, fixed_dt=0.05, launch_offset=1e-4):
    verts, faces, areas, context = _parallel_triangle_geometry()
    return solve_charged_surface_cascade_3d(
        (incident,), {incident.name: charge}, response, context, verts, faces, areas,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=launch_offset, fixed_dt=fixed_dt, max_steps=200,
        max_bounces=max_bounces, device="cpu")


def test_perfect_absorber_closes_one_response_with_exact_charge_ledger():
    _, _, _, context = _parallel_triangle_geometry()
    incident = _impact("electron", 0, 2.5e7, context, [-1.0, 0.0, 0.0])
    result = _solve(incident, -1, PerfectAbsorberChargedSurfaceResponse3D())

    assert result.completed
    assert len(result.transfers) == 1
    assert result.flights_by_bounce == ()
    assert np.isclose(
        result.initial_incident_charge_rate_c_s,
        result.deposited_charge_rate_c_s, rtol=3e-16)
    assert result.escaped_charge_rate_c_s == 0.0
    assert result.unresolved_charge_rate_c_s == 0.0
    assert result.relative_charge_balance_error < 3e-16


def test_closed_specular_cavity_keeps_cap_remainder_in_charge_ledger():
    _, _, _, context = _parallel_triangle_geometry()
    incident = _impact("electron", 0, 2.5e7, context, [-1.0, 0.0, 0.0])
    result = _solve(incident, -1, _PerfectSpecularResponse(), max_bounces=4)

    assert not result.completed
    assert len(result.transfers) == 4
    assert len(result.flights_by_bounce) == 4
    assert result.deposited_charge_rate_c_s == 0.0
    assert result.escaped_charge_rate_c_s == 0.0
    assert result.unresolved_incident[0].event_face.tolist() == [0]
    assert np.isclose(
        result.unresolved_charge_rate_c_s,
        result.initial_incident_charge_rate_c_s, rtol=3e-16)
    assert result.relative_charge_balance_error < 3e-16
    assert all(item.relative_charge_balance_error < 5e-16 for item in result.transfers)


def test_one_lambertian_electron_per_ion_closes_charge_and_q1_projection():
    verts, faces, _, context = _parallel_triangle_geometry()
    incident_rate = 2.5e7
    incident = _impact("Ar+", 0, incident_rate, context, [-1.0, 0.0, 0.0])
    response = _LambertianElectronPerIonResponse()
    result = _solve(incident, 1, response)

    emitted = response.last_outgoing
    assert emitted is not None
    assert np.isclose(np.sum(emitted.event_rate_s), incident_rate, rtol=2e-16)
    cosine = np.einsum(
        "rc,rc->r", emitted.event_velocity_sqrt_eV,
        context.face_gas_normal[emitted.source_face]) / np.sqrt(response.energy_eV)
    assert np.all(cosine > 0.0)
    assert np.isclose(np.mean(cosine), 2.0 / 3.0, rtol=2e-3)
    assert result.completed
    assert np.isclose(result.deposited_charge_rate_c_s + result.escaped_charge_rate_c_s,
                      ECHARGE * incident_rate, rtol=5e-15)
    projected = lump_triangle_sheet_charge_3d(
        (2, 2, 2), verts, faces, result.face_current_density_a_m2,
        grid_origin=(0.0, 0.0, 0.0), grid_spacing=1.0,
        coordinate_length_unit_m=1e-6)
    assert np.isclose(np.sum(projected), result.deposited_charge_rate_c_s, rtol=5e-15)
    assert result.relative_charge_balance_error < 5e-15
    assert all(item.relative_charge_balance_error < 5e-16 for item in result.transfers)


def test_reimpact_face_and_charge_are_stable_under_timestep_and_offset_refinement():
    _, _, _, context = _parallel_triangle_geometry()
    incident = _impact("electron", 0, 2.5e7, context, [-1.0, 0.0, 0.0])
    coarse = _solve(
        incident, -1, _PerfectSpecularResponse(), max_bounces=1,
        fixed_dt=0.05, launch_offset=1e-4)
    refined = _solve(
        incident, -1, _PerfectSpecularResponse(), max_bounces=1,
        fixed_dt=0.025, launch_offset=5e-5)

    assert np.array_equal(
        coarse.unresolved_incident[0].event_face,
        refined.unresolved_incident[0].event_face)
    assert np.allclose(
        coarse.unresolved_incident[0].event_energy_eV,
        refined.unresolved_incident[0].event_energy_eV, rtol=0.0, atol=2e-6)
    assert coarse.unresolved_charge_rate_c_s == refined.unresolved_charge_rate_c_s


def test_lambertian_surface_quadrature_refines_the_landing_escape_partition():
    _, _, _, context = _parallel_triangle_geometry()
    incident = _impact("Ar+", 0, 2.5e7, context, [-1.0, 0.0, 0.0])
    level_8 = _solve(incident, 1, _LambertianElectronPerIonResponse(256))
    level_9 = _solve(incident, 1, _LambertianElectronPerIonResponse(512))

    assert np.isclose(
        level_8.deposited_charge_rate_c_s,
        level_9.deposited_charge_rate_c_s, rtol=3e-15)
    assert np.isclose(
        level_8.escaped_charge_rate_c_s,
        level_9.escaped_charge_rate_c_s, rtol=3e-15)
    assert level_8.relative_charge_balance_error < 5e-15
    assert level_9.relative_charge_balance_error < 5e-15


def _reflection_model(*, probability=0.95, exponent=3.0, retention=0.9):
    names = (
        "grazing_reflection_probability", "angular_exponent",
        "energy_retention_fraction")
    return GrazingSpecularIonReflection3D(
        material_id="Si", ion_species_name="Ar+",
        grazing_reflection_probability=probability,
        angular_exponent=exponent, energy_retention_fraction=retention,
        parameter_evidence={
            name: ParameterEvidence("manufactured reflection gate", "analytic")
            for name in names},
        parameter_bounds={name: (0.0, 8.0) for name in names})


def test_grazing_specular_reflection_closes_particle_charge_and_energy_ledgers():
    _, _, _, context = _parallel_triangle_geometry()
    context = ChargedSurfaceContext3D(
        context.face_area_m2, context.face_gas_normal, np.array(["Si", "Si"]))
    direction = np.array([-0.1, 0.0, -np.sqrt(0.99)])
    incident_rate = 2.5e7
    incident = _impact("Ar+", 0, incident_rate, context, direction, energy_eV=100.0)
    incident = FaceResolvedEnergeticFlux(
        incident.name, incident.face_count, incident.event_face,
        incident.event_flux_m2_s, incident.event_energy_eV, [0.1],
        event_position=incident.event_position,
        event_incident_direction=incident.event_incident_direction)
    model = _reflection_model()
    transfer = model.evaluate((incident,), {"Ar+": 1}, context)

    probability = model.reflection_probability([0.1])[0]
    reflected = transfer.outgoing[0]
    reflected_rate = float(np.sum(reflected.event_rate_s))
    absorbed_rate = incident_rate - reflected_rate
    reflected_energy_rate = reflected_rate * 90.0
    assert np.isclose(reflected_rate, probability * incident_rate, rtol=2e-16)
    assert np.isclose(absorbed_rate + reflected_rate, incident_rate, rtol=2e-16)
    assert np.isclose(
        transfer.deposited_charge_rate_c_s / ECHARGE, absorbed_rate, rtol=3e-15)
    assert np.isclose(
        transfer.outgoing_kinetic_energy_rate_eV_s, reflected_energy_rate, rtol=3e-16)
    assert np.isclose(
        transfer.deposited_kinetic_energy_rate_eV_s + reflected_energy_rate,
        100.0 * incident_rate, rtol=3e-16)
    assert transfer.relative_charge_balance_error < 5e-15
    assert transfer.relative_kinetic_energy_balance_error == 0.0
    normal = context.face_gas_normal[0]
    reflected_direction = reflected.event_velocity_sqrt_eV[0] / np.sqrt(90.0)
    assert np.allclose(
        reflected_direction,
        direction - 2.0 * np.dot(direction, normal) * normal, rtol=0.0, atol=2e-16)
    assert np.allclose(
        reflected_direction - np.dot(reflected_direction, normal) * normal,
        direction - np.dot(direction, normal) * normal, rtol=0.0, atol=2e-16)


def test_grazing_reflection_refuses_inconsistent_angle_direction_lineage():
    _, _, _, context = _parallel_triangle_geometry()
    context = ChargedSurfaceContext3D(
        context.face_area_m2, context.face_gas_normal, np.array(["Si", "Si"]))
    incident = _impact(
        "Ar+", 0, 2.5e7, context, [-1.0, 0.0, 0.0], energy_eV=100.0)
    inconsistent = FaceResolvedEnergeticFlux(
        incident.name, incident.face_count, incident.event_face,
        incident.event_flux_m2_s, incident.event_energy_eV, [0.1],
        event_position=incident.event_position,
        event_incident_direction=incident.event_incident_direction)

    with pytest.raises(ValueError, match="inconsistent"):
        _reflection_model().evaluate((inconsistent,), {"Ar+": 1}, context)


def test_grazing_reflection_creates_floor_corner_flux_and_off_switch_removes_it():
    # One vertical wall at x=0 and a floor at z=0. A grazing downward ion reflects from the wall
    # and reaches the floor close to their shared corner, the transport precursor of microtrenching.
    verts = np.array([
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.5)
    normals = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    context = ChargedSurfaceContext3D(
        areas * 1e-12, normals, np.array(["Si", "Si"]))
    direction = np.array([-0.1, 0.0, -np.sqrt(0.99)])
    incident = FaceResolvedEnergeticFlux(
        "Ar+", 2, [0], [2.0e7 / context.face_area_m2[0]], [100.0], [0.1],
        event_position=[[0.0, 0.2, 0.8]], event_incident_direction=[direction])
    reflected = _reflection_model(probability=1.0).evaluate(
        (incident,), {"Ar+": 1}, context)
    absorbed = _reflection_model(probability=0.0).evaluate(
        (incident,), {"Ar+": 1}, context)

    flight = trace_charged_surface_events_field_3d(
        reflected.outgoing, verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=1e-4, fixed_dt=0.01, max_steps=500, device="cpu")[0]
    assert flight.termination.tolist() == [1]
    assert flight.hit_face.tolist() == [1]
    assert flight.incident.event_position[0, 0] < 0.1
    assert absorbed.outgoing == ()
    incident_rate_c_s = ECHARGE * 2.0e7
    assert absorbed.deposited_charge_rate_c_s == incident_rate_c_s


def test_reflected_reimpacts_are_in_the_chemistry_facing_surface_flux():
    verts = np.array([
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.5)
    normals = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    context = ChargedSurfaceContext3D(
        areas * 1e-12, normals, np.array(["Si", "Si"]))
    direction = np.array([-0.1, 0.0, -np.sqrt(0.99)])
    incident = FaceResolvedEnergeticFlux(
        "Ar+", 2, [0], [2.0e7 / context.face_area_m2[0]], [100.0], [0.1],
        event_position=[[0.0, 0.2, 0.8]], event_incident_direction=[direction])
    cascade = solve_charged_surface_cascade_3d(
        (incident,), {"Ar+": 1}, _reflection_model(), context, verts, faces, areas,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=1e-4, fixed_dt=0.01, max_steps=500,
        max_bounces=16, device="cpu")
    primary = BoundaryTransport3DResult(
        SurfaceFluxes({}, (incident,)), {"Ar+": 1.0}, {"Ar+": 0.0}, {"Ar+": 0.0},
        "manufactured first hit", ("no surface reflection or neutral re-emission",))
    effective = augment_transport_with_charged_reimpacts_3d(primary, cascade)
    impacts, = effective.surface_fluxes.energetic_fluxes

    assert cascade.completed
    assert len(impacts.event_face) >= 2
    assert impacts.event_face[0] == 0
    assert 1 in impacts.event_face[1:]
    assert impacts.flux_m2_s[1] > 0.0
    assert np.all(impacts.event_energy_eV[1:] < impacts.event_energy_eV[0])
    assert "charged_surface_reimpact_cascade" in effective.transport_model


def test_reflected_flight_segment_is_reciprocal_under_path_reversal():
    verts = np.array([
        [0.25, 0.0, 0.25], [0.25, 1.0, 0.25], [0.25, 0.0, 1.0],
        [0.25, 0.0, 0.25], [1.0, 0.0, 0.25], [0.25, 1.0, 0.25]])
    faces = np.array([[0, 1, 2], [3, 4, 5]])
    areas = np.full(2, 0.375)
    normals = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    context = ChargedSurfaceContext3D(
        areas * 1e-12, normals, np.array(["Si", "Si"]))
    direction = np.array([-0.1, 0.0, -np.sqrt(0.99)])
    source = np.array([0.25, 0.2, 0.8])
    incident = FaceResolvedEnergeticFlux(
        "Ar+", 2, [0], [2.0e7 / context.face_area_m2[0]], [100.0], [0.1],
        event_position=[source], event_incident_direction=[direction])
    reflected, = _reflection_model(probability=1.0).evaluate(
        (incident,), {"Ar+": 1}, context).outgoing
    launch_offset = 1e-5
    forward, = trace_charged_surface_events_field_3d(
        (reflected,), verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=launch_offset, fixed_dt=0.005, max_steps=1000, device="cpu")
    reverse = type(reflected)(
        reflected.name, reflected.charge_number, reflected.face_count, [1],
        reflected.event_rate_s, forward.incident.event_position,
        -reflected.event_velocity_sqrt_eV)
    backward, = trace_charged_surface_events_field_3d(
        (reverse,), verts, faces, areas, normals,
        nodal_potential_v=np.zeros((3, 3, 3)), potential_origin=(0.0, 0.0, 0.0),
        potential_spacing=0.5, mesh_length_unit_m=1e-6,
        launch_offset=launch_offset, fixed_dt=0.005, max_steps=1000, device="cpu")

    assert forward.hit_face.tolist() == [1]
    assert backward.hit_face.tolist() == [0]
    # Two outward launch offsets accumulate geometrically; at this grazing angle their pathwise
    # displacement is bounded by 12 launch offsets and vanishes under offset refinement.
    assert np.linalg.norm(backward.incident.event_position[0] - source) < 12 * launch_offset
