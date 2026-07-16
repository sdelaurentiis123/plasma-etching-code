import numpy as np

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.boundary_transport_3d import BoundaryTransport3DResult
from petch.charged_surface_response_3d import ChargedSurfaceContext3D
from petch.chlorine_poly_si import HwangGiapisClSiMechanism
from petch.feature_step_3d import (
    _face_material_ids, _surface_gas_normals, advance_feature_step_3d,
    make_rectangular_trench_geometry_3d,
)
from petch.hwang_giapis_scatter_3d import (
    HwangGiapisSiO2ForwardScatter3D,
    apply_hwang_giapis_forward_scatter_to_transport_3d,
)
from petch.surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes
from petch.threed import extract_mesh_3d


def _corner_mesh():
    vertices = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
        [1.0, 0.0, 1.0], [1.0, 1.0, 1.0],
    ])
    faces = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]])
    areas = np.full(4, 0.5)
    normals = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0],
                        [-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    material = np.array([2, 2, 1, 1])
    return vertices, faces, areas, normals, material


def test_hwang_giapis_scatter_probability_and_hard_sphere_energy_are_analytic():
    model = HwangGiapisSiO2ForwardScatter3D(2)
    cosine = np.cos(np.deg2rad([30.0, 45.0, 60.0, 90.0]))
    assert np.allclose(model.scattering_probability(cosine), [0.0, 0.0, 1/3, 1.0])
    assert np.allclose(model.energy_retention_fraction(cosine)[1:], [0.0, 0.25, 1.0])


def test_neutralized_sio2_scatter_hits_sidewall_and_closes_particle_energy_ledgers():
    vertices, faces, areas, normals, material = _corner_mesh()
    angle = np.deg2rad(60.0)
    direction = np.array([np.sin(angle), 0.0, -np.cos(angle)])
    incident = FaceResolvedEnergeticFlux(
        "Cl+", 4, [0], [2e20], [100.0], [np.cos(angle)],
        event_position=[[0.2, 0.25, 0.0]], event_incident_direction=[direction])
    transport = BoundaryTransport3DResult(
        SurfaceFluxes({}, (incident,)), {"Cl+": 1.0}, {"Cl+": 0.0},
        {"Cl+": 0.0}, "manufactured", ())
    context = ChargedSurfaceContext3D(areas * 1e-12, normals, material)

    augmented, result = apply_hwang_giapis_forward_scatter_to_transport_3d(
        transport, HwangGiapisSiO2ForwardScatter3D(2), context,
        vertices, faces, areas, domain_minimum=[0.0, 0.0, 0.0],
        domain_maximum=[1.2, 1.0, 1.2], mesh_length_unit_m=1e-6,
        launch_offset=1e-7)

    scattered = next(
        item for item in augmented.surface_fluxes.energetic_fluxes
        if item.name == "Cl_fast_neutral")
    assert result.flight.landed_rate_s == result.scattered_rate_s
    assert result.flight.escaped_rate_s == 0.0
    assert scattered.event_face[0] in (2, 3)
    assert np.isclose(scattered.event_energy_eV[0], 25.0)
    assert result.relative_surface_particle_balance_error < 1e-14
    assert result.relative_surface_energy_balance_error < 1e-14
    assert "neutralized_sio2_forward_scatter" in augmented.transport_model


def test_forward_scatter_off_below_critical_angle_adds_empty_neutral_measure():
    vertices, faces, areas, normals, material = _corner_mesh()
    angle = np.deg2rad(30.0)
    direction = np.array([np.sin(angle), 0.0, -np.cos(angle)])
    incident = FaceResolvedEnergeticFlux(
        "Cl+", 4, [0], [2e20], [100.0], [np.cos(angle)],
        event_position=[[0.2, 0.25, 0.0]], event_incident_direction=[direction])
    transport = BoundaryTransport3DResult(
        SurfaceFluxes({}, (incident,)), {"Cl+": 1.0}, {"Cl+": 0.0},
        {"Cl+": 0.0}, "manufactured", ())
    context = ChargedSurfaceContext3D(areas * 1e-12, normals, material)

    augmented, result = apply_hwang_giapis_forward_scatter_to_transport_3d(
        transport, HwangGiapisSiO2ForwardScatter3D(2), context,
        vertices, faces, areas, domain_minimum=[0.0, 0.0, 0.0],
        domain_maximum=[1.2, 1.0, 1.2], mesh_length_unit_m=1e-6)

    assert result.scattered_rate_s == 0.0
    assert result.flight.relative_particle_balance_error == 0.0
    assert augmented.surface_fluxes.energetic_fluxes[-1].event_face.size == 0


def test_common_feature_step_routes_scattered_fast_neutral_to_same_chemistry_engine():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=1.0, cell_length=0.4, domain_height=2.0, dx=0.1,
        opening_width=0.4, mask_thickness=0.3, substrate_top=1.0,
        etched_depth=0.2, substrate_material_id=1, mask_material_id=2)
    vertices, faces, centroids, areas = extract_mesh_3d(geometry.phi, geometry.dx)
    normals = _surface_gas_normals(vertices, faces, centroids, geometry)
    material = _face_material_ids(centroids, geometry)
    source_face = int(np.flatnonzero(material == 2)[0])
    normal = normals[source_face]
    trial_axis = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial_axis, normal)) > 0.8:
        trial_axis = np.array([0.0, 1.0, 0.0])
    tangent = np.cross(normal, trial_axis)
    tangent /= np.linalg.norm(tangent)
    angle = np.deg2rad(60.0)
    incident_direction = np.sin(angle) * tangent - np.cos(angle) * normal
    ion_events = FaceResolvedEnergeticFlux(
        "Cl+", len(faces), [source_face], [1e20], [100.0], [np.cos(angle)],
        event_position=[centroids[source_face]],
        event_incident_direction=[incident_direction])
    transport = BoundaryTransport3DResult(
        SurfaceFluxes({}, (ion_events,)), {"Cl+": 1.0}, {"Cl+": 0.0},
        {"Cl+": 0.0}, "manufactured", ())
    boundary = PlasmaBoundaryState((SpeciesBoundaryState(
        "Cl+", 1, 35.45, 1e20, [[0.0, 0.0, 10.0]], [1.0]),),
        reference_plane_m=1.8e-6)

    result = advance_feature_step_3d(
        geometry, boundary, {"Cl+": "energetic_bombardment"},
        HwangGiapisClSiMechanism(), etchable_material_ids=(1,), duration_s=0.0,
        source_bounds=(0.0, 1.0, 0.0, 0.4), source_z=1.8,
        precomputed_transport=transport,
        neutral_forward_scatter=HwangGiapisSiO2ForwardScatter3D(2),
        neutral_forward_scatter_options={"launch_offset": 1e-6},
        reinitialize=False, transport_device="cpu")

    names = {item.name for item in result.transport.surface_fluxes.energetic_fluxes}
    assert names == {"Cl+", "Cl_fast_neutral"}
    assert result.neutral_forward_scatter is not None
    assert result.diagnostics["neutral_forward_scatter_applied"] is True
    assert result.diagnostics["neutral_forward_scatter_particle_balance_error"] < 1e-14
