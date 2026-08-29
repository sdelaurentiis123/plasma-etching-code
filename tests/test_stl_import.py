"""Gates for STL import: reader, exact SDF rasterization, axisymmetric routing."""

import itertools
from types import SimpleNamespace

import numpy as np
import pytest

from petch.boundary_state import PlasmaBoundaryState, SpeciesBoundaryState
from petch.feature_step_3d import advance_feature_step_3d
from petch.stl_import import (
    StlMesh,
    assign_materials_by_z,
    build_feature_geometry_from_stl,
    diagnose_mesh,
    drop_degenerate_faces,
    extract_axisymmetric_profile,
    rasterize_signed_distance,
    read_stl,
    revolved_stl_mesh,
    to_axisymmetric_profile,
    write_stl,
)
from petch.stl_import import _point_triangle_distance_and_solid_angle, _weld
from petch.surface_kinetics import (
    MechanismValidity, SiO2SurfaceState, SurfaceMaterialExchange,
)
from petch.threed import extract_mesh_3d

_BOX_FACES = (
    (0, 3, 2), (0, 1, 3), (4, 6, 7), (4, 7, 5), (0, 5, 1), (0, 4, 5),
    (2, 3, 7), (2, 7, 6), (0, 2, 6), (0, 6, 4), (1, 5, 7), (1, 7, 3))


def _box_mesh(x_range, y_range, z_range):
    """Axis-aligned closed box with outward winding."""
    corners = np.asarray(
        list(itertools.product(x_range, y_range, z_range)), dtype=float)
    triangles = np.asarray(
        [[corners[i], corners[j], corners[k]] for i, j, k in _BOX_FACES])
    vertices, faces = _weld(triangles)
    return StlMesh(vertices, faces)


def _cylinder(radius=1.0, height=2.0, n_theta=64):
    return revolved_stl_mesh(
        np.array([0.0, height]), np.full(2, radius), n_theta=n_theta)


def test_revolved_mesh_is_watertight_outward_and_encloses_the_exact_volume():
    mesh = _cylinder(n_theta=256)
    report = diagnose_mesh(mesh)

    assert report.is_watertight and report.consistently_oriented
    assert report.outward_oriented
    assert report.failure_reason() is None
    # A regular n-gon prism encloses (n/2) R^2 sin(2 pi / n) h, exactly.
    exact = 0.5 * 256 * 1.0 ** 2 * np.sin(2.0 * np.pi / 256) * 2.0
    assert report.signed_volume == pytest.approx(exact, rel=1e-12)


def test_exact_point_triangle_distance_is_a_tight_lower_bound():
    """The closest-point cascade must never overshoot a densely sampled
    reference, and must be tight to the sampling resolution."""
    rng = np.random.default_rng(0)
    tri = rng.normal(size=(6, 3, 3))
    points = 1.5 * rng.normal(size=(150, 3))
    distance, _ = _point_triangle_distance_and_solid_angle(points, tri)

    axis = np.linspace(0.0, 1.0, 200)
    u, v = np.meshgrid(axis, axis, indexing="ij")
    inside = (u + v) <= 1.0
    u, v = u[inside], v[inside]
    reference = np.full(len(points), np.inf)
    for triangle in tri:
        sampled = (triangle[0] + u[:, None] * (triangle[1] - triangle[0])
                   + v[:, None] * (triangle[2] - triangle[0]))
        reference = np.minimum(reference, np.sqrt(np.sum(
            (points[:, None, :] - sampled[None]) ** 2, axis=2)).min(axis=1))

    assert np.all(distance <= reference + 1e-12)
    assert np.max(reference - distance) < 1e-2


def test_cylinder_sdf_is_exact_to_the_tessellation():
    """Inside a convex body the signed distance is min over the bounding
    planes, so the only admissible error is the STL's own faceting."""
    mesh = _cylinder(n_theta=64)
    phi, report = rasterize_signed_distance(mesh, dx=0.1)
    axes = [report.origin[axis] + np.arange(report.shape[axis]) * report.dx
            for axis in range(3)]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    analytic = np.minimum.reduce(
        [1.0 - np.sqrt(x ** 2 + y ** 2), z, 2.0 - z])
    faceting = 1.0 - np.cos(np.pi / 64)

    assert np.all((phi > 0.0) == (analytic > 0.0))
    interior = analytic > 0.0
    assert np.max(np.abs(phi - analytic)[interior]) <= faceting + 1e-12
    # No node may land exactly on the surface: phi == 0 is neither solid nor
    # gas and breaks marching cubes (the default grid is half-cell staggered).
    assert not np.any(phi == 0.0)
    assert report.max_winding_residual < 1e-9


def test_cylinder_marching_cubes_surface_matches_the_analytic_radius():
    mesh = _cylinder(n_theta=64)
    phi, report = rasterize_signed_distance(mesh, dx=0.1)
    _, _, centroids, _ = extract_mesh_3d(phi, report.dx)
    world = centroids + np.asarray(report.origin)
    band = (world[:, 2] > 0.4) & (world[:, 2] < 1.6)
    radius = np.hypot(world[band, 0], world[band, 1])
    wall = radius > 0.7

    assert wall.sum() > 100
    assert np.mean(radius[wall]) == pytest.approx(1.0, abs=1e-2)
    assert np.max(np.abs(radius[wall] - 1.0)) < 0.1 * report.dx


def test_axisymmetry_detector_separates_a_cylinder_from_a_square_prism():
    cylinder = extract_axisymmetric_profile(_cylinder(n_theta=64), n_levels=16)
    prism = extract_axisymmetric_profile(
        _box_mesh((-0.5, 0.5), (-0.5, 0.5), (0.0, 2.0)), n_levels=12)

    assert cylinder.is_axisymmetric
    assert cylinder.relative_deviation <= 1.0 - np.cos(np.pi / 64)
    assert np.allclose(cylinder.r, 1.0, atol=1.0 - np.cos(np.pi / 64))
    assert np.all(cylinder.z > 0.0) and np.all(cylinder.z < 2.0)

    # A square prism is a body of four azimuthal samples: its own faceting
    # bound would excuse it, so the absolute out-of-roundness cap must refuse.
    assert not prism.is_axisymmetric
    assert prism.relative_deviation > 0.1
    with pytest.raises(ValueError, match="not axisymmetric"):
        to_axisymmetric_profile(prism)


def test_tapered_hole_profile_matches_the_analytic_frustum():
    mesh = revolved_stl_mesh(
        np.array([0.0, 3.0]), np.array([1.0, 0.4]), n_theta=96)
    report = extract_axisymmetric_profile(mesh, n_levels=24)
    analytic = 1.0 + (0.4 - 1.0) * report.z / 3.0

    assert report.is_axisymmetric
    assert np.max(np.abs(report.r - analytic)) <= report.facet_bound
    profile = to_axisymmetric_profile(report)
    assert np.all(np.diff(profile.z) > 0.0) and np.all(profile.r > 0.0)


def test_extracted_profile_survives_a_revolve_round_trip():
    original = revolved_stl_mesh(
        np.array([0.0, 1.0, 3.0]), np.array([0.8, 0.5, 0.5]), n_theta=96)
    first = extract_axisymmetric_profile(original, n_levels=32)
    second = extract_axisymmetric_profile(
        revolved_stl_mesh(first.z, first.r, n_theta=96), n_levels=32)

    assert np.max(np.abs(np.interp(second.z, first.z, first.r) - second.r)) < 5e-3


def test_ascii_and_binary_readers_agree(tmp_path):
    """Bitwise on float32-representable coordinates; to float32 epsilon
    otherwise, since binary STL stores single precision by specification."""
    dyadic = _box_mesh((-0.5, 0.5), (-0.5, 0.5), (0.0, 2.0))
    binary = read_stl(write_stl(tmp_path / "d.stl", dyadic, binary=True))
    ascii_mesh = read_stl(write_stl(tmp_path / "d.txt", dyadic, binary=False))

    assert np.array_equal(binary.vertices, ascii_mesh.vertices)
    assert np.array_equal(binary.faces, ascii_mesh.faces)
    assert binary.signed_volume == ascii_mesh.signed_volume == 2.0

    general = _cylinder(n_theta=32)
    from_binary = read_stl(write_stl(tmp_path / "g.stl", general, binary=True))
    from_ascii = read_stl(write_stl(tmp_path / "g.txt", general, binary=False))
    assert np.max(np.abs(np.sort(from_binary.vertices, axis=0)
                         - np.sort(from_ascii.vertices, axis=0))) < 1e-6
    # Format detection must not rely on the leading keyword: a binary file whose
    # 80-byte header begins with "solid" still reads as binary.
    assert len(from_binary.faces) == len(general.faces)


def test_non_watertight_stl_is_rejected_with_a_diagnostic():
    holed = StlMesh(_box_mesh((-0.5, 0.5), (-0.5, 0.5), (0.0, 2.0)).vertices,
                    _box_mesh((-0.5, 0.5), (-0.5, 0.5), (0.0, 2.0)).faces[:-1])
    report = diagnose_mesh(holed)

    assert not report.is_watertight
    assert report.n_boundary_edges == 3
    assert "not watertight" in report.failure_reason()
    with pytest.raises(ValueError, match="not watertight"):
        rasterize_signed_distance(holed, dx=0.2)


def test_isolated_zero_area_exporter_facet_requires_receipted_repair():
    clean = _box_mesh((-0.5, 0.5), (-0.5, 0.5), (0.0, 2.0))
    point = len(clean.vertices)
    raw = StlMesh(
        np.concatenate((clean.vertices, [[0.5, 0.5, 2.0]]), axis=0),
        np.concatenate((clean.faces, [[point, point, point]]), axis=0),
    )

    diagnosis = diagnose_mesh(raw)
    assert diagnosis.n_degenerate_faces == 1
    assert "zero-area" in diagnosis.failure_reason()
    with pytest.raises(ValueError, match="zero-area"):
        rasterize_signed_distance(raw, dx=0.2)

    repaired, receipt = drop_degenerate_faces(raw)
    assert diagnose_mesh(repaired).failure_reason() is None
    assert receipt["removed_face_indices"] == [len(raw.faces) - 1]
    assert receipt["relative_volume_change"] <= 128.0 * np.finfo(float).eps
    assert repaired.signed_volume == clean.signed_volume


def test_material_layers_reconstruct_the_solid_sign_exactly():
    mesh = _box_mesh((-1.0, 1.0), (-1.0, 1.0), (0.0, 2.0))
    phi, report = rasterize_signed_distance(mesh, dx=0.2)
    material, levelsets = assign_materials_by_z(
        phi, dx=report.dx, origin=report.origin,
        layers=[(-10.0, 1.0, 1), (1.0, 10.0, 2)])
    union = np.maximum.reduce(list(levelsets.values()))

    # The invariant FeatureGeometry3D enforces on material level sets.
    assert np.all((union >= 0.0) == (phi >= 0.0))
    assert set(np.unique(material)) == {0, 1, 2}
    assert np.all(material[phi <= 0.0] == 0)
    z = report.origin[2] + np.arange(report.shape[2]) * report.dx
    lower = np.broadcast_to(z < 1.0, phi.shape)
    assert np.all(material[(phi > 0.0) & lower] == 1)
    assert np.all(material[(phi > 0.0) & ~lower] == 2)


def test_exterior_solid_region_inverts_the_body():
    mesh = _cylinder(radius=0.6, height=2.0, n_theta=48)
    void, _ = rasterize_signed_distance(mesh, dx=0.15, solid_region="exterior")
    solid, _ = rasterize_signed_distance(mesh, dx=0.15, solid_region="interior")

    assert np.array_equal(void, -solid)
    assert np.count_nonzero(void > 0.0) > np.count_nonzero(solid > 0.0)


def test_min_feature_cells_reports_an_underresolved_wall():
    """A wall thinner than a couple of cells must be reported, not silently
    imported: the level set cannot represent it."""
    thin = _box_mesh((-0.06, 0.06), (-2.0, 2.0), (0.0, 1.0))
    _, report = rasterize_signed_distance(thin, dx=0.1, padding_cells=3.0)

    assert report.min_feature_cells <= 2.0
    assert "refine dx" in report.resolution_warning()
    coarse = _box_mesh((-1.0, 1.0), (-1.0, 1.0), (0.0, 2.0))
    _, wide = rasterize_signed_distance(coarse, dx=0.1)
    assert wide.resolution_warning() is None


class _ManufacturedDirectionalEtch:
    """Test-only nonnegative directional removal law for the round-trip gate.

    It replaces no numerical operator: it hands the public engine one declared
    velocity proportional to the incident energetic flux, which is enough to
    certify that an STL-built geometry is accepted and advanced.
    """

    density_m3 = 1.0e28
    etch_velocity_m_s = 4.0e-8

    @staticmethod
    def initial_state(shape=()):
        return SiO2SurfaceState.bare(shape)

    def advance(self, state, fluxes, duration_s):
        population = next(item for item in fluxes.energetic_fluxes
                          if item.name == "etch+")
        incident = np.asarray(population.flux_m2_s, dtype=float)
        peak = float(np.max(incident)) if incident.size else 0.0
        etch = (np.zeros(incident.shape) if peak == 0.0
                else self.etch_velocity_m_s * incident / peak)
        removed = etch * self.density_m3 * float(duration_s)
        return SimpleNamespace(
            state=state,
            etch_velocity_m_s=etch,
            normal_growth_velocity_m_s=np.zeros(incident.shape),
            material_exchange=SurfaceMaterialExchange(
                removed_units_m2={"solid_unit": removed},
                outgoing_units_m2={},
                unresolved_units_m2={"solid_unit": removed},
                deposited_units_m2={"solid_unit": np.zeros(incident.shape)},
                known_limitations=("manufactured STL round-trip gate",)),
            product_populations=(),
            validity=MechanismValidity(
                within_declared_scope=True, reasons=(),
                unsupported_neutral_species=(),
                known_model_form_omissions=("manufactured directional law",),
                parameter_evidence_supports_prediction=False,
                nonpredictive_parameters=("manufactured_normal_velocity",)))


def test_stl_trench_round_trips_through_one_engine_step():
    """End to end: STL void -> level set -> materials -> one advanced step."""
    dx = 0.05
    origin = (-0.325, -0.075, -0.025)
    shape = (13, 3, 15)
    trench = _box_mesh((-0.1, 0.1), (-1.0, 1.0), (0.2, 5.0))
    geometry, report = build_feature_geometry_from_stl(
        trench, dx=dx, mesh_length_unit_m=1e-6, shape=shape, origin=origin,
        solid_region="exterior", solid_ceiling=0.6,
        layers=[(-10.0, 0.45, 1), (0.45, 10.0, 2)])

    assert geometry.phi.shape == shape
    assert np.any(geometry.phi > 0.0) and np.any(geometry.phi < 0.0)
    assert set(np.unique(geometry.material_id)) == {0, 1, 2}
    assert report.solid_region == "exterior"

    # Mesh coordinates run index * dx from zero; the STL-frame origin lives in
    # geometry.mesh_origin_m, so the source plane in SI is origin + extent.
    source_z = (shape[2] - 1) * dx
    reference_plane_m = geometry.mesh_origin_m[2] + source_z * 1.0e-6
    boundary = PlasmaBoundaryState(
        (SpeciesBoundaryState("etch+", 1, 40.0, 2.2e21, [[0.0, 0.0, 10.0]], [1.0]),),
        reference_plane_m=reference_plane_m)
    result = advance_feature_step_3d(
        geometry, boundary, {"etch+": "energetic_bombardment"},
        # Only the lower material etches; the mask above 0.45 stays pinned,
        # which is the point of importing stacked materials from one STL.
        _ManufacturedDirectionalEtch(), etchable_material_ids=(1,),
        duration_s=0.05,
        # Padded a hair beyond the domain, as the engine's own gates do: a
        # source rectangle flush with the mesh edge trips the projected-area
        # conservation guard on float32 mesh vertices.
        source_bounds=(-0.1 * dx, (shape[0] - 1 + 0.1) * dx,
                       -0.1 * dx, (shape[1] - 1 + 0.1) * dx),
        source_z=source_z,
        n_position=4, seed=11, cfl_number=0.25, transport_device="cpu",
        ballistic_transport="face_gather", ballistic_face_quadrature_points=1,
        reinitialization_method="cr2",
        topology_change_policy="continue_gas_cavity")

    assert result.geometry.phi.shape == shape
    assert np.any(result.geometry.phi > 0.0) and np.any(result.geometry.phi < 0.0)
