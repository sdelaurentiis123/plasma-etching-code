import subprocess
import sys

import numpy as np
import pytest

from petch.feature_geometry_backend_3d import (
    FeatureGeometryBackend3D,
    UniformFeatureGeometryBackend3D,
)
from petch.feature_geometry_state_3d import (
    FeatureGeometry3D as StateFeatureGeometry3D,
    face_material_ids_3d,
)
from petch.feature_step_3d import (
    FeatureGeometry3D,
    _extract_uniform_surface_arrays,
    _face_material_ids,
    make_rectangular_trench_geometry_3d,
)
from petch.threed import (
    _cancel_duplicate_marching_cubes_faces_3d,
    _certify_extracted_surface_topology_3d,
    extract_mesh_3d,
)


def test_feature_step_reexports_exact_dependency_neutral_authorities():
    assert FeatureGeometry3D is StateFeatureGeometry3D
    assert _face_material_ids is face_material_ids_3d


@pytest.mark.parametrize("module_order", [
    ("petch.feature_geometry_backend_3d", "petch.feature_step_3d"),
    ("petch.feature_step_3d", "petch.feature_geometry_backend_3d"),
    ("petch", "petch.feature_geometry_backend_3d", "petch.feature_step_3d"),
])
def test_geometry_modules_are_import_order_independent(module_order):
    source_root = str(__file__).rsplit("/tests/", 1)[0] + "/src"
    code = (
        "import importlib, sys\n"
        f"sys.path.insert(0, {source_root!r})\n"
        f"order = {module_order!r}\n"
        "for name in order: importlib.import_module(name)\n"
        "from petch.feature_geometry_state_3d import "
        "FeatureGeometry3D as state_geometry, face_material_ids_3d\n"
        "from petch.feature_step_3d import "
        "FeatureGeometry3D as step_geometry, _face_material_ids\n"
        "assert state_geometry is step_geometry\n"
        "assert face_material_ids_3d is _face_material_ids\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr


def _plane_geometry(*, origin=(0.0, 0.0, 0.0)):
    dx = 0.25
    shape = (5, 4, 7)
    x, y, z = (np.arange(size) * dx for size in shape)
    _, _, Z = np.meshgrid(x, y, z, indexing="ij")
    phi = 0.625 - Z
    material = np.where(phi >= 0.0, 1, 0)
    return FeatureGeometry3D(
        phi, material, dx, 1e-6, origin,
        material_levelsets={1: phi})


def _trench_geometry():
    base = make_rectangular_trench_geometry_3d(
        cell_width=0.04,
        cell_length=0.02,
        domain_height=0.10,
        dx=0.005,
        opening_width=0.02,
        mask_thickness=0.02,
        substrate_top=0.05,
        etched_depth=0.01,
    )
    return FeatureGeometry3D(
        base.phi, base.material_id, base.dx, base.mesh_length_unit_m,
        (1e-6, 2e-6, 3e-6), material_levelsets=base.material_levelsets)


def test_uniform_backend_reports_physical_extent_and_periodic_metadata():
    geometry = _trench_geometry()
    backend = UniformFeatureGeometryBackend3D(geometry, periodic_axes=(1, 0))
    expected_mesh = tuple((np.asarray(geometry.phi.shape) - 1) * geometry.dx)
    expected_m = tuple(value * geometry.mesh_length_unit_m for value in expected_mesh)

    assert isinstance(backend, FeatureGeometryBackend3D)
    assert backend.backend_kind == "uniform_dense_reference_v1"
    assert backend.shape == geometry.phi.shape
    assert backend.finest_spacing_mesh_units == geometry.dx
    assert backend.finest_spacing_m == geometry.dx * geometry.mesh_length_unit_m
    assert np.allclose(backend.domain_extent_mesh_units, expected_mesh)
    assert np.allclose(backend.domain_extent_m, expected_m)
    assert backend.domain_bounds_m[0] == geometry.mesh_origin_m
    assert np.allclose(
        backend.domain_bounds_m[1], np.asarray(geometry.mesh_origin_m) + expected_m)
    assert backend.periodic_metadata.axes == (0, 1)
    assert backend.periodic_metadata.lengths_mesh_units == (
        expected_mesh[0], expected_mesh[1], None)
    assert backend.periodic_metadata.lengths_m == (expected_m[0], expected_m[1], None)
    assert backend.periodic_metadata.duplicate_endpoint_planes == (True, True, False)


def test_uniform_backend_samples_mesh_and_si_signed_distance():
    geometry = _plane_geometry(origin=(2e-6, -1e-6, 4e-6))
    backend = UniformFeatureGeometryBackend3D(geometry)
    points_mesh = np.array([
        [0.13, 0.42, 0.10],
        [0.77, 0.51, 0.625],
        [0.40, 0.20, 1.20],
    ])
    expected_mesh = 0.625 - points_mesh[:, 2]
    points_m = (np.asarray(geometry.mesh_origin_m)[None, :]
                + points_mesh * geometry.mesh_length_unit_m)

    assert np.allclose(
        backend.sample_signed_distance_mesh(points_mesh), expected_mesh, atol=2e-16)
    assert np.allclose(
        backend.sample_signed_distance_mesh(points_mesh, material_id=1),
        expected_mesh, atol=2e-16)
    assert np.allclose(
        backend.sample_signed_distance_m(points_m),
        expected_mesh * geometry.mesh_length_unit_m, atol=2e-22)


def test_uniform_backend_material_owner_is_material_local_and_gas_is_zero():
    geometry = _trench_geometry()
    backend = UniformFeatureGeometryBackend3D(geometry, periodic_axes=(0, 1))
    points_mesh = np.array([
        [0.020, 0.010, 0.030],  # substrate below the etched floor
        [0.005, 0.010, 0.060],  # mask outside the opening
        [0.020, 0.010, 0.060],  # gas inside the opening
    ])
    points_m = (np.asarray(geometry.mesh_origin_m)[None, :]
                + points_mesh * geometry.mesh_length_unit_m)

    assert np.array_equal(
        backend.sample_material_owner_mesh(points_mesh), np.array([1, 2, 0]))
    assert np.array_equal(
        backend.sample_material_owner_m(points_m), np.array([1, 2, 0]))


def test_uniform_backend_surface_extraction_is_exact_reference_parity():
    geometry = _trench_geometry()
    backend = UniformFeatureGeometryBackend3D(geometry, periodic_axes=(0, 1))
    expected = extract_mesh_3d(geometry.phi, geometry.dx)
    surface = backend.extract_surface()

    assert np.array_equal(surface.vertices_mesh, expected[0])
    assert np.array_equal(surface.faces, expected[1])
    assert np.array_equal(surface.centroids_mesh, expected[2])
    assert np.array_equal(surface.areas_mesh2, expected[3])
    assert np.array_equal(
        surface.face_material_id, _face_material_ids(expected[2], geometry))
    assert np.allclose(
        surface.vertices_m,
        np.asarray(geometry.mesh_origin_m)[None, :]
        + expected[0] * geometry.mesh_length_unit_m)
    assert np.allclose(
        surface.areas_m2, expected[3] * geometry.mesh_length_unit_m ** 2)
    assert not surface.vertices_mesh.flags.writeable
    assert not surface.faces.flags.writeable


def test_extract_mesh_retains_positive_slivers_needed_for_topological_closure():
    # One barely positive interior node creates a closed subcell octahedron.  Every face is much
    # smaller than the former float32-scaled area floor, but all eight faces are representable and
    # topologically necessary.  Removing them converts a closed scalar-field contour into a hole.
    phi = -np.ones((3, 3, 3))
    phi[1, 1, 1] = 1.0e-3

    vertices, faces, _centroids, areas = extract_mesh_3d(phi, 1.0)
    edge = np.sort(
        faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    _unique, count = np.unique(edge, axis=0, return_counts=True)
    former_floor = 32.0 * np.finfo(np.float32).eps

    assert len(faces) == 8
    assert np.all((areas > 0.0) & (areas < former_floor))
    assert np.array_equal(count, np.full_like(count, 2))
    _certify_extracted_surface_topology_3d(vertices, faces, (2.0, 2.0, 2.0))


def test_extracted_surface_topology_refuses_an_interior_hole():
    vertices = np.asarray([
        [0.8, 0.8, 0.8],
        [1.2, 0.8, 0.8],
        [1.0, 1.2, 0.8],
        [1.0, 1.0, 1.2],
    ])
    # A closed tetrahedron with one missing face has three unmatched interior edges.
    faces = np.asarray([[0, 1, 2], [0, 3, 1], [1, 3, 2]])

    with pytest.raises(RuntimeError, match="3 unmatched interior edges"):
        _certify_extracted_surface_topology_3d(
            vertices, faces, (2.0, 2.0, 2.0))


def test_marching_cubes_duplicate_faces_reduce_as_an_oriented_surface_chain():
    vertices = np.asarray([
        [0.5, 0.5, 0.5],
        [1.0, 0.5, 0.5],
        [0.5, 1.0, 0.5],
        [0.5, 0.5, 1.0],
    ])
    faces = np.asarray([
        [0, 1, 2],
        [0, 2, 1],  # opposite-winding duplicate pair: both cancel
        [0, 1, 3],
        [0, 1, 3],  # same-winding duplicate pair: one remains
    ])

    reduced = _cancel_duplicate_marching_cubes_faces_3d(vertices, faces)

    assert np.array_equal(reduced, np.asarray([[0, 1, 3]]))


def test_uniform_surface_bridge_preserves_exact_legacy_writable_array_contract():
    geometry = _trench_geometry()
    expected = (*extract_mesh_3d(geometry.phi, geometry.dx),)
    expected = (*expected, _face_material_ids(expected[2], geometry))
    observed = _extract_uniform_surface_arrays(geometry)

    assert len(observed) == len(expected) == 5
    for actual, reference in zip(observed, expected):
        assert np.array_equal(actual, reference)
        assert actual.dtype == reference.dtype
        assert actual.shape == reference.shape
        assert actual.flags.c_contiguous
        assert actual.flags.writeable


def test_face_material_authority_preserves_layered_and_fallback_owner_paths():
    layered = _trench_geometry()
    layered_centroids = extract_mesh_3d(layered.phi, layered.dx)[2]
    layered_owner = face_material_ids_3d(layered_centroids, layered)
    assert np.array_equal(
        layered_owner, _face_material_ids(layered_centroids, layered))
    assert set(np.unique(layered_owner)) == {1, 2}

    plane = _plane_geometry()
    fallback = FeatureGeometry3D(
        plane.phi, plane.material_id, plane.dx, plane.mesh_length_unit_m,
        plane.mesh_origin_m, material_levelsets=None)
    fallback_centroids = extract_mesh_3d(fallback.phi, fallback.dx)[2]
    fallback_owner = face_material_ids_3d(fallback_centroids, fallback)
    assert np.array_equal(
        fallback_owner, _face_material_ids(fallback_centroids, fallback))
    assert np.array_equal(fallback_owner, np.ones(len(fallback_centroids), dtype=int))


def test_periodic_sampling_wraps_while_nonperiodic_sampling_refuses():
    geometry = _plane_geometry()
    periodic = UniformFeatureGeometryBackend3D(geometry, periodic_axes=(0, 1))
    extent = np.asarray(periodic.domain_extent_mesh_units)
    reference = np.array([[0.10, 0.20, 0.50]])
    shifted = reference + np.array([[extent[0], -2.0 * extent[1], 0.0]])

    assert np.array_equal(
        periodic.sample_signed_distance_mesh(reference),
        periodic.sample_signed_distance_mesh(shifted))
    with pytest.raises(ValueError, match="outside the domain"):
        UniformFeatureGeometryBackend3D(geometry).sample_signed_distance_mesh(shifted)


def test_backend_and_surface_fingerprints_are_deterministic_and_sensitive():
    geometry = _plane_geometry()
    copied = FeatureGeometry3D(
        geometry.phi.copy(), geometry.material_id.copy(), geometry.dx,
        geometry.mesh_length_unit_m, geometry.mesh_origin_m,
        material_levelsets={1: geometry.material_levelsets[1].copy()})
    first = UniformFeatureGeometryBackend3D(geometry, periodic_axes=(0, 1))
    second = UniformFeatureGeometryBackend3D(copied, periodic_axes=(1, 0))

    assert first.fingerprint == first.fingerprint
    assert first.fingerprint == second.fingerprint
    assert first.extract_surface().fingerprint == second.extract_surface().fingerprint
    assert first.fingerprint != UniformFeatureGeometryBackend3D(
        copied, periodic_axes=(0,)).fingerprint

    shifted_origin = FeatureGeometry3D(
        copied.phi, copied.material_id, copied.dx, copied.mesh_length_unit_m,
        (1e-12, 0.0, 0.0), material_levelsets=copied.material_levelsets)
    assert first.fingerprint != UniformFeatureGeometryBackend3D(
        shifted_origin, periodic_axes=(0, 1)).fingerprint


@pytest.mark.parametrize("axes", [(True,), (0, 0), (3,), (-1,)])
def test_uniform_backend_refuses_invalid_periodic_axes(axes):
    with pytest.raises(ValueError, match="periodic axes"):
        UniformFeatureGeometryBackend3D(_plane_geometry(), periodic_axes=axes)


def test_extract_mesh_survives_grid_aligned_near_zero_plateau():
    """Production w80 moving-Cr failure, step 80: a shielded, grid-aligned
    TiO2 top leaves float-noise phi (1e-19..3e-7) at mesh nodes under the
    retreating Cr cap corner; marching cubes then emitted degenerate
    vertex-on-node triangles and the watertight certification failed with
    "3 unmatched interior edges".  The vertex-on-surface guard must extract
    this exact captured field cleanly."""
    from pathlib import Path
    import numpy as np
    from petch.threed import extract_mesh_3d

    data = np.load(Path(__file__).parent / "data"
                   / "w80_step80_grid_aligned_zero_phi.npz")
    phi = np.asarray(data["phi"], dtype=float)
    dx = float(data["dx"])
    assert int(np.sum(np.abs(phi) < 1.0e-4 * dx)) >= 30  # the noise plateau
    verts, faces, centroids, areas = extract_mesh_3d(phi, dx)
    assert len(faces) > 1000
    assert np.all(areas > 0.0)
