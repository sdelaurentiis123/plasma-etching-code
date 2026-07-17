import numpy as np
import pytest

from petch.surface_mesh_3d import TriangleSurface3D


def _single_triangle(*, material=1, periodic_lengths=(None, None, None)):
    return TriangleSurface3D(
        vertices=np.asarray([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]),
        faces=np.asarray([[0, 1, 2]]),
        face_material_id=np.asarray([material]),
        periodic_lengths=periodic_lengths,
    )


@pytest.mark.parametrize(
    ("point", "expected_distance", "expected_closest"),
    [
        ([0.25, 0.25, 2.0], 2.0, [0.25, 0.25, 0.0]),
        ([0.75, 0.75, 0.0], np.sqrt(0.125), [0.5, 0.5, 0.0]),
        ([-1.0, -1.0, 0.0], np.sqrt(2.0), [0.0, 0.0, 0.0]),
    ],
)
def test_exact_distance_handles_plane_edge_and_vertex(point, expected_distance, expected_closest):
    result = _single_triangle().nearest_brute_force(point)
    np.testing.assert_allclose(result.distance, [expected_distance], rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(result.closest_point, [expected_closest], rtol=0.0, atol=1e-14)
    np.testing.assert_array_equal(result.face_index, [0])


def test_periodic_nearest_image_is_seam_equivalent():
    surface = TriangleSurface3D(
        vertices=np.asarray([
            [0.05, 0.0, 0.0],
            [0.05, 1.0, 0.0],
            [0.05, 0.0, 1.0],
        ]),
        faces=np.asarray([[0, 1, 2]]),
        face_material_id=np.asarray([1]),
        periodic_lengths=(1.0, None, None),
    )
    result = surface.nearest_brute_force(np.asarray([
        [0.98, 0.2, 0.2],
        [-0.02, 0.2, 0.2],
        [1.98, 0.2, 0.2],
    ]))
    np.testing.assert_allclose(result.distance, 0.07, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(result.periodic_shift, [[1.0, 0.0, 0.0]] * 3)
    np.testing.assert_allclose(result.closest_point[:, 0], 1.05)


def test_material_restricted_query_never_falls_back_to_nearer_other_material():
    surface = TriangleSurface3D(
        vertices=np.asarray([
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [10.0, 0.0, 0.0], [10.0, 1.0, 0.0], [10.0, 0.0, 1.0],
        ]),
        faces=np.asarray([[0, 1, 2], [3, 4, 5]]),
        face_material_id=np.asarray([1, 2]),
    )
    unrestricted = surface.nearest_brute_force([0.1, 0.2, 0.2])
    restricted = surface.nearest_brute_force(
        [0.1, 0.2, 0.2], material_id=2)
    np.testing.assert_array_equal(unrestricted.face_index, [0])
    np.testing.assert_allclose(unrestricted.distance, [0.1])
    np.testing.assert_array_equal(restricted.face_index, [1])
    np.testing.assert_allclose(restricted.distance, [9.9])
    with pytest.raises(ValueError, match="no material 3"):
        surface.nearest([0.1, 0.2, 0.2], material_id=3)


def test_fingerprint_is_layout_independent_order_sensitive_and_ties_use_face_order():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0],
    ])
    faces = np.asarray([[0, 1, 2], [3, 4, 5]])
    material = np.asarray([1, 2])
    first = TriangleSurface3D(vertices, faces, material)
    layout_copy = TriangleSurface3D(
        np.asfortranarray(vertices.astype(np.float32)),
        np.asfortranarray(faces.astype(np.int32)),
        material.astype(np.int32),
    )
    reordered = TriangleSurface3D(vertices, faces[::-1], material[::-1])
    assert first.fingerprint == layout_copy.fingerprint
    assert first.fingerprint != reordered.fingerprint
    result = first.nearest_brute_force([0.2, 0.2, 1.0])
    indexed = first.nearest([0.2, 0.2, 1.0])
    np.testing.assert_array_equal(result.face_index, [0])
    np.testing.assert_array_equal(indexed.face_index, result.face_index)
    np.testing.assert_array_equal(indexed.periodic_shift, result.periodic_shift)


def test_surface_and_result_arrays_are_immutable_and_areas_are_exact():
    surface = _single_triangle()
    np.testing.assert_allclose(surface.face_area, [0.5], rtol=0.0, atol=0.0)
    for array in (
            surface.vertices, surface.faces, surface.face_material_id,
            surface.triangles, surface.face_area, surface.face_centroid,
            surface.face_radius):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.flat[0] = 4
    result = surface.nearest([0.2, 0.2, 0.5], maximum_distance=1.0)
    for array in (
            result.distance, result.face_index, result.closest_point,
            result.periodic_shift, result.candidate_count, result.found):
        assert not array.flags.writeable


@pytest.mark.parametrize(
    "vertices",
    [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    ],
)
def test_degenerate_triangles_are_refused(vertices):
    with pytest.raises(ValueError, match="nondegenerate"):
        TriangleSurface3D(
            vertices=np.asarray(vertices),
            faces=np.asarray([[0, 1, 2]]),
            face_material_id=np.asarray([1]),
        )


def test_bounded_candidates_match_brute_force_on_deterministic_manufactured_mesh():
    vertices = []
    faces = []
    material = []
    for ix in range(4):
        for iy in range(3):
            base = len(vertices)
            vertices.extend((
                (2.0 * ix, 2.0 * iy, 0.0),
                (2.0 * ix + 0.8, 2.0 * iy, 0.0),
                (2.0 * ix, 2.0 * iy + 0.8, 0.0),
            ))
            faces.append((base, base + 1, base + 2))
            material.append(1 + (ix + iy) % 2)
    surface = TriangleSurface3D(
        np.asarray(vertices), np.asarray(faces), np.asarray(material))
    points = np.asarray([
        [2.0 * ix + 0.2, 2.0 * iy + 0.2, 0.25]
        for ix in range(4) for iy in range(3)
    ])
    bounded = surface.nearest(points, maximum_distance=0.5)
    brute = surface.nearest_brute_force(points, maximum_distance=0.5)
    np.testing.assert_allclose(bounded.distance, brute.distance, rtol=0.0, atol=1e-15)
    np.testing.assert_array_equal(bounded.face_index, brute.face_index)
    np.testing.assert_allclose(bounded.closest_point, brute.closest_point, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(bounded.periodic_shift, brute.periodic_shift, rtol=0.0, atol=0.0)
    assert np.all(bounded.candidate_count < brute.candidate_count)
    assert np.all(bounded.found)


def test_bounded_query_reports_certified_absence_without_crossing_the_bound():
    result = _single_triangle().nearest([0.2, 0.2, 2.0], maximum_distance=0.25)
    np.testing.assert_array_equal(result.face_index, [-1])
    assert np.isinf(result.distance[0])
    assert np.all(np.isnan(result.closest_point[0]))
    assert not result.found[0]


def test_indexed_periodic_queries_match_brute_force_face_image_and_tie_authority():
    vertices = []
    faces = []
    material = []
    for ix in range(7):
        for iy in range(5):
            x = 0.08 + 0.13 * ix
            y = 0.06 + 0.18 * iy
            z = 0.03 * np.sin(ix + 2.0 * iy)
            base = len(vertices)
            vertices.extend((
                (x, y, z),
                (x + 0.07, y, z + 0.01),
                (x, y + 0.09, z - 0.015),
            ))
            faces.append((base, base + 1, base + 2))
            material.append(1 + (ix + iy) % 2)
    surface = TriangleSurface3D(
        np.asarray(vertices), np.asarray(faces), np.asarray(material),
        periodic_lengths=(1.0, 1.0, None))
    generator = np.random.default_rng(47382)
    points = generator.uniform(
        [-1.4, -1.2, -0.2], [2.6, 2.4, 0.2], size=(160, 3))
    for material_id in (None, 1, 2):
        for maximum_distance in (None, 0.045, 0.18):
            indexed = surface.nearest(
                points, material_id=material_id,
                maximum_distance=maximum_distance)
            reference = surface.nearest_brute_force(
                points, material_id=material_id,
                maximum_distance=maximum_distance)
            np.testing.assert_array_equal(indexed.face_index, reference.face_index)
            np.testing.assert_array_equal(indexed.periodic_shift, reference.periodic_shift)
            np.testing.assert_allclose(
                indexed.distance, reference.distance, rtol=0.0, atol=2e-15)
            np.testing.assert_allclose(
                indexed.closest_point, reference.closest_point,
                rtol=0.0, atol=2e-15)
            assert np.all(indexed.candidate_count <= reference.candidate_count)

    # The point is exactly halfway between the -1 and primary images.  Image
    # authority is the lexicographically smaller shift, not tree traversal.
    tie_surface = TriangleSurface3D(
        np.asarray([
            [0.5, 0.0, 0.0], [0.5, 0.2, 0.0], [0.5, 0.0, 0.2],
        ]),
        np.asarray([[0, 1, 2]]), np.asarray([1]),
        periodic_lengths=(1.0, None, None))
    tie = tie_surface.nearest([0.0, 0.05, 0.05])
    brute_tie = tie_surface.nearest_brute_force([0.0, 0.05, 0.05])
    np.testing.assert_array_equal(tie.face_index, brute_tie.face_index)
    np.testing.assert_array_equal(tie.periodic_shift, [[-1.0, 0.0, 0.0]])


def test_indexed_periodic_centroid_neighbors_match_exhaustive_physical_reduction():
    vertices = []
    faces = []
    for face in range(12):
        x = 0.04 + 0.075 * face
        base = len(vertices)
        vertices.extend((
            (x, 0.1, 0.0), (x + 0.03, 0.1, 0.0),
            (x, 0.14, 0.01),
        ))
        faces.append((base, base + 1, base + 2))
    surface = TriangleSurface3D(
        np.asarray(vertices), np.asarray(faces), np.ones(12, dtype=int),
        periodic_lengths=(1.0, None, None))
    points = np.asarray([
        [0.0, 0.12, 0.0], [0.5, 0.13, 0.02], [0.99, 0.11, -0.01],
    ])
    result = surface.nearest_face_centroids(points, count=5, material_id=1)
    shifts = surface._periodic_shifts
    for row, point in enumerate(points):
        image = surface.face_centroid[:, None, :] + shifts[None, :, :]
        distance = np.linalg.norm(image - point[None, None, :], axis=2)
        image_choice = np.argmin(distance, axis=1)
        physical_distance = distance[np.arange(len(surface.faces)), image_choice]
        order = np.lexsort((np.arange(len(surface.faces)), physical_distance))[:5]
        np.testing.assert_array_equal(result.face_index[row], order)
        np.testing.assert_allclose(
            result.distance[row], physical_distance[order], rtol=0.0, atol=2e-15)
        np.testing.assert_array_equal(
            result.periodic_shift[row], shifts[image_choice[order]])


def test_certified_image_candidate_rows_match_exhaustive_sphere_bound():
    vertices = []
    faces = []
    for face in range(9):
        x = 0.05 + 0.1 * face
        base = len(vertices)
        vertices.extend((
            (x, 0.1, 0.0), (x + 0.04, 0.1, 0.0),
            (x, 0.16, 0.01),
        ))
        faces.append((base, base + 1, base + 2))
    surface = TriangleSurface3D(
        np.asarray(vertices), np.asarray(faces), np.ones(9, dtype=int),
        periodic_lengths=(1.0, None, None))
    points = np.asarray([[0.01, 0.13, 0.0], [0.52, 0.12, 0.02]])
    query_radius = np.asarray([0.12, 0.09])
    result = surface.certified_triangle_image_candidates(
        points, query_radius=query_radius, material_id=1)
    image_centroid = (
        surface.face_centroid[:, None, :] + surface._periodic_shifts[None, :, :]
    ).reshape(-1, 3)
    image_face = np.repeat(np.arange(len(surface.faces)), len(surface._periodic_shifts))
    image_shift = np.tile(surface._periodic_shifts, (len(surface.faces), 1))
    image_radius = np.repeat(surface.face_radius, len(surface._periodic_shifts))
    for row, point in enumerate(points):
        distance = np.linalg.norm(image_centroid - point[None, :], axis=1)
        expected = np.flatnonzero(
            distance <= query_radius[row] + image_radius + 1e-14)
        order = np.lexsort((
            image_shift[expected, 2], image_shift[expected, 1],
            image_shift[expected, 0], image_face[expected]))
        expected = expected[order]
        start, stop = result.row_offsets[row:row + 2]
        np.testing.assert_array_equal(
            result.face_index[start:stop], image_face[expected])
        np.testing.assert_array_equal(
            result.periodic_shift[start:stop], image_shift[expected])
        np.testing.assert_allclose(
            result.centroid_distance[start:stop], distance[expected],
            rtol=0.0, atol=2e-15)


def test_periodic_float32_endpoint_roundoff_is_canonical_but_real_excursion_refuses():
    epsilon = 5.0e-10
    rounded = TriangleSurface3D(
        np.asarray([
            [-epsilon, 0.0, 0.0], [1.0 + epsilon, 0.0, 0.0],
            [0.0, 0.5, 0.0],
        ]),
        np.asarray([[0, 1, 2]]), np.asarray([1]),
        periodic_lengths=(1.0, None, None))
    np.testing.assert_array_equal(rounded.vertices[:, 0], [0.0, 1.0, 0.0])

    with pytest.raises(ValueError, match="primary cell"):
        TriangleSurface3D(
            np.asarray([
                [-1e-3, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.5, 0.0],
            ]),
            np.asarray([[0, 1, 2]]), np.asarray([1]),
            periodic_lengths=(1.0, None, None))
