import numpy as np

from petch.surface_charge_remap_3d import remap_surface_charge_3d


def _plane_mesh(cells, z=0.0):
    coordinate = np.linspace(0.0, 1.0, int(cells) + 1)
    x, y = np.meshgrid(coordinate, coordinate, indexing="ij")
    vertices = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, float(z))))
    faces = []
    stride = int(cells) + 1
    for i in range(int(cells)):
        for j in range(int(cells)):
            lower = i * stride + j
            faces.extend(((lower, lower + stride, lower + stride + 1),
                          (lower, lower + stride + 1, lower + 1)))
    return vertices, np.asarray(faces, dtype=int)


def _centroids(vertices, faces):
    return np.asarray(vertices)[np.asarray(faces)].mean(axis=1)


def _sigma(point):
    point = np.asarray(point)
    return 2.0e-6 + 0.4e-6 * np.sin(2.0 * np.pi * point[:, 0]) * np.cos(
        2.0 * np.pi * point[:, 1])


def test_zero_motion_is_bitwise_identical_and_closes_signed_ledger():
    vertices, faces = _plane_mesh(3)
    sigma = np.linspace(-3.0e-6, 4.0e-6, len(faces))
    result = remap_surface_charge_3d(
        vertices, faces, sigma, np.ones(len(faces), dtype=int), np.zeros(len(faces)),
        vertices.copy(), faces.copy(), np.ones(len(faces), dtype=int),
        mesh_length_unit_m=1.0)

    assert np.array_equal(result.sigma_c_per_m2, sigma)
    assert np.array_equal(result.removed_charge_by_old_face_c, np.zeros(len(faces)))
    assert result.removed_positive_charge_c == 0.0
    assert result.removed_negative_charge_c == 0.0
    assert result.relative_charge_balance_error == 0.0


def test_translating_material_surface_carries_charge_with_exact_conservation():
    old_vertices, old_faces = _plane_mesh(4, z=0.0)
    new_vertices, new_faces = _plane_mesh(4, z=0.125)
    sigma = _sigma(_centroids(old_vertices, old_faces))
    material = np.ones(len(old_faces), dtype=int)
    result = remap_surface_charge_3d(
        old_vertices, old_faces, sigma, material, np.full(len(old_faces), 0.125),
        new_vertices, new_faces, material, mesh_length_unit_m=1.0,
        neighbor_count=1, maximum_distance=0.13)

    assert np.allclose(result.sigma_c_per_m2, sigma, rtol=2e-16, atol=0.0)
    assert result.removed_net_charge_c == 0.0
    assert np.isclose(result.face_charge_c.sum(), sigma.mean(), rtol=2e-16)
    assert result.relative_charge_balance_error < 2e-15


def test_uniform_recession_removes_all_charge_and_itemizes_both_signs():
    old_vertices, old_faces = _plane_mesh(3, z=0.2)
    new_vertices, new_faces = _plane_mesh(3, z=0.1)
    sigma = np.linspace(-2.0e-6, 3.0e-6, len(old_faces))
    material = np.ones(len(old_faces), dtype=int)
    result = remap_surface_charge_3d(
        old_vertices, old_faces, sigma, material, np.full(len(old_faces), -0.1),
        new_vertices, new_faces, material, mesh_length_unit_m=1.0)

    old_charge = sigma / len(old_faces)
    assert np.array_equal(result.sigma_c_per_m2, np.zeros(len(new_faces)))
    assert np.allclose(result.removed_charge_by_old_face_c, old_charge, rtol=3e-16)
    assert np.isclose(result.removed_positive_charge_c,
                      np.maximum(old_charge, 0.0).sum(), rtol=2e-16)
    assert np.isclose(result.removed_negative_charge_c,
                      np.maximum(-old_charge, 0.0).sum(), rtol=2e-16)
    assert result.relative_charge_balance_error < 3e-16


def test_advancing_plane_remap_converges_under_surface_refinement():
    errors = []
    for cells in (4, 8, 16):
        old_vertices, old_faces = _plane_mesh(cells, z=0.0)
        new_vertices, new_faces = _plane_mesh(2 * cells, z=0.1)
        old_sigma = _sigma(_centroids(old_vertices, old_faces))
        result = remap_surface_charge_3d(
            old_vertices, old_faces, old_sigma, np.ones(len(old_faces), dtype=int),
            np.full(len(old_faces), 0.1), new_vertices, new_faces,
            np.ones(len(new_faces), dtype=int), mesh_length_unit_m=1.0,
            neighbor_count=4, maximum_distance=0.3)
        exact = _sigma(_centroids(new_vertices, new_faces))
        errors.append(float(np.sqrt(np.mean((result.sigma_c_per_m2 - exact) ** 2))))
        assert result.relative_charge_balance_error < 2e-15

    order = np.log2(np.asarray(errors[:-1]) / np.asarray(errors[1:]))
    assert np.all(order > 0.8)
