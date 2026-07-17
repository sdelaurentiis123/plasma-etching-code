import numpy as np
import pytest

from petch.surface_mesh_3d import TriangleSurface3D
from petch.surface_partitioned_overlap_3d import (
    build_partitioned_surface_overlap_transfer_3d,
)


def _patch(axis, coordinate, first_range, second_range, *, sign, material):
    tangent = [item for item in range(3) if item != axis]
    vertex = np.zeros((4, 3), dtype=float)
    vertex[:, axis] = coordinate
    vertex[:, tangent[0]] = [first_range[0], first_range[1], first_range[0], first_range[1]]
    vertex[:, tangent[1]] = [
        second_range[0], second_range[0], second_range[1], second_range[1]]
    face = np.asarray([[0, 1, 3], [0, 3, 2]], dtype=int)
    triangle = vertex[face[0]]
    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
    if np.sign(normal[axis]) != int(sign):
        face = face[:, [0, 2, 1]]
    return vertex, face, np.full(2, int(material), dtype=int)


def _surface(patches, *, periodic_lengths=(None, None, None)):
    vertices = []
    faces = []
    material = []
    for vertex, face, owner in patches:
        offset = len(vertices)
        vertices.extend(vertex)
        faces.extend(face + offset)
        material.extend(owner)
    return TriangleSurface3D(
        np.asarray(vertices), np.asarray(faces), np.asarray(material),
        periodic_lengths=periodic_lengths)


def _three_sided_trench(*, floor=0.0, left=0.0, right=1.0):
    return _surface([
        _patch(2, floor, (0.0, 1.0), (0.0, 1.0), sign=1, material=1),
        _patch(0, left, (0.0, 1.0), (0.0, 1.0), sign=1, material=2),
        _patch(0, right, (0.0, 1.0), (0.0, 1.0), sign=-1, material=2),
    ], periodic_lengths=(None, 1.0, None))


def test_floor_and_two_walls_transfer_as_one_conservative_operator():
    old = _three_sided_trench()
    new = _three_sided_trench(floor=0.01, left=0.01, right=0.99)
    transfer = build_partitioned_surface_overlap_transfer_3d(
        old, new, maximum_normal_distance=0.02)

    assert transfer.old_patch_count == 3
    assert transfer.new_patch_count == 3
    assert transfer.positive_patch_pair_count == 3
    np.testing.assert_allclose(
        transfer.old_covered_area, old.face_area, rtol=0.0, atol=3e-15)
    np.testing.assert_allclose(
        transfer.new_covered_area, new.face_area, rtol=0.0, atol=3e-15)
    extensive = transfer.apply_extensive(
        np.where(old.face_material_id == 1, 2.0, 7.0),
        newly_exposed_density=0.0)
    np.testing.assert_allclose(
        extensive.values[new.face_material_id == 1], 2.0, atol=3e-14)
    np.testing.assert_allclose(
        extensive.values[new.face_material_id == 2], 7.0, atol=3e-14)
    assert extensive.maximum_relative_balance_error < 2e-14
    assert transfer.geometry_receipt["row_column_area_residual"] == 0.0


def test_piecewise_planar_partial_area_uses_explicit_removed_and_exposed_closures():
    old = _surface([
        _patch(2, 0.0, (0.0, 1.0), (0.0, 1.0), sign=1, material=1),
        _patch(0, 0.0, (0.0, 1.0), (0.0, 1.0), sign=1, material=2),
    ])
    new = _surface([
        _patch(2, 0.01, (0.25, 1.25), (0.0, 1.0), sign=1, material=1),
        _patch(0, 0.01, (0.0, 1.0), (0.0, 1.0), sign=1, material=2),
    ])
    transfer = build_partitioned_surface_overlap_transfer_3d(
        old, new, maximum_normal_distance=0.02)
    assert np.sum(transfer.old_uncovered_area) == pytest.approx(0.25, abs=4e-15)
    assert np.sum(transfer.new_uncovered_area) == pytest.approx(0.25, abs=4e-15)
    result = transfer.apply_extensive(
        np.where(old.face_material_id == 1, 2.0, 5.0),
        newly_exposed_density=np.where(new.face_material_id == 1, 3.0, 5.0))
    ledger = result.material_ledger[1]
    assert ledger["removed_inventory"] == pytest.approx(0.5, abs=4e-15)
    assert ledger["newly_exposed_inventory"] == pytest.approx(0.75, abs=4e-15)
    assert ledger["new_inventory"] == pytest.approx(2.25, abs=4e-15)


def test_patch_normal_to_periodic_axis_uses_nearest_cell_image_once():
    old = _surface([
        _patch(0, 0.98, (0.0, 1.0), (0.0, 1.0), sign=1, material=1),
    ], periodic_lengths=(1.0, None, None))
    new = _surface([
        _patch(0, 0.02, (0.0, 1.0), (0.0, 1.0), sign=1, material=1),
    ], periodic_lengths=(1.0, None, None))
    transfer = build_partitioned_surface_overlap_transfer_3d(
        old, new, maximum_normal_distance=0.05)
    assert np.sum(transfer.overlap_area) == pytest.approx(1.0, abs=3e-15)
    assert transfer.patch_receipts[0]["normal_image_shift"] == pytest.approx(1.0)
    np.testing.assert_allclose(transfer.old_covered_area, old.face_area, atol=3e-15)
    np.testing.assert_allclose(transfer.new_covered_area, new.face_area, atol=3e-15)


def test_nearby_rotating_surface_refuses_instead_of_declaring_wholesale_replacement():
    old = _surface([
        _patch(2, 0.0, (0.0, 1.0), (0.0, 1.0), sign=1, material=1),
    ])
    vertices = np.array(old.vertices, copy=True)
    vertices[[1, 2], 2] += 0.01
    rotated = TriangleSurface3D(vertices, old.faces, old.face_material_id)
    with pytest.raises(ValueError, match="nearby nonparallel surface"):
        build_partitioned_surface_overlap_transfer_3d(
            old, rotated, maximum_normal_distance=0.05)
