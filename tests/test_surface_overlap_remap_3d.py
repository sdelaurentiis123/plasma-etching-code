import time

import numpy as np
import pytest

from petch.surface_mesh_3d import TriangleSurface3D
from petch.surface_overlap_remap_3d import build_surface_overlap_transfer_3d
from petch.surface_transfer_3d import build_surface_transfer_3d


def _grid_surface(
        x_coordinates, y_coordinates, *, z=0.0, flip_diagonal=False,
        material=None, periodic_lengths=(None, None, None)):
    x = np.asarray(x_coordinates, dtype=float)
    y = np.asarray(y_coordinates, dtype=float)
    vertices = np.asarray([(xx, yy, z) for xx in x for yy in y])
    faces = []
    face_material = []
    for ix in range(len(x) - 1):
        for iy in range(len(y) - 1):
            lower = ix * len(y) + iy
            upper = (ix + 1) * len(y) + iy
            if flip_diagonal:
                faces.extend(((lower, upper, lower + 1),
                              (upper, upper + 1, lower + 1)))
            else:
                faces.extend(((lower, upper, upper + 1),
                              (lower, upper + 1, lower + 1)))
            selected = 1 if material is None else int(material(ix, iy))
            face_material.extend((selected, selected))
    return TriangleSurface3D(
        vertices, np.asarray(faces), np.asarray(face_material),
        periodic_lengths=periodic_lengths)


def _square_center_fan(*, z=0.0, material=1):
    vertices = np.asarray([
        [0.0, 0.0, z], [1.0, 0.0, z], [1.0, 1.0, z],
        [0.0, 1.0, z], [0.5, 0.5, z],
    ])
    faces = np.asarray([
        [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4],
    ])
    return TriangleSurface3D(
        vertices, faces, np.full(4, int(material), dtype=int))


def _build(old, new, *, maximum_normal_distance=0.02):
    return build_surface_overlap_transfer_3d(
        old, new, projection_axis=2, orientation_sign=1,
        maximum_normal_distance=maximum_normal_distance)


def test_noop_is_bitwise_and_receipts_are_immutable_and_replayable():
    surface = _grid_surface([0.0, 1.0], [0.0, 1.0])
    first = _build(surface, surface, maximum_normal_distance=0.0)
    replay = _build(surface, surface, maximum_normal_distance=0.0)
    np.testing.assert_array_equal(first.row_offsets, [0, 1, 2])
    np.testing.assert_array_equal(first.old_face_index, [0, 1])
    np.testing.assert_array_equal(first.overlap_area, surface.face_area)
    assert first.fingerprint == replay.fingerprint
    extensive = first.apply_extensive(
        [2.0, 7.0], newly_exposed_density=0.0)
    intensive = first.apply_intensive([0.2, 0.8])
    explicit_fill = first.apply_intensive([0.2, 0.8], uncovered_fill=0.0)
    np.testing.assert_array_equal(extensive.values, [2.0, 7.0])
    np.testing.assert_array_equal(intensive.values, [0.2, 0.8])
    np.testing.assert_array_equal(explicit_fill.values, intensive.values)
    assert explicit_fill.application_fingerprint != intensive.application_fingerprint
    assert extensive.maximum_relative_balance_error == 0.0
    assert first.geometry_receipt["row_column_area_residual"] == 0.0
    for array in (
            first.row_offsets, first.old_face_index, first.overlap_area,
            first.old_covered_area, first.new_covered_area,
            first.old_uncovered_area, first.new_uncovered_area,
            extensive.values, intensive.values):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.flat[0] = 9
    with pytest.raises(TypeError):
        first.geometry_receipt["operator"] = "changed"
    with pytest.raises(TypeError):
        extensive.material_ledger[1]["new_inventory"] = 0.0


def test_retriangulation_and_subdivision_close_area_and_inventory_ledgers():
    old = _grid_surface([0.0, 1.0], [0.0, 1.0])
    new = _square_center_fan(z=0.01)
    transfer = _build(old, new)
    np.testing.assert_allclose(
        transfer.old_covered_area, old.face_area, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(
        transfer.new_covered_area, new.face_area, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(transfer.old_uncovered_area, 0.0, atol=2e-15)
    np.testing.assert_allclose(transfer.new_uncovered_area, 0.0, atol=2e-15)
    extensive = transfer.apply_extensive(
        [2.0, 4.0], newly_exposed_density=0.0)
    ledger = extensive.material_ledger[1]
    assert abs(ledger["old_balance_residual"]) < 2e-15
    assert abs(ledger["new_balance_residual"]) < 2e-15
    assert ledger["removed_inventory"] == pytest.approx(0.0, abs=2e-15)
    assert ledger["newly_exposed_inventory"] == pytest.approx(0.0, abs=2e-15)
    constant = transfer.apply_intensive([0.375, 0.375])
    np.testing.assert_allclose(constant.values, 0.375, rtol=0.0, atol=2e-15)


def test_partial_overlap_itemizes_removed_and_newly_exposed_closures():
    old = _grid_surface([0.0, 1.0], [0.0, 1.0])
    new = _grid_surface([0.25, 1.25], [0.0, 1.0], z=0.01,
                        flip_diagonal=True)
    transfer = _build(old, new)
    assert np.sum(transfer.overlap_area) == pytest.approx(0.75, abs=3e-15)
    assert np.sum(transfer.old_uncovered_area) == pytest.approx(0.25, abs=3e-15)
    assert np.sum(transfer.new_uncovered_area) == pytest.approx(0.25, abs=3e-15)
    with pytest.raises(TypeError):
        transfer.apply_extensive([2.0, 2.0])
    extensive = transfer.apply_extensive(
        [2.0, 2.0], newly_exposed_density=3.0)
    ledger = extensive.material_ledger[1]
    assert ledger["old_inventory"] == pytest.approx(2.0, abs=3e-15)
    assert ledger["retained_inventory"] == pytest.approx(1.5, abs=3e-15)
    assert ledger["removed_inventory"] == pytest.approx(0.5, abs=3e-15)
    assert ledger["newly_exposed_inventory"] == pytest.approx(0.75, abs=3e-15)
    assert ledger["new_inventory"] == pytest.approx(2.25, abs=3e-15)
    assert extensive.maximum_relative_balance_error < 2e-15
    with pytest.raises(ValueError, match="requires uncovered_fill"):
        transfer.apply_intensive([1.0, 2.0])
    intensive = transfer.apply_intensive([1.0, 2.0], uncovered_fill=4.0)
    assert np.all(intensive.values >= 1.0)
    assert np.all(intensive.values <= 4.0)


def test_shifted_periodic_retriangulation_counts_fundamental_area_once():
    old = _grid_surface(
        [0.0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.5, 1.0],
        periodic_lengths=(1.0, 1.0, None))
    new = _grid_surface(
        [0.0, 0.125, 0.375, 0.625, 0.875, 1.0],
        [0.0, 0.25, 0.75, 1.0], z=0.005, flip_diagonal=True,
        periodic_lengths=(1.0, 1.0, None))
    transfer = _build(old, new)
    assert np.sum(transfer.overlap_area) == pytest.approx(1.0, abs=8e-15)
    np.testing.assert_allclose(
        transfer.old_covered_area, old.face_area, rtol=0.0, atol=5e-15)
    np.testing.assert_allclose(
        transfer.new_covered_area, new.face_area, rtol=0.0, atol=5e-15)
    assert np.all(transfer.old_covered_area <= old.face_area + 5e-15)
    assert np.all(transfer.new_covered_area <= new.face_area + 5e-15)
    constant = transfer.apply_extensive(
        np.full(len(old.faces), 2.5), newly_exposed_density=0.0)
    np.testing.assert_allclose(constant.values, 2.5, rtol=0.0, atol=2e-13)
    assert constant.maximum_relative_balance_error < 2e-14


def test_mixed_material_overlap_never_crosses_junction():
    material = lambda ix, iy: 1 if ix == 0 else 2
    old = _grid_surface([0.0, 0.5, 1.0], [0.0, 0.5, 1.0], material=material)
    new = _grid_surface(
        [0.0, 0.5, 1.0], [0.0, 0.25, 0.75, 1.0], z=0.005,
        flip_diagonal=True, material=material)
    transfer = _build(old, new)
    row = transfer.row_face_index
    np.testing.assert_array_equal(
        old.face_material_id[transfer.old_face_index],
        new.face_material_id[row])
    values = np.where(old.face_material_id == 1, 3.0, 9.0)
    applied = transfer.apply_intensive(values)
    np.testing.assert_allclose(
        applied.values[new.face_material_id == 1], 3.0, atol=3e-14)
    np.testing.assert_allclose(
        applied.values[new.face_material_id == 2], 9.0, atol=3e-14)


def test_linear_intensive_field_converges_under_opposite_diagonal_retriangulation():
    errors = []
    for cell_count in (4, 8, 16):
        coordinate = np.linspace(0.0, 1.0, cell_count + 1)
        old = _grid_surface(coordinate, coordinate)
        new = _grid_surface(
            coordinate, coordinate, z=0.002, flip_diagonal=True)
        transfer = _build(old, new)
        old_value = (
            0.3 + old.face_centroid[:, 0]
            + 2.0 * old.face_centroid[:, 1])
        exact = (
            0.3 + new.face_centroid[:, 0]
            + 2.0 * new.face_centroid[:, 1])
        predicted = transfer.apply_intensive(old_value).values
        errors.append(float(np.sqrt(np.mean((predicted - exact) ** 2))))
    assert errors[1] < 0.6 * errors[0]
    assert errors[2] < 0.6 * errors[1]


def test_nonparallel_nonplanar_and_large_normal_motion_refuse():
    old = _grid_surface([0.0, 1.0], [0.0, 1.0])
    tilted_vertices = np.array(old.vertices, copy=True)
    tilted_vertices[-1, 2] = 0.1
    tilted = TriangleSurface3D(
        tilted_vertices, old.faces, old.face_material_id)
    with pytest.raises(ValueError, match="oriented parallel|noncoplanar"):
        _build(old, tilted, maximum_normal_distance=0.2)
    far = _grid_surface([0.0, 1.0], [0.0, 1.0], z=0.5)
    with pytest.raises(ValueError, match="large normal motion"):
        _build(old, far, maximum_normal_distance=1.0)


def test_planar_2048_face_runtime_is_bounded_against_indexed_knn_reference():
    coordinate = np.linspace(0.0, 1.0, 33)
    old = _grid_surface(
        coordinate, coordinate, periodic_lengths=(1.0, 1.0, None))
    new = _grid_surface(
        coordinate, coordinate, z=0.005, flip_diagonal=True,
        periodic_lengths=(1.0, 1.0, None))

    def knn_workload():
        return build_surface_transfer_3d(
            old, new, neighbor_count=4, maximum_distance=0.02)

    def overlap_workload():
        return _build(old, new)

    # Warm JIT and both candidate paths before interleaved timings.
    knn_workload()
    overlap_workload()
    knn_wall = []
    overlap_wall = []
    for _ in range(3):
        started = time.perf_counter()
        knn_workload()
        knn_wall.append(time.perf_counter() - started)
        started = time.perf_counter()
        overlap_workload()
        overlap_wall.append(time.perf_counter() - started)
    ratio = float(np.median(overlap_wall) / np.median(knn_wall))
    assert ratio <= 1.5, (
        f"planar overlap is {ratio:.3f}x indexed KNN; "
        f"knn={knn_wall}, overlap={overlap_wall}")
