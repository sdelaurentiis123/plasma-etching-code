from dataclasses import replace
import time

import numpy as np
import pytest

from petch.surface_mesh_3d import TriangleSurface3D
from petch.surface_transfer_3d import (
    _normalize_nonnegative_weights,
    build_surface_transfer_3d,
)


class _MixedRemapState:
    def __init__(self, inventory, fraction):
        self.inventory = np.asarray(inventory, dtype=float)
        self.fraction = np.asarray(fraction, dtype=float)

    def conservative_surface_fields(self):
        return {"inventory": self.inventory, "fraction": self.fraction}

    def conservative_surface_upper_bounds(self):
        return {"inventory": None, "fraction": 1.0}

    def surface_field_remap_modes(self):
        return {"inventory": "conservative", "fraction": "intensive"}

    def with_conservative_surface_fields(self, fields):
        return type(self)(fields["inventory"], fields["fraction"])


def _periodic_grid_pair(cell_count):
    coordinates = np.linspace(0.0, 1.0, int(cell_count) + 1)
    vertices = np.asarray([
        (x, y, 0.0) for x in coordinates for y in coordinates])
    faces = []
    for ix in range(int(cell_count)):
        for iy in range(int(cell_count)):
            lower = ix * (int(cell_count) + 1) + iy
            upper = (ix + 1) * (int(cell_count) + 1) + iy
            faces.extend(((lower, upper, upper + 1),
                          (lower, upper + 1, lower + 1)))
    faces = np.asarray(faces, dtype=int)
    material = np.where(
        np.mean(vertices[faces], axis=1)[:, 0] < 0.5, 1, 2)
    old = TriangleSurface3D(
        vertices, faces, material, periodic_lengths=(1.0, 1.0, None))
    moved = vertices.copy()
    moved[:, 2] = (
        0.01
        + 0.002 * np.sin(2.0 * np.pi * moved[:, 0])
        * np.sin(2.0 * np.pi * moved[:, 1]))
    new = TriangleSurface3D(
        moved, faces, material, periodic_lengths=(1.0, 1.0, None))
    return old, new


def _square_surface(z=0.0, *, material=(1, 1), periodic_lengths=(None, None, None)):
    return TriangleSurface3D(
        np.asarray([
            [0.0, 0.0, z], [1.0, 0.0, z],
            [1.0, 1.0, z], [0.0, 1.0, z],
        ]),
        np.asarray([[0, 1, 2], [0, 2, 3]]),
        np.asarray(material),
        periodic_lengths=periodic_lengths,
    )


def test_identity_has_exact_one_face_rows_and_bitwise_intensive_replay():
    surface = _square_surface()
    transfer = build_surface_transfer_3d(
        surface, surface, neighbor_count=4, maximum_distance=0.1)
    values = np.asarray([0.125, 17.0])
    applied = transfer.apply_intensive(values)
    np.testing.assert_array_equal(applied.values, values)
    np.testing.assert_array_equal(transfer.row_offsets, [0, 1, 2])
    np.testing.assert_array_equal(transfer.old_face_index, [0, 1])
    np.testing.assert_array_equal(transfer.weight, [1.0, 1.0])
    assert transfer.maximum_exact_surface_distance == 0.0
    assert transfer.maximum_nearest_centroid_distance == 0.0


def test_unrelated_plane_retriangulation_uses_exact_surface_not_centroid_distance():
    old = TriangleSurface3D(
        np.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]),
        np.asarray([[0, 1, 2]]), np.asarray([1]))
    new = TriangleSurface3D(
        np.asarray([
            [0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.4, 0.0],
            [1.6, 0.0, 0.0], [2.0, 0.0, 0.0], [1.6, 0.4, 0.0],
        ]),
        np.asarray([[0, 1, 2], [3, 4, 5]]), np.asarray([1, 1]))
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=4, maximum_distance=1e-12)
    assert transfer.maximum_exact_surface_distance <= 1e-15
    assert transfer.maximum_nearest_centroid_distance > 0.5
    np.testing.assert_array_equal(
        transfer.apply_intensive([3.25]).values, [3.25, 3.25])


def test_periodic_seam_uses_one_nearest_image_per_physical_face():
    old = TriangleSurface3D(
        np.asarray([[0.0, 0.98, 0.0], [0.1, 0.98, 0.0], [0.0, 1.0, 0.0]]),
        np.asarray([[0, 1, 2]]), np.asarray([1]),
        periodic_lengths=(None, 1.0, None))
    new = TriangleSurface3D(
        np.asarray([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.02, 0.0]]),
        np.asarray([[0, 1, 2]]), np.asarray([1]),
        periodic_lengths=(None, 1.0, None))
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=9, maximum_distance=0.014)
    np.testing.assert_array_equal(transfer.old_face_index, [0])
    np.testing.assert_allclose(transfer.source_periodic_shift, [[0.0, -1.0, 0.0]])
    assert transfer.maximum_exact_surface_distance < 0.014
    assert transfer.row_offsets[-1] == 1


def test_material_junction_refuses_instead_of_borrowing_nearby_other_material():
    old = TriangleSurface3D(
        np.asarray([
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [10.0, 0.0, 0.0], [10.0, 1.0, 0.0], [10.0, 0.0, 1.0],
        ]),
        np.asarray([[0, 1, 2], [3, 4, 5]]), np.asarray([1, 2]))
    new = TriangleSurface3D(
        np.asarray([
            [10.0, 0.0, 0.0], [10.0, 1.0, 0.0], [10.0, 0.0, 1.0],
            [10.0, 0.0, 0.1], [10.0, 1.0, 0.1], [10.0, 0.0, 1.1],
        ]),
        np.asarray([[0, 1, 2], [3, 4, 5]]), np.asarray([1, 2]))
    with pytest.raises(ValueError, match="borrowing across materials is forbidden"):
        build_surface_transfer_3d(
            old, new, neighbor_count=2, maximum_distance=0.2)


def test_rows_close_are_material_local_and_sparse_arrays_are_immutable():
    old = _square_surface(material=(1, 2))
    new = _square_surface(z=0.05, material=(1, 2))
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=2, maximum_distance=0.1)
    np.testing.assert_allclose(transfer.row_sum, 1.0, rtol=0.0, atol=8e-16)
    for row in range(len(new.faces)):
        start, stop = transfer.row_offsets[row:row + 2]
        source = transfer.old_face_index[start:stop]
        assert len(np.unique(source)) == len(source)
        assert np.all(old.face_material_id[source] == new.face_material_id[row])
    for array in (
            transfer.row_offsets, transfer.old_face_index, transfer.weight,
            transfer.source_periodic_shift, transfer.centroid_distance,
            transfer.row_sum):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.flat[0] = 9


def test_row_closure_preserves_a_tiny_positive_tail_without_negative_roundoff():
    # Exact binary64 witness captured from the failure class. The historical
    # final-entry closure produces -2.22e-16 even though every raw weight is
    # positive.
    raw = np.array([
        80513463366.02206,
        758073282.2955991,
        5.924681595697332e-16,
    ])
    historical = raw / float(np.sum(raw))
    assert 1.0 - float(np.sum(historical[:-1])) < 0.0

    closed = _normalize_nonnegative_weights(raw)

    assert np.all(closed >= 0.0)
    assert closed[-1] > 0.0
    assert float(np.sum(closed)) == 1.0
    np.testing.assert_allclose(
        closed, historical, rtol=3.0e-16, atol=0.0)


def test_row_closure_is_nonnegative_across_extreme_dynamic_range():
    generator = np.random.default_rng(20260821)
    for count in range(1, 9):
        raw = 10.0 ** generator.uniform(-250.0, 250.0, (500, count))
        for row in raw:
            closed = _normalize_nonnegative_weights(row)
            assert np.all(np.isfinite(closed))
            assert np.all(closed >= 0.0)
            assert float(np.sum(closed)) == pytest.approx(
                1.0, rel=0.0, abs=8e-16)


def test_sparse_contract_names_the_exact_rejected_invariant():
    old = _square_surface()
    new = _square_surface(z=0.05)
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=2, maximum_distance=0.1)
    bad_row_sum = transfer.row_sum.copy()
    bad_row_sum[1] = np.nextafter(1.0, 2.0, dtype=np.float64)
    bad_row_sum[1] = np.nextafter(bad_row_sum[1], 2.0, dtype=np.float64)
    bad_row_sum[1] = np.nextafter(bad_row_sum[1], 2.0, dtype=np.float64)
    bad_row_sum[1] = np.nextafter(bad_row_sum[1], 2.0, dtype=np.float64)
    with pytest.raises(
            ValueError,
            match=(
                r"invalid sparse surface-transfer weights: "
                r"row_sum_not_unit count=1, first_row=1, .*"
                r"absolute_tolerance=8e-16")):
        replace(transfer, row_sum=bad_row_sum)

    bad_weight = transfer.weight.copy()
    bad_weight[0] = np.nan
    with pytest.raises(
            ValueError,
            match=(
                r"invalid sparse surface-transfer weights: "
                r"nonfinite_weights count=1, first_entry=0")):
        replace(transfer, weight=bad_weight)


def test_deterministic_replay_and_order_sensitive_surface_contract():
    old = _square_surface()
    new = _square_surface(z=0.05)
    first = build_surface_transfer_3d(
        old, new, neighbor_count=2, maximum_distance=0.1)
    replay = build_surface_transfer_3d(
        TriangleSurface3D(
            np.asfortranarray(old.vertices),
            np.asfortranarray(old.faces.astype(np.int32)),
            old.face_material_id.astype(np.int32)),
        TriangleSurface3D(
            np.asfortranarray(new.vertices),
            np.asfortranarray(new.faces.astype(np.int32)),
            new.face_material_id.astype(np.int32)),
        neighbor_count=2, maximum_distance=0.1)
    assert first.fingerprint == replay.fingerprint
    np.testing.assert_array_equal(first.row_offsets, replay.row_offsets)
    np.testing.assert_array_equal(first.old_face_index, replay.old_face_index)
    np.testing.assert_array_equal(first.weight, replay.weight)


def test_exact_surface_distance_gate_refuses_remote_same_material():
    old = _square_surface()
    new = _square_surface(z=1.0)
    with pytest.raises(ValueError, match="surface transfer distance"):
        build_surface_transfer_3d(
            old, new, neighbor_count=2, maximum_distance=0.2)


def test_intensive_and_extensive_application_metadata_are_explicit():
    old = _square_surface()
    new = TriangleSurface3D(
        np.asarray([
            [0.0, 0.0, 0.05], [0.5, 0.0, 0.05], [0.5, 1.0, 0.05],
            [0.0, 1.0, 0.05], [1.0, 0.0, 0.05], [1.0, 1.0, 0.05],
        ]),
        np.asarray([[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]]),
        np.ones(4, dtype=int))
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=2, maximum_distance=0.1)
    intensive = transfer.apply_intensive([0.2, 0.8], lower_bound=0.0, upper_bound=1.0)
    extensive = transfer.apply_extensive([2.0, 4.0])
    assert intensive.semantics == "intensive"
    assert intensive.metadata["area_integral_preserved"] is False
    assert extensive.semantics == "extensive"
    assert extensive.metadata["area_integral_preserved"] is True
    assert extensive.maximum_relative_integral_error < 2e-15
    assert not intensive.values.flags.writeable
    assert not extensive.values.flags.writeable
    with pytest.raises(TypeError):
        extensive.metadata["area_integral_preserved"] = False
    with pytest.raises(TypeError):
        extensive.material_integrals[1]["new_area_integral"] = 0.0
    assert np.isclose(
        np.dot([2.0, 4.0], old.face_area),
        np.dot(extensive.values, new.face_area), rtol=2e-15)


def test_sparse_transfer_matches_legacy_manufactured_surface_state():
    from petch.feature_step_3d import conservative_remap_surface_state
    from petch.surface_kinetics import SiO2SurfaceState

    old = _square_surface(material=(1, 1))
    new = _square_surface(z=-0.05, material=(1, 1))
    state = SiO2SurfaceState(
        [0.2, 0.8], [1e18, 3e18], [2e17, 6e17],
        [0.1, 0.4], [0.25, 0.75])
    legacy, _ = conservative_remap_surface_state(
        state,
        old.face_centroid, old.face_area, old.face_material_id,
        new.face_centroid, new.face_area, new.face_material_id,
        dx=0.05, mesh_length_unit_m=1.0,
        neighbor_count=2, maximum_distance=0.1,
        old_triangles=old.triangles)
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=2, maximum_distance=0.1)
    for name, before in state.conservative_surface_fields().items():
        upper = state.conservative_surface_upper_bounds()[name]
        applied = transfer.apply_extensive(before, upper_bound=upper)
        np.testing.assert_allclose(
            applied.values,
            legacy.conservative_surface_fields()[name],
            rtol=2e-15, atol=0.0)


def test_indexed_shared_transfer_matches_brute_force_mixed_state_on_periodic_surface():
    old, new = _periodic_grid_pair(8)
    coordinate = np.linspace(0.0, 1.0, len(old.faces))
    state = _MixedRemapState(
        2.0e18 * (1.0 + coordinate),
        0.15 + 0.7 * coordinate)
    transfer = build_surface_transfer_3d(
        old, new, neighbor_count=4, maximum_distance=0.05)
    shared = {
        "inventory": transfer.apply_extensive(state.inventory).values,
        "fraction": transfer.apply_intensive(
            state.fraction, lower_bound=0.0, upper_bound=1.0).values,
    }
    brute_raw = {
        name: np.empty(len(new.faces), dtype=float)
        for name in state.conservative_surface_fields()
    }
    shifts = old._periodic_shifts
    for new_face, point in enumerate(new.face_centroid):
        material = int(new.face_material_id[new_face])
        old_index = np.flatnonzero(old.face_material_id == material)
        image = old.face_centroid[old_index, None, :] + shifts[None, :, :]
        image_distance = np.linalg.norm(image - point[None, None, :], axis=2)
        image_choice = np.argmin(image_distance, axis=1)
        distance = image_distance[np.arange(len(old_index)), image_choice]
        order = np.lexsort((old_index, distance))[:4]
        source = old_index[order]
        selected_distance = distance[order]
        selected_shift = shifts[image_choice[order]]
        floor = 64.0 * np.finfo(float).eps * max(
            0.05, float(np.max(np.abs(old.face_centroid[old_index]))),
            float(np.max(np.abs(point))), 1.0)
        if selected_distance[0] <= floor:
            source = source[:1]
            selected_distance = selected_distance[:1]
            selected_shift = selected_shift[:1]
            weight = np.ones(1)
        else:
            weight = old.face_area[source] / np.maximum(
                selected_distance ** 2, floor ** 2)
            weight /= np.sum(weight)
            weight[-1] = 1.0 - np.sum(weight[:-1])
        start, stop = transfer.row_offsets[new_face:new_face + 2]
        np.testing.assert_array_equal(
            transfer.old_face_index[start:stop], source)
        np.testing.assert_array_equal(
            transfer.source_periodic_shift[start:stop], selected_shift)
        np.testing.assert_allclose(
            transfer.centroid_distance[start:stop], selected_distance,
            rtol=0.0, atol=2e-15)
        np.testing.assert_allclose(
            transfer.weight[start:stop], weight, rtol=0.0, atol=2e-15)
        for name, value in state.conservative_surface_fields().items():
            brute_raw[name][new_face] = float(np.dot(weight, value[source]))

    expected_inventory = brute_raw["inventory"].copy()
    for material in (1, 2):
        old_selected = old.face_material_id == material
        new_selected = new.face_material_id == material
        target = float(np.dot(
            state.inventory[old_selected], old.face_area[old_selected]))
        current = float(np.dot(
            expected_inventory[new_selected], new.face_area[new_selected]))
        expected_inventory[new_selected] *= target / current
    np.testing.assert_allclose(
        shared["inventory"], expected_inventory, rtol=3e-14, atol=0.0)
    np.testing.assert_allclose(
        shared["fraction"], brute_raw["fraction"], rtol=3e-14, atol=0.0)


def test_indexed_shared_transfer_meets_2048_face_legacy_runtime_gate():
    from petch.feature_step_3d import conservative_remap_surface_state

    old, new = _periodic_grid_pair(32)
    coordinate = np.linspace(0.0, 1.0, len(old.faces))
    state = _MixedRemapState(
        2.0e18 * (1.0 + coordinate),
        0.15 + 0.7 * coordinate)
    legacy_arguments = dict(
        dx=0.01, mesh_length_unit_m=1.0,
        neighbor_count=4, maximum_distance=0.05,
        old_triangles=old.triangles,
        periodic_lengths=old.periodic_lengths)

    def legacy_workload():
        return conservative_remap_surface_state(
            state,
            old.face_centroid, old.face_area, old.face_material_id,
            new.face_centroid, new.face_area, new.face_material_id,
            **legacy_arguments)

    def shared_workload():
        transfer = build_surface_transfer_3d(
            old, new, neighbor_count=4, maximum_distance=0.05)
        transfer.apply_extensive(state.inventory)
        transfer.apply_intensive(
            state.fraction, lower_bound=0.0, upper_bound=1.0)
        return transfer

    # Warm both code paths, then interleave three measurements so host load
    # does not systematically favor either implementation.
    legacy_workload()
    shared_workload()
    legacy_wall = []
    shared_wall = []
    for _ in range(3):
        started = time.perf_counter()
        legacy_workload()
        legacy_wall.append(time.perf_counter() - started)
        started = time.perf_counter()
        shared_workload()
        shared_wall.append(time.perf_counter() - started)
    ratio = float(np.median(shared_wall) / np.median(legacy_wall))
    assert ratio <= 1.5, (
        f"indexed shared transfer is {ratio:.3f}x legacy on 2048 faces; "
        f"legacy={legacy_wall}, shared={shared_wall}")
