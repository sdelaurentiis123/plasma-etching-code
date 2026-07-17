import numpy as np
import pytest

from petch.surface_common_refinement_3d import (
    build_surface_common_refinement_transfer_3d,
)
from petch.surface_mesh_3d import TriangleSurface3D
from petch.surface_partitioned_overlap_3d import (
    build_partitioned_surface_overlap_transfer_3d,
)


def _square(*, z=0.0, x_shift=0.0, slope=0.0, reverse=False, periodic=None):
    vertex = np.asarray([
        [x_shift + 0.0, 0.0, z + slope * 0.0],
        [x_shift + 1.0, 0.0, z + slope * 1.0],
        [x_shift + 0.0, 1.0, z + slope * 0.0],
        [x_shift + 1.0, 1.0, z + slope * 1.0],
    ])
    face = np.asarray([[0, 1, 3], [0, 3, 2]], dtype=int)
    if reverse:
        face = face[:, [0, 2, 1]]
    return TriangleSurface3D(
        vertex, face, np.ones(2, dtype=int),
        periodic_lengths=(None, None, None) if periodic is None else periodic)


def _x_square(coordinate, *, periodic_length):
    vertex = np.asarray([
        [coordinate, 0.0, 0.0], [coordinate, 1.0, 0.0],
        [coordinate, 0.0, 1.0], [coordinate, 1.0, 1.0],
    ])
    face = np.asarray([[0, 1, 3], [0, 3, 2]], dtype=int)
    return TriangleSurface3D(
        vertex, face, np.ones(2, dtype=int),
        periodic_lengths=(periodic_length, None, None))


def test_identity_is_bitwise_and_inventory_exact():
    surface = _square()
    transfer = build_surface_common_refinement_transfer_3d(
        surface, surface, maximum_normal_distance=0.0)
    density = np.asarray([2.0, 7.0])
    result = transfer.apply_extensive(density, newly_exposed_density=3.0)

    np.testing.assert_array_equal(result.values, density)
    np.testing.assert_array_equal(transfer.old_covered_area, surface.face_area)
    np.testing.assert_array_equal(transfer.new_covered_area, surface.face_area)
    assert result.maximum_relative_balance_error == 0.0
    assert transfer.geometry_receipt["old_matched_area_fraction"] == 1.0
    assert transfer.geometry_receipt["new_matched_area_fraction"] == 1.0


def test_parallel_translation_matches_planar_exact_authority():
    old = _square()
    new = _square(z=0.01, x_shift=0.25)
    common = build_surface_common_refinement_transfer_3d(
        old, new, maximum_normal_distance=0.02)
    exact = build_partitioned_surface_overlap_transfer_3d(
        old, new, maximum_normal_distance=0.02)

    np.testing.assert_allclose(
        common.old_covered_area, exact.old_covered_area, rtol=0.0, atol=4e-15)
    np.testing.assert_allclose(
        common.new_covered_area, exact.new_covered_area, rtol=0.0, atol=4e-15)
    common_result = common.apply_extensive(
        np.asarray([2.0, 2.0]), newly_exposed_density=3.0)
    exact_result = exact.apply_extensive(
        np.asarray([2.0, 2.0]), newly_exposed_density=3.0)
    np.testing.assert_allclose(
        common_result.values, exact_result.values, rtol=0.0, atol=4e-14)


def test_tilted_new_surface_closes_retained_and_fresh_area_separately():
    old = _square()
    new = _square(z=0.01, slope=0.05)
    transfer = build_surface_common_refinement_transfer_3d(
        old, new, maximum_normal_distance=0.1)
    result = transfer.apply_extensive(
        np.full(2, 2.0), newly_exposed_density=np.full(2, 3.0))

    ledger = result.material_ledger[1]
    assert ledger["old_inventory"] == pytest.approx(2.0, abs=3e-14)
    assert ledger["removed_inventory"] == pytest.approx(0.0, abs=3e-14)
    assert ledger["newly_exposed_inventory"] > 0.0
    assert ledger["new_balance_residual"] == pytest.approx(0.0, abs=3e-14)
    assert result.maximum_relative_balance_error < 3e-14
    assert transfer.geometry_receipt["old_matched_area_fraction"] == pytest.approx(
        1.0, abs=3e-14)
    assert transfer.geometry_receipt["new_matched_area_fraction"] < 1.0


def test_opposite_new_sheet_is_created_not_a_teleported_predecessor():
    old = _square()
    new = _square(z=0.01, reverse=True)
    transfer = build_surface_common_refinement_transfer_3d(
        old, new, maximum_normal_distance=0.02)
    result = transfer.apply_extensive(
        np.full(2, 2.0), newly_exposed_density=np.full(2, 3.0))

    assert transfer.overlap_area.size == 0
    np.testing.assert_allclose(result.values, 3.0, rtol=0.0, atol=0.0)
    ledger = result.material_ledger[1]
    assert ledger["removed_inventory"] == pytest.approx(2.0)
    assert ledger["newly_exposed_inventory"] == pytest.approx(3.0)
    assert ledger["retained_inventory"] == 0.0


def test_periodic_normal_translation_uses_one_nearest_image():
    old = _x_square(0.98, periodic_length=1.0)
    new = _x_square(0.02, periodic_length=1.0)
    transfer = build_surface_common_refinement_transfer_3d(
        old, new, maximum_normal_distance=0.05)

    np.testing.assert_allclose(transfer.old_covered_area, old.face_area, atol=3e-15)
    np.testing.assert_allclose(transfer.new_covered_area, new.face_area, atol=3e-15)
    assert transfer.geometry_receipt["total_overlap_area"] == pytest.approx(1.0)
