import numpy as np
import pytest

from petch.block_levelset_3d import (
    ACTIVE_BLOCK,
    build_block_sparse_levelset_3d,
)
from petch.feature_step_3d import (
    FeatureGeometry3D,
    make_rectangular_trench_geometry_3d,
)


def _periodic_plane():
    dx = 0.125
    shape = (9, 5, 25)
    x, y, z = (np.arange(size) * dx for size in shape)
    _, _, Z = np.meshgrid(x, y, z, indexing="ij")
    phi = 1.5 - Z
    material = np.where(phi >= 0.0, 1, 0)
    return FeatureGeometry3D(
        phi, material, dx, 1e-6, material_levelsets={1: phi})


def test_sparse_blocks_cover_partial_boundaries_and_leave_unambiguous_far_labels():
    geometry = _periodic_plane()
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(3, 3, 5), band_width_cells=2,
        periodic_axes=(0, 1))

    assert sparse.block_grid_shape == (3, 2, 5)
    assert len(sparse.blocks) == 30
    assert sparse.blocks[-1].cell_shape == (2, 1, 4)
    assert sparse.active_blocks
    assert sparse.far_blocks
    assert all(block.label != ACTIVE_BLOCK for block in sparse.far_blocks)
    assert all(block.combined_phi is None for block in sparse.far_blocks)
    assert tuple(block.key for block in sparse.blocks) == tuple(sorted(
        block.key for block in sparse.blocks))


def test_exact_interface_band_and_truncated_far_view_are_consistent():
    geometry = _periodic_plane()
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(4, 2, 4), band_width_cells=2,
        periodic_axes=(0, 1))
    combined, owner, materials = sparse.to_truncated_dense()

    interface_cell = np.zeros(tuple(np.asarray(geometry.phi.shape) - 1), dtype=bool)
    for axis in range(3):
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        lower[axis] = slice(0, -1)
        upper[axis] = slice(1, None)
        changed = (
            (geometry.phi[tuple(lower)] >= 0.0)
            != (geometry.phi[tuple(upper)] >= 0.0))
        # Collapse node-oriented changes onto the common cell domain.
        trim = tuple(slice(0, size) for size in interface_cell.shape)
        interface_cell |= changed[trim]
    nodes = np.zeros(geometry.phi.shape, dtype=bool)
    for corner in np.ndindex((2, 2, 2)):
        selector = tuple(
            slice(corner[axis], corner[axis] + interface_cell.shape[axis])
            for axis in range(3))
        nodes[selector] |= interface_cell

    assert np.array_equal(combined[nodes], geometry.phi[nodes])
    assert np.array_equal(owner[nodes], geometry.material_id[nodes])
    assert np.array_equal(materials[1][nodes], geometry.phi[nodes])
    active_node = np.zeros(geometry.phi.shape, dtype=bool)
    for block in sparse.active_blocks:
        selector = tuple(slice(left, right + 1)
                         for left, right in zip(block.cell_start, block.cell_stop))
        active_node[selector] = True
    assert np.all(np.abs(combined[~active_node]) == sparse.truncation_distance)


def test_periodic_neighbor_and_halo_indices_use_one_canonical_endpoint():
    geometry = _periodic_plane()
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(4, 2, 4), band_width_cells=2,
        periodic_axes=(0, 1))
    packed = sparse.pack_active_halos(1)

    first_key = tuple(packed.keys[0])
    first_index = next(index for index, block in enumerate(sparse.blocks)
                       if block.key == first_key)
    x_minus = sparse.neighbor_index[first_index, 0, 0]
    assert x_minus >= 0
    assert sparse.blocks[x_minus].key[0] == sparse.block_grid_shape[0] - 1
    valid_index = packed.global_node_index[packed.valid_node]
    assert np.all(valid_index >= 0)
    assert np.all(valid_index[:, 0] < geometry.phi.shape[0] - 1)
    assert np.all(valid_index[:, 1] < geometry.phi.shape[1] - 1)
    assert np.all(valid_index[:, 2] < geometry.phi.shape[2])
    for location, accepted in np.ndenumerate(packed.valid_node):
        if not accepted:
            continue
        source = tuple(packed.global_node_index[location])
        assert packed.combined_phi[location] == sparse.to_truncated_dense()[0][source]


def test_sparse_fingerprint_and_packing_are_deterministic_and_read_only():
    geometry = _periodic_plane()
    first = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(4, 2, 4), band_width_cells=2,
        periodic_axes=(0, 1))
    second = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(4, 2, 4), band_width_cells=2,
        periodic_axes=(1, 0))

    assert first.fingerprint == second.fingerprint
    assert np.array_equal(first.neighbor_index, second.neighbor_index)
    packed = first.pack_active_halos(1)
    assert not packed.combined_phi.flags.writeable
    assert not packed.global_node_index.flags.writeable
    with pytest.raises(ValueError):
        packed.combined_phi.flat[0] = 1.0


def test_indexed_halo_pool_exactly_gathers_the_replicated_reference():
    geometry = _periodic_plane()
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(4, 2, 4), band_width_cells=2,
        periodic_axes=(0, 1))
    replicated = sparse.pack_active_halos(1)
    indexed = sparse.pack_active_indexed_halos(1)
    valid = indexed.valid_node
    gather = indexed.node_index[valid]

    assert np.array_equal(valid, replicated.valid_node)
    assert np.array_equal(
        indexed.global_node_index[gather],
        replicated.global_node_index[replicated.valid_node])
    assert np.array_equal(
        indexed.combined_phi[gather],
        replicated.combined_phi[replicated.valid_node])
    assert np.array_equal(
        indexed.material_owner[gather],
        replicated.material_owner[replicated.valid_node])
    assert all(np.array_equal(
        indexed.material_phi[material_id][gather],
        replicated.material_phi[material_id][replicated.valid_node])
        for material_id in sparse.material_ids)
    assert tuple(map(tuple, indexed.global_node_index)) == tuple(sorted(
        map(tuple, indexed.global_node_index)))


def test_multimaterial_krueger_band_preserves_junction_and_reports_storage():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.13, cell_length=0.02, domain_height=2.8,
        dx=0.01, opening_width=0.09, mask_thickness=0.85,
        substrate_top=1.8, etched_depth=0.0)
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(8, 8, 8), band_width_cells=8,
        periodic_axes=(0, 1))
    combined, owner, materials = sparse.to_truncated_dense()
    near = np.abs(geometry.phi) <= 2.0 * geometry.dx

    assert sparse.material_ids == (1, 2)
    assert np.array_equal(combined[near], geometry.phi[near])
    assert np.array_equal(owner[near], geometry.material_id[near])
    assert all(np.array_equal(
        materials[material_id][near], geometry.material_levelsets[material_id][near])
        for material_id in sparse.material_ids)
    receipt = sparse.storage_receipt(halo_cells=1)
    assert receipt["active_block_count"] < receipt["block_count"]
    assert receipt["core_memory_reduction"] > 1.0


def test_five_nm_krueger_roundoff_periodicity_is_canonicalized_and_priced():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.13, cell_length=0.02, domain_height=2.8,
        dx=0.005, opening_width=0.09, mask_thickness=0.85,
        substrate_top=1.8, etched_depth=0.0)
    sparse = build_block_sparse_levelset_3d(
        geometry, block_cell_shape=(8, 4, 16), band_width_cells=8,
        periodic_axes=(0, 1))
    combined, owner, materials = sparse.to_truncated_dense()

    assert 0.0 < sparse.periodic_endpoint_max_abs_difference
    assert (sparse.periodic_endpoint_max_abs_difference / geometry.dx) < 1.0e-9
    assert np.array_equal(combined[0], combined[-1])
    assert np.array_equal(owner[0], owner[-1])
    assert all(np.array_equal(field[0], field[-1])
               for field in materials.values())


def test_sparse_builder_refuses_false_periodicity_and_missing_material_authority():
    geometry = _periodic_plane()
    roundoff_phi = geometry.phi.copy()
    roundoff_phi[-1] += 0.5e-9 * geometry.dx
    roundoff = FeatureGeometry3D(
        roundoff_phi, np.where(roundoff_phi >= 0.0, 1, 0), geometry.dx,
        geometry.mesh_length_unit_m, material_levelsets={1: roundoff_phi})
    canonical = build_block_sparse_levelset_3d(roundoff, periodic_axes=(0,))
    assert canonical.periodic_endpoint_max_abs_difference == pytest.approx(
        0.5e-9 * geometry.dx)

    broken_phi = geometry.phi.copy()
    broken_phi[-1] += 0.1
    broken = FeatureGeometry3D(
        broken_phi, np.where(broken_phi >= 0.0, 1, 0), geometry.dx,
        geometry.mesh_length_unit_m, material_levelsets={1: broken_phi})
    with pytest.raises(ValueError, match="differs across its duplicate endpoint"):
        build_block_sparse_levelset_3d(broken, periodic_axes=(0,))

    multi = make_rectangular_trench_geometry_3d(
        cell_width=0.04, cell_length=0.02, domain_height=0.10,
        dx=0.01, opening_width=0.02, mask_thickness=0.02,
        substrate_top=0.05, etched_depth=0.0)
    unauthoritative = FeatureGeometry3D(
        multi.phi, multi.material_id, multi.dx, multi.mesh_length_unit_m,
        material_levelsets=None)
    with pytest.raises(ValueError, match="requires material level sets"):
        build_block_sparse_levelset_3d(unauthoritative)
