"""Deterministic fixed-spacing block-sparse narrow-band geometry.

This module is deliberately standalone.  It classifies and packs an existing
uniform :class:`FeatureGeometry3D` without changing the live profile operator.
Every interface and material junction plus a declared buffer remains exact;
far blocks carry only an unambiguous gas/material label and reconstruct a
truncated signed distance.  The representation is the correctness precursor
to sparse evolution and one-level 2:1 AMR.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import product
from types import MappingProxyType

import numpy as np

from .feature_geometry_state_3d import FeatureGeometry3D


ACTIVE_BLOCK = -1
FAR_GAS = 0
_SCHEMA = b"petch-block-sparse-levelset-3d-v1"
PERIODIC_ENDPOINT_TOLERANCE_CELLS = 1.0e-9


def _readonly(value, dtype=None):
    array = np.asarray(value, dtype=dtype).copy()
    array.setflags(write=False)
    return array


def _as_three_positive_ints(value, name):
    raw = tuple(value)
    if (len(raw) != 3
            or any(isinstance(item, (bool, np.bool_)) for item in raw)
            or any(int(item) != item or int(item) <= 0 for item in raw)):
        raise ValueError(f"{name} must contain three positive integers")
    return tuple(int(item) for item in raw)


def _periodic_axes(value):
    raw = tuple(value)
    if any(isinstance(axis, (bool, np.bool_)) or int(axis) != axis for axis in raw):
        raise ValueError("periodic axes must be unique integer axes")
    axes = tuple(sorted(int(axis) for axis in raw))
    if len(set(axes)) != len(axes) or any(axis < 0 or axis >= 3 for axis in axes):
        raise ValueError("periodic axes must be unique integer axes")
    return axes


def _shift_mask(mask, offset, periodic_axes):
    """Shift a canonical-node mask, wrapping only declared physical axes."""
    output = np.asarray(mask, dtype=bool)
    for axis, delta in enumerate(offset):
        if delta == 0:
            continue
        if axis in periodic_axes:
            output = np.roll(output, int(delta), axis=axis)
            continue
        shifted = np.zeros_like(output)
        source = [slice(None)] * 3
        target = [slice(None)] * 3
        if delta > 0:
            source[axis] = slice(0, -delta)
            target[axis] = slice(delta, None)
        else:
            source[axis] = slice(-delta, None)
            target[axis] = slice(0, delta)
        shifted[tuple(target)] = output[tuple(source)]
        output = shifted
    return output


def _dilate(mask, iterations, periodic_axes):
    current = np.asarray(mask, dtype=bool).copy()
    offsets = tuple(product((-1, 0, 1), repeat=3))
    for _ in range(int(iterations)):
        expanded = np.zeros_like(current)
        for offset in offsets:
            expanded |= _shift_mask(current, offset, periodic_axes)
        current = expanded
    return current


def _canonical_view(array, periodic_axes):
    selector = tuple(
        slice(0, size - 1) if axis in periodic_axes else slice(None)
        for axis, size in enumerate(array.shape))
    return np.asarray(array)[selector]


def _expand_periodic_endpoints(canonical, full_shape, periodic_axes):
    index = []
    for axis, size in enumerate(full_shape):
        values = np.arange(size, dtype=int)
        if axis in periodic_axes:
            values[-1] = 0
        index.append(values)
    return np.asarray(canonical)[np.ix_(*index)]


def _interface_seed(labels, phi, periodic_axes):
    labels = np.asarray(labels, dtype=int)
    seed = np.abs(np.asarray(phi, dtype=float)) <= (
        64.0 * np.finfo(float).eps
        * max(float(np.max(np.abs(phi), initial=0.0)), 1.0))
    for axis in range(3):
        if labels.shape[axis] <= 1:
            continue
        lower = [slice(None)] * 3
        upper = [slice(None)] * 3
        lower[axis] = slice(0, -1)
        upper[axis] = slice(1, None)
        changed = labels[tuple(lower)] != labels[tuple(upper)]
        seed[tuple(lower)] |= changed
        seed[tuple(upper)] |= changed
        if axis in periodic_axes:
            first = [slice(None)] * 3
            last = [slice(None)] * 3
            first[axis] = 0
            last[axis] = -1
            wrapped = labels[tuple(first)] != labels[tuple(last)]
            seed[tuple(first)] |= wrapped
            seed[tuple(last)] |= wrapped
    return seed


def _digest_array(digest, name, value, dtype):
    array = np.ascontiguousarray(value, dtype=dtype)
    encoded = str(name).encode("utf-8")
    digest.update(np.asarray([len(encoded)], dtype="<u8").tobytes())
    digest.update(encoded)
    digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
    digest.update(array.tobytes())


@dataclass(frozen=True)
class SparseLevelSetBlock3D:
    """One lexicographically indexed block of physical cells and their nodes."""

    key: tuple[int, int, int]
    cell_start: tuple[int, int, int]
    cell_stop: tuple[int, int, int]
    label: int
    combined_phi: np.ndarray | None = None
    material_phi: object | None = None
    material_owner: np.ndarray | None = None

    def __post_init__(self):
        key = tuple(int(value) for value in self.key)
        start = tuple(int(value) for value in self.cell_start)
        stop = tuple(int(value) for value in self.cell_stop)
        label = int(self.label)
        if (len(key) != 3 or len(start) != 3 or len(stop) != 3
                or any(value < 0 for value in key + start)
                or any(right <= left for left, right in zip(start, stop))
                or label < ACTIVE_BLOCK):
            raise ValueError("invalid sparse level-set block metadata")
        expected = tuple(right - left + 1 for left, right in zip(start, stop))
        if label == ACTIVE_BLOCK:
            combined = _readonly(self.combined_phi, float)
            owner = _readonly(self.material_owner, int)
            fields = {
                int(material_id): _readonly(field, float)
                for material_id, field in dict(self.material_phi or {}).items()
            }
            if (combined.shape != expected or owner.shape != expected or not fields
                    or any(material_id <= 0 or field.shape != expected
                           for material_id, field in fields.items())
                    or np.any(~np.isfinite(combined))
                    or any(np.any(~np.isfinite(field)) for field in fields.values())
                    or np.any(owner < 0)):
                raise ValueError("invalid active sparse block payload")
            object.__setattr__(self, "combined_phi", combined)
            object.__setattr__(self, "material_phi", MappingProxyType(fields))
            object.__setattr__(self, "material_owner", owner)
        elif any(value is not None for value in (
                self.combined_phi, self.material_phi, self.material_owner)):
            raise ValueError("far sparse blocks may carry only one uniform label")
        object.__setattr__(self, "key", key)
        object.__setattr__(self, "cell_start", start)
        object.__setattr__(self, "cell_stop", stop)
        object.__setattr__(self, "label", label)

    @property
    def active(self):
        return self.label == ACTIVE_BLOCK

    @property
    def cell_shape(self):
        return tuple(
            right - left for left, right in zip(self.cell_start, self.cell_stop))

    @property
    def node_shape(self):
        return tuple(value + 1 for value in self.cell_shape)


@dataclass(frozen=True)
class PackedSparseLevelSetBlocks3D:
    """Fixed-shape active-block payload with deterministic stencil halos."""

    keys: np.ndarray
    cell_start: np.ndarray
    cell_stop: np.ndarray
    global_node_index: np.ndarray
    valid_node: np.ndarray
    combined_phi: np.ndarray
    material_phi: object
    material_owner: np.ndarray
    halo_cells: int

    def __post_init__(self):
        keys = _readonly(self.keys, np.int32)
        start = _readonly(self.cell_start, np.int32)
        stop = _readonly(self.cell_stop, np.int32)
        index = _readonly(self.global_node_index, np.int32)
        valid = _readonly(self.valid_node, bool)
        combined = _readonly(self.combined_phi, float)
        owner = _readonly(self.material_owner, int)
        fields = {
            int(material_id): _readonly(field, float)
            for material_id, field in dict(self.material_phi).items()
        }
        count = len(keys)
        local_shape = valid.shape[1:]
        if (keys.shape != (count, 3) or start.shape != (count, 3)
                or stop.shape != (count, 3)
                or index.shape != (count, *local_shape, 3)
                or combined.shape != valid.shape or owner.shape != valid.shape
                or any(field.shape != valid.shape for field in fields.values())
                or not fields or int(self.halo_cells) < 0
                or np.any(valid & ~np.isfinite(combined))
                or any(np.any(valid & ~np.isfinite(field)) for field in fields.values())):
            raise ValueError("invalid packed sparse-block halo payload")
        for array in (keys, start, stop, index, valid, combined, owner):
            array.setflags(write=False)
        object.__setattr__(self, "keys", keys)
        object.__setattr__(self, "cell_start", start)
        object.__setattr__(self, "cell_stop", stop)
        object.__setattr__(self, "global_node_index", index)
        object.__setattr__(self, "valid_node", valid)
        object.__setattr__(self, "combined_phi", combined)
        object.__setattr__(self, "material_phi", MappingProxyType(fields))
        object.__setattr__(self, "material_owner", owner)
        object.__setattr__(self, "halo_cells", int(self.halo_cells))


@dataclass(frozen=True)
class IndexedSparseLevelSetBlocks3D:
    """Unique canonical node pool plus fixed block-to-node stencil indices."""

    keys: np.ndarray
    cell_start: np.ndarray
    cell_stop: np.ndarray
    node_index: np.ndarray
    global_node_index: np.ndarray
    combined_phi: np.ndarray
    material_phi: object
    material_owner: np.ndarray
    halo_cells: int

    def __post_init__(self):
        keys = _readonly(self.keys, np.int32)
        start = _readonly(self.cell_start, np.int32)
        stop = _readonly(self.cell_stop, np.int32)
        node_index = _readonly(self.node_index, np.int32)
        global_index = _readonly(self.global_node_index, np.int32)
        combined = _readonly(self.combined_phi, float)
        owner = _readonly(self.material_owner, int)
        fields = {
            int(material_id): _readonly(field, float)
            for material_id, field in dict(self.material_phi).items()}
        block_count = len(keys)
        node_count = len(global_index)
        if (keys.shape != (block_count, 3) or start.shape != (block_count, 3)
                or stop.shape != (block_count, 3)
                or node_index.shape[0] != block_count
                or global_index.shape != (node_count, 3)
                or combined.shape != (node_count,) or owner.shape != (node_count,)
                or any(field.shape != (node_count,) for field in fields.values())
                or not fields or int(self.halo_cells) < 0
                or np.any(node_index < -1) or np.any(node_index >= node_count)
                or np.any(~np.isfinite(combined))
                or any(np.any(~np.isfinite(field)) for field in fields.values())):
            raise ValueError("invalid indexed sparse-block payload")
        object.__setattr__(self, "keys", keys)
        object.__setattr__(self, "cell_start", start)
        object.__setattr__(self, "cell_stop", stop)
        object.__setattr__(self, "node_index", node_index)
        object.__setattr__(self, "global_node_index", global_index)
        object.__setattr__(self, "combined_phi", combined)
        object.__setattr__(self, "material_phi", MappingProxyType(fields))
        object.__setattr__(self, "material_owner", owner)
        object.__setattr__(self, "halo_cells", int(self.halo_cells))

    @property
    def valid_node(self):
        return self.node_index >= 0


@dataclass(frozen=True)
class BlockSparseLevelSet3D:
    """Immutable active/far block decomposition at one global spacing."""

    shape: tuple[int, int, int]
    dx: float
    mesh_length_unit_m: float
    mesh_origin_m: tuple[float, float, float]
    material_ids: tuple[int, ...]
    block_cell_shape: tuple[int, int, int]
    band_width_cells: int
    periodic_axes: tuple[int, ...]
    periodic_endpoint_max_abs_difference: float
    blocks: tuple[SparseLevelSetBlock3D, ...]
    neighbor_index: np.ndarray

    def __post_init__(self):
        shape = tuple(int(value) for value in self.shape)
        block_shape = _as_three_positive_ints(
            self.block_cell_shape, "block_cell_shape")
        axes = _periodic_axes(self.periodic_axes)
        materials = tuple(sorted(int(value) for value in self.material_ids))
        blocks = tuple(self.blocks)
        neighbor = _readonly(self.neighbor_index, np.int32)
        if (len(shape) != 3 or any(value < 2 for value in shape)
                or not np.isfinite(self.dx) or self.dx <= 0.0
                or not np.isfinite(self.mesh_length_unit_m)
                or self.mesh_length_unit_m <= 0.0
                or len(self.mesh_origin_m) != 3
                or np.any(~np.isfinite(self.mesh_origin_m))
                or not materials or any(value <= 0 for value in materials)
                or int(self.band_width_cells) != self.band_width_cells
                or int(self.band_width_cells) < 1
                or not np.isfinite(self.periodic_endpoint_max_abs_difference)
                or self.periodic_endpoint_max_abs_difference < 0.0
                or not blocks or neighbor.shape != (len(blocks), 3, 2)):
            raise ValueError("invalid block-sparse level-set metadata")
        keys = tuple(block.key for block in blocks)
        if keys != tuple(sorted(keys)) or len(set(keys)) != len(keys):
            raise ValueError("sparse blocks must have unique lexicographic keys")
        expected_count = tuple(
            int(np.ceil((size - 1) / width))
            for size, width in zip(shape, block_shape))
        if len(blocks) != int(np.prod(expected_count)):
            raise ValueError("sparse block table does not cover the physical domain")
        if any(
                block.cell_start != tuple(key[axis] * block_shape[axis]
                                          for axis in range(3))
                or block.cell_stop != tuple(min(
                    (key[axis] + 1) * block_shape[axis], shape[axis] - 1)
                    for axis in range(3))
                for key, block in zip(keys, blocks)):
            raise ValueError("sparse block bounds do not match their keys")
        neighbor.setflags(write=False)
        object.__setattr__(self, "shape", shape)
        object.__setattr__(self, "block_cell_shape", block_shape)
        object.__setattr__(self, "periodic_axes", axes)
        object.__setattr__(self, "material_ids", materials)
        object.__setattr__(self, "band_width_cells", int(self.band_width_cells))
        object.__setattr__(
            self, "periodic_endpoint_max_abs_difference",
            float(self.periodic_endpoint_max_abs_difference))
        object.__setattr__(self, "blocks", blocks)
        object.__setattr__(self, "neighbor_index", neighbor)
        object.__setattr__(
            self, "mesh_origin_m", tuple(float(value) for value in self.mesh_origin_m))

    @property
    def block_grid_shape(self):
        return tuple(
            int(np.ceil((size - 1) / width))
            for size, width in zip(self.shape, self.block_cell_shape))

    @property
    def active_blocks(self):
        return tuple(block for block in self.blocks if block.active)

    @property
    def far_blocks(self):
        return tuple(block for block in self.blocks if not block.active)

    @property
    def truncation_distance(self):
        return float(self.band_width_cells * self.dx)

    @property
    def fingerprint(self):
        digest = sha256()
        digest.update(_SCHEMA)
        _digest_array(digest, "shape", self.shape, "<i8")
        _digest_array(digest, "metadata", (
            self.dx, self.mesh_length_unit_m, *self.mesh_origin_m,
            self.band_width_cells, self.periodic_endpoint_max_abs_difference), "<f8")
        _digest_array(digest, "block_shape", self.block_cell_shape, "<i8")
        _digest_array(digest, "periodic", self.periodic_axes, "<i8")
        _digest_array(digest, "materials", self.material_ids, "<i8")
        _digest_array(digest, "neighbors", self.neighbor_index, "<i8")
        for block in self.blocks:
            _digest_array(digest, "block", (
                *block.key, *block.cell_start, *block.cell_stop, block.label), "<i8")
            if block.active:
                _digest_array(digest, "combined", block.combined_phi, "<f8")
                _digest_array(digest, "owner", block.material_owner, "<i8")
                for material_id in self.material_ids:
                    _digest_array(
                        digest, f"material-{material_id}",
                        block.material_phi[material_id], "<f8")
        return digest.hexdigest()

    def to_truncated_dense(self):
        """Return a development-only dense view of exact band and labeled far blocks."""
        threshold = self.truncation_distance
        combined = np.full(self.shape, np.nan, dtype=float)
        owner = np.full(self.shape, -1, dtype=int)
        materials = {
            material_id: np.full(self.shape, np.nan, dtype=float)
            for material_id in self.material_ids}
        for block in self.far_blocks:
            selector = tuple(
                slice(left, right + 1)
                for left, right in zip(block.cell_start, block.cell_stop))
            combined[selector] = threshold if block.label > 0 else -threshold
            owner[selector] = block.label
            for material_id in self.material_ids:
                materials[material_id][selector] = (
                    threshold if block.label == material_id else -threshold)
        # Exact active payload has authority on shared block-boundary nodes.
        for block in self.active_blocks:
            selector = tuple(
                slice(left, right + 1)
                for left, right in zip(block.cell_start, block.cell_stop))
            combined[selector] = block.combined_phi
            owner[selector] = block.material_owner
            for material_id in self.material_ids:
                materials[material_id][selector] = block.material_phi[material_id]
        if (np.any(~np.isfinite(combined)) or np.any(owner < 0)
                or any(np.any(~np.isfinite(field)) for field in materials.values())):
            raise RuntimeError("sparse block table left an uncovered diagnostic node")
        return combined, owner, MappingProxyType(materials)

    def pack_active_halos(self, halo_cells=1):
        """Pack fixed-shape exact/truncated halos around every active block."""
        if (isinstance(halo_cells, (bool, np.bool_))
                or int(halo_cells) != halo_cells or int(halo_cells) < 0
                or int(halo_cells) > self.band_width_cells):
            raise ValueError("halo_cells must be an integer inside the stored band")
        halo = int(halo_cells)
        active = self.active_blocks
        local_shape = tuple(value + 1 + 2 * halo for value in self.block_cell_shape)
        payload_shape = (len(active), *local_shape)
        keys = np.asarray([block.key for block in active], dtype=np.int32)
        start = np.asarray([block.cell_start for block in active], dtype=np.int32)
        stop = np.asarray([block.cell_stop for block in active], dtype=np.int32)
        index = np.full((*payload_shape, 3), -1, dtype=np.int32)
        valid = np.zeros(payload_shape, dtype=bool)
        combined = np.full(payload_shape, np.nan, dtype=float)
        owner = np.full(payload_shape, -1, dtype=int)
        fields = {
            material_id: np.full(payload_shape, np.nan, dtype=float)
            for material_id in self.material_ids}
        interval = np.asarray(self.shape) - 1
        exact_node = {}
        for block in active:
            for local in np.ndindex(block.node_shape):
                node = np.asarray(block.cell_start) + np.asarray(local)
                for axis in self.periodic_axes:
                    node[axis] %= interval[axis]
                key = tuple(int(value) for value in node)
                value = (
                    float(block.combined_phi[local]),
                    int(block.material_owner[local]),
                    tuple(float(block.material_phi[material_id][local])
                          for material_id in self.material_ids),
                )
                previous = exact_node.setdefault(key, value)
                if previous != value:
                    raise RuntimeError(
                        "shared active-block node has inconsistent exact payload")
        block_by_key = {block.key: block for block in self.blocks}
        for block_index, block in enumerate(active):
            cell_shape = np.asarray(block.cell_shape)
            for local in np.ndindex(local_shape):
                offset = np.asarray(local) - halo
                if np.any(offset > cell_shape + halo):
                    continue
                global_node = np.asarray(block.cell_start) + offset
                accepted = True
                for axis in range(3):
                    if axis in self.periodic_axes:
                        global_node[axis] %= interval[axis]
                    elif global_node[axis] < 0 or global_node[axis] >= self.shape[axis]:
                        accepted = False
                        break
                if not accepted:
                    continue
                destination = (block_index, *local)
                source = tuple(int(value) for value in global_node)
                valid[destination] = True
                index[destination] = global_node
                if source in exact_node:
                    exact_combined, exact_owner, exact_material = exact_node[source]
                    combined[destination] = exact_combined
                    owner[destination] = exact_owner
                    for material_index, material_id in enumerate(self.material_ids):
                        fields[material_id][destination] = exact_material[material_index]
                    continue
                cell = np.minimum(global_node, interval - 1)
                key = tuple(int(cell[axis] // self.block_cell_shape[axis])
                            for axis in range(3))
                far = block_by_key[key]
                if far.active:
                    raise RuntimeError("active halo node is missing its exact sparse payload")
                combined[destination] = (
                    self.truncation_distance if far.label > 0
                    else -self.truncation_distance)
                owner[destination] = far.label
                for material_id in self.material_ids:
                    fields[material_id][destination] = (
                        self.truncation_distance
                        if far.label == material_id else -self.truncation_distance)
        return PackedSparseLevelSetBlocks3D(
            keys, start, stop, index, valid, combined, fields, owner, halo)

    def pack_active_indexed_halos(self, halo_cells=1):
        """Pack halos through one sorted unique-node pool instead of replicated values."""
        if (isinstance(halo_cells, (bool, np.bool_))
                or int(halo_cells) != halo_cells or int(halo_cells) < 0
                or int(halo_cells) > self.band_width_cells):
            raise ValueError("halo_cells must be an integer inside the stored band")
        halo = int(halo_cells)
        active = self.active_blocks
        local_shape = tuple(value + 1 + 2 * halo for value in self.block_cell_shape)
        keys = np.asarray([block.key for block in active], dtype=np.int32)
        start = np.asarray([block.cell_start for block in active], dtype=np.int32)
        stop = np.asarray([block.cell_stop for block in active], dtype=np.int32)
        interval = np.asarray(self.shape) - 1

        exact_node = {}
        for block in active:
            for local in np.ndindex(block.node_shape):
                node = np.asarray(block.cell_start) + np.asarray(local)
                for axis in self.periodic_axes:
                    node[axis] %= interval[axis]
                key = tuple(int(value) for value in node)
                value = (
                    float(block.combined_phi[local]),
                    int(block.material_owner[local]),
                    tuple(float(block.material_phi[material_id][local])
                          for material_id in self.material_ids),
                )
                previous = exact_node.setdefault(key, value)
                if previous != value:
                    raise RuntimeError(
                        "shared active-block node has inconsistent exact payload")
        block_by_key = {block.key: block for block in self.blocks}

        def canonical_source(block, local):
            offset = np.asarray(local) - halo
            if np.any(offset > np.asarray(block.cell_shape) + halo):
                return None
            node = np.asarray(block.cell_start) + offset
            for axis in range(3):
                if axis in self.periodic_axes:
                    node[axis] %= interval[axis]
                elif node[axis] < 0 or node[axis] >= self.shape[axis]:
                    return None
            return tuple(int(value) for value in node)

        sources = set()
        for block in active:
            for local in np.ndindex(local_shape):
                source = canonical_source(block, local)
                if source is not None:
                    sources.add(source)
        ordered_source = tuple(sorted(sources))
        source_index = {source: index for index, source in enumerate(ordered_source)}
        node_index = np.full((len(active), *local_shape), -1, dtype=np.int32)
        for block_index, block in enumerate(active):
            for local in np.ndindex(local_shape):
                source = canonical_source(block, local)
                if source is not None:
                    node_index[(block_index, *local)] = source_index[source]

        combined = np.empty(len(ordered_source), dtype=float)
        owner = np.empty(len(ordered_source), dtype=int)
        fields = {
            material_id: np.empty(len(ordered_source), dtype=float)
            for material_id in self.material_ids}
        for pool_index, source in enumerate(ordered_source):
            if source in exact_node:
                value_combined, value_owner, value_material = exact_node[source]
            else:
                cell = np.minimum(np.asarray(source), interval - 1)
                block_key = tuple(
                    int(cell[axis] // self.block_cell_shape[axis])
                    for axis in range(3))
                far = block_by_key[block_key]
                if far.active:
                    raise RuntimeError("active indexed node lacks exact sparse payload")
                value_combined = (
                    self.truncation_distance if far.label > 0
                    else -self.truncation_distance)
                value_owner = far.label
                value_material = tuple(
                    self.truncation_distance
                    if far.label == material_id else -self.truncation_distance
                    for material_id in self.material_ids)
            combined[pool_index] = value_combined
            owner[pool_index] = value_owner
            for material_index, material_id in enumerate(self.material_ids):
                fields[material_id][pool_index] = value_material[material_index]
        return IndexedSparseLevelSetBlocks3D(
            keys, start, stop, node_index,
            np.asarray(ordered_source, dtype=np.int32),
            combined, fields, owner, halo)

    def storage_receipt(self, *, halo_cells=None):
        """Report stored bytes and node fractions without claiming kernel speed."""
        dense_field_count = 2 + len(self.material_ids)  # owner + combined + each material
        dense_bytes = int(np.prod(self.shape)) * (
            8 * (1 + len(self.material_ids)) + 8)
        active_core_nodes = sum(int(np.prod(block.node_shape))
                                for block in self.active_blocks)
        sparse_bytes = active_core_nodes * (
            8 * (1 + len(self.material_ids)) + 8)
        sparse_bytes += len(self.blocks) * (10 * 8 + 3 * 2 * 4)
        payload = {
            "dense_node_count": int(np.prod(self.shape)),
            "dense_field_count": dense_field_count,
            "dense_bytes": dense_bytes,
            "block_count": len(self.blocks),
            "active_block_count": len(self.active_blocks),
            "far_block_count": len(self.far_blocks),
            "active_core_node_count_with_shared_boundaries": active_core_nodes,
            "sparse_core_bytes": int(sparse_bytes),
            "core_memory_reduction": float(dense_bytes / sparse_bytes),
            "active_block_fraction": float(len(self.active_blocks) / len(self.blocks)),
            "periodic_endpoint_max_abs_difference": (
                self.periodic_endpoint_max_abs_difference),
            "periodic_endpoint_max_difference_cells": float(
                self.periodic_endpoint_max_abs_difference / self.dx),
        }
        if halo_cells is not None:
            packed = self.pack_active_halos(halo_cells)
            halo_bytes = (
                packed.combined_phi.nbytes + packed.material_owner.nbytes
                + packed.valid_node.nbytes + packed.global_node_index.nbytes
                + sum(field.nbytes for field in packed.material_phi.values()))
            payload["packed_halo_cells"] = int(halo_cells)
            payload["packed_halo_bytes"] = int(halo_bytes)
            payload["memory_reduction_with_packed_halo"] = float(
                dense_bytes / (sparse_bytes + halo_bytes))
            indexed = self.pack_active_indexed_halos(halo_cells)
            indexed_bytes = (
                indexed.node_index.nbytes + indexed.global_node_index.nbytes
                + indexed.combined_phi.nbytes + indexed.material_owner.nbytes
                + sum(field.nbytes for field in indexed.material_phi.values()))
            payload["indexed_unique_node_count"] = len(indexed.global_node_index)
            payload["indexed_halo_bytes"] = int(indexed_bytes)
            payload["indexed_halo_only_memory_reduction"] = float(
                dense_bytes / indexed_bytes)
            payload["memory_reduction_with_indexed_halo"] = float(
                dense_bytes / (sparse_bytes + indexed_bytes))
        return payload


def _canonicalize_periodic_geometry(geometry, periodic_axes, material_fields):
    combined = np.asarray(geometry.phi, dtype=float).copy()
    owner = np.asarray(geometry.material_id, dtype=int).copy()
    fields = {
        int(material_id): np.asarray(field, dtype=float).copy()
        for material_id, field in material_fields.items()}
    maximum_difference = 0.0
    for axis in periodic_axes:
        first = [slice(None)] * 3
        last = [slice(None)] * 3
        first[axis] = 0
        last[axis] = -1
        if not np.array_equal(owner[tuple(first)], owner[tuple(last)]):
            raise ValueError(
                f"periodic axis {axis} changes material ownership across its duplicate endpoint")
        for name, array in (("combined", combined), *(
                (f"material-{material_id}", field)
                for material_id, field in sorted(fields.items()))):
            first_value = array[tuple(first)]
            last_value = array[tuple(last)]
            difference = float(np.max(np.abs(first_value - last_value), initial=0.0))
            scale = max(
                float(np.max(np.abs(first_value), initial=0.0)),
                float(np.max(np.abs(last_value), initial=0.0)),
                float(geometry.dx), 1.0)
            tolerance = max(
                512.0 * np.finfo(float).eps * scale,
                PERIODIC_ENDPOINT_TOLERANCE_CELLS * float(geometry.dx))
            if difference > tolerance:
                raise ValueError(
                    f"periodic axis {axis} field {name} differs across its duplicate "
                    f"endpoint by {difference:g} > {tolerance:g}")
            maximum_difference = max(maximum_difference, difference)
            array[tuple(last)] = first_value
        owner[tuple(last)] = owner[tuple(first)]
    return combined, owner, fields, maximum_difference


def _neighbor_table(block_grid_shape, periodic_axes):
    keys = tuple(np.ndindex(block_grid_shape))
    lookup = {key: index for index, key in enumerate(keys)}
    neighbor = np.full((len(keys), 3, 2), -1, dtype=np.int32)
    for index, key in enumerate(keys):
        for axis in range(3):
            for side_index, delta in enumerate((-1, 1)):
                adjacent = list(key)
                adjacent[axis] += delta
                if axis in periodic_axes:
                    adjacent[axis] %= block_grid_shape[axis]
                elif not 0 <= adjacent[axis] < block_grid_shape[axis]:
                    continue
                neighbor[index, axis, side_index] = lookup[tuple(adjacent)]
    return neighbor


def build_block_sparse_levelset_3d(
        geometry, *, block_cell_shape=(8, 8, 8), band_width_cells=8,
        periodic_axes=()):
    """Classify one dense reference geometry into exact-band and labeled far blocks."""
    if not isinstance(geometry, FeatureGeometry3D):
        raise TypeError("block-sparse construction requires FeatureGeometry3D")
    block_shape = _as_three_positive_ints(block_cell_shape, "block_cell_shape")
    axes = _periodic_axes(periodic_axes)
    if (isinstance(band_width_cells, (bool, np.bool_))
            or int(band_width_cells) != band_width_cells
            or int(band_width_cells) < 1):
        raise ValueError("band_width_cells must be a positive integer")
    band_width = int(band_width_cells)
    material_ids = tuple(
        int(value) for value in np.unique(geometry.material_id) if int(value) > 0)
    if geometry.material_levelsets is None:
        if len(material_ids) != 1:
            raise ValueError("multi-material sparse geometry requires material level sets")
        material_fields = {material_ids[0]: np.asarray(geometry.phi, dtype=float)}
    else:
        material_fields = {
            int(material_id): np.asarray(field, dtype=float)
            for material_id, field in geometry.material_levelsets.items()}
        if tuple(sorted(material_fields)) != tuple(sorted(material_ids)):
            raise ValueError("material level sets do not match represented solid materials")
    canonical_phi, canonical_owner, material_fields, periodic_difference = (
        _canonicalize_periodic_geometry(geometry, axes, material_fields))

    full_label = np.where(
        canonical_phi >= 0.0, canonical_owner, FAR_GAS)
    canonical_label = _canonical_view(full_label, axes)
    physical_phi = _canonical_view(canonical_phi, axes)
    seed = _interface_seed(canonical_label, physical_phi, axes)
    active_canonical = _dilate(seed, band_width, axes)
    active_node = _expand_periodic_endpoints(
        active_canonical, geometry.phi.shape, axes)

    intervals = np.asarray(geometry.phi.shape) - 1
    grid_shape = tuple(int(np.ceil(intervals[axis] / block_shape[axis]))
                       for axis in range(3))
    blocks = []
    for key in np.ndindex(grid_shape):
        start = tuple(key[axis] * block_shape[axis] for axis in range(3))
        stop = tuple(min((key[axis] + 1) * block_shape[axis], intervals[axis])
                     for axis in range(3))
        selector = tuple(slice(left, right + 1) for left, right in zip(start, stop))
        labels = np.asarray(full_label[selector], dtype=int)
        active = bool(np.any(active_node[selector]))
        unique = np.unique(labels)
        if not active and len(unique) != 1:
            active = True
        if active:
            blocks.append(SparseLevelSetBlock3D(
                key, start, stop, ACTIVE_BLOCK,
                np.asarray(canonical_phi[selector]),
                {material_id: field[selector]
                 for material_id, field in material_fields.items()},
                labels,
            ))
        else:
            blocks.append(SparseLevelSetBlock3D(
                key, start, stop, int(unique[0])))
    return BlockSparseLevelSet3D(
        tuple(int(value) for value in geometry.phi.shape),
        float(geometry.dx), float(geometry.mesh_length_unit_m),
        tuple(float(value) for value in geometry.mesh_origin_m),
        tuple(sorted(material_fields)), block_shape, band_width, axes,
        periodic_difference, tuple(blocks), _neighbor_table(grid_shape, axes))
