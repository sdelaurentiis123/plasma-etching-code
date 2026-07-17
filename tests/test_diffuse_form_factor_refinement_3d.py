import numpy as np
import pytest

import petch.diffuse_form_factor_refinement_3d as refinement
from petch.boundary_transport_3d import (
    _diffuse_form_factor_ray_sample_block_3d,
    _diffuse_form_factor_ray_samples_3d,
    estimate_diffuse_form_factors_3d,
    trace_diffuse_form_factor_events_float64_3d,
)
from petch.neutral_radiosity_3d import DiffuseFormFactors3D


def _mixed_escape_geometry():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.8, 0.8, 0.0], [0.9, 0.8, 0.0], [0.8, 0.9, 0.0],
        [0.25, 0.25, 0.4], [0.75, 0.25, 0.4],
        [0.75, 0.75, 0.4], [0.25, 0.75, 0.4],
    ])
    faces = np.asarray([
        [0, 1, 2], [3, 4, 5], [6, 8, 7], [6, 9, 8],
    ])
    centroids = vertices[faces].mean(axis=1)
    normals = np.asarray([
        [0.0, 0.0, 1.0], [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0], [0.0, 0.0, -1.0],
    ])
    return vertices, faces, centroids, normals


def _options(*, periodic=False):
    return {
        "domain_size": (1.0, 1.0, 1.0),
        "periodic_lateral": periodic,
        "ray_offset": 1.0e-5,
        "source_sampling": "triangle_area",
        "visibility_mode": "cellwise_certified",
        "maximum_visibility_wraps": 64,
        "maximum_visibility_replay_wraps": 4096,
        "device": "cpu",
    }


def _estimate(geometry, rays, seed, *, periodic=False):
    vertices, faces, centroids, normals = geometry
    return estimate_diffuse_form_factors_3d(
        vertices, faces, centroids, normals,
        rays_per_face=rays, seed=seed, **_options(periodic=periodic))


def _identity(geometry, factors, seed, epoch="manufactured", *, periodic=False):
    vertices, faces, centroids, normals = geometry
    return refinement.diffuse_form_factor_operator_identity_3d(
        vertices, faces, centroids, normals,
        rays_per_face=factors.rays_per_face, seed=seed,
        operator_identity={"epoch": epoch}, form_factors=factors,
        **_options(periodic=periodic))


def _refine(geometry, base, identity, *, selected=(0,), refined=16,
            epoch="manufactured"):
    return refinement.refine_diffuse_form_factor_rows_nested_3d(
        base, identity, *geometry, selected_source_faces=selected,
        refined_rays_per_face=refined,
        expected_operator_identity={"epoch": epoch})


def _dense(factors):
    result = np.zeros((factors.face_count, factors.face_count))
    result[factors.source_face, factors.target_face] = factors.transfer_fraction
    return result


def test_selected_sample_block_is_exact_slice_of_global_nested_rule():
    geometry = _mixed_escape_geometry()
    vertices, faces, centroids, normals = geometry
    full_source, full_origin, full_direction = _diffuse_form_factor_ray_samples_3d(
        vertices, faces, centroids, normals, rays_per_face=16, seed=29,
        ray_offset=1.0e-5, source_sampling="triangle_area")
    source, origin, direction = _diffuse_form_factor_ray_sample_block_3d(
        vertices, faces, centroids, normals, source_faces=[0, 2],
        sobol_index_start=8, sobol_index_stop=16, seed=29,
        ray_offset=1.0e-5, source_sampling="triangle_area")
    full_origin = full_origin.reshape(4, 16, 3)
    full_direction = full_direction.reshape(4, 16, 3)

    assert np.array_equal(source, np.repeat([0, 2], 8))
    assert np.array_equal(origin.reshape(2, 8, 3), full_origin[[0, 2], 8:])
    assert np.array_equal(direction.reshape(2, 8, 3), full_direction[[0, 2], 8:])
    assert np.array_equal(full_source.reshape(4, 16)[0], np.zeros(16, dtype=int))


@pytest.mark.parametrize("seed", (3, 5))
def test_selected_refinement_matches_global_rows_and_preserves_untouched_rows(seed):
    geometry = _mixed_escape_geometry()
    base = _estimate(geometry, 8, seed)
    fine = _estimate(geometry, 16, seed)
    receipt = _refine(geometry, base, _identity(geometry, base, seed))

    assert np.array_equal(_dense(receipt.form_factors)[0], _dense(fine)[0])
    assert receipt.form_factors.escape_fraction[0] == fine.escape_fraction[0]
    assert np.array_equal(_dense(receipt.form_factors)[1:], _dense(base)[1:])
    assert np.array_equal(
        receipt.form_factors.escape_fraction[1:], base.escape_fraction[1:])
    assert np.array_equal(receipt.row_ray_count, [16, 8, 8, 8])
    assert receipt.seed == seed
    assert receipt.traced_ray_count == 8
    classified = receipt.escape_ray_count + np.bincount(
        receipt.form_factors.source_face, weights=receipt.transfer_ray_count,
        minlength=4).astype(np.int64)
    assert np.array_equal(classified, receipt.row_ray_count)
    assert 0 < receipt.open_escape_count < receipt.traced_ray_count


def _periodic_shared_edge_geometry():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
        [0.1, 0.1, 0.5], [0.2, 0.1, 0.5],
        [0.2, 0.2, 0.5], [0.1, 0.2, 0.5],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6]])
    normals = np.asarray([
        [0.0, 0.0, 1.0], [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0], [0.0, 0.0, -1.0],
    ])
    centroids = vertices[faces].mean(axis=1)
    _source, origin, direction = _diffuse_form_factor_ray_sample_block_3d(
        vertices, faces, centroids, normals, source_faces=[0],
        sobol_index_start=8, sobol_index_stop=16, seed=1,
        ray_offset=1.0e-5, source_sampling="triangle_area")
    shared_ray = 2
    distance = (0.5 - origin[shared_ray, 2]) / direction[shared_ray, 2]
    impact = np.mod(
        origin[shared_ray, :2] + distance * direction[shared_ray, :2], 1.0)
    half_width = 0.08
    vertices[4:8] = np.asarray([
        [impact[0] - half_width, impact[1] - half_width, 0.5],
        [impact[0] + half_width, impact[1] - half_width, 0.5],
        [impact[0] + half_width, impact[1] + half_width, 0.5],
        [impact[0] - half_width, impact[1] + half_width, 0.5],
    ])
    return vertices, faces, vertices[faces].mean(axis=1), normals, shared_ray


def test_periodic_shared_edge_row_matches_global_certified_refinement():
    vertices, faces, centroids, normals, shared_ray = (
        _periodic_shared_edge_geometry())
    geometry = vertices, faces, centroids, normals
    _source, origin, direction = _diffuse_form_factor_ray_sample_block_3d(
        *geometry, source_faces=[0], sobol_index_start=8,
        sobol_index_stop=16, seed=1, ray_offset=1.0e-5,
        source_sampling="triangle_area")
    reference = trace_diffuse_form_factor_events_float64_3d(
        origin[shared_ray:shared_ray + 1], direction[shared_ray:shared_ray + 1],
        vertices, faces, normals, domain_size=(1.0, 1.0, 1.0),
        periodic_lateral=True, maximum_wraps=64)
    assert reference.hit_face[0] in (2, 3)
    assert reference.wrap_count[0] >= 1

    base = _estimate(geometry, 8, 1, periodic=True)
    fine = _estimate(geometry, 16, 1, periodic=True)
    receipt = _refine(
        geometry, base, _identity(geometry, base, 1, epoch="edge", periodic=True),
        epoch="edge")
    assert np.array_equal(_dense(receipt.form_factors)[0], _dense(fine)[0])
    assert receipt.form_factors.escape_fraction[0] == fine.escape_fraction[0]
    assert np.array_equal(_dense(receipt.form_factors)[1:], _dense(base)[1:])
    assert receipt.maximum_wrap_count >= 1


def test_receipt_is_immutable_and_binds_counts_mesh_seed_levels_and_sources():
    geometry = _mixed_escape_geometry()
    base = _estimate(geometry, 8, 3)
    identity = _identity(geometry, base, 3)
    receipt = _refine(geometry, base, identity)

    assert len(receipt.sha256) == len(identity.sha256) == 64
    assert receipt.construction_identity["base_operator_identity_sha256"] == (
        identity.sha256)
    assert receipt.construction_identity["sobol_index_interval"] == (8, 16)
    assert receipt.construction_identity["selected_source_face"] == (0,)
    assert receipt.construction_identity["count_merge"].startswith("integer")
    assert not receipt.row_ray_count.flags.writeable
    assert not receipt.transfer_ray_count.flags.writeable
    with pytest.raises((TypeError, ValueError)):
        receipt.construction_identity["seed"] = 99
    with pytest.raises(ValueError):
        receipt.row_ray_count[0] = 8


@pytest.mark.parametrize(
    "selected,refined,cap,match",
    [
        ((0,), 8, 0.25, "increasing powers"),
        ((0,), 12, 0.25, "increasing powers"),
        ((0, 0), 16, 0.25, "duplicate or out of range"),
        ((4,), 16, 0.25, "duplicate or out of range"),
        ((0, 1), 16, 0.25, "25% cap"),
        ((0,), 16, 0.5, "25% cap"),
    ],
)
def test_refinement_refuses_invalid_levels_sources_and_fraction(
        selected, refined, cap, match):
    geometry = _mixed_escape_geometry()
    base = _estimate(geometry, 8, 3)
    identity = _identity(geometry, base, 3)
    with pytest.raises(refinement.DiffuseFormFactorRefinementRefusal, match=match):
        refinement.refine_diffuse_form_factor_rows_nested_3d(
            base, identity, *geometry, selected_source_faces=selected,
            refined_rays_per_face=refined,
            expected_operator_identity={"epoch": "manufactured"},
            selected_source_fraction_cap=cap)


def test_refinement_refuses_operator_mesh_and_base_sample_mismatch():
    geometry = _mixed_escape_geometry()
    base = _estimate(geometry, 8, 3)
    identity = _identity(geometry, base, 3)
    with pytest.raises(
            refinement.DiffuseFormFactorRefinementRefusal,
            match="caller operator identity"):
        _refine(geometry, base, identity, epoch="another-operator")

    changed = list(geometry)
    changed[0] = geometry[0].copy()
    changed[0][0, 0] += 1.0e-6
    changed[2] = changed[0][changed[1]].mean(axis=1)
    with pytest.raises(
            refinement.DiffuseFormFactorRefinementRefusal,
            match="mesh or sampling operator"):
        _refine(tuple(changed), base, identity)

    other_base = _estimate(geometry, 8, 5)
    with pytest.raises(
            refinement.DiffuseFormFactorRefinementRefusal,
            match="mesh or sampling operator"):
        _refine(geometry, other_base, identity)


def test_refinement_refuses_noninteger_base_and_incomplete_replay(monkeypatch):
    geometry = _mixed_escape_geometry()
    invalid = DiffuseFormFactors3D(
        4, [0], [2], [0.3], [0.7, 1.0, 1.0, 1.0], 8)
    invalid_identity = _identity(geometry, invalid, 3)
    with pytest.raises(
            refinement.DiffuseFormFactorRefinementRefusal,
            match="integer ray counts"):
        _refine(geometry, invalid, invalid_identity)

    base = _estimate(geometry, 8, 3)
    identity = _identity(geometry, base, 3)

    def refuse(*_args, **_kwargs):
        raise RuntimeError("exact replay encountered periodic-wrap budget exhaustion")

    monkeypatch.setattr(
        refinement, "trace_diffuse_form_factor_events_cellwise_certified_3d",
        refuse)
    with pytest.raises(RuntimeError, match="exact replay.*budget exhaustion"):
        _refine(geometry, base, identity)
