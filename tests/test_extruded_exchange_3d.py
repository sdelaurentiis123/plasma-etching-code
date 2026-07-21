import numpy as np
import pytest

from petch.extruded_exchange_3d import build_extruded_triangle_exchange_3d


def _extrude_segments(segments, normals, *, y_planes=(0.0, 1.0)):
    vertices = []
    faces = []
    face_normals = []
    for segment, normal in zip(segments, normals):
        for lower, upper in zip(y_planes[:-1], y_planes[1:]):
            start, stop = segment
            offset = len(vertices)
            vertices.extend([
                [start[0], lower, start[1]],
                [stop[0], lower, stop[1]],
                [start[0], upper, start[1]],
                [stop[0], upper, stop[1]],
            ])
            faces.extend([
                [offset, offset + 1, offset + 2],
                [offset + 1, offset + 3, offset + 2],
            ])
            face_normals.extend([
                [normal[0], 0.0, normal[1]],
                [normal[0], 0.0, normal[1]],
            ])
    return np.asarray(vertices), np.asarray(faces), np.asarray(face_normals)


def _rectangle():
    segments = np.array([
        [[0.0, 0.0], [1.0, 0.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [[0.0, 1.0], [1.0, 1.0]],
        [[0.0, 0.0], [0.0, 1.0]],
    ])
    normals = np.array([
        [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0],
    ])
    return segments, normals


def test_extruded_rectangle_groups_triangle_strips_and_closes_enclosure():
    segments, normals = _rectangle()
    vertices, faces, face_normals = _extrude_segments(
        segments, normals, y_planes=(0.0, 0.4, 1.0))
    result = build_extruded_triangle_exchange_3d(
        vertices, faces, face_normals, extrusion_length=1.0,
        exchange_relative_tolerance=1e-8)
    assert result.group_count == 4
    assert result.strip_count == 2
    assert result.face_count == 16
    assert result.maximum_group_area_relative_error < 2e-15
    assert result.maximum_area_reciprocity_error < 2e-15
    assert np.allclose(result.group_area, 1.0, rtol=0.0, atol=2e-15)
    assert np.allclose(
        result.line_exchange.transfer_fraction.sum(axis=1), 1.0,
        rtol=0.0, atol=2e-14)
    assert np.allclose(result.form_factors.escape_fraction, 0.0, atol=2e-14)


def test_triangle_repartition_does_not_change_group_operator():
    segments, normals = _rectangle()
    coarse = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 1.0)),
        extrusion_length=1.0)
    split = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 0.2, 0.6, 1.0)),
        extrusion_length=1.0)
    assert np.allclose(
        coarse.line_exchange.exchange_length,
        split.line_exchange.exchange_length, rtol=0.0, atol=2e-15)
    assert coarse.line_exchange.escape_fraction == pytest.approx(
        split.line_exchange.escape_fraction, abs=2e-15)


def test_expanded_face_operator_preserves_area_reciprocity():
    segments, normals = _rectangle()
    result = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 0.25, 1.0)),
        extrusion_length=1.0)
    factors = result.form_factors
    lookup = {
        (int(source), int(target)): float(fraction)
        for source, target, fraction in zip(
            factors.source_face, factors.target_face, factors.transfer_fraction)
    }
    for (source, target), fraction in lookup.items():
        assert result.face_area[source] * fraction == pytest.approx(
            result.face_area[target] * lookup[(target, source)], abs=2e-15)


def test_nonextruded_triangle_refuses_instead_of_projecting():
    segments, normals = _rectangle()
    vertices, faces, face_normals = _extrude_segments(segments, normals)
    vertices[0, 0] += 0.02
    with pytest.raises(ValueError, match="not an extrusion"):
        build_extruded_triangle_exchange_3d(
            vertices, faces, face_normals, extrusion_length=1.0,
            geometry_tolerance=1e-9)


def test_extrusion_field_certification_is_area_weighted_and_fail_closed():
    segments, normals = _rectangle()
    result = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 0.25, 1.0)),
        extrusion_length=1.0)
    value = 2.0 + result.face_group_index.astype(float)
    mean = result.certify_face_field(
        value, relative_tolerance=0.0, absolute_tolerance=0.0)
    assert np.array_equal(mean, np.arange(2.0, 6.0))
    first_group = np.flatnonzero(result.face_group_index == 0)
    first_strip = first_group[
        result.face_strip_index[first_group] == result.face_strip_index[first_group[0]]]
    value[first_strip] += 1e-3
    with pytest.raises(ValueError, match="violates the declared extrusion invariance"):
        result.certify_face_field(
            value, relative_tolerance=1e-6, absolute_tolerance=1e-9)


def test_triangle_partition_variation_is_averaged_without_hiding_strip_variation():
    segments, normals = _rectangle()
    result = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 0.5, 1.0)),
        extrusion_length=1.0)
    value = 2.0 + result.face_group_index.astype(float)
    member = np.flatnonzero(
        (result.face_group_index == 0) & (result.face_strip_index == 0))
    assert member.size == 2
    value[member[0]] -= 0.2
    value[member[1]] += 0.2
    mean = result.certify_face_field(
        value, relative_tolerance=0.0, absolute_tolerance=1e-14)
    assert mean[0] == pytest.approx(2.0, abs=1e-15)


def test_geometry_extrusion_projection_flattens_noise_and_guards_3d():
    from petch.extruded_exchange_3d import project_geometry_to_extrusion
    from petch.feature_geometry_state_3d import FeatureGeometry3D

    x = np.arange(8) * 0.01
    z = np.arange(10) * 0.01
    phi0 = (0.05 - z)[None, None, :] + 0.0 * x[:, None, None]
    phi = np.repeat(np.repeat(phi0, 3, axis=1), 8, axis=0)[:8, :3, :]
    rng = np.random.default_rng(11)
    noise = 1e-7 * rng.standard_normal(phi.shape)
    noisy = phi + noise
    levels = {1: noisy.copy()}
    material = np.where(noisy >= 0.0, 1, 0)
    geometry = FeatureGeometry3D(noisy, material, 0.01, 1e-6, material_levelsets=levels)
    projected, deviation = project_geometry_to_extrusion(geometry)
    assert 0.0 < deviation < 1e-6
    projected_phi = np.asarray(projected.phi)
    assert np.max(np.abs(projected_phi - projected_phi[:, :1, :])) == 0.0
    union = np.asarray(projected.material_levelsets[1])
    assert not np.any((union >= 0.0) != (projected_phi >= 0.0))

    coarse = phi + 0.5 * 0.01 * np.sin(
        np.arange(3)[None, :, None] * 2.0)  # genuine 3-D variation
    material3d = np.where(coarse >= 0.0, 1, 0)
    geometry3d = FeatureGeometry3D(
        coarse, material3d, 0.01, 1e-6, material_levelsets={1: coarse.copy()})
    with pytest.raises(ValueError, match="exceeds the projection guard"):
        project_geometry_to_extrusion(geometry3d)


def test_field_projection_floors_shadowed_group_scale_but_guards_bright_groups():
    segments, normals = _rectangle()
    result = build_extruded_triangle_exchange_3d(
        *_extrude_segments(segments, normals, y_planes=(0.0, 0.5, 1.0)),
        extrusion_length=1.0)
    value = np.full(result.face_count, 100.0)
    member = np.flatnonzero(
        (result.face_group_index == 0) & (result.face_strip_index == 0))
    other = np.flatnonzero(
        (result.face_group_index == 0) & (result.face_strip_index == 1))
    assert member.size and other.size
    # Deeply shadowed group: ~1e-4 of peak with 50% strip-to-strip spread passes under
    # the floored scale (absolute noise is negligible against the field).
    value[member] = 0.012
    value[other] = 0.008
    mean, variation = result.project_face_field(value, relative_guard=0.05)
    assert variation < 0.05
    # The same 50% spread at significant flux still refuses.
    value[member] = 150.0
    value[other] = 50.0
    with pytest.raises(ValueError, match="mean-field projection guard"):
        result.project_face_field(value, relative_guard=0.05)
