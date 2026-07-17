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
    value[0] += 1e-3
    with pytest.raises(ValueError, match="violates the declared extrusion invariance"):
        result.certify_face_field(
            value, relative_tolerance=1e-6, absolute_tolerance=1e-9)
