"""Adversarial parity checks for the diffuse form-factor visibility fast path.

These tests compare the legacy float32 Warp classifier with the cell-by-cell float64
authority on deliberately small hard-visibility cases.  Known differences are asserted as
diagnostics rather than being hidden behind a looser geometric tolerance.
"""

import numpy as np
import pytest

import petch.boundary_transport_3d as boundary_transport_3d


def _trace_pair(
        verts, faces, normals, origin, direction, *, domain=(1.0, 1.0, 1.0),
        periodic=False, maximum_wraps=1024):
    origin = np.asarray(origin, dtype=float)
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    domain = np.asarray(domain, dtype=float)
    fast = boundary_transport_3d._trace_diffuse_form_factor_events_warp_3d(
        np.asarray(verts, dtype=float), np.asarray(faces, dtype=int), origin,
        direction, domain, periodic, "cpu")
    reference = boundary_transport_3d.trace_diffuse_form_factor_events_float64_3d(
        origin, direction, np.asarray(verts, dtype=float), np.asarray(faces, dtype=int),
        np.asarray(normals, dtype=float), domain_size=domain,
        periodic_lateral=periodic, maximum_wraps=maximum_wraps)
    return fast, reference


def _shared_square():
    verts = np.asarray([
        [0.0, 0.0, 0.5],
        [1.0, 0.0, 0.5],
        [1.0, 1.0, 0.5],
        [0.0, 1.0, 0.5],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]])
    normals = np.tile([0.0, 0.0, -1.0], (2, 1))
    return verts, faces, normals


def test_fast_visibility_matches_float64_on_exact_shared_edge_at_normal_incidence():
    verts, faces, normals = _shared_square()
    fast, reference = _trace_pair(
        verts, faces, normals, [[0.5, 0.5, 0.1]], [[0.0, 0.0, 1.0]])

    assert np.array_equal(reference.termination, [1])
    assert fast[0] in (0, 1)
    assert reference.hit_face[0] in (0, 1)
    assert fast[0] == reference.hit_face[0]


def test_fast_visibility_matches_float64_on_both_sides_of_shared_edge():
    verts, faces, normals = _shared_square()
    edge_offset = 2.0e-5
    origin = np.asarray([
        [0.5, 0.5 - edge_offset, 0.1],
        [0.5, 0.5 + edge_offset, 0.1],
    ])
    direction = np.tile([0.0, 0.0, 1.0], (2, 1))
    fast, reference = _trace_pair(verts, faces, normals, origin, direction)

    assert np.array_equal(reference.termination, [1, 1])
    assert np.array_equal(fast, reference.hit_face)
    assert set(fast.tolist()) == {0, 1}


def test_production_grazing_shared_edge_is_cross_examined_against_float64():
    # These normal/direction values are the production shared-edge regression preserved in
    # test_charged_surface_cascade_3d.py.  Construct a square whose diagonal passes through the
    # impact point, then launch along that nearly tangent direction toward the diagonal.
    normal = np.asarray([-0.7071067690849304, 0.0, 0.7071067690849304])
    normal /= np.linalg.norm(normal)
    direction = np.asarray([
        -0.042131607648075275, -0.9982234525891459, -0.04212916255984473])
    direction /= np.linalg.norm(direction)
    along_edge = np.asarray([0.0, 1.0, 0.0])
    across_edge = np.cross(normal, along_edge)
    across_edge /= np.linalg.norm(across_edge)
    impact = np.asarray([0.5, 0.5, 0.8])
    half_width = 0.4
    verts = np.asarray([
        impact - half_width * along_edge - half_width * across_edge,
        impact + half_width * along_edge - half_width * across_edge,
        impact + half_width * along_edge + half_width * across_edge,
        impact - half_width * along_edge + half_width * across_edge,
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]])
    normals = np.repeat(normal[None, :], 2, axis=0)
    origin = impact - 0.25 * direction
    fast, reference = _trace_pair(
        verts, faces, normals, origin[None, :], direction[None, :],
        domain=(1.0, 1.0, 1.5))

    assert np.array_equal(reference.termination, [1])
    assert reference.hit_face[0] in (0, 1)
    # Keep exact fast/reference ownership visible.  If float32 misses this near-tangent seam,
    # this strict assertion fails rather than reclassifying the ray as an escape.
    assert fast[0] == reference.hit_face[0]


def test_fast_visibility_matches_float64_after_one_periodic_wrap():
    verts = np.asarray([
        [0.0, 0.0, 0.3],
        [0.35, 0.0, 0.3],
        [0.35, 1.0, 0.3],
        [0.0, 1.0, 0.3],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]])
    normals = np.tile([0.0, 0.0, -1.0], (2, 1))
    fast, reference = _trace_pair(
        verts, faces, normals, [[0.9, 0.25, 0.1]], [[1.0, 0.0, 1.0]],
        periodic=True, maximum_wraps=8)

    assert np.array_equal(reference.termination, [1])
    assert reference.wrap_count[0] == 1
    assert fast[0] == reference.hit_face[0]


def test_fast_open_escape_matches_float64_open_top_event():
    verts = np.asarray([
        [0.0, 0.0, 0.3],
        [0.2, 0.0, 0.3],
        [0.0, 0.2, 0.3],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    fast, reference = _trace_pair(
        verts, faces, normals, [[0.8, 0.8, 0.2]], [[0.0, 0.0, 1.0]],
        periodic=True)

    assert np.array_equal(fast, [-1])
    assert np.array_equal(reference.termination, [2])
    assert np.array_equal(reference.hit_face, [-1])
    assert reference.hit_position[0, 2] == 1.0


def test_fast_miss_code_cannot_distinguish_periodic_exhaustion_from_escape():
    # The fast API has only a face index, so -1 conflates a certified open escape with exhausting
    # its fixed internal wrap loop.  The float64 authority preserves this as refusal event 3.
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    fast, reference = _trace_pair(
        verts, faces, normals, [[0.5, 0.5, 0.5]], [[1.0, 0.0, 0.0]],
        periodic=True, maximum_wraps=3)

    assert np.array_equal(fast, [-1])
    assert np.array_equal(reference.termination, [3])
    assert np.array_equal(reference.hit_face, [-1])
    assert reference.wrap_count[0] == 4


def test_replay_hardened_path_recovers_a_fast_shared_edge_miss(monkeypatch):
    verts, faces, normals = _shared_square()
    monkeypatch.setattr(
        boundary_transport_3d,
        "_trace_diffuse_form_factor_events_warp_3d",
        lambda *args, **kwargs: np.asarray([-1], dtype=int))

    result = boundary_transport_3d.trace_diffuse_form_factor_events_replay_hardened_3d(
        [[0.5, 0.5, 0.1]], [[0.0, 0.0, 1.0]], verts, faces, normals,
        domain_size=(1.0, 1.0, 1.0), device="cpu")

    assert np.array_equal(result.termination, [1])
    assert result.hit_face[0] in (0, 1)
    assert result.replay_count == result.replay_eligible_count == 1
    assert result.recovered_hit_count == 1
    assert result.open_escape_count == 0


def test_replay_hardened_path_refuses_periodic_wrap_exhaustion():
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])

    with pytest.raises(RuntimeError, match="periodic-wrap budget"):
        boundary_transport_3d.trace_diffuse_form_factor_events_replay_hardened_3d(
            [[0.5, 0.5, 0.5]], [[1.0, 0.0, 0.0]], verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
            maximum_wraps=3, device="cpu")


def test_cellwise_certified_path_recovers_a_float32_exhaustion(monkeypatch):
    verts, faces, normals = _shared_square()
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [4]))

    result = (
        boundary_transport_3d
        .trace_diffuse_form_factor_events_cellwise_certified_3d(
            [[0.5, 0.5, 0.1]], [[0.0, 0.0, 1.0]], verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), maximum_wraps=8, device="cpu"))

    assert np.array_equal(result.termination, [1])
    assert result.hit_face[0] in (0, 1)
    assert result.replay_count == 1
    assert result.recovered_hit_count == 1
    assert result.maximum_wrap_count == 0


def test_cellwise_certified_path_refuses_incomplete_exact_replay(monkeypatch):
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [4]))

    with pytest.raises(RuntimeError, match="exact replay.*periodic-wrap budget"):
        (boundary_transport_3d
         .trace_diffuse_form_factor_events_cellwise_certified_3d(
             [[0.5, 0.5, 0.5]], [[1.0, 0.0, 0.0]], verts, faces, normals,
             domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
             maximum_wraps=3, device="cpu"))


def test_cellwise_certified_path_uses_a_separate_bounded_exact_replay_horizon(
        monkeypatch):
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [2]))
    direction = np.asarray([[1.0, 0.0, 0.2]])
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)

    result = (
        boundary_transport_3d
        .trace_diffuse_form_factor_events_cellwise_certified_3d(
            [[0.1, 0.5, 0.5]], direction, verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
            maximum_wraps=1, maximum_exact_replay_wraps=8, device="cpu"))

    assert np.array_equal(result.termination, [2])
    assert np.array_equal(result.hit_face, [-1])
    assert result.replay_count == 1
    assert result.maximum_wrap_count == 2
    assert result.derived_horizon_extension_count == 0
    assert result.initial_maximum_wraps == 1
    assert result.final_maximum_wraps == 8


def test_cellwise_certified_path_derives_open_top_horizon_for_grazing_ray(
        monkeypatch):
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [4]))
    direction = np.asarray([[1.0, 0.0, 1.0e-3]])
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)

    result = (
        boundary_transport_3d
        .trace_diffuse_form_factor_events_cellwise_certified_3d(
            [[0.1, 0.5, 0.5]], direction, verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
            maximum_wraps=3, device="cpu"))

    assert np.array_equal(result.termination, [2])
    assert np.array_equal(result.hit_face, [-1])
    assert result.replay_count == 1
    assert result.derived_horizon_extension_count == 1
    assert result.maximum_wrap_count == 500
    assert result.initial_maximum_wraps == 3
    assert result.final_maximum_wraps == 504 + 64 + 504 // 16

    # Preserve the exact production lineage that exposed the former guessed 1024-wrap cap.  Its
    # true float64 trace escaped after 1113 wraps; the conservative geometric authority is 1117.
    production_origin = np.asarray([
        0.020309919291093343, 0.007027431235114469, 2.6508378871010123])
    production_direction = np.asarray([
        0.9363883881477437, -0.3509498266169299, 0.0033174899574309746])
    assert boundary_transport_3d._derived_vertical_domain_wrap_horizon_3d(
        production_origin, production_direction, np.asarray([0.13, 0.02, 2.8])) == 1117


def test_cellwise_certified_path_derives_lower_domain_horizon_for_grazing_ray(
        monkeypatch):
    verts, faces, _ = _shared_square()
    normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [4]))
    direction = np.asarray([[1.0, 0.0, -1.0e-3]])
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)

    result = (
        boundary_transport_3d
        .trace_diffuse_form_factor_events_cellwise_certified_3d(
            [[0.1, 0.5, 0.6]], direction, verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
            maximum_wraps=3, device="cpu"))

    assert np.array_equal(result.termination, [1])
    assert result.hit_face[0] in (0, 1)
    assert result.replay_count == 1
    assert result.derived_horizon_extension_count == 1
    assert result.maximum_wrap_count == 100
    assert result.initial_maximum_wraps == 3
    assert result.final_maximum_wraps == 604 + 64 + 604 // 16

    # Preserve the second production lineage.  Its exact replay hit a real face after 1110 wraps;
    # the conservative lower-domain integrity horizon is finite and therefore recoverable.
    production_origin = np.asarray([
        0.108438417364833, 0.004037909820028365, 2.134689932283348])
    production_direction = np.asarray([
        -0.00391740694902086, -0.9999561259207863, -0.008508828138368412])
    horizon = boundary_transport_3d._derived_vertical_domain_wrap_horizon_3d(
        production_origin, production_direction, np.asarray([0.13, 0.02, 2.8]))
    assert horizon == 12556
    assert horizon > 1110


def test_float64_authority_exposes_solid_facing_intersection_as_refusal():
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    reference = boundary_transport_3d.trace_diffuse_form_factor_events_float64_3d(
        [[0.5, 0.5, 0.1]], [[0.0, 0.0, 1.0]], verts, faces,
        solid_facing_normals, domain_size=(1.0, 1.0, 1.0))

    assert np.array_equal(reference.termination, [4])
    assert reference.hit_face[0] in (0, 1)
    assert reference.hit_cosine[0] == -1.0

    with pytest.raises(RuntimeError, match="solid-facing hard intersection"):
        boundary_transport_3d.trace_diffuse_form_factor_events_replay_hardened_3d(
            [[0.5, 0.5, 0.1]], [[0.0, 0.0, 1.0]], verts, faces,
            solid_facing_normals, domain_size=(1.0, 1.0, 1.0), device="cpu")


def test_cellwise_warp_candidate_matches_explicit_float64_event_codes():
    verts, faces, normals = _shared_square()
    origin = np.asarray([
        [0.5, 0.5, 0.1],
        [0.25, 0.75, 0.1],
        [1.25, 0.5, 0.1],
    ])
    direction = np.tile([0.0, 0.0, 1.0], (3, 1))
    candidate = boundary_transport_3d.trace_diffuse_form_factor_events_warp_cellwise_3d(
        origin, direction, verts, faces, normals,
        domain_size=(2.0, 2.0, 1.0), device="cpu")
    authority = boundary_transport_3d.trace_diffuse_form_factor_events_float64_3d(
        origin, direction, verts, faces, normals,
        domain_size=(2.0, 2.0, 1.0))

    assert np.array_equal(candidate.termination, authority.termination)
    assert np.array_equal(candidate.hit_face, authority.hit_face)

    solid_normals = -normals
    solid_candidate = (
        boundary_transport_3d.trace_diffuse_form_factor_events_warp_cellwise_3d(
            origin[:1], direction[:1], verts, faces, solid_normals,
            domain_size=(2.0, 2.0, 1.0), device="cpu"))
    assert np.array_equal(solid_candidate.termination, [4])


def test_cellwise_certified_recovers_solid_facing_replay_from_source_relaunch():
    # A quadrature point shifted by the declared normal offset can land across an adjacent
    # facet at a sharp concave fold: the exact replay then sees a back-face and refuses.
    # The bounded recovery relaunches only that ray from its support-face centroid origin
    # with the direction unchanged; the receipt records the count and distance.
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    origin = np.asarray([[0.5, 0.5, 0.1]])
    direction = np.asarray([[0.0, 0.0, 1.0]])
    relaunch = np.asarray([[0.5, 0.5, 0.6]])

    with pytest.raises(RuntimeError, match="solid-facing hard intersection"):
        boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
            origin, direction, verts, faces, solid_facing_normals,
            domain_size=(1.0, 1.0, 1.0), device="cpu")

    events = boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
        origin, direction, verts, faces, solid_facing_normals,
        domain_size=(1.0, 1.0, 1.0), source_relaunch_origin=relaunch, device="cpu")
    assert np.array_equal(events.termination, [2])
    assert events.open_escape_count == 1
    assert events.source_relaunch_count == 1
    assert events.maximum_source_relaunch_distance == pytest.approx(0.5)


def test_cellwise_certified_second_solid_facing_refusal_remains_fatal():
    # The relaunch is a single bounded correction of the launch point, not a search:
    # if the support-centroid origin also produces a solid-facing intersection, the
    # operator still refuses rather than guessing.
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    origin = np.asarray([[0.5, 0.5, 0.1]])
    direction = np.asarray([[0.0, 0.0, 1.0]])
    still_solid = np.asarray([[0.5, 0.5, 0.05]])

    with pytest.raises(RuntimeError, match="solid-facing hard intersection"):
        boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
            origin, direction, verts, faces, solid_facing_normals,
            domain_size=(1.0, 1.0, 1.0), source_relaunch_origin=still_solid,
            device="cpu")


def test_cellwise_certified_overlap_skip_admits_grazing_exit_by_perpendicular_depth():
    # The interpenetration overlap is a thin slab: a grazing ray can travel many cells
    # ALONG the skin while its launch sits only a fraction of a cell BEHIND the crossed
    # face's plane.  Admission is therefore by perpendicular depth, not path length.
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    origin = np.asarray([[0.1, 0.5, 0.4999]])
    direction = np.asarray([[1.0, 0.0, 1e-3]])
    direction /= np.linalg.norm(direction)

    events = boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
        origin, direction, verts, faces, solid_facing_normals,
        domain_size=(1.0, 1.0, 1.0), overlap_skip_depth_limit=1e-3, device="cpu")
    assert np.array_equal(events.termination, [2])
    assert events.overlap_skip_count == 1
    # perpendicular depth is the 1e-4 skin, far below the ~0.1 path length to the crossing
    assert events.maximum_overlap_skip_depth == pytest.approx(1e-4, rel=1e-6)


def test_cellwise_certified_overlap_skip_refuses_deep_back_face():
    # A back face whose plane lies deeper behind the launch than the declared limit is
    # not an interpenetration artifact; the operator still refuses.
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    origin = np.asarray([[0.5, 0.5, 0.1]])
    direction = np.asarray([[0.0, 0.0, 1.0]])

    with pytest.raises(RuntimeError, match="solid-facing hard intersection"):
        boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
            origin, direction, verts, faces, solid_facing_normals,
            domain_size=(1.0, 1.0, 1.0), overlap_skip_depth_limit=0.05, device="cpu")


def test_cellwise_certified_overlap_exit_admits_by_authority_surface_distance():
    # When the caller supplies the launch's distance from the authority (trilinear)
    # surface, admission uses it directly: a launch within the declared discretization
    # skin exits the artifact even when the crossed plane lies deep along the ray
    # (extended contact slabs traversed at grazing incidence).
    verts, faces, _ = _shared_square()
    solid_facing_normals = np.tile([0.0, 0.0, 1.0], (2, 1))
    origin = np.asarray([[0.5, 0.5, 0.1]])
    direction = np.asarray([[0.0, 0.0, 1.0]])

    events = boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
        origin, direction, verts, faces, solid_facing_normals,
        domain_size=(1.0, 1.0, 1.0), overlap_skip_depth_limit=0.05,
        overlap_exit_authority_distance=np.asarray([0.01]), device="cpu")
    assert np.array_equal(events.termination, [2])
    assert events.overlap_skip_count == 1

    with pytest.raises(RuntimeError, match="solid-facing hard intersection"):
        boundary_transport_3d.trace_diffuse_form_factor_events_cellwise_certified_3d(
            origin, direction, verts, faces, solid_facing_normals,
            domain_size=(1.0, 1.0, 1.0), overlap_skip_depth_limit=0.05,
            overlap_exit_authority_distance=np.asarray([0.2]), device="cpu")


def test_cellwise_certified_retries_exhaustion_when_horizon_undercounts(monkeypatch):
    # The derived horizon's whole-cell seam count carries only a small tie allowance, so a
    # near-grazing ray can exhaust a budget at or above its derived bound by a handful of
    # crossings.  The retry budget is the proof bound plus an explicit margin; the ray must
    # then classify.  (Production lineage: 10 nm base step 362, one ray of 14656 exhausted
    # 1024 wraps with a derived horizon just below the budget.)
    verts = np.asarray([
        [0.0, 0.0, 0.25],
        [0.1, 0.0, 0.25],
        [0.0, 0.1, 0.25],
    ])
    faces = np.asarray([[0, 1, 2]])
    normals = np.asarray([[0.0, 0.0, -1.0]])
    monkeypatch.setattr(
        boundary_transport_3d,
        "trace_diffuse_form_factor_events_warp_cellwise_3d",
        lambda *args, **kwargs:
        boundary_transport_3d.DiffuseFormFactorEventsWarpCellwise3D(
            [-1], [3], [4]))
    # Undercount the horizon: report fewer wraps than the ray truly needs (~500).
    monkeypatch.setattr(
        boundary_transport_3d,
        "_derived_vertical_domain_wrap_horizon_3d",
        lambda origin, direction, domain: 470)
    direction = np.asarray([[1.0, 0.0, 1.0e-3]])
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)

    result = (
        boundary_transport_3d
        .trace_diffuse_form_factor_events_cellwise_certified_3d(
            [[0.1, 0.5, 0.5]], direction, verts, faces, normals,
            domain_size=(1.0, 1.0, 1.0), periodic_lateral=True,
            maximum_wraps=480, device="cpu"))

    # Old behavior refused: horizon (470) <= budget (480) meant no retry.  The margin
    # retry (470 -> 480 + 64 + 30) lets the true ~500-wrap escape classify.
    assert np.array_equal(result.termination, [2])
    assert result.derived_horizon_extension_count == 1
    assert result.final_maximum_wraps == 480 + 64 + 480 // 16
