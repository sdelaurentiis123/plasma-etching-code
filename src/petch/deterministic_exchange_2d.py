"""Deterministic diffuse exchange for surfaces extruded normal to a 2-D section.

For an infinitely extruded diffuse surface, Hottel's crossed-string relation gives the
unobstructed exchange length between two line elements exactly.  This module combines that
relation with deterministic visibility refinement when third line elements may block only part
of the exchange.  It is an authority/cross-check for line-trench mean profiles; it is not a
replacement for the genuinely three-dimensional hard-visibility operator.

The primary stored quantity is the symmetric exchange length ``H_ij``.  Form factors are derived
as ``F_ij = H_ij / L_i``, so reciprocity is true by construction rather than repaired after the
fact.  Open-boundary escape is the unassigned part of each source row.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import heapq

import numpy as np

from .neutral_radiosity_3d import DiffuseFormFactors3D


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


def _cross_2d(left, right):
    return left[..., 0] * right[..., 1] - left[..., 1] * right[..., 0]


def _point_on_segment(point, start, stop, tolerance):
    direction = stop - start
    offset = point - start
    scale = max(float(np.linalg.norm(direction)), 1.0)
    if abs(float(_cross_2d(direction, offset))) > tolerance * scale:
        return False
    return bool(
        np.dot(offset, direction) >= -tolerance * scale
        and np.dot(point - stop, direction) <= tolerance * scale)


def _proper_segment_intersection(first_start, first_stop, second_start, second_stop,
                                 tolerance):
    """Whether a blocker crosses the open interior of the first segment."""
    first_direction = first_stop - first_start
    second_direction = second_stop - second_start
    denominator = float(_cross_2d(first_direction, second_direction))
    scale = max(
        float(np.linalg.norm(first_direction)),
        float(np.linalg.norm(second_direction)), 1.0)
    if abs(denominator) <= tolerance * scale:
        # Collinear overlap blocks only when it occupies a nonzero open interval of the
        # connector. Merely sharing a surface endpoint is not an obstruction.
        if abs(float(_cross_2d(second_start - first_start, first_direction))) > tolerance * scale:
            return False
        length2 = float(np.dot(first_direction, first_direction))
        if length2 <= tolerance * tolerance:
            return False
        lo, hi = sorted((
            float(np.dot(second_start - first_start, first_direction) / length2),
            float(np.dot(second_stop - first_start, first_direction) / length2)))
        return max(lo, tolerance) < min(hi, 1.0 - tolerance)
    delta = second_start - first_start
    first_parameter = float(_cross_2d(delta, second_direction) / denominator)
    second_parameter = float(_cross_2d(delta, first_direction) / denominator)
    return bool(
        tolerance < first_parameter < 1.0 - tolerance
        and -tolerance <= second_parameter <= 1.0 + tolerance)


def _convex_hull(points):
    unique = sorted(set(map(tuple, np.asarray(points, dtype=float))))
    if len(unique) <= 1:
        return np.asarray(unique, dtype=float)

    def half(sequence):
        output = []
        for point in sequence:
            point = np.asarray(point, dtype=float)
            while (len(output) >= 2
                   and _cross_2d(output[-1] - output[-2], point - output[-1]) <= 0.0):
                output.pop()
            output.append(point)
        return output

    lower = half(unique)
    upper = half(reversed(unique))
    return np.asarray(lower[:-1] + upper[:-1], dtype=float)


def _point_in_convex_polygon(point, polygon, tolerance):
    if len(polygon) < 3:
        return any(_point_on_segment(point, polygon[index], polygon[(index + 1) % len(polygon)],
                                     tolerance)
                   for index in range(len(polygon)))
    sign = []
    for index in range(len(polygon)):
        value = float(_cross_2d(
            polygon[(index + 1) % len(polygon)] - polygon[index],
            point - polygon[index]))
        if abs(value) > tolerance:
            sign.append(np.sign(value))
    return not sign or all(value == sign[0] for value in sign)


def _point_strictly_in_convex_polygon(point, polygon, tolerance):
    if len(polygon) < 3:
        return False
    values = np.asarray([
        _cross_2d(
            polygon[(index + 1) % len(polygon)] - polygon[index],
            point - polygon[index])
        for index in range(len(polygon))], dtype=float)
    return bool(
        np.all(values > tolerance) or np.all(values < -tolerance))


def _transverse_segment_intersection(first_start, first_stop, second_start, second_stop,
                                     tolerance):
    first_direction = first_stop - first_start
    second_direction = second_stop - second_start
    denominator = float(_cross_2d(first_direction, second_direction))
    scale = max(
        float(np.linalg.norm(first_direction)),
        float(np.linalg.norm(second_direction)), 1.0)
    if abs(denominator) <= tolerance * scale:
        return False
    delta = second_start - first_start
    first_parameter = float(_cross_2d(delta, second_direction) / denominator)
    second_parameter = float(_cross_2d(delta, first_direction) / denominator)
    return bool(
        tolerance < first_parameter < 1.0 - tolerance
        and tolerance < second_parameter < 1.0 - tolerance)


def _blocker_may_cross_connector_hull(hull, lower, upper, blocker_start, blocker_stop,
                                      tolerance):
    """Whether a blocker can cross any connector contained by a precomputed hull."""
    blocker_lower = np.minimum(blocker_start, blocker_stop)
    blocker_upper = np.maximum(blocker_start, blocker_stop)
    if np.any(blocker_upper < lower) or np.any(blocker_lower > upper):
        return False
    midpoint = 0.5 * (blocker_start + blocker_stop)
    if (_point_strictly_in_convex_polygon(blocker_start, hull, tolerance)
            or _point_strictly_in_convex_polygon(blocker_stop, hull, tolerance)
            or _point_strictly_in_convex_polygon(midpoint, hull, tolerance)):
        return True
    return any(
        _transverse_segment_intersection(
            blocker_start, blocker_stop, hull[index], hull[(index + 1) % len(hull)], tolerance)
        for index in range(len(hull)))


def _candidate_blockers(first, second, segments, excluded, tolerance, *, candidates=None,
                        segment_lower=None, segment_upper=None):
    """Return the blockers that can intersect the convex connector family.

    Descendant exchange cells inherit their parent's candidate set.  A blocker excluded by a
    parent hull cannot re-enter a child hull, which makes adaptive shadow refinement scale with
    the local occluders instead of repeatedly scanning the entire surface.
    """
    hull = _convex_hull([first[0], first[1], second[0], second[1]])
    if len(hull) == 0:
        return ()
    lower = np.min(hull, axis=0) - tolerance
    upper = np.max(hull, axis=0) + tolerance
    if segment_lower is None:
        segment_lower = np.min(segments, axis=1)
    if segment_upper is None:
        segment_upper = np.max(segments, axis=1)
    if candidates is None:
        candidate = np.arange(len(segments), dtype=int)
    else:
        candidate = np.asarray(candidates, dtype=int)
    if candidate.size:
        overlap = np.all(
            (segment_upper[candidate] >= lower)
            & (segment_lower[candidate] <= upper), axis=1)
        candidate = candidate[overlap]
    if excluded:
        candidate = candidate[
            ~np.isin(candidate, np.fromiter(excluded, dtype=int))]
    if candidate.size == 0:
        return ()

    # Test every candidate against the same connector hull in one vectorized batch.  This is the
    # dominant production path: a trench contains O(N^2) surface pairs, so invoking Python once
    # per possible blocker would turn an inexpensive geometry build into an accidental O(N^3)
    # interpreter loop.
    blocker_start = segments[candidate, 0]
    blocker_stop = segments[candidate, 1]
    point = np.stack((blocker_start, blocker_stop, 0.5 * (blocker_start + blocker_stop)), axis=1)
    inside = np.zeros(candidate.size, dtype=bool)
    if len(hull) >= 3:
        edge_start = hull
        edge_direction = np.roll(hull, -1, axis=0) - hull
        offset = point[:, :, None, :] - edge_start[None, None, :, :]
        signed = (
            edge_direction[None, None, :, 0] * offset[..., 1]
            - edge_direction[None, None, :, 1] * offset[..., 0])
        point_inside = (
            np.all(signed > tolerance, axis=2)
            | np.all(signed < -tolerance, axis=2))
        inside = np.any(point_inside, axis=1)

    edge_start = hull
    edge_direction = np.roll(hull, -1, axis=0) - hull
    blocker_direction = blocker_stop - blocker_start
    denominator = (
        blocker_direction[:, None, 0] * edge_direction[None, :, 1]
        - blocker_direction[:, None, 1] * edge_direction[None, :, 0])
    scale = np.maximum(
        np.maximum(
            np.linalg.norm(blocker_direction, axis=1)[:, None],
            np.linalg.norm(edge_direction, axis=1)[None, :]),
        1.0)
    transverse = np.abs(denominator) > tolerance * scale
    safe_denominator = np.where(transverse, denominator, 1.0)
    delta = edge_start[None, :, :] - blocker_start[:, None, :]
    blocker_parameter = (
        delta[..., 0] * edge_direction[None, :, 1]
        - delta[..., 1] * edge_direction[None, :, 0]) / safe_denominator
    edge_parameter = (
        delta[..., 0] * blocker_direction[:, None, 1]
        - delta[..., 1] * blocker_direction[:, None, 0]) / safe_denominator
    crosses_boundary = np.any(
        transverse
        & (blocker_parameter > tolerance)
        & (blocker_parameter < 1.0 - tolerance)
        & (edge_parameter > tolerance)
        & (edge_parameter < 1.0 - tolerance),
        axis=1)
    return tuple(map(int, candidate[inside | crosses_boundary]))


def _connection_visible(start, stop, blockers, excluded, tolerance, candidates=None):
    return _connection_blocker(
        start, stop, blockers, excluded, tolerance, candidates) is None


def _connection_blocker(start, stop, blockers, excluded, tolerance, candidates=None):
    iterator = range(len(blockers)) if candidates is None else candidates
    for index in iterator:
        index = int(index)
        if index in excluded:
            continue
        blocker = blockers[index]
        if _proper_segment_intersection(start, stop, blocker[0], blocker[1], tolerance):
            return int(index)
    return None


def unobstructed_crossed_string_exchange_2d(first, second):
    """Return the exact symmetric exchange length for two unobstructed line elements."""
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if first.shape != (2, 2) or second.shape != (2, 2):
        raise ValueError("crossed-string elements must each contain two 2-D endpoints")
    distances = (
        np.linalg.norm(first[0] - second[1])
        + np.linalg.norm(first[1] - second[0])
        - np.linalg.norm(first[0] - second[0])
        - np.linalg.norm(first[1] - second[1]))
    # Endpoint orientation changes only the sign of the named crossed/uncrossed strings.
    return 0.5 * abs(float(distances))


def _subdivide(segment, count):
    parameter = np.linspace(0.0, 1.0, count + 1)
    points = ((1.0 - parameter[:, None]) * segment[0]
              + parameter[:, None] * segment[1])
    return np.stack((points[:-1], points[1:]), axis=1)


def _facing(first, second, first_normal, second_normal, tolerance):
    direction = np.mean(second, axis=0) - np.mean(first, axis=0)
    distance = float(np.linalg.norm(direction))
    if distance <= tolerance:
        return False
    unit = direction / distance
    return bool(
        np.dot(first_normal, unit) > tolerance
        and np.dot(second_normal, -unit) > tolerance)


def _classify_exchange_cell(first, second, segments, excluded, tolerance, *, candidates=None,
                            segment_lower=None, segment_upper=None):
    possible = _candidate_blockers(
        first, second, segments, excluded, tolerance, candidates=candidates,
        segment_lower=segment_lower, segment_upper=segment_upper)
    if not possible:
        return "visible", 1.0, possible

    # Use interior Gauss abscissae. Endpoint strings have zero measure in the exchange
    # integral and, in a closed polygon, lie exactly along adjacent opaque walls.
    offset = np.sqrt(3.0 / 5.0) / 2.0
    parameters = (0.5 - offset, 0.5, 0.5 + offset)
    blockers = []
    for first_parameter in parameters:
        first_point = ((1.0 - first_parameter) * first[0]
                       + first_parameter * first[1])
        for second_parameter in parameters:
            second_point = ((1.0 - second_parameter) * second[0]
                            + second_parameter * second[1])
            blockers.append(_connection_blocker(
                first_point, second_point, segments, excluded, tolerance, possible))
    if all(blocker is None for blocker in blockers):
        return "visible", 1.0, possible
    if blockers[0] is not None and all(blocker == blockers[0] for blocker in blockers):
        # A single straight opaque segment intercepts all extreme and interior connector
        # samples. On this convex connector cell it separates the two elements completely.
        return "blocked", 0.0, possible
    return (
        "mixed", float(sum(blocker is None for blocker in blockers) / len(blockers)),
        possible)


def _adaptive_pair_exchange(segments, normals, first_index, second_index, *,
                            relative_tolerance, absolute_tolerance,
                            minimum_refinement_level, maximum_refinement_level,
                            geometry_tolerance, segment_lower, segment_upper):
    first = segments[first_index]
    second = segments[second_index]
    if not _facing(
            first, second, normals[first_index], normals[second_index],
            geometry_tolerance):
        return 0.0, 0.0, 0
    unobstructed = unobstructed_crossed_string_exchange_2d(first, second)
    if unobstructed <= absolute_tolerance:
        return unobstructed, 0.0, 0

    excluded = {first_index, second_index}
    status, visible_fraction, possible = _classify_exchange_cell(
        first, second, segments, excluded, geometry_tolerance,
        segment_lower=segment_lower, segment_upper=segment_upper)
    if status == "visible":
        return unobstructed, 0.0, 0
    if status == "blocked" and minimum_refinement_level == 0:
        return 0.0, 0.0, 0

    # Mixed cells carry a worst-case error equal to their entire unobstructed exchange.
    # Refine the largest such cell first. Resolved children contribute exact crossed-string
    # exchange, so the remaining heap sum is a conservative unresolved-shadow budget.
    serial = 0
    estimate = visible_fraction * unobstructed
    unresolved = unobstructed
    heap = [(
        -unobstructed, serial, first, second, 0, visible_fraction, status, possible)]
    maximum_used = 0
    threshold = absolute_tolerance + relative_tolerance * unobstructed
    while heap and (unresolved > threshold or maximum_used < minimum_refinement_level):
        (negative_exchange, _, first_cell, second_cell, depth, fraction, cell_status,
         cell_candidates) = heapq.heappop(heap)
        cell_exchange = -negative_exchange
        estimate -= fraction * cell_exchange
        unresolved -= cell_exchange
        if depth >= maximum_refinement_level:
            raise RuntimeError(
                "deterministic line exchange did not close its shadow-boundary error budget; "
                f"pair=({first_index},{second_index}), unresolved={unresolved + cell_exchange:.6g}")
        first_children = _subdivide(first_cell, 2)
        second_children = _subdivide(second_cell, 2)
        child_depth = depth + 1
        maximum_used = max(maximum_used, child_depth)
        for first_child in first_children:
            for second_child in second_children:
                if not _facing(
                        first_child, second_child,
                        normals[first_index], normals[second_index], geometry_tolerance):
                    continue
                child_exchange = unobstructed_crossed_string_exchange_2d(
                    first_child, second_child)
                if child_exchange <= 0.0:
                    continue
                child_status, child_fraction, child_candidates = _classify_exchange_cell(
                    first_child, second_child, segments, excluded, geometry_tolerance,
                    candidates=cell_candidates, segment_lower=segment_lower,
                    segment_upper=segment_upper)
                if child_status == "visible":
                    estimate += child_exchange
                elif child_status == "blocked" and child_depth >= minimum_refinement_level:
                    continue
                else:
                    estimate += child_fraction * child_exchange
                    unresolved += child_exchange
                    serial += 1
                    heapq.heappush(heap, (
                        -child_exchange, serial, first_child, second_child,
                        child_depth, child_fraction, child_status, child_candidates))
    return float(estimate), float(max(unresolved, 0.0)), int(maximum_used)


@dataclass(frozen=True)
class DeterministicLineExchange2D:
    """Reciprocal open-boundary exchange operator per unit extrusion depth."""

    segments: np.ndarray
    gas_normals: np.ndarray
    exchange_length: np.ndarray
    transfer_fraction: np.ndarray
    escape_fraction: np.ndarray
    refinement_level: np.ndarray
    estimated_absolute_error: np.ndarray
    relative_tolerance: float
    absolute_tolerance: float
    maximum_refinement_level: int
    fingerprint: str

    def __post_init__(self):
        segments = _readonly(self.segments)
        normals = _readonly(self.gas_normals)
        exchange = _readonly(self.exchange_length)
        transfer = _readonly(self.transfer_fraction)
        escape = _readonly(self.escape_fraction)
        level = _readonly(self.refinement_level, int)
        error = _readonly(self.estimated_absolute_error)
        count = len(segments)
        if (segments.shape != (count, 2, 2) or normals.shape != (count, 2)
                or exchange.shape != (count, count) or transfer.shape != (count, count)
                or escape.shape != (count,) or level.shape != (count, count)
                or error.shape != (count, count)
                or np.any(~np.isfinite(exchange)) or np.any(exchange < 0.0)
                or np.any(~np.isfinite(transfer)) or np.any(transfer < 0.0)
                or np.any(~np.isfinite(escape)) or np.any(escape < 0.0)
                or not np.allclose(exchange, exchange.T, rtol=0.0, atol=5e-13)
                or not np.allclose(
                    np.sum(transfer, axis=1) + escape, 1.0,
                    rtol=0.0, atol=5e-13)):
            raise ValueError("invalid deterministic line-exchange result")
        for name, value in (
                ("segments", segments), ("gas_normals", normals),
                ("exchange_length", exchange), ("transfer_fraction", transfer),
                ("escape_fraction", escape), ("refinement_level", level),
                ("estimated_absolute_error", error)):
            object.__setattr__(self, name, value)

    @property
    def segment_length(self):
        return np.linalg.norm(self.segments[:, 1] - self.segments[:, 0], axis=1)

    @property
    def maximum_estimated_absolute_error(self):
        return float(np.max(self.estimated_absolute_error, initial=0.0))

    def as_diffuse_form_factors_3d(self):
        """Adapt the deterministic matrix to the common radiosity consumer contract.

        ``rays_per_face=1`` is a compatibility cardinality, not a sampling claim.  The
        construction method and convergence receipt remain authoritative on this object.
        """
        source, target = np.nonzero(self.transfer_fraction > 0.0)
        return DiffuseFormFactors3D(
            len(self.segments), source, target,
            self.transfer_fraction[source, target], self.escape_fraction,
            rays_per_face=1)


def build_deterministic_line_exchange_2d(
        segments, gas_normals, *, relative_tolerance=1.0e-5,
        absolute_tolerance=1.0e-12, minimum_refinement_level=2,
        maximum_refinement_level=18, geometry_tolerance=1.0e-12):
    """Build a reciprocal crossed-string operator with deterministic blocking refinement."""
    segment = np.asarray(segments, dtype=float)
    normal = np.asarray(gas_normals, dtype=float)
    count = len(segment)
    length = np.linalg.norm(segment[:, 1] - segment[:, 0], axis=1)
    normal_length = np.linalg.norm(normal, axis=1)
    relative_tolerance = float(relative_tolerance)
    absolute_tolerance = float(absolute_tolerance)
    minimum = int(minimum_refinement_level)
    maximum = int(maximum_refinement_level)
    geometry_tolerance = float(geometry_tolerance)
    if (segment.shape != (count, 2, 2) or count == 0 or normal.shape != (count, 2)
            or np.any(~np.isfinite(segment)) or np.any(~np.isfinite(normal))
            or np.any(length <= 0.0) or np.any(normal_length <= 0.0)
            or relative_tolerance <= 0.0 or absolute_tolerance < 0.0
            or minimum < 0 or maximum < minimum or geometry_tolerance <= 0.0):
        raise ValueError("invalid deterministic line-exchange inputs")
    normal = normal / normal_length[:, None]
    segment_lower = np.min(segment, axis=1)
    segment_upper = np.max(segment, axis=1)

    exchange = np.zeros((count, count), dtype=float)
    level_used = np.zeros((count, count), dtype=int)
    error = np.zeros((count, count), dtype=float)
    for first in range(count):
        for second in range(first + 1, count):
            current, uncertainty, level = _adaptive_pair_exchange(
                segment, normal, first, second,
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                minimum_refinement_level=minimum,
                maximum_refinement_level=maximum,
                geometry_tolerance=geometry_tolerance,
                segment_lower=segment_lower, segment_upper=segment_upper)
            exchange[first, second] = exchange[second, first] = current
            level_used[first, second] = level_used[second, first] = level
            error[first, second] = error[second, first] = uncertainty

    transfer = exchange / length[:, None]
    outgoing = np.sum(transfer, axis=1)
    closure_tolerance = np.maximum(
        absolute_tolerance / length,
        16.0 * relative_tolerance)
    if np.any(outgoing > 1.0 + closure_tolerance):
        worst = int(np.argmax(outgoing - 1.0))
        raise RuntimeError(
            f"deterministic line exchange exceeds unit row closure on segment {worst}")
    escape = np.maximum(1.0 - outgoing, 0.0)
    # Roundoff-scale closure excess is assigned to escape, never redistributed among receivers.
    row_total = outgoing + escape
    escape += 1.0 - row_total

    digest = sha256()
    digest.update(b"petch.deterministic-line-exchange-2d.v1\0")
    for name, value, dtype in (
            ("segments", segment, "<f8"), ("gas_normals", normal, "<f8"),
            ("exchange_length", exchange, "<f8"),
            ("refinement_level", level_used, "<i8")):
        array = np.ascontiguousarray(value, dtype=dtype)
        digest.update(name.encode("ascii") + b"\0")
        digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
        digest.update(array.tobytes())
    digest.update(np.asarray([
        relative_tolerance, absolute_tolerance, minimum, maximum,
        geometry_tolerance], dtype="<f8").tobytes())
    return DeterministicLineExchange2D(
        segments=segment, gas_normals=normal, exchange_length=exchange,
        transfer_fraction=transfer, escape_fraction=escape,
        refinement_level=level_used, estimated_absolute_error=error,
        relative_tolerance=relative_tolerance, absolute_tolerance=absolute_tolerance,
        maximum_refinement_level=maximum, fingerprint=digest.hexdigest())
