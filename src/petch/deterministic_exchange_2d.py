"""Deterministic diffuse exchange for surfaces extruded normal to a 2-D section.

For an infinitely extruded diffuse surface, Hottel's crossed-string relation gives the
unobstructed exchange length between two line elements exactly.  Obstructed pairs are
resolved by the ``analytic_occlusion`` method: for each source point, every visibility
transition along the target is a projective image of a blocker endpoint (or a facing-clip
root), each elementary target interval is classified once by the exact connector predicate,
and visible intervals contribute the closed-form two-dimensional point-to-segment factor.
Only the outer source integral is numerical -- certified adaptive Simpson panels split at
the analytic visibility events -- so shadow boundaries are never located by bisection and
grazing shadows cost the same as transversal ones.  This module is an authority/cross-check
for line-trench mean profiles; it is not a replacement for the genuinely three-dimensional
hard-visibility operator.

The older ``adaptive_refinement`` method (two-dimensional cell bisection with Gauss-sampled
visibility classification) remains as the per-pair fallback and cross-check.  Its receipt
covers unresolved mixed cells but NOT blocked slivers that evade all Gauss samples inside a
cell it classifies visible, so on grazing-shadow geometry it overcounts by up to a few
tenths of a percent (see the module tests, which document the gap against dense references);
the analytic method has no such term.

The primary stored quantity is the symmetric exchange length ``H_ij``.  Form factors are derived
as ``F_ij = H_ij / L_i``, so reciprocity is true by construction rather than repaired after the
fact.  Open-boundary escape is the unassigned part of each source row.

The geometry kernels operate on scalar floats (with a vectorized candidate prefilter for the
full-surface scan).  A production trench section visits hundreds of thousands of evaluation
points, so per-point work must not pay array-dispatch overhead on two-element vectors.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from math import sqrt
import heapq
import multiprocessing
import os
import sys

import numpy as np

from .neutral_radiosity_3d import DiffuseFormFactors3D

# Interior Gauss abscissae for connector visibility sampling.  Endpoint strings have zero
# measure in the exchange integral and, in a closed polygon, lie exactly along adjacent
# opaque walls.
_GAUSS_OFFSET = sqrt(3.0 / 5.0) / 2.0
_GAUSS_PARAMETERS = (0.5 - _GAUSS_OFFSET, 0.5, 0.5 + _GAUSS_OFFSET)


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


def _endpoint_hull(fax, fay, fbx, fby, sax, say, sbx, sby):
    """Convex hull (counterclockwise monotone chain) of the four cell endpoints."""
    return _hull_points(((fax, fay), (fbx, fby), (sax, say), (sbx, sby)))


def _hull_points(points):
    unique = sorted(set(points))
    if len(unique) <= 2:
        return unique

    def half(sequence):
        output = []
        for point in sequence:
            while len(output) >= 2:
                bx, by = output[-2]
                cx, cy = output[-1]
                if (cx - bx) * (point[1] - cy) - (cy - by) * (point[0] - cx) <= 0.0:
                    output.pop()
                else:
                    break
            output.append(point)
        return output

    lower = half(unique)
    upper = half(unique[::-1])
    return lower[:-1] + upper[:-1]


def _proper_connector_hit(px, py, qx, qy, bx0, by0, bx1, by1, tolerance):
    """Whether a blocker crosses the open interior of the connector (p, q)."""
    fdx = qx - px
    fdy = qy - py
    sdx = bx1 - bx0
    sdy = by1 - by0
    denominator = fdx * sdy - fdy * sdx
    scale = max(sqrt(fdx * fdx + fdy * fdy), sqrt(sdx * sdx + sdy * sdy), 1.0)
    if abs(denominator) <= tolerance * scale:
        # Collinear overlap blocks only when it occupies a nonzero open interval of the
        # connector. Merely sharing a surface endpoint is not an obstruction.
        if abs((bx0 - px) * fdy - (by0 - py) * fdx) > tolerance * scale:
            return False
        length2 = fdx * fdx + fdy * fdy
        if length2 <= tolerance * tolerance:
            return False
        start = ((bx0 - px) * fdx + (by0 - py) * fdy) / length2
        stop = ((bx1 - px) * fdx + (by1 - py) * fdy) / length2
        low, high = (start, stop) if start <= stop else (stop, start)
        return max(low, tolerance) < min(high, 1.0 - tolerance)
    dx0 = bx0 - px
    dy0 = by0 - py
    first_parameter = (dx0 * sdy - dy0 * sdx) / denominator
    if not tolerance < first_parameter < 1.0 - tolerance:
        return False
    second_parameter = (dx0 * fdy - dy0 * fdx) / denominator
    return -tolerance <= second_parameter <= 1.0 + tolerance


def _transverse_edge_hit(bx0, by0, bx1, by1, ex0, ey0, ex1, ey1, tolerance):
    fdx = bx1 - bx0
    fdy = by1 - by0
    sdx = ex1 - ex0
    sdy = ey1 - ey0
    denominator = fdx * sdy - fdy * sdx
    scale = max(sqrt(fdx * fdx + fdy * fdy), sqrt(sdx * sdx + sdy * sdy), 1.0)
    if abs(denominator) <= tolerance * scale:
        return False
    dx0 = ex0 - bx0
    dy0 = ey0 - by0
    first_parameter = (dx0 * sdy - dy0 * sdx) / denominator
    if not tolerance < first_parameter < 1.0 - tolerance:
        return False
    second_parameter = (dx0 * fdy - dy0 * fdx) / denominator
    return tolerance < second_parameter < 1.0 - tolerance


def _blocker_intersects_hull(hull, bx0, by0, bx1, by1, tolerance):
    """Strict-interior or transverse-boundary test of one blocker against a hull."""
    if len(hull) >= 3:
        mx = 0.5 * (bx0 + bx1)
        my = 0.5 * (by0 + by1)
        for px, py in ((bx0, by0), (bx1, by1), (mx, my)):
            positive = True
            negative = True
            for index in range(len(hull)):
                ex0, ey0 = hull[index]
                ex1, ey1 = hull[(index + 1) % len(hull)]
                signed = (ex1 - ex0) * (py - ey0) - (ey1 - ey0) * (px - ex0)
                if signed <= tolerance:
                    positive = False
                if signed >= -tolerance:
                    negative = False
                if not positive and not negative:
                    break
            if positive or negative:
                return True
    for index in range(len(hull)):
        ex0, ey0 = hull[index]
        ex1, ey1 = hull[(index + 1) % len(hull)]
        if _transverse_edge_hit(bx0, by0, bx1, by1, ex0, ey0, ex1, ey1, tolerance):
            return True
    return False


def _candidate_blockers(cell_first, cell_second, geometry, tolerance, candidates):
    """Return the blockers that can intersect the convex connector family.

    Descendant exchange cells inherit their parent's candidate set.  A blocker excluded by a
    parent hull cannot re-enter a child hull, which makes adaptive shadow refinement scale with
    the local occluders instead of repeatedly scanning the entire surface.  The full-surface
    scan (``candidates is None``) is vectorized; inherited sets are small and stay scalar.
    """
    fax, fay, fbx, fby = cell_first
    sax, say, sbx, sby = cell_second
    hull = _endpoint_hull(fax, fay, fbx, fby, sax, say, sbx, sby)
    if len(hull) == 0:
        return ()
    return _segment_hull_candidates(hull, geometry, tolerance, candidates)


def _segment_hull_candidates(hull, geometry, tolerance, candidates):
    lower_x = min(point[0] for point in hull) - tolerance
    lower_y = min(point[1] for point in hull) - tolerance
    upper_x = max(point[0] for point in hull) + tolerance
    upper_y = max(point[1] for point in hull) + tolerance
    segments_f = geometry.segments_f
    bounds_f = geometry.bounds_f
    if candidates is None:
        overlap = (
            (geometry.segment_upper[:, 0] >= lower_x)
            & (geometry.segment_upper[:, 1] >= lower_y)
            & (geometry.segment_lower[:, 0] <= upper_x)
            & (geometry.segment_lower[:, 1] <= upper_y))
        overlap[geometry.excluded_first] = False
        overlap[geometry.excluded_second] = False
        pool = np.flatnonzero(overlap)
    else:
        pool = candidates
    result = []
    for index in pool:
        index = int(index)
        if candidates is not None:
            if index == geometry.excluded_first or index == geometry.excluded_second:
                continue
            blx, bly, bux, buy = bounds_f[index]
            if bux < lower_x or buy < lower_y or blx > upper_x or bly > upper_y:
                continue
        bx0, by0, bx1, by1 = segments_f[index]
        if _blocker_intersects_hull(hull, bx0, by0, bx1, by1, tolerance):
            result.append(index)
    return tuple(result)


def _classify_exchange_cell(cell_first, cell_second, geometry, tolerance, candidates):
    possible = _candidate_blockers(cell_first, cell_second, geometry, tolerance, candidates)
    if not possible:
        return "visible", 1.0, possible

    fax, fay, fbx, fby = cell_first
    sax, say, sbx, sby = cell_second
    segments_f = geometry.segments_f
    first_blocker = None
    uniform = True
    visible_count = 0
    blocked_count = 0
    for first_parameter in _GAUSS_PARAMETERS:
        px = (1.0 - first_parameter) * fax + first_parameter * fbx
        py = (1.0 - first_parameter) * fay + first_parameter * fby
        for second_parameter in _GAUSS_PARAMETERS:
            qx = (1.0 - second_parameter) * sax + second_parameter * sbx
            qy = (1.0 - second_parameter) * say + second_parameter * sby
            hit = None
            for index in possible:
                bx0, by0, bx1, by1 = segments_f[index]
                if _proper_connector_hit(px, py, qx, qy, bx0, by0, bx1, by1, tolerance):
                    hit = index
                    break
            if hit is None:
                visible_count += 1
                uniform = False
            else:
                blocked_count += 1
                if first_blocker is None:
                    if visible_count:
                        uniform = False
                    first_blocker = hit
                elif hit != first_blocker:
                    uniform = False
    if blocked_count == 0:
        return "visible", 1.0, possible
    if uniform and first_blocker is not None and visible_count == 0:
        # A single straight opaque segment intercepts all extreme and interior connector
        # samples. On this convex connector cell it separates the two elements completely.
        return "blocked", 0.0, possible
    return "mixed", float(visible_count / (visible_count + blocked_count)), possible


def _facing_cells(cell_first, cell_second, first_normal, second_normal, tolerance):
    fax, fay, fbx, fby = cell_first
    sax, say, sbx, sby = cell_second
    dx = 0.5 * (sax + sbx) - 0.5 * (fax + fbx)
    dy = 0.5 * (say + sby) - 0.5 * (fay + fby)
    distance = sqrt(dx * dx + dy * dy)
    if distance <= tolerance:
        return False
    ux = dx / distance
    uy = dy / distance
    return (first_normal[0] * ux + first_normal[1] * uy > tolerance
            and second_normal[0] * (-ux) + second_normal[1] * (-uy) > tolerance)


def _crossed_string_cells(cell_first, cell_second):
    fax, fay, fbx, fby = cell_first
    sax, say, sbx, sby = cell_second
    dx = fax - sbx
    dy = fay - sby
    total = sqrt(dx * dx + dy * dy)
    dx = fbx - sax
    dy = fby - say
    total += sqrt(dx * dx + dy * dy)
    dx = fax - sax
    dy = fay - say
    total -= sqrt(dx * dx + dy * dy)
    dx = fbx - sbx
    dy = fby - sby
    total -= sqrt(dx * dx + dy * dy)
    # Endpoint orientation changes only the sign of the named crossed/uncrossed strings.
    return 0.5 * abs(total)


def _split_cell(cell):
    ax, ay, bx, by = cell
    mx = (1.0 - 0.5) * ax + 0.5 * bx
    my = (1.0 - 0.5) * ay + 0.5 * by
    return ((ax, ay, mx, my), (mx, my, bx, by))


class _PairGeometry:
    """Shared per-build scalar geometry plus the vectorized prefilter arrays."""

    __slots__ = ("segments_f", "bounds_f", "segment_lower", "segment_upper",
                 "excluded_first", "excluded_second")

    def __init__(self, segments, segment_lower, segment_upper):
        self.segments_f = [
            (float(item[0][0]), float(item[0][1]), float(item[1][0]), float(item[1][1]))
            for item in segments]
        self.bounds_f = [
            (float(low[0]), float(low[1]), float(high[0]), float(high[1]))
            for low, high in zip(segment_lower, segment_upper)]
        self.segment_lower = segment_lower
        self.segment_upper = segment_upper
        self.excluded_first = -1
        self.excluded_second = -1


def _adaptive_pair_exchange(geometry, normals_f, first_index, second_index, *,
                            relative_tolerance, absolute_tolerance,
                            minimum_refinement_level, maximum_refinement_level,
                            geometry_tolerance):
    first = geometry.segments_f[first_index]
    second = geometry.segments_f[second_index]
    first_normal = normals_f[first_index]
    second_normal = normals_f[second_index]
    if not _facing_cells(first, second, first_normal, second_normal, geometry_tolerance):
        return 0.0, 0.0, 0
    unobstructed = _crossed_string_cells(first, second)
    if unobstructed <= absolute_tolerance:
        return unobstructed, 0.0, 0

    geometry.excluded_first = first_index
    geometry.excluded_second = second_index
    status, visible_fraction, possible = _classify_exchange_cell(
        first, second, geometry, geometry_tolerance, None)
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
        first_children = _split_cell(first_cell)
        second_children = _split_cell(second_cell)
        child_depth = depth + 1
        maximum_used = max(maximum_used, child_depth)
        for first_child in first_children:
            for second_child in second_children:
                if not _facing_cells(
                        first_child, second_child, first_normal, second_normal,
                        geometry_tolerance):
                    continue
                child_exchange = _crossed_string_cells(first_child, second_child)
                if child_exchange <= 0.0:
                    continue
                child_status, child_fraction, child_candidates = _classify_exchange_cell(
                    first_child, second_child, geometry, geometry_tolerance,
                    cell_candidates)
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


def _point_segment_exchange(px, py, source_normal, target, target_normal, blockers_f,
                            tolerance):
    """Exact per-source-point exchange factor to the visible parts of the target segment.

    Every visibility transition along the target parameter is a projective image of a
    blocker endpoint (or a facing-clip root), so the parameter axis is cut at those points
    and each elementary interval is classified once, at its midpoint, by the authoritative
    connector predicate.  Visible intervals contribute the exact two-dimensional
    point-to-segment factor 0.5 * |sin(theta_hi) - sin(theta_lo)|.
    """
    tax, tay, tbx, tby = target
    dbx = tbx - tax
    dby = tby - tay
    target_length = sqrt(dbx * dbx + dby * dby)
    if abs(dbx * (py - tay) - dby * (px - tax)) <= tolerance * target_length:
        return 0.0
    npx, npy = source_normal
    ntx, nty = target_normal
    cuts = [0.0, 1.0]
    # Facing clips are linear in the target parameter.
    for constant, slope in (
            (npx * (tax - px) + npy * (tay - py), npx * dbx + npy * dby),
            (ntx * (px - tax) + nty * (py - tay), -(ntx * dbx + nty * dby))):
        if slope != 0.0:
            root = -constant / slope
            if 0.0 < root < 1.0:
                cuts.append(root)
    rx = tax - px
    ry = tay - py
    for bx0, by0, bx1, by1 in blockers_f:
        for ex, ey in ((bx0, by0), (bx1, by1)):
            dx = ex - px
            dy = ey - py
            denominator = dx * dby - dy * dbx
            if denominator == 0.0:
                continue
            forward = (rx * dby - ry * dbx) / denominator
            if forward <= 0.0:
                continue
            root = (rx * dy - ry * dx) / denominator
            if 0.0 < root < 1.0:
                cuts.append(root)
    cuts.sort()
    factor = 0.0
    for index in range(len(cuts) - 1):
        low = cuts[index]
        high = cuts[index + 1]
        if high - low <= 1.0e-14:
            continue
        middle = 0.5 * (low + high)
        qx = tax + middle * dbx
        qy = tay + middle * dby
        dx = qx - px
        dy = qy - py
        if npx * dx + npy * dy <= 0.0:
            continue
        if ntx * (px - qx) + nty * (py - qy) <= 0.0:
            continue
        blocked = False
        for blocker in blockers_f:
            if _proper_connector_hit(
                    px, py, qx, qy, blocker[0], blocker[1], blocker[2], blocker[3],
                    tolerance):
                blocked = True
                break
        if blocked:
            continue
        lx = tax + low * dbx - px
        ly = tay + low * dby - py
        hx = tax + high * dbx - px
        hy = tay + high * dby - py
        sine_low = (npx * ly - npy * lx) / sqrt(lx * lx + ly * ly)
        sine_high = (npx * hy - npy * hx) / sqrt(hx * hx + hy * hy)
        factor += 0.5 * abs(sine_high - sine_low)
    return factor


def _panel_events(first, second, first_normal, second_normal, blockers_f):
    """Source parameters where the target visibility-interval structure can change.

    Interval births and deaths happen exactly when the source point becomes collinear with
    a blocker endpoint and a target endpoint, or when a facing clip crosses a target
    endpoint; each such condition is linear in the source parameter.
    """
    fax, fay, fbx, fby = first
    dax = fbx - fax
    day = fby - fay
    sax, say, sbx, sby = second
    npx, npy = first_normal
    ntx, nty = second_normal
    events = []

    def add_root(constant, slope):
        if slope != 0.0:
            root = -constant / slope
            if 0.0 < root < 1.0:
                events.append(root)

    for qx, qy in ((sax, say), (sbx, sby)):
        # cross(e - p(s), q - p(s)) = 0 is linear in s for every blocker endpoint e.
        for blocker in blockers_f:
            for ex, ey in ((blocker[0], blocker[1]), (blocker[2], blocker[3])):
                constant = (ex - fax) * (qy - fay) - (ey - fay) * (qx - fax)
                slope = ((ey - qy) * dax - (ex - qx) * day)
                add_root(constant, slope)
        add_root(npx * (qx - fax) + npy * (qy - fay), -(npx * dax + npy * day))
        add_root(ntx * (fax - qx) + nty * (fay - qy), ntx * dax + nty * day)
    return sorted(set(events))


def _analytic_pair_exchange(geometry, normals_f, first_index, second_index, *,
                            relative_tolerance, absolute_tolerance,
                            minimum_refinement_level, maximum_refinement_level,
                            geometry_tolerance):
    """Deterministic obstructed exchange: exact inner occlusion, adaptive outer quadrature.

    The inner integral over the target is evaluated exactly per source point (projective
    shadow intervals classified by the connector predicate).  Only the outer integral over
    the source runs adaptive Simpson panels split at the analytic visibility events, so
    the shadow boundary never has to be located by bisection.  Exhausting the quadrature
    depth falls back to the conservative adaptive shadow refinement.
    """
    first = geometry.segments_f[first_index]
    second = geometry.segments_f[second_index]
    first_normal = normals_f[first_index]
    second_normal = normals_f[second_index]
    if not _facing_cells(first, second, first_normal, second_normal, geometry_tolerance):
        return 0.0, 0.0, 0
    unobstructed = _crossed_string_cells(first, second)
    if unobstructed <= absolute_tolerance:
        return unobstructed, 0.0, 0
    geometry.excluded_first = first_index
    geometry.excluded_second = second_index
    candidates = _candidate_blockers(first, second, geometry, geometry_tolerance, None)
    if not candidates:
        return unobstructed, 0.0, 0
    blockers_f = [geometry.segments_f[index] for index in candidates]

    fax, fay, fbx, fby = first
    dax = fbx - fax
    day = fby - fay
    source_length = sqrt(dax * dax + day * day)

    def integrand(parameter):
        return _point_segment_exchange(
            fax + parameter * dax, fay + parameter * day, first_normal, second,
            second_normal, blockers_f, geometry_tolerance)

    threshold = absolute_tolerance + relative_tolerance * unobstructed
    panels = [0.0] + _panel_events(
        first, second, first_normal, second_normal, blockers_f) + [1.0]
    panel_budget = threshold / (source_length * max(len(panels) - 1, 1))
    total = 0.0
    error_estimate = 0.0
    deepest = 0
    for index in range(len(panels) - 1):
        low = panels[index]
        high = panels[index + 1]
        if high - low <= 1.0e-14:
            continue
        stack = [(low, high, integrand(low), integrand(0.5 * (low + high)),
                  integrand(high), 0)]
        while stack:
            left, right, f_left, f_middle, f_right, depth = stack.pop()
            width = right - left
            coarse = width / 6.0 * (f_left + 4.0 * f_middle + f_right)
            middle = 0.5 * (left + right)
            f_lq = integrand(0.5 * (left + middle))
            f_rq = integrand(0.5 * (middle + right))
            fine = width / 12.0 * (
                f_left + 4.0 * f_lq + 2.0 * f_middle + 4.0 * f_rq + f_right)
            difference = abs(fine - coarse)
            deepest = max(deepest, depth + 1)
            if (difference <= panel_budget * width / max(high - low, 1.0e-300)
                    and depth >= minimum_refinement_level):
                total += fine
                error_estimate += difference
                continue
            if depth >= maximum_refinement_level:
                # The outer quadrature refused to certify; only the shadow-refinement
                # bound is safe for this pair.
                return _adaptive_pair_exchange(
                    geometry, normals_f, first_index, second_index,
                    relative_tolerance=relative_tolerance,
                    absolute_tolerance=absolute_tolerance,
                    minimum_refinement_level=minimum_refinement_level,
                    maximum_refinement_level=maximum_refinement_level,
                    geometry_tolerance=geometry_tolerance)
            stack.append((left, middle, f_left, f_lq, f_middle, depth + 1))
            stack.append((middle, right, f_middle, f_rq, f_right, depth + 1))
    exchange = min(max(total * source_length, 0.0), unobstructed)
    return float(exchange), float(error_estimate * source_length), int(deepest)


_PAIR_METHODS = ("analytic_occlusion", "adaptive_refinement")


def _pair_kernel(method):
    return _analytic_pair_exchange if method == "analytic_occlusion" else _adaptive_pair_exchange


def _pair_rows_task(payload):
    """Compute the adaptive exchange for a batch of matrix rows in a worker process.

    Each pair is computed independently by the identical serial kernel, so the assembled
    matrix does not depend on the worker count or on scheduling order.
    """
    (segment, segment_lower, segment_upper, normals_f, rows, count, method,
     relative_tolerance, absolute_tolerance, minimum, maximum,
     geometry_tolerance) = payload
    geometry = _PairGeometry(segment, segment_lower, segment_upper)
    kernel = _pair_kernel(method)
    results = []
    for first in rows:
        for second in range(first + 1, count):
            results.append((first, second) + kernel(
                geometry, normals_f, first, second,
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                minimum_refinement_level=minimum,
                maximum_refinement_level=maximum,
                geometry_tolerance=geometry_tolerance))
    return results


_PAIR_POOL = {}


def _pair_pool(workers):
    pool = _PAIR_POOL.get(workers)
    if pool is None:
        context = multiprocessing.get_context(
            "fork" if sys.platform.startswith("linux") else "spawn")
        pool = context.Pool(processes=workers)
        _PAIR_POOL[workers] = pool
    return pool


def _configured_worker_count():
    """Opt-in process parallelism for the pair loop; performance-only, output-identical."""
    raw = os.environ.get("PETCH_DETERMINISTIC_EXCHANGE_WORKERS", "").strip()
    if not raw:
        return 1
    workers = int(raw)
    if workers < 1:
        raise ValueError("PETCH_DETERMINISTIC_EXCHANGE_WORKERS must be a positive integer")
    if multiprocessing.current_process().daemon:
        return 1
    return workers


def unobstructed_crossed_string_exchange_2d(first, second):
    """Return the exact symmetric exchange length for two unobstructed line elements."""
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if first.shape != (2, 2) or second.shape != (2, 2):
        raise ValueError("crossed-string elements must each contain two 2-D endpoints")
    return _crossed_string_cells(
        (float(first[0][0]), float(first[0][1]), float(first[1][0]), float(first[1][1])),
        (float(second[0][0]), float(second[0][1]), float(second[1][0]), float(second[1][1])))


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
    method: str
    fingerprint: str

    def __post_init__(self):
        if self.method not in _PAIR_METHODS:
            raise ValueError("unknown deterministic line-exchange method")
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
        segments, gas_normals, *, method="analytic_occlusion", relative_tolerance=1.0e-5,
        absolute_tolerance=1.0e-12, minimum_refinement_level=2,
        maximum_refinement_level=18, geometry_tolerance=1.0e-12):
    """Build a reciprocal crossed-string operator with deterministic blocking resolution.

    ``analytic_occlusion`` resolves blocking exactly in the inner (target) integral via
    projective shadow intervals and runs only a one-dimensional certified quadrature over
    the source, falling back per pair to the conservative adaptive shadow refinement if
    that quadrature cannot certify; ``adaptive_refinement`` uses the two-dimensional
    shadow refinement for every pair and serves as the independent cross-check.  The
    tolerances govern both the quadrature certification and the refinement.
    """
    if method not in _PAIR_METHODS:
        raise ValueError("unknown deterministic line-exchange method")
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
    geometry = _PairGeometry(segment, segment_lower, segment_upper)
    normals_f = [(float(item[0]), float(item[1])) for item in normal]

    exchange = np.zeros((count, count), dtype=float)
    level_used = np.zeros((count, count), dtype=int)
    error = np.zeros((count, count), dtype=float)
    workers = _configured_worker_count()
    if workers > 1 and count >= 64:
        # Interleaved row batches balance the triangular pair workload across processes.
        batches = [range(start, count, 4 * workers) for start in range(4 * workers)]
        payloads = [(
            segment, segment_lower, segment_upper, normals_f, tuple(rows), count, method,
            relative_tolerance, absolute_tolerance, minimum, maximum,
            geometry_tolerance) for rows in batches if len(rows)]
        for results in _pair_pool(workers).map(_pair_rows_task, payloads):
            for first, second, current, uncertainty, level in results:
                exchange[first, second] = exchange[second, first] = current
                level_used[first, second] = level_used[second, first] = level
                error[first, second] = error[second, first] = uncertainty
    else:
        kernel = _pair_kernel(method)
        for first in range(count):
            for second in range(first + 1, count):
                current, uncertainty, level = kernel(
                    geometry, normals_f, first, second,
                    relative_tolerance=relative_tolerance,
                    absolute_tolerance=absolute_tolerance,
                    minimum_refinement_level=minimum,
                    maximum_refinement_level=maximum,
                    geometry_tolerance=geometry_tolerance)
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
    digest.update(b"petch.deterministic-line-exchange-2d.v2\0")
    digest.update(method.encode("ascii") + b"\0")
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
        maximum_refinement_level=maximum, method=method,
        fingerprint=digest.hexdigest())
