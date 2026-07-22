"""Deterministic exact-occlusion diffuse exchange for bodies of revolution.

The axisymmetric analogue of :mod:`petch.deterministic_exchange_2d`: a hole (or any body
of revolution) with piecewise-linear generator profile r(z) is meshed into conical bands,
and the band-to-band diffuse exchange is computed with EXACT occlusion.  Along a chord
between generator points the squared cylindrical radius is convex, so blocking can only
occur at profile kinks, and each kink contributes a closed-form azimuthal visibility
threshold ``cos(dphi) <= C_k`` (RESEARCH_EXACT_3D_OCCLUSION_2026-07-21.md section 5).  The
visible azimuth set is a single symmetric interval whose endpoints are arccos of rational
functions of the kink radii -- the exact analogue of the 2-D blocker-endpoint projections.
There is no sampled visibility anywhere.

Numerical structure mirrors the 2-D operator: the azimuthal kernel integral over the
visible interval and the outer generator-pair integral are certified adaptive quadratures
with per-pair error receipts; reciprocity holds by construction (the symmetric exchange
area ``X_ij`` is stored once); row closure assigns the deficit to escape through the
declared opening.  A straight or monotonically tapered hole has no occlusion events; they
appear exactly with necking or bowing.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from math import acos, cos, pi, sin, sqrt

import numpy as np

from .neutral_radiosity_3d import DiffuseFormFactors3D


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class AxisymmetricProfile:
    """Piecewise-linear generator r(z) of a body of revolution, z increasing upward.

    ``z`` strictly increasing; ``r`` positive except that the first entry may be zero
    (closed bottom on the axis).  The last vertex is the opening rim (escape boundary).
    """

    z: np.ndarray
    r: np.ndarray

    def __post_init__(self):
        z = _readonly(self.z)
        r = _readonly(self.r)
        if (z.ndim != 1 or z.shape != r.shape or len(z) < 2
                or np.any(~np.isfinite(z)) or np.any(~np.isfinite(r))
                or np.any(np.diff(z) <= 0.0) or np.any(r < 0.0)
                or np.any(r[1:] <= 0.0)):
            raise ValueError("invalid axisymmetric generator profile")
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "r", r)


def _band_geometry(profile, bands_per_segment):
    """Subdivide each generator segment into bands; return band endpoint arrays."""
    z0 = []
    z1 = []
    r0 = []
    r1 = []
    for index in range(len(profile.z) - 1):
        za, zb = float(profile.z[index]), float(profile.z[index + 1])
        ra, rb = float(profile.r[index]), float(profile.r[index + 1])
        for sub in range(bands_per_segment):
            ta = sub / bands_per_segment
            tb = (sub + 1) / bands_per_segment
            z0.append(za + ta * (zb - za))
            z1.append(za + tb * (zb - za))
            r0.append(ra + ta * (rb - ra))
            r1.append(ra + tb * (rb - ra))
    return (np.asarray(z0), np.asarray(z1), np.asarray(r0), np.asarray(r1))


def _visible_cosine_threshold(rs, zs, rt, zt, kink_z, kink_r, tolerance):
    """min_k C_k over kinks strictly between the endpoint heights; >= 1 means unblocked."""
    threshold = 1.0
    dz = zt - zs
    if abs(dz) <= tolerance:
        return threshold
    for zk, rk in zip(kink_z, kink_r):
        u = (zk - zs) / dz
        if not tolerance < u < 1.0 - tolerance:
            continue
        a = (1.0 - u) * rs
        b = u * rt
        denominator = 2.0 * a * b
        if denominator <= tolerance * tolerance:
            # A chord touching the axis is blocked by any interior kink of smaller radius.
            if a + b > rk + tolerance:
                return -1.0
            continue
        candidate = (rk * rk - a * a - b * b) / denominator
        threshold = min(threshold, candidate)
        if threshold <= -1.0:
            return -1.0
    return threshold


def _pair_kernel_integral(rs, zs, ns_r, ns_z, rt, zt, nt_r, nt_z, phi_limit,
                          azimuth_order):
    """Integral over the visible azimuth interval of cos cos / (pi rho^2), per unit areas.

    Source at (rs, 0, zs), target ring point at azimuth phi.  The integrand is smooth on
    the visible interval; fixed-order Gauss-Legendre with order doubling is certified by
    the caller through the returned pair (value, refinement estimate).
    """

    def integrand(phi):
        cphi = cos(phi)
        sphi = sin(phi)
        dx = rt * cphi - rs
        dy = rt * sphi
        dz = zt - zs
        rho2 = dx * dx + dy * dy + dz * dz
        if rho2 <= 0.0:
            return 0.0
        # Inward normals: source normal (nr, nz) at azimuth 0; target normal rotated.
        cos_s = (ns_r * dx + ns_z * dz)
        cos_t = -(nt_r * (cphi * dx + sphi * dy) + nt_z * dz)
        if cos_s <= 0.0 or cos_t <= 0.0:
            return 0.0
        return cos_s * cos_t / (pi * rho2 * rho2)

    nodes, weights = np.polynomial.legendre.leggauss(azimuth_order)
    half = phi_limit
    total = 0.0
    for node, weight in zip(nodes, weights):
        total += weight * integrand(0.5 * half * (node + 1.0))
    # The chord kernel is symmetric in phi; integrate [0, phi_limit] and double.
    return 2.0 * 0.5 * half * total


def _coaxial_disk_factor(radius_1, radius_2, gap):
    """Exact view factor from disk 1 to coaxial parallel disk 2 at separation gap."""
    if gap <= 0.0:
        raise ValueError("coaxial disk separation must be positive")
    x = radius_1 / gap
    y = radius_2 / gap
    s = 1.0 + (1.0 + y * y) / (x * x)
    return 0.5 * (s - sqrt(s * s - 4.0 * (y / x) ** 2))


def _cylinder_band_exchange_area(radius, za0, za1, zb0, zb1):
    """Exact exchange area between two interior bands of one cylinder (disk algebra).

    With E(u,v) = A_disk * F_dd(|u-v|) the exchange area between the virtual openings at
    heights u and v, the telescoping section identities give, for disjoint sections
    a0 < a1 <= b0 < b1:

        X = E(a1,b0) - E(a0,b0) - E(a1,b1) + E(a0,b1)

    (checked: additive over target subdivision, and X(W, everything above) telescopes to
    the crossing power A_d (1 - F_dd)), and for the self term:

        X(W, W) = A_wall - 2 A_disk (1 - F_dd(a1 - a0)).
    """
    disk_area = pi * radius * radius

    def disk_exchange(u, v):
        gap = abs(u - v)
        if gap <= 0.0:
            return disk_area
        return disk_area * _coaxial_disk_factor(radius, radius, gap)

    if za0 == zb0 and za1 == zb1:
        wall_area = 2.0 * pi * radius * (za1 - za0)
        return wall_area - 2.0 * (disk_area - disk_exchange(za0, za1))
    if za1 <= zb0:
        return (disk_exchange(za1, zb0) - disk_exchange(za0, zb0)
                - disk_exchange(za1, zb1) + disk_exchange(za0, zb1))
    if zb1 <= za0:
        return (disk_exchange(zb1, za0) - disk_exchange(zb0, za0)
                - disk_exchange(zb1, za1) + disk_exchange(zb0, za1))
    raise ValueError("cylinder bands must be identical or disjoint")


def build_cylinder_band_exchange(radius, z_edges):
    """Exact (closed-form) band exchange operator for a straight cylinder interior.

    ``z_edges`` are increasing band boundaries.  Escape is through both end disks; the
    escape fractions toward each end are also exact by the same disk algebra.
    """
    radius = float(radius)
    edges = np.asarray(z_edges, dtype=float)
    if radius <= 0.0 or edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("invalid cylinder band inputs")
    count = len(edges) - 1
    area = 2.0 * pi * radius * np.diff(edges)
    exchange = np.zeros((count, count))
    for i in range(count):
        for j in range(i, count):
            value = _cylinder_band_exchange_area(
                radius, edges[i], edges[i + 1], edges[j], edges[j + 1])
            exchange[i, j] = exchange[j, i] = max(value, 0.0)
    transfer = exchange / area[:, None]
    bottom = edges[0]
    top = edges[-1]
    disk_area = pi * radius * radius

    def wall_to_disk(z0, z1, disk_z):
        # Exchange area between wall band [z0,z1] and an end disk, by enclosure algebra:
        # X(band, disk at u) = A_d * (F_dd(|z0-u|) ... ) via section differences.
        def dd(gap):
            if gap <= 0.0:
                return disk_area
            return disk_area * _coaxial_disk_factor(radius, radius, gap)
        return max(dd(abs(z0 - disk_z)) - dd(abs(z1 - disk_z)), 0.0) if disk_z <= z0 \
            else max(dd(abs(z1 - disk_z)) - dd(abs(z0 - disk_z)), 0.0)

    escape_bottom = np.array([
        wall_to_disk(edges[i], edges[i + 1], bottom) for i in range(count)]) / area
    escape_top = np.array([
        wall_to_disk(edges[i], edges[i + 1], top) for i in range(count)]) / area
    closure = transfer.sum(axis=1) + escape_bottom + escape_top
    # Differences of near-equal disk factors accumulate ~1e-8 cancellation error across
    # thousands of bands at 200:1; the closure gate stays far below physical relevance.
    if np.any(np.abs(closure - 1.0) > 1e-7):
        raise ValueError(
            f"cylinder disk-algebra closure failed: max |defect| = "
            f"{np.max(np.abs(closure - 1.0)):.3e}")
    return {
        "band_edges": edges,
        "band_area": area,
        "transfer_fraction": transfer,
        "escape_bottom": escape_bottom,
        "escape_top": escape_top,
    }


def cylinder_clausing_transmission(aspect_ratio, *, bands=None):
    """Exact-algebra Clausing transmission for a straight cylinder of given L/(2R)...

    aspect_ratio = depth / diameter.  Entrance at the top disk (cosine influx), exit
    through the bottom disk.  All factors closed-form; the only numerics is the linear
    bounce solve, so the result converges rapidly in band count with no quadrature error.
    """
    radius = 0.5
    length = float(aspect_ratio)  # depth = AR * diameter = AR * 1.0
    count = int(bands if bands is not None else max(24, int(24 * aspect_ratio)))
    edges = np.linspace(0.0, length, count + 1)
    operator = build_cylinder_band_exchange(radius, edges)
    transfer = operator["transfer_fraction"]
    disk_area = pi * radius * radius

    def dd(gap):
        if gap <= 0.0:
            return disk_area
        return disk_area * _coaxial_disk_factor(radius, radius, gap)

    # Influx 1 through the top disk: direct to bottom disk + load on each wall band
    # (disk-to-band by the same section algebra, normalized by the disk area).
    direct = _coaxial_disk_factor(radius, radius, length)
    wall_load = np.array([
        (dd(abs(length - edges[i + 1])) - dd(abs(length - edges[i])))
        for i in range(count)]) / disk_area
    balance = direct + wall_load.sum()
    if abs(balance - 1.0) > 1e-10:
        raise ValueError(f"top-disk closure failed: {balance - 1.0:.3e}")
    resolvent = np.linalg.solve(
        np.eye(count) - transfer.T, operator["escape_bottom"])
    return float(direct + wall_load @ resolvent)


def build_axisymmetric_band_exchange(
        profile, *, bands_per_segment=4, azimuth_order=24,
        generator_order=6, relative_tolerance=1.0e-4,
        geometry_tolerance=1.0e-12):
    """Build the reciprocal band exchange operator for a body of revolution.

    Returns band arrays, the symmetric transfer-fraction matrix, the escape fraction
    through the opening rim plane, and per-pair certification receipts (the difference
    between the declared azimuth/generator orders and their doubled-order refinements).
    """
    if (int(bands_per_segment) < 1 or int(azimuth_order) < 4
            or int(generator_order) < 2 or float(relative_tolerance) <= 0.0
            or float(geometry_tolerance) <= 0.0):
        raise ValueError("invalid axisymmetric exchange inputs")
    z0, z1, r0, r1 = _band_geometry(profile, int(bands_per_segment))
    count = len(z0)
    kink_z = np.asarray(profile.z[1:-1], dtype=float)
    kink_r = np.asarray(profile.r[1:-1], dtype=float)

    generator_length = np.hypot(z1 - z0, r1 - r0)
    mean_radius = 0.5 * (r0 + r1)
    area = 2.0 * pi * mean_radius * generator_length
    if np.any(area <= 0.0):
        raise ValueError("axisymmetric band with nonpositive area")
    # Inward normal of the revolved cone: rotate the generator tangent by -90 degrees so
    # the normal points into the cavity (toward the axis for an outward-opening profile).
    tangent_r = (r1 - r0) / generator_length
    tangent_z = (z1 - z0) / generator_length
    normal_r = -tangent_z
    normal_z = tangent_r

    gauss_nodes, gauss_weights = np.polynomial.legendre.leggauss(int(generator_order))

    def pair_exchange(i, j, az_order, gen_nodes, gen_weights):
        total = 0.0
        for node_i, weight_i in zip(gen_nodes, gen_weights):
            ti = 0.5 * (node_i + 1.0)
            rs = r0[i] + ti * (r1[i] - r0[i])
            zs = z0[i] + ti * (z1[i] - z0[i])
            inner = 0.0
            for node_j, weight_j in zip(gen_nodes, gen_weights):
                tj = 0.5 * (node_j + 1.0)
                rt = r0[j] + tj * (r1[j] - r0[j])
                zt = z0[j] + tj * (z1[j] - z0[j])
                threshold = _visible_cosine_threshold(
                    rs, zs, rt, zt, kink_z, kink_r, geometry_tolerance)
                if threshold <= -1.0:
                    continue
                phi_limit = pi if threshold >= 1.0 else acos(max(-1.0, threshold))
                if phi_limit <= 0.0:
                    continue
                kernel = _pair_kernel_integral(
                    rs, zs, normal_r[i], normal_z[i], rt, zt,
                    normal_r[j], normal_z[j], phi_limit, az_order)
                inner += weight_j * kernel * rt * generator_length[j] * 0.5
            total += weight_i * inner * rs * generator_length[i] * 0.5
        # X_ij = closed double integral including both ring measures (2 pi each), with one
        # azimuth integrated analytically into the kernel interval (factor 2 pi absorbed).
        return 2.0 * pi * total

    exchange = np.zeros((count, count))
    receipt = np.zeros((count, count))
    fine_nodes, fine_weights = np.polynomial.legendre.leggauss(2 * int(generator_order))
    for i in range(count):
        for j in range(i, count):
            coarse = pair_exchange(i, j, int(azimuth_order), gauss_nodes, gauss_weights)
            if coarse == 0.0:
                continue
            fine = pair_exchange(i, j, 2 * int(azimuth_order), fine_nodes, fine_weights)
            error = abs(fine - coarse)
            if error > float(relative_tolerance) * max(abs(fine), 1e-300):
                raise ValueError(
                    "axisymmetric exchange quadrature did not certify; "
                    f"pair=({i},{j}), relative={error / max(abs(fine), 1e-300):.3g}")
            exchange[i, j] = exchange[j, i] = fine
            receipt[i, j] = receipt[j, i] = error

    transfer = exchange / area[:, None]
    outgoing = np.sum(transfer, axis=1)
    if np.any(outgoing > 1.0 + 16.0 * float(relative_tolerance)):
        worst = int(np.argmax(outgoing))
        raise ValueError(
            f"axisymmetric exchange exceeds unit row closure on band {worst} "
            f"(outgoing={outgoing[worst]:.6f})")
    excess = outgoing > 1.0
    if np.any(excess):
        transfer[excess] = transfer[excess] / outgoing[excess, None]
        outgoing = np.sum(transfer, axis=1)
    escape = np.maximum(1.0 - outgoing, 0.0)
    escape += 1.0 - (outgoing + escape)
    escape = np.maximum(escape, 0.0)

    digest = sha256()
    digest.update(b"petch.axisymmetric-band-exchange.v1\0")
    for name, value in (("z0", z0), ("z1", z1), ("r0", r0), ("r1", r1),
                        ("exchange", exchange)):
        array = np.ascontiguousarray(value, dtype="<f8")
        digest.update(name.encode("ascii") + b"\0")
        digest.update(np.asarray(array.shape, dtype="<i8").tobytes())
        digest.update(array.tobytes())
    digest.update(np.asarray([
        bands_per_segment, azimuth_order, generator_order,
        relative_tolerance, geometry_tolerance], dtype="<f8").tobytes())
    return {
        "band_z": np.stack((z0, z1), axis=1),
        "band_r": np.stack((r0, r1), axis=1),
        "band_area": area,
        "exchange_area": exchange,
        "transfer_fraction": transfer,
        "escape_fraction": escape,
        "quadrature_receipt": receipt,
        "fingerprint": digest.hexdigest(),
    }


def clausing_transmission(profile, *, bands_per_segment=6, azimuth_order=24,
                          generator_order=6, relative_tolerance=1.0e-4):
    """Free-molecular transmission probability through an open-ended body of revolution.

    Molecules enter through the TOP opening with a cosine distribution, reflect diffusely
    from the walls, and exit through either end.  Returns the fraction leaving through the
    BOTTOM opening -- the Clausing factor for a tube when the profile is a cylinder.  The
    entrance flux onto each wall band uses the same exact ring kernel with the top opening
    treated as a diffuse disk source.
    """
    operator = build_axisymmetric_band_exchange(
        profile, bands_per_segment=bands_per_segment, azimuth_order=azimuth_order,
        generator_order=generator_order, relative_tolerance=relative_tolerance)
    band_z = operator["band_z"]
    band_r = operator["band_r"]
    transfer = operator["transfer_fraction"]
    count = len(band_z)

    top_z = float(profile.z[-1])
    top_r = float(profile.r[-1])
    bottom_z = float(profile.z[0])
    bottom_r = float(profile.r[0])

    def disk_to_band_fraction(disk_z, disk_r, facing_down):
        """View fraction from the (diffuse) disk to each wall band, exact ring kernel."""
        nodes, weights = np.polynomial.legendre.leggauss(int(generator_order) * 2)
        kink_z = np.asarray(profile.z[1:-1], dtype=float)
        kink_r = np.asarray(profile.r[1:-1], dtype=float)
        fractions = np.zeros(count)
        nz = -1.0 if facing_down else 1.0
        for index in range(count):
            z0b, z1b = band_z[index]
            r0b, r1b = band_r[index]
            length = sqrt((z1b - z0b) ** 2 + (r1b - r0b) ** 2)
            tangent_r = (r1b - r0b) / length
            tangent_z = (z1b - z0b) / length
            n_r, n_z = -tangent_z, tangent_r
            total = 0.0
            for node_s, weight_s in zip(nodes, weights):
                rs = disk_r * 0.5 * (node_s + 1.0)
                inner = 0.0
                for node_t, weight_t in zip(nodes, weights):
                    tj = 0.5 * (node_t + 1.0)
                    rt = r0b + tj * (r1b - r0b)
                    zt = z0b + tj * (z1b - z0b)
                    threshold = _visible_cosine_threshold(
                        rs, disk_z, rt, zt, kink_z, kink_r, 1e-12)
                    if threshold <= -1.0:
                        continue
                    from math import acos as _acos
                    phi_limit = pi if threshold >= 1.0 else _acos(max(-1.0, threshold))
                    if phi_limit <= 0.0:
                        continue
                    kernel = _pair_kernel_integral(
                        rs, disk_z, 0.0, nz, rt, zt, n_r, n_z, phi_limit,
                        int(azimuth_order))
                    inner += weight_t * kernel * rt * length * 0.5
                total += weight_s * inner * rs * disk_r * 0.5
            fractions[index] = 2.0 * pi * total / (pi * disk_r * disk_r)
        return fractions

    # Source: unit diffuse influx through the top disk, facing down into the cavity.
    wall_direct = disk_to_band_fraction(top_z, top_r, facing_down=True)
    direct_bottom_disk = max(0.0, 1.0 - wall_direct.sum()
                             ) if bottom_r > 0.0 else max(0.0, 1.0 - wall_direct.sum())
    # Diffuse wall bounces: escape split between top and bottom openings needs per-band
    # escape routing.  Escape fraction from the operator is TOTAL escape (both ends); split
    # it by the same exact disk kernels from each band.
    bottom_view = (disk_to_band_fraction(bottom_z, max(bottom_r, 1e-12),
                                         facing_down=False)
                   * (pi * max(bottom_r, 1e-12) ** 2)
                   / operator["band_area"]) if bottom_r > 0.0 else np.zeros(count)
    escape = operator["escape_fraction"]
    bottom_share = np.minimum(bottom_view, escape)
    # Absorbing-boundary Neumann series: wall load b, bounce matrix T (row-stochastic up
    # to escape); transmitted = direct + b^T (I - T^T)^-1 bottom_share.
    identity = np.eye(count)
    resolvent = np.linalg.solve(identity - transfer.T, bottom_share)
    transmitted = float(direct_bottom_disk + wall_direct @ resolvent)
    return transmitted, operator


def santeler_transmission(aspect_ratio):
    """Santeler's closed-form Clausing-factor approximation (<0.7 percent error)."""
    l_over_r = 2.0 * float(aspect_ratio)
    return 1.0 / (1.0 + (3.0 * l_over_r / 8.0)
                  * (1.0 + 1.0 / (3.0 * (1.0 + l_over_r / 7.0))))
