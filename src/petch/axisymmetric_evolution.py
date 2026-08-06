"""Axisymmetric (r, z) profile evolution driver for high-aspect-ratio holes.

Phase-2 of ``HOLE_STUDY_PLAN_2026-08-05.md``.  Phase 1 characterised transport
through a *frozen* hole; this module lets the profile move, coupling the same
validated pieces:

* geometry -- a single-valued generator ``r(depth)`` on a fixed axial band grid
  plus a moving floor plane (see *Representation* below);
* thermal neutrals -- the exact closed-form cylinder band-exchange enclosure
  (:func:`petch.axisymmetric_exchange_3d.build_cylinder_band_exchange`) closed
  with a floor disk and a mouth aperture, solved for multi-bounce diffuse
  transport at the mechanism's own per-species reaction probability;
* energetic particles -- the deterministic specular cascade, production
  reaction rule verbatim (``react = clip(0.9*kress(cos), 0, 1)``, continuing
  weight ``1 - react``, Eq. 2.34 retention), binned per band and per incidence
  cosine so the surface model evaluates its nonlinear laws per event;
* thermalised return (E8) -- cascade weight that loses specular retention is
  reborn into the neutral ledger at the band where it thermalised, at a
  *declared* fluorocarbon fraction (unpublished for every reactor in the
  corpus; swept, never fitted);
* chemistry -- the same :class:`petch.mixed_layer_mechanism.MixedLayerMechanism`
  the trench validation uses, one face per band plus one for the floor.

Representation
--------------
The front is tracked directly as a single-valued generator ``r(depth)`` on a
fixed axial grid, with the floor as a separate scalar depth.  Chosen over a
level set in (r, z) because:

1. it *is* the input the validated transport already consumes
   (:class:`~petch.axisymmetric_exchange_3d.AxisymmetricProfile`), so no
   conversion, reinitialisation or marching-cubes step sits between the moving
   front and the operator that was benchmarked against Clausing;
2. band areas, and hence the volume-vs-ledger conservation check, are exact
   closed-form annuli rather than a reconstructed iso-surface;
3. the etch front of a HAR hole is single-valued in ``r(depth)`` except under
   true undercut.

**Declared representation limit:** a re-entrant (overhanging) wall, where one
depth carries two radii, cannot be represented.  Bowing -- a larger radius at
mid-depth -- is single-valued and *is* representable.  Undercut is out of scope
for v1 and the driver has no way to detect it, so it is a declared limitation,
not a guarded one.

**Declared transport envelope:** both transport channels are exact for a
straight cylinder.  As the profile evolves the driver measures its own
straightness (max relative radius deviation across exposed bands) and evaluates
the operators on the area-weighted mean radius.  Past a declared tolerance the
run *stops with a receipt* rather than extrapolating an operator outside the
regime it was gated in -- the general body-of-revolution operator does not yet
certify its own self-pair quadrature (``HOLE_STUDY_RESULTS_2026-08-05.md`` §2),
so there is nothing honest to fall back to.  The stop depth is itself a
reported observable: it is where taper stops being a perturbation.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import pi

import numpy as np

from .axisymmetric_exchange_3d import (
    _coaxial_disk_factor,
    build_cylinder_band_exchange,
)
from .mixed_layer_mechanism import MixedLayerSurfaceState
from .surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes

__all__ = [
    "HoleEnclosure",
    "HoleGeometry",
    "HoleEvolutionState",
    "build_hole_enclosure",
    "solve_diffuse_hole_delivery",
    "cascade_hole_delivery",
    "advance_hole_step",
    "evolve_hole",
]

#: Eq. 2.34 retention constants, ``boundary_transport_3d`` verbatim.
_E_THERMAL_SPECULAR_EV = 100.0
_THETA_CRITICAL_DEG = 70.0
#: Production angular form (Kress B = 9.3), ``boundary_transport_3d`` verbatim.
_KRESS_B = 9.3


def _kress(cosine):
    return np.maximum((1.0 + _KRESS_B * (1.0 - cosine ** 2)) * cosine, 0.0)


# --- geometry -----------------------------------------------------------------


@dataclass(frozen=True)
class HoleGeometry:
    """Discrete axisymmetric hole.

    ``radius[i]`` is the wall radius of band ``i``, which spans depths
    ``[i*band_height, (i+1)*band_height]`` below the mouth plane.  The last
    exposed band is truncated at ``floor_depth``.  Lengths are metres.
    """

    band_height: float
    radius: np.ndarray
    floor_depth: float
    floor_radius: float

    def __post_init__(self):
        radius = np.asarray(self.radius, dtype=float).copy()
        if (self.band_height <= 0.0 or radius.ndim != 1 or radius.size < 1
                or np.any(~np.isfinite(radius)) or np.any(radius <= 0.0)
                or not np.isfinite(self.floor_depth) or self.floor_depth <= 0.0
                or not np.isfinite(self.floor_radius) or self.floor_radius <= 0.0):
            raise ValueError("invalid hole geometry")
        expected = int(np.ceil(self.floor_depth / self.band_height - 1e-12))
        if radius.size != max(expected, 1):
            raise ValueError(
                f"radius array holds {radius.size} bands, floor depth implies {expected}")
        radius.setflags(write=False)
        object.__setattr__(self, "radius", radius)

    @classmethod
    def straight(cls, radius, depth, band_height):
        """A straight-walled hole of the given radius and depth."""
        count = max(int(np.ceil(float(depth) / float(band_height) - 1e-12)), 1)
        return cls(float(band_height), np.full(count, float(radius)),
                   float(depth), float(radius))

    @property
    def band_count(self):
        return int(self.radius.size)

    @property
    def band_lower_depth(self):
        """Depth of the top edge of each band (towards the mouth)."""
        return self.band_height * np.arange(self.band_count, dtype=float)

    @property
    def band_upper_depth(self):
        """Depth of the bottom edge of each band, truncated at the floor."""
        edges = self.band_height * (np.arange(self.band_count, dtype=float) + 1.0)
        return np.minimum(edges, self.floor_depth)

    @property
    def band_exposed_height(self):
        return np.maximum(self.band_upper_depth - self.band_lower_depth, 0.0)

    @property
    def band_area(self):
        """Lateral area of each exposed wall band (m^2)."""
        return 2.0 * pi * self.radius * self.band_exposed_height

    @property
    def floor_area(self):
        return pi * self.floor_radius ** 2

    @property
    def mean_radius(self):
        area = self.band_area
        total = area.sum()
        if total <= 0.0:
            return float(self.floor_radius)
        return float((self.radius * area).sum() / total)

    @property
    def straightness_deviation(self):
        """Max relative radius deviation across exposed bands and the floor rim."""
        mean = self.mean_radius
        values = np.concatenate([self.radius, [self.floor_radius]])
        return float(np.max(np.abs(values - mean)) / mean)

    @property
    def aspect_ratio(self):
        return float(self.floor_depth / (2.0 * self.mean_radius))

    def solid_volume_removed(self, initial):
        """Solid volume (m^3) removed relative to ``initial`` (same grid)."""
        def swept(geometry):
            wall = float(
                (pi * geometry.radius ** 2 * geometry.band_exposed_height).sum())
            return wall
        return swept(self) - swept(initial)


# --- diffuse (thermal) transport ----------------------------------------------


@dataclass(frozen=True)
class HoleEnclosure:
    """Closed-form exchange factors for a straight cylinder with a floor disk.

    Faces are ordered ``[band 0 .. band N-1, floor]``; band 0 is at the mouth.
    ``mouth_to_face`` is the fraction of a cosine-distributed influx through the
    mouth aperture that arrives directly at each face.  ``face_to_face[j, i]``
    is the fraction of diffuse emission from ``j`` arriving at ``i``;
    ``face_to_mouth`` is the escaping remainder.
    """

    radius: float
    depth: float
    area: np.ndarray
    face_to_face: np.ndarray
    face_to_mouth: np.ndarray
    mouth_to_face: np.ndarray
    mouth_direct_escape: float
    closure_residual: float


def build_hole_enclosure(radius, depth, band_edges_depth):
    """Assemble the exact cylinder enclosure: wall bands + floor + mouth aperture.

    ``band_edges_depth`` are increasing depths from the mouth (0) to the floor.
    Every factor is closed form -- the wall-wall and wall-disk blocks come from
    the gated :func:`build_cylinder_band_exchange` algebra, the disk-disk terms
    from the coaxial-disk factor, and floor-to-wall by reciprocity.
    """
    radius = float(radius)
    depth = float(depth)
    edges = np.asarray(band_edges_depth, dtype=float)
    if (radius <= 0.0 or depth <= 0.0 or edges.ndim != 1 or edges.size < 2
            or np.any(np.diff(edges) <= 0.0)
            or abs(edges[0]) > 1e-12 * depth or abs(edges[-1] - depth) > 1e-9 * depth):
        raise ValueError("invalid hole enclosure grid")
    count = edges.size - 1
    # The exchange operator works in height above the floor; our grid is depth
    # from the mouth.  Reverse so band 0 (mouth) is the topmost strip.
    z_edges = np.sort(depth - edges)
    operator = build_cylinder_band_exchange(radius, z_edges)
    # Row/column k of the operator is the strip [z_k, z_k+1]; band i (depth
    # order) is strip (count-1-i).
    order = np.arange(count)[::-1]
    wall_area = operator["band_area"][order]
    wall_to_wall = operator["transfer_fraction"][np.ix_(order, order)]
    wall_to_floor = operator["escape_bottom"][order]
    wall_to_mouth = operator["escape_top"][order]

    disk_area = pi * radius * radius

    def disk_disk(gap):
        if gap <= 0.0:
            return 1.0
        return _coaxial_disk_factor(radius, radius, gap)

    floor_to_mouth = disk_disk(depth)
    # Reciprocity: A_floor F(floor->band) = A_band F(band->floor).
    floor_to_wall = wall_to_floor * wall_area / disk_area
    # Cosine influx through the mouth: same section algebra as the Clausing
    # solve -- direct to the floor plus the load deposited on each band.
    mouth_direct_floor = disk_disk(depth)
    lower = np.minimum(edges[:-1], depth)
    upper = np.minimum(edges[1:], depth)
    mouth_to_wall = np.array(
        [disk_disk(lower[i]) - disk_disk(upper[i]) for i in range(count)])

    area = np.concatenate([wall_area, [disk_area]])
    face_to_face = np.zeros((count + 1, count + 1))
    face_to_face[:count, :count] = wall_to_wall
    face_to_face[:count, count] = wall_to_floor
    face_to_face[count, :count] = floor_to_wall
    face_to_mouth = np.concatenate([wall_to_mouth, [floor_to_mouth]])
    mouth_to_face = np.concatenate([mouth_to_wall, [mouth_direct_floor]])

    closure = np.max(np.abs(face_to_face.sum(axis=1) + face_to_mouth - 1.0))
    influx_closure = abs(mouth_to_face.sum() - 1.0)
    residual = float(max(closure, influx_closure))
    if residual > 1e-7:
        raise ValueError(
            f"hole enclosure closure failed: max |defect| = {residual:.3e}")
    return HoleEnclosure(
        radius=radius, depth=depth, area=area, face_to_face=face_to_face,
        face_to_mouth=face_to_mouth, mouth_to_face=mouth_to_face,
        mouth_direct_escape=0.0, closure_residual=residual)


def solve_diffuse_hole_delivery(enclosure, sticking, mouth_rate, born_rate=None):
    """Multi-bounce diffuse delivery through the enclosure.

    ``sticking`` is the per-face loss probability per collision (the mechanism's
    own ``neutral_reaction_probability``), ``mouth_rate`` the total entering
    rate (particles/s) and ``born_rate`` an optional per-face birth rate, used
    for the thermalised (E8) return: a reborn particle is *emitted* from its
    band rather than deposited there.

    Returns ``(arrival_rate, absorbed_rate, escaped_rate)`` where the arrival
    rate is per face and dimensionally particles/s (divide by area for the flux
    density the surface model consumes).
    """
    stick = np.clip(np.asarray(sticking, dtype=float), 0.0, 1.0)
    faces = enclosure.area.size
    if stick.shape != (faces,):
        raise ValueError("sticking must be one value per face")
    born = (np.zeros(faces) if born_rate is None
            else np.asarray(born_rate, dtype=float))
    if born.shape != (faces,) or np.any(born < 0.0):
        raise ValueError("invalid born rate")
    direct = float(mouth_rate) * enclosure.mouth_to_face + born @ enclosure.face_to_face
    reflect = enclosure.face_to_face.T * (1.0 - stick)[None, :]
    arrival = np.linalg.solve(np.eye(faces) - reflect, direct)
    absorbed = stick * arrival
    escaped = float(
        ((1.0 - stick) * arrival + born) @ enclosure.face_to_mouth)
    return arrival, absorbed, escaped


# --- energetic (specular cascade) transport -----------------------------------


def cascade_hole_delivery(radius, depth, band_edges_depth, iadf, energy_eV, *,
                          max_bounces=8, n_polar=192, n_azimuth=64,
                          n_radial=24, n_cosine_bins=8, minimum_weight=1e-4):
    """Deterministic specular cascade in a straight cylinder, resolved per band.

    Same construction as ``scripts/hole_study_phase1.cascade_delivery`` (a
    specular reflection off a cylinder preserves polar angle and impact
    parameter, so a ray's strikes are equally spaced in depth and the cascade is
    exact algebra per ray), extended to report per-band *arrival* rate resolved
    in incidence cosine -- what the surface model needs -- rather than only the
    reacting share.

    Returns a dict of per-face arrival histograms normalised to unit flux
    entering the mouth aperture.  ``wall_hist[i, k]`` is the arriving weight at
    band ``i`` in cosine bin ``k``; ``floor_hist[k]`` the same at the floor.
    """
    depth = float(depth)
    radius = float(radius)
    edges = np.asarray(band_edges_depth, dtype=float)
    count = edges.size - 1
    bins = int(n_cosine_bins)
    cos_edges = np.linspace(0.0, 1.0, bins + 1)
    cos_centres = 0.5 * (cos_edges[:-1] + cos_edges[1:])

    nodes, weights = np.polynomial.legendre.leggauss(int(n_radial))
    entry_r = np.sqrt(0.5 * (nodes + 1.0) * radius ** 2)
    entry_w = weights / weights.sum()
    phi = (np.arange(int(n_azimuth)) + 0.5) * np.pi / int(n_azimuth)
    phi_w = np.full(phi.shape, 1.0 / int(n_azimuth))
    polar_deg, polar_w = iadf.polar_quadrature(float(energy_eV), n_polar=int(n_polar))
    polar = np.deg2rad(polar_deg)

    wall_hist = np.zeros((count, bins))
    floor_hist = np.zeros(bins)
    thermalised = np.zeros(count)
    absorbed = np.zeros(count)
    direct_bottom = 0.0
    cascaded_bottom = 0.0
    generations = 0

    rr = entry_r[:, None]
    ww = entry_w[:, None] * phi_w[None, :]
    p_dot_d = -rr * np.cos(phi)[None, :]
    impact = np.abs(rr * np.sin(phi)[None, :])
    half_chord = np.sqrt(np.maximum(radius ** 2 - impact ** 2, 0.0))
    first_transverse = p_dot_d + half_chord
    chord = 2.0 * half_chord
    cos_geometry = half_chord / radius
    cos_critical = np.cos(np.deg2rad(_THETA_CRITICAL_DEG))

    for angle, mass in zip(polar, polar_w):
        if mass <= 0.0:
            continue
        tangent = np.tan(angle)
        axial_cos = float(np.cos(angle))
        floor_bin = min(int(np.digitize(axial_cos, cos_edges) - 1), bins - 1)
        if tangent <= 0.0:
            direct_bottom += float(mass)
            floor_hist[floor_bin] += float(mass)
            continue
        sin_theta = np.sin(angle)
        cosine = np.clip(sin_theta * cos_geometry, 0.0, 1.0)
        wall_bin = np.clip(np.digitize(cosine, cos_edges) - 1, 0, bins - 1)
        react = np.clip(0.9 * _kress(cosine), 0.0, 1.0)
        continue_weight = 1.0 - react
        retained_full = (cosine < cos_critical) & (float(energy_eV) > _E_THERMAL_SPECULAR_EV)

        z_first = first_transverse / tangent
        step = chord / tangent
        weight = np.full(rr.shape, 1.0) * ww * float(mass)
        direct = z_first >= depth
        direct_bottom += float(weight[direct].sum())
        floor_hist[floor_bin] += float(weight[direct].sum())
        weight = np.where(direct, 0.0, weight)
        z_hit = np.where(direct, np.inf, z_first)

        for bounce in range(int(max_bounces)):
            alive = (weight > 0.0) & np.isfinite(z_hit) & (z_hit < depth)
            if not np.any(alive):
                break
            generations = max(generations, bounce + 1)
            band = np.clip(np.digitize(z_hit[alive], edges) - 1, 0, count - 1)
            np.add.at(wall_hist, (band, wall_bin[alive]), weight[alive])
            np.add.at(absorbed, band, weight[alive] * react[alive])
            surviving = weight[alive] * continue_weight[alive]
            keep = retained_full[alive] & (surviving > minimum_weight * weight[alive])
            np.add.at(thermalised, band, np.where(keep, 0.0, surviving))
            new_weight = np.zeros_like(weight)
            new_weight[alive] = np.where(keep, surviving, 0.0)
            weight = new_weight
            z_next = np.where(np.isfinite(z_hit), z_hit + step, np.inf)
            arrived = (weight > 0.0) & (z_next >= depth)
            landed = float(weight[arrived].sum())
            cascaded_bottom += landed
            floor_hist[floor_bin] += landed
            weight = np.where(arrived, 0.0, weight)
            z_hit = np.where(arrived, np.inf, z_next)
        remaining = weight > 0.0
        if np.any(remaining):
            band = np.clip(np.digitize(np.where(np.isfinite(z_hit), z_hit, depth),
                                       edges) - 1, 0, count - 1)
            np.add.at(thermalised, band[remaining], weight[remaining])

    total_bottom = float(direct_bottom + cascaded_bottom)
    return {
        "cosine_bin_centres": cos_centres,
        "wall_hist": wall_hist,
        "floor_hist": floor_hist,
        "thermalised_per_band": thermalised,
        "absorbed_per_band": absorbed,
        "direct_bottom": float(direct_bottom),
        "cascaded_bottom": float(cascaded_bottom),
        "total_bottom": total_bottom,
        "bounce_generations": int(generations),
        # Weight budget: every entering particle either reaches the floor, is
        # consumed by the reacting share at a wall strike, or thermalises
        # (including truncation at the bounce cap).
        "closure_residual": float(
            total_bottom + absorbed.sum() + thermalised.sum() - 1.0),
    }


# --- evolution ----------------------------------------------------------------


@dataclass(frozen=True)
class HoleEvolutionState:
    """Geometry plus surface chemistry state (bands first, floor last)."""

    geometry: HoleGeometry
    surface: MixedLayerSurfaceState
    time_s: float = 0.0


def _grow_state(surface, count):
    """Extend a per-face state to ``count`` bands + floor, new bands bare."""
    fields = ("n_c_film", "n_f_film", "n_si", "n_o", "n_c", "n_f",
              "n_xl_film", "n_act", "removed_formula_units_m2")
    values = [np.asarray(getattr(surface, name), dtype=float) for name in fields]
    current = values[0].size - 1
    if count == current:
        return surface
    grown = []
    for value in values:
        bands, floor = value[:-1], value[-1:]
        pad = np.zeros(count - current)
        grown.append(np.concatenate([bands, pad, floor]))
    return MixedLayerSurfaceState(*grown)


def _build_fluxes(mechanism, state, enclosure, geometry, mouth_flux_m2_s,
                  cascade, ion_flux_m2_s, energy_eV, thermalized_return_fraction,
                  thermalized_stoichiometry):
    """Assemble per-face :class:`SurfaceFluxes` from the two transport channels."""
    area = enclosure.area
    mouth_area = pi * enclosure.radius ** 2
    probability = mechanism.neutral_reaction_probability(state)
    # E8: thermalised cascade weight is reborn as radicals at its band, at a
    # declared fluorocarbon share (unpublished; swept, never fitted).
    born_total = np.zeros(area.size)
    born_total[:-1] = (np.asarray(cascade["thermalised_per_band"], dtype=float)
                       * float(ion_flux_m2_s) * mouth_area
                       * float(thermalized_return_fraction))

    neutral = {}
    escaped = {}
    balance = 0.0
    for name, flux in mouth_flux_m2_s.items():
        stick = np.asarray(probability.get(name, 0.0), dtype=float)
        stick = np.broadcast_to(stick, area.shape).astype(float)
        born = born_total * float(thermalized_stoichiometry.get(name, 0.0))
        entering = float(flux) * mouth_area
        arrival, absorbed, escape = solve_diffuse_hole_delivery(
            enclosure, stick, entering, born_rate=born)
        neutral[name] = arrival / area
        escaped[name] = escape
        supply = entering + float(born.sum())
        if supply > 0.0:
            balance = max(balance, abs(
                float(absorbed.sum()) + escape - supply) / supply)

    scale = float(ion_flux_m2_s) * mouth_area
    faces = area.size
    centres = cascade["cosine_bin_centres"]
    hist = np.vstack([cascade["wall_hist"], cascade["floor_hist"][None, :]])
    rate = hist * scale
    face_index, bin_index = np.nonzero(rate > 0.0)
    events = FaceResolvedEnergeticFlux(
        "ions", faces, face_index,
        rate[face_index, bin_index] / area[face_index],
        np.full(face_index.shape, float(energy_eV)),
        np.clip(centres[bin_index], 1e-6, 1.0))
    return SurfaceFluxes(neutral, (events,)), escaped, born_total, balance


def advance_hole_step(state, mechanism, *, mouth_flux_m2_s, ion_flux_m2_s,
                      iadf, energy_eV, dt_s, thermalized_return_fraction=0.0,
                      thermalized_stoichiometry=None, max_bounces=8,
                      cascade_kwargs=None):
    """Advance the hole one step; returns ``(new_state, record)``."""
    geometry = state.geometry
    edges = np.concatenate([geometry.band_lower_depth, [geometry.floor_depth]])
    enclosure = build_hole_enclosure(geometry.mean_radius, geometry.floor_depth, edges)
    cascade = cascade_hole_delivery(
        geometry.mean_radius, geometry.floor_depth, edges, iadf, energy_eV,
        max_bounces=max_bounces, **(cascade_kwargs or {}))
    fluxes, escaped, born, neutral_balance = _build_fluxes(
        mechanism, state.surface, enclosure, geometry, mouth_flux_m2_s, cascade,
        ion_flux_m2_s, energy_eV, thermalized_return_fraction,
        thermalized_stoichiometry or {})
    result = mechanism.advance(state.surface, fluxes, float(dt_s), strict=False)

    etch = np.asarray(result.etch_velocity_m_s, dtype=float)
    growth = np.asarray(result.normal_growth_velocity_m_s, dtype=float)
    normal = (etch - growth) * float(dt_s)
    radius = geometry.radius + normal[:-1]
    floor_depth = geometry.floor_depth + float(normal[-1])
    # The floor disk spans the wall radius at its own depth.
    floor_radius = float(radius[-1])
    if np.any(radius <= 0.0) or floor_depth <= 0.0:
        raise ValueError("hole geometry collapsed (clogged or inverted)")

    count = max(int(np.ceil(floor_depth / geometry.band_height - 1e-12)), 1)
    if count > radius.size:
        radius = np.concatenate([radius, np.full(count - radius.size, floor_radius)])
    moved = HoleGeometry(geometry.band_height, radius, floor_depth, floor_radius)
    surface = _grow_state(result.state, count)

    removed_units = np.asarray(result.removed_bare_formula_units_m2, dtype=float)
    record = {
        "time_s": float(state.time_s + dt_s),
        "dt_s": float(dt_s),
        "floor_depth_m": float(floor_depth),
        "aspect_ratio": float(moved.aspect_ratio),
        "mean_radius_m": float(moved.mean_radius),
        "straightness_deviation": float(moved.straightness_deviation),
        "mouth_radius_m": float(radius[0]),
        "floor_etch_velocity_m_s": float(etch[-1]),
        "floor_growth_velocity_m_s": float(growth[-1]),
        "max_wall_growth_velocity_m_s": float(growth[:-1].max()),
        "enclosure_closure_residual": float(enclosure.closure_residual),
        "neutral_balance_residual": float(neutral_balance),
        "cascade_closure_residual": float(cascade["closure_residual"]),
        "cascade_bottom_delivery": float(cascade["total_bottom"]),
        "cascade_bounce_generations": int(cascade["bounce_generations"]),
        "thermalised_born_rate_s": float(born.sum()),
        "neutral_escape_rate_s": {k: float(v) for k, v in escaped.items()},
        "removed_units_m2": removed_units,
        "band_area_m2": enclosure.area,
    }
    return HoleEvolutionState(moved, surface, float(state.time_s + dt_s)), record


def evolve_hole(state, mechanism, *, mouth_flux_m2_s, ion_flux_m2_s, iadf,
                energy_eV, duration_s, dt_s, straightness_tolerance=0.02,
                max_steps=100000, **step_kwargs):
    """Time-step the hole, stopping on the declared straight-wall envelope.

    Returns ``(final_state, records, stop_reason)``.  ``stop_reason`` is
    ``"duration"`` for a completed run or
    ``"profile_left_straight_wall_envelope"`` when the measured straightness
    deviation exceeds the declared tolerance -- the point past which the exact
    cylinder operators are outside the regime they were gated in.
    """
    records = []
    reason = "duration"
    elapsed = 0.0
    initial_volume = state.geometry
    for _ in range(int(max_steps)):
        if elapsed >= float(duration_s) - 1e-15:
            break
        step = min(float(dt_s), float(duration_s) - elapsed)
        state, record = advance_hole_step(
            state, mechanism, mouth_flux_m2_s=mouth_flux_m2_s,
            ion_flux_m2_s=ion_flux_m2_s, iadf=iadf, energy_eV=energy_eV,
            dt_s=step, **step_kwargs)
        elapsed = state.time_s
        record["solid_volume_removed_m3"] = float(
            state.geometry.solid_volume_removed(initial_volume))
        records.append(record)
        if record["straightness_deviation"] > float(straightness_tolerance):
            reason = "profile_left_straight_wall_envelope"
            break
    return state, records, reason
