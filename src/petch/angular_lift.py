"""Axisymmetric angular lift: recover a polar measure from a planar marginal.

A published IEAD from an axisymmetric reactor model (or a plane-resolving
analyser) supplies the ion angle *projected into one plane* — the signed angle
whose tangent is ``v_x / v_z``.  Lifting that measurement to the 3-D angular
distribution a feature-scale transport solver needs is not a relabelling: the
published measure has already been marginalized over the unresolved transverse
component, and the lift must invert that marginalization.

Geometry.  For an ion with polar angle ``theta`` and azimuth ``phi`` about the
surface normal, the projected planar angle ``theta_x`` obeys

    tan(theta_x) = tan(theta) * cos(phi)                                    (1)

exactly.  With ``phi`` uniform (the axisymmetric closure), a shell of ions at
fixed ``theta`` projects onto ``|theta_x| <= theta`` with the arcsine law

    F_theta(a) = P(|theta_x| <= a) = (2/pi) * arcsin(tan a / tan theta)     (2)

so the planar marginal ``m`` and the polar measure ``p`` are related by the
Abel-type integral equation

    m(a) = integral_a^inf  p(theta) / (pi * sqrt(tan^2 theta - tan^2 a))
                            * sec^2(a) d(theta)                             (3)

whose small-angle form is the classical Abel transform pair.  This module
inverts (3) numerically by *onion peeling*: the outermost bin can only be fed
by the outermost shell, which fixes that shell's mass; its computed spill into
every inner bin is subtracted and the peel proceeds inward.  The kernel is the
analytic bin mass of (2), so the discretization error comes only from placing
each shell at its bin centre, not from the projection law.

Two consequences are worth stating because they are exact and independent of
the distribution's shape:

* **The second-moment identity.** Squaring (1) and averaging over uniform
  azimuth gives ``E[tan^2 theta_x] = E[tan^2 theta] / 2`` for *any*
  axisymmetric measure — Gaussian or not.  The polar width therefore exceeds
  the planar width by exactly sqrt(2) in rms.  Identifying the published planar
  angle with the polar angle (the pre-2026-08-02 closure) discards that factor.
* **Support is preserved.** The largest projected angle equals the largest
  polar angle, since a shell projects onto its own radius at ``phi = 0``.  The
  lift therefore reweights the published angular support; it never extends it.

The inversion is applied to the angular *marginal* and realized as a per-bin
weight multiplier, so nodes keep their own energies and angles.  This is a
declared limitation: it lifts the angular distribution uniformly in energy,
which is what a digitized point cloud (a few nodes per energy row) can support.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "shell_projection_mass",
    "invert_planar_marginal",
    "axisymmetric_polar_weights",
]


def shell_projection_mass(polar_deg, lower_deg, upper_deg):
    """Mass a shell at ``polar_deg`` projects into folded bin [lower, upper].

    Implements the arcsine projection law (2).  Angles in degrees; the bin is
    on the folded axis ``|theta_x| >= 0``.  Broadcasting follows numpy rules.
    """
    polar = np.deg2rad(np.asarray(polar_deg, dtype=float))
    tan_polar = np.tan(polar)
    with np.errstate(divide="ignore", invalid="ignore"):
        def cumulative(edge):
            ratio = np.tan(np.deg2rad(np.asarray(edge, dtype=float))) / tan_polar
            return (2.0 / np.pi) * np.arcsin(np.clip(ratio, 0.0, 1.0))

        mass = cumulative(upper_deg) - cumulative(lower_deg)
    return np.where(np.isfinite(mass), mass, 0.0)


def invert_planar_marginal(edges_deg, planar_mass):
    """Onion-peel the folded planar histogram into a polar shell measure.

    ``edges_deg`` are ascending folded bin edges starting at 0; ``planar_mass``
    holds the mass in each bin.  Returns ``(shell_deg, shell_mass,
    diagnostics)`` with one shell per bin, placed at the bin centre.

    The peel is exact back-substitution on a triangular system (a shell feeds
    only bins at or below its own radius), so no least-squares or
    regularization parameter enters.  Abel inversion amplifies noise in the
    outer bins; any resulting negative shell mass is clamped to zero and the
    clamped fraction reported rather than silently absorbed.
    """
    edges = np.asarray(edges_deg, dtype=float)
    mass = np.asarray(planar_mass, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("edges_deg must be a 1-D array of at least two edges")
    if np.any(np.diff(edges) <= 0.0) or edges[0] < 0.0:
        raise ValueError("edges_deg must ascend from a non-negative first edge")
    if mass.shape != (edges.size - 1,):
        raise ValueError("planar_mass must hold one entry per bin")

    centres = 0.5 * (edges[:-1] + edges[1:])
    residual = mass.copy()
    shell = np.zeros_like(mass)
    for j in range(mass.size - 1, -1, -1):
        own = float(shell_projection_mass(centres[j], edges[j], edges[j + 1]))
        if own <= 0.0:
            continue
        shell[j] = residual[j] / own
        if j:
            residual[:j] -= shell[j] * shell_projection_mass(
                centres[j], edges[:j], edges[1:j + 1])

    negative = float(-shell[shell < 0.0].sum())
    clamped = np.clip(shell, 0.0, None)
    total = float(clamped.sum())
    source_total = float(mass.sum())
    if total > 0.0 and source_total > 0.0:
        clamped = clamped * (source_total / total)
    diagnostics = {
        "bin_count": int(mass.size),
        "clamped_negative_mass": negative,
        "clamped_negative_fraction": (
            negative / abs(source_total) if source_total else 0.0),
        "planar_total": source_total,
        "polar_total": float(clamped.sum()),
    }
    return centres, clamped, diagnostics


def axisymmetric_polar_weights(signed_angle_deg, weight, *, bin_deg=0.25,
                               max_clamped_fraction=0.05):
    """Per-node weights realizing the axisymmetric lift of a planar IEAD.

    Nodes keep their published angle (used as the polar angle magnitude, whose
    support the lift preserves) and are reweighted by the ratio of inverted
    polar mass to planar mass in their folded angular bin.  Returns
    ``(polar_deg, lifted_weight, diagnostics)`` with ``lifted_weight`` summing
    to the input total.

    ``bin_deg`` must match the angular resolution the source measure actually
    carries.  Abel inversion amplifies noise like the reciprocal of the bin
    width, so peeling a digitized point cloud finer than its own grid does not
    buy resolution — it manufactures negative shells.  On the Krüger Figure-4
    digitization the clamped fraction is exactly zero at the source grid of
    0.25 degrees and reaches 2.3 (i.e. the inversion has diverged) at 0.10.
    Exceeding ``max_clamped_fraction`` is therefore refused rather than
    silently clamped.
    """
    angle = np.abs(np.asarray(signed_angle_deg, dtype=float))
    mass = np.asarray(weight, dtype=float)
    if angle.shape != mass.shape:
        raise ValueError("signed_angle_deg and weight must share a shape")
    if not np.isfinite(bin_deg) or bin_deg <= 0.0:
        raise ValueError("bin_deg must be positive and finite")
    if np.any(angle >= 90.0):
        raise ValueError("planar angles must lie strictly inside +/-90 degrees")
    total = float(mass.sum())
    if angle.size == 0 or total <= 0.0 or float(angle.max()) <= 0.0:
        return angle, mass, {"bin_count": 0, "clamped_negative_fraction": 0.0,
                             "planar_total": total, "polar_total": total,
                             "degenerate": True}

    edge_count = int(np.ceil(float(angle.max()) / float(bin_deg))) + 1
    edges = np.arange(edge_count, dtype=float) * float(bin_deg)
    edges[-1] = max(edges[-1], float(angle.max()) * (1.0 + 1e-12))
    index = np.clip(np.searchsorted(edges, angle, side="right") - 1,
                    0, edges.size - 2)
    planar_mass = np.bincount(index, weights=mass, minlength=edges.size - 1)

    _, shell_mass, diagnostics = invert_planar_marginal(edges, planar_mass)
    clamped_fraction = float(diagnostics["clamped_negative_fraction"])
    if (max_clamped_fraction is not None
            and clamped_fraction > float(max_clamped_fraction)):
        raise ValueError(
            f"axisymmetric lift diverged: {clamped_fraction:.3g} of the mass "
            f"peeled negative at bin_deg={bin_deg:g}, above the "
            f"{float(max_clamped_fraction):.3g} guard — the source measure does "
            "not carry that angular resolution")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(planar_mass > 0.0, shell_mass / planar_mass, 0.0)
    lifted = mass * ratio[index]
    lifted_total = float(lifted.sum())
    if lifted_total > 0.0:
        lifted = lifted * (total / lifted_total)
    diagnostics = dict(diagnostics)
    diagnostics["bin_deg"] = float(bin_deg)
    diagnostics["degenerate"] = False
    return angle, lifted, diagnostics
