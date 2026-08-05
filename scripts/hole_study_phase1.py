"""Phase-1 of the HAR hole study: geometry routing, transport reference, delivery.

Series 1-3 of ``HOLE_STUDY_PLAN_2026-08-05.md``.  Everything here is
frozen-geometry and deterministic -- no Monte Carlo, no evolution, no fitting.

Series 1  transport reference: Clausing transmission vs aspect ratio through
          three independent paths (exact band algebra, the general
          body-of-revolution operator fed by the STL-derived generator, and
          Santeler's published closed form).
Series 2  ion delivery: cone acceptance of the two-component beam vs aspect
          ratio, swept over the declared tail-fraction band.
Series 3  the cascade: a deterministic specular ray cascade inside the
          cylinder, using the production reaction rule verbatim
          (``split_grazing_ion_reflection``: react = clip(0.9*kress(cos),0,1),
          continuing weight = 1 - react, Eq. 2.34 retention).  Returns bottom
          delivery split into direct ions and cascaded hot neutrals, plus the
          sidewall energy deposition profile vs depth.

Geometry convention throughout: hole diameter 1, radius R = 0.5, depth = AR,
axis +z, entrance at the top (z = depth), etch front at z = 0.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from petch.axisymmetric_exchange_3d import (  # noqa: E402
    AxisymmetricProfile,
    clausing_transmission,
    cylinder_clausing_transmission,
    santeler_transmission,
)
from petch.iadf_two_component import (  # noqa: E402
    acceptance_half_angle_deg,
    kim_2025_reference_iadf,
)
from petch.stl_import import (  # noqa: E402
    diagnose_mesh,
    extract_axisymmetric_profile,
    read_stl,
    revolved_stl_mesh,
    to_axisymmetric_profile,
    write_stl,
)

RADIUS = 0.5
#: Reference ion energy for the sweep.  Well above the Eq. 2.34 specular
#: threshold (100 eV), so retention is full at grazing incidence and the energy
#: enters only through the beam width sigma ~ E**-0.5.  Chosen inside the
#: measurement's stated 1.4-2.0 keV band (Kim 2025 conditions).
REFERENCE_ENERGY_EV = 1500.0
#: Production cascade bounce cap (boundary_transport_3d.split_grazing_ion_reflection
#: default ``max_bounces=8``).  Reported alongside the uncapped physics limit.
PRODUCTION_BOUNCE_CAP = 8
TAIL_FRACTIONS = (0.0, 0.35, 0.50, 0.65)
ASPECT_RATIOS = (50.0, 100.0, 150.0, 200.0)


# --- geometry -----------------------------------------------------------------


def build_hole_stl(aspect_ratio, path, *, n_theta=192, n_axial=64):
    """Write a straight-walled hole of the given aspect ratio as a binary STL."""
    depth = float(aspect_ratio)
    z = np.linspace(0.0, depth, int(n_axial))
    r = np.full(z.shape, RADIUS)
    mesh = revolved_stl_mesh(z, r, n_theta=int(n_theta))
    write_stl(path, mesh, binary=True, name="har_hole")
    return mesh


def geometry_receipt(aspect_ratio, path):
    """Full import chain: STL on disk -> diagnostics -> axisymmetry -> profile."""
    build_hole_stl(aspect_ratio, path)
    mesh = read_stl(path)
    diagnostics = diagnose_mesh(mesh)
    report = extract_axisymmetric_profile(mesh)
    profile = to_axisymmetric_profile(report)
    deviation = float(report.relative_deviation)
    return {
        "aspect_ratio": float(aspect_ratio),
        "triangles": int(len(mesh.triangles)),
        "watertight": bool(diagnostics.is_watertight),
        "consistently_oriented": bool(diagnostics.consistently_oriented),
        "relative_deviation": deviation,
        "facet_bound": float(report.facet_bound),
        "is_axisymmetric": bool(report.is_axisymmetric),
        "routed_to": "axisymmetric" if report.is_axisymmetric else "full_3d",
    }, profile


# --- series 1: transport reference --------------------------------------------


def general_operator_probe(aspect_ratio, band_length, azimuth_order,
                           generator_order):
    """Try the general body-of-revolution operator; report its certification.

    The straight-hole numbers come from the exact-algebra path (the gated
    route).  This probe records whether the *general* operator -- the one a
    tapered profile would need -- certifies its own quadrature, and at what
    residual.  It is never bypassed: a failure is reported, not tolerated.
    """
    count = int(round(float(aspect_ratio) / float(band_length))) + 1
    z = np.linspace(0.0, float(aspect_ratio), count)
    profile = AxisymmetricProfile(z, np.full(z.shape, RADIUS))
    try:
        value = clausing_transmission(
            profile, bands_per_segment=1, azimuth_order=int(azimuth_order),
            generator_order=int(generator_order))
        return {"certified": True, "value": float(value), "residual": None}
    except ValueError as error:
        message = str(error)
        residual = None
        if "relative=" in message:
            residual = float(message.rsplit("relative=", 1)[1].split()[0])
        return {"certified": False, "value": None, "residual": residual,
                "message": message}


def transport_reference(aspect_ratios, workdir):
    rows = []
    for aspect in aspect_ratios:
        receipt, _profile = geometry_receipt(
            aspect, Path(workdir) / f"hole_ar{int(aspect)}.stl")
        exact = cylinder_clausing_transmission(aspect)
        santeler = santeler_transmission(aspect)
        rows.append({
            **receipt,
            "clausing_exact_algebra": float(exact),
            "santeler_closed_form": float(santeler),
            "exact_vs_santeler_relative": float(abs(exact - santeler) / santeler),
        })
    return rows


# --- series 2: cone acceptance -------------------------------------------------


def cone_acceptance(aspect_ratios, tail_fractions, energy_eV):
    rows = []
    for fraction in tail_fractions:
        iadf = kim_2025_reference_iadf(tail_fraction=fraction)
        for aspect in aspect_ratios:
            alpha = float(acceptance_half_angle_deg(aspect))
            cone = float(iadf.acceptance_fraction_cone(alpha, energy_eV))
            rows.append({
                "tail_fraction": float(fraction),
                "aspect_ratio": float(aspect),
                "acceptance_half_angle_deg": alpha,
                "cone_acceptance": cone,
                "sidewall_share": 1.0 - cone,
            })
    return rows


# --- series 3: deterministic specular cascade ---------------------------------


def _kress(cosine):
    """Production angular form, ``boundary_transport_3d`` line 775 verbatim."""
    return np.maximum((1.0 + 9.3 * (1.0 - cosine ** 2)) * cosine, 0.0)


def cascade_delivery(aspect_ratio, iadf, energy_eV, *, n_polar=192,
                     n_azimuth=64, n_radial=24, max_bounces=PRODUCTION_BOUNCE_CAP,
                     depth_bins=100, minimum_weight=1e-4):
    """Deterministic ray cascade in the cylinder.

    A specular reflection off a cylinder wall preserves the polar angle and the
    impact parameter, so successive wall strikes of one ray are equally spaced
    in depth -- the cascade is exact algebra per ray, with quadrature only over
    the entry disk and the beam's angular measure.

    Returns bottom delivery (direct + cascaded), the thermalised share, and the
    sidewall energy deposition profile vs depth.
    """
    depth = float(aspect_ratio)
    # Entry position: area measure over the disk, Gauss-Legendre in s^2.
    nodes, weights = np.polynomial.legendre.leggauss(int(n_radial))
    s_squared = 0.5 * (nodes + 1.0) * RADIUS ** 2
    entry_r = np.sqrt(s_squared)
    entry_w = weights / weights.sum()
    # Direction azimuth relative to the entry radius: uniform, midpoint rule
    # over [0, pi) (the other half is the mirror image).
    phi = (np.arange(int(n_azimuth)) + 0.5) * np.pi / int(n_azimuth)
    phi_w = np.full(phi.shape, 1.0 / int(n_azimuth))
    # Polar angle: the beam's own measure.
    polar_deg, polar_w = iadf.polar_quadrature(energy_eV, n_polar=int(n_polar))
    polar = np.deg2rad(polar_deg)

    bins = np.linspace(0.0, depth, int(depth_bins) + 1)
    wall_energy = np.zeros(int(depth_bins))
    wall_rate = np.zeros(int(depth_bins))
    direct_bottom = 0.0
    cascaded_bottom = 0.0
    cascaded_bottom_energy = 0.0
    thermalised = 0.0
    generations = 0

    # Broadcast the (entry x azimuth) transverse geometry once.
    rr = entry_r[:, None]
    ww = (entry_w[:, None] * phi_w[None, :])
    p_dot_d = -rr * np.cos(phi)[None, :]
    impact = np.abs(rr * np.sin(phi)[None, :])
    half_chord = np.sqrt(np.maximum(RADIUS ** 2 - impact ** 2, 0.0))
    first_transverse = p_dot_d + half_chord      # >= 0 for a point inside
    chord = 2.0 * half_chord
    # Wall incidence cosine (from the wall normal) for unit transverse speed.
    cos_geometry = half_chord / RADIUS

    for angle, mass in zip(polar, polar_w):
        if mass <= 0.0:
            continue
        tangent = np.tan(angle)
        if tangent <= 0.0:
            direct_bottom += float(mass)
            continue
        sin_theta = np.sin(angle)
        cosine = np.clip(sin_theta * cos_geometry, 0.0, 1.0)
        react = np.clip(0.9 * _kress(cosine), 0.0, 1.0)
        continue_weight = 1.0 - react
        # Eq. 2.34: at these energies (>100 eV) and grazing incidence
        # (cos < cos(70 deg)) the continuing particle keeps full energy.
        grazing = cosine < np.cos(np.deg2rad(70.0))
        retained_full = grazing & (energy_eV > 100.0)
        # Rays whose first strike is beyond the bottom leave directly.
        z_first = first_transverse / tangent
        dz = chord / tangent
        weight = np.full(rr.shape, 1.0) * ww * float(mass)
        direct = z_first >= depth
        direct_bottom += float(weight[direct].sum())
        weight = np.where(direct, 0.0, weight)
        z_hit = np.where(direct, np.inf, z_first)
        # Depth is measured downward from the entrance.
        for bounce in range(int(max_bounces)):
            alive = (weight > 0.0) & np.isfinite(z_hit) & (z_hit < depth)
            if not np.any(alive):
                break
            generations = max(generations, bounce + 1)
            depth_from_top = z_hit[alive]
            # Energy deposited into the wall by the reacting share.
            index = np.clip(
                np.digitize(depth - depth_from_top, bins) - 1, 0,
                int(depth_bins) - 1)
            deposited = weight[alive] * react[alive]
            np.add.at(wall_rate, index, deposited)
            np.add.at(wall_energy, index, deposited * energy_eV)
            surviving = weight[alive] * continue_weight[alive]
            # Particles that do not keep specular energy leave the cascade.
            keep = retained_full[alive] & (
                surviving > minimum_weight * weight[alive])
            thermalised += float(surviving[~keep].sum())
            new_weight = np.zeros_like(weight)
            new_weight[alive] = np.where(keep, surviving, 0.0)
            weight = new_weight
            z_next = np.where(np.isfinite(z_hit), z_hit + dz, np.inf)
            arrived = (weight > 0.0) & (z_next >= depth)
            cascaded_bottom += float(weight[arrived].sum())
            cascaded_bottom_energy += float(
                (weight[arrived] * energy_eV).sum())
            weight = np.where(arrived, 0.0, weight)
            z_hit = np.where(arrived, np.inf, z_next)
        thermalised += float(weight[weight > 0.0].sum())

    centres = 0.5 * (bins[:-1] + bins[1:])
    return {
        "aspect_ratio": depth,
        "max_bounces": int(max_bounces),
        "direct_bottom": direct_bottom,
        "cascaded_bottom": cascaded_bottom,
        "total_bottom": direct_bottom + cascaded_bottom,
        "cascaded_share_of_bottom": float(
            cascaded_bottom / max(direct_bottom + cascaded_bottom, 1e-300)),
        "thermalised_share": thermalised,
        "bounce_generations": int(generations),
        # Wall deposition indexed by HEIGHT ABOVE THE ETCH FRONT (bin 0 is at
        # the front, the last bin at the entrance), normalised to unit entering
        # flux; the energy profile carries the same weighting times the beam
        # energy.
        "wall_rate_profile": wall_rate.tolist(),
        "wall_energy_profile": wall_energy.tolist(),
        "wall_height_bin_centres": centres.tolist(),
        "wall_total": float(wall_rate.sum()),
        "closure_residual": float(
            direct_bottom + cascaded_bottom + thermalised
            + float(wall_rate.sum()) - 1.0),
    }


# --- series 4: ARDE of the total energetic delivery ---------------------------

#: Aspect ratios for the ARDE curve.  The low end matters: the trench arc found
#: that the narrow published beam plus the cascade produces *anti*-ARDE (floor
#: recession rising with depth) over AR 0-4, which is unphysical, so the sign of
#: this trend is a gate on whether a beam model may be quoted at all.
ARDE_ASPECT_RATIOS = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 100.0, 150.0, 200.0)


def arde_sweep(iadf, energy_eV, aspect_ratios=ARDE_ASPECT_RATIOS, **kwargs):
    """Total energetic delivery to the etch front vs aspect ratio.

    The hole's etch-front area is aspect-ratio independent, so delivery per
    unit front area is proportional to the delivered *fraction*: a monotone
    decline is ARDE, any rise is anti-ARDE and disqualifies the configuration.
    """
    rows = []
    for aspect in aspect_ratios:
        record = cascade_delivery(aspect, iadf, energy_eV, **kwargs)
        rows.append({
            "aspect_ratio": float(aspect),
            "direct_bottom": record["direct_bottom"],
            "cascaded_bottom": record["cascaded_bottom"],
            "total_bottom": record["total_bottom"],
            "cascaded_share_of_bottom": record["cascaded_share_of_bottom"],
        })
    totals = np.array([row["total_bottom"] for row in rows])
    increments = np.diff(totals)
    return {
        "rows": rows,
        "monotone_declining": bool(np.all(increments <= 1e-12)),
        "max_rise": float(increments.max()) if increments.size else 0.0,
        "total_ratio_first_to_last": float(totals[-1] / totals[0]),
        "arde_verdict": ("ARDE" if np.all(increments <= 1e-12)
                         else "ANTI-ARDE (not usable)"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/curated/hole_study")
    parser.add_argument("--energy-eV", type=float, default=REFERENCE_ENERGY_EV)
    parser.add_argument("--quick", action="store_true",
                        help="coarse quadrature for a smoke pass")
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    workdir = out / "stl"
    workdir.mkdir(exist_ok=True)

    aspects = ASPECT_RATIOS if not args.quick else (50.0, 200.0)
    fractions = TAIL_FRACTIONS if not args.quick else (0.0, 0.65)
    quad = dict(n_polar=48, n_azimuth=24, n_radial=12) if args.quick else {}

    print("series 1: transport reference", flush=True)
    series1 = transport_reference(aspects, workdir)
    for row in series1:
        print(f"  AR {row['aspect_ratio']:6.0f}  exact {row['clausing_exact_algebra']:.6f}"
              f"  santeler {row['santeler_closed_form']:.6f}"
              f"  rel {row['exact_vs_santeler_relative']:.2e}"
              f"  dev {row['relative_deviation']:.2e} -> {row['routed_to']}",
              flush=True)

    print("series 1b: general body-of-revolution operator certification", flush=True)
    probes = []
    for band, az, gen in ((0.20, 32, 8), (0.10, 32, 8), (0.10, 48, 16),
                          (0.05, 48, 16)):
        probe = general_operator_probe(10.0, band, az, gen)
        probe.update({"aspect_ratio": 10.0, "band_length": band,
                      "azimuth_order": az, "generator_order": gen})
        probes.append(probe)
        print(f"  band {band:.2f} az {az} gen {gen}: "
              f"{'CERTIFIED ' + format(probe['value'], '.6f') if probe['certified'] else 'refused, residual ' + format(probe['residual'], '.2e')}",
              flush=True)

    print("series 2: cone acceptance", flush=True)
    series2 = cone_acceptance(aspects, fractions, args.energy_eV)
    for row in series2:
        if row["aspect_ratio"] == max(aspects):
            print(f"  tail {row['tail_fraction']:.2f}  AR {row['aspect_ratio']:.0f}"
                  f"  cone {row['cone_acceptance']:.4f}"
                  f"  sidewall {row['sidewall_share']:.4f}", flush=True)

    print("series 3: cascade", flush=True)
    series3 = []
    for fraction in fractions:
        iadf = kim_2025_reference_iadf(tail_fraction=fraction)
        for aspect in aspects:
            for cap in (PRODUCTION_BOUNCE_CAP, 64):
                record = cascade_delivery(
                    aspect, iadf, args.energy_eV, max_bounces=cap, **quad)
                record["tail_fraction"] = float(fraction)
                series3.append(record)
                print(f"  tail {fraction:.2f} AR {aspect:6.0f} cap {cap:3d}"
                      f"  direct {record['direct_bottom']:.5f}"
                      f"  cascaded {record['cascaded_bottom']:.5f}"
                      f"  total {record['total_bottom']:.5f}"
                      f"  gens {record['bounce_generations']:2d}"
                      f"  closure {record['closure_residual']:+.2e}", flush=True)

    print("series 4: ARDE of total energetic delivery", flush=True)
    series4 = {}
    for fraction in fractions:
        iadf = kim_2025_reference_iadf(tail_fraction=fraction)
        sweep = arde_sweep(iadf, args.energy_eV, **quad)
        series4[f"tail_{fraction:.2f}"] = sweep
        trend = "  ".join(
            f"{row['aspect_ratio']:.0f}:{row['total_bottom']:.4f}"
            for row in sweep["rows"])
        print(f"  tail {fraction:.2f}  {sweep['arde_verdict']}"
              f"  last/first {sweep['total_ratio_first_to_last']:.4f}"
              f"  max_rise {sweep['max_rise']:+.2e}", flush=True)
        print(f"    {trend}", flush=True)

    payload = {
        "reference_energy_eV": float(args.energy_eV),
        "radius": RADIUS,
        "series1_transport_reference": series1,
        "series1b_general_operator_certification": probes,
        "series2_cone_acceptance": series2,
        "series3_cascade": series3,
        "series4_arde": series4,
    }
    (out / "phase1.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote {out / 'phase1.json'}")


if __name__ == "__main__":
    main()
