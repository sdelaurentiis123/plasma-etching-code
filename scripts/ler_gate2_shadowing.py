"""LER Gate 2 — the Constantoudis ion-shadowing transfer rule, graded exactly.

Published rule (RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29 sec. 2.2, from
Constantoudis, Kokkoris & Gogolides, JM3 12(4), 041310 (2013) and the authors'
SPIE Newsroom summary):

  * substrate sigma_S vs resist sigma_R is linear with slope ~ 0.5 above a
    threshold;
  * threshold sigma_R* ~ xi_R / (c tan theta_R) with c ~ 2.0-2.5;
  * equivalently LER reduction requires (sigma_R/xi_R) tan theta_R > 1/c.

Their stated mechanism is pure ion shadowing by the rough, tapered resist
sidewall -- the substrate edge inherits the envelope of the resist edge. That is
the same geometry our exchange operator computes exactly, so the rule is
gradeable with a static calculation over a synthesised rough-edge ensemble: no
fitting, no digitisation, no feature evolution.

Geometry (all lengths nm, angles deg):

  * resist edge displacement u(y) at the mask top, synthesised Palasantzas;
  * tapered sidewall at theta_R from horizontal, so the wall surface sits at
    x_wall(y, z) = u(y) + (h - z) cot(theta_R) for 0 <= z <= h -- the foot
    protrudes into the open space by h cot(theta_R);
  * ions arrive with transverse tangents (tx, ty) ~ N(0, s^2), s = tan(sigma_th),
    the same tangent-plane Gaussian the two-component IADF uses per component;
  * a ray leaving substrate point (x, y) clears the mask iff for every height z
        x + z tx >= u(y + z ty) + (h - z) cot(theta_R),
    so the blocking value is B(y; tx, ty) = max_z [u(y + z ty) + (h - z) cot -
    z tx] and the received flux fraction at (x, y) is the weighted CDF of B.
    The etched edge is the level-quantile of B -- exact, no bisection.

Vertical ions (s -> 0) give B = u(y) + h cot, i.e. sigma_S = sigma_R and no
transfer at all; the y-coupling that produces the published rule requires the
angular spread, and the threshold is where the roughness slope sigma_R/xi_R
beats the wall taper cot(theta_R) as seen along a tilted ray. That makes c a
measure of the effective angular spread, which this script measures.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from petch.ler_metrology import EdgeStatistics, synthesize_edge_nm  # noqa: E402


def _gauss_hermite_tangents(sigma_theta_deg: float, n_nodes: int):
    """Probabilists' Gauss-Hermite nodes/weights for N(0, s^2) tangents."""
    nodes, weights = np.polynomial.hermite_e.hermegauss(int(n_nodes))
    weights = weights / np.sqrt(2.0 * np.pi)
    scale = np.tan(np.deg2rad(float(sigma_theta_deg)))
    return nodes * scale, weights


def _interp_periodic(u_nm, y_query_nm, spacing_nm):
    """Linear interpolation of the periodic edge at arbitrary y."""
    n = u_nm.size
    position = y_query_nm / spacing_nm
    lower = np.floor(position).astype(np.int64)
    frac = position - lower
    lower = np.mod(lower, n)
    upper = np.mod(lower + 1, n)
    return u_nm[lower] * (1.0 - frac) + u_nm[upper] * frac


def shadowed_substrate_edge_nm(u_nm, *, spacing_nm: float, mask_height_nm: float,
                               sidewall_angle_deg: float,
                               sigma_theta_deg: float, n_nodes: int = 11,
                               n_height: int = 25, level: float = 0.5):
    """Etched-edge displacement from exact ray shadowing of a rough taper.

    Returns the edge at the requested flux level (0.5 = half of open field).
    """
    u = np.asarray(u_nm, dtype=float)
    cot = 1.0 / np.tan(np.deg2rad(float(sidewall_angle_deg)))
    tangents, weights = _gauss_hermite_tangents(sigma_theta_deg, n_nodes)
    tx = tangents[:, None]
    ty = tangents[None, :]
    weight = (weights[:, None] * weights[None, :]).ravel()
    heights = np.linspace(0.0, float(mask_height_nm), int(n_height))

    y = np.arange(u.size, dtype=float) * spacing_nm
    # (Ny, Ndir, Nz) blocking surface; max over z gives the binding constraint.
    y_samples = (y[:, None, None]
                 + heights[None, None, :] * ty.ravel()[None, :, None])
    wall = _interp_periodic(u, y_samples, spacing_nm)
    wall = wall + (mask_height_nm - heights)[None, None, :] * cot
    wall = wall - heights[None, None, :] * tx.ravel()[None, :, None]
    blocking = wall.max(axis=2)                                   # (Ny, Ndir)

    order = np.argsort(blocking, axis=1)
    sorted_blocking = np.take_along_axis(blocking, order, axis=1)
    sorted_weight = weight[order]
    cumulative = np.cumsum(sorted_weight, axis=1)
    cumulative /= cumulative[:, -1][:, None]
    index = (cumulative >= float(level)).argmax(axis=1)
    edge = sorted_blocking[np.arange(u.size), index]
    return edge - edge.mean()


def transfer_sweep(sigma_values_nm, *, correlation_length_nm: float,
                   roughness_exponent: float, sigma_theta_deg: float,
                   sidewall_angle_deg: float, mask_height_nm: float,
                   n_points: int, spacing_nm: float, seeds):
    """Ensemble sigma_in -> sigma_out over the declared roughness sweep."""
    out = []
    for sigma_in in sigma_values_nm:
        realized_in, realized_out = [], []
        for seed in seeds:
            statistics = EdgeStatistics(float(sigma_in),
                                        float(correlation_length_nm),
                                        float(roughness_exponent))
            edge = synthesize_edge_nm(statistics, n_points=n_points,
                                      spacing_nm=spacing_nm, seed=int(seed))
            etched = shadowed_substrate_edge_nm(
                edge, spacing_nm=spacing_nm, mask_height_nm=mask_height_nm,
                sidewall_angle_deg=sidewall_angle_deg,
                sigma_theta_deg=sigma_theta_deg)
            realized_in.append(float(np.std(edge)))
            realized_out.append(float(np.std(etched)))
        out.append({
            "sigma_declared_nm": float(sigma_in),
            "sigma_in_nm": float(np.mean(realized_in)),
            "sigma_out_nm": float(np.mean(realized_out)),
            "ratio": float(np.mean(realized_out) / np.mean(realized_in)),
        })
    return out


def fit_slope_and_threshold(sweep, *, ratio_cut: float = 0.9):
    """Slope of the reducing branch and its crossing with sigma_out = sigma_in.

    The published form is sigma_out = m sigma_in + b above a threshold, with
    sigma_out = sigma_in below it; the threshold is where the two lines meet.
    """
    sigma_in = np.array([row["sigma_in_nm"] for row in sweep])
    sigma_out = np.array([row["sigma_out_nm"] for row in sweep])
    ratio = sigma_out / sigma_in
    overall_slope = float(np.polyfit(sigma_in, sigma_out, 1)[0])
    reducing = ratio < ratio_cut
    if reducing.sum() < 2:
        # No reducing branch: the transfer is a rigid copy over the whole sweep.
        return {"slope": float("nan"), "intercept_nm": float("nan"),
                "threshold_nm": float("nan"), "n_reducing": int(reducing.sum()),
                "overall_slope": overall_slope,
                "min_ratio": float(ratio.min())}
    slope, intercept = np.polyfit(sigma_in[reducing], sigma_out[reducing], 1)
    threshold = intercept / (1.0 - slope) if slope < 1.0 else float("nan")
    return {"slope": float(slope), "intercept_nm": float(intercept),
            "threshold_nm": float(threshold), "n_reducing": int(reducing.sum()),
            "overall_slope": overall_slope, "min_ratio": float(ratio.min())}


def measured_c(threshold_nm, *, correlation_length_nm, sidewall_angle_deg):
    """c := xi / (sigma* tan theta_R) -- the published rule solved for c."""
    if not np.isfinite(threshold_nm) or threshold_nm <= 0.0:
        return float("nan")
    return float(correlation_length_nm
                 / (threshold_nm * np.tan(np.deg2rad(sidewall_angle_deg))))


def main():
    # Constantoudis Fig. 4 conditions: xi = 30 nm, alpha = 0.6, resist 150 nm,
    # sidewall 86.2 deg (angle from horizontal; the convention check in the
    # research doc rejects the from-vertical reading as unphysical).
    base = dict(correlation_length_nm=30.0, roughness_exponent=0.6,
                sidewall_angle_deg=86.2, mask_height_nm=150.0,
                n_points=1024, spacing_nm=1.0, seeds=(11, 12, 13, 14))
    sigma_values = (0.15, 0.3, 0.6, 1.2, 2.5, 5.0, 9.0)

    report = {"conditions": {k: v for k, v in base.items() if k != "seeds"},
              "seeds": list(base["seeds"]),
              "sigma_declared_nm": list(sigma_values),
              "beam_sweep": [], "xi_scaling": [], "taper_scaling": []}

    for sigma_theta in (2.0, 10.0, 25.0, 45.0):
        sweep = transfer_sweep(sigma_values, sigma_theta_deg=sigma_theta, **base)
        fit = fit_slope_and_threshold(sweep)
        report["beam_sweep"].append({
            "sigma_theta_deg": sigma_theta, "sweep": sweep, "fit": fit,
            "c": measured_c(fit["threshold_nm"],
                            correlation_length_nm=base["correlation_length_nm"],
                            sidewall_angle_deg=base["sidewall_angle_deg"]),
        })
        print(f"beam sigma_theta={sigma_theta:5.1f} deg  slope={fit['slope']:.3f}"
              f"  threshold={fit['threshold_nm']:.3f} nm"
              f"  c={report['beam_sweep'][-1]['c']:.2f}", flush=True)

    # Structural claim 1: threshold proportional to xi at fixed taper and beam.
    for xi in (15.0, 30.0, 60.0):
        conditions = dict(base, correlation_length_nm=xi)
        sweep = transfer_sweep(sigma_values, sigma_theta_deg=25.0, **conditions)
        fit = fit_slope_and_threshold(sweep)
        report["xi_scaling"].append({
            "correlation_length_nm": xi, "fit": fit,
            "c": measured_c(fit["threshold_nm"], correlation_length_nm=xi,
                            sidewall_angle_deg=base["sidewall_angle_deg"]),
            "sweep": sweep})
        print(f"xi={xi:5.1f} nm  threshold={fit['threshold_nm']:.3f} nm"
              f"  slope={fit['slope']:.3f}", flush=True)

    # Structural claim 2: threshold proportional to cot(theta_R).
    for angle in (84.0, 86.2, 88.0):
        conditions = dict(base, sidewall_angle_deg=angle)
        sweep = transfer_sweep(sigma_values, sigma_theta_deg=25.0, **conditions)
        fit = fit_slope_and_threshold(sweep)
        report["taper_scaling"].append({
            "sidewall_angle_deg": angle, "fit": fit,
            "c": measured_c(fit["threshold_nm"],
                            correlation_length_nm=base["correlation_length_nm"],
                            sidewall_angle_deg=angle),
            "sweep": sweep})
        print(f"theta_R={angle:5.1f} deg  threshold={fit['threshold_nm']:.3f} nm"
              f"  slope={fit['slope']:.3f}", flush=True)

    destination = Path("results/curated/ler_demonstrator")
    destination.mkdir(parents=True, exist_ok=True)
    with open(destination / "gate2_shadowing.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {destination / 'gate2_shadowing.json'}")


if __name__ == "__main__":
    main()
