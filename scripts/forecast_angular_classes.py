"""Forecast the depth-rate change from the Appendix-B angular classes.

The three previously-unimplemented angular markers (SiO2(s)+Ar+ class 1,
SiO2CF(s)+Ar+ class 2, AC(s)+Ar+ class 1) multiply the oxide and mask ion
kernels.  Both classes are unity at normal incidence, so the question is what
the *distribution* of incidence angles over the etch front does to the
depth-setting channels.

This script measures that distribution from an archived checkpoint's actual
etch-front geometry, folds in the ion angular distribution the boundary
delivers, and reports the flux-weighted change in the two oxide kernels --
i.e. the predicted depth-rate change -- before any box run is spent.

Usage:  python scripts/forecast_angular_classes.py [checkpoint.npz]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "src")

from petch.mixed_layer import (  # noqa: E402
    _angular_chemical_sputter,
    _angular_physical_sputter,
)

DEFAULT_CKPT = (
    "results/curated/mixed_layer_feature_v1/ml18-depxl-12s/checkpoint.npz")

# Krueger Fig-4 digitised beam, planar sigma (repo digitisation band
# [0.822, 0.860] deg); the axisymmetric lift makes the polar rms sqrt(2) wider
# (RESULTS_ANGULAR_CONVERGENCE_P0_2026-08-02, P1a).
BEAM_SIGMA_PLANAR_DEG = 0.8334


def etch_front_normals(path: str):
    """Surface-normal tilt from vertical, per cell, over the etch front.

    phi > 0 is solid.  The interface normal is grad(phi); its angle from the
    -z axis is the local surface tilt, which equals the ion incidence angle
    for a perfectly vertical beam.
    """
    data = np.load(path)
    phi = np.asarray(data["phi"], dtype=float)
    gz, gy, gx = np.gradient(phi)  # index order (z, y, x) per the engine layout
    mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    band = (np.abs(phi) < 1.0) & (mag > 1e-12)
    if not band.any():
        raise SystemExit("no interface cells found")
    # cosine between the surface normal and the vertical beam axis
    cos_tilt = np.abs(gz[band]) / mag[band]
    return np.clip(cos_tilt, 0.0, 1.0), phi, band


def flux_weighted_factor(cos_face, sigma_deg, factor_fn, n_beam=4001):
    """Average an angular class over the beam spread on each face.

    Each face of tilt alpha sees ions arriving within +-a few sigma of vertical,
    so the incidence angle is alpha convolved with the beam.  Returns the
    per-face flux-weighted factor and the areal-cosine weight of that face.
    """
    alpha = np.degrees(np.arccos(cos_face))
    # beam sample (planar signed angle), Gaussian weights
    span = 4.0 * sigma_deg
    beam = np.linspace(-span, span, n_beam)
    wts = np.exp(-0.5 * (beam / sigma_deg) ** 2)
    wts /= wts.sum()
    inc = np.abs(alpha[:, None] + beam[None, :])
    inc = np.clip(inc, 0.0, 90.0)
    cos_inc = np.cos(np.radians(inc))
    fac = factor_fn(cos_inc)
    # areal flux onto the face scales with cos(incidence)
    areal = cos_inc
    num = (fac * areal * wts[None, :]).sum(axis=1)
    den = (areal * wts[None, :]).sum(axis=1)
    return num / np.maximum(den, 1e-300), den


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CKPT
    cos_face, phi, band = etch_front_normals(path)

    # Restrict to the etch front: cells in the lower half of the feature, i.e.
    # the horizontal-ish surfaces that set depth (tilt < 45 deg from vertical
    # beam means a near-horizontal floor).
    floor = cos_face > np.cos(np.radians(45.0))
    wall = ~floor

    print(f"checkpoint          : {path}")
    print(f"interface cells     : {cos_face.size} "
          f"(floor-like {floor.sum()}, wall-like {wall.sum()})")
    print(f"beam sigma (planar) : {BEAM_SIGMA_PLANAR_DEG} deg")
    print()

    rows = []
    for label, mask in (("floor (depth-setting)", floor),
                        ("wall (lateral)", wall),
                        ("all interface", np.ones_like(floor, dtype=bool))):
        if not mask.any():
            continue
        c = cos_face[mask]
        f1, w = flux_weighted_factor(c, BEAM_SIGMA_PLANAR_DEG,
                                     _angular_physical_sputter)
        f2, _ = flux_weighted_factor(c, BEAM_SIGMA_PLANAR_DEG,
                                     _angular_chemical_sputter)
        # weight faces by the ion flux they intercept
        f1_bar = float((f1 * w).sum() / w.sum())
        f2_bar = float((f2 * w).sum() / w.sum())
        tilt = np.degrees(np.arccos(c))
        rows.append((label, float(np.median(tilt)), f1_bar, f2_bar))

    print(f"{'region':<24}{'median tilt':>12}{'class1 (bare)':>16}"
          f"{'class2 (complex)':>18}")
    for label, tilt, f1_bar, f2_bar in rows:
        print(f"{label:<24}{tilt:>11.2f}d{f1_bar:>16.4f}{f2_bar:>18.4f}")
    print()
    print("Both classes are 1.000 before this change, so the numbers above are")
    print("the multiplicative change in each kernel.")
    print()

    floor_row = next(r for r in rows if r[0].startswith("floor"))
    _, _, f1_floor, f2_floor = floor_row
    # SiO2 removal = bare (class 1) + complex (class 2).  At the Krueger base
    # the complex channel dominates the floor; report both bounds and the
    # complex-dominated estimate that sets depth.
    print(f"predicted depth-rate change (complex-dominated floor): "
          f"{100.0 * (f2_floor - 1.0):+.1f}%")
    print(f"predicted depth-rate change (bare-dominated floor):    "
          f"{100.0 * (f1_floor - 1.0):+.1f}%")
    # Measured split of SiO2 removal between the two channels at the Krueger
    # base, film-free floor (complex 1.652e20 vs bare 1.671e20 -> 49.7 % / 50.3 %).
    complex_share = 0.497
    combined = complex_share * f2_floor + (1.0 - complex_share) * f1_floor
    print(f"predicted depth-rate change (measured 49.7/50.3 split):"
          f" {100.0 * (combined - 1.0):+.1f}%")
    print()
    for base, name in ((825.0, "experimental target"),):
        lo, hi = 0.95 * base, 1.05 * base
        print(f"{name}: {base:.0f} nm, gate band [{lo:.0f}, {hi:.0f}]")


if __name__ == "__main__":
    main()
