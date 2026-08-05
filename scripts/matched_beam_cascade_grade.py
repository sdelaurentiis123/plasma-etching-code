#!/usr/bin/env python3
"""First-principles audit of the reflection cascade's aspect-ratio behaviour.

WHY NOT A DIRECT GRADE AGAINST HUANG'S NUMBERS
----------------------------------------------
petch produces anti-ARDE on trenches (074068e: floor recession RISES 27% from
AR 0 to 4; total energetic delivery falls 3% to AR 16) where Huang's MCFPM
produces the physical trend (etch rate -80% by AR 40).  The obvious move is to
grade our cascade against his published fluxes.  That move fails on a bound
that can be derived before running anything, so it is documented here rather
than attempted:

  For ANY fixed angular distribution and straight walls, the direct-ion
  acceptance of a slot of aspect ratio A with entry position uniform across the
  opening is F(A) = E[max(0, 1 - A|tan b|)].  For a wide beam this is
  asymptotically c/A, so the ratio between two aspect ratios is bounded below
  by their inverse ratio:  F(53)/F(13) >= 13/53 = 0.245.

  Huang's published ion decay (huang_thesis.txt L5405-5407, oxide AR 0 -> 40
  under a PR of AR 13, i.e. total AR 13 -> 53) is
  0.3e15/2.0e15 = 0.150 -- BELOW that bound.

  A numerical bisection over beam width confirms it: the ratio rails at 0.238
  at sigma = 30 deg and never reaches 0.150.

His decay is therefore not fixed-beam straight-wall shadowing, and his own text
says why -- the profile evolves ("the etch front evolving from flat to a
tapered profile", L5443-5446), the PR erodes (L5586-5596: "As the mask erodes,
the conductance limit of neutrals ... relaxes"), the feature is a via rather
than our trench, and every etch-front strike including re-arrivals increments
his flux count (L5399-5402).  Grading our straight-wall frozen-geometry scan
against those numbers would be comparing different observables.  This is the
same conclusion the funnelling pass reached (3489dbe) and it is restated here
with a derivation rather than an observation.

WHAT REPLACES IT: THE GOVERNING LAW, DERIVED
--------------------------------------------
In a straight trench with specular walls, a reflection off a vertical wall
flips the lateral velocity and preserves the axial one, so the polar angle b is
CONSERVED along the cascade.  Eq. 2.34 (L2336-2341, E_ts=100 eV, E_c=10 eV,
theta_c=70 deg) then returns full energy for every bounce with b < 20 deg, so
energy is conserved too.  The only attrition is the leftover-rule reaction
probability per bounce, and with

    bounces to traverse aspect ratio AR:  n = AR * tan(b)
    per-bounce loss:                      r = p0 * kress(90 deg - b)
    small angle:                          kress -> (1 + B) * sin(b)

the survival of a ray to the etch front is

    S(b, AR) = (1 - r)^n  ->  exp( -p0 (1 + B) * AR * b^2 )        [b in rad]

with p0 = 0.9 and B = 9.3, i.e. exp(-9.3 * AR * b^2).  Verified against the
exact product form below: 0.3% at b = 0.5 deg, 1.0% at b = 1 deg (the
small-angle form degrades above ~2 deg, where kress is no longer linear).

The consequence is the entire result: **cascade attenuation is QUADRATIC in
beam width**.  AR enters linearly, b enters squared.  A beam twice as wide
attenuates four times as fast per unit aspect ratio.  Nothing about this is
fitted and it needs no reference simulation to state.

So anti-ARDE is not a cascade defect -- it is the correct behaviour of a
specular cascade under a beam narrow enough that AR*sigma^2 stays small.  The
Krueger IEAD (planar sigma 0.833 deg = 0.0145 rad) gives 9.3*40*0.0145^2 =
0.078 at AR 40, i.e. 92% survival: essentially lossless funnelling.  The
physical trend requires angular width, and the measured collisional tail is
where that width comes from.
"""
import argparse
import json
from pathlib import Path

import numpy as np

# numpy>=2 renamed trapz; the repo carries the same compat alias.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---- Eq. 2.34 constants, verbatim (huang_thesis.txt L2337-2341) -------------
E_TS, E_C, THETA_C = 100.0, 10.0, 70.0
# Production film-sputter law (Krueger Appendix B polymer row, angle class 1).
KRESS_B, REACT_P0 = 9.3, 0.9

# Beam widths under test (planar sigma, degrees).
KRUEGER_PLANAR_SIGMA_DEG = 0.8334      # P1a-lifted digitised IEAD
# Kim 2025 measured two-component IADF (iadf_two_component.py provenance):
# core T_perp 0.044 eV, tail T_perp 0.57 eV -> width ratio sqrt(0.57/0.044).
KIM_CORE_SIGMA_DEG = 0.380
KIM_TAIL_WIDTH_RATIO = float(np.sqrt(0.57 / 0.044))

# Huang Fig. 6.6(a)/6.7 anchors (L5405-5449, L5410-5417) -- context recorded,
# NOT used as a grading target (see module docstring).
HUANG = {
    "ion_ar0": 2.0e15, "ion_ar40": 0.3e15,
    "hot_ar0": 3.1e15, "hot_ar4": 8.0e15, "hot_ar40": 1.1e15,
    "pr_ar": 13.0,
    "pr_sweep_ion_ar0": 4.2e15, "pr_sweep_ion_ar20": 1.3e15,
}


def kress(cos_t):
    return np.maximum((1.0 + KRESS_B * (1.0 - cos_t ** 2)) * cos_t, 0.0)


def retained_energy(E, theta_deg):
    """Eq. 2.34 exactly as boundary_transport_3d.py:782-790 applies it."""
    E = np.asarray(E, float)
    th = np.asarray(theta_deg, float)
    interp = E * ((th - THETA_C) / (90.0 - THETA_C)) * ((E - E_C) / (E_TS - E_C))
    return np.maximum(np.where((E > E_TS) & (th > THETA_C), E,
                               np.where((E < E_C) | (th < THETA_C), 0.0, interp)),
                      0.0)


def survival_exact(beta_deg, aspect_ratio, energy_eV=1500.0):
    """Exact bounce-product survival of one ray (specular, straight trench)."""
    beta_deg = float(abs(beta_deg))
    if beta_deg <= 0.0:
        return 1.0, energy_eV
    n = aspect_ratio * np.tan(np.radians(beta_deg))
    theta_wall = 90.0 - beta_deg
    cos_wall = float(np.cos(np.radians(theta_wall)))
    r = float(np.clip(REACT_P0 * kress(cos_wall), 0.0, 1.0))
    E = float(energy_eV)
    if theta_wall > THETA_C and E > E_TS:
        # Eq. 2.34 returns the incident energy unchanged in this regime, so the
        # bounce loop is a fixed point: energy is conserved for every bounce.
        return (1.0 - r) ** n, E
    whole = int(np.floor(n))
    for _ in range(min(whole, 4096)):
        E = float(retained_energy(E, theta_wall))
        if E <= E_C:
            return 0.0, 0.0
    return (1.0 - r) ** n, E


def survival_small_angle(beta_deg, aspect_ratio):
    b = np.radians(np.abs(beta_deg))
    return float(np.exp(-REACT_P0 * (1.0 + KRESS_B) * aspect_ratio * b ** 2))


def beam_survival(sigma_deg, aspect_ratio, *, tail_fraction=0.0,
                  tail_width_ratio=1.0, n_nodes=4001, max_sigma=8.0):
    """Flux-weighted survival of a (optionally two-component) Gaussian beam.

    Returns the fraction of entering energetic flux still alive at the etch
    front, i.e. the cascade's aspect-ratio transfer function.
    """
    widths = [(1.0 - tail_fraction, sigma_deg)]
    if tail_fraction > 0.0:
        widths.append((tail_fraction, sigma_deg * tail_width_ratio))
    total = 0.0
    for weight, sig in widths:
        if weight <= 0.0:
            continue
        grid = np.linspace(-max_sigma * sig, max_sigma * sig, n_nodes)
        pdf = np.exp(-0.5 * (grid / sig) ** 2)
        pdf /= _trapz(pdf, grid)
        surv = np.array([survival_exact(b, aspect_ratio)[0] for b in grid])
        total += weight * float(_trapz(pdf * surv, grid))
    return total


def acceptance_fraction(aspect_ratio, sigma_deg, n_rays, seed):
    """Direct-hit fraction: F(A) = E[max(0, 1 - A|tan b|)] by sampling."""
    rng = np.random.default_rng(seed)
    beta = rng.normal(0.0, sigma_deg, n_rays)
    x0 = rng.uniform(0.0, 1.0, n_rays)
    avail = np.where(beta >= 0.0, 1.0 - x0, x0)
    need = np.abs(np.tan(np.radians(beta))) * aspect_ratio
    return float(np.mean(need <= avail))


def bound_demonstration(n_rays=200_000, seed=11):
    """Show Huang's ion ratio lies below the fixed-beam straight-wall bound."""
    lo, hi = HUANG["pr_ar"], HUANG["pr_ar"] + 40.0
    rows = []
    for sigma in (0.8334, 2.0, 5.0, 10.0, 30.0):
        f_lo = acceptance_fraction(lo, sigma, n_rays, seed)
        f_hi = acceptance_fraction(hi, sigma, n_rays, seed)
        rows.append({"sigma_deg": sigma, "F_ar13": f_lo, "F_ar53": f_hi,
                     "ratio": (f_hi / f_lo) if f_lo > 0 else None})
    return {
        "huang_published_ratio": HUANG["ion_ar40"] / HUANG["ion_ar0"],
        "wide_beam_bound": lo / hi,
        "sweep": rows,
        "verdict": "Huang's 0.150 lies below the 13/53 = 0.245 wide-beam bound; "
                   "his ion decay is not fixed-beam straight-wall shadowing "
                   "(tapering + PR erosion + via geometry + re-arrival counting).",
    }


def run(args):
    ars = [1.0, 2.0, 4.0, 8.0, 16.0, 40.0, 100.0, 200.0]
    out = {
        "law": {
            "form": "S(b,AR) = (1-p0*kress(90-b))^(AR*tan b) -> "
                    "exp(-p0*(1+B)*AR*b^2), b in radians",
            "p0": REACT_P0, "B": KRESS_B,
            "derivation": "specular wall reflection conserves polar angle and "
                          "(via Eq. 2.34, b<20deg) energy; only the leftover-"
                          "rule reaction attenuates",
            "validation": [],
        },
        "bound_demonstration": bound_demonstration(),
        "beam_transfer": {},
    }
    for beta in (0.5, 1.0, 2.0, 5.0):
        for ar in (10.0, 40.0):
            exact = survival_exact(beta, ar)[0]
            approx = survival_small_angle(beta, ar)
            out["law"]["validation"].append({
                "beta_deg": beta, "aspect_ratio": ar, "exact": exact,
                "small_angle": approx,
                "rel_error": abs(approx - exact) / max(exact, 1e-12)})

    beams = {
        "krueger_digitised_narrow": dict(sigma_deg=KRUEGER_PLANAR_SIGMA_DEG,
                                         tail_fraction=0.0),
        "kim_core_only": dict(sigma_deg=KIM_CORE_SIGMA_DEG, tail_fraction=0.0),
        "kim_two_component_f0.35": dict(sigma_deg=KIM_CORE_SIGMA_DEG,
                                        tail_fraction=0.35,
                                        tail_width_ratio=KIM_TAIL_WIDTH_RATIO),
        "kim_two_component_f0.50": dict(sigma_deg=KIM_CORE_SIGMA_DEG,
                                        tail_fraction=0.50,
                                        tail_width_ratio=KIM_TAIL_WIDTH_RATIO),
        "kim_two_component_f0.65": dict(sigma_deg=KIM_CORE_SIGMA_DEG,
                                        tail_fraction=0.65,
                                        tail_width_ratio=KIM_TAIL_WIDTH_RATIO),
    }
    for label, cfg in beams.items():
        rows = []
        for ar in ars:
            rows.append({"aspect_ratio": ar,
                         "cascade_survival": beam_survival(aspect_ratio=ar, **cfg)})
        base = rows[0]["cascade_survival"]
        out["beam_transfer"][label] = {
            "config": cfg, "rows": rows,
            "survival_ar40_over_ar1": rows[5]["cascade_survival"] / base,
            "survival_ar200_over_ar1": rows[-1]["cascade_survival"] / base,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=Path("results/curated/cascade_matched_beam/grade.json"))
    args = ap.parse_args()
    out = run(args)

    print("LAW VALIDATION  S = exp(-9.3*AR*b^2)")
    print(f"  {'beta':>6} {'AR':>6} {'exact':>10} {'small-angle':>12} {'rel err':>9}")
    for r in out["law"]["validation"]:
        print(f"  {r['beta_deg']:6.1f} {r['aspect_ratio']:6.0f} {r['exact']:10.5f}"
              f" {r['small_angle']:12.5f} {r['rel_error']:8.1%}")

    bd = out["bound_demonstration"]
    print(f"\nBOUND: Huang ion ratio {bd['huang_published_ratio']:.3f} vs "
          f"wide-beam bound {bd['wide_beam_bound']:.3f}")
    for r in bd["sweep"]:
        print(f"  sigma {r['sigma_deg']:6.2f} deg  F(13)={r['F_ar13']:.4f} "
              f"F(53)={r['F_ar53']:.4f}  ratio={r['ratio']:.3f}")

    print("\nCASCADE TRANSFER (fraction of entering energetic flux alive at front)")
    hdr = f"  {'beam':<30}" + "".join(f"{a:>9.0f}" for a in
                                      [1, 2, 4, 8, 16, 40, 100, 200])
    print(hdr.replace("       1", "     AR1"))
    for label, blk in out["beam_transfer"].items():
        row = "".join(f"{r['cascade_survival']:9.4f}" for r in blk["rows"])
        print(f"  {label:<30}{row}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
