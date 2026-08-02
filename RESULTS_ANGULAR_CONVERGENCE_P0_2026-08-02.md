# P0 — angular-convergence harness results (2026-08-02)

Executes S0 of `RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md` / P0 of
`ROADMAP_PERFECTION_2026-07-29.md`: **build the ruler before building physics.**

Harness: `scripts/angular_convergence_harness.py` (plots:
`scripts/plot_angular_convergence.py`). Static ideal trench geometry, the
production deterministic face gather, **no chemistry, no surface evolution, no
box** — every number below is pure transport. Data and figures:
`results/curated/angular_convergence_p0/`.

Mesh: straight trench, width 1, depth = AR, floor banded into 200 cells across
the opening, each sidewall into 60 depth bands (120 for the AR-9 mouth study),
mask top on both sides, periodic laterally (trench array), source plane one
width above the mask. Face count is independent of AR, so the ablation isolates
*angular* resolution from spatial resolution.

## Gate — analytic view factor (passed)

A single collimated direction at polar angle θ into a straight trench must
illuminate a floor fraction `1 - AR·tan θ`.

| θ | observed | analytic | abs error |
|---|---|---|---|
| 0.5° | 0.737671 | 0.738194 | 5.2e-04 |
| 1.0° | 0.475214 | 0.476348 | 1.1e-03 |
| 1.5° | 0.214959 | 0.214422 | 5.4e-04 |

All within the 1/200 = 5e-3 floor discretization. (The gate *earned its keep*:
the first harness build used a single wide floor face, which resolves shadowing
at two centroids and gave errors up to 0.21 — a spatial artifact that would have
contaminated every angular conclusion.)

---

## EXP A — what angular quadrature buys

### A1. Gauss–Hermite transverse order (virtual-sheath beam, σ = 0.148°)

| AR | GH 3 | GH 5 | GH 9 | GH 17 | 9→17 change |
|---|---|---|---|---|---|
| 30 | 0.0416 | 0.0480 | 0.0523 | 0.0547 | +4.6% |
| 100 | 0.1387 | 0.1593 | 0.1734 | 0.1816 | +4.7% |
| 200 | 0.2513 | 0.2984 | 0.3295 | 0.3438 | +4.3% |

(wall flux / mouth flux)

**σ_θ is exactly invariant at 0.1478° for every order**, while the transported
wall flux moves **+37% at AR 200** across the sweep and is **still not converged
at order 17**. This is the research doc's pathology #2 confirmed
quantitatively: Gauss–Hermite converges *moments*, and deep features are
decided by the *support*. A moment-converging quadrature is the wrong
discretization for HAR, and no order of it is "converged" in the sense that
matters.

### A2. Azimuthal order (digitized-IEAD axisymmetric lift)

| AR | az 4 | az 8 | az 16 | az 32 | az 64 | 16→64 change |
|---|---|---|---|---|---|---|
| 30 | 0.2420 | 0.2224 | 0.2181 | 0.2170 | 0.2168 | −0.6% |
| 100 | 0.6231 | 0.5615 | 0.5399 | 0.5372 | 0.5363 | −0.7% |
| 200 | 0.8085 | 0.7483 | 0.7164 | 0.7073 | 0.7056 | −1.5% |

**The production azimuthal order (16) is converged to ~1.5% even at AR 200.**
This *refutes* the natural reading that the standing "1.5× azimuth deficit" is a
node-count problem — more azimuthal nodes buy essentially nothing. See EXP C for
what the deficit actually is.

### A3. Polar bin width (digitized IEAD)

| AR | exact digitized | 0.25° | 0.50° | 1.00° |
|---|---|---|---|---|
| 30 | 0.2181 | 0.2181 | 0.2181 | 0.2182 |
| 100 | 0.5399 | 0.5399 | 0.5445 | 0.5662 |
| 200 | 0.7164 | 0.7164 | 0.7372 | 0.8256 |

Two findings. (1) Binning at 0.25° is **bitwise identical** to the exact
digitized measure — confirming the source grid *is* 0.25°, so the production bin
is free, and **refinement below 0.25° is impossible on this data at any cost**.
(2) Coarsening to 1.0° moves the AR-200 floor delivery by −38%, so angular bin
width is a first-order control at AR ≥ 100 — the 0.25° data grid is the hard
ceiling on AR-200 fidelity from this boundary, exactly as the research doc
warned.

---

## EXP B — the mouth question at AR 9 (the ml13 geometry)

Standing hypothesis: the ml13 mouth-opening gap (24.8 vs 45 nm) is a lost
wide-angle tail. Four beams, matched quadrature, same geometry:

| beam | planar σ | floor frac | wall frac | mouth-band ratio vs narrow |
|---|---|---|---|---|
| collisionless core only (petch virtual sheath) | 0.148° | 0.9814 | 0.0186 | 1.00 |
| measured core + collision tail (Kim 2025) | 0.428° | 0.9510 | 0.0490 | 2.63 |
| **Krüger digitized IEAD, production lift (the ml13 path)** | 0.589° | 0.9341 | 0.0658 | **3.53** |
| Krüger digitized IEAD, closure-corrected | 0.833° | 0.9069 | 0.0930 | **4.99** |

**Finding 1 — the premise was wrong, and this is the important correction.**
The ml13 feature validation does **not** run the narrow collisionless beam. It
consumes the digitized Krüger IEAD directly
(`build_krueger_2024_development_boundary` → `load_krueger_2024_digitized_iead`),
which **already contains the collisional tail**. So the mouth gap cannot be
explained by "the tail is missing from the boundary" — at AR 9 the ml13 path
already delivered 3.5× the mouth-region wall flux of a tail-free beam. The S0
question "is the tail missing, or does transport lose it?" resolves to: **the
tail is present; transport loses part of it.** How much is EXP C.

**Finding 2.** The measured bi-Gaussian (built from Kim 2025's core 0.044 eV /
tail 0.57 eV at 3465 eV, tail fraction 0.65) sits *between* the narrow beam and
the digitized IEAD. The digitized HPEM IEAD is wider than the measured
distribution — consistent with the research doc's open note that the digitized
σ_tail (0.946°) exceeds the measured one (0.520°), plausibly log-colorbar
digitization over-weighting the far tail.

---

## EXP C — the azimuthal closure error (the actual deficit)

The production lift takes the published **planar signed angle** of the Krüger
IEAD and uses it **directly as the polar angle** of an axisymmetric 3-D beam
(`development_species`, `azimuthal_closure="axisymmetric_uniform"`). For any
axisymmetric beam with Gaussian transverse components, the polar-angle rms
exceeds the planar-marginal rms by exactly √2. Measured on the actual data:

```
published signed-angle (planar) sigma = 0.8334 deg
lifted planar sigma (production)      = 0.5893 deg
ratio                                 = 1.4141   (sqrt2 = 1.41421)
```

**The production lift discards exactly √2 of the published angular width.** This
is an analytic, exactly-reproducible closure error, not an empirical fudge.

Transport consequence (production lift vs the same IEAD with polar angles scaled
by √2 so the lifted planar marginal reproduces the published width):

| AR | wall frac (production) | wall frac (corrected) | wall ratio | floor ratio |
|---|---|---|---|---|
| **9** | 0.0658 | 0.0930 | **×1.414** | ×0.971 |
| 30 | 0.2181 | 0.2994 | ×1.373 | ×0.896 |
| 100 | 0.5399 | 0.6324 | ×1.171 | ×0.798 |
| 200 | 0.7164 | 0.7871 | ×1.099 | ×0.749 |

**This is the headline result.** At AR 9 — the ml13 geometry — correcting the
closure multiplies sidewall flux by **1.414**, and the standing mouth-gap
hypothesis called for ≈1.5×. The long-quoted "1.5× azimuth deficit" is, to
within its own uncertainty, **this √2 closure error**, now identified in closed
form and reproducible in one line.

The error is *not* AR-independent: on the walls it decays with AR (geometry
saturates — at AR 200 nearly everything already hits a wall), while on the
**floor** it grows monotonically (×0.971 → ×0.749). At AR 200 the closure error
alone inflates bottom delivery by **33%** (1/0.749).

---

## Verdict for the roadmap

1. **P0's hypothesis is resolved, and differently than expected.** The ml13
   mouth gap is not a missing physical tail — the tail was in the boundary all
   along. It is an **azimuthal-closure error of exactly √2** in lifting a
   published planar IEAD to 3-D. Cheap to state, cheap to test, and it lands
   precisely on the magnitude the mouth residual needed.
2. **Azimuthal node count is exonerated** (converged to 1.5% at 16 nodes) —
   spending effort there would have been wasted.
3. **Gauss–Hermite is the wrong representation for HAR** and is not converged at
   any tested order; the S1 replacement is justified by measurement, not taste.
4. **0.25° is a hard data ceiling** on the Krüger boundary. AR 200 needs 0.057°,
   which this source physically cannot supply — S1's gate must anchor on the
   Kim 2025 measured IADF, as the research doc already argued.

### Recommended P1 spec (revised by these measurements)

- **P1a (do first, hours not days): fix the azimuthal closure.** Infer the polar
  density from the published planar marginal instead of identifying them —
  exactly the √2 scaling for Gaussian transverse components, and an Abel-type
  inversion in general. Gate: the lifted beam's planar marginal reproduces the
  published σ within the repo's own digitization band (0.822–0.860°). Then
  re-run the ml13 base case: preregistered expectation is the mouth opens
  toward 45 nm without touching a single chemistry constant. This is the
  cheapest open shot at the last standing feature-scale residual.
- **P1b: the two-component `IonAngularEnergyDistribution`** as specified in S1
  (explicit (E, θ, φ), core + tail, analytic erf acceptance, T⊥ derived from
  T_gas), gated against Kim 2025 rather than HPEM output.
- **P1c: retire Gauss–Hermite** for the ion transverse representation; A1 shows
  it cannot converge the observable that matters.
- **Angular AMR (P5/S3) is confirmed cheap**: at AR 200 azimuth is already
  converged at 16 nodes, so refinement effort belongs entirely in the polar
  variable inside the acceptance cone — which is what S3 proposed.

---

## P1a landed (2026-08-02)

`src/petch/angular_lift.py` replaces the identity closure with an **onion-peel
Abel inversion** of the published planar marginal, wired into
`Krueger2024DigitizedIEAD.development_species` and receipted in species
provenance (`three_dimensional_polar_inversion*`). Gates:
`tests/test_angular_lift.py` (8).

| quantity | before | after |
|---|---|---|
| lifted planar rms (production path) | 0.5893° | **0.8233°** (band 0.822–0.860) |
| published / lifted deficit | ×1.4141 | **×1.0124** |
| E[tan²θ] / E[tan²θ_x] | 1.000 | **2.000 (exact, 0 ulp)** |

Three results worth recording beyond the fix itself.

1. **The √2 is not a Gaussian result.** Squaring `tan θ_x = tan θ · cos φ` and
   averaging over uniform azimuth gives `E[tan²θ] = E[tan²θ]/2` for *any*
   axisymmetric measure. The P0 text called it "for Gaussian transverse
   components"; it is in fact shape-independent, which is why the digitized
   (non-Gaussian) data reproduced 1.4141 so precisely. The discrete lift
   honours it to 0 ulp because a uniform azimuthal ring averages cos² to ½
   exactly for order ≥ 3 — so the gate is machine-precision, not statistical.
2. **The inversion is exact but ill-posed, and 0.25° is not a convention —
   it is the resolution the data carries.** Onion peeling is exact
   back-substitution (round-trip reproduces the planar histogram to 2.8e-17,
   no regularization parameter). But peeling *finer* than the source grid
   diverges into negative shells: clamped mass fraction is 0.00 at 0.25°,
   0.28 at 0.20°, and 2.3 (i.e. divergence) at 0.10°. The lift now **refuses**
   above a 5% clamp guard rather than returning a manufactured answer. This
   independently confirms A3's finding that 0.25° is the hard data ceiling —
   from the inversion side rather than the transport side.
3. **The correction is ×1.397 in width, not the full ×1.414.** The
   EXP C diagnostic scaled all polar angles by √2 uniformly; the true
   inversion redistributes weight shape-dependently and lands at 0.8233 rather
   than 0.8334 (1.2% low, inside the digitization band). Expected AR-9 sidewall
   gain is therefore ≈×1.40, marginally below the diagnostic's ×1.414.

**Preregistered ml16.** Pilot flags (add `--resume` for restart):

```
python scripts/krueger_2024_trench_pilot.py \
  --dx-um 0.01 --radiosity-backend deterministic_extruded_2d \
  --transport-device cuda:0 --surface-state-remap-backend common_refinement \
  --topology-change-policy continue_gas_cavity --surface-model mixed_layer \
  --max-wall-s 86400 --duration-s 60 \
  --mixed-layer-volatilization-yield 1.0 \
  --grazing-ion-reflection literature_v1 --output <OUT>
```

**Scope note the runner must resolve before launching.** These flags at current
HEAD do *not* reproduce ml13's chemistry: the closure batch (`82c223b`,
`23fe6de`) moved the constants to the Table-6.5 converged set with the
activated-SiO₂ and de-crosslink channels — that is the **ml15** configuration.
So as written this is "ml15 constants + corrected lift". Running "ml13
constants + corrected lift" (the literal preregistration) requires reverting
those two commits' constant changes first. Both are defensible experiments and
the choice is a physics decision, not a flag: ml13-constants isolates the lift
against the config of record, ml15-constants tests whether the corrected wall
flux was the missing ingredient that made the verbatim-complete set
underperform. Recommend ml15-constants first, since the P0 finding predicts the
lift supplies exactly the mouth-region flux whose absence ml15 was blamed on.

### Declared limitations of this pass

- Ideal straight-wall geometry; no evolving profile, no chemistry, no
  redeposition. These are transport ratios, not predicted CDs.
- The √2 correction implemented in EXP C is a **diagnostic scaling**; the
  principled inversion landed in P1a above and does differ in the far tail
  (0.8233° vs the diagnostic's 0.8334°).
- P1a lifts the angular **marginal** with a per-bin weight multiplier, so the
  inversion is uniform in energy. The digitization is a point cloud with 1–4
  angle nodes per energy row, which cannot support a per-energy peel; an
  energy-resolved lift needs the S1 two-component distribution, not this data.
- The AR-9 "mouth band" is defined as the upper quarter of the sidewall; the
  ratio is insensitive to that choice at the reported precision but the
  definition is arbitrary.
- Floor discretization (200 bands) sets a ~5e-3 floor on all fractions.
