# Falsifying the wall-slope verdict before implementing it (2026-08-04)

`RESEARCH_LIP_CERTAINTY_2026-08-04.md` (committed 2d57a86) concluded that the
remaining 15 % lip imbalance is a **wall-slope** effect: ion removal at a near-vertical
wall carries `cos θ` twice, so net velocity is monotone in the tilt `α` from
vertical, the balance zero-crosses at `α ≈ 9-11°`, and our evolution runs pinch
because the simulated shoulder is "too vertical".  The proposed fix was to
initialise the mask shoulder from the digitised SEM.

That verdict was checked before any of it was implemented.  **It does not
survive.**  The measurements below are the reason nothing was changed.

## (a) Our simulated shoulder is not vertical

Wall angle from the local slope of the half-aperture, per depth band, measured
on the archived checkpoints through `scripts/regrade_neck_metrics.py`:

| depth below mask top | ml16a | ml16b | ml13 |
|---|---|---|---|
| 0-50 nm | **17.31°** | **17.69°** | 14.28° |
| 50-100 nm | **10.56°** | 8.45° | 11.99° |
| 100-150 nm | 2.20° | 0.30° | 4.64° |
| 150-200 nm | 0.93° | 3.75° | 5.90° |
| 200-270 nm | 9.29° | 8.93° | 1.72° |

The premise is refuted directly: the top shoulder sits at 10-17° from vertical,
*inside* the claimed 9-11° "open" regime, and closes anyway.  The near-zero
readings at 100-200 nm are where the profile passes its own minimum — and a
smooth minimum has zero slope by definition.  That is the same tautology the
research pass correctly identified in the probe's Gaussian bump, now recognised
in our evolved geometry too: `α → 0` at a pinch is a consequence of pinching,
not its cause.

## (b) Krüger's own neck sits at 1.9° and stays open

Same measurement on the digitised Fig. 7 profiles (`tmp/mouth_profiles/`):

| depth | MCFPM sim | SEM |
|---|---|---|
| 0-50 nm | 7.29° | 44.69° (corner clip artefact) |
| 50-100 nm | 9.59° | 11.02° |
| 100-150 nm | 5.96° | 9.08° |
| 150-200 nm | 3.46° | 6.24° |
| **200-270 nm** | **1.86°** | 20.13° (past the SEM minimum) |
| mean 0-850 nm | 1.21° | 0.81° |

His simulated neck — the converged 38.8 nm minimum at 271 nm — sits on a wall
tilted **1.86°** from vertical, and it holds for the full 60 s.  If net velocity
zero-crossed at 9-11°, that neck would close.  **The pure-angle theory cannot be
the mechanism that holds 45 nm open.**

## (c) Angle is a real modifier, but there is no zero crossing

`scripts/mouth_equilibrium_probe.py --angle-deg` (new) replaces the Gaussian
with a straight wedge of prescribed half-angle, choosing the apex aperture so
that the aperture at the measurement band midpoint is 45 nm for *every* angle —
so angle varies at fixed local width, which the aperture sweep could not do.
Band: 15-55 nm above the apex, apex 250 nm below the mask top, dx = 10 nm.

| prescribed α | realised α | band aperture | net velocity | mean cos θ |
|---|---|---|---|---|
| 1.0° | 1.00° | 45.00 nm | −0.00352 nm/s | 0.0271 |
| 3.0° | 3.00° | 45.00 nm | −0.00243 nm/s | 0.0561 |
| 8.0° | 8.01° | 45.00 nm | −0.00037 nm/s | 0.1400 |
| 12.0° | 11.98° | 45.01 nm | −0.00004 nm/s | 0.2070 |
| 15.0° | 14.90° | 44.99 nm | −0.00001 nm/s | 0.2564 |

Angle dependence is real and strong — 350× over 1→15°, the right order for the
`sin²α (1 + B cos²α)` scaling the research pass derived.  But the sign never
changes: the rate decays asymptotically toward zero and every geometry closes.
Angle *modulates* the imbalance; it does not create an equilibrium, so there is
no `α ≈ 9-11°` crossing to initialise the shoulder onto.

Structural note from the same sweep: substrate etch is **identically zero on
every mask face at every angle** — correct, since Krüger's a-C mask is
sputter-armoured and his `hm` target is the full 850 nm — so the entire lip
balance is carried by the fluorocarbon film, and `net` here is exactly
−(film growth).

## The finding that replaces it: the defect is depth-localised

Mean closure rate per side, inferred from the 90 nm initial opening over 60 s,
ours (ml16a) against his MCFPM:

| depth band | our aperture | his aperture | our rate | his rate | ratio |
|---|---|---|---|---|---|
| 0-50 nm | 50.2 | 86.3 | 0.332 nm/s | 0.031 nm/s | **10.8×** |
| 50-100 nm | 25.1 | 70.1 | 0.541 | 0.166 | **3.26×** |
| 100-150 nm | 13.1 | 57.3 | 0.641 | 0.273 | 2.35× |
| 150-200 nm | 15.0 | 47.9 | 0.625 | 0.351 | 1.78× |
| 200-270 nm | 29.4 | 42.1 | 0.505 | 0.399 | 1.27× |
| 270-400 nm | 53.2 | 48.0 | 0.307 | 0.350 | 0.88× |
| 400-600 nm | 69.8 | 71.8 | 0.168 | 0.152 | 1.11× |

Below ~270 nm we reproduce his closure rate to within ±12 %.  The excess switches
on above 270 nm and rises monotonically to **10.8× at the mask top**.

This signature is itself a strong discriminator.  **Anything depth-uniform is
excluded**: sticking magnitudes, the film sputter law, the per-cell-vs-per-atom
asymmetry (rank 2 of the research pass), the mask AC Kress factor (rank 3) — all
would bias every band alike, and the deep bands are already correct.  What
remains are quantities whose delivery is maximal at the mask top and shadowed
with depth:

1. **The O-radical channel.**  Krüger names it as the controlling lever
   verbatim ("necking and clogging can ultimately be controlled by the reaction
   probability of the O based polymer etching"), and its flux is maximal at the
   unshadowed top.  In our runs it is nearly inert: ml16a → ml16b raises
   `p_ox` 0.0423 → 0.0628 (+48 %) and moves the neck only 11.1 → 12.3 nm (+11 %)
   and top CD 62.3 → 62.7 nm (+0.6 %).  If the O channel carried the ~19.5 % of
   the deposition budget his Table-I fluxes imply, a +48 % change would cut the
   closure rate by roughly two thirds.  It does not, so the channel is throttled
   somewhere between flux and film carbon.
2. **The deposition view factor at the top lip**, which is the other quantity
   that peaks at the surface and decays with shadowing.

## Verdict: nothing implemented

- Option (i) — SEM-initialised shoulder — **rejected**: the mechanism it was
  meant to exploit is falsified by (a) and (b).  (Krüger's initial mask is a
  vertical 90 nm opening, so an SEM-derived initial shoulder would also not be
  faithful to the run being reproduced.)
- Option (ii) — Cho/Schaepkens-bounded angular law — **rejected as the mouth
  fix**, and it points the wrong way: at the lip's incidence (83°), `B = 9.3`
  gives an angular factor of 1.24 while a measured-bounded `B ≈ 1.7` (peak/normal
  = 1.3) gives 0.33.  Adopting the measured bound would *reduce* lip removal ~3.8×
  and close the top faster.  The B = 9.3 over-peaking is still a real
  [VERIFY] item (its source, Kress 1999, is Ar⁺ on Cu(111) at 50-250 eV) but it
  is not this defect, and it must not be changed while it is load-bearing
  elsewhere.
- Option (iii) — per-atom asymmetry — **not implemented**, as directed, and now
  additionally excluded by the depth signature.  The receipt it needs remains a
  blanket-film QCM comparison (Chae, Vitale & Sawin 2003).

## Methodology caveat on the probe (new, and it bounds the earlier run)

The frozen-geometry probe reports absolute rates ~10× below what the evolution
runs require: at the 45 nm geometry its fastest band (0-50 nm) reads
−0.034 nm/s, while ml16a's top band closed at 0.33 nm/s.  On a frozen geometry
the film piles up indefinitely instead of being converted into boundary motion,
so the relaxed state is not the evolution's quasi-steady state.  **Relative**
comparisons at equal relaxation (the angle sweep, the depth profile) remain
valid; the earlier run's *absolute* claims — "removal is 85.5 % of deposition",
"no equilibrium aperture exists" — are bounded by this and should be re-derived
once the relaxation is anchored to a representative film thickness.

## Next decisive test

Decompose the *top-band* budget into deposition, O-removal and ion-removal per
face, and compare O-removal / deposition against the 0.195 implied by his
Table-I fluxes.  A ratio far below 0.195 confirms a throttled O channel — the
single mechanism that matches both the depth signature and Krüger's own
statement about what controls necking.  This needs one probe extension, no box.

## Outcome of that test (2026-08-04)

Run and answered in `RESULTS_O_CHANNEL_2026-08-04.md`.  The O channel *was*
throttled — by 2.69×, from a per-cell probability applied per atom — and the fix
is a one-line removal of a spurious composition factor, with no constant
changed.  But the budget arithmetic that the corrected ratio makes possible also
bounds it: the fix closes **15 %** of the top-band excess (10.8× → 9.4×), and
Krüger's top-band removal must be ~91 % of his deposition against our ~20 %.
The remaining ~0.71 × deposition of unexplained removal sits at the grazing
ion-removal law on the 86-89° lip faces — the [VERIFY] this document already
carried.  The graded box run should wait for that audit.
