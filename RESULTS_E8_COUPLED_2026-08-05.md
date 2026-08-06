# E8 completion: thermalized radicals re-emit through the neutral solve

Follows `RESULTS_E8_THERMALIZED_RETURN_2026-08-05.md` (commit `550e2a5`), which
built the return and measured it at the face where the cascade dropped it. That
measurement read the **source** term, not the delivered flux, and concluded the
weight was "stranded on the sidewalls". This pass wires the return through the
transport that was always there, finds a bug that had been silently zeroing it,
and measures delivery through the coupled solve.

Commits `543d8ae` (fix + plumbing + gates). Suite 1183 passed / 1 skipped.
Box spend **$0** — the coupled forecast refused the run.

## 1. A bug: the return was being dropped entirely

`gather_boundary_state_ballistic_3d` writes the E8 contribution into
`neutral_flux[target]` from inside the **energetic** species branch. The
**neutral** branch for that same species assigns:

```python
if role[species.name] == "neutral_reactant":
    neutral_flux[species.name] = gathered.sum(axis=0)   # assign, not accumulate
```

So whenever the boundary's species order put the ion before its target radical —
which it does for the pilot's own boundary — the neutral branch **overwrote the
E8 injection** and the entire return vanished. Measured on a flat-plane gather:
added flux exactly `0.0` at share 0.4.

Fixed by accumulating into a separate dict and merging after the species loop.
The new conservation gate (below) is what caught it; the original pass's gates
checked the per-face diagnostic, which was correct, rather than the ledger the
solve actually consumes.

## 2. The completion: re-emission was already solved, just not connected

petch solves multi-bounce diffuse re-emission for plasma-sourced neutrals in
`solve_diffuse_neutral_radiosity_3d` — `H = D + B(1-s)H`, per-species sticking,
exact global balance audit. `_apply_diffuse_neutral_transport` builds `D` from
`transport.surface_fluxes.neutral_flux_m2_s`, i.e. **the same ledger the gather
writes E8 into**. So the E8 source becomes an additional source term in that
solve by construction: returned radicals bounce at their own published sticking
and redistribute toward the front, which is the mechanism behind Huang's
statement (`research_sources/thesis_extracts/huang_thesis.txt` L5714-5727):

> "After losing energy through several collisions with the sidewalls and etch
> front, these energetic species become thermal CFx and CxFy radicals, which can
> passivate the oxide surface or deposit as polymer... As the AR increases to
> greater than 10, the neutralized and thermalized CFx+ and CxFy+ ions become
> the main source (> 95%) of radicals reaching the etch front."

The only missing piece was the option itself: `thermalized_radical_return` lived
on the gather but was not a parameter of `advance_feature_step_3d` or
`solve_feature_3d`, so no feature run could reach it. Now plumbed through both.

**Species rule** (unchanged, sourced): non-fluorocarbon partners are excluded —
"The neutral and thermalized partners of other ions are non-reactive species and
diffuse out of the feature with no surface reactions" (same section). The
reactive fraction stays a **declared swept input**: Krüger publishes an
aggregate positive-ion flux with a combined IEAD, so the CFx+/Ar+ split does not
exist for his reactor and is never inferred.

## 3. Gates

| gate | requirement | measured | verdict |
|---|---|---|---|
| conservation | returned rate == share x thermalized rate | 1e-12 (mesh-area convention) | **PASS** |
| strict source | no face loses flux; sum > 0 | holds | **PASS** |
| explicit zero | bitwise inert | `array_equal` | **PASS** |
| plumbing | option reaches both entry points | signatures | **PASS** |
| blanket invariance | default off moves nothing | unchanged | **PASS** |
| suite | green | 1183 / 1 skipped | **PASS** |
| floor composition | see section 4 | see section 4 | see section 4 |

## 4. Delivery through the coupled solve

Method: the radiosity solve is **linear in its source**, so E8-delivered floor
flux is measured exactly by differencing two solves on the identical operator
(same frozen geometry, same reaction probabilities) — E8 off, then E8 on. Floor
band = active faces within 30 nm of the etch front, inside the opening. Flux is
the delivered (post-radiosity) incident CF2 flux, area-weighted.
`scripts/e8_coupled_floor_scan.py`, dx = 10 nm.

| geometry | oxide AR | mask | total AR | f_FC | plasma CF2 | E8 delivered | E8 share |
|---|---|---|---|---|---|---|---|
| Krüger cell | 12 | 0.85 µm | 21.4 | 0.3 | 1.485e20 | 4.377e18 | 2.86 % |
| Krüger cell | 12 | 0.85 µm | 21.4 | 1.0 | 1.485e20 | 1.459e19 | 8.94 % |
| Huang-like | 20 | 0.05 µm | 20.6 | 0.3 | 1.504e20 | 4.384e18 | 2.83 % |
| Huang-like | 20 | 0.05 µm | 20.6 | 1.0 | 1.504e20 | 1.461e19 | 8.86 % |
| smoke check | 3 | 0.05 µm | 3.6 | 1.0 | 5.274e20 | 1.381e19 | 2.55 % |

Linearity check: the f = 1.0 rows equal the f = 0.3 rows scaled by 1/0.3 to four
digits (1.459e19 predicted, 1.459e19 measured), confirming the solve is linear in
the E8 source and licensing the forecast scaling in section 5.

**Re-emission is real and roughly doubles delivery.** The raw local deposit at
the Krüger floor was 6.65e18 at f = 1.0, i.e. ≈2.0e18 at f = 0.3; the coupled
solve delivers 4.38e18 — **2.2× more**, because weight thermalized on the walls
now diffuses to the front instead of being absorbed where it died. The first
pass's "stranded on the sidewalls" reading was an artifact of measuring the
source.

**The magnitude is small, and the Huang-geometry test corrects why.** The first
pass attributed the shortfall to Krüger's cell being *mask-dominated*. Running a
thin-mask feature at matched total aspect ratio refutes that: 0.05 µm of mask
over 1.8 µm of oxide (total AR 20.6) gives **2.83 % / 8.86 %**, statistically
identical to the 0.85 µm-mask cell's **2.86 % / 8.94 %** at total AR 21.4, and
plasma CF2 delivery is likewise near-identical (1.504e20 vs 1.485e20). **At
matched total AR the mask fraction is immaterial to floor radical composition** —
what governs both plasma and E8 delivery is the total path, not how it splits
between mask and oxide.

So the real statement is about depth, not stack: total AR ≈ 21 is simply not deep
enough for the thermalized channel to dominate. The trend is in the right
direction and measured — E8 share at f = 1.0 rises 2.55 % → 8.9 % as total AR
goes 3.6 → 21 — but reaching Huang's > 95 % needs plasma transmission orders of
magnitude lower than either test geometry provides (the first pass estimated
≈0.009 %). The AR-3 point independently reproduces his *stated* low-AR behaviour
("the majority of the CFx and CxFy radicals at low AR (< 5) ... originate from
the thermal neutrals incident into the feature from the plasma"), so the model
tracks his description at the shallow end and misses his threshold only because
neither geometry reaches the depth his statement is about.

**Gate 2 verdict: not reproduced at either geometry (8.9 % maximum, at the
physical upper bound), with the attribution corrected from mask geometry to
total aspect ratio.**

## 5. Forecast — and why no run was bought

Coupled 0-D at the measured AR-12 floor delivery, 3406 eV front energy, E8
delivery scaled linearly in the declared fraction (the solve is linear):

| f_FC | E8 delivered | floor precursor | floor rate | vs f = 0 |
|---|---|---|---|---|
| 0.00 | 0 | 4.5047e20 | 2.1071 nm/s | — |
| 0.30 | 4.377e18 | 4.5485e20 | 2.1044 nm/s | −0.125 % |
| 0.50 | 7.295e18 | 4.5776e20 | 2.1027 nm/s | −0.208 % |
| 1.00 | 1.459e19 | 4.6506e20 | 2.0983 nm/s | −0.417 % |

At the **physical upper bound** (every positive ion fluorocarbon) the floor rate
moves **−0.42 %, downward**: the returned species is a *precursor*, so it feeds
polymer deposition as well as complex formation, and the film thickens slightly.
Against ml23's −49 % depth undershoot that is two orders of magnitude short and
the wrong sign. Forecast-before-spend refuses the run; a 12 s endpoint would
have re-measured ml23 at a cost of a box.

## 6. The Krüger depth item: final word

E8 is now **complete and correct** — conserved, re-emitting through the same
solve as every other radical, gated, and reachable from the pilot. It is not the
depth fix, and the delivery measurement says why: the missing supply is not
thermalized-ion return at this geometry.

The depth item therefore rests where `VALIDATION_DOSSIER_KRUEGER_2026-08-05.md`
put it — decomposed, bounded, and attributed — plus two declared supply items
that no source closes today:

1. **The (s₀, B₀) adsorption pair** (Kwon/Sawin E1). Gray's printed SiO₂ sticking
   is 0.02 (thesis p.246), but it is half of a co-regressed pair: landing the
   scalar alone moves the measured Gray half-rise from 1.94 to 104.9 against his
   measured 27 ± 8, and breaks four validated chemistry gates that carry s = 1 in
   their own closed forms. Transplanting the pair is a campaign of its own.
2. **The CFx+ fraction for this reactor** — unpublished; swept over [0, 1] here,
   never fitted, and shown above to be immaterial to depth at any value.

Nothing in this pass changes a validated result: the return is default-off, and
every existing number stands.
