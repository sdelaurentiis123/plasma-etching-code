# The muted O-radical channel: a per-cell probability applied per atom (2026-08-04)

`RESULTS_WALL_SLOPE_FALSIFICATION_2026-08-04.md` (committed `3b7f104`) localised
the remaining trench defect: our closure rate matches Krüger's to within ±12 %
below 270 nm and diverges monotonically toward the surface, reaching **10.8× at
the mask top**.  That depth signature excludes every depth-uniform suspect and
points at a channel whose delivery is maximal at the unshadowed top.  The
surviving candidate was the O-radical polymer etch — the lever Krüger names
verbatim as controlling necking, and one that behaved as if it were inert in our
runs: raising `p_ox` by 48 % (ml16a → ml16b) moved the neck by only 11 %.

This pass measured the channel and found the cause.

## The target is exact, not approximate

Krüger's two relevant rows are both **per-cell**: a deposition event adds one
polymer cell, and `O(g) + P(s) → products` at `p_ox = 0.0423` removes one.  Their
ratio at his Table-I fluxes is therefore a geometry-free property of the
mechanism — both are thermal fluxes, so view factors cancel:

```
p_ox * J_O                       0.0423 * 7.70e20     3.2571e19
------------------------  =  ----------------------  =  ---------  =  0.1953
sum_i p_dep,i * J_i          0.1*4.4e20 + 0.1*9.4e20   1.6680e20
                             + 0.1*8.4e19 + 0.03*6.8e20
```

Reproducing 0.1953 to four digits from the published fluxes and probabilities
confirms the reading: **one O collision with polymer removes one polymer unit.**

## What our model did instead

`src/petch/mixed_layer.py` computed the film-oxidation carbon loss as

```python
ox_c = params.oxidation_probability * fluxes.oxygen_flux * theta_film * x_c
f_per_ox_c = _guarded_ratio(state.n_f_film, state.n_c_film, 2.0)
ox_f = ox_c * f_per_ox_c
```

The `f_per_ox_c` line already implements per-cell removal — its own comment says
"each oxidized C carries along the local film F/C ratio".  Multiplying the carbon
by `x_c = n_c/(n_c + n_f)` **as well** double-counts the composition: `p_ox` is a
per-collision-with-polymer probability that already subsumes which atom the
oxygen lands on.  The two factors cancel exactly in the atom count, so the term
removed precisely **one atom per collision** instead of one polymer unit:

```
ox_c * (1 + F/C) = p_ox * J_O * theta * x_c * (1 + n_f/n_c) = p_ox * J_O * theta
```

## Candidate verdicts

| candidate | verdict | evidence |
|---|---|---|
| (a) `theta_film` gating | **not the cause** | the pinch regime is thick film: the ml16a checkpoint's mask faces carry 1.173e21 atoms/m² = 15.6 nm, so `theta_film = 1.000000` to machine zero |
| (b) site competition / availability clamp | **not the cause** | `scale_film` only engages when losses exceed availability; at 15.6 nm of film the reservoir is four orders above the per-step loss |
| (c) O transport delivery to the top | **not the cause** | the O channel is thermal and the top faces are unshadowed; the throttle is present at full source flux in the 0-D budget, which has no transport at all |
| (d) **normalisation vs the published per-collision row** | **CONFIRMED** | `x_c = 0.371 ± 0.004` across every filmed face of the ml16a checkpoint (F/C = 1.690), giving a uniform **1/x_c = 2.69× throttle** wherever film exists |

The composition is remarkably uniform — `x_c ∈ [0.370, 0.381]` over 750 mask
faces — so this is a clean multiplicative throttle, not a state-dependent one.

## The fix

Drop the spurious factor; gate on exposed polymer alone, as the published row
states:

```python
ox_c = params.oxidation_probability * fluxes.oxygen_flux * theta_film
```

No constant was changed, nothing was imported, no mechanism was added.

## Gate numbers

| quantity | before | after | target |
|---|---|---|---|
| O-removal / deposition (per-cell, thick film, Table-I fluxes) | 0.0726 | **0.1953** | 0.1953 (0.02 % off) |
| carbon removed per reactive O collision | `x_c` = 0.371 | **1.0** | 1 (one polymer unit) |
| atoms removed per reactive O collision | 1.00 | **2.69** | 1 + F/C = 2.69 |
| response of removal to +48 % `p_ox` | throttled | **+48 %, exactly proportional** | proportional |

Gates live in `tests/test_o_channel_budget.py` (5 tests).  Element ledgers still
close below 1e-9 at the stronger rate.

## Test repaired, not weakened

`test_oxygen_thins_film_and_moves_clog_boundary` asserted that precursor
1.0e21 clogs at `J_O` = 2e21 and is rescued at 8e21.  With O 2.69× stronger,
2e21 now rescues it — the clog boundary moved up by the same factor in precursor
flux, which is the corrected physics, not a regression.  The probe point was
moved to precursor 2.0e21, where the identical assertion holds (clogs at 2e21,
rescued at 8e21).  The physics asserted is unchanged: a clog boundary exists and
oxygen moves it.

Measured boundary after the fix:

| precursor | `J_O` = 0 | 5e20 | 2e21 | 8e21 | 2e22 |
|---|---|---|---|---|---|
| 1.0e21 | clog | clog | etch | etch | etch |
| 2.0e21 | clog | clog | **clog** | **etch** | etch |

## What this does and does not claim — it closes 15 % of the top-band gap

The channel is now correct, but it is **not sufficient**, and the budget
arithmetic says so before any run does.  Normalising to deposition at the top
band, with ion removal there at the measured 88.7° incidence
(`cos²θ (1 + B sin²θ) = 0.0053`):

| quantity (× deposition) | pre-fix | post-fix |
|---|---|---|
| O removal | 0.073 | **0.195** |
| ion removal at 88.7° | 0.005 | 0.005 |
| net closure | 0.922 | **0.800** |

Krüger's top-band rate is 10.8× smaller than our pre-fix rate, which implies his
net closure is 0.085 × deposition — i.e. **his removal is ~91 % of his
deposition at the mask top, against our ~20 % after this fix.**  So:

- top-band excess: 10.8× → **9.4×** (target ≤ 1.3×)
- **fraction of the gap closed: 15 %**
- still unexplained: **~0.71 × deposition** of removal at the mask top

That residual is far too large for any remaining thermal-channel bookkeeping and
sits exactly where the falsification pass left its open [VERIFY]: the **angular
ion-removal law at near-grazing incidence**.  Our lip faces sit at 86–89°, where
`cos²θ (1 + B sin²θ)` collapses by 190× and delivers 0.5 % of deposition, while a
per-particle voxel code lets grazing ions and their reflected hot neutrals keep
sputtering the lip.  The measured above-cosine data (You 2023) could not be
tested against our neck geometry because those faces span only 78–80°; testing it
needs a geometry spanning 30–70°, which the new `--angle-deg` probe mode can now
build.

The frozen-geometry probe's absolute rates run ~10× below evolution (methodology
caveat in the falsification doc), so the graded box run remains the arbiter of
the closure ratio.  **It should not be spent on this fix alone** — on this
arithmetic it would return ~9× and read as a null result.  The next free step is
the grazing-incidence removal audit, not a run.
