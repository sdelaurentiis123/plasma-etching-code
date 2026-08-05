# Morning report — 2026-08-05

Overnight state of the Krüger base validation.  Everything below is committed
and pushed; no box is left running.

## Where the validation stands

| gate | target | best measured | verdict |
|---|---|---|---|
| mask remaining | 850 nm | **850.2** (ml19) | **PASS**, exact |
| mask aperture | 45 nm (`w_m`) | **50.9** equilibrium (ml19) | near-miss, +13 % |
| constriction depth | 200 (SEM) / 271 (sim) | 170 nm (ml19) | −15 % vs SEM |
| closure/etch, t ≥ 8 s | 0.0310 | 0.0257 / 0.0153 / 0.0160 | **in band** |
| closure/etch, t = 1–4 s | 0.0310 | 0.0569 (ml18) | 1.83×, above band |
| etch depth at 60 s | 825 ± 5 % | ~1066 extrapolated (ml19) | **MISS**, +29 % |

The mouth arc closed at the mechanism level yesterday: ml19 is the first run
in the campaign whose aperture **equilibrates** (50.92 nm mean over
t = 34–46 s, drift −17 pm/s) instead of sealing.  Every earlier configuration
sealed monotonically.

## What was established overnight

**1. The crosslink partner count is published — the last mouth `[VERIFY]` is
closed.**  Krüger et al., *JVST A* **42**, 043008 (2024), sec. III
(`tmp/pdfs/krueger-2024.txt` L388-390):

> "...based on the number of available bonds (three in the example in Fig. 5).
> For example, CF2 would have a maximum of two crosslinks and CF3 would have a
> maximum of a single crosslink."

Two worked examples pin the rule uniquely — `available = 4n − m − 2(n−1)` for
C_nF_m, reproducing CF → 3, CF₂ → 2, CF₃ → 1 — so it is applied to the
published stoichiometry and passed to the layer as the deposition-weighted
mean.  Nothing is fitted.  Forecast effect at the audited lip delivery:
x_xl 0.703 → **0.793**, film growth 0.831 → **0.673 nm/s**, i.e. the lip
closure driver moves **1.95× → 1.58×** Krüger's 0.427 nm/s per side.
Gated by 7 tests including the worked examples themselves.  Commit `914bb8d`.

A composition reading of the same rule was implemented first, then **reverted
and re-derived from the published examples** — the interim version flattered
the result without a tabulated source, which is the definition of fitting.

**2. The depth gate is an energetic-delivery question, not a chemistry
question — and this is the sharpest result of the night.**  Scanning the
delivery plane (`scripts/floor_delivery_scan.py`), oxide recession is

- **exactly proportional to ion flux** (7.54 / 15.10 / 22.78 nm/s at 1× / 2× / 3×), and
- **flat to under 1 % across a 50× range of neutral delivery**.

So the floor is a pure energetic observable at these conditions.  ml19's
measured 21.9 nm/s (AR ≈ 11.7) and 16.7 nm/s (AR ≈ 32) therefore imply
**2.9× and 2.2× the source ion flux** arriving at the floor, against
Krüger's implied **1.8×**.  Our funnelled cascade over-delivers to the floor
by 1.2–1.6×, and that ratio *is* the depth overshoot.

**No chemistry constant can close this gate** without moving the blanket rate
that is not in question — a depth-uniform change scales every delivery
equally.  The next test is the cascade's floor delivery against Huang's
published funnelling curve (his total energetic flux falls *below* incident by
AR 40; ours is still 2.2× at AR 32).

**3. Aspect ratios were being quoted wrong.**  The etch front sits below an
850 nm mask, so ml19 spans AR ≈ 11.7 → 32, not the "1 → 15" used in earlier
notes.  Conclusions are unchanged; the comparison to funnelling data must use
the corrected values.

**4. Chang & Sawin ∠=2 roll-off downgraded from `[VERIFY]`** to a quantified
approximation — the digitized curve tracks our class-2 form within 0.065
absolute.  No code change needed.

## Runs and spend

| run | purpose | result |
|---|---|---|
| ml19 (yesterday) | 60 s endpoint, crosslink fix | mouth equilibrates 50.9; stopped t = 46.2 s on a mask sliver |
| ml20 (overnight) | 12 s confirmation of published partner counts | see below |

Box spend overnight: one RTX 3090 at $0.113/h for the ml20 confirmation.
Preregistered gate for ml20: closure/etch per window against Krüger's 0.0310
and against ml18's 0.0569 / 0.0495 / 0.0367, with the forecast predicting a
~19 % improvement and an aperture at t = 12 s above ml18's 64.89 nm.

## Blockers and next moves, in order

1. **Floor energetic delivery** — the depth gate.  Measure the cascade's
   ion + hot-neutral flux at the floor versus depth from the ml19 checkpoint
   and grade against Huang's funnelling numbers.  Free, receipts-based.
2. **Sub-resolution sliver reabsorption** — 60 s runs at dx = 10 nm stop
   around t = 46 s when a 2-cell fragment detaches from the armoured mask
   (a remeshing artifact; the mask cannot physically fragment).  Needed before
   any endpoint run can complete.  Precedent exists in the interior-gas
   nucleation guard: restore pre-step values for sub-resolution components
   with a receipt.
3. **Scorecard** once the base gates settle — the 8-condition O₂ and power
   sweeps are the real transfer test.

## Standing doctrine observed

Every number committed tonight carries a quoted source line.  One value was
implemented, caught without a tabulation, and reverted before the source was
found — recorded here because the revert is as much a result as the fix.
