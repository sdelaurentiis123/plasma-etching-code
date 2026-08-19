# Pause-state handoff — 2026-08-06

> Historical record. For the current build and scientific state, read
> `CURRENT_STATE.md` and `CONTINUATION_STATE_2026-08-18.md` first.

Boxes: all campaign instances destroyed (only the user's own `exp015n-spectate`
remains — NOT ours, an external spawner the user should investigate). Suite at
this handoff: 1240 passed, 1 skipped. Everything below is committed and pushed
to origin/codex/validation-first-multiphysics.

## Where the science stands — correction

The direct-beam and depth-impossibility claims in the original pause handoff
were retracted after a full Figure-4 pixel audit. The authoritative correction
is `RESULTS_DEPTH_IDENTIFIABILITY_2026-08-06.md`.

- **Species-agnostic validation was false.** The default end-to-end mechanism
  returned the same `0.380584 SiO2/ion` for F+, CF+, CF2+, and CF3+ at 1000 eV
  because energetic ion identity was discarded. Karahashi instead measures a
  strong species ladder (`0.3232`, `0.6751`, `1.1957`, `1.4703`) and CF3+
  reaches `1.8736` at 1500 eV. The former “4.7% independent validation” gate
  was removed.
- **The endpoint comparison is numerically unchanged.** Mask thickness,
  endpoint power ratios, neck, mouth history, and the simulated `346.833 nm`
  versus measured `825 nm` MISS are unchanged.
- **The impossibility proof was false.** `2.521 SiO2/wafer-ion` is a
  run-average lower-bound normalization, not a yield sustained at the floor.
  Takada's independent stable-C5F8/Ar+ beam experiment reports `2.4–2.5` at
  900 eV and ratio 1 (a retained journal/conference source discrepancy),
  proving that Karahashi's rounded `1.5` pure-CF3+ value is not a universal
  surface ceiling.
- **Absolute depth is underidentified.** Krüger publishes an aggregate ion
  flux/IEAD but no positive-ion species composition and no stable C4F6 wafer
  flux. Both missing variables have order-one effects in direct surface
  measurements. Exact 825 nm prediction is not authorized from the published
  boundary.
- **A direct C4F6 diagnostic now strengthens that verdict.** Kim et al. 2021
  measured an undissociated C4F6 mass signal and a multi-species CFx+/CxFy+
  ion distribution at another CCP's powered electrode. Its count rates cannot
  calibrate Krüger, but it confirms the omitted boundary classes are physical,
  not speculative.
- **The code now makes the boundary explicit.** An opt-in Karahashi
  species-resolved table reproduces the measured ladder only inside measured
  energy/angle support and leaves the default aggregate path bitwise
  unchanged. Takada Figure 3 is separately archived, digitized, visually
  audited, and labeled as a C5F8 analog—not a C4F6 law.

## Deliverables shelf (all committed)

- Hole study phases 1+2 (first coupled 200:1 rate-ARDE anywhere; tail changes
  ARDE sign; E8 immaterial-negative). Phase-3 blocker: analytic self-pair
  kernel (tapered profiles), then mask layer in the axisym driver.
- LER demonstrator (chain closed, exact null control; static-shadowing 1:1
  negative with mechanism => dynamic mask-erosion rung preregistered).
- Validation dossier + benchmark certification (the skeptic-facing records).
- Field map (RESEARCH_EXTREME_AR_FIELD): 200:1 = undemonstrated whitespace;
  ranked 10-gap list; our transport matches Lam's published numbers.
- Literature library: research_sources/LIBRARY.md + library/ (94 sources, 858
  claim rows, reverse constant->source index, quarantine table) + 32 full-text
  extracts. Conventions: fetch=>extract+entry same commit; bibkey provenance.

## Open decisions for the user

1. Do **not** relaunch the former blanket-anchored board as a prediction. There
   is no published blanket datum; any ion normalization is a fit to 825 nm.
   First obtain or predict species-resolved ions and stable C4F6 at the wafer.
2. Send the Krueger email (drafted in RESEARCH_F_SUPPLY_BAND §6) — may yield
   his true fluxes.
3. Rebuild the public HTML artifact from the certification (old page is stale
   both directions).
4. Sean's reactor model: integrate as a boundary provider and grade its
   species-resolved outputs. For this depth case it must be extended from SF6
   to C4F6 and validated against ion-composition and stable-parent diagnostics.
5. Kill the exp*-spectate spawner on the vast account (not ours, ~$0.31/hr).

## Standing rules (unchanged)

Zero undeclared knobs; look-it-up-first with verbatim quotes (memory:
agents-must-look-it-up); full-text only, no abstract numbers; forecast before
box spend; make_box_archive.sh for ALL deploys (two-layer partner-content
scan); one physics change per run; destroy own boxes; never touch
exp*-spectate; [VERIFY] stays declared without a source.
