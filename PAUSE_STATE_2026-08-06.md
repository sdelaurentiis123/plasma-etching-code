# Pause-state handoff — 2026-08-06

Boxes: all campaign instances destroyed (only the user's own `exp015n-spectate`
remains — NOT ours, an external spawner the user should investigate). Suite at
last commit: 1202 passed, 1 skipped. Everything below is committed and pushed
to origin/codex/validation-first-multiphysics.

## Where the science stands

- **Model verified against direct measurement.** Per-ion yields match the
  Karahashi single-species beam experiment to 4.7% (chemical channel) / 1.27x
  (physical) at 1 keV; angular bound and sqrt-E scaling independently
  corroborated. Gates pin this (tests/test_rate_gap_supply_bound.py).
- **Krueger validation (VALIDATION_DOSSIER + BENCHMARK_CERTIFICATION):** mask
  exact (850.2/850) at all six conditions; power-transfer ratios PASS at
  endpoint (r4/6=0.909, r8/6=1.036); necking exact; mouth passes through 45nm
  at t=44s then over-closes to 39.8 (band edge); depth 347 vs 825 = MISS with
  the proof that his own published inputs cannot reach his own figure (blanket
  rate demands ~2.5-3 units/ion vs measured ceiling 1.5 => his Table-I ion flux
  understated ~2x). Karahashi 1keV grade supersedes the 350eV two-channel
  decomposition (reconciliation note in RESULTS_RATE_GAP_CLOSURE).
- **IN FLIGHT WHEN PAUSED (relaunch to finish): the blanket-anchored
  prediction** — calibrate boundary ion flux on his PUBLISHED BLANKET rate
  (independent observable, yields stay measured), then the six-condition board
  as out-of-sample tests. Fork was killed after step 1 of the base run; boxes
  destroyed; the implementation may be partially in the tree/commits of that
  fork — check git log for its commits before redoing. Expectation: depth
  ~2x up toward band, closure ratio drops toward Krueger's, clog/ordering may
  co-recover. This is the "predict depth honestly" move.

## Deliverables shelf (all committed)

- Hole study phases 1+2 (first coupled 200:1 rate-ARDE anywhere; tail changes
  ARDE sign; E8 immaterial-negative). Phase-3 blocker: analytic self-pair
  kernel (tapered profiles), then mask layer in the axisym driver.
- LER demonstrator (chain closed, exact null control; static-shadowing 1:1
  negative with mechanism => dynamic mask-erosion rung preregistered).
- Validation dossier + benchmark certification (the skeptic-facing records).
- Field map (RESEARCH_EXTREME_AR_FIELD): 200:1 = undemonstrated whitespace;
  ranked 10-gap list; our transport matches Lam's published numbers.
- Literature library: research_sources/LIBRARY.md + library/ (92 sources, 844
  claim rows, reverse constant->source index, quarantine table) + 28 full-text
  extracts. Conventions: fetch=>extract+entry same commit; bibkey provenance.

## Open decisions for the user

1. Relaunch the blanket-anchored board (the depth answer) — ~$3, an evening.
2. Send the Krueger email (drafted in RESEARCH_F_SUPPLY_BAND §6) — may yield
   his true fluxes.
3. Rebuild the public HTML artifact from the certification (old page is stale
   both directions).
4. Sean's reactor model: reply drafted (in session); integrate as boundary
   provider; his sources -> library; knobs->profile SF6 loop = first joint
   deliverable.
5. Kill the exp*-spectate spawner on the vast account (not ours, ~$0.31/hr).

## Standing rules (unchanged)

Zero undeclared knobs; look-it-up-first with verbatim quotes (memory:
agents-must-look-it-up); full-text only, no abstract numbers; forecast before
box spend; make_box_archive.sh for ALL deploys (two-layer partner-content
scan); one physics change per run; destroy own boxes; never touch
exp*-spectate; [VERIFY] stays declared without a source.
