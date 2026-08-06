# Validation dossier — Krüger 2024 Ar/C₄F₆/O₂ → SiO₂ trench

**Status: campaign CLOSED 2026-08-05.** Written for a reader who intends to
disbelieve it. Every number below is reproducible from a committed artifact;
every commit hash is on `origin/codex/validation-first-multiphysics`.

Reference: J. Krüger, S. Huang *et al.*, *JVST A* **42**, 043008 (2024) and the
accompanying PhD thesis. Target feature: 90 nm line opening in an 850 nm
amorphous-carbon mask over SiO₂, 60 s etch, Ar/C₄F₆/O₂ CCP.

---

## 1. The claim

petch reproduces this feature **with zero constants fitted to it**. Every
chemistry number is lifted from the source's published mechanism (Appendix B /
Table 6.5) or from an independent measurement, each carrying a provenance
string. No knob was tuned to move a gate. Where a number could not be sourced,
the gate was left missed and the gap declared rather than closed by fitting —
§5 is the complete ledger of those.

The comparison is *hostile to us by construction*: the source's own five tuned
constants were fitted by his optimizer against these very metrics, ours were
not. The honest framing is therefore "how far does a zero-added-knob transcription
of a published mechanism get, and where does it break" — not "we match".

## 2. Final gate table

| gate | target | petch | verdict |
|---|---|---|---|
| mask remaining | 850 nm | **850.2 nm** (ml19) | **PASS — exact** |
| mask constriction (`w_m`) | 45 nm | **50.9 nm, equilibrium** (drift −17 pm/s over 12 s) | **+13 %**, residual bounded (§5.1) |
| aperture at 270 nm depth | 41.1 nm (his sim) | **41.8 nm** (ml16a) | **PASS — 1.7 %** |
| aperture at 350 nm depth | 50.1 nm (his sim) | 57.1 nm | +14 % |
| closure/etch, t ≥ 8 s | 0.0310 | 0.0257 / 0.0153 / 0.0160 | **in band** |
| closure/etch, t = 1–4 s | 0.0310 | 0.0569 (1.83×) | above band, §5.1 |
| clog @ O₂ 0.5 | sealed | **0.000 — sealed** | **PASS — exact** |
| necking absent @ O₂ 2.5 | open | open | **PASS** |
| etch depth, reference-energy conditions | 825 ± 5 % | in band (ml6 811, ml9a 791, ml13 852) | **PASS** |
| etch depth, feature keV conditions | 825 ± 5 % | ~1066 extrapolated (+29 %) | **MISS — attributed, §3** |

Sources: `RESULTS_ML19_ENDPOINT_2026-08-05.md` (`fabc2d8`),
`RESULTS_NECK_REGRADE_2026-08-04.md` (`fcc98eb`),
`RESULTS_ML18_DEPOSITION_CROSSLINK_2026-08-04.md` (`391ee43`),
scorecard 1 (`d6068f3`), campaigns 2–9 in
`MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md`.

**The mouth result is the campaign's headline.** Every run before `3a931b1`
sealed monotonically (11.1 nm, 7.2 nm). ml19 is the first to hold a stable
aperture — the mechanism, not the number, is what changed.

**The neck-band result matters more than the scalar.** At the depth where the
SEM and his simulation actually neck, petch reads 41.8 against his 41.1. The
defect is confined to the top ~150 nm of mask, which the single-scalar
`mask_opening` metric could not distinguish from a size error
(`RESULTS_NECK_REGRADE`).

## 3. The depth miss, attributed

Depth is in band at reference-energy conditions and **+29 % at the multi-keV
energies actually measured at the etch front (3406 eV, mean cos 0.768)**. Four
independent routes converge on one term — the ion-energy channel magnitude,
which petch inherits from Appendix B row-for-row:

1. **The cascade audit** (`fd61bb7`, `RESULTS_CASCADE_FIRST_PRINCIPLES_2026-08-05.md`):
   coupled floor rate measured ~1.4× high, with transport exonerated —
   Eq. 2.34 transcription verified exact against `huang_thesis.txt` L2336-2341,
   and delivery to the front measured flat to −1 % across AR 0–16.
2. **The joint model solve** (`953f01c`/`f322cb4`,
   `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`): the oxide assignment uniquely
   selected by Gray's measured beam curve (peak-normalisation, dynamic range
   0.210 against a measured 0.20–0.30) undershoots the depth requirement by
   **1.24–1.35×**. The assignment the measurements select is the one the depth
   gate rejects.
3. **The source author states it himself.** Krüger thesis L4884-4888, verbatim:

   > "These trends indicate that above a certain threshold energy the etch
   > progression and the mask removal process are not ion starved, but rather
   > limited by neutral gas transport. To some degree this trend is reproduced
   > by the simulations where etch depth does increase with increasing low
   > frequency power however the rate of increase is substantially sublinear.
   > These outcomes indicate that **the effect of ion energy (for example in
   > sputter yield or related processes) might be overestimated in the
   > mechanism.**"

4. **Direct row-vs-published-yield, per channel** (`74eb2fa`,
   `RESULTS_ABSOLUTE_YIELD_2026-08-05.md`). Gray publishes *absolute* yields,
   not only the ratios the joint solve graded, so each channel's magnitude can
   be read separately at 350 eV:

   | quantity | petch | Gray 1993 | ratio |
   |---|---|---|---|
   | floor, F/Ar⁺ → 0 (bare row alone) | 0.341 | **0.28** | **1.22× too strong** |
   | plateau, F-saturated (complex row) | 0.390 | **1.10** | **0.35× — 2.8× too weak** |

   The floor is the bare row read directly — `0.0852 · (350−70)/(140−70) = 0.341`
   — with no coverage or transport interpretation in between.

**This decomposes the bound rather than merely attributing it.** The two ion
channels are wrong in *opposite* directions, which is why every single-parameter
attempt over-corrected something else, and why the model is structurally stuck
ion-limited: with the chemical channel 2.8× too weak it can never dominate the
fluorine-free sputter channel, so neutral-limited behaviour is not reachable by
re-weighting at all. Together the two also reproduce the dynamic-range miss
exactly (measured 0.255 vs petch 0.873 — floor high and plateau low, from both
sides).

**Receipted, not applied.** Both are absolute measurements at a single energy
(350 eV) against a feature running at 3406 eV, so they cannot be landed as
magnitude corrections until the 350 eV → keV scaling is settled — which is
precisely the open ZBL-vs-published-linear question
(`RESULTS_LIMITING_REGIME` §3). Landing a 350 eV correction under an unsettled
keV scaling would be fitting at the reference point and hoping. What they settle
is the *shape* of the departure the depth gate needs, with a published number on
each channel — making a future departure from Appendix B's magnitudes a beam
measurement rather than a tuned constant.

Consistent with this, petch's oxide chemistry measures **exactly ion-limited**
at feature energies — a 50× radical sweep moves the rate 0.7 % while fluorine
coverage collapses 12× (`RESULTS_LIMITING_REGIME_2026-08-05.md`, `830e5c5`) —
where the real process is neutral-transport-limited. An ion-limited chemistry
fed by an AR-flat energetic supply cannot produce ARDE, which is why the missing
ARDE and the depth overshoot are one defect, not two.

**Stop rule.** The joint solve enumerated the full discrete sourced space
(64 combinations, 4 axes) and returned an **empty survivor set**: the required
depth factor 0.735–0.812 falls in the *gap between adjacent published options*
(0.628 | 0.868). Closing it means departing from the published mechanism with a
receipt — i.e. new physics with its own validation — not further transcription.
The campaign therefore closes with the bound declared rather than fitted around.

## 4. Method: four defects found by measurement, not tuning

Each was found because a gate refused, and each is a permanent engine fix that
transfers to every future chemistry.

| defect | receipt | fix |
|---|---|---|
| **√2 axisymmetric lift** — published planar IEAD angle used directly as the polar angle, making the beam 1.41× too narrow | lifted planar σ 0.5893° vs published 0.8334°, ratio **1.4141** (√2 = 1.41421) | Abel/onion-peel inversion of the planar marginal → 0.8233°, inside the digitisation band; deficit ×1.4141 → **×1.0124** (`6e97ef3`) |
| **O-channel normalisation** — a per-cell probability applied per atom | O-removal/deposition **0.0726** against the mechanism's geometry-free **0.1953** (both thermal, view factors cancel) | one line; ratio reproduces 0.1953 to **0.02 %** (`63cfefa`) |
| **crosslink inversion** — creation put on the ion channel, which collapses ~200× on a vertical lip; the source creates crosslinks *at deposition* and *breaks* them with ions | lip x_xl 0.163 measured where the source's lip is the most crosslinked surface | zero new constants (his own Table 6.2 stoichiometry): lip x_xl 0.176 → **0.703**, lip growth 1.759 → **0.831 nm/s** (`3a931b1`), confirmed at feature scale 5.07× → **1.83×** (`391ee43`) |
| **probe transport flag** — the diagnostic itself discarded ~97 % of thermal neutrals while keeping the ion beam | delivered/source 0.0281 → **1.0021** after the fix | conservation gates added so the class cannot recur silently (`cbbd2d6`); all probe-derived numbers from 08-02/04 were retracted |

Discipline notes, because they are the actual product: the **√2 was found by a
convergence harness built to prove the opposite hypothesis**; the wall-slope
theory was **falsified before implementation** (his neck holds at 1.86°, our
shoulder already sits at 10–17°); a "unity cumulative cap" was implemented and
then **killed by the source's own arithmetic** (his 13.75 nm/s requires 3.15
formula units per incident ion — a unity cap allows 1.00 = 4.36 nm/s); and the
angular-convention change was **refuted by the coupled forecast** after passing
in beam mode. Roughly $1 of compute was spent across the entire campaign because
forecast-before-spend killed four would-be runs.

The forecaster earned its status twice: predicted 1.95× / measured 1.83× (ml18),
predicted 0.810 / measured 0.811 (ml20). Balance questions are now answered in
minutes on frozen geometry rather than hours on a GPU.

## 5. Declared-open ledger

Nothing here is fitted around. Each item states what would close it.

**5.1 Mouth residual (+13 %).** The per-material crosslink bond multiplicity is
now sourced — *JVST A* **42**, 043008 (2024) §IV: "CF₂ would have a maximum of
two crosslinks and CF₃ would have a maximum of a single crosslink" — and landed
(`914bb8d`, `924cc04`), moving the lip driver 1.95× → 1.58×. What remains
unpublished is the **crosslinking probability** itself; the MCFPM User's Manual
does not exist publicly (every cpseg page link-extracted,
`RESEARCH_VERIFY_HUNT_2026-08-05.md`).

**5.2 Depth (+29 % at feature energies) — the single named path to closing it.**
§3 decomposes the bound per channel against published absolute yields: the
physical-sputter row is **1.22× too strong** and the chemically-enhanced row
**2.8× too weak** (`74eb2fa`). Closing depth therefore requires, in order:
(i) settling the 350 eV → keV energy scaling — the open ZBL-vs-published-linear
question (`RESULTS_LIMITING_REGIME` §3), which is what licenses carrying a beam
measurement to feature conditions; (ii) landing both channel magnitudes from the
Gray/Kwon beam data under that scaling, with the N1/N2 gates as acceptance; and
(iii) the Kwon/Sawin site-limited adsorption element (E1, and §5.3's `s`), which
is what makes the chemical channel able to become the limiting one. All three
have published sources and preregistered gates
(`RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`); none is a free parameter.
This is a *departure from Appendix B with a receipt* — new physics carrying its
own validation — which is why it belongs to a future campaign rather than this
one.

**5.3 Adsorption coefficient `s`.** Gray's measured half-rise (F/Ar⁺ ≈ 27 ± 8)
implies `s ≈ 0.06` on bare SiO₂ against petch's 1.0. Recorded `[VERIFY]` —
inverted from a published curve, not read from a table, therefore not adopted.

**5.4 Angular form provenance.** Krüger's class-1 citation (Kress 1999) is a
molecular-dynamics study of **Ar⁺ on Cu(111)** at 50–250 eV; its peak/normal of
4.17 is 3× the only in-chemistry measurements (Cho 2000: 1.30; Schaepkens 1998:
1.33) and it cannot reach the lineage's own stated 60° peak at any `B`
(`f′=0 ⟹ cos²θ = (1+B)/(3B)`, bounded by 54.74°). Oxide/mask rows are now bounded
at B = 1.7 (`830e5c5`); the polymer row keeps 9.3 so validated lip results are
untouched. The measured FC-film curve (Barklund & Blom 1992, peak 1.448 at 65°)
uniquely selects the yield reading and is the standing candidate for that row.

**5.5 Class-2 roll-off.** Downgraded from `[VERIFY]` to a quantified
approximation: `min(1, cosθ/cos45°)` agrees with the digitised Chang curve to
**≤ 0.065 absolute**, worst at 70°.

**5.6 Not modelled in this pipeline.** Feature charging (module exists,
notching-validated, not wired into the trench pipeline and unvalidated at deep
AR — gate ladder in `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`);
thermalised-FC return to the radical ledger (element E8, specified with its gate,
not implemented); sputtered-product redeposition beyond the three p = 0.01 rows.

## 6. What transfers

**Engine and instruments — all chemistry-agnostic, all gated:** the √2 lift and
the two-component IADF (`b79e968`); the exact MCFPM reflection cascade with the
derived attenuation law `S = exp(−p₀(1+B)·AR·β²)`, validated to 1 % against the
exact bounce product; the axisymmetric hole operator (0.656 % against Clausing
at 200:1); the STL importer (`02e4a62`); the frozen-geometry forecaster; the
angular-convergence harness; the neck-metric regrader; the preflight harness;
the choke-point strip symmetrisation; sliver dissolution; the `bincount` segment
reduction (4.8×, bitwise).

**Mechanism forms — reusable physics, not this system's numbers:**
deposition-driven crosslinking with ion breaking (the passivation-hardening
mechanism of *any* polymerising chemistry); the element-resolved two-reservoir
mixed layer with exact ledgers; per-event ion chemistry; the critical-wall-angle
lip criterion; site-turnover blocking (gated at > 1390× against the source rows).

**System-specific — the deck, and only the deck:** all Ar/C₄F₆/O₂→SiO₂
constants, now data with provenance (`chemistry_deck.py`, `6499734`/`e97fc95`,
bitwise-gated against the previous hardcoded construction). A new chemistry is a
deck plus a validation campaign, not a code change.

**The one engine-matched compensation, declared:** ml13-class configurations
carried the paper-set O-etch constant (0.0628 vs the converged 0.0423). It is one
number, declared, against the source's own five optimizer-fitted values.

## 7. Reproduction

```
python scripts/grade_ml19_endpoint.py     # §2 gate table
python scripts/grade_ml18_crosslink.py    # closure/etch windows
python scripts/regrade_neck_metrics.py    # neck triple on archived checkpoints
python scripts/gate_n1_beam.py            # Gray 1993 beam gates
python scripts/ion_channel_model_solve.py # §3 joint solve matrix
```

Archived runs: `results/curated/mixed_layer_feature_v1/`. Source extracts and
digitised figure data: `research_sources/` (`555b7fc`). Suite: **1174 passed,
1 skipped**.

---

## Addendum — §5.2 executed (2026-08-05, `79cb426`)

The path this dossier named for closing depth was run the same day. Outcome:
**the regime defect closed, the depth error changed sign.**

Krüger's Appendix B carries `n = 1` (linear, Sigmund) on exactly two rows — the
SiO₂ rows that set the etch — against 416 rows at `n = 0.5` in the same table.
Gray tested that form against √E on this exact system and rejected it (MIT
thesis 1993, Fig. 5-2, p. 161). Replacing those two rows with Gray's own
measured laws, absolutely calibrated and free of any reference-energy
normalisation:

- physical sputter at the measured 3406 eV front: **4.060 → 0.752 per ion**
  (the channel was **5.4× too strong**);
- rate response to a 100× radical shadow: **< 1 % → 0.010×** — the model is
  neutral-limited for the first time;
- early-span rate-vs-depth slope: **+0.7530 → −0.0033** — the anti-ARDE
  pathology is gone;
- aperture behaviour unchanged (86.75 vs 86.00 nm at t = 1 s) — no cross-channel
  leakage;
- depth: **+29 % over → −49 % under**.

The ion-yield question is therefore **settled** and should not be reopened: the
magnitudes now come from the primary source's own absolute measurements, and
both remaining sourced options bracket the requirement rather than meeting it.

Depth is relocated, with a number attached. Gray defines β_e as SiF₄ removed
*from fluorine-saturated regions* per ion; at the front energy that is **2.99**
against the **≈3.15 units/ion** Krüger's blanket arithmetic implies (95 %). The
chemical channel carries the observed etch **at fluorine saturation** — ours is
F-starved. Depth is now a fluorine-delivery-and-coverage question at the etch
front, and the largest named candidate is a *transport* feature: thermalised-ion
return to the radical ledger (E8), since Huang measures **> 95 %** of the floor's
radicals above AR 10 to be thermalised CF_x⁺, a population petch does not
deliver.

Declared-open ledger, updated: `_THERMAL_F_STICKING = 1.0` against Gray's
printed **0.02** and Krüger's absence of any thermal-F-on-bare-oxide row —
sourced, not landed (it would have confounded this grading, and its sign
reduces F further).

---

## Addendum, 2026-08-06 — the 60 s endpoints, measured

Every depth number in this dossier above was either measured under a superseded
mechanism or extrapolated. The final mechanism's own 60 s endpoint had never
been reached: `ml23` was cut at t = 2.80 s and projected linearly, and the
same-day scorecard stalled at t = 1.948 s on a box whose warp build exposed no
CUDA device. All six conditions have now been run to **t = 60.000 s** on a
GPU-verified box (`RESULTS_SCORECARD_ENDPOINT_2026-08-06.md`).

| gate | target | measured endpoint | verdict |
|---|---|---|---|
| mask remaining (all six conditions) | 850 +/- 2 % | 850.20 - 850.40 | **PASS** |
| r(4/6) | [0.84, 0.94] | **0.909** | **PASS** |
| r(8/6) | [0.97, 1.06] | **1.036** | **PASS** |
| necking absent at O2 2.5 | open | 49.60 nm | **PASS** |
| etch depth | 825 +/- 5 % | **346.8 nm** | **MISS (-58.0 %)** |
| mask constriction | 45 +/- 5 | **39.82 nm** | **MISS** (0.18 nm below band) |
| clog at O2 0.5 | sealed | 13.08 nm | **MISS** |
| O2 depth rank max at 1.5 | 1.5 highest | 2.5 highest | **MISS** |

**What the endpoints add to the record.**

The power-transfer result is the substantive gain: scorecard-1 missed both
ratios under the ml9a mechanism (0.672, 1.085) and attributed the misses to the
over-narrowed mouth. Under the final mechanism both land inside their published
bands **at the published endpoint**, within 0.02 of the matched-time reading, so
the transfer behaviour across process conditions is certified rather than
indicated. Mask survival is exact at every oxygen ratio and both powers.

The depth miss is larger than the projection: −58.0 % measured against the
−49 % extrapolated. The projection was optimistic because it drew a straight
line through a still-decelerating span; the aspect-ratio dependence the Gray
laws restored bends the real trajectory below it. This does not change the
attribution — it sharpens the magnitude the decomposed channel bound has to
account for.

The mouth result changed character and should be read with the depth result,
not beside it. Pre-Gray the aperture equilibrated at 50.9 nm (drift −17 pm/s);
it now passes through 45 nm at t ≈ 44 s and ends at 39.82 nm. The run-average
closure/etch is 0.0723 against Krüger's 0.0310. Depth and mouth are one
statement — the etch is too slow relative to lip closure — and the two former
"independent" residuals collapse into the single open channel-magnitude item
already named in §5.2.

Clog at O2 0.5 is directionally preserved (3x narrower than base, still closing
at the endpoint) but does not seal inside 60 s, and is graded MISS on the
published criterion rather than reinterpreted. The 1.0/1.5 oxygen inversion the
matched-time pass recorded at 0.60 nm is 6.2 nm at the endpoint, so it is a real
feature of the mechanism.

No mechanism change was made in this pass. The dossier's §5.2 path — the two
receipted per-channel magnitudes, held pending the 350 eV -> keV scaling — is
unchanged and remains the single named route to closing depth.
