# Anchoring the two ion channels to Gray's absolute yields

Follow-on to `RESULTS_ABSOLUTE_YIELD_2026-08-05.md`, which measured the two SiO2
ion channels against Gray 1993 and found them wrong in opposite directions.
This pass forecasts what happens if each channel's MAGNITUDE is re-anchored to
that measurement, each row keeping its own published energy FORM.

## 1. A hypothesis of mine, refuted first

Before the anchoring test I proposed that the energy FORM was the root cause:
that the canonical Steinbruchel `sqrt(E) - sqrt(Eth)` law (transfer factor 4.8x
from 350 eV to 3406 eV, against 9.9x for ZBL and 11.9x for the Appendix-B
linear form) would simultaneously fix the depth overshoot, invert the
bare/complex ratio, and restore neutral limitation.

**Forecast, and it is refuted.**  Applying `sqrt(E)` to both rows:

| observable | current | sqrt(E) | measured / needed |
|---|---|---|---|
| Gray floor (350 eV) | 0.341 | **0.254** | 0.28 |
| Gray plateau | 0.390 | 0.318 | 1.10 |
| Gray dynamic range | 0.873 | 0.800 | 0.255 |
| coupled depth factor | 1.000 | **0.322** | 0.735-0.812 |
| neutral response (0.05x radicals) | 1.017 | 0.967 | << 1 |
| ARDE (AR16/AR0) | 0.972 | 0.969 | < 1 |

It improves the bare-row magnitude and nothing else: the plateau stays ~3.5x
low, the depth *overshoots* the correction by more than 2x, and the regime does
not move.  The reasoning error was mine — both rows carry similar thresholds, so
changing the exponent rescales them **together** and the bare/complex ratio is
preserved (1.06 after, 0.96 before).  The predicted inversion cannot happen.
Caught for free, before any run.

## 2. What the data actually supports: magnitude anchoring

Gray publishes two absolutes at 350 eV — floor **0.28** (F-free, physical
sputter alone) and plateau **1.10** (F-saturated).  Those fix one scale per
channel, each row keeping its own published energy form:

    bare row    x 0.821   (0.341 -> 0.280)
    complex row x 2.82    (0.390 -> 1.095)

| observable | current | Gray-anchored | measured / needed |
|---|---|---|---|
| Gray floor | 0.341 | **0.280** | 0.28 (anchor) |
| Gray plateau | 0.390 | **1.095** | 1.10 (anchor) |
| **Gray dynamic range** | 0.873 | **0.256** | **0.255 — free check, passes to 0.4%** |
| neutral response, 0.05x radicals | 1.017 (flat) | **0.832** | responsive |
| neutral response, 0.3x radicals | — | 0.932 | responsive |
| coupled depth factor | 1.000 | 1.045 | 0.735-0.812 **FAIL** |
| ARDE (AR16/AR0) | 0.972 | 0.968 | << 1 **FAIL** |

Two anchors, three measurements: **the dynamic range is not fitted, it follows,
and it lands within 0.4% of Gray's measured value.**  That is an independent
validation of the two-channel structure petch carries — the structure is right,
the magnitudes were wrong.

And it moves the regime for the first time in the campaign: starving radicals
20x now costs 17% of the etch rate, where every previous configuration was flat
to 1.7% (`RESULTS_LIMITING_REGIME` §2).  The model is no longer structurally
ion-limited.

## 3. What it does NOT fix, and where that leaves the two open gates

Depth stays at 1.045 (needs 0.735-0.812) and ARDE stays at 0.968.  Taken with
everything else this campaign has eliminated, the map is now complete:

| candidate | status |
|---|---|
| transport / cascade | exonerated by the matched-beam audit (`fd61bb7`) |
| angular assignments | joint solve: no sourced combination reaches the depth band (`f322cb4`) |
| grid resolution | acquitted, 0.05% converged floor rate (`ab5b880`) |
| chemistry rows | closed; every Appendix-B row transcribed and audited |
| **channel magnitudes** | **anchorable to Gray, fixes the REGIME, does not fix depth or ARDE (this pass)** |
| ion-energy channel magnitude | the residual depth term; Krueger flags it himself (L4884-4888) |
| **neutral delivery vs depth** | **the residual ARDE term — see below** |

The ARDE number localises sharply now that the chemistry responds to radicals.
Our thermal delivery to the floor falls only **12%** from AR 0 to AR 16
(`results/curated/cascade_funnelling/scan.json`), while Huang's published ladder
has N/I falling **~75x** down a real feature (L5430-5478).  A 17%-per-20x
chemical sensitivity against a 12% delivery change can only produce ~1% of rate
change — which is exactly the 0.968 measured.  **The missing ARDE is a neutral
CONDUCTANCE defect, not a chemistry defect**, and the named unimplemented item
sits precisely there: E8, the thermalised fluorocarbon ions that Huang measures
as >95% of the floor radical flux above AR 10 and that petch's cascade currently
discards from the ledger entirely (`RESULTS_ANGULAR_CONVENTION` §8).

## 4. Status: forecast only, nothing landed

The anchoring is **not applied**.  Two reasons, both doctrine:

1. The anchors are measured at **350 eV** against a feature front at **3406 eV**.
   Carrying them requires the energy form, and §1 just demonstrated how badly a
   form change can behave when its consequences are not forecast coupled.  The
   form question is out for sourcing.
2. Re-anchoring is a **departure from Appendix B's magnitudes**.  It is a
   receipted departure — Gray is a published beam measurement and the dynamic
   range check passes independently — but it should land as one deliberate,
   graded change with the energy form settled, not as a scale factor bolted onto
   an unsettled scaling.

What the pass establishes is that the departure is *warranted and testable*: the
two-channel structure survives an independent check, the regime defect is a
magnitude error rather than a missing mechanism, and the two remaining gates
have distinct, named owners — ion-energy magnitude for depth, neutral
conductance (E8) for ARDE.
