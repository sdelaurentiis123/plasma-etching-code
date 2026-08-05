# The depth gate is an energetic-delivery question, not a chemistry question

Receipts: `scripts/floor_delivery_scan.py`,
`scripts/forecast_bond_multiplicity.py`,
`results/curated/mixed_layer_feature_v1/ml19-depxl-60s/audit.json`.

## Correction to the aspect ratios quoted so far

The etch front sits below an 850 nm mask, so ml19's *first* record is already
AR ≈ 11.7 (929 nm below the mask top through a 79 nm aperture) and its last is
AR ≈ 32.  Earlier notes that labelled these "AR 1" and "AR 15" used the etched
depth alone.  Nothing downstream of the label changes, but the comparison to
Huang's funnelling curve has to be made at the corrected values.

## The floor is exactly ion-limited, and linear

Oxide recession relaxed to steady state over the delivery plane (nm/s):

| ion / source | neutral 0.02 | 0.05 | 0.10 | 0.30 | 0.60 | 1.00 |
|---|---|---|---|---|---|---|
| 1.0 | 7.59 | 7.58 | 7.58 | 7.56 | 7.55 | 7.54 |
| 2.0 | 15.18 | 15.18 | 15.17 | 15.15 | 15.12 | 15.10 |
| 3.0 | 22.78 | 22.77 | 22.76 | 22.73 | 22.70 | 21.42 |

Across a **50× range of neutral delivery the rate moves under 1 %**, while it
is **exactly proportional to ion flux**.  The floor carries no F-starvation
sensitivity at these conditions: it is a pure energetic-delivery observable.

## What ml19's measured rates therefore imply

| quantity | rate | implied energetic flux / source ion flux |
|---|---|---|
| ml19 at AR ≈ 11.7 | 21.91 nm/s | **2.9×** |
| ml19 at AR ≈ 32 | 16.71 nm/s | **2.2×** |
| Krüger run average | 13.75 nm/s | **1.8×** |
| 0-D blanket, full delivery | 7.54 nm/s | 1.0× (by construction) |

A shadowed floor cannot receive more *neutral* flux than the open field, and
the scan shows neutrals are irrelevant anyway.  So the feature's floor is
being supplied 2.2–2.9× the source ion flux — the funnelled cascade (grazing
ions converted to hot neutrals that re-arrive at the floor).  Krüger's own
depth implies 1.8× through the same mechanism, so **our energetic delivery to
the floor runs 1.2–1.6× above his**, and that ratio is the entire depth
overshoot (measured 1.215× on the late rate).

**Consequence for the fix list.** No chemistry constant can close the depth
gate without breaking the ion-linearity the scan just established — a change
that slows the floor by 1.2–1.6× would slow it by the same factor at every
delivery, including the blanket rate that is not in question.  The depth gate
is owned by the cascade's floor delivery, which is testable against Huang's
published funnelling curve (ion 2.0 → 0.3 ×10¹⁵ and hot neutrals 3.1 → 8.0 →
1.1 ×10¹⁵ over AR 0 → 4 → 40): his total energetic flux *falls below* the
incident ion flux by AR 40, where ours is still at 2.2× at AR 32.

## Crosslink partner count — now sourced, no longer `[VERIFY]`

Krüger et al., *JVST A* **42**, 043008 (2024), sec. III
(`tmp/pdfs/krueger-2024.txt` L388-390), verbatim:

> "...based on the number of available bonds (three in the example in Fig. 5).
> For example, CF2 would have a maximum of two crosslinks and CF3 would have a
> maximum of a single crosslink."

Thesis sec. 2.2.3 (`tmp/pdfs/krueger_thesis.txt` L2475-2476) states the same
basis without the worked examples.  Those two examples pin the rule uniquely:

    available = 4n - m - 2(n - 1)   for C_nF_m

reproducing CF → 3, CF2 → 2, CF3 → 1 exactly.  Applied to the published
stoichiometry and passed to the layer as the deposition-weighted mean.
Nothing is fitted.

Effect at the audited lip delivery: x_xl 0.703 → **0.793**, film growth
0.831 → **0.673 nm/s**, i.e. **1.95× → 1.58×** Krüger's 0.427 nm/s per-side
closure.  The mouth residual shrinks by the same factor it predicted.

## Still open

- Floor energetic delivery (above) — the depth gate.
- The ∠=2 roll-off *shape* stays `[VERIFY]`: Huang (`huang_thesis.txt`
  L2295-2296) and Huard (`huard_chad_phd_thesis.txt` L3020) state only
  "unity for normal incidence and angles up to 45°, with a monotonic roll-off
  to zero probability at grazing incidence"; no source reprints the table.
- 60 s runs at dx = 10 nm stop on a 2-cell mask sliver (remeshing artifact,
  ml19 at t = 46.2 s); a sub-resolution fragment reabsorption rule is needed
  before an endpoint run can complete.
