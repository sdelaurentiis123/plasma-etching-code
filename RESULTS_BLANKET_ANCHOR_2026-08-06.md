# Historical ion-flux anchor — bound retracted 2026-08-06

> **RETRACTED BOUND.** The former `2.40x` “measurement-only” ion-flux lower
> bound depended on treating Karahashi's rounded `1.5 SiO2/ion` CF3+ value as a
> universal ceiling and on applying a final-geometry delivery fraction to the
> full 60 s history. Both steps are invalid. See
> `RESULTS_DEPTH_IDENTIFIABILITY_2026-08-06.md`.

The mandate was to anchor a boundary-flux normalization on Krüger's "blanket
etch rate" of 13.75 nm/s, on the reading that his Table-I ion flux is
understated ~2x against his own blanket data.

**There is no blanket datum.** 13.75 nm/s is not a published observable: it is
`825 nm / 60 s`, the feature depth target divided by the etch time, and every
one of the twelve places the repo uses it derives it that way
(`scripts/forecast_bond_multiplicity.py:39` is explicit:
`KRUEGER_AVG = 825.0 / 60.0`). The string "13.75" appears in **zero** of the 24
archived source extracts, and no blanket or unpatterned etch rate appears
anywhere in the thesis, the paper, the OSTI manuscript, or the 2022 companion.

Anchoring a boundary on it is therefore not a blanket calibration. It is a fit
to the feature depth itself. The `ion_flux_normalization` implementation is
retained as an explicit sensitivity/target-fit axis, defaulting to bitwise
identity, but it has no independently measured calibration value.

## 1. The arithmetic, step by step

Every number carries its source.

| # | step | value | source |
|---|---|---|---|
| 1 | SiO2 formula-unit density | 2.2e28 m^-3 | deck |
| 2 | experimental depth / time | 825 nm / 60 s | Krüger Fig. 7(b) SEM; thesis Table 6.4 target |
| 3 | depth-integrated average removal flux | **3.025e20 units m^-2 s^-1** | (2) x (1) |
| 4 | published ion flux | **1.2e20 m^-2 s^-1** | `krueger-2024.txt` L298, Table I, "Ions 1.2 x 10^16 cm^-2s^-1" |
| 5 | target normalized by wafer-plane ions | **2.521** | (3) / (4), a run-average lower-bound normalization |
| 6 | simulated 346.833 nm on the same basis | **1.060** | endpoint depth normalized by (1), (4), and 60 s |
| 7 | unresolved effective gap | **1.461** | (5) - (6) |
| 8 | final-geometry ion delivery diagnostic | 0.70 | cascade-funnelling scan; not a 60 s history average |
| 9 | counterfactual target if 0.70 held for all 60 s | **3.601 per delivered floor ion** | (5) / (8), diagnostic only |

There is no hard ion-flux lower bound in this table. Karahashi Figure 4 reaches
`1.8736` for CF3+ at 1500 eV, and Takada measures `2.5` for stable C5F8/Ar+
co-incidence at 900 eV and ratio 1. The latter is a different molecule and not
a C4F6 law, but it proves the old ceiling premise false.

## 2. The cross-check that refutes the "understated 2x" premise

The published ion flux is **not** anomalous for its reactor class:

| quantity | Krüger | Huang (`huang_thesis.txt` L5270-5280) |
|---|---|---|
| total ion flux | 1.2e16 cm^-2s^-1 | ~7.7e15 cm^-2s^-1 |
| reactive-FC / ion ratio | 25.8 | 30 (stated base case) |

Krüger's ion flux is 1.56x Huang's in absolute terms and his neutral-to-ion
ratio sits just below Huang's — both consistent with a higher-power reactor.
Nothing in the peer literature marks the published value as low. The
"understated 2x" premise has no independent support, and the normalization
implemented here is therefore a **declared calibration**, not an inference.

## 3. Why the trilemma failed

The former trilemma omitted two variables:

1. Krüger's positive ions are published only as an aggregate population even
   though reactive-ion yield depends strongly on ion identity.
2. Table I omits stable parent C4F6 even though direct beam experiments show
   that a surviving fluorocarbon molecule can add an order-one,
   non-monotone ion-assisted channel.

The ion flux is HPEM output rather than a direct measurement, so scaling it is
a legitimate sensitivity. It is not the only unmeasured variable and is not
singled out by a valid surface-yield ceiling.

## 4. What was implemented

`ion_flux_normalization` on both the base and transfer boundary builders,
default 1.0 and inert at the identity by construction. It scales the aggregate
positive-ion flux only; every neutral flux is bitwise unchanged, and a scaled
flux carries a new first-class evidence kind, `declared_calibration`, so no
audit can mistake it for a measured or model-published value.

One bug was caught by direct wiring verification before any box time: the
Figure-16 oxygen table supplies its own fluxes and discarded the normalization,
so the O2 conditions would have run uncalibrated while the power conditions ran
calibrated — the two halves of the out-of-sample scorecard graded under
different boundaries. Fixed and gated for all five transfer conditions.

## 5. Historical target-fit forecast

Frozen-geometry floor rates, ion-only against uniform scaling
(`scripts/forecast_ion_flux_normalization.py`):

| scale | ion-only | ratio | uniform | ratio |
|---|---|---|---|---|
| 1.0 | 4.100 nm/s | 1.000 | 4.100 | 1.000 |
| 2.0 | 7.238 | 1.766 | 8.199 | 2.000 |
| 2.4 | 8.457 | 2.063 | 9.839 | 2.400 |
| 2.8 | 9.663 | 2.357 | 11.479 | 2.800 |
| 4.0 | 13.110 | 3.198 | 15.398 | 3.756 |

Ion-only is the mode implemented, for a physical reason: uniform scaling raises
deposition in lockstep with removal, so the closure/etch ratio (measured 2.33x
Krüger's) does not improve, while ion-only raises removal against a fixed
depositor flux and so acts on the depth gate and the closure ratio together.

The measured 60 s depth of 346.8 nm needs 2.379x to reach 825 nm, and the
frozen-geometry response reaches that rate ratio near an ion scale of 2.83x.
This remains a sensitivity result only. A coupled run at that scale would fit
the target and cannot be advertised as a prediction or a blanket-anchored
calibration.

Note this also corrects a supply-bound reading recorded earlier: scaling the
*yield magnitude* saturates by 4x (10% at 8x), but scaling the *boundary flux*
is near-linear at the floor. Both are true; they are different levers, and only
the second has authority over the rate. Neither identifies the missing
reactive-ion mixture or stable-parent channel.
