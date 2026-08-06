# The ion-flux anchor: the premise corrected, the bound derived (2026-08-06)

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

Anchoring a boundary on it is therefore not a blanket calibration — it is a
calibration on the feature depth itself. That is worth doing, but it must be
labelled for what it is, and this document does that.

## 1. The arithmetic, step by step

Every number carries its source.

| # | step | value | source |
|---|---|---|---|
| 1 | SiO2 formula-unit density | 2.2e28 m^-3 | deck |
| 2 | experimental depth / time | 825 nm / 60 s | Krüger Fig. 7(b) SEM; thesis Table 6.4 target |
| 3 | required removal flux at the floor | **3.025e20 units m^-2 s^-1** | (2) x (1) |
| 4 | published ion flux | **1.2e20 m^-2 s^-1** | `krueger-2024.txt` L298, Table I, "Ions 1.2 x 10^16 cm^-2s^-1" |
| 5 | required units per incident ion | **2.521** | (3) / (4) |
| 6 | measured per-ion ceiling, >=1 keV, saturating | **1.5 molecules/ion** | Karahashi 2007 Fig. 3, CF3+ (`karahashi_2007_sio2_cfx_ionbeam.txt` L118-127) |
| 7 | hard lower bound on floor ion flux | 2.017e20 m^-2 s^-1 = **1.68x** published | (3) / (6) |
| 8 | ion delivery to an AR~21 floor | 0.70 | cascade funnelling scan |
| 9 | hard lower bound at the wafer plane | 2.881e20 = **2.40x** published | (7) / (8) |

Step 9 uses only measured quantities — an SEM depth, a beam-measured yield
ceiling, and a computed-then-benchmarked delivery fraction. It contains no
petch parameter.

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

## 3. The trilemma

Steps 1-9 are each independently sourced, and they cannot all hold with the
published boundary. Exactly one of the following must be false:

1. the wafer-plane ion flux is the published 1.2e20 m^-2 s^-1;
2. the per-ion removal ceiling is the beam-measured 1.5 molecules/ion
   (corroborated by two independent experiments, Gray 1993 and Karahashi 2004,
   which agree with petch's channels to 4.7%);
3. the experiment removed 825 nm of SiO2 in 60 s at these boundary conditions.

The source reaches (3) while holding (1) only by fitting `ps,SiO2` over the
optimizer range 0.0-0.3 (`krueger_thesis_2024.txt` L4570-4583) inside an energy
law that both beam experiments contradict. Under the measured sqrt(E) law the
same row yields 0.752 units/ion where his converged value yields 4.06.

Of the three, (1) is the only one with no measurement behind it: the ion flux is
HPEM model output (`evidence_type: HPEM_simulation`) from a reactor model the
thesis itself treats as ground truth without experimental validation. That is
why the normalization implemented here acts on the ion flux and never on the
yields.

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

## 5. Forecast

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

The measured 60 s depth of 346.8 nm needs 2.379x to reach 825 nm, which the
frozen-geometry response reaches at **2.83x**. The coupled run carries mouth
feedback the frozen estimate cannot, so the realized scale is expected to be
lower; the run measures it.

Note this also corrects a supply-bound reading recorded earlier: scaling the
*yield magnitude* saturates by 4x (10% at 8x), but scaling the *boundary flux*
is near-linear at the floor. Both are true; they are different levers, and only
the second has authority over the rate.
