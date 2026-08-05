# ml18: deposition-driven crosslinking at feature scale — the early transient closes

Confirmation run for `3a931b1` (deposition-driven crosslinking, the lip
inversion).  Graded per `RESULTS_EARLY_TRANSIENT_2026-08-04.md`: **12 s, not
60 s** (88 % of the defect was complete by t = 8 s), on `closure/etch` by
window against Krüger's run-average **0.0310** — never on aperture at 60 s.

Run: `ml18-depxl-12s`, dx = 10 nm, 83 steps to t = 12.000 s, RTX 3090,
clean `git archive` deploy (`fresh_fraction` verified present post-extraction),
box destroyed.  Artefacts: `results/curated/mixed_layer_feature_v1/ml18-depxl-12s/`.
Reproduce the grade: `python scripts/grade_ml18_crosslink.py`.

## Preregistered gate: closure/etch by window

| window | ml18 | × Krüger | ml16a baseline | × Krüger | improvement |
|---|---|---|---|---|---|
| 1–4 s | **0.0569** | 1.83× | 0.1572 | 5.07× | **2.76×** |
| 4–8 s | **0.0495** | 1.59× | 0.1091 | 3.52× | **2.21×** |
| 8–12 s | **0.0367** | 1.18× | 0.0526 | 1.69× | **1.43×** |

The declared pass band was `0.031–0.1`.  **All three windows land inside it**,
and the third is within 18 % of his run-average.  The 5.1× / 3.5× / 1.7×
early-transient signature that defined the defect is gone.

## Matched-time trajectory

| t (s) | opening ml18 | baseline | Δ | depth ml18 | baseline | Δ |
|---|---|---|---|---|---|---|
| 1 | 86.00 | 80.96 | +5.04 | 23.3 | 22.7 | +2.6 % |
| 2 | 83.35 | 73.14 | +10.21 | 45.3 | 44.7 | +1.4 % |
| 4 | 78.54 | 60.61 | +17.93 | 88.8 | 87.4 | +1.6 % |
| 8 | 70.32 | 44.79 | +25.53 | 172.0 | 159.9 | +7.5 % |
| 12 | **64.89** | 38.25 | **+26.64** | **246.0** | 222.1 | **+10.7 %** |

Closure budget spent by t = 12 s: **49 %** of Krüger's full-run 51.2 nm, against
the baseline's **101 %**.  The baseline reached his *final* neck value five
times early and kept going; ml18 is on a trajectory that leaves room for the
remaining 48 s.

Depth rises 10.7 % without any change to the floor chemistry — the wider mouth
delivers more flux, the coupling the neck regrade (`fcc98eb`) predicted when it
showed depth is inherited from the throttled aperture.  Mask remains
**850.21 nm** (target 850, unchanged and exact).

The 0-D probe predicted this within its own error bar: it forecast the lip
going from 4.1× to **1.95×** too fast, and the feature run delivered **1.83×**
in the matching window.  That is an independent validation of the frozen-geometry
probe as a forecasting instrument, not just a diagnostic one.

## Residual, and the integer that bounds it

The first window remains 1.83× his ratio.  Per
`RESULTS_LIP_CROSSLINK_2026-08-04.md` §5 the implemented form is the *minimal*
transcription — one bond per deposition event, steady state `x = (1+m)/(2+m)`
with `m = 1`, giving `x = 2/3`.  His module allows a per-material maximum bond
count that the thesis never tabulates (it appears only as "3 in the example
depicted in Figure 2.2"):

| m | lip x_xl | est. first-window ratio |
|---|---|---|
| 1 (implemented) | 0.667 → measured 0.703 | **1.83×** (measured) |
| 3 (his figure's example) | 0.80 | ≈ 1.4× |
| 8 | 0.90 | ≈ 1.0× |

The residual is bounded by that one undetermined integer, interpolated on the
doc's own 0-D family (x = 0.703 → 0.831 nm/s, x = 0.9 → 0.427 nm/s, Krüger
0.427).  Adopting m > 1 without a tabulated source would be fitting, so it
stays `[VERIFY]`.

## What this does not yet show

`neck_cd_nm = 27.78 at 1090 nm below mask top` is the *etch-front* taper at the
trench bottom (floor at 1.55 µm, mask exit at 1.81 µm), not a mask neck — at
t = 12 s only 246 nm is etched, so the run carries no mask-neck observable
comparable to his 60 s 38.8 nm @ 271 nm.  A 60 s run is still required for the
endpoint metrics (neck CD, neck depth, final depth 825 ± 5 %), and its
preregistered expectation is now: aperture trajectory continuing from 64.89 at
t = 12 with closure/etch held near 0.031–0.037.
