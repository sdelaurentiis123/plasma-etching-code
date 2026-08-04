# Neck regrade — Top CD / Neck CD / neck depth on the archived runs (2026-08-04)

Campaign 5 closed on an ambiguity: `mask_opening_nm` minimises over the mask
band and returns one scalar, so a run that is too narrow and a run that necks
at the wrong depth score identically.  The frozen-geometry probe
(`RESULTS_MOUTH_EQUILIBRIUM_PROBE_2026-08-02.md`) measured closure 30-50x
faster at the mask top than at the 200-250 nm band where the SEM and MCFPM
neck, which predicts exactly that confusion.  The pilot now emits the
community triple (`top_cd_nm`, `neck_cd_nm`, `neck_z_um`,
`neck_depth_from_mask_top_nm`) plus a coarse `aperture_profile`, and
`scripts/regrade_neck_metrics.py` replays archived checkpoints through it.

## Targets (digitised from Fig. 7, `tmp/mouth_profiles/`)

| source | neck CD | neck depth below mask top |
|---|---|---|
| Krüger MCFPM simulation | **38.8 nm** | **271 nm** |
| SEM experiment | **39.0 nm** | **200 nm** |

Digitisation caveat, verified this pass: the experimental trace dives to
27.0 nm at 26 nm depth and rebounds to 90.6 nm by 40 nm — a mask-corner
clipping artifact in the top ~50 nm.  Excluding the top 50 nm reproduces the
39.0 nm / 200 nm figure quoted in `RESEARCH_MOUTH_MECHANISM_KRUEGER_2026-08-02`.
The simulated trace has no such artifact (38.8 nm / 271 nm with or without the
exclusion).  The experimental profile then **re-opens** below its neck
(61.2 nm at 230, 83.5 at 270, 89.9 at 350) while the simulated one does not
(41.6 / 41.1 / 50.1) — Krüger's own unreproduced-taper caveat.

## Regrade

| run | top CD | neck CD | neck depth | mask-band neck | 200-270 nm band | etch depth |
|---|---|---|---|---|---|---|
| ml13-base-cascade | 86.8 | 18.7 | 1558 nm (in trench) | 24.8 @ 228 nm | 24.8-32.1 | 852.1 |
| ml16a-verbatim-lift | 62.3 | 11.1 | 130 nm | 11.1 @ 130 nm | 19.5-39.1 | 590.5 |
| ml16b-ml13c-lift | 62.7 | 12.3 | 120 nm | 12.3 @ 120 nm | 18.2-36.0 | 638.8 |

**The legacy scalar was not comparing like with like.** ml13's 24.8 nm is a
*mid-mask* neck at 228 nm; ml16a/b's 11.1/12.3 nm are *top-region* pinches at
120-130 nm.  The apparent "24.8 → 11.1 regression" conflates a size change with
a location change.  Caveat on ml13's apparent location agreement: its mask
eroded 132 nm (`mask_top_z_um` 2.518 vs the intact 2.650), so measured from the
*original* mask plane its neck sits at 360 nm, below the target band.  The
ml16 runs hold the mask at 2.650 exactly (armoured a-C, 850 nm target met).

## Aperture profile vs depth — where the error actually lives

| depth below mask top | ml13 | ml16a | ml16b | Krüger sim | SEM |
|---|---|---|---|---|---|
| 50 nm | 64.7 | 34.0 | 28.8 | 80.1 | 85.0 |
| 100 nm | 42.1 | 15.9 | 13.0 | 62.5 | 66.2 |
| 130 nm | 38.4 | 11.1 | 13.7 | 56.7 | 57.4 |
| 200 nm | 27.3 | 19.5 | 18.2 | 45.0 | 39.9 |
| 230 nm | 24.8 | 29.5 | 25.8 | 41.6 | 61.2 |
| 270 nm | 32.1 | **41.8** | 38.7 | **41.1** | 83.5 |
| 350 nm | 53.6 | 57.1 | 55.0 | 50.1 | 89.9 |

Two structural findings, neither visible in the legacy metric:

1. **The defect is confined to the top ~150 nm of the mask.** At 270 nm ml16a
   reads 41.8 nm against the simulated 41.1 nm, and at 350 nm 57.1 vs 50.1 —
   agreement at the few-nm level in exactly the band where Krüger necks.  The
   miss is a factor 3-4 at 100-130 nm and vanishes by 270 nm.
2. **Our aperture profile has the wrong shape, not merely the wrong size.**
   Krüger's simulated profile narrows monotonically from 88 nm at the top to
   its neck; every petch run instead forms a constriction high in the mask
   (100-130 nm) and then *widens* below it.  The corrected sqrt-2 lift
   (ml16a/b) sharpens that high constriction — top CD falls 86.8 → 62.3 —
   which is the mechanism behind campaign 5's apparent regression.

## What this licenses

- Grade the mouth on `neck_cd_nm` **and** `neck_depth_from_mask_top_nm`; a
  single aperture cannot separate the two failure modes and previously did not.
- Aim the next fix at the **top-of-mask lip balance** (probe: closure
  30-50x faster there; here: the entire aperture error).  The deep-mask
  chemistry already reproduces the target band and must not be disturbed.
- Krüger's simulated profile, not the SEM's, is the shape target above 200 nm:
  the SEM's re-opening below its neck is unreproduced by his model too.

Artifacts: `scripts/regrade_neck_metrics.py`,
`tests/test_neck_metrics.py`, pilot `final_metrics` gains the triple plus
`aperture_profile`.
