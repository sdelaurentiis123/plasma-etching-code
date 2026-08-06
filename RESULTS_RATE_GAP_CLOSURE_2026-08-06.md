# Historical rate-gap investigation — corrected 2026-08-06

> **RETRACTED CONCLUSION.** This document formerly claimed that Karahashi's
> rounded 1000 eV CF3+ yield (`1.5 SiO2/ion`) independently validated petch's
> neutral-assisted complex-channel magnitude to 4.7% and supplied a universal
> physical ceiling. Both claims were false. The authoritative correction is
> `RESULTS_DEPTH_IDENTIFIABILITY_2026-08-06.md`.

## What remains valid

The measured endpoint remains `346.833 nm` versus an `825 nm` target, a
`2.3787x` rate/depth gap under the published aggregate boundary.

The incumbent neutral-assisted complex channel is supply-bounded at the
final-geometry floor. Scaling only its Gray coefficient gives:

| complex-yield scale | floor rate | versus base |
|---:|---:|---:|
| 1x | 4.100 nm/s | 1.00 |
| 2x | 4.390 nm/s | 1.07 |
| 4x | 4.518 nm/s | 1.10 |
| 8x | 4.518 nm/s | 1.10 |

That result, gated in `tests/test_rate_gap_supply_bound.py`, rejects one narrow
proposal: multiplying the existing complex-channel coefficient cannot be the
sole fix. It does not bound missing reactive-ion, stable-parent,
impact-fragmentation, or molecule/ion co-incidence channels.

The following diagnostics also remain valid within their scopes:

- atomic F over the sourced band barely moves the final-geometry frozen-floor
  rate;
- the existing activated-site population self-limits;
- the implemented post-wall-collision thermalized-ion return is small;
- ballistic transport passes its analytic and Lam AR-50 checks.

## What failed

At 1000 eV the old end-to-end mechanism returned the same
`0.380584 SiO2/ion` for F+, CF+, CF2+, and CF3+ because energetic ion identity
was discarded. Comparing two internal formula values (`0.381` and `1.570`) to
Karahashi's F+ and CF3+ markers was therefore not an end-to-end validation.

The visually audited Figure-4 data instead give:

| ion at 1000 eV | measured SiO2/ion |
|---|---:|
| F+ | 0.3232 |
| CF+ | 0.6751 |
| CF2+ | 1.1957 |
| CF3+ | 1.4703 |

CF3+ reaches `1.8736` at 1500 eV. The 1000 eV rounded value was never a hard
ceiling, and a pure-ion radical-free experiment cannot bound molecule-assisted
plasma removal.

Takada's archived radical-free C5F8/Ar+ experiment directly measures that
missing class: `~1.2` at 400 eV and `2.5` at 900 eV, both at molecule/ion ratio
1. It cannot be transplanted to C4F6, but it invalidates the former
impossibility argument.

## Corrected target normalization

`825 nm / 60 s`, SiO2 density `2.2e28 m^-3`, and published aggregate wafer ion
flux `1.2e20 m^-2 s^-1` give `2.5208 SiO2` per **wafer-plane ion** as a
run-average lower-bound normalization. It is not a yield “sustained at the
floor.” Applying the final diagnostic ion delivery of 0.70 to all 60 s would
instead give a counterfactual `3.6012` per delivered floor ion; the evolving
feature requires a time integral.

## Corrected verdict

The rate gap is not closed and is not proven impossible. Krüger publishes:

- an aggregate positive-ion flux and combined IEAD, but no ion-species
  composition;
- no stable C4F6 wafer flux;
- no direct C4F6 molecule/ion surface law.

Those omissions are exactly the variables to which direct measurements show
order-one sensitivity. The current result is therefore an underidentified
boundary/mechanism problem. The code now includes an opt-in, support-refusing
Karahashi species-resolved closure and a separate Takada analog evidence table;
neither is silently applied to the Krüger feature.
