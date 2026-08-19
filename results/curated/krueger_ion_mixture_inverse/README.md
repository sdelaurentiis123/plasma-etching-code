# Krüger ion-mixture inverse and feature-prefix gate

Status: **the published IEAD cannot identify the ion mixture, and the first
independently measured C4F6 mixture band does not close the feature-prefix
rate**.

## What changed

The feature engine now accepts an explicit fractional composition for
Krüger's aggregate `ions` population. The total ion flux, combined IEAD,
neutral fluxes and mask mechanism stay unchanged; only the atom inventory
coincident with ion impacts changes. Fractions may sum below one, with the
remainder retained as non-incorporating ions. Every such run remains a
declared sensitivity rather than a predictive Krüger boundary.

The inverse audit also corrects a tempting visual misread. Figure 16's
approximately 250 eV energy ladder is the committed digitizer's 250 eV binning
operator, not a set of species peaks. Krüger's thesis reports only about a
60 eV O+ to CF3+ mean-energy shift in a related waveform-tailoring condition.
That is 0.24 of one digitization bin. More fundamentally, the paper publishes
one combined IEAD per power and no species kernels, so the published
composition-contrast operator has exact rank zero.

## Source-conditioned feature result

Three 10 nm, 0.25 s feature runs used Benck et al.'s mass-resolved C4F6/Ar
measurements at 10 mTorr. No Krüger depth or flux multiplier selected the
mixtures.

| composition sensitivity | depth at 0.25 s | prefix rate | error vs 13.75 nm/s scale |
|---|---:|---:|---:|
| Benck 50% C4F6, unresolved ions inert | 2.9535 nm | 11.8139 nm/s | -14.08% |
| Benck 50% C4F6, all unresolved ions assigned CF3+ | 3.0363 nm | 12.1452 nm/s | -11.67% |
| Benck 75% C4F6, all unresolved ions assigned CF3+ | 3.1355 nm | 12.5420 nm/s | -8.79% |

The F-rich assignments are upper sensitivities, not measurements: Benck's
unresolved current includes other fragments and etch-product ions. The
experiment is also an ICP with no oxygen and a grounded diagnostic surface,
not Krüger's high-power, multifrequency C4F6/Ar/O2 CCP. Even so, the entire
source-conditioned band misses the required run-average rate before deep
feature transport has developed.

Pure-CF2+ and pure-CF3+ endpoints still bracket the required rate, as recorded
in `../guo_krueger_ion_identity_envelope/README.md`, but this audit shows that
those extreme endpoints are not supported by the closest quantitative C4F6
composition measurement.

## Verdict

The code can now propagate any proposed mixture and can prove that the
published plot does not uniquely determine one. It has **not** recovered a
Krüger mixture or predicted 825 nm. A defensible closure now requires either:

1. the original species-resolved HPEM/PCMCM wafer output from Krüger/Kushner;
   or
2. a C4F6/Ar/O2 reactor model first graded against Benck's mass-resolved
   composition and pressure response, then applied to Krüger without using the
   825 nm endpoint.

Machine-readable receipts are in `audit.json` and `feature_prefix.json`.
