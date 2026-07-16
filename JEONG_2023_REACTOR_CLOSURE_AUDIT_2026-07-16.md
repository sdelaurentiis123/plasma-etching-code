# Jeong 2023 reactor-to-feature closure audit

Date: 2026-07-16

## Outcome

The current feature engine is not the limiting component in the Jeong density sweep. The
public experiment does not provide the species-resolved ion flux, IEAD, complete radical
wall-flux vector, or hot-neutral boundary needed to drive the feature model predictively.
The historical all-Ar/Bohm boundary was a development proxy, not a reactor model.

No full moving-profile matrix is earned from the present evidence. This is a deliberate
bounded stop, not a numerical failure.

Machine-readable evidence:
`results/jeong_2023_reactor_closure_audit/audit.json`.

Static diagnostic:
`results/jeong_2023_reactor_closure_audit/closure_audit.png`.

## What was implemented

1. `Jeong2023IonBoundaryClosure` now permits either:
   - an explicit positive-ion density mixture with species-specific Bohm fluxes; or
   - species-resolved fluxes supplied directly by a diagnostic or validated reactor model.
2. The Jeong runner classifies every positive ion as energetic bombardment and reports
   species-resolved boundary flux, energy, hit, and escape diagnostics.
3. The existing grazing-reflection law remains single-species. A multispecies Jeong boundary
   must disable it until each ion has a declared material-response law; applying the Ar law to
   fluorocarbon ions silently is forbidden.
4. The reduced Huang--Kushner mechanism now includes the published 5--70 eV direct
   fluorocarbon-ion polymer-deposition channel on oxide-fluorocarbon complexes. The channel is
   species-selective, conservative, and included in the material ledger.
5. The Jeong summary now compares archived implementation checksums to the current operator.
   Historical runs are automatically labeled `historical_stale_operator` rather than presented
   as current validation.

## Quantitative findings

### Historical result is useful but stale

The archived 200 nm campaign used the older implicit substrate-polymer projection:

- experimental energy-sweep endpoint gain: 830.58 nm;
- historical prediction: 863.51 nm;
- experimental density-sweep endpoint gain: 451.47 nm;
- historical prediction: 1176.45 nm.

The energy trend was encouraging. The density response was much too steep. Because the chemistry
and boundary operator have since changed, these numbers are historical evidence only.

### Current source-backed additions are insufficient

At the current zero-dimensional anchor, the one legal energetic-response scale is
`1.526120494`, replaying the 1223.163 nm calibration depth to numerical precision.
This is not yet a frozen profile calibration.

The strongest bounded low-energy test moved an intentionally generous 20% of total ion flux to
15 eV `CF+`. Published activation plus direct ion deposition reduced the predicted density-sweep
gain:

- baseline current 0-D gain: 999.02 nm;
- with 20% low-energy `CF+`: 774.96 nm;
- experiment: 451.47 nm.

Thus the new channel supplies only 40.9% of the required reduction, even at a synthetic fraction
far larger than the present transport operator produces. A full profile rerun is not justified
by this mechanism alone.

The separate collisionless virtual-sheath audit produced no 5--70 eV population and at most
6.68% density-dependent yield flattening. It is also insufficient.

### Inverting the experimental curve diagnoses the missing response

Using the three 200 nm density-sweep depths as development data, the current model requires
effective ion-response multipliers of:

| Electron density (m^-3) | Required multiplier | Shape-equivalent mass if low point is 40 amu |
| --- | ---: | ---: |
| 1.1e15 | 1.2889 | 40.0 amu |
| 1.9e15 | 1.0543 | 59.8 amu |
| 3.1e15 | 0.7367 | 122.5 amu |

The corresponding effective response scales approximately as
`Gamma_effective proportional to n_e^0.463`, rather than the declared linear proxy.
This is an inversion of the scored data, not a model or validation result. It shows that the
missing boundary/surface feedback has ample leverage, but does not identify its cause.

## First-principles interpretation

Jeong et al. measured electron density and self-bias. They explicitly could not use an RFEA
because fluorocarbon film contaminated the analyzer. Their `Gamma_ion proportional to n_e`
axis assumes the positive-ion Bohm speeds are similar despite possible masses from Ar to
`C3F5+`. Figure 6 reports selected *volume densities* from a global model, not wall fluxes.

That distinction matters. A feature code needs flux-weighted species and joint energy-angle
distributions at a physical reference plane. Converting selected volume densities to
`n v_thermal / 4` does not reproduce a global model's wall-loss solution, and representing all
positive ions as Ar removes fluorocarbon-ion fragmentation and polymer-delivery physics.

State-of-the-art HPEM/MCFPM calculations supply the feature model with reactor-computed fluxes
and IEADs, neutralize ions at their first surface collision, and continue their hot-neutral
partners through multiple collisions. They also show that the fraction of fluorocarbon ions can
materially alter deep-feature passivation and etch rate. Petch can now consume that boundary;
the Jeong publication does not expose it.

## Decision and next step

Do not fit the held-out Jeong density sweep and do not spend another long profile matrix on the
current guessed boundary. Retain the new engine capabilities and proceed in either of two ways:

1. obtain species-resolved flux/IEAD output for the Jeong discharge from a validated reactor
   model or diagnostics, then rerun the frozen validation; or
2. continue experimental validation on benchmarks whose source boundary is sufficiently
   specified, while treating Jeong as an open reactor-boundary target.

The two failed 90/180-step current-operator anchor pilots were conservative-remap preflights:
the source-faithful operator etched too far per coarse profile step, so the engine refused at
steps 3 and 6. They do not indicate missing equilibrium or solver failure. A future earned
profile run needs at least a finer schedule or adaptive retry, followed by timestep refinement.

## Primary sources

- Jeong et al. 2023:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10222222/
- Huang and Kushner 2019:
  https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
- Li et al. 2004 ion composition/IED measurements:
  https://cpseg.eecs.umich.edu/pub/articles/jvsta_22_500_2004.pdf
- Kim et al. 2021 capacitively coupled C4F8/Ar ion mass/energy measurements:
  https://doi.org/10.3390/coatings11080993
