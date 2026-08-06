# E8: thermalized cascade weight as a radical source — built, gated, and bounded

Commit `550e2a5`. Suite 1180 passed / 1 skipped. Box spend $0 (no run bought —
the forecast refused it).

## What was built

The reflection cascade previously *discarded* the weight that drops below the
Eq-2.34 cutoffs, recording it only as a scalar `thermalized_rate` diagnostic.
Huang says that weight is a radical source (sec. 6.4.3,
`research_sources/thesis_extracts/huang_thesis.txt` L5714-5727, verbatim):

> "After losing energy through several collisions with the sidewalls and etch
> front, these energetic species become thermal CFx and CxFy radicals, which
> can passivate the oxide surface or deposit as polymer. In the base case, the
> majority of the CFx and CxFy radicals at low AR (< 5) ... originate from the
> thermal neutrals incident into the feature from the plasma. As the AR
> increases to greater than 10, the neutralized and thermalized CFx+ and CxFy+
> ions become the main source (> 95%) of radicals reaching the etch front."

and fixes the species rule by exclusion (same section):

> "The neutral and thermalized partners of other ions are non-reactive species
> and diffuse out of the feature with no surface reactions (only scattering at
> the surface). There is little surface passivation or polymer deposition
> resulting from those thermalized species."

The cascade now accumulates thermalized weight **per face**, and the gather
returns it to the neutral ledger through an explicit
`thermalized_radical_return` mapping. Escaped weight (continues but leaves the
feature) is not returned.

**The reactive fraction is a declared caller input, never inferred.** Krüger
publishes only an aggregate positive-ion flux with a combined IEAD — his Table
6.1 lists neutrals only (C3F4, C2F3, CF, CF2, CF3, O), and the paper says "The
combined IEAD of all positive ion species". The CFx+/Ar+ split does not exist
for this reactor, so it is a swept input, not a constant. Huang's own base case
is 30% for *his* reactor (L5737) and is not transplanted.

## Gate results

| gate | requirement | measured | verdict |
|---|---|---|---|
| 1a per-face ledger | == scalar diagnostic | 1e-12 | **PASS** |
| 1b weight accounting | nothing vanishes | bound holds, all faces >= 0 | **PASS** |
| 1c deposition | faces only, finite | 6/6 faces, finite | **PASS** |
| 1d locality | normal incidence thermalizes on the struck face | 100% on face 4 | **PASS** |
| 3 blanket 0-D | bitwise unchanged | default off — no result moves | **PASS** |
| 4 Gray N1 beam gates | stay in band | unchanged | **PASS** |
| 5 suite | green | 1180 / 1 skipped | **PASS** |
| **2 floor composition** | **>90% of floor radicals thermalized at AR>=10** | **1.46% at AR 12** | **FAIL** |

## Why gate 2 fails — three measured reasons, none a bug

Frozen straight-trench gathers, Krüger cell (0.09 um opening, 0.85 um mask),
dx = 10 nm, floor band at the etch front:

| oxide AR | E8 source (m^-2 s^-1) | plasma neutrals | E8 share |
|---|---|---|---|
| 2 | 6.792e18 | 5.949e20 | 1.13 % |
| 12 | 6.654e18 | 4.505e20 | 1.46 % |

**(a) The floor's E8 source is pinned at 10% of the direct ion flux by the
leftover rule.** At normal incidence `kress(1) = 1`, so `react = 0.9` and the
continuing weight is exactly `0.1`; the Eq-2.34 angle test fails (incidence is
not beyond 70 deg), so retained energy is zero and the whole 0.1 thermalizes on
the face it struck. Predicted 0.1 x 6.95e19 = 6.95e18 against the measured
6.65e18 — the number is the rule, not a free quantity.

**(b) Local deposition strands the rest on the sidewalls.** Most cascade weight
thermalizes at grazing collisions on the walls, and this implementation deposits
it where it dies. Huang's >95% requires those radicals to *diffuse* after
thermalizing — they carry sticking coefficients of order 0.02-0.1, so they
bounce many times and redistribute toward the floor. That is a thermal
re-emission coupling into the neutral transport solve, not a bookkeeping
change. **This is the E8 completion item**, and it is what the >95% claim rests
on.

**(c) Krüger's geometry is mask-dominated, so his AR threshold does not
transfer.** At oxide AR 12 the total path is 0.85 um mask + 1.08 um oxide over a
0.09 um opening (total AR ~22), and measured plasma-neutral transmission to the
floor is ~12% of the Table 6.1 sum — high enough that plasma neutrals still
supply 98.5% of floor radicals. Huang's >95% is stated for his feature, where
neutral transmission is far lower. For E8 to reach 95% here, plasma delivery
would have to fall to 3.5e17, i.e. ~0.009% transmission — a much deeper feature.
**Transplanting his AR-10 threshold to this geometry was the specification's
error, not the model's.**

## Forecast — and why no run was bought

Coupled 0-D at the measured AR-12 floor delivery (3406 eV front energy), sweeping
the declared fluorocarbon-ion fraction:

| f_FC | floor precursor flux | floor rate |
|---|---|---|
| 0.00 | 4.5047e20 | 2.730 nm/s |
| 0.30 | 4.5247e20 | 2.729 nm/s |
| 1.00 | 4.5713e20 | 2.726 nm/s |

Even at f = 1 (every ion fluorocarbon — the physical upper bound) the floor rate
moves **0.15%, and slightly downward**: the returned carbon is a precursor, so it
also thickens the film. Against ml23's depth undershoot of −49%, E8 with local
deposition is three orders of magnitude short of the gate. Forecast-before-spend
therefore refuses the run; a 60 s endpoint would have re-measured ml23.

## Status

E8 is **implemented, conserved, gated, and default-off** — a correct piece of
bookkeeping that is inert until its transport half exists. The depth undershoot
owner is unchanged and now sharper: it is not the *existence* of a thermalized
radical source but the *delivery* of one, which requires thermal re-emission.

Declared-open, in priority order:

1. **Thermal re-emission of thermalized radicals** (the E8 completion) — the
   only path by which the source reaches the etch front, and the mechanism
   behind Huang's >95%.
2. **The (s0, B0) adsorption pair** (Kwon/Sawin E1). Gray's printed SiO2
   sticking is 0.02 (thesis p.246: "setting s,=0.2 and 0.02 for the cases of
   silicon and SiO, etching respectively"), but it is half of a co-regressed
   pair: landing the scalar alone moves the measured Gray half-rise from 1.94 to
   104.9 against his measured 27 +/- 8, and breaks four validated chemistry
   gates carrying s = 1 in their own closed forms. Left at 1.0 with a gate
   pinning the choice.
3. The CFx+ fraction for this reactor — unpublished; swept, never fitted.
