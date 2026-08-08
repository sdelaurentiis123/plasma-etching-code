# Native reactor-to-feature build status - 2026-08-08

## Bottom line

The native reactor program has advanced from a design document to a
conservation-audited zero-dimensional implementation with executable primary
measurement boards. It is **not yet a predictive knobs-to-wafer-boundary
model**, and it has **not matched the Krueger absolute depth**.

The authoritative Krueger endpoint remains `346.833 nm` simulated versus
`825 nm` measured at 60 s. The old website's apparent depth pass came from
canceling defects and is superseded. Krueger does not publish the
species-resolved positive-ion composition or stable-C4F6 wafer flux needed to
identify that feature boundary. No reactor or surface parameter may be chosen
from the `825 nm` feature target.

## Build ladder

| rung | current state | evidence boundary |
|---|---|---|
| conserved 0-D infrastructure | **implemented and verified** | open-system species, charge, pressure, wall, exhaust, and power bookkeeping; structural correctness is not experimental validation |
| argon source-model rung | **implemented** | source reproduction and transport checks; not a complete equipment-specific predictor |
| chlorine particle balance | **implemented, fail-closed** | solves five heavy-species balances, quasineutrality, and pressure-controller exhaust at supplied `Te`; atom/current ledgers close |
| chlorine chemistry | **partly upgraded** | Lee--Lieberman replay plus evaluated ionization/attachment support, eight state-resolved Hamilton dissociation channels, and separate Kemaneci 36-row forward / 44-row COMSOL nonelastic replays |
| neutral and wall transport | **implemented with applicability gates** | exact cylindrical Robin roots and state-dependent conditioned-wall providers; missing absolute uncertainties and reactor-wall transfer prevent predictive status |
| electron-energy ledger | **operator implemented; model incomplete** | fixed event losses, attachment energy moments, and exact elastic energy transfer are distinct; several collision channels/support tails remain open |
| direct molecular-Cl2 swarm board | **implemented from primary measurement** | 52 pure-Cl2 drift, effective-ionization, and longitudinal-diffusion markers validate a collision/EEDF solver only |
| Lam Alliance state boards | **implemented as held-out observables** | 62 electron-temperature, 27 volume-average electron-density, and 38 chlorine-dissociation markers; these do not measure absorbed power or wafer flux |
| collision arrays and EEDF solve | **open** | no complete native Kawaguchi arrays and no no-fit Boltzmann/Monte Carlo replay against the 52-marker board yet |
| absorbed-power-to-`Te` closure | **open** | generator or forward-minus-reflected power cannot silently become absorbed plasma power |
| charged transport, sheath, IED/IAD, wafer boundary | **open for prediction** | source-model ion mobilities exist, but public evaluated species-resolved transport and a validated equipment closure do not |
| feature coupling | **interfaces exist; predictive chain absent** | nothing currently supports a knobs -> fluxes -> profile/depth claim |

## What the two downloaded papers closed

The official Gonzalez-Magana--de Urquijo 2018 PDF was hash checked, rendered,
and visually audited. Tables A1--A3 contain 52 native pure-Cl2 measurements:

- 23 electron-drift-velocity points;
- 21 effective-ionization points;
- eight density-normalized longitudinal-diffusion points.

The committed transcription preserves the printed values, units, typical
uncertainty ranges, method, pressure and temperature ranges, source pages, and
the paper's own body/table range conflict. The digitizer hash-gates both the
publisher PDF and the exact 300-dpi page renders. Only derived data and audit
metadata are version controlled.

The Kawaguchi et al. 2020 full paper establishes the provenance of its 20-row
collision set. It retains the Gregorio molecular-Cl2 momentum-transfer row
below `1 eV`, modifies it above `1 eV` using integrated differential
measurements, and explicitly applies `2.5x`, `4.0x`, `1.05x`, and `0.85x`
changes to other channel families. That makes it a valuable
**swarm-calibrated candidate provider**, not an untouched fundamental deck.
The article does not supply the native numeric cross-section arrays, so a
no-fit local reproduction remains open. Digitizing its plots would be a
declared lossy fallback.

The unavailable Griffin et al. 1995 full paper is an atomic-Cl elastic
calculation. Its absence keeps the atomic-Cl elastic row calculation-grade and
fail-closed, but it does not block the new direct molecular-Cl2 swarm board.

## Validation state

The complete reactor test surface passes: **209 passed** on 2026-08-08. The
repository-wide suite passes **1,688 passed, 1 skipped** in 12m25s. The new
board alone passes four independent replay, hash, unit, count, and evidence
boundary tests. The repository contains **132 source entries**. ReactorLab was
audited read-only and rejected as a boundary provider at its audited commit;
no ReactorLab source, fit, or claimed validation was integrated.

These results verify implementation invariants and primary-data fidelity. They
do not upgrade an unmeasured reactor boundary into a depth prediction.

## Shortest defensible path to depth prediction

1. Obtain native Kawaguchi cross-section arrays from the authors or a usable
   primary repository; keep plot digitization as a last-resort sensitivity.
2. Implement and independently verify a two-term Boltzmann or Monte Carlo EEDF
   solver, then preregister a no-fit grade against all 52 direct swarm markers.
3. Land only the collision channels that pass, complete their electron-energy
   moments and uncertainty/support gates, and solve `Te` from independently
   measured or validated absorbed power.
4. Grade predicted Lam electron temperature, density, dissociation, and a
   mass-resolved ion-flux board without using feature depth.
5. Add equipment-specific charged transport, sheath, species-resolved IED/IAD,
   neutral/radical delivery, and covariance to the feature-plane boundary.
6. Predict depths on multiple held-out chemistries/reactors. Krueger may be one
   test only if its missing boundary is independently measured or recovered;
   it must never calibrate the reactor chain.

That sequence is slower than multiplying an unknown flux by a fitted scale,
but it is the shortest route whose final nanometers mean anything physically.
