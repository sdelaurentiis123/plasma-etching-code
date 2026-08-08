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
| electron-energy ledger | **two-term density-coupled operator implemented; model incomplete** | fixed event losses, attachment energy moments, exact elastic transfer, and the isotropic Hagelaar--Pitchford electron--electron Landau term are distinct; e-ion, anisotropic Coulomb momentum, excited states, and several collision channels remain open |
| direct molecular-Cl2 swarm board | **implemented from primary measurement** | 52 pure-Cl2 mean-arrival-time drift, effective-ionization, and spatiotemporal longitudinal-diffusion markers validate a collision/EEDF solver only |
| Lam Alliance state boards | **implemented as held-out observables** | 62 electron-temperature, 27 volume-average electron-density, and 38 chlorine-dissociation markers; these do not measure absorbed power or wafer flux |
| user collision-deck ingestion | **implemented, rights-safe, and hash-gated** | parses local BOLSIG/LXCat process arrays without packaging third-party bytes; multi-term readiness fails closed unless an elastic angular closure is declared |
| transport-definition-safe swarm grade | **implemented** | refuses flux mobility in place of pulsed-Townsend mean-arrival-time drift and scalar diffusion in place of spatiotemporal longitudinal diffusion |
| deterministic nonconservative swarm/EEDF solve | **two-term temporal-growth rung implemented; multi-term board open** | the conservative eigen-root and density-dependent isotropic e-e fixed point pass manufactured and independent numerical gates, but mean-arrival-time drift, longitudinal diffusion, and steady-state Townsend still require the preregistered multi-term hierarchy |
| absorbed-power-to-EEPF closure | **six-condition physical sensitivity board runs; prediction gate open** | coupled eight-equation reference maps an explicit absorbed-power sensitivity to EEPF, six species, and mass-resolved axial ion fluxes; the complete Hamilton/atomic-Cl/e-e board finishes in 98 s, while outer implicit derivatives, measured absorbed power, and complete excited-state chemistry remain open |
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

The complete rung-by-rung execution and claim plan is recorded in
`FULL_STACK_EXECUTION_PLAN_2026-08-08.md`.

## Validation state

The current reactor and Malyshev test surface passes: **277 passed** on
2026-08-08. The
repository-wide suite passes **1,696 passed, 1 skipped** in 16m25s. A final
source-definition/angular-closure correction made after that run collected was
separately replayed on the exact affected surface: **12 passed**. The swarm
board, collision parser, and strict grade carry independent replay, hash, unit,
count, definition, corruption, angular-closure, and evidence-boundary tests.
The repository contains **132 source entries**. ReactorLab was audited
read-only and rejected as a boundary provider at its audited commit; no
ReactorLab source, fit, or claimed validation was integrated.

These results verify implementation invariants and primary-data fidelity. They
do not upgrade an unmeasured reactor boundary into a depth prediction.

## Latest density-coupling rung

The exact isotropic electron--electron Fokker--Planck term from
Hagelaar--Pitchford equations 34--41 is now coupled to temporal electron
growth through a deterministic fixed point. Its piecewise-constant moments
have analytic batch JVP/VJP operators; the discrete operator conserves
particles to roundoff and its energy defect falls below `0.3%` under the
declared refinement test. Solve-local continuation reduced a representative
coupled eigen solve from `479` growth-root evaluations to `11` without
changing the converged state outside its `1e-7` tolerance. The complete
six-condition Lam sensitivity now takes `98 s`; the held-out 30%-absorbed,
500-W point takes `35.7 s` instead of running beyond 13 minutes.

This new physics does **not** close the Lam state board. At the held-out point,
the previous Hamilton/atomic-Cl replay gave `E/N=1091.2 Td`, electron-density
error `-15.1%`, energy-proxy error `+1.3%`, and Cl2 proxy error `-16.9`
percentage points. Adding isotropic e-e collisions gives `1077.9 Td`,
`-18.5%`, `+7.2%`, and `-17.1` points. Across the full sensitivity it lowers
the axial ion flux by roughly 1--3%. It is therefore retained as a physical
density-coupling sensitivity, not selected as a correction and not promoted
to wafer-flux or feature-depth support.

## First reactor-diagnostic-conditioned equipment transfer

A single constant source-to-plasma fraction was inferred from only the 300 W
Malyshev volume-average electron-density marker over the pre-existing
`0.30--0.50` sensitivity bracket. The deterministic inverse solve converged to
`0.3571645`, or `107.15 W` effective absorbed power. With that fraction frozen,
the untouched 500 W density is predicted within `-1.01%` and the global axial
positive-ion flux is `5.708e19 m^-2 s^-1`.

This is a successful power-scaling transfer but not a direct absorbed-power
measurement: the paper reports forward power into the matching network and no
density uncertainty. The independent 500 W residuals remain `+7.70%` for the
mean-energy temperature proxy and `-18.71` percentage points for the Cl2
proxy. The result therefore separates the next work cleanly: the constant
equipment-power scale is adequate for the density trend, while chemistry,
wall state, diagnostic forward models, and spatial wafer transfer remain
open. No temperature, dissociation, held-out measurement, or feature depth
selected the fraction, and every wafer/depth support flag remains false.

## Shortest defensible path to depth prediction

1. Load and hash a user-downloaded Gregorio/SIGLO molecular-Cl2 deck from
   LXCat; pursue native Kawaguchi arrays separately and keep plot digitization
   as a declared last-resort sensitivity.
2. Complete and independently verify a deterministic, differentiable,
   nonconservative multi-term Boltzmann solver with density-gradient and
   steady-state-Townsend modes, then preregister a no-fit grade against all 52
   direct swarm markers.
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
