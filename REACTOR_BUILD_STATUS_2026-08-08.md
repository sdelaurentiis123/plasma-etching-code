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
| neutral and wall transport | **implemented with applicability gates** | exact cylindrical Robin roots, Chapman--Enskog Cl-in-Cl2 diffusion, and a bounded differentiable coverage response replace the former unphysical raw-ratio exponential extrapolation; ratio/temperature transfer remains sensitivity evidence |
| electron-energy ledger | **two-term density-coupled operator implemented; model incomplete** | fixed event losses, attachment energy moments, exact elastic transfer, and the isotropic Hagelaar--Pitchford electron--electron Landau term are distinct; e-ion, anisotropic Coulomb momentum, excited states, and several collision channels remain open |
| direct molecular-Cl2 swarm board | **implemented from primary measurement** | 52 pure-Cl2 mean-arrival-time drift, effective-ionization, and spatiotemporal longitudinal-diffusion markers validate a collision/EEDF solver only |
| Lam Alliance state boards | **implemented as held-out observables** | 62 electron-temperature, 27 volume-average electron-density, and 38 chlorine-dissociation markers; these do not measure absorbed power or wafer flux |
| user collision-deck ingestion | **implemented, rights-safe, and hash-gated** | parses local BOLSIG/LXCat process arrays without packaging third-party bytes; multi-term readiness fails closed unless an elastic angular closure is declared |
| transport-definition-safe swarm grade | **implemented** | refuses flux mobility in place of pulsed-Townsend mean-arrival-time drift and scalar diffusion in place of spatiotemporal longitudinal diffusion |
| deterministic nonconservative swarm/EEDF solve | **two-term temporal-growth rung implemented; multi-term board open** | the conservative eigen-root and density-dependent isotropic e-e fixed point pass manufactured and independent numerical gates, but mean-arrival-time drift, longitudinal diffusion, and steady-state Townsend still require the preregistered multi-term hierarchy |
| absorbed-power-to-EEPF closure | **six-condition atomic-Cl board closes; prediction gate open** | coupled eight-equation reference maps absorbed-power sensitivity to EEPF, six species, and mass-resolved axial ion fluxes; the obsolete molecular-only deck now fails explicitly because atomic-Cl electron collisions are mandatory |
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

The repository-wide suite passes **1,748 passed, 1 skipped** in 12m23s after
the Eq.-11, state-dependent transport-density, bounded-wall, and regenerated-
receipt changes. The swarm board, collision parser, and strict grade carry
independent replay, hash, unit, count, definition, corruption, angular-
closure, and evidence-boundary tests. The repository contains **132 source
entries**. ReactorLab was audited
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

The first version of that board compared Malyshev's measured Eq. 11 quantity
to the wrong model observable, `nCl2/(nCl2+nCl)`, and incorrectly constrained
the powered neutral particle density to the gauge density. The corrected
closure enforces `nCl2+nCl/2=nCl2^0`, evaluates neutral diffusion and ion
mobility at the actual powered particle density, and compares only
`100*nCl2/(nCl2+nCl/2)` to the measurement. It also uses the paper's printed
333 K initial state and source-parameterized Chapman--Enskog diffusion.

At the held-out 30%-absorbed, 500 W point, the Hamilton/atomic-Cl replay gives
`E/N=880.5 Td`, electron-density error `-15.9%`, non-equivalent energy-proxy
error `-1.3%`, and Eq.-11 Cl2 error `-4.8` percentage points (`-6.6%`
relative), which passes the source's about `+/-25%` absolute-density accuracy
band. Adding isotropic e-e collisions gives `865.7 Td`, `-19.2%`, `+4.6%`,
and `-5.2` points; it remains a negative density sensitivity. The corrected
non-Coulomb board is converged from 415 to 813 cells with a largest physical
change of `0.0348%`.

## First reactor-diagnostic-conditioned equipment transfer

A single constant source-to-plasma fraction was inferred from only the 300 W
Malyshev volume-average electron-density marker over the pre-existing
`0.30--0.50` sensitivity bracket. The deterministic inverse solve converged to
`0.3646468`, or `109.39 W` effective absorbed power. With that fraction frozen,
the untouched 500 W density is predicted within `+1.76%` and the global axial
positive-ion flux is `5.779e19 m^-2 s^-1`.

This is a successful power-scaling transfer but not a direct absorbed-power
measurement: the paper reports forward power into the matching network and no
density uncertainty. The independent 500 W Eq.-11 Cl2 residual is `-8.96`
percentage points, or `-12.28%` relative, and passes the source's about
`+/-25%` absolute-density accuracy statement. That accuracy is not relabeled
as a statistical sigma. The mean-energy temperature proxy is `+4.98%`, but it
is not the OES forward observable and receives no validation label. The
bounded wall ratio/temperature transfer, unmeasured absorbed power, omitted
5% rare-gas collision/ion channels, and lack of a local wafer boundary keep
reactor-prediction, wafer-flux, and depth support false. No temperature,
dissociation, held-out measurement, or feature depth selected the fraction.

## Shortest defensible path to depth prediction

1. Retain the hash-gated user-downloaded SIGLO molecular-Cl2 deck and evaluated
   Hamilton/atomic-Cl upgrade; pursue native Kawaguchi arrays separately and
   keep plot digitization as a declared last-resort sensitivity.
2. Complete and independently verify a deterministic, differentiable,
   nonconservative multi-term Boltzmann solver with density-gradient and
   steady-state-Townsend modes, then preregister a no-fit grade against all 52
   direct swarm markers.
3. Land only the collision channels that pass, complete their electron-energy
   moments and uncertainty/support gates, and solve `Te` from independently
   measured or validated absorbed power.
4. Add the declared 5% diagnostic-gas collision/ion channels and a real OES
   forward observable, then extend the Lam grade across power, pressure, and
   both gaps without using feature depth.
5. Add equipment-specific charged transport, sheath, species-resolved IED/IAD,
   neutral/radical delivery, and covariance to the feature-plane boundary.
6. Predict depths on multiple held-out chemistries/reactors. Krueger may be one
   test only if its missing boundary is independently measured or recovered;
   it must never calibrate the reactor chain.

That sequence is slower than multiplying an unknown flux by a fitted scale,
but it is the shortest route whose final nanometers mean anything physically.
