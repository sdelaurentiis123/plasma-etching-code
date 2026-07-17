# Surface-state / neutral-radiosity coupling plan

Date: 2026-07-17
Status: bounded diagnostic authorized; engine integration held behind the gates below
Scope: one reusable common-engine repair, not a Krueger-only chemistry retune

Implementation checkpoint, 2026-07-17: SC0 cache replay and the SC1 standalone manufactured operator
now pass. `src/petch/surface_radiosity_coupling_3d.py` contains the immutable cache and transactional
embedded integrator; `22` focused plus adjacent radiosity/material-router tests pass. It is deliberately
not wired into profile evolution. SC2 still requires nested form-factor-ray refinement before SC3 may
begin. Later the same day, the finite-area/replicated estimator, cell-by-cell float64 visibility
authority, exact-periodic certified fast path, and exact-overlap fixed-physical-patch scorer passed
their gates. The real-checkpoint visibility gate is now closed: all `14,664` eight-ray events match
float64 exactly, including paths with `320` wraps. This is numerical infrastructure, not SC2 closure;
one replicated patch-scored shortest-horizon receipt remains open.

## 1. Why this task exists

The current feature step computes direct neutral transport and diffuse radiosity once, then advances
the surface mechanism over the full profile timestep. That is exact only when surface reaction
probabilities remain effectively constant during the step.

The repaired Krueger endpoint checkpoint disproved that assumption. At only one-sixteenth of the
next declared profile step (`0.00779536 s`), the exact material mechanism produced finite oxide
removal and a maximum local reaction-probability change of about `0.0596`. Geometry moved only about
`0.021 dx`, so the first under-resolved coupling is surface state to neutral redistribution—not the
level-set move.

This is a useful engine finding. Polymer coverage changes neutral sticking; sticking changes how
unreacted neutrals bounce through a cavity; the redistributed flux changes polymer and oxide
chemistry again. Holding that feedback fixed for a whole outer step can give the wrong accumulated
etch rate even when transport, chemistry, and geometry each conserve correctly in isolation.

This dependency is not a new petch invention. Ertl and Selberherr's feature-scale formulation states
explicitly that flux distributions depend on surface coverages through effective sticking, making
the coverage/flux problem recursive; they solve it iteratively and limit each surface advance to a
small fraction of a grid spacing ([ECS Transactions 23, 61--68
(2009)](https://www.iue.tuwien.ac.at/pdf/ib_2009/CP2009_Ertl_2.pdf)). A later three-dimensional
fluorocarbon model likewise evolves coverage equations from neutral/ion flux and assumes a
pseudo-steady coverage only when adsorption/desorption is demonstrably faster than surface evolution
([Journal of Computational Electronics 22,
2023](https://doi.org/10.1007/s10825-023-02068-y)). The diagnostic therefore tests the condition under
which the cheap pseudo-steady/lagged approximation is legal; it does not add an exotic mechanism.

## 2. First-principles decomposition

On one fixed surface mesh, for each diffuse neutral species,

```text
q(z) = d + F^T diag(1 - p(z)) q(z)
dz/dt = f(z, q(z))
```

where:

- `d` is direct boundary-to-surface flux;
- `F` is the geometry-only diffuse form-factor operator;
- `z` is conservative surface state;
- `p(z)` is the state-dependent reaction probability;
- `q(z)` is total incident flux after diffuse redistribution;
- `f` is the existing exact material mechanism.

While geometry and the boundary are fixed, caching `d` and `F` is physically exact. Freezing `q` is
not. Therefore the lasting repair is to trace direct particles once, estimate the form factors once,
and cheaply re-solve radiosity as chemistry evolves.

## 3. Stable large task and adaptive subtasks

The large task is stable: **make time-dependent surface chemistry and diffuse neutral transport one
converged operator inside the common feature engine.** The subtasks may change when evidence lands.

| Subtask | Cheapest adequate experiment | Pass gate | Failure branch |
| --- | --- | --- | --- |
| SC0 cache identity | same checkpoint/state through ordinary and cached paths | facewise flux agreement for every species; full mesh/config/operator fingerprint match | repair cache construction; no chemistry run |
| SC1 manufactured DAE | one-face, state-independent, two-face cavity, and nonmonotone-sticking cases | embedded step-halving convergence, exact accepted ledgers, rejected trials contribute zero | repair integrator/error norm |
| SC2 bounded Krueger diagnostic | fixed R1.9 geometry; horizon ladder `dt/16` to `dt`; R17/R19 mechanisms | direction survives coupling-step and nested ray refinement; every state/inventory/displacement gate passes | classify coupling, sampling, or stale-geometry limit |
| SC3 common-engine module | explicit opt-in branch in `advance_feature_step_3d` | bitwise zero-step/no-coupling parity; deterministic replay; one direct trace and one form-factor build per outer step | keep diagnostic-only |
| SC4 short profile truth | current uniform 10/5 nm, about `0.5 s` | timestep-refined profile/ledger agreement and expected chemistry direction | repair operator before AMR/endpoints |
| SC5 acceleration | profile call counts and wall-time only after SC4 | same accepted trajectory within tolerance; measured speedup | retain simple implementation |

No subtask is allowed to loosen the calibration or held-out contract. A diagnostic failure closes an
implementation hypothesis; it does not authorize a new fitted parameter.

### Current earned boundary

- SC0: the repaired periodic production call and cached replay agree face-by-face and hash-by-hash on
  all `1,833` checkpoint faces; direct transport and form factors are each built once.
- SC1: manufactured open-surface, constant-probability, reflective-cavity, nonmonotone-probability,
  deterministic replay, identity refusal, minimum-step refusal, incomplete routing, product refusal,
  and zero-duration gates pass. A `4x` accepted-substep refinement contracts the measured cavity error
  by about `5.15x` against the finer reference.
- SC2 partial: all 8/16/32-ray shortest-horizon runs pass and preserve the R17/R19 integrated oxide
  direction. Full-state receipts show why this is not yet engine-integration authority: 16 -> 32
  changes integrated oxide removal by at most `0.024%` and the paired effect by `0.0066%`, but the
  maximum local fluorocarbon deposition/growth L1 change is about `18%` and its worst face about
  `58%`. `FORM_FACTOR_ACCURACY_AUDIT_2026-07-17.md` defines FF0--FF5: independent replicated
  scrambled-QMC uncertainty, fixed-physical-patch gates, visibility/reciprocity certification, and
  adaptive/hierarchical refinement before SC3 wiring. The common foundations for the first three now
  pass manufactured tests; no real-checkpoint uncertainty or production ray count is claimed yet.

The mechanism API still does not expose neutral-species-to-material stoichiometry. The module therefore
certifies the radiosity projectile ledger and all four material exchange ledgers independently. It does
not manufacture a combined cross-ledger claim.

## 4. Required integrator

Endpoint reaction-probability change is a safety limiter, not a numerical error estimate. A small
change can be amplified by a reflective cavity, and a large change on a dark face can be irrelevant.
It can also miss nonmonotone motion.

Use embedded step-doubling:

```text
at state z_n and candidate step h
    solve q(z_n)
    coarse: advance chemistry once by h

    fine:   advance chemistry by h/2
            solve q at the midpoint state
            advance chemistry by h/2

    compare coarse and fine conservative outputs
    if all norms and safety gates pass:
        accept the fine path
    else:
        reject without touching state/ledgers and retry at h/2
```

Do not Richardson-extrapolate the state. The fine mechanism result preserves its declared bounds,
positivity, and conservation; an extrapolated state may not.

The error comparison must cover:

- every conservative surface-state increment;
- removed, outgoing, unresolved, and deposited inventories for every named inventory;
- per-face time-integrated recession and growth;
- area-integrated oxide and mask observables;
- radiosity source/reacted/escaped balance;
- exact material ledgers;
- the outer frozen-geometry displacement limit.

Local `max |Δp|` remains a hard cap. It is never the sole acceptance rule.

## 5. Cache contract

The cached operator is valid only for one exact geometry/boundary/operator epoch. Its identity binds:

- vertices, faces, areas, normals, active-face map, and material IDs;
- periodic domain and seam convention;
- direct boundary state and sampling configuration;
- form-factor seed, ray count, offsets, and visibility operator;
- engine/operator epoch.

Cache certification is facewise for every neutral species, not merely an integrated flux comparison.
Perturbing any bound input must invalidate or refuse the cache. Instrumentation must prove one direct
transport evaluation and one form-factor construction per accepted fixed-geometry outer step.

## 6. Conservation and declared limits

Every accepted radiosity solve retains source = reacted + escaped. Every accepted material step
retains removed = outgoing + unresolved, by inventory and face. Rejected trials contribute exactly
zero. Final-minus-initial stored inventory must reconcile with the accumulated exchange on the frozen
surface.

The first engine version refuses:

- product-emitting mechanisms whose products require transport/redeposition;
- incomplete active-face routing;
- a cache whose geometry or operator identity changed;
- a chemistry step below the declared minimum timestep;
- a frozen horizon whose predicted gross surface displacement exceeds `0.05 dx`.

Species-resolved neutral consumption is not yet exposed by the mechanism API. Until it is, the engine
may claim exact radiosity closure and exact material closure separately, but not a combined
neutral-to-material cross-ledger.

## 7. Runtime discipline

The bounded diagnostic has a hard wall-time limit and stops at the first failed horizon. A timeout
means the current controller or cold radiosity implementation is inefficient; it does not prove the
physical coupling impossible.

Only after correctness passes may implementation work reuse sparse matrix structure, preconditioner
state, or GMRES warm starts. Only after the short profile truth test passes may another long endpoint
run be purchased. No GPU and no held-out observation are needed for SC0--SC3.

## 8. Promotion statement

Successful SC0--SC2 establishes a causal direction and an implementation contract. It does not by
itself validate Krueger or permit a `60 s` rerun. SC3 plus SC4 are the authority gates. The first
held-out claim remains downstream of a frozen, timestep/grid/sample-refined base calibration on the
same final operator epoch.

## 9. 2026-07-17 checkpoint decision

The first sealed replicated checkpoint screen reached the coupling operator and completed under its
five-minute authority limit. It passes exact material/radiosity closure and the integrated profile
directions, but it is not promoted: `dt_next/16` permits roughly `0.73 dx` gross local motion and
20/40 nm fluorocarbon/radiosity patches remain sampling dominated. Stage B and live feature-step
wiring are held. The earned successor is a shorter common frozen horizon plus selected-source nested
ray allocation using the persisted direct/form-factor artifacts. No additional ballistic transport,
global ray escalation, moving profile, or held-out observation is needed for that diagnosis.
