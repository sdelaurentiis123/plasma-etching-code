# Local history reconciliation

Audit date: 2026-07-12. Scope: all local commits after `origin/main` (`49da05c`) through branch commit
`9f419ee`, plus the working changes subsequently committed as `bb65f34`. No remote operation was used.

## Repository topology

- One worktree exists: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`.
- `main` is one linear chain 114 commits ahead of `origin/main` at `109badb`.
- `codex/unified-engine-root-fixes` starts at `109badb` and adds the geometry/numerics corrections through
  `9f419ee`, followed by the architecture reconciliation at `bb65f34`.
- There are no stashes and no second local branch containing competing solver work. The same-timestamp
  commit pairs are parent/child pairs, not forks.
- `git fsck --unreachable --no-reflogs` found one 3,161-byte blob: an earlier Jeon README body. It is not
  executable code or a lost commit.

The two Codex clients therefore wrote serially into one history. There is no alternate branch to merge.
The reconciliation problem is semantic: rapid successive implementations created multiple solver
lineages and allowed evidence from one lineage to be described as evidence for another.

## What changed after the last remote checkpoint

The committed range is approximately 25,599 insertions and 47 deletions across 80 files. It is primarily
additive. No pre-existing test file was changed in the committed range; the additional tests are new.
The pre-existing implementation files materially changed only in:

- `charging_backward.py`: experimental 2-D adjoint/current-balance development;
- `threed.py`: mesh consistency, reinitialization, periodic transport wiring, and neutral-radiosity hooks;
- `__init__.py`: exports for the new contracts.

The chronological work falls into seven causal groups:

1. **2-D charging diagnosis:** open-wafer gates, current pooling, forward/adjoint reciprocity, physical-face
   launch, exit-state weights, support-complete proposals, uncertainty-aware bidirectional estimates, and
   boundary-fitted nodal experiments.
2. **Common plasma boundary:** finite-transit sheath, immutable dimensional species distributions, analytic
   and tabulated densities, and source-independent forward/adjoint adapters.
3. **Compatible electrostatics:** 2-D nodal experiments, then independent Q1 3-D Poisson, triangle-charge
   projection, field transport, replicated current confidence, and safeguarded dielectric current balance.
4. **Common feature evolution:** dimensional surface kinetics, face-resolved energetic events, conservative
   state remap, additive material level sets, CR-2 redistancing, and geometry-dependent charging rebuilds.
5. **Chemistry and data:** versioned interaction tables, checksum-pinned Si-Cl2-Ar+ data, reduced SiO2
   state, Bosch wafer data, Krueger transfer data, and Jeon depth/control extraction with preregistered
   calibration/held-out splits.
6. **Neutral transport:** conservative diffuse form factors, state/material-dependent radiosity, periodic
   transport, and deterministic face gathering.
7. **Jeon numerical ladder and root fixes:** the initial ladder exposed material ownership and subcell
   interface defects; the branch corrected those defects and then invalidated/reset the old ladder rather
   than preserving a favorable obsolete score.

## The three lineages that must not be conflated

| Lineage | Entry point | Purpose | Status |
|---|---|---|---|
| Legacy product/demo | `api.py -> threed.py::run_etch_3d` | Preserve SF6/O2, Bosch, cryo, ALE, ViennaPS parity and historical demonstrations | Runnable; contains calibrated closures and at least one region-specific redeposition suppression; not the common product engine |
| 2-D charging research | `charging_backward.py`, `boundary_transport.py`, `charging_nodal_fixed_point.py` | Derive and falsify adjoint/current-balance estimators cheaply | Experimental; several analytic gates pass, but converged HG/reference/product claims do not close simultaneously |
| Common 3-D engine | `feature_step_3d.py`, `boundary_transport_3d.py`, `charging_coupled_3d.py`, `surface_kinetics.py` | One dimensional boundary/transport/electrostatics/state/interface contract | Authoritative migration target; partial, validity-reporting, and not yet the public product API |

The common 3-D engine may reuse independently gated numerical primitives from `threed.py` (mesh extraction,
velocity extension, level-set advection/redistancing). That does not make `run_etch_3d` the common engine.

## Rollback decision

Do **not** reset to the last remote checkpoint or to an arbitrary commit from yesterday:

- the history is linear, so rollback does not recover a hidden competing branch;
- the range is overwhelmingly additive and contains the dimensional contracts needed for the target;
- the complete current suite passes 277 tests with one skip;
- the branch fixes concrete geometry and conservation defects discovered by the new experimental ladder;
- rolling back would restore stronger-looking Jeon numbers that were generated before those defects were
  understood and would delete the evidence that invalidated them.

Preserve `origin/main` as the immutable pre-campaign reference and preserve every current local commit.
If a subsystem later fails review, revert or replace that bounded subsystem on this branch with a new
local commit; do not rewrite the entire history.

## Reconciliation rules

1. `feature_step_3d` is the only target for new feature-scale governing physics.
2. The legacy path remains compatibility-only until each mechanism is re-earned through the common
   boundary, transport, material-state, remap, validity, and profile contracts.
3. The 2-D charging modules remain research/reference tools. Their fitted or benchmark-emulation source
   laws may not enter production defaults.
4. No benchmark name, aspect-ratio band, expected answer, or surface region may select a governing formula.
5. Unknown reaction probabilities are explicit, provenance-bearing calibration inputs with uncertainty;
   they are not described as first-principles predictions.
6. A benchmark promotes physics only when the experiment supplies enough boundary/material information to
   test that physics. Code parity and published simulations remain separate evidence classes.
7. Every migration uses one causal observable, analytic/numerical invariants, a preregistered external
   metric, and a compute ceiling. Full regression protects compatibility but does not promote evidence.

## Next bounded work

1. Introduce a versioned common-engine case/result schema beside the compatibility API. Do not silently
   redirect `Process.run()`.
2. Close the species/material mass loop: emitted etch products, transport, sticking/reaction, and mask
   evolution through the same material mechanism interface. Do not port the legacy notch-foot suppression.
3. Add sourced ion-reflection interaction tables with energy/angle/material/roughness support and refusal.
4. Re-run a time-resolved profile campaign and held-out geometry/process transfer through the common engine.
5. Enable charging in a product gate only when charging-on/off changes an experimental profile above the
   combined numerical and measurement error budget.

The roadmap target remains the same: a trustworthy forward feature solver first; data-efficient
calibration and held-out transfer next; verified gradients and inverse recipe design after the forward
physics earns those gates.
