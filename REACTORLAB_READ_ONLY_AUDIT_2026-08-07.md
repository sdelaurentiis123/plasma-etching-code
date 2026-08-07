# ReactorLab read-only audit — 2026-08-07

## Decision

**Do not integrate ReactorLab as a reactor boundary provider.**

The repository is useful as a map to primary literature and as an example of
reaction-ledger and serialization structure. Its present numerical calibration,
uncertainty claims, low-pressure branch handling, sheath reduction, and
feature-scale boundary are not strong enough for the validation-first
multiphysics certification standard in this repository.

No ReactorLab source was copied, modified, committed, or pushed. Any future use
of a paper named by ReactorLab must start from the primary paper and an
independent extraction in `research_sources/`.

## Frozen source identity

- Repository: `git@github-personal:seanfoleydesigns/reactorlab.git`
- Branch: `main`
- Commit: `92c87be2d0871c8756dd9a3193ac5a1111f78841`
- Tree: `cc6435cfa8dd5ae274698af26f841e639ba0be1e`
- Tracked worktree state after audit: clean

This identity is recorded so the verdict applies only to the audited state.

## What was checked

The audit read the project overview, validity statement, physics questions,
decisions, and public interface; inspected the global balance solver, chemistry
loader and rate laws, transport and geometry, calibration objective and
published result, sheath/IEDF/IADF reduction, output schema, petch adapter, and
tests; rendered and inspected the published calibration plots; and executed the
test suite locally without editing tracked files.

The fast test selection produced **171 passes, 1 failure, 37 deselections**.
The complete suite independently reproduced the same result with
**208 passes, 1 failure** in 23 minutes 22 seconds. The failure was
`TestLiebermanLimits.test_te_falls_with_pressure`: at 3 mTorr the reported
steady state sits on the imposed 30 eV electron-temperature ceiling and
electron-density floor while the ODE and Newton routes disagree by
approximately `3.58e9 %`. This is a plasma-off/ignition-branch or convergence
classification failure, not a small tolerance miss.

## Calibration audit

The published calibration result has SHA-256
`a4917916c507798baf6d2f6106d7dc9e74d685e53892d4c88b52066ada39fd80`.
The documentation describes five calibrated parameters, while the executable
parameter table contains six:

1. absorbed-power fraction,
2. fluorine wall loss,
3. oxygen wall loss,
4. dissociation scale,
5. ion-ion recombination scale,
6. attachment scale.

Those six parameters are fit to 15 scalar values from only five unique
conditions, with three observables per condition. Two neutral-flux observables
are model-assisted rather than direct. The optimum places the ion-ion
recombination scale exactly on its lower bound.

Recalculation from the published central values gives:

| Quantity | Calibration-set MAPE |
|---|---:|
| Direct measured ion flux | 54.84% |
| Model-assisted F flux | 13.12% |
| Model-assisted O flux | 62.84% |
| All 15 scalar values | 43.60% |

For the nominal held-out condition, predicted/measured central ratios are
`0.7131` for F, `0.7077` for O, and `0.4755` for ion flux, or 36.79% MAPE.
The corresponding reported 95% geometric span factors are approximately
22.65, 216.4, and 20.27. These bands are too broad to certify a quantitative
knobs-to-flux boundary. They also come from a local covariance/pseudoinverse
construction that extends fitted variables beyond their declared physical
bounds; they are not posterior-predictive intervals.

The held-out target is excluded from the scalar objective, but it is computed
and retained during every Latin-hypercube evaluation. This is useful development
monitoring, not an operationally blind test. The raw sample trace is not present
in the audited published tree.

An independent pure-SF6 comparison is also not closed: the reported ion-flux
prediction/measurement ratios are approximately 0.516, 0.175, and 0.246 at
5, 25, and 75 mTorr.

## Physics blockers

- Absorbed power is a fitted fraction of generator power, and its optimum is
  about 0.996 rather than a measured plasma absorption boundary.
- Reactor dimensions are assumed for the Lam TCP 9400 validation rather than
  carried as measured machine geometry with uncertainty.
- The 5% argon actinometer used in the experiment is omitted from the modeled
  gas composition.
- Twenty-nine cross-family ion-ion recombination channels reuse one Pateau
  coefficient, then apply a fitted global scale. This is not species-resolved
  reaction physics.
- Several oxygen excitation reactions terminate in a `lost` pseudo-species,
  creating a bulk atomic sink instead of explicitly returning or quenching the
  excited atom. The stoichiometric bookkeeping is visible, but the tracked
  physical system is not atom-conservative.
- All positive ions are reduced to one flux-weighted mean mass and one
  parametric IEDF/IADF. The IADF includes an engineering transverse-heating
  fraction of 0.1.
- The implemented sheath remains collisionless when its own mean-free-path
  diagnostic enters the collisional regime; it emits a warning rather than
  changing models or refusing the prediction.
- The feature adapter reduces the reactor result to three high-level flux
  channels, one effective ion, one parametric IEDF, and one angular width. It
  cannot carry a species-resolved atomic boundary with covariance and
  applicability gates.

These limitations do not make the project valueless. They mean its current
outputs must not be presented as a validated reactor-to-feature prediction or
used to repair an absolute-depth discrepancy.

## Useful parts, with strict boundary

The reaction tables point to primary sources worth obtaining independently:

- Kokkoris et al., *J. Phys. D* **42**, 055209 (2009),
  DOI `10.1088/0022-3727/42/5/055209`;
- Lallement et al., *Plasma Sources Sci. Technol.* **18**, 025001 (2009),
  DOI `10.1088/0963-0252/18/2/025001`;
- Pateau et al., *J. Vac. Sci. Technol. A* **32** (2014),
  DOI `10.1116/1.4853675`;
- Panagopoulos and Economou, *J. Appl. Phys.* **85**, 3435–3443 (1999);
- the Belen/Gomez reactor and etch datasets identified in the project notes.

The native implementation may reuse concepts—explicit stoichiometric ledgers,
source-bound reaction metadata, deterministic serialization—but not ReactorLab
code, fitted constants, extracted text, or validation claims. Each adopted
mechanism must be rebuilt from a visually checked primary source, conserve
elements and charge over the declared system, expose its validity domain, and
pass an independent held-out gate before it can supply feature-scale depth.

## Consequence for the depth program

ReactorLab does not close the Krüger depth because it models a different
chemistry and does not provide the missing species-resolved fluorocarbon
boundary for that reactor. The correct path remains:

1. keep the published-input Krüger miss visible;
2. infer or measure reactor boundaries only from independent reactor-scale
   observables;
3. carry species, energy, angle, covariance, and regime validity into the
   feature solver;
4. validate absolute depth on multiple chemistries without retuning the
   surface laws to each trench.

This audit therefore changes no existing certification result. It prevents an
insufficiently identified reactor fit from being used to manufacture depth
agreement.
