# Numerical geometry, multi-fidelity calibration, and AMR plan

Date: 2026-07-17
Status: active plan; bounded topology/remap/refinement/CUDA prerequisites are closing before any new endpoint

Related documents:

- [`VIENNAPS_CODEBASE_AUDIT_2026-07-17.md`](VIENNAPS_CODEBASE_AUDIT_2026-07-17.md)
- [`KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md`](KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md)
- [`PHYSICS_FIRST_UNIFIED_ENGINE_2026-07-17.md`](PHYSICS_FIRST_UNIFIED_ENGINE_2026-07-17.md)

## 1. Three different problems

Three issues have repeatedly been conflated:

| Problem | Question | Correct tool |
| --- | --- | --- |
| Surface correspondence | Where did old surface state go after remeshing? | closest-surface and conservative overlap transfer |
| Spatial resolution | Where are small cells actually needed? | sparse narrow bands and AMR |
| Parameter identification | Which uncertain physical closures fit calibration data? | multi-fidelity trust-region/Bayesian design |

AMR makes one simulation cheaper. Better calibration reduces the number of simulations. Robust
remapping makes both scientifically meaningful.

## 2. Current Krüger remap incident

At 56.482184 s on the 5 nm development run, timestep reduction by 8x left the apparent remap jump
essentially unchanged at about 10.2 nm. The cause was not physical velocity. The 20 nm periodic
lateral dimension had been passed to transport and topology but not to the surface-state neighbor
search. Faces adjacent through the periodic seam therefore looked separated in ordinary Euclidean
coordinates.

The narrow repair now under regression uses:

- exact point-to-triangle distance after KD-tree candidate culling;
- periodic image lookup for both certification and interpolation;
- unchanged material-local conservative and bounded state reconstruction;
- continued refusal for a truly distant surface;
- manufactured tests for retriangulation and a periodic seam.

The mixed-operator endpoint is development-only. An authoritative fine-grid confirmation must start
at `t=0` with one checksum-bound remap operator.

## 3. Production remapping ladder

### R0: present first-order transfer

Material-local inverse-distance interpolation followed by an exact area-integral correction. It is
bounded and conservative, but correspondence is based on neighboring samples rather than geometric
overlap.

### R1: accelerated closest-surface certification

Replace temporary KD-tree candidate culling with a reusable triangle BVH/AABB structure. The gate
measures physical surface distance, including periodic images, not centroid displacement.

### R2: conservative triangle-overlap remapping

Project old/new local patches into a common tangent chart, intersect their polygons, and transfer each
extensive field by overlap area. Intensive coverages are reconstructed with bounded monotone weights.
Uncovered new surface receives the declared newly-exposed state; disappearing old area is assigned to
the material-removal ledger.

### R3: refine/coarsen transfer for AMR

Use the same overlap weights to prolong/restrict surface state across refinement levels. Extensive
inventories conserve; intensive fields remain bounded; material identity never averages across a
junction.

## 4. Remap verification matrix

| Test | Exact requirement |
| --- | --- |
| Identical mesh, zero motion | bitwise-identical state |
| Same surface, different triangulation | identical integrals; convergent local field |
| Periodic seam crossing | wrapped result equals translated interior result |
| Translating plane | exact extensive conservation and declared intensive transport |
| Uniformly receding surface | removed inventory equals geometric ledger |
| Growing surface | retained state rides material; new area uses declared closure |
| Curved sphere/cylinder | refinement order reported |
| Material junction | no cross-material borrowing |
| Discontinuous surface field | no new extrema; L1 convergence under refinement |
| Repeated refine/coarsen cycle | bounded accumulated diffusion and round-trip error |
| True topology change | structured refusal or explicit topology-event handler |
| CPU/GPU implementation | common answer within declared roundoff/error tolerance |

The current global conservation correction alone is not enough. A remap can conserve total inventory
while putting it on the wrong patch; local convergence and junction tests are therefore mandatory.

## 5. AMR is a cost method, not an accuracy waiver

A uniform 2.5 nm, 60 s Krüger run is not the next test. The efficient evidence ladder is:

```text
initial state                          late 56.48 s state
     |                                       |
0.5 s at 20/10/5 nm              frozen-rate + 0.1--0.5 s bursts
     +-------------------+-------------------+
                         |
          locate the grid-sensitive operator/region
                         |
            refine only where the estimator demands
                         |
         AMR must reproduce uniform 5 nm short tests
                         |
             one authoritative long AMR/fine run
```

The late checkpoint is essential: a shallow open feature does not exercise narrow access, corner
curvature, polymer necking, or material-junction remapping.

### Refinement indicators

Candidate indicators are:

- signed-distance band around every material interface;
- curvature and under-resolved gap width;
- material junctions and topology-event proximity;
- gradients in surface coverage, inventory, normal velocity, and charge;
- changes in visibility/access probability;
- Poisson residual and electric-field gradient when charging is active.

### Mesh contract

The first true AMR implementation should be block-structured and 2:1 balanced. It must carry:

- material-local signed-distance fields;
- authoritative material ownership;
- conservative prolongation/restriction of volume and surface inventories;
- periodic neighbor relationships;
- consistent marching-surface extraction across coarse/fine boundaries;
- checkpoint/restart and deterministic mesh refinement decisions;
- CPU/GPU parity and explicit work/error reporting.

Sparse narrow-band storage with one global `dx` should be benchmarked first. It reduces inactive
volume work and is simpler, but it is not AMR. Surface transport/radiosity cost still scales with the
resolved surface mesh, so AMR must also coarsen smooth, low-gradient surface regions to produce a
large end-to-end gain.

### AMR acceptance gates

- exact material and surface ledgers at every refine/coarsen event;
- no new extrema in bounded intensive state;
- identical topology classification to the uniform reference;
- short-burst depth/opening increments within the uniform 5 nm uncertainty;
- surface-velocity area-L1/RMS and worst-feature error reported;
- timestep/refinement schedule invariance;
- lower wall time or memory at a matched error;
- no experimental parameter changes between uniform and AMR runs.

## 6. Parameter calibration is not a rectangular grid search

The current Krüger problem has two fitted closures and two base observables:

```text
crosslinked-film fraction  ----mostly----> mask opening
oxide-yield scale          ----mostly----> etch depth
                       with measured cross-coupling
```

A local response matrix already improves dramatically on a parameter grid:

\[
\Delta \mathbf y \approx J\,\Delta\boldsymbol\theta,
\qquad
\mathbf y=(\text{opening},\text{depth}),
\quad
\boldsymbol\theta=(f_{crosslink},s_{yield}).
\]

The ban on frozen-map root solvers applies to the nonsmooth, sample-biased charging equilibrium. It
does not prohibit a safeguarded response solve for two smooth, full-profile calibration observables.

## 7. Multi-fidelity trust-region calibration

The lasting low-dimensional method should model the fine response as

\[
\mathbf y_f(\boldsymbol\theta)
=\mathbf y_c(\boldsymbol\theta)+\mathbf d(\boldsymbol\theta),
\]

where `c` is a cheap/coarse or reduced operator and `d` is a paired coarse-to-fine discrepancy.

One iteration is:

1. fit/update a local linear or quadratic response from calibration-only runs;
2. solve the bounded step inside a trust region;
3. evaluate the proposed point cheaply;
4. promote only informative points to the fine operator;
5. compare predicted versus observed improvement;
6. accept/shrink/grow the trust region;
7. stop when base residual, parameter uncertainty, numerical uncertainty, and identifiability gates
   pass.

For Krüger, the coarse and fine runs must share the same physical mechanism and differ only in
declared numerics. The final parameters and final base confirmation belong to the authoritative fine
or certified-AMR operator. A coarse result cannot silently inherit fine-grid calibration authority.

### Calibration gates

- only declared calibration observations enter the objective;
- held-out outcomes remain unread and checksum-sealed;
- parameter bounds and transformations are physical;
- response/Jacobian condition and parameter correlation are reported;
- paired numerical/sample uncertainty accompanies every response;
- lack of identifiability causes refusal or a prior/data request, not another free parameter;
- the final point receives an independent full-operator confirmation;
- parameter stability under refinement is reported rather than calibrated away indefinitely.

## 8. When to use other optimization methods

| Situation | Preferred method | Reason |
| --- | --- | --- |
| 2--4 parameters, smooth local response | bounded multi-fidelity trust region/secant | interpretable and sample-efficient |
| 5--15 expensive uncertain parameters | Gaussian-process Bayesian optimization | selects informative simulations and carries uncertainty |
| Many correlated parameters | active subspaces/sensitivity screening first | avoids fitting unidentifiable directions |
| Smooth differentiable PDE subsystem | adjoint gradient | cost weakly dependent on parameter count |
| Ray hits/topology/discrete events | derivative-free local model or likelihood-free method | naive gradients are unreliable |
| Large validated simulation corpus | learned surrogate with exact verification | fast design mode, not authority |

Bayesian optimization is not automatically better for two parameters. Its overhead and prior choices
are unnecessary when a well-conditioned local response already exists. Conversely, a rectangular grid
becomes indefensible once each evaluation costs hours.

## 9. Calibration edge as a product capability

A useful engine should expose calibration as a first-class, evidence-aware operation:

```text
calibration specification
  - parameters, bounds, priors, provenance
  - calibration observables and uncertainty
  - sealed held-out set
  - exact/cheap operators and costs
             |
      experiment selector
             |
  response + discrepancy model
             |
  identifiability/error report
             |
      frozen reveal manifest
```

The deliverable is not only a fitted parameter vector. It includes posterior/trust-region uncertainty,
correlations, numerical discrepancy, sensitivity, applicability, exact file checksums, and the blind
prediction receipt. This is a practical advantage over ad hoc manual tuning even before a learned
foundation model exists.

## 10. Vienna comparison program

ViennaPS is scored as a separate black-box backend on identical overlapping equations. The comparison
has four columns:

| Dimension | Vienna question | Petch question | Superset gate |
| --- | --- | --- | --- |
| Accuracy | Does Vienna match analytic/manufactured truth? | Does petch? | petch no worse within uncertainty |
| Runtime | Cost at matched error | Cost at matched error | petch meets release budget or documents gap |
| Reliability | Does it finish unattended? | Does petch recover/refuse correctly? | no silent loss; comparable completion rate |
| Physics | Which channels exist? | Which are unified and verified? | every overlap plus evidence-backed additions |

Current honest positioning:

- Vienna leads in production geometry, sparse level sets, model/API ergonomics, backend breadth, and CI;
- petch leads architecturally in dimensional boundaries, charging, response lineage, conservative
  inventories, stochastic modes, uncertainty, and calibration/validation contracts;
- neither statement proves better experimental prediction;
- the held-out league decides that claim.

## 11. Current bounded execution order

1. preserve and classify the completed mixed-operator endpoint as development evidence — complete;
2. run the initial paired 10/5 nm short burst — complete; late topology-equivalent evidence remains open;
3. quantify the 10-to-5 nm difference — complete for the initial state: global depth/opening/mask
   thickness are close, while local maximum/top width differs by 71.235%;
4. explicitly select and short-test the indexed/common-refinement remap in the Krüger worker —
   complete at 10 nm: common refinement is the production candidate after a paired two-step gate;
   one bounded 5 nm confirmation remains before the backend is frozen;
5. implement the calibration trust-region/discrepancy artifact without reading held-outs — complete;
   the real readiness receipt blocks a proposal until remap selection and a current-epoch fine anchor;
6. benchmark sparse narrow-band storage — complete with a bounded Krüger-specific no-go: at 900 nm
   depth the best safe 5 nm core/indexed reductions are 1.415x/1.148x and the optimistic work ceiling
   is 1.921x, below the predeclared 3x/2x gates;
7. hold fixed-`dx` sparse evolution and one-level Krüger volume AMR; the existing exact extruded
   diffuse path has now been costed on the real 10 nm base surface and is `7.75%` slower in the paired
   local timing, so it remains a deterministic reference rather than replacing production transport;
8. remove repeated immutable full-3-D event validation/yield work, then revisit AMR only on a larger
   genuinely volume-dominated geometry or after a cost model earns it;
9. perform one clean authoritative fine/AMR base confirmation;
10. freeze parameters, run held-outs once, and compare the overlapping cases with ViennaPS.

The corrected unified-device CUDA preflight takes 7.780 s per warmed 5 nm step. Chemistry/material
routing (29.5%), ballistic transport (28.8%), diffuse exchange (15.3%), and legacy remap (13.7%)
dominate; level-set redistance is only 0.8%. `KRUEGER_2024_CUDA_PROFILE_REPORT_2026-07-17.md`
therefore sets the speed work order. A one-step arithmetic 60 s projection is 5.19 h and is not an
endpoint prediction; no blind long run follows from it.

This sequence attacks both costs: fewer runs and cheaper runs, while keeping the experimental judge
independent.
