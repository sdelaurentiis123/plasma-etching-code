# Differentiable Transport Research — 2026-07-21

**Subject:** gradients through the exact deterministic 2-D diffuse-exchange operator
(`src/petch/deterministic_exchange_2d.py`, method `analytic_occlusion`) and onward through the
radiosity solve, surface chemistry, and level-set evolution — feasibility, calculus, literature
map, staged build plan, framework choice.

**Status:** research memo, not committed. Companion to the `Differentiability` row of
`END_STATE_VERIFICATION.md` (currently **open**: "hard hit/escape indicators still move
discontinuously with field and geometry for a finite sample set … the gradient gate requires a
verified discontinuity treatment"). Today's operator removes exactly that blocker for the 2-D
neutral-transport lane: the hit/miss indicator has been replaced by closed-form pieces with
analytically enumerated seams.

---

## 0. Executive verdict

The `analytic_occlusion` operator is, in the precise sense used by the differentiable-rendering
literature, an operator whose **boundary (silhouette) integrals are computed analytically instead
of sampled**. The three known families of discontinuity treatments — edge sampling
([Li, Aittala, Durand, Lehtinen 2018](https://dl.acm.org/doi/10.1145/3272127.3275109) /
[redner](https://github.com/BachiLi/redner)), reparameterization
([Loubet, Holzschuch, Jakob 2019](https://dl.acm.org/doi/10.1145/3355089.3356510)), and
warped-area sampling ([Bangaru, Li, Durand 2020](https://dl.acm.org/doi/10.1145/3414685.3417833))
— exist **because** those authors could not enumerate visibility boundaries in general 3-D scenes
and had to estimate boundary terms by Monte Carlo. In our 2-D extruded setting the boundaries
*are* enumerated (projective images of blocker endpoints; linear-in-source collinearity events),
so the Reynolds/Leibniz boundary terms can be evaluated (or shown to vanish) in closed form.
The closest true prior art is not any 2018+ paper but
[Arvo 1994, "The Irradiance Jacobian for Partially Occluded Polyhedral Sources"](https://dl.acm.org/doi/10.1145/192161.192250)
— analytic derivatives of irradiance including occlusion-boundary terms — and the compiler
formalization of exactly our situation is
[Teg — Bangaru, Michel et al. 2021, "Systematically Differentiating Parametric Discontinuities"](https://dl.acm.org/doi/10.1145/3450626.3459775).
Our operator is a hand-instantiated Teg-style program.

**Key analytical finding (Section 1):** for this operator the per-source integrand is
*continuous* across every visibility event (shadow intervals are born/die at zero width, and cuts
enter/leave through the target endpoints), so the outer Leibniz boundary terms **vanish
identically** in generic configurations. `H_ij` is C^1 in vertex coordinates away from a
measure-zero set of degenerate configurations (exact collinearity of a blocker with the source
segment, grazing tangency sustained over finite source measure). Differentiating the closed-form
per-interval formulas with frozen panel/branch structure gives the exact gradient — no Dirac
terms, no edge sampling, no reparameterization needed.

**First milestone (Stage A, Section 5):** forward-mode (dual-number) derivative of
`build_deterministic_line_exchange_2d` w.r.t. one vertex coordinate, spike-tested against
Richardson-extrapolated central FD on trench geometries including grazing configs. Estimated
scope: one focused session; the kernels are deliberately scalar float math, which makes a
~100-line dual-number overlay feasible without touching the certified quadrature logic.

**First commercial value (Stage B/C):** tangent-mode (forward) sensitivities of endpoint profile
observables w.r.t. 2–10 chemistry/process calibration knobs, validated against the deterministic
FD baselines we already have. This unlocks gradient calibration (calibrate-N → predict-N+1) with
Gauss-Newton instead of black-box search, and — because the trajectory is deterministic — every
gradient claim is falsifiable against FD to near machine precision.

---

## 1. The calculus: dH_ij/d(parameter) for our operator

### 1.1 Structure of the operator

Per pair (i, j), with source segment parameterized by s ∈ [0,1] and target by t ∈ [0,1]:

```
H_ij = L_i ∫₀¹ g(s; p) ds
g(s; p) = Σ_{k ∈ visible intervals} ½ | sin θ(t_{k+1}(s,p)) − sin θ(t_k(s,p)) |
```

where p is any parameter (vertex coordinate, and later any knob that moves geometry). The cut
points t_k are, from `_point_segment_exchange`:

- fixed endpoints 0, 1 (independent of p);
- facing-clip roots: `t = −constant/slope`, both affine in the coordinates → rational
  (degree-1/degree-1) in p — a **Möbius/projective function**, smooth wherever slope ≠ 0;
- projective images of blocker endpoints: `t = (r × d_e)/(d_e × d_b)` — again a ratio of
  bilinear forms in the coordinates, smooth away from the pole `denominator = 0` (the pole is
  exactly the configuration where the endpoint image leaves through infinity, i.e. the cut has
  already exited [0,1]).

The outer panel boundaries s_m (from `_panel_events`) are roots of conditions **linear in s**
with affine-in-p coefficients, so `s_m(p) = −c(p)/m(p)` is likewise projective and smooth away
from `m = 0`.

### 1.2 Leibniz/Reynolds applied

The general rule ([Leibniz integral rule / Reynolds transport theorem](https://en.wikipedia.org/wiki/Reynolds_transport_theorem);
1-D form in [Cimbala's note](https://www.me.psu.edu/cimbala/me521web_Fall_2007/Lectures/Leibniz_Theorem_and_RTT.pdf)) for
`H(p) = ∫ g(s;p) ds` with breakpoints s_m(p):

```
dH/dp = Σ_m ∫_{s_m}^{s_{m+1}} ∂g/∂p ds  +  Σ_m [ g(s_m⁻;p) − g(s_m⁺;p) ] · ds_m/dp
```

Two layers of "discontinuity" must be examined:

**Inner layer (target integral).** Already done symbolically. The inner integral was evaluated
in closed form per visible interval, so differentiating it is just the chain rule through the
closed-form endpoint functions t_k(s,p):

```
∂/∂p [ ½|sin θ(t_{k+1}) − sin θ(t_k)| ] = ± ½ [ cosθ·θ'(t_{k+1})·∂t_{k+1}/∂p − (same at t_k) ]
        + direct dependence of sinθ on p at fixed t
```

There is **no Dirac term** at interval birth/death: a shadow interval is born with
t_{k+1} = t_k (zero width), where its contribution and its first derivative in the
birth-parameter are both continuous (the ½|sinθ−sinθ| factor vanishes to first order in the
interval width times the sine slope — value continuous; the derivative acquires only a kink,
i.e. a jump in the *second* derivative, at the event itself). This is the structural payoff of
having integrated the discontinuous visibility indicator *exactly* per piece: the Teg analysis
([Bangaru, Michel et al. 2021](https://dl.acm.org/doi/10.1145/3450626.3459775); follow-up
[Michel et al. 2024, "Distributions for Compositionally Differentiating Parametric Discontinuities"](https://dl.acm.org/doi/10.1145/3649843))
shows the delta contributions arise from differentiating *under* an integral sign across a moving
discontinuity of the *integrand*; once the inner integral is symbolic, the "delta" has already
been absorbed into the endpoint chain-rule terms.

**Outer layer (source integral).** The question is whether g(s;p) jumps across a panel event
s_m. Enumerate the event types (`_panel_events`):

1. *Interior shadow-interval birth/death* (source point crosses the extension line of a
   blocker): the two endpoint images coincide then separate — interval width → 0 continuously.
   g continuous.
2. *Cut entering/leaving through a target endpoint* (collinearity of source point, blocker
   endpoint, target endpoint): just past the event, the newly blocked (or revealed) sliver near
   the target endpoint has width → 0 continuously. g continuous.
3. *Facing-clip root crossing a target endpoint*: same zero-width argument. g continuous.
4. *Projective pole* (blocker endpoint image → ±∞): the cut has left [0,1] before the pole is
   reached; g locally independent of that cut. g continuous.

Hence `g(s_m⁻) = g(s_m⁺)` at every event: **the outer boundary terms vanish**, and

```
dH_ij/dp = L_i Σ_m ∫_{s_m(p)}^{s_{m+1}(p)} ∂g/∂p (s;p) ds   (+ dL_i/dp · ∫ g ds when p moves segment i's endpoints)
```

with ∂g/∂p available in closed form per panel. The panels still must be split at s_m because
∂g/∂s and ∂g/∂p are only piecewise-smooth — the same certified adaptive Simpson machinery
applies, now integrating the (integrand, derivative-integrand) pair jointly so the certificate
covers both.

### 1.3 When is H differentiable at all?

- **Generic configurations:** H is C^1 (above), and piecewise-analytic in p; the kink set (where
  some event coincides with a degeneracy) has measure zero in knob space.
- **Degenerate/grazing configurations** (one-sided derivatives only, Clarke subdifferential):
  - a blocker edge exactly collinear with the source segment (event condition degenerates:
    slope = 0 and constant = 0 simultaneously — the event holds for *all* s, so a whole panel
    sits at a visibility crease);
  - a blocker endpoint exactly touching the target segment or a connector endpoint (the
    tolerance-banded predicates `_proper_connector_hit` / `_facing_cells` change branch);
  - coincident segment endpoints in `_crossed_string_cells` (the |·| of the crossed-string sum
    and the four sqrt distances are non-smooth when a distance → 0 — for closed polygonal
    surfaces adjacent segments share endpoints, so the *unobstructed* term's derivative at the
    shared-endpoint pair must be taken as a directional derivative; it is finite because the
    crossed-string combination cancels the singular direction, but AD through `sqrt(0)` must be
    guarded).
- **Code-level non-smoothness to freeze or guard for AD** (these are properties of the
  implementation, not the mathematics):
  - `min(max(total·L, 0), unobstructed)` clamp (line ~578), escape clamping and excess-row
    renormalization at the end of `build_deterministic_line_exchange_2d` — all active only
    within tolerance bands; freeze the branch during a derivative epoch;
  - tolerance-band predicates (`> tolerance` tests) — classify once at the primal, freeze;
  - the `adaptive_refinement` **fallback** path — a different algorithm whose value agrees within
    tolerance but whose derivative path differs; a derivative epoch must either forbid fallback
    (raise) or FD that pair;
  - adaptive panel subdivision — freeze the panel tree from the primal pass and differentiate
    the *frozen discretization* (standard "differentiate-the-discretization" discipline; the
    certificate should be extended to the derivative integrand, see Stage A).

This is exactly the "differentiate integrals whose discontinuity locations are differentiable
functions of the parameter" regime formalized for rendering by the
[SIGGRAPH 2020 course "Physics-Based Differentiable Rendering: A Comprehensive Introduction"
(Zhao, Jakob, Li)](https://shuangz.com/courses/pbdr-course-sg20/downloads/pbdr-course-sg20-notes.pdf),
which derives the interior+boundary split from the Reynolds transport theorem — except our
boundary terms are provably zero and our interior terms are closed-form.

---

## 2. Differentiable-rendering literature mapped onto our structure

| Concept | Source | What it does | Transfer to petch |
|---|---|---|---|
| Edge/silhouette sampling | [Li et al. 2018, redner](https://dl.acm.org/doi/10.1145/3272127.3275109); [GitHub](https://github.com/BachiLi/redner) | Monte Carlo estimate of the boundary integral over silhouette edges | We *enumerate* the "silhouette" events analytically; their boundary integral is our (vanishing) jump-terms sum. Their importance-sampling machinery is unnecessary in 2-D. Their formulation is the correctness yardstick: any AD of our operator must agree with the interior+boundary decomposition. |
| Reparameterization | [Loubet, Holzschuch, Jakob 2019](https://dl.acm.org/doi/10.1145/3355089.3356510) ([PDF](https://rgl.s3.eu-central-1.amazonaws.com/media/papers/Loubet2019Reparameterizing.pdf)) | Change of variables so discontinuities do not move w.r.t. samples | Conceptually what our per-panel parameterization already achieves *exactly*: within a panel, the visibility structure is fixed, so the integrand is smooth in (s, p) jointly. Panels-split-at-events *is* the exact reparameterization. |
| Warped-area sampling | [Bangaru, Li, Durand 2020](https://dl.acm.org/doi/10.1145/3414685.3417833) | Divergence theorem converts boundary integral → area integral, sampled unbiasedly | Not needed (boundary terms vanish/are closed-form); relevant only if we later go to true 3-D hard visibility where enumeration fails — then WAS or [projective sampling (Zhang et al. 2023)](https://dl.acm.org/doi/10.1145/3618385) are the fallbacks. |
| Path-space differentiable rendering | [Zhang, Miller, Yan, Gkioulekas, Zhao 2020](https://shuangz.com/projects/psdr-sg20/); [boundary-integral estimation follow-up 2022](https://dl.acm.org/doi/10.1145/3528223.3530080) | Differential path integrals with separated interior + boundary path integrals | The multi-bounce analogue: our radiosity Neumann series is the diffuse path integral. Because our single-bounce operator's derivative is exact, the multi-bounce derivative follows by differentiating the linear solve (Section 3.1) — no boundary *path* sampling needed. |
| Teg: compiler for parametric discontinuities | [Bangaru, Michel et al. 2021](https://dl.acm.org/doi/10.1145/3450626.3459775); [Michel et al. 2024](https://dl.acm.org/doi/10.1145/3649843) | Language-level Leibniz rule; Dirac deltas at parametric discontinuities; proved closure under differentiation | The formal semantics for what Section 1 does by hand. Their correctness theorems justify freeze-branch-then-chain-rule. If we ever want machine-checked derivative code for the kernels, Teg's semantics is the spec. |
| Analytic occlusion Jacobians | [Arvo 1994, Irradiance Jacobian](https://dl.acm.org/doi/10.1145/192161.192250) | Closed-form spatial gradients of irradiance from partially occluded polygonal luminaires, incl. occlusion-boundary motion | **Closest prior art.** Same idea one dimension up: differentiate a clipped analytic form factor through the motion of occlusion boundaries. Validates the whole approach; also warns that the bookkeeping (which vertex generates which boundary) is the hard part — which our `_panel_events`/cuts machinery already does. |
| Differentiable radiosity | [Gaussian-surfel adapted radiosity, 2025 (arXiv 2509.18497)](https://arxiv.org/pdf/2509.18497); classical analytic form factors [Baum et al. 1989](https://dl.acm.org/doi/10.1145/74333.74367) | Recent inverse-rendering radiosity differentiates smoothed/soft visibility; classical work computed analytic (but undifferentiated) form factors | No prior art found that differentiates an *exact occluded* radiosity operator through visibility-boundary motion — Arvo 1994 is irradiance (one bounce, one receiver point). A short methods paper is plausibly available here ("exactly differentiable occluded radiosity in 2-D"), which also matters for the Resona differentiable-engine pitch framing. |
| Frameworks | [Dr.Jit / Mitsuba 3 (Jakob et al. 2022)](https://dl.acm.org/doi/abs/10.1145/3528223.3530099) | JIT megakernel tracing for differentiable MC rendering | Evidence that the rendering world needed heavy compiler machinery because their operators are *sampled*; ours is closed-form and small — see Section 6 for why we do not need this. |

---

## 3. Through the solve chain

### 3.1 Radiosity linear solve — trivial (implicit differentiation)

With `(I − Γ F) B = E` (F from `transfer_fraction`, Γ reflectances/sticking complements):
`dB = (I − ΓF)⁻¹ (dΓ F B + Γ dF B + dE)`. Adjoint: one transposed solve per objective,
independent of parameter count. This is textbook implicit differentiation
([Blondel et al. 2022, "Efficient and Modular Implicit Differentiation" (jaxopt)](https://arxiv.org/abs/2105.15183);
[docs](https://jaxopt.github.io/stable/implicit_diff.html)): differentiate the optimality/fixed-point
condition, never the iterative solver. The same recipe covers any charging fixed point later
(the END_STATE gate's "implicit differentiation of the fixed point" clause), *provided* the
fixed-point map itself is differentiable — which is what today's operator delivers for the
neutral-transport factor.

### 3.2 Surface-chemistry step

Coverage steady states (mixed-layer/site-balance algebra): differentiate the residual, IFT again
— one small linear solve per face, closed form for 1–3 site species. Time-stepped stiff ODEs:
discrete adjoint of the actual integrator (implicit Euler/BDF step: differentiate the step's
Newton solve by IFT). Reference implementation pattern:
[PETSc TSAdjoint (Zhang et al., arXiv 1912.07696)](https://arxiv.org/pdf/1912.07696). Chemistry
nonsmoothness (max(0, ·) rate clamps, regime switches like etch-stop when Fs → 0) are genuine
kinks: freeze the active set within a derivative epoch and report one-sided derivatives at
switch points — same discipline as the transport branches.

### 3.3 Level-set evolution over many steps

- **Continuous adjoint / shape-derivative lineage:** level-set inverse problems begin with
  [Santosa 1996](https://www.esaim-cocv.org/articles/cocv/abs/1996/01/cocv-Vol1.17/cocv-Vol1.17.html)
  and the shape-calculus of Sokolowski & Zolésio (Introduction to Shape Optimization, 1992);
  the canonical HJ-adjoint practice is
  [Allaire, Jouve, Toader 2004, "Structural optimization using sensitivity analysis and a level-set method", JCP 194:363-393](https://ui.adsabs.harvard.edu/abs/2004JCoPh.194..363A/abstract):
  shape derivative via adjoint → normal velocity → HJ transport. This gives *first-order shape
  gradients of endpoint functionals* and is the right formalism for Stage D inverse design.
- **Discrete adjoint:** differentiate the actual upwind/WENO update. Recent practice:
  [discrete-adjoint level-set topology optimization (arXiv 2205.09807)](https://arxiv.org/pdf/2205.09807).
  Discrete beats continuous here because our observables are extracted from the discrete φ
  (profile metrics), and discrete adjoints are exactly consistent with the FD baselines we
  validate against.
- **Kinks/shocks caveat:** HJ viscosity solutions develop corners (trench bottom corners,
  merging fronts). Linearized/adjoint consistency at shocks is the classic hazard —
  [Giles & Ulbrich 2010, SINUM 48:882-904 (Part 1)](https://epubs.siam.org/doi/10.1137/080727464)
  and Part 2 (905-921) prove convergence of discrete adjoints only for particular schemes, and
  [Ulbrich's shift-differentiability calculus](https://epubs.siam.org/doi/10.1137/S0363012900370764)
  is the honest sensitivity notion across moving discontinuities. Practical consequence: expect
  clean gradients while the front is smooth and graceful degradation (mesh-dependent smearing of
  the adjoint at corners), not blow-up — dissipative upwinding regularizes the adjoint the same
  way it regularizes the primal.
- **Reinitialization** (our `reinit` = fast marching via skfmm; 3-D `reinit_fsm/narrow/cr2`):
  known to perturb the interface and hinder gradient correctness. Literature practice:
  (a) *constrained reinitialization* that preserves the zero level set
  ([Hartmann et al., JCP 2010](https://www.ljll.fr/~frey/papers/distance/Hartmann%20D.,%20The%20constrained%20reinitialization%20equation%20for%20level%20set%20methods.pdf));
  (b) *gradient-preserving reinitialization*
  ([arXiv 1504.02064, J. Sci. Comput. 2017](https://arxiv.org/abs/1504.02064));
  (c) avoid reinit entirely via distance-regularization energy
  ([variational distance-regularized level set, CMAME 2017](https://www.sciencedirect.com/science/article/abs/pii/S0045782516314062))
  or RBF/global parameterizations
  ([extended level-set method](https://www.sciencedirect.com/science/article/pii/S0021999106002968)).
  **Recommended petch policy:** fast marching is non-differentiable (upwind ordering branches);
  do NOT differentiate through it. Since exact-reinit preserves the zero level set to scheme
  order, treat reinit in the adjoint as the identity on interface perturbations (project the
  incoming adjoint onto the interface-normal component) — this is the de facto
  topology-optimization practice — and *validate this approximation against FD* (Stage C gate).
  If FD disagrees, fall back to a few steps of differentiable PDE-reinit (Sussman-Smereka-Osher
  relaxation with smoothed sign) used *only inside derivative epochs*, or to constrained reinit
  (a)/(b).
- **Conservative remap** (`feature_step_3d` area-conservative bounded remap; 2-D analogue):
  piecewise-linear in geometry with branchy stencil assignment. Same discipline: freeze the
  overlap stencil at the primal, differentiate the (then-linear) transfer weights. The remap's
  structured refusals (topology events) define natural epoch boundaries: gradients are only
  claimed within a fixed-topology window (Section 4).
- **Memory:** checkpointed adjoints are solved technology —
  [Griewank & Walther, revolve, ACM TOMS 26(1), 2000](https://dl.acm.org/doi/10.1145/347837.347846)
  (logarithmic memory, provably optimal), multistage variants in
  [PETSc TSAdjoint](https://arxiv.org/pdf/1912.07696) and
  [optimal multistage checkpointing (arXiv 2106.13879)](https://arxiv.org/pdf/2106.13879).
  For petch 2-D trajectories (10²–10³ steps, per-step state = φ + coverages + fluxes, MBs), a
  checkpoint-every-step store is affordable; revolve is a later luxury.

### 3.4 Whole-chain assembly

Per time step: geometry → exchange matrix H → radiosity solve → fluxes → chemistry rates →
normal velocity → HJ advect (→ occasional reinit/remap). Every arrow above now has a defined
derivative: closed-form (H), IFT (solves), discrete chain rule (advect), frozen-projection
(reinit/remap). The full trajectory derivative is a standard discrete adjoint/tangent
composition over steps with per-epoch frozen branch structure.

---

## 4. The chaos question

**Why this system should be (and empirically is) sensitivity-friendly.** Etch front motion is a
*dissipative, monotone* front propagation: normal flow with non-negative etch speed is a
contraction on front perturbations (viscosity-solution comparison principle — the same property
that makes the primal robust). There is no feedback loop with positive Lyapunov exponents in the
neutral-transport + chemistry chain; the radiosity fixed point is a contraction (spectral radius
of ΓF < 1 via escape); coverage dynamics relax to stable quasi-steady states. The observed
smooth 2–3 nm endpoint responses to 0.01–0.02 knob perturbations are consistent with a
non-chaotic parameter-to-endpoint map. This is qualitatively unlike turbulent-flow sensitivity,
where adjoints diverge exponentially.

**Where trouble can still arise (and what it looks like):**
1. *Topology events* (pinch-off, breakthrough, mask collapse, remap refusals): genuine
   discontinuities of the endpoint map. Policy: gradients are defined within fixed-topology
   epochs; calibration windows must not straddle an event (the remap's structured refusal
   already detects them — reuse it as the gradient-epoch guard).
2. *Grid-crossing staircase*: FD of a discretized front can look noisy at sub-cell amplitude
   even when the underlying map is smooth. The discrete adjoint/tangent of the *frozen
   discretization* does not see this noise; compare against Richardson FD at step sizes above
   the staircase floor.
3. *Corner/shock adjoint smearing* (Section 3.3, Giles-Ulbrich): degraded accuracy near
   profile corners, not divergence.
4. *Charging coupling later*: the charging fixed point with ill-conditioned response maps
   (condition numbers 10⁴–10⁵ per END_STATE_VERIFICATION) is the one place where
   long-horizon sensitivities might misbehave once coupled.

**Fallbacks if long-horizon sensitivities ever misbehave:** the chaotic-sensitivity toolbox —
[Least-Squares Shadowing (Wang, Hu, Blonigan 2014, arXiv 1204.0159 lineage)](https://arxiv.org/pdf/1204.0159),
[NILSS (Ni & Wang 2017, JCP)](https://www.sciencedirect.com/science/article/abs/pii/S0021999117304783)
([arXiv 1611.00880](https://arxiv.org/pdf/1611.00880)),
[adjoint NILSS (Ni, JCP 2018)](https://www.sciencedirect.com/science/article/abs/pii/S0021999117305739),
[FD-NILSS (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0021999119304115),
and [multiple-shooting shadowing](https://www.sciencedirect.com/science/article/abs/pii/S002199911730791X).
NILSS cost scales with the number of positive Lyapunov exponents — for petch that count is
almost certainly zero, which is itself a cheap diagnostic to run (tangent-growth test over a
trajectory). The pragmatic first fallback is simpler: truncated-horizon/windowed adjoints
(calibrate on segments, chain via checkpointed restarts), standard in PDE-constrained
optimization and directly supported by revolve-style checkpointing.

---

## 5. Staged plan for petch

**Stage A — exchange-operator gradients alone (spike, first milestone).**
Forward-mode dual numbers over the scalar kernels of `deterministic_exchange_2d.py`:
`_crossed_string_cells`, `_point_segment_exchange` cut formulas, `_panel_events` roots, Simpson
accumulation with the panel tree and every branch/classification **frozen from the primal pass**.
Deliverable: `dH/d(vertex coordinate)` for all pairs; gate test `tests/test_exchange_gradients.py`
comparing against Richardson central FD (deterministic ⇒ agreement to ~1e-8 relative away from
degenerate configs; one-sided FD agreement at constructed grazing configs). Also extend the
Simpson certificate to the derivative integrand (integrate the (g, ∂g/∂p) pair; accept a panel
only when both converge). Risks: sqrt(0) at shared endpoints (guard), fallback path (forbid in
derivative epochs). *Effort: small — the module is 832 lines of scalar math and was explicitly
written to be scalar (module docstring).*

**Stage B — one-step gradients through transport + chemistry.**
Chain Stage A into the radiosity solve (adjoint/tangent of the linear system, Section 3.1) and
one chemistry evaluation (IFT on coverage steady state). Parameters: 2–4 process knobs
(e.g. cal_F-analogue neutral flux scale, sticking/recombination, ion yield factor). Gate: FD
agreement on d(per-face etch rate)/d(knob) for a static trench. *This already produces the
sensitivity map partners ask for ("which knob moves the bottom rate").*

**Stage C — full-trajectory tangent/adjoint w.r.t. 2 calibration knobs. ← first commercial value.**
For ≤10 knobs, run **tangent (forward) mode** — no reverse sweep, no checkpointing, no reinit
adjoint subtleties: propagate (φ̇, coverage-dot) alongside the primal, with reinit handled by
re-projecting φ̇ onto the interface (Section 3.3 policy) and topology-epoch guards active.
Gate: d(profile observables: depth, CD, bow, taper)/d(knob₁, knob₂) vs the existing
deterministic FD baselines over a full etch trajectory. **Value:** gradient calibration —
Gauss-Newton on (observables − experiment) over 2–10 knobs converges in a handful of forward
runs vs dozens for black-box search, gives local identifiability (J^T J spectrum) for free, and
makes calibrate-N-predict-N+1 a least-squares story instead of a search story. This is the first
stage a design partner pays for, and it does not require the adjoint at all.

**Stage D — geometry/shape gradients (adjoint required).**
Many parameters (per-vertex, mask-shape, or level-set-parameterized initial geometry) ⇒ reverse
mode: discrete adjoint over the step chain with checkpointing
([revolve](https://dl.acm.org/doi/10.1145/347837.347846) if memory ever binds), shape-derivative
assembly in the Allaire-Jouve-Toader pattern for inverse design ("what mask/recipe yields this
profile"). This is the differentiable-engine moonshot deliverable; it reuses Stage A-C parts
plus a reverse sweep and the reinit-adjoint projection validated in Stage C.

Ordering rationale: A→B→C is strictly increasing scope with an FD gate at each step, and C is
reachable without solving the two genuinely research-grade items (adjoint-through-reinit
correctness at corners; boundary terms for true 3-D visibility), which are deferred to D and to
the 3-D lane respectively.

---

## 6. Framework: hand-rolled vs JAX/PyTorch vs Enzyme

- **Hand-rolled forward-mode (recommended for Stages A–C).** The kernels are deliberately
  scalar Python floats (module docstring: per-point work must not pay array dispatch on
  two-element vectors); the operator's derivative structure is closed-form; the branch-freezing
  discipline must be explicit anyway. A ~100-line dual-number class (value, partials) threaded
  through the kernel functions — or a source-copied `_grad` variant of each kernel — keeps the
  certified-quadrature logic byte-identical and testable against the primal. Complex-step
  differentiation is a good *secondary validator* on branch-frozen paths (matches the
  END_STATE gate's "finite-difference/complex-step agreement" language) but cannot be primary
  because of `abs`/branches. Adjoints for the linear solves are 5-line transposed solves.
- **JAX/PyTorch rewrite: not now.** The operator is control-flow-heavy (heapq refinement,
  data-dependent cut lists, per-pair dynamic candidate sets) — hostile to jit/vmap without a
  padded-static restructure; PyTorch adds per-op overhead on scalar math. A JAX port pays off
  only where arrays dominate: the *chemistry + level-set inner loops* (Stage C/D) are
  numpy-array code and are the natural JAX or hand-tangent candidates. If a rewrite is ever
  done, do it operator-by-operator with the hand-rolled gradients as the oracle. jaxopt-style
  [implicit-diff decorators](https://jaxopt.github.io/stable/implicit_diff.html) are the model
  for the solve wrappers either way.
- **Enzyme: only if kernels move to C/C++/Rust.** [Enzyme (Moses & Churavy, NeurIPS 2020)](https://arxiv.org/abs/2010.01709)
  ([enzyme.mit.edu](https://enzyme.mit.edu/)) differentiates LLVM IR — it does not apply to
  CPython float math. If the pair loop is ever ported to a compiled kernel for speed (plausible,
  given the multiprocessing pool already exists), Enzyme becomes the zero-maintenance way to
  keep gradients in sync; [Enzyme-JAX](https://github.com/EnzymeAD/Enzyme-JAX/blob/main/README.md)
  would then bridge it into a JAX-side chemistry/level-set stack. Not on the critical path.
- **Teg** ([Bangaru/Michel 2021](https://dl.acm.org/doi/10.1145/3450626.3459775)) is a research
  DSL, not production tooling — use its *semantics* as the spec for kernel derivatives, not its
  implementation.

---

## 7. Source index

Differentiable rendering / discontinuities: [Li et al. 2018 (redner)](https://dl.acm.org/doi/10.1145/3272127.3275109) · [redner GitHub](https://github.com/BachiLi/redner) · [Loubet et al. 2019](https://dl.acm.org/doi/10.1145/3355089.3356510) · [Bangaru et al. 2020 WAS](https://dl.acm.org/doi/10.1145/3414685.3417833) · [Zhang et al. 2020 PSDR](https://shuangz.com/projects/psdr-sg20/) · [PSDR boundary-integral estimation 2022](https://dl.acm.org/doi/10.1145/3528223.3530080) · [Teg 2021](https://dl.acm.org/doi/10.1145/3450626.3459775) · [Distributions follow-up 2024](https://dl.acm.org/doi/10.1145/3649843) · [PBDR course notes (Zhao/Jakob/Li 2020)](https://shuangz.com/courses/pbdr-course-sg20/downloads/pbdr-course-sg20-notes.pdf) · [PBDR survey 2025 (arXiv 2504.01402)](https://arxiv.org/pdf/2504.01402) · [projective sampling 2023](https://dl.acm.org/doi/10.1145/3618385) · [Dr.Jit 2022](https://dl.acm.org/doi/abs/10.1145/3528223.3530099).
Radiosity/analytic: [Arvo 1994 Irradiance Jacobian](https://dl.acm.org/doi/10.1145/192161.192250) · [Baum et al. 1989 analytic form factors](https://dl.acm.org/doi/10.1145/74333.74367) · [Gaussian-surfel radiosity 2025](https://arxiv.org/pdf/2509.18497).
Calculus: [Reynolds transport theorem](https://en.wikipedia.org/wiki/Reynolds_transport_theorem) · [Leibniz/RTT note (Cimbala)](https://www.me.psu.edu/cimbala/me521web_Fall_2007/Lectures/Leibniz_Theorem_and_RTT.pdf).
Level set / adjoint: [Santosa 1996](https://www.esaim-cocv.org/articles/cocv/abs/1996/01/cocv-Vol1.17/cocv-Vol1.17.html) · [Allaire-Jouve-Toader 2004](https://ui.adsabs.harvard.edu/abs/2004JCoPh.194..363A/abstract) · [discrete-adjoint level-set TO (arXiv 2205.09807)](https://arxiv.org/pdf/2205.09807) · [Hartmann constrained reinit](https://www.ljll.fr/~frey/papers/distance/Hartmann%20D.,%20The%20constrained%20reinitialization%20equation%20for%20level%20set%20methods.pdf) · [gradient-preserving reinit (arXiv 1504.02064)](https://arxiv.org/abs/1504.02064) · [distance-regularized LS TO](https://www.sciencedirect.com/science/article/abs/pii/S0045782516314062) · [extended LS method](https://www.sciencedirect.com/science/article/pii/S0021999106002968).
Shocks/adjoint consistency: [Giles & Ulbrich Part 1](https://epubs.siam.org/doi/10.1137/080727464) · [Ulbrich shift-differentiability](https://epubs.siam.org/doi/10.1137/S0363012900370764).
Chaos/shadowing: [LSS (arXiv 1204.0159)](https://arxiv.org/pdf/1204.0159) · [NILSS](https://www.sciencedirect.com/science/article/abs/pii/S0021999117304783) · [adjoint NILSS](https://www.sciencedirect.com/science/article/abs/pii/S0021999117305739) · [FD-NILSS](https://www.sciencedirect.com/science/article/abs/pii/S0021999119304115) · [multiple-shooting shadowing](https://www.sciencedirect.com/science/article/abs/pii/S002199911730791X).
Implicit diff / adjoint infrastructure: [Blondel et al. 2022 (arXiv 2105.15183)](https://arxiv.org/abs/2105.15183) · [jaxopt implicit diff](https://jaxopt.github.io/stable/implicit_diff.html) · [revolve (Griewank & Walther)](https://dl.acm.org/doi/10.1145/347837.347846) · [PETSc TSAdjoint](https://arxiv.org/pdf/1912.07696) · [optimal multistage checkpointing](https://arxiv.org/pdf/2106.13879).
AD tooling: [Enzyme (NeurIPS 2020)](https://arxiv.org/abs/2010.01709) · [enzyme.mit.edu](https://enzyme.mit.edu/) · [Enzyme-JAX](https://github.com/EnzymeAD/Enzyme-JAX/blob/main/README.md).
Adjacent (litho, for the commercial framing): [gradient-based inverse lithography / ILT](https://arxiv.org/html/2408.08969v1) · [hierarchical ILT 2025](https://doi.org/10.3390/mi16070798) — inverse lithography is already gradient-native; a differentiable etch engine is the missing sibling.
