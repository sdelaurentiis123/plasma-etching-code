# Research — event-resolved (atom-resolved) surface chemistry for the mixed layer

Date: 2026-07-26 (drafted 2026-07-25). Author: Opus research agent (parallel to the
overnight per-event chemistry build). **Do not commit.**

Scope: answer the five research questions in `OVERNIGHT_PLAN_2026-07-25.md` §Research —
how the established feature-scale codes couple *per-particle* energy/angle to surface
reaction probability, the published Jensen/compression bias for IEAD-driven yields, the
numerical patterns for atom-resolved surface kinetics at scale, the hot-neutral yield-law
question, and a verdict on petch's planned formulation (state frozen per substep +
per-atom rates summed per face).

Locators that were not opened to the primary text in this session are marked **[VERIFY]**.
Journal/volume/page anchors that were confirmed through search result metadata are called
out as such; the exact in-text quotes behind them still warrant a read before they enter a
citable claim.

---

## 0. The one-sentence finding

Every established feature-scale etch code — Kushner's MCFPM, the Zhang/Economou →
Krüger fluorocarbon detailed surface model, Kuboi's voxel-slab, and ViennaPS's ray
tracer — evaluates the **nonlinear yield at each particle's own (E, θ) and only then
accumulates over particles**. None of them evaluates a yield at a flux-weighted mean
(E, cos). The "yield-at-summed-moments" order that ml8 exposed is not an approximation
any of the reference codes make; it is a compression petch's adapter introduced at the
`ModuleFluxes` boundary. petch's planned move (per-event yield, segment-summed per face)
is the **correct and standard** order of operations. Confirmed.

Notably petch already *contains* the correct pattern: `FaceResolvedEnergeticFlux.
yield_rate_m2_s` (src/petch/surface_kinetics.py:339–358) computes
`event_yield = yield_law.evaluate(event_energy_eV, event_cosine_incidence)` then
`np.bincount(event_face, weights=event_flux * event_yield)`. That is per-event-then-sum.
The mixed-layer adapter (`_module_fluxes`, mixed_layer_mechanism.py:343–431) throws that
away by collapsing each population to `mean_energy`/`mean_cosine` before the module's
stopping-table and sputter laws run. The fix is to push the surface-kinetics pattern that
already exists through the mixed-layer channels.

---

## 1. Per-particle → surface reaction coupling in the reference codes

### 1.1 Kushner MCFPM (Monte Carlo Feature Profile Model)

Order of operations is **pseudo-particle by pseudo-particle**:

- Gas-phase ions and neutrals are represented by pseudo-particles ("computational
  particles"), each launched toward the feature with an energy and angle **sampled from
  the ion/neutral energy-angle distributions (EADs/IEADs)** supplied by the reactor-scale
  Plasma Chemistry Monte Carlo Model. The EAD is the sampling distribution; individual
  particles carry individual (E, θ). (Huard PhD thesis, U. Michigan; Huang et al. 2019.)
- A particle's trajectory is traced (ballistic for neutrals, sheath-accelerated /
  optionally re-emitted for ions) until it **strikes a mesh cell**. At the impact, a
  reaction is chosen by Monte Carlo **from probability arrays indexed by the striking
  particle's energy and angle and by the identity/state of the struck cell**. The chosen
  reaction mutates that cell (etch = remove cell, deposit = add cell, mix/implant =
  change cell identity). ("When a pseudoparticle hits the surface, a reaction is chosen
  based on probability arrays using Monte Carlo techniques, with the selected reaction
  changing the identities of mesh cells, or cells being added or removed" — Zhang/Economou
  detailed-surface-model lineage; the same mechanism is inherited by MCFPM.)
- **State coupling is per impact and sequential**: the cell state the next particle sees
  already reflects the reaction the previous particle caused. There is no averaging of
  energy across the arriving population before the yield/branching is evaluated — the
  branching *is* the per-event yield evaluation.

This is the strongest confirmation of the plan: the authoritative kinetic code does not
merely "sum per-event yields," it evaluates the reaction **against the live, mutating
surface state one particle at a time**. That is a stricter ordering than petch's planned
"state frozen per substep." See §5 for why the frozen-substep relaxation is nonetheless
the right engineering choice for a deterministic (non-Monte-Carlo) integrator.

Locators:
- C. Huard, *Nano-Scale Feature Profile Modeling of Plasma Material Processing*, PhD
  thesis, University of Michigan (Kushner group). cpseg.eecs.umich.edu/pub/theses/
  huard_chad_phd_thesis.pdf — **[VERIFY]** exact chapter on the per-particle reaction
  selection loop and the EAD-sampling step.
- S. Huang, S. Shim, S. Nam, M. J. Kushner, "Plasma etching of high aspect ratio features
  in SiO2 using Ar/C4F8/O2", *J. Vac. Sci. Technol. A* **37**, 031304 (2019).
  cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf — confirmed via search
  metadata; **[VERIFY]** in-text hot-neutral treatment (see §4). (Direct fetch failed on a
  TLS cert error this session; open on-box.)

### 1.2 Zhang/Economou → Krüger fluorocarbon detailed surface model

petch's `KRUEGER_2024_*` constant blocks (mixed_layer_mechanism.py:458–485) trace to the
detailed fluorocarbon SiO2/Si surface model in the Zhang & Economou / Kushner lineage,
carried into Florian Krüger's more recent work.

- D. Zhang & M. J. Economou, "SiO2 and Si etching in fluorocarbon plasmas: A detailed
  surface model coupled with a complete plasma and profile simulator" (and the companion
  "Etching of SiO2 and Si in fluorocarbon plasmas: A detailed surface model accounting for
  etching and deposition"). The surface is a **fluorocarbon film over a mixed/reactive
  layer**; SiO2 is removed by energetic ions *in the presence of* the FC film, exactly the
  two-reservoir structure petch's `mixed_layer` implements. Reaction selection when a
  pseudoparticle hits is per-particle Monte Carlo against the local layer state.
  academia.edu/16909363 (403 this session — **[VERIFY]** on-box); ADS 2000PhDT........80Z
  (Zhang PhD).
- Florian Krüger is in the Kushner-group lineage; the "Krüger 2024" thesis/constants in
  petch (mixing efficiency, crosslink channel, substrate-dependent deposition
  probabilities ~100× higher on polymer than bare oxide) are **[VERIFY]** — confirm the
  exact document (Krüger PhD thesis, U. Michigan, ~2024) and appendix that publishes the
  complex-formation and deposition-probability tables the code lifts verbatim
  (mixed_layer_mechanism.py:467–485). A Krüger 2024 *Phys. Plasmas* 31, 033508 paper on
  voltage-waveform HAR SiO2 etching exists and is the same author, but is a reactor-scale
  paper, not the surface-constants source.

Key point for the plan: the crosslink term the plan wants per-atom —
`crosslink (E_i − E_iface_i)` — is intrinsically per-event. `E_iface` (the local
crosslinking/activation threshold) is compared to **each ion's** deposited energy; a
flux-weighted-mean ion energy cannot represent the fraction of the population above
threshold when the distribution straddles it. This is the same threshold-straddling
failure as the sqrt-yield Jensen effect (§2) and is why the compression "poisons" faces
that see a bimodal primary+reflected population.

### 1.3 Kuboi voxel-slab model (Sony)

- N. Kuboi et al., "Insights into different etching properties of continuous wave and
  atomic layer etching processes for SiO2 and Si3N4 films using voxel-slab model",
  *J. Vac. Sci. Technol. A* **37**, 051004 (2019). Also the 3D voxel damage-distribution
  papers, *JVSTA* **33**, 061308 (2015). (Confirmed via search metadata; **[VERIFY]**
  in-text.)
- Structure: surface = a **C–F polymer layer** on top of a **reactive layer divided into
  ~0.25 nm lattice-scale slabs**. Ion penetration deposits energy across slabs; each
  slab's state evolves. This is the same two-layer + depth-resolved-deposition physics as
  petch's `mixed_layer` + ZBL stopping table.
- Directly load-bearing finding for ml8: Kuboi reports that **using a monoenergetic ion
  energy improves controllability of both the polymer-layer and reactive-layer
  thickness** relative to a broad IED — i.e. the *width* of the ion energy distribution
  materially changes the surface-layer evolution, not just its mean. That is an
  experimental/simulation statement that the distribution cannot be replaced by its mean
  without changing the answer — the physical counterpart of the Jensen bias.
- Kuboi also uses a **statistical-ensemble** method (IOP, 10.35848/1347-4065/acbebb) to
  carry per-site distributions over large patterns rather than mean fields — same
  philosophy as retaining events instead of moments.

### 1.4 ViennaPS ray-level yield accumulation

- ViennaPS (TU Wien IUE; Manstetten/Quell/Rodgers/… lineage; JOSS/SoftwareX
  "ViennaPS: A flexible framework for semiconductor process simulation",
  ScienceDirect S2352711025004194; docs viennatools.github.io/ViennaPS).
- Mechanism: a top-down Monte Carlo **ray tracer** launches rays with directions/energies
  sampled from source distributions. On intersection with a surface element, the particle
  model's `surfaceCollision` deposits a **per-ray contribution** into that element's rate
  accumulator; the ray may reflect (carrying reduced/redirected energy) or terminate.
  Yield/sticking laws (e.g. sputter yield vs energy and angle) are applied **at the
  collision, per ray**, to that ray's own (E, θ) — then accumulated. After tracing,
  `SurfaceModel::calculateVelocities` maps accumulated rates → normal velocity. ("Particles
  and their impacts must be accumulated to calculate the surface rate at a single
  location" — ViennaPS GPU raytracing thesis, Riedel 2026, repositum.tuwien.at 226294.)
- So ViennaPS is **per-ray-yield-then-accumulate**, identical in spirit to petch's
  `np.bincount(event_face, weights=event_flux * event_yield)`. The accumulation is over the
  nonlinear yield already evaluated per ray, not a yield of an accumulated flux. **[VERIFY]**
  the exact SurfaceModel/particle API names against the current release; the ordering claim
  is robust across versions.

### 1.5 Established order of operations — verdict on Q1

**Correct order: evaluate the (nonlinear) yield per event at that event's own (E, θ)
against the current surface state, then segment-sum per face.** Evaluating a yield at
summed/mean moments (`Y(⟨E⟩, ⟨cos⟩)`) is what every reference code avoids and what ml8
proved harmful. petch's `FaceResolvedEnergeticFlux.yield_rate_m2_s` already does it right;
the mixed-layer adapter must stop compressing.

---

## 2. Published quantification of the Jensen / compression bias for IEAD yields

### 2.1 The literature anchor for petch's sqrt-E remark

petch's `chemistry.py:10–34` docstring states the effect: yield
`Y = A·max(√E − √E_th, 0)`, `√E` concave ⇒ `⟨Y(E)⟩ < Y(⟨E⟩)` (Jensen), so evaluating at
the mean **over-estimates** the yield, and the effect "grows with distribution width and
is sharply nonlinear when the low horn nears E_th."

The literature anchor is **Steinbrüchel's universal square-root law**:

- C. Steinbrüchel, "Universal energy dependence of physical and ion-enhanced chemical etch
  yields at low ion energy", *Appl. Phys. Lett.* **55**(19), 1960–1962 (1989),
  DOI 10.1063/1.102336. Establishes that physical and ion-enhanced chemical etch yields are
  **linear in √E down to a threshold E_th** for Ar/reactive-ion sputtering of Si and SiO2
  and for ion-beam-enhanced chemical etching — i.e. `Y ∝ (√E − √E_th)`. (Confirmed via
  search metadata.) This is the concave law whose curvature produces the Jensen gap; every
  petch yield of the `SteinbruchelYield` form (surface_kinetics.py) inherits it.
- The `viennaps` angular factors and Steinbrüchel √E law together mean petch has **two**
  nonlinear compositions — `f(cos)` and `√E` — so the compression bias is a joint-moment
  bias, not a 1-D one. Collapsing to `(⟨E⟩, ⟨cos⟩)` discards the covariance between energy
  and angle within a face's population, which for the primary-ion + specular-hot-neutral
  mixture is strongly negative (hot neutrals arrive oblique and, after inelastic
  reflection, cooler) — precisely the heterogeneity ml8 flagged.

### 2.2 Directly-published bias quantifications

- **Bimodal RF-sheath IED (arcsine "horns").** The physically correct low-frequency-bias
  sheath IED is bimodal with horns at ±ΔE/2 (arcsine density), not Gaussian. petch's
  `_ied_yield` bimodal branch (chemistry.py:30–34) implements it and cites **Kawamura et
  al. 1999** [VERIFY exact: E. Kawamura, V. Vahedi, M. A. Lieberman, C. K. Birdsall,
  "Ion energy distributions in RF sheaths; review, analysis and simulation", *Plasma
  Sources Sci. Technol.* 8, R45 (1999)]. Because a large fraction of the flux sits in the
  low horn, `⟨√E⟩` is materially below `√⟨E⟩` whenever the low horn approaches E_th — the
  bias is maximized exactly in the threshold-straddling regime, and it is the *shape*, not
  the mean, that sets it. This is the closest published quantification: same-mean IEDs of
  different width give different integrated yields.
- **Krüger 2024 voltage-waveform HAR SiO2** (*Phys. Plasmas* 31, 033508, 2024) and the
  dual-frequency CCP HAR SiO2 studies (MDPI *Materials* 16, 3820, 2023; *Nanomaterials* 12,
  4457, 2022) all report that etch outcome depends on the **IED shape at fixed mean energy**
  (narrow high-frequency IED reduces corner faceting/sputtering) — an experimental
  statement that `Y(⟨E⟩)` is insufficient. **[VERIFY]** whether any of these tabulates the
  mean-vs-integrated yield gap numerically; most report profile consequences, not the bare
  Jensen number.
- **Direct mean-vs-distribution yield comparison paper:** I did **not** find a single paper
  whose headline result is "yield-at-mean minus distribution-integrated yield = X%". The
  effect is treated as folklore-grade obvious (it is Jensen on Steinbrüchel's √E) and is
  usually *avoided by construction* (all the §1 codes integrate over the distribution), so
  nobody publishes the error of the approximation they don't make. **The honest framing for
  petch: the anchor is Steinbrüchel (concave law) + Kawamura (the real bimodal width) +
  Jensen's inequality; the magnitude is петch's own ml-series measurement, which is a
  legitimate contribution rather than a re-derivation.** Petch's two-atom Jensen gate
  (plan §3) is the right way to pin the number in-repo.

Verdict on Q2: the concave-law anchor is **Steinbrüchel 1989 APL 55, 1960**; the
distribution-width-matters anchor is **Kawamura 1999 (bimodal IED)** plus Kuboi's
monoenergetic-vs-broad controllability result (§1.3). A clean published "mean vs
distribution yield error" number for FC/SiO2 does not appear to exist — petch measuring it
is novel, not redundant.

---

## 3. Numerical patterns for atom-resolved surface kinetics at scale

### 3.1 The segment-sum / scatter-add formulation

The operation is a **segmented reduction** keyed by face id:
`rate[face] = Σ_{atoms on face} flux_atom · term(E_atom, cos_atom, state[face])`.

- **CPU (numpy):** two idioms.
  - `np.bincount(event_face, weights=w, minlength=F)` — what petch already uses; fastest
    for a single weighted sum, returns a dense length-F array. Use one bincount per additive
    term (attenuation-weighted flux, film-sputter contribution, Kress term, crosslink
    term, each layer channel). Each term is `weights = event_flux · term_i(E,cos,state)`.
  - `np.add.at(out, event_face, w)` — unbuffered scatter-add; correct for **repeated
    indices** (many atoms → one face) where naive `out[event_face] += w` silently drops
    collisions. Slower than bincount but works for multi-column / non-scalar accumulation.
    petch's adapter already uses `np.add.at` for the flux/energy/cosine moment sums
    (mixed_layer_mechanism.py:361–367) — the same call sites become the per-term yield
    sums.
  - Prefer **bincount per term** over `add.at` when each term reduces to one scalar column;
    it is measurably faster and immutable-friendly (petch caches frozen bincount results).
- **GPU:** the same segmented reduction. Options, in rough order of determinism:
  - **Sort-by-face + segmented reduction** (Modern GPU `segreduce`, CUB
    `DeviceSegmentedReduce`): deterministic, order-independent, the gold standard for
    reproducibility. petch already does a "cupy edge-sort" in the tracer (see MEMORY
    speed-sota-roadmap) — the same sort-then-segment pattern applies to atom→face.
  - **Atomic scatter-add** (`atomicAdd`, `cupyx.scatter_add`, `jax.ops.segment_sum`,
    `torch_scatter.scatter_add`): simplest, but **non-associative float atomics make it
    run-to-run nondeterministic** (the Locality-Aware AD paper, arXiv 2509.00406, and the
    pytorch_scatter docs both flag this). Contention is low here (atoms per face ~O(100),
    faces ~9k) so performance is fine, but for petch's receipt/determinism discipline
    (1e-9 ledger closure, byte-reproducible gates) **prefer sort+segment, or sum in a fixed
    order, so the ledger closes deterministically.**
  - `jax.ops.segment_sum` / `scatter_add` are the differentiable path — relevant to the
    Route-B differentiable-charging moonshot in MEMORY, and worth keeping the formulation
    compatible with (all per-atom terms are smooth functions of (E,cos) except the
    threshold `max(·,0)` kinks, which are subgradient-friendly).

### 3.2 Operator ordering when many atoms share a face in one timestep

Two admissible orderings:

1. **State frozen per substep (petch's plan, "Jacobi"):** all atoms on a face see the
   *same* frozen `state[face]` from the start of the substep; their per-atom rates are
   summed; the summed rate advances the state once. Order-independent, vectorizable,
   deterministic, trivially the same result regardless of atom permutation. **Bias is
   O(dt)**: the state does not relax between atoms within the substep, so a term that
   saturates (e.g. film sputter cannot remove more film than exists, F-layer depletion)
   can be over-driven if the substep delivers more than the reservoir holds. Controlled by
   the substep cap `default_max_step_s` and the overdraw/representational-floor clamps
   (mixed_layer_mechanism.py:305–311) — **but the clamp must be applied to the summed rate,
   not per atom** (see §5 pitfalls).
2. **Sequential per atom ("Gauss-Seidel", the MCFPM ordering):** each atom sees the state
   left by the previous atom. Exact saturation handling (no over-draw possible), but
   **order-dependent** (needs a canonical atom order for reproducibility), serial (kills
   vectorization), and for a deterministic integrator introduces a spurious dependence on
   an arbitrary atom ordering that has no physical meaning (the real ordering is the random
   arrival sequence MCFPM samples; a deterministic code has no such sequence).

**Stability:** the frozen-substep (Jacobi) scheme is stable provided the substep obeys a
CFL-like bound — no reservoir may be driven negative within one substep. petch already has
the machinery: `n_steps = ceil(duration/cap)` subdivides, and the "materially negative"
guard (line 308) catches integrator failure. The added requirement for atoms is that the
**cap must bound the summed per-face draw**, i.e. the many-atom sum can be up to ~100× a
single event's draw, so the effective per-substep depletion is larger — the cap may need to
shrink, or the saturating terms need a per-substep analytic clamp (1 − exp(−Σrate·dt))
rather than a linear Σrate·dt. Recommend the **exponential/rate-limited update for
depleting reservoirs** so that no ordering and no atom count can overshoot, which removes
the CFL sensitivity entirely for those channels.

Verdict on Q3: use **bincount-per-term segment sums** (CPU) / **sort+segmented-reduce**
(GPU) for determinism; adopt **state-frozen-per-substep (Jacobi)**; apply saturation
clamps to the **summed** face rate; and use **rate-limited (exponential) updates for any
depleting reservoir** so atom count/ordering cannot overshoot.

---

## 4. Hot-neutral vs ion yield laws (reflected energetic neutrals)

Question: are Krüger's/MCFPM's yield laws charge-agnostic? Do the same (E, θ) yield laws
apply to specularly-reflected hot neutrals as to primary ions, on FC film and SiO2?

**Answer: yes, the surface yield laws are charge-agnostic — they depend on the projectile's
mass, energy, and angle, not its charge state.** Support:

- **Huang et al. 2019 (JVSTA 37, 031304)** — the HAR SiO2 Ar/C4F8/O2 MCFPM study — treats
  ions that strike a surface and **reflect as "hot neutrals": they neutralize on impact but
  retain most of their energy and, being now uncharged, travel ballistically deeper into
  the feature where they participate in the same reactions as ions** (physical/chemical
  sputtering, mixing, implantation) governed by the **same energy/angle yield functions**.
  In HAR features the hot-neutral channel is essential — it is how energy reaches the
  sidewall/bottom after the sheath-directed ion is lost to a wall. (Confirmed via search
  metadata; **[VERIFY]** the exact statement and whether Huang uses a reduced reflection
  coefficient or energy-loss factor at neutralization — direct fetch failed on TLS this
  session.)
- **Physical basis:** sputtering/knock-on and ion-assisted bond breaking are momentum/
  energy-transfer processes (the Steinbrüchel √E law, ZBL stopping). Charge is neutralized
  by Auger/resonant processes in the top ~few Å and does **not** enter the collision cascade
  energetics. So a 200 eV Ar hot neutral and a 200 eV Ar+ deposit the same energy via the
  same ZBL table petch already uses (`zbl-deposited-in-layer`, mixed_layer_mechanism.py:166)
  and drive the same yield. This is standard and is why MCFPM/ViennaPS apply one yield law
  to "energetic particles" irrespective of charge.
- **Caveats the plan should record (charge is NOT irrelevant everywhere):**
  1. **Transport differs by charge, yield does not.** The neutral is not steered by the
     sheath/local charging field; only its *arrival* (E, θ) distribution differs. So the
     correct decomposition is: charge-dependent **transport** produces two populations with
     different (E, θ) spectra hitting a face; a **single charge-agnostic yield law**
     evaluates each event. This is exactly petch's atom-array design — different events,
     same `yield_law`. 
  2. **Charge-exchange / reflection energy loss.** A specular reflection is rarely elastic:
     the reflected neutral carries a fraction (energy-loss/accommodation coefficient) of the
     incident energy, and the reflected-angle distribution is not a perfect mirror. That is
     a **transport-side** parameter (sets the hot-neutral E, θ), not a yield-law change. The
     plan's "reflection literature_v1" (ml9) is where this coefficient lives — keep it in the
     adapter's event-generation, not in the yield law.
  3. **Chemistry threshold reuse.** Whether a hot CFx neutral still *chemisorbs/deposits*
     like a thermal radical is a separate question from its *sputter* yield — a hot CFx may
     sputter (energetic channel) rather than stick (thermal channel). petch already
     separates these (energetic populations → `mixed_layer_step` ion channels; thermal
     radicals → `precursor_flux`/chemisorption). Route reflected **hot** CFx into the
     energetic channel, not the sticking channel, or it will double-count as both a
     depositor and a sputterer. **[VERIFY]** how Krüger/Huang route reflected polyatomic
     neutrals specifically (monatomic Ar hot neutral is unambiguous; CFx is not).

Verdict on Q4: **the yield laws are charge-agnostic** — apply the same `yield_law.evaluate
(E, cos)` to primary-ion events and reflected-hot-neutral events. The charge dependence
lives entirely in **transport/event generation** (which populations, with which (E, θ), and
energy loss on reflection), which is precisely the boundary petch's atom-array design draws.
Guard the polyatomic-neutral routing (energetic vs thermal channel) against double-counting.

---

## 5. VERDICT — recommended atom-resolved formulation for petch's mixed layer

**The plan is confirmed against the literature.** "State frozen per substep + per-atom rates
summed per face" is the correct, standard, and (for a deterministic integrator) the *right*
formulation. Specifics and guardrails:

### 5.1 Confirmed
- **Per-event yield then segment-sum** is the established order (MCFPM per-particle, ViennaPS
  per-ray, Kuboi per-slab, Zhang/Krüger per-impact). Yield-at-mean-moments is the error, and
  it is unique to petch's current adapter compression — remove it. (§1)
- petch **already implements the correct pattern** at
  `FaceResolvedEnergeticFlux.yield_rate_m2_s`; the work is to route the mixed-layer's
  per-atom ion terms (attenuation `exp(−d_FC/λ(E_i))`, film sputter, Kress, crosslink
  `(E_i − E_iface)`, layer channels) through the same bincount-per-term reduction instead of
  through `mean_energy`/`mean_cosine`. (§1.1, §3)
- **State frozen per substep (Jacobi)** is preferable to sequential-per-atom for a
  deterministic code: order-independent, vectorizable, byte-reproducible, and it does not
  invent a physically-meaningless atom ordering (MCFPM's sequential update is physical only
  because its order *is* the sampled random arrival sequence — a deterministic solver has
  none). (§3.2)
- **Charge-agnostic yield law** for ions and reflected hot neutrals; charge lives in
  transport/event generation. (§4)

### 5.2 Pitfalls to gate against (these are the ways the atom sum can still be wrong)

1. **Saturation clamps must be applied to the summed face rate, NOT per atom.** If each
   atom's film-sputter is independently clamped to "≤ film present," then N atoms can each
   remove up to the *whole* film → N× over-removal. Correct: sum the unclamped per-atom
   removal rates, *then* clamp the sum against the reservoir (or use the rate-limited
   exponential update). The plan's ledger-closure gate (<1e-9) will catch a per-atom clamp
   because it double-counts removed film. **This is the single most likely bug.**
2. **Overdraw scaling / normalization order.** If the module rescales rates when the total
   requested draw exceeds available inventory (overdraw protection), that rescale must
   operate on the **summed** face rate and be applied **uniformly to every atom's
   contribution** (a single per-face scale factor), so the ledger and the state update use
   the same scaled numbers. Rescaling per atom, or rescaling the state update but not the
   emitted-product ledger, breaks conservation. Order: (a) sum per-atom rates per face,
   (b) compute one per-face overdraw scale, (c) apply it to the summed rate *and* to the
   product/exchange ledger identically.
3. **Attenuation is per atom and state-coupled — do not pull `d_FC` out of the sum.**
   `exp(−d_FC/λ(E_i))` depends on both the current film thickness `d_FC` (face state, frozen
   for the substep — fine) and the atom energy `E_i` through `λ(E_i)` (per atom). Because
   `λ(E)` and the sputter yield are both nonlinear in E, this term is a genuine per-atom
   evaluation; the whole motivation of the refactor is that `⟨exp(−d/λ(E))·Y(E)⟩ ≠
   exp(−d/λ(⟨E⟩))·Y(⟨E⟩)`. Keep `d_FC` frozen per substep (Jacobi), vary E per atom.
4. **Crosslink threshold straddle.** `(E_i − E_iface)` with `max(·,0)` is exactly the
   threshold term where compression is worst (a bimodal primary+reflected population can
   have `⟨E⟩ > E_iface` while half the flux is below it, or vice versa). Per-atom evaluation
   is mandatory here; a two-atom gate straddling `E_iface` should show the largest
   correction of any channel. Good canary.
5. **Determinism of the sum.** On GPU, use sort+segmented-reduce (or a fixed summation
   order), not float atomicAdd, so the <1e-9 ledger-closure and single-atom==scalar (1e-12)
   gates stay reproducible. (§3.1)
6. **Empty-face / zero-atom faces.** bincount with `minlength=face_count` and the
   `ion_flux>0` guards already handle faces with no ion atoms; keep the scalar (no-atom)
   path byte-identical when the atom arrays are absent (plan §1 "scalar path unchanged when
   absent") — the single-atom==scalar 1e-12 gate is the correct guard.
7. **Cap may need to shrink.** The many-atom summed draw per substep is larger than a single
   event's; either shrink `default_max_step_s` or (better) switch depleting-reservoir
   channels to rate-limited `(1 − exp(−Σr·dt))` updates so no atom count/ordering can
   overshoot, removing CFL sensitivity. (§3.2)

### 5.3 One structural recommendation beyond the plan
Do the per-atom reduction **inside `mixed_layer.step`** (or a thin per-event kernel it
calls), taking the sparse `(event_face, event_flux, E, cos)` arrays directly, rather than
pre-reducing in the adapter. The adapter's job becomes *only* event assembly (build atom
arrays from `FaceResolvedEnergeticFlux` events and broadcast `EnergeticFlux` spectrum rows
per face — plan §4), with **no moment compression anywhere**. That keeps the single
"yield-at-mean" temptation out of the module boundary permanently and mirrors how MCFPM/
ViennaPS keep the (E, θ) all the way to the reaction site. The current `_module_fluxes`
mean-energy/mean-cosine computation (mixed_layer_mechanism.py:385–388) should be **deleted**
on the atom path, not merely bypassed, so it cannot be reintroduced.

---

## Reference list (verify locators before citing)

1. C. Steinbrüchel, *Appl. Phys. Lett.* **55**(19), 1960–1962 (1989), DOI 10.1063/1.102336 —
   universal √E − √E_th etch-yield law. **Jensen anchor.** (search-confirmed)
2. S. Huang, S. Shim, S. Nam, M. J. Kushner, *J. Vac. Sci. Technol. A* **37**, 031304 (2019)
   — HAR SiO2 Ar/C4F8/O2 MCFPM; hot-neutral (reflected-ion) treatment. **[VERIFY]** in-text
   (TLS fetch failed). cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
3. C. Huard, PhD thesis, U. Michigan (Kushner group) — MCFPM per-particle reaction loop.
   **[VERIFY]** cpseg.eecs.umich.edu/pub/theses/huard_chad_phd_thesis.pdf
4. D. Zhang & M. J. Economou (also Zhang PhD, ADS 2000PhDT........80Z) — detailed FC SiO2/Si
   surface model, two-layer film+reactive, per-impact Monte Carlo reaction. **[VERIFY]**
5. F. Krüger — "Krüger 2024" surface constants in petch (complex-formation, deposition,
   crosslink, mixing). **[VERIFY]** the exact thesis/appendix (Kushner-group lineage, ~2024);
   distinct from Krüger et al. *Phys. Plasmas* **31**, 033508 (2024) (reactor-scale waveform).
6. N. Kuboi et al., *J. Vac. Sci. Technol. A* **37**, 051004 (2019); **33**, 061308 (2015);
   IOP 10.35848/1347-4065/acbebb (statistical-ensemble) — voxel-slab, monoenergetic-vs-broad
   controllability. **[VERIFY]** in-text.
7. ViennaPS — "ViennaPS: A flexible framework for semiconductor process simulation"
   (ScienceDirect S2352711025004194); X. Riedel, "GPU-Accelerated Ray-Tracing for Particle
   Simulations in ViennaPS", TU Wien 2026 (repositum 226294) — per-ray flux accumulation.
   **[VERIFY]** exact SurfaceModel/particle API.
8. E. Kawamura, V. Vahedi, M. A. Lieberman, C. K. Birdsall, *Plasma Sources Sci. Technol.*
   **8**, R45 (1999) — bimodal RF-sheath IED (arcsine horns). **[VERIFY]** exact cite;
   petch chemistry.py cites "Kawamura 1999".
9. Locality-Aware AD on GPU (arXiv 2509.00406); Modern GPU segmented reduction
   (moderngpu.github.io/segreduce); pytorch_scatter (github rusty1s/pytorch_scatter) —
   scatter-add nondeterminism and sort+segment determinism. (search-confirmed)

## In-repo cross-references (verified this session by reading the source)
- src/petch/surface_kinetics.py:339–358 — `FaceResolvedEnergeticFlux.yield_rate_m2_s`:
  the correct per-event-then-bincount pattern already exists.
- src/petch/surface_kinetics.py:242–277 — `EnergeticFlux.mean_yield`: integrates the yield
  over the population's own (E,cos) spectrum weights (also correct — it does NOT collapse to
  a mean energy first; `np.dot(weight, yield_law.evaluate(energy, cosine))`). The bug is not
  here either.
- src/petch/mixed_layer_mechanism.py:356–388 — `_module_fluxes`: **the compression site.**
  Collapses every population to `mean_energy`/`mean_cosine` before the module runs. This is
  what ml8 indicted; delete on the atom path.
- src/petch/mixed_layer_mechanism.py:305–311 — `representational_floor`: the negative-
  reservoir guard; ensure atom-summed clamps route through the summed rate, not per atom.
- src/petch/chemistry.py:10–34 — `_ied_yield`: petch's own statement of the Jensen effect
  and its bimodal-IED integration; the mixed layer should reuse this discipline.
