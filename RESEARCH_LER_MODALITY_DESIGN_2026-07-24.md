# LER / PSD-Roughness Modality — Research + Design

**Date:** 2026-07-24 · **Status:** design memo, not committed. Companion to
`RESEARCH_LER_EXPERIMENTAL_SOURCES_2026-07-21.md` (validation targets) and
`PROGRAM_ROADMAP_2026-07-24.md` Stage D.

**One-line thesis:** petch's continuum engine is smooth by construction; roughness is
fluctuation physics it cannot spontaneously produce. The right first product is not a
new stochastic solver but a **differentiable roughness *transfer function*** — propagate a
*measured* mask-edge PSD through the deterministic etch via the engine's own boundary
sensitivity, and add a physically-scaled intrinsic-noise term for the roughness the etch
itself injects. This composes exactly with the planned dual-number differentiability
(`RESEARCH_DIFFERENTIABLE_TRANSPORT_2026-07-21.md`) and honors the no-fitted-knobs doctrine:
`T(f)` is read from the engine, `PSD_in(f)` is measured resist, the intrinsic term is scaled
from countable fluctuation sources. **Nobody has a PSD-validated roughness engine** (§5 of the
sources doc) — this is an open lane, and doing it *blind* on a published transfer measurement
would be a first.

---

## 0. What "LER" means here, and why the engine can't currently make it

An etched line edge is a 1-D height signal `y(x)` along the line (top-down LER) or a 2-D
sidewall field `y(x,z)` (sidewall roughness). The community-standard description
(Mack 2018; Constantoudis 2018; imec roughness protocol) is the **self-affine triple plus the
full PSD**:

- **σ** — RMS edge deviation (unbiased 3σ LWR ≈ 3–5 nm post-litho EUV, dropping 20–50% post-etch).
- **ξ** — correlation length (~10–30 nm resist; grows through etch).
- **α (H)** — roughness/Hurst exponent, from the PSD log-log slope: **α = (slope_PSD − 1)/2**
  ([Constantoudis, Patsis, Leunissen, Gogolides, JVST B 22(4), 1974 (2004)](https://pubs.aip.org/avs/jvb/article-abstract/22/4/1974/589603/)).
- **PSD(f)** — the power spectral density, `f` from ~0.5/µm (box length) to ~50/µm (sampling),
  with `PSD(0)` the low-frequency plateau (∝ σ²ξ) and a high-frequency roll-off set by α.

petch's fields are C¹-smooth: the level-set `φ`, the radiosity flux, the mixed-layer element
ledgers are all continuous by construction, and the `analytic_occlusion` operator was purpose-built
to remove the last hard indicator (`RESEARCH_DIFFERENTIABLE_TRANSPORT` §0). There is no term in
the velocity that carries lateral fluctuation. Roughness must be **injected** (a mask-edge PSD
boundary condition) and **generated** (intrinsic etch stochastics), then **transported** by the
existing deterministic engine. The design question is *how* to add the injection + generation
without breaking smoothness, differentiability, GPU-shapeability, or the exact ledgers.

---

## 1. Literature map

### 1.1 The metrology substrate (must-match, else we validate against artifacts)
Fully covered in `RESEARCH_LER_EXPERIMENTAL_SOURCES` §1. Load-bearing facts for this design:
pre-2017 σ is +30–50% noise-biased; PSDs below ~10 nm wavelength are white-noise-floor dominated;
the **imec roughness protocol** ([Lorusso, JM3 17(4), 041009 (2018)](https://www.spiedigitallibrary.org/journals/journal-of-micro-nanopatterning-materials-and-metrology/volume-17/issue-04/041009/))
sets box ≥2 µm, mandatory noise-floor subtraction (Mack), PSD + autocorrelation reporting. **We
must ship a metrology-twin** (sample the simulated edge like a CD-SEM, add noise, run the Mack
estimator) before any comparison is meaningful.

### 1.2 Constantoudis / Gogolides / Kokkoris lineage (NCSR Demokritos) — the measure-AND-model group
This is the closest existing template to what petch wants, and it splits into **two model classes**
that map cleanly onto our architecture options:

- **Full 3-D stochastic/cellular Monte Carlo etch simulator** — cellular surface representation,
  MC mass transport + reaction kinetics, follows macroscopic profile evolution *and* roughening.
  Validated against their own σ/ξ/α before/after (this is architecture (c), the kMC-skin analog).
  Core LER-transfer paper: **Constantoudis, Kokkoris, Xydi, Gogolides, Pargon, Martin, "Line-edge-
  roughness transfer during plasma etching: modeling approaches and comparison with experimental
  results," JM3 8(4), 043004 (2009)**
  ([SPIE](https://www.spiedigitallibrary.org/journals/journal-of-micro-nanolithography-mems-and-moems/volume-8/issue-4/043004/) ·
  [HAL](https://dumas.ccsd.cnrs.fr/FMNT/hal-00399796v1)).
- **Geometrical (non-MC) transfer model** — **Constantoudis, Kokkoris, Gogolides, "Three-dimensional
  geometrical modeling of plasma transfer effects on line edge roughness," JM3 12(4), 041310 (2013)**
  ([SPIE](https://www.spiedigitallibrary.org/journals/journal-of-micro-nanolithography-mems-and-moems/volume-12/issue-4/041310/)).
  A *geometric* representation of the resist sidewall is transferred anisotropically and reproduces
  the main experimental trends **without a full stochastic solver**. This is a coarse ancestor of our
  proposed architecture (a): they transfer the geometry deterministically; petch improves it by
  transferring through the *actual physics* and reading the transfer as a differentiable sensitivity.
- **PSD/HHCF analysis methods** — the group *invented* the σ/ξ/α self-affine LER description
  ([Constantoudis, JVST B 22(4), 1974 (2004)](https://pubs.aip.org/avs/jvb/article-abstract/22/4/1974/589603/);
  review [Constantoudis, Papavieros, JM3 17(4), 041014 (2018)](https://www.spiedigitallibrary.org/journals/journal-of-micro-nanolithography-mems-and-moems/volume-17/issue-4/041014/)).
  Practical estimator note we must implement: **PSD method wins for α, HHCF wins for ξ** in the
  noisy/finite-N regime — so the metrology-twin should report *both*.
- Sidewall anisotropy: their "is the resist sidewall isotropic or anisotropic?" work measured
  **anisotropic vertical/horizontal correlation** — tells us how to seed 3-D sidewall noise (vertical
  striation direction ≠ line direction). `RESEARCH_LER_EXPERIMENTAL_SOURCES` §3.

### 1.3 Kushner group roughness work
- **Wang, Wang, Biolsi, Kushner, "Scaling of atomic layer etching of SiO2 in fluorocarbon plasmas:
  Transient etching and surface roughness," JVST A 39, 033003 (2021)**
  ([AIP](https://pubs.aip.org/avs/jva/article/39/3/033003/1079734/) · [OSTI](https://www.osti.gov/pages/biblio/1853526)).
  MCFPM prediction: **statistical cycle-to-cycle polymer residue → spatially-varying EPC → ALE
  roughness**. This is the mechanistic anchor for the *polymer-cluster discreteness* fluctuation
  source (§2). Simulation-only, no experimental PSD — the gap petch can fill.
- Krüger et al. JVST A 42, 043008 (2024) and waveform-tailoring companions: **verified to contain
  NO roughness data** (`RESEARCH_LER_EXPERIMENTAL_SOURCES` §2). The Krüger *blind protocol* transfers;
  the Krüger *dataset* does not. We reuse the freeze→sealed-holdout→reveal methodology only.

### 1.4 Sandia / imec / Leti resist-LER-transfer studies
- **imec ADI→AEI EUV PSD pairs** — single-step transfer functions on real stacks; the industry design
  rule *"litho owns low-frequency roughness, etch owns high-frequency (grows ξ)"* is the qualitative
  shape any `T(f)` must reproduce. Fractilia/Intel SPIE 2024: **material ranking flips between ADI and
  AEI** → you cannot predict AEI roughness without the etch model. (`RESEARCH_LER_EXPERIMENTAL_SOURCES` §2/§4).
- **CEA-Leti / Pargon gate-patterning series** — unbiased PSDs at *every* stack step (resist→BARC→
  hard-mask→poly-Si) across ≥5 pre-treatments. The richest per-step measured transfer functions
  (`RESEARCH_LER_EXPERIMENTAL_SOURCES` §2, ranked target #1).

### 1.5 kMC / etch-front-roughening theory (the intrinsic-noise anchor)
- **Drotar, Zhao, Lu, Wang, "Mechanisms for plasma and reactive ion etch-front roughening,"
  Phys. Rev. B 61, 3012 (2000)**
  ([APS](https://link.aps.org/doi/10.1103/PhysRevB.61.3012) ·
  [RG PDF](https://www.researchgate.net/publication/235508642_)).
  (2+1)-D flux-redistribution MC: etch-front roughening with re-emission gives **universal KPZ
  exponents α ≈ β ≈ z ≈ 1**, matching experiment. This fixes the *scaling* of the stochastic-forcing
  term — the intrinsic-roughness PSD slope and its temporal growth law are not free; they are the
  KPZ/re-emission class. Enumerated roughening mechanisms: stochastic noise, shadowing, etchant
  re-emission, micromasking, ion scattering.
- **Kuboi et al., "Insights into different etching properties of CW and ALE for SiO2 and Si3N4 using
  voxel-slab model," JVST A 37(5), 051004 (2019)**
  ([AIP](https://pubs.aip.org/avs/jva/article-abstract/37/5/051004/245895/)).
  Two-layer (C-F polymer + reactive layer) **voxel slab** on a continuum; residual-F and polymer-F
  statistics drive roughness. This is the reference design for architecture (c) — a discrete skin on
  a continuum bulk — and its two-layer decomposition **mirrors petch's own mixed-layer film/reaction
  split**, so a kMC skin would inherit the ledger structure rather than fight it.
- Martin & Cunge JVST B 26, 1281 (2008) (`sources` §2): **plasma can SMOOTH**; roughening needs
  micromasking. Hard constraint — any intrinsic term that grows roughness unconditionally is wrong.

### 1.6 Resist-side stochastics (how the input PSD is characterized / where it comes from)
- **Gallatin, "Resist blur and line edge roughness" (invited)**
  ([RG](https://www.researchgate.net/publication/286141743_)) and **Naulleau & Gallatin, "Line-edge
  roughness transfer function," Appl. Opt. 42(17) (2003)** ([LBNL PDF](https://goldberg.lbl.gov/papers/Naulleau_AO_42(17)2003.pdf)) —
  established the **frequency-domain LER transfer function**: `PSD_out(f) = |T(f)|² PSD_in(f) + PSD_intrinsic(f)`.
  This is literally the equation petch's modality implements; petch's contribution is computing `T(f)`
  from *first-principles etch physics* rather than a fitted Gaussian blur.
- **Gallatin, depth-dependence of LER ∝ photoacid diffusion length, JVST B 20(6), 2927 (2002)**
  ([AIP](https://pubs.aip.org/avs/jvb/article-abstract/20/6/2927/1074454/)); Hansen / IRM shot-noise
  stochastics: LER floor from **Poisson photon counting, σ_N/N = 1/√N**
  ([Hansen, JVST B 35(6), 061602 (2017)](https://pubs.aip.org/avs/jvb/article/35/6/061602/591639/)).
  These characterize the *input* edge — they give us the shape (σ/ξ/α) and physical origin of the
  resist PSD we inject, and the Poisson-counting logic reused for the etch-side ion-shot-noise term.

**Positioning:** Demokritos validated *trends* (geometrical + MC), Sawin validated *angle physics on
blanket surfaces* (dead lineage), Kushner has the *machinery* but validates *profiles not spectra*,
Gallatin/Naulleau gave the *transfer-function algebra* but for resist blur not etch physics. **No one
computes the etch `T(f)` from a validated first-principles transport+chemistry engine and predicts a
held-out post-etch PSD blind.** That is petch's opening, and the differentiability roadmap makes the
sensitivity computation nearly free.

---

## 2. Fluctuation sources, ranked for fluorocarbon SiO2/Si etch

Ranking = expected contribution to *post-etch sidewall* PSD in our chemistry. Each source is
parameterized from a measurable quantity — **no free knobs** (petch doctrine): every entry is
"derive / measure / declare as fab-measurable."

| Rank | Source | Physical scale | Parameterization (no knob) | Literature anchor |
|---|---|---|---|---|
| **1** | **Incoming mask-edge roughness (PSD transfer)** | Full band; dominates low-f (< ~1/50nm) | `PSD_in(f)` = *measured* resist PSD (σ,ξ,α); transferred by engine `T(f)`. Zero fitted content. | Naulleau/Gallatin AO 42 (2003); imec ADI→AEI; Pargon per-step PSDs |
| **2** | **Polymer cluster discreteness** | ξ set by CFx cluster/island size ~1–5 nm; sets high-f floor & micromasking | Film areal density `n_C_film,n_F_film` already in `mixed_layer.py`; cluster area from monolayer density `_MONOLAYER_AREAL_M2` → discreteness variance = 1/(cluster count per correlation cell). Measurable: FC film thickness vs bias (Oehrlein), CF2 sticking (Graves/Coburn) — already Stage-B constants | Wang/Kushner JVST A 39, 033003 (2021); Kuboi voxel-slab JVST A 37 (2019) |
| **3** | **Ion shot noise (flux statistics per correlation area)** | Poisson: `σ_N/N = 1/√N`, N = ions per correlation cell per etch time; ~atomic-to-nm amplitude | N = (ion flux Γ_i)·(cell area A_c)·(dwell t). Γ_i is an engine boundary input; A_c = (ξ_intrinsic)²; **no knob** — pure counting on quantities the engine already carries | Hansen JVST B 35 (2017) counting logic; Drotar/Wang PRB 61 (2000) flux-redistribution |
| **4** | **Stochastic sticking / reaction events** | Per-site Bernoulli on sticking `s_p`, oxidation `p_ox` (already constants in `MixedLayerParams`) | Variance = s(1−s)/N_sites per cell; `s_p=0.0842`, `p_ox=0.0628` already anchored (Krüger). Turns the *mean* rate the engine uses into a per-cell *variance* — reuses existing constants | `mixed_layer.py` constants; Wang/Kushner EPC-variance mechanism |
| **5** | **Redeposition granularity** | Discrete sputter-product landing; contributes micromasking bumps, ~1–3 nm | Redep flux already computed (`surface_product_redeposition_3d.py`); granularity variance = Poisson on landed-product count per cell. Gated OFF when redep flux → 0 (respects Martin-Cunge "no unconditional roughening") | Martin & Cunge JVST B 26 (2008); Nakazaki/Ono two-mode roughening |

**Reading:** Source 1 is the whole ballgame for the *first* validation (top-down LER transfer through
etch — it needs only `T(f)` + a measured input, no intrinsic model). Sources 2–5 are the additive
`PSD_intrinsic(f)` term — the roughness the etch *adds* (the high-f/ξ-growth part of the imec design
rule). They share one structural trick: **each converts a mean quantity the deterministic engine
already computes (flux, film density, sticking probability, redep flux) into a per-correlation-cell
variance via a counting/Bernoulli argument.** No new physical constant is introduced; the variance is
a function of quantities already in the ledgers. That is what keeps the modality knob-free.

---

## 3. Architecture options

Constraints to satisfy simultaneously: observables span **~1/µm to ~1/nm** (needs lateral resolution
≤ ~1 nm along the edge, but only in the roughness dimension — the profile dimension stays at the
engine's native 5–10 nm); **σ ~1–3 nm**; **GPU-shapeable**; **composes with exact ledgers**; **honors
differentiability**; **no fitted knobs**.

### (a) Perturbation / linear-response transfer function  ← RECOMMENDED SPINE
Treat the mask edge as a nominal straight boundary plus a small perturbation `δy_in(x) = Σ_k a_k e^{ikx}`.
Because petch is (or is becoming) differentiable, the map from boundary perturbation to final etched-edge
displacement is a **linear operator `T`** in the small-amplitude limit:
`δy_out(x) = ∫ T(x,x') δy_in(x') dx'`, translationally invariant along a straight line ⇒ diagonal in
Fourier: `δŷ_out(k) = T̂(k) δŷ_in(k)`. Then the measured transfer algebra falls straight out:
**`PSD_out(k) = |T̂(k)|² · PSD_in(k) + PSD_intrinsic(k)`** — exactly the Naulleau/Gallatin form, with
`T̂(k)` computed from the engine instead of fitted.

Two ways to get `T̂(k)`: (i) **finite-difference probe** — perturb the mask with a single-wavenumber
sinusoid of small amplitude, run the deterministic etch, read the output amplitude at the same `k`;
sweep `k`. Works *today*, no differentiability required. (ii) **dual-number / adjoint** — once
`RESEARCH_DIFFERENTIABLE_TRANSPORT` lands, `T̂(k)` is the Jacobian of the boundary→edge map, one solve.

- **Pros:** the first validation observable *is* a transfer function, so this option validates the exact
  thing the datasets measure; knob-free (`T` from engine, `PSD_in` measured); cheapest (deterministic
  solve + linear sweep, embarrassingly parallel over `k` on GPU); composes trivially with exact ledgers
  (chemistry runs unchanged, we only perturb the boundary); it is the *native* consumer of the planned
  differentiability — the same Jacobian that does inverse design gives `T(f)` for free; strictly better
  than the Demokritos 2013 *geometrical* transfer model because it transfers through real physics.
- **Cons:** linear ⇒ valid only for small σ/ξ (roughness ≪ CD, true for CMOS LER: σ~2nm ≪ CD~20nm, so
  fine; breaks for wiggling/line-collapse); cannot itself *generate* the intrinsic etch roughness
  (needs (b) for `PSD_intrinsic`); assumes translational invariance (fine for straight lines, needs
  block-Toeplitz treatment for corners/2-D layouts — deferred).

### (b) Stochastic forcing (physically-scaled noise per correlation cell)  ← RECOMMENDED, supplies `PSD_intrinsic`
Add a random field `η(x,z,t)` to the level-set normal velocity, with amplitude and correlation set by
the §2 sources (ion shot noise, polymer/redep discreteness) — *not* a fitted noise level. Per
correlation cell of area `A_c`, velocity variance = (per-event Δ)²·(event count variance)/(cell·dt),
where event counts are Poisson/Bernoulli on the fluxes and probabilities the engine already carries.
The **KPZ/re-emission scaling (Drotar/Wang, α≈β≈z≈1)** fixes how this forcing grows and correlates —
the spectral slope and temporal growth are predictions, not fits.

- **Pros:** generates the etch-*added* roughness (the piece (a) structurally cannot); scaling is
  literature-fixed (KPZ), so still knob-free; naturally produces ξ-growth and the high-f content; can be
  run as a *single* stochastic realization or as an *ensemble* to estimate `PSD_intrinsic` directly.
- **Cons:** each realization needs lateral resolution ≤ ξ_intrinsic (~1 nm) → expensive in 3-D unless
  restricted to a thin edge band; stochastic velocity fights the level-set smoothness/reinit (needs
  a band-limited noise + careful reinit cadence); ensemble averaging for a PSD is N_realizations×
  cost; **not differentiable** as-is (breaks the dual-number lane — so keep it *additive and separable*
  from (a), never inside the transferred path).

### (c) Hybrid discrete-continuum (kMC skin on continuum bulk, Kuboi-style)
A thin voxel/kMC slab (Kuboi voxel-slab; Kushner MCFPM skin) resolves discrete surface events in the
top few nm; the continuum level-set carries the bulk. The skin's two layers (C-F polymer + reactive)
**map 1:1 onto petch's existing mixed-layer film/reaction split**, so the ledgers transfer.

- **Pros:** most physically complete — captures micromasking, striations, the nonlinear
  wiggling/line-collapse regime (a) and (b) miss; directly comparable to Kushner/Kuboi/Demokritos MC.
- **Cons:** heaviest build (a whole kMC engine + continuum↔discrete coupling); slowest; hardest to make
  differentiable; biggest surface area for new knobs to sneak in. **Overkill for the first validation**
  (top-down PSD transfer), which is linear-regime physics.

### Recommendation
**Ship (a) as the spine + (b) as the additive intrinsic-noise source; defer (c).**

`PSD_out(f) = |T̂(f)|²·PSD_in(f) + PSD_intrinsic(f)` where **`T̂(f)` from (a)** (finite-difference probe
now, dual-number later) and **`PSD_intrinsic(f)` from (b)** (KPZ-scaled counting noise, ensemble or
analytic). This decomposition is not a compromise — it *is* the physical structure of the observable
(Naulleau/Gallatin) and it keeps the differentiable path (a) clean by pushing the non-differentiable
stochastic part (b) into a separable additive term. (c) becomes the roadmap escape hatch for the day a
dataset forces us into the nonlinear striation/micromasking regime, and its Kuboi two-layer structure
means it can be grafted onto the existing mixed-layer ledgers when that day comes. GPU-shape: (a) is a
parallel-over-`k` batch of the existing solver; (b) is a per-cell RNG on the edge band — both native.

---

## 4. Validation ladder (mapped to `RESEARCH_LER_EXPERIMENTAL_SOURCES` datasets)

Sequencing follows the sources doc §6 "practical sequencing," made concrete for the (a)+(b) architecture.

**Rung 0 — Metrology-twin (prerequisite, no physics claim).** Sample a *known synthetic* edge
(prescribed σ/ξ/α) as a CD-SEM would (imec protocol: box ≥2 µm, sampling interval, add white-noise
floor), run the PSD *and* HHCF estimators with Mack noise-floor subtraction, recover the input σ/ξ/α.
**Gate:** recover synthetic σ/ξ/α within the estimator's own bias band; PSD-method α and HHCF ξ agree
with ground truth. Boundary data: none needed (self-test). *Without this rung, every downstream number
is uninterpretable.*

**Rung 1 — Reproduce Demokritos+Leti JM3 8, 043004 (2009) (sources target #2).** Purpose-built
model-vs-experiment σ/ξ/α before/after in three scenarios: trim-only, transfer-only, trim+transfer.
This is the *first PSD-transfer reproduction*. Boundary data the dataset provides: input resist σ/ξ/α
(→ `PSD_in`), etch chemistry/step, output σ/ξ/α. **First blind-reproducible observable: does petch's
`|T̂(f)|²` reproduce the measured before→after σ/ξ/α in the transfer-only scenario?** Preregister:
(i) σ decreases on transfer (T<1 at high f), (ii) T→1 as f→0, (iii) ξ grows. **Gate:** all three
qualitative signs correct AND σ_out within a declared band (e.g. ±20%) *without tuning `T`* (it comes
from the engine). Digitize PSDs from figures [VERIFY exact σ/ξ/α values against the paper before gating].

**Rung 2 — Blind protocol on CEA-Leti/Pargon per-step series (sources target #1).** Unbiased PSDs at
every gate-stack step across ≥5 pre-treatments. Freeze the engine, seal a held-out subset of
treatments/steps, calibrate any *declared* boundary inputs (fluxes/IADF for that recipe) on the
non-held-out subset, then **blind-predict the held-out per-step output PSDs**. Boundary data: per-step
input PSD + recipe. **Gate (Krüger-style):** predicted held-out `PSD_out(f)` inside a preregistered
band across the frequency window, revealed once. Weakness to declare: Cl2/HBr Si chemistry (matches our
Si capability, not SiO2/FC) — so run Rung 2 on Si, keep FC for Rung 3.

**Rung 3 — Intrinsic-noise validation, decoupled from transfer.** Use Oehrlein C4F8/Ar blanket
resist-roughening (sources target #4, *our* FC chemistry, no transport confound) to test the (b) term
alone: the **"frozen spectrum, growing amplitude"** signature (spectrum set in ~seconds, then amplitude
grows) is a sharp falsifiable test of the noise statistics. Calibrate `PSD_intrinsic` scaling at one
condition, predict spectrum evolution at others. Cross-check the KPZ β≈1 growth law. **Gate:** predicted
amplitude-growth exponent + fixed spectral shape match at held-out conditions.

**Rung 4 — external 3-D geometry deliverable.** On an imported 3-D geometry, report
sidewall σ/ξ/α with the metrology-twin; spot-check absolute sidewall σ/ξ against **photonics
waveguide-sidewall PSDs** (arXiv 2105.11477 / 2501.11590, FC-etched Si/SiN, σ~0.5–3nm, ξ~30–100nm) and
**Goldfarb striation phenomenology**. This is a demonstration, not a blind gate (no matched input PSD),
but it exercises the full 3-D path through `feature_step_3d`.

Rungs are strictly ordered: 0 gates everything; 1 (reproduce) precedes 2 (blind); 3 constrains (b)
independently of (a) so a Rung-2 miss can be localized to transfer vs intrinsic.

---

## 5. First build slice (v1, ~1–2 sessions) — a falsifiable PSD prediction

**Goal:** the smallest thing that outputs `PSD_out(f)` for a real trench and can be checked against
Rung-0 + the start of Rung-1. This is **architecture (a), 2-D, finite-difference probe, no new physics,
no differentiability dependency, no intrinsic-noise term yet.**

**Scope:** a straight line/trench in the *existing* engine (`feature_step_3d` in a quasi-2-D / single
transverse slice configuration, or the 2-D `deterministic_exchange_2d` lane), mask edge treated as the
perturbable boundary.

**Build steps:**
1. **Metrology-twin module** (`src/petch/ler_metrology.py`, new) — synthetic self-affine edge generator
   (σ,ξ,α), CD-SEM sampler (imec box/interval + white-noise floor), PSD + HHCF estimators with Mack
   noise-floor subtraction, σ/ξ/α extraction (α = (slope−1)/2). **Passes Rung 0 on its own.** This is
   the piece with the most standalone value and zero engine coupling — build and gate it first.
2. **Transfer-function probe harness** (`src/petch/ler_transfer.py`, new) — given the engine's nominal
   trench solve, apply a single-wavenumber sinusoidal mask-edge perturbation of small amplitude `ε`
   (a few % of CD), run the deterministic etch, extract the output-edge displacement amplitude at that
   `k` → `T̂(k)`. Sweep `k` across the metrology window (0.5–50/µm), batched. Verify linearity by
   halving `ε` (amplitude should scale, phase stable) — a built-in knob-free self-check.
3. **Compose** `PSD_out(k) = |T̂(k)|²·PSD_in(k)` with a *measured* resist `PSD_in` (start with the
   Demokritos/Leti Rung-1 input PSD, [VERIFY digitized values]); output the predicted post-etch PSD and
   its σ/ξ/α via the twin.
4. **Falsifiable prediction:** the three preregistered transfer signs (T<1 high-f, T→1 low-f, ξ grows)
   plus a σ_out number. If the engine's deterministic etch does *not* produce a high-f roll-off in
   `T̂(k)` — i.e. it transfers all frequencies equally — that is a real, publishable negative result
   about the engine's lateral coupling, and it is falsified by the first probe sweep.

**Explicitly deferred from v1:** the intrinsic-noise term (b) — so v1 predicts transfer-only, correct
for the Rung-1 "transfer-only" scenario and honest about it; full 3-D sidewall (v1 is a single slice);
dual-number `T̂` (finite-difference is enough to prove the observable); block-Toeplitz for corners.

**Why this is the right first slice:** it produces a *number the datasets also produce* (`PSD_out` /
σ,ξ,α) from a knob-free computation (`T̂` from the engine, `PSD_in` measured), reuses the entire
existing deterministic solver untouched, has an internal linearity self-check, and its two modules
(metrology-twin, transfer-probe) are each independently useful and independently gateable. It sets up
Rung-1 reproduction directly and leaves a clean seam for (b) and dual-number `T̂` to slot in.

---

## 6. Open questions / risks to carry into the build

- **Translational-invariance assumption in (a):** valid for straight lines; an imported STL geometry may
  have corners/curvature → `T̂(k)` diagonal breaks, need block operator. Fine for v1 (straight line),
  flag for Rung 4.
- **Level-set + stochastic forcing (b) numerics:** band-limited noise vs reinit cadence is unsolved;
  prototype on a 1-D edge before 3-D. KPZ scaling gives the *target*, not the discretization.
- **Metrology self-consistency:** our simulated edge has no SEM noise; the datasets' *unbiased* PSDs
  already subtracted it. Decide once: compare in the unbiased domain (cleaner) and only forward-model SEM
  noise when comparing to *raw* (pre-2017) data. Prefer unbiased-domain comparison throughout.
- **[VERIFY] all digitized PSD values** from Demokritos/Leti JM3 8, 043004 (2009) and Pargon figures
  before any gate — the sources doc flags these as figure-digitized, not tabulated.
- **Doctrine watch:** the one place a knob could sneak in is `PSD_intrinsic` amplitude in (b). Guard it:
  the amplitude *must* be a counting/Bernoulli function of engine-carried fluxes/densities, checked
  against the KPZ growth law, never a free scalar.

---

### Appendix — locator status
Verified via web this session: Constantoudis JVST B 22(4) 1974 (2004) α-formula; Constantoudis/Kokkoris/
Gogolides JM3 12(4) 041310 (2013) geometrical model; Kuboi JVST A 37(5) 051004 (2019) voxel-slab; Wang/
Kushner JVST A 39 033003 (2021) ALE roughness; Drotar/Wang PRB 61 3012 (2000) KPZ α≈β≈z≈1; Naulleau/
Gallatin AO 42(17) (2003) transfer function; Hansen JVST B 35(6) 061602 (2017) shot-noise. Inherited
verified from `RESEARCH_LER_EXPERIMENTAL_SOURCES_2026-07-21.md`: all §1–§6 metrology + dataset locators.
Not independently re-verified here: exact numeric σ/ξ/α in the Demokritos/Leti and Pargon datasets —
marked [VERIFY] at point of use.
