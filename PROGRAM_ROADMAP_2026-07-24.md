# petch program roadmap — full-physics, any-etch, atomic-accuracy (staged 2026-07-24)

Doctrine (standing): no fitted knobs — derive, measure, or declare as fab-measurable
chemical constants; every capability lands through preregistered gates; everything
GPU-shaped; validated on real experiments or it doesn't count.

## Where the program stands

Exact deterministic exchange transport (analytic occlusion, receipted 1e-9);
blind-validated Krüger pipeline (freeze → sealed held-out → reveal, executed once);
K24-DEKNOB-1: yield-scale knob retired, power saturation emergent parameter-free;
mixed-layer element-ledger chemistry (6 constants, all anchored) wired through the
3-D engine behind `--surface-model mixed_layer`.

## Stages

### Stage A — mixed-layer at feature scale (NOW, days)
1. Fix the common_refinement sliver-face edge at near-zero surface velocity
   (probes running; fix + regression test).
2. Krüger base case on GPU at 10 nm: anchor k_v (one constant, base only,
   declared procedure), read the EMERGENT opening (the retired knob's old job).
3. If opening lands: rerun the 8-condition scorecard under mixed_layer —
   same declared criteria as K24-DEKNOB-1. If not: the residual localizes a
   missing mechanism (facet sputter at grazing incidence is the expected one);
   document, extend, re-gate.
Gate: scorecard section appended to KNOB_RETIREMENT_STUDY (or successor doc).

### Stage B — constants from beam data + spectrum fidelity (next)
1. Sawin/Yin 2008 sealed blind campaign (declared #1): calibrate one slice,
   blind-predict energies/angles/flux-ratios. Retires the yield-law constants
   as fitted quantities entirely.
2. Per-event ion spectrum integration in the mixed-layer adapter (removes the
   declared mean-energy compression; cached-table lookups over events).
3. Constants sweep from the thesis-mining catalog: CF2 sticking (Graves/Coburn),
   L_mix vs E (Humbird MD), film thickness vs bias (Oehrlein) — each becomes a
   pinned test asserting the model reproduces the measurement.
Gate: each constant either measured-and-matched or listed as fab-measurement target.

### Stage C — charging + deep AR (unlocks AR>20 claims)
1. GPU charging tracer (Route A: the MC charged-transport port; plan exists in
   CHARGING_PHYSICS_PLAN).
2. Sub-degree IADF handling (the de Boer AR>20 floor-collapse frontier).
3. 5 nm confirmation fix (guard scales with dx) + axis-graded AMR for holes.
Gate: de Boer high-AR floor + notching gates under one operator.

### Stage D — LER/PSD stochastic modality (design NOW, build after B)
The continuum engine is smooth by construction; LER is fluctuation physics.
Design doc in flight (Opus agent): seed measured mask-edge PSDs, add discrete
fluctuation sources (polymer clusters, ion shot noise), validate PSD *transfer
functions* against the Leti/Pargon series. Nobody has a PSD-validated roughness
simulation — open lane.
Gate: reproduce one published PSD-transfer measurement blind.

### Stage E — Tier-1 reactor model (design NOW, build parallel to B/C)
Global CCP/ICP model translating recipes (W, sccm, mTorr) into the fluxes+IADF
boundary petch consumes. Unlocks every industry dataset that publishes recipes
without boundaries (TEL HARC profiles, 3D-NAND charging, CD-SAXS). Design doc in
flight (Opus agent).
Gate: reproduce the published Krüger HPEM wafer fluxes from his recipe within a
declared band; then one TEL-recipe blind attempt.

Reactor tiers are a ladder to a full HPEM-class stack, not a substitute:
Tier-1 (0-D global, above) → Tier-2 (1-D/2-D axisymmetric fluid + EM power
deposition, GPU-native and differentiable from day one) → Tier-3 (kinetic
electron/ion MC modules, reusing the feature-side GPU tracing machinery).
Identical boundary contract at every tier; each tier gates against the one
below plus measured densities before it earns trust. Tier-1-first is a
measurement, not a compromise: Gate 1 against Krüger's known boundary
separates input error (cross sections, wall coefficients — shared by all
tiers) from dimensionality error (what Tier-2 buys), and that decides how
much full-solver to build.

### Stage F — differentiability + speed endgame
1. CUDA port of the analytic-occlusion exchange (10-50x expected).
2. Dual-number forward mode through transport+chemistry (boundary terms vanish
   identically per RESEARCH_DIFFERENTIABLE_TRANSPORT) → gradient calibration,
   then inverse design (target profile → recipe).
3. Chemistry expansion in mixed-layer form: Cl2/HBr poly-Si (Sawin DB has the
   tables), SiN, ALE module (Lam 2017 synergy tables), cryo-HF (new condensation
   physics — the one place new mechanism code is unavoidable).

## Dependency spine

A → (B, and A's scorecard gates D's build) ; B → C claims at depth ; E is
independent until its gate, then feeds every fab dataset ; F.1/F.2 can start
any time compute allows; F.3 rides on B's beam-data machinery.

## Standing constraints

Boxes: provision fast, kill when done. Blind credit only via fresh sealed
campaigns. Every stage ends in a committed gate artifact, not a vibe.
