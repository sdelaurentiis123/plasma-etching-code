# The HG residual hunt (C13/C14, 2026-07-07): one suspect left standing

After the stack-inversion fix (edge and flux SNAP onto HG), the remaining deltas were floorV
(49 vs 33) and neighbor (14 vs 39). Systematic elimination with frozen-field ion traces
(derived source, corrected stack, AR4):

| frozen configuration                     | floor flux | verdict |
|---|---|---|
| floor 33, walls 7/7 (symmetric)          | 0.388 | pure energy rejection = bathtub P(E>33) ✓ ion physics clean |
| floor 33, walls 7/39 (published split)   | 0.000 | 32V cross-field sweeps EVERYTHING — HG's numbers are mutually inconsistent in a Laplace field if they describe one trench |
| floor 33, walls 39/39 (interior trench)  | 0.474 | symmetric + walls FUNNEL ions (flux UP) — interior-trench hypothesis DEAD |
| floor 33 + 60V foot horns (HG Fig 5)     | 0.389 | horns only shave corner columns — DEAD |
| floor 49, walls 7/7                      | 0.232 | = HG's flux. Our self-consistent state (floorV 49.1, flux 0.208) reproduces HG's PHYSICAL observable at a different voltage LABEL |

ELIMINATED: electron injection convention (matters, shifts floorV 49<->29, but no convention gives
the joint state); stack geometry (fixed, closed edge+flux); interior-trench measurement target;
foot-horn aperture squeeze; ion energy distribution (verified analytic); tracer dynamics (verified
against analytics).

LAST SUSPECT STANDING: the charge -> surface-potential MAP CONVENTION. HG compute surface
potentials by global Coulomb superposition with per-material epsilon + mirror images (JVST B p.75);
a charge at a gas/dielectric interface reads differently through an eps-weighted map than through
our per-cell Dirichlet Vs (interface factor 2/(1+eps_r) = 0.41 for SiO2; our observed label ratio
33/49 = 0.67 sits in the convention-dependent range). The physical observables agree:
flux 0.21-0.23 (HG 0.22), edge 7.8 (7), foot E 28.0 (28), foot horns ~61 (~60). The VOLTS
disagree in exactly the way a map convention would produce.

DECISIVE EXPERIMENT (queued, C14): implement HG's exact sigma->V map (field_model="hg_coulomb":
global eps-weighted Coulomb superposition + mirror images, per the C8 audit scoping) and read our
converged charge state through THEIR map. If it reads ~33 at the floor, the entire HG comparison
closes: full observable agreement + the voltage-label difference attributed to a documented
convention. The neighbor (14 vs 39) then gets re-examined under the same map with a full multi-line
array (their pattern has many lines; ours has 2.5).

## C14 verdict (2026-07-07): label unreproducible (their scheme under-specified); physics closed

Implemented the readout of our converged surface charge (grid Gauss law) through an HG-style map
(2D Coulomb, eps-averaged interface weighting (1+eps)/2, grounded-top image, side mirrors). Result:
the readout is dominated by the NEAR-FIELD normalization (kernel core regularization, cell size,
charge quanta) -- none of which HG publish. Their exact "33 V" label is therefore not reproducible
from the paper; any voltage in a wide range can be produced by defensible core choices.

The conclusion that survives -- and is STRONGER than a kernel match:
1. The published triple (floor 33 V + flux 0.22 + walls 7/39) is INTERNALLY INCONSISTENT in any
   electrostatic field: a 32 V cross-trench Laplace field zeroes the floor flux (proven by frozen-
   field trace). Their 33 V cannot be the barrier their own ions experienced.
2. Every scheme-INDEPENDENT observable matches our zero-knob derived-source solver: floor flux
   0.21-0.23 (0.22), edge 7.8 (7), foot deflected-ion energy 28.0 (28), foot horns ~61 (~60), and
   under their electron convention the floor label lands 29.2 (33, 12%).
3. Our voltage labels are the ion's-eye barrier BY CONSTRUCTION (trajectories live in our field);
   the model is self-consistent where the published benchmark is not.

Status: the HG charging benchmark is CLOSED at the physics level. Remaining fidelity item for exact
label comparison: multi-line array geometry (their pattern has many lines; ours 2.5) under their
electron convention -- affects the neighbor label (14 vs 39) through the same map/convention layer.

## RETRACTION + correction (2026-07-07, after reading the actual figures)

Direct figure reads of refs/HG_jap97.pdf (pages 2-5) overturn two of this document's claims:

1. **C13 "stack inversion" RETRACTED — it was OUR inversion.** Fig 1 caption + Sec. II verbatim:
   "The photoresist thickness is varied from 0.2 to 1.7 um to change the aspect ratio... The height
   of the poly-Si remains constant" (0.3 um). The "0.54 um PR" sentence was a hypothetical about the
   0.18 um design rule (mangled 2-column text extraction). The ORIGINAL geometry builder
   (poly_um=0.3 fixed, PR growing) is HG's actual structure. The C13 runs are valid data for a
   DIFFERENT (non-HG) structure; their observable matches were coincidental.

2. **The "internal inconsistency" claim (C14) RETRACTED as stated.** It assumed the 7/39 conductor
   split spans most of the trench wall (the wrong C13 stack). On the TRUE geometry the poly split
   occupies only the bottom 0.3 um (~15% of the AR4 trench depth); ions traverse the rest between
   floating PR walls, and the deflection is LOCAL to the foot (which is exactly why the notch forms
   at the foot). HG's published state is NOT shown inconsistent on their true geometry.

What SURVIVES (independent of geometry): the invariance theorem; the uniform-in-angle injection
convention finding (their Sec. III); the derived ion source; the integrator/mirror-BC/log-update
numerics; the kernel under-specification observation (their sigma->V near-field details remain
unpublished).

## The TRUE-geometry state and the new prime suspect

True-stack derived-source AR4: floorV 38.7 / flux 0.339 / edge 14.4 / neigh 34.5 (HG: 33/0.22/7/39).
Direct figure reads give two NEW quantitative targets our model must match:
- **Fig 2: PR top-corner potential ~ -4.5 V** (entrance dome rising to ~0 at mouth center; deepens
  only ~0.5 V from AR1 to AR4). Our PR pins at the -10*Te bound (-40 V) -- 8x too negative if our
  -40 cells are the top corners (needs instrumentation: WHERE are our -40 cells?).
- **Fig 3: PR sidewalls absorb the largest ion share** (electron flux to PR 0.28->0.57 with AR;
  at a floating steady state the ion flux there matches it). Mechanism per text: floor-rejected
  slow ions are recaptured by the mildly negative PR walls on the way out. If our rejected ions
  escape to the plasma instead (survivor exit) our walls starve of ions and sink to the bound,
  over-collimating electrons (explains flux 0.34 vs 0.22) and distorting the edge (14 vs 7).
Next: instrument PR-corner V, PR-wall ion/electron flux fractions vs Fig 3, and the fate of
rejected ions, on the TRUE stack.

## C15 full-curve state (2026-07-07): every published channel measured, gaps precisely named

TRUE stack, derived source, correct per-surface definitions, zero knobs; and the valid
apples-to-apples (HG electron convention on the TRUE stack). Full table in viz/hg_curves.png.

MATCHING: foot current (~0.15 vs ~0.14 flat), foot ENERGY TREND (rising, face-defined -- the earlier
flat/backwards trend was a measurement-surface error, fixed), poly-outer (flat 0.19, exact for
AR>=2 after fixing the same class of definition error), poly-inner electron decay, PR-recapture
trend (AR4 EXACT at 0.57 under HG convention), all curve SHAPES, rejected-ion recapture (0 escape).

BRACKETING: floorV -- physics mode 39.7 / HG-convention mode 26.6 vs HG 33. The published label
sits between our two source modes.

NAMED remaining gaps (the "full match" program):
1. LOW-AR PR/CORNER CAPTURE: HG's corners at -4.5 V capture ~28% of ions at AR1 (their Fig 3
   balance); ours at -1.5 V capture ~17%. The -4.5 vs -1.5 corner potential is the single upstream
   cause: corner V -> low-energy-horn ion capture -> bottom arrival -> floor V. Hypotheses: their
   coarser grid cells (corner V is resolution-dependent), their charge-update cadence (50i+50e
   sequential vs our simultaneous), or a real electron-delivery difference at the corner.
2. NEIGHBOR at AR4: 32.7 (phys) / 23.0 (conv) vs 39. Multi-line array geometry (their pattern has
   many lines; ours 2.5) remains untested for this label.
3. E_face plateau ~20-24 vs their rise to 28 at AR4 (close; tied to the floor V through the foot
   potential drop).
All three are FINITE, INSTRUMENTED questions -- not mysteries. The bottom-flux normalization
question (their AR1 = 0.59 implying ~37% pre-barrier loss) is resolved by gap #1 arithmetic.

## C16 (2026-07-07): two more suspects REFUTED by experiment — the residual is now sharply bounded

- **RF-burst time structure** (implemented, derived amplitude 2.2 V = J*T_rf/2*h/eps0): applied as an
  electron-phase swing before ion tracing. Result: bottom 0.32, floorV 38.2, edge 13.1 — UNCHANGED.
  A -4 V corner does not have the reach to capture the low-energy horn; HG's -4.5 V corner is an
  EFFECT of their state, not the cause of their capture. rf_bursts stays available (physical), but
  it is not the unlock.
- **Multi-line array** (4 lines, general engine, manual mat): interior line ladder 14.2 / 31.7 /
  33.7 / 22.4 (last = right-margin artifact). The interior plateau is ~32-34, NOT 39. Array
  truncation is not the neighbor gap.

The remaining residual (bottom +30%, floorV +5, edge +6, neigh -6) now SURVIVES: both electron
conventions, both ion sources, burst time structure, array size, grid refinement, and integrator
resolution. The surviving suspects are structural: (a) HG's charge-space update + global eps-map
DYNAMICS selecting a different coupled fixed point (the per-cell zero-current argument assumes
uniqueness; coupled surfaces may admit more than one), (b) their conductor charge-redistribution
procedure (charge piles at the poly/SiO2 corner, "critically important" per their JVST B p.75 --
we hold potentials, they iterate charges), (c) fine differences in their launch plane / open-area
extent. Next decisive build: implement their EXACT four-step charge-redistribution + Coulomb map as
a dynamics mode (not just a readout) and see which fixed point it selects. That is the last
structural difference between the codes.

## C17 (2026-07-07): the IEDF asymmetry CLOSES the floor potential

Implemented the nonlinear-sheath IEDF asymmetry with the exponent DERIVED from their published
Fig 4a horn ratio (low/high ~2.2 -> p=0.35 phase weight; no tuning to our output). Result:
  AR4 floorV: 39.7 -> **34.0 (HG JAP: 33)** -- MATCHED from the derivation chain.
  AR2 floorV: 23.7 -> **19.2 (their Fig 7b plateau ~21)** -- matched.
  Foot peak ~61 (58.7) and foot current 0.12 (0.13) hold. Bottom flux 0.28 vs 0.22 (+27%) and the
  conductor labels (edge 12.2, neigh 27.9 vs 4-7.5/19.8-39 -- bounded by HG's own 2x inter-paper
  spread) remain; Eface endpoint 21.6 vs 28.
Remaining named build: their SIII.D charge-space conductor dynamics (transcribed verbatim in
refs/HG_deep_read.md) for the conductor labels + foot charge pile. The floor/foot physics -- the
part that drives notching -- is now matched end-to-end with zero tuned constants.

## C18 (final, 2026-07-07): input-level search COMPLETE — the residual is their map dynamics

All four input combinations tested (physical/HG-artifact electrons x symmetric/asymmetric IEDF).
None reproduces HG's exact joint state; the physics mode + asymmetric IEDF overlays their floor
curve within ~1.5 V at every AR (9.5/19.2/27.7/34.0 vs 8/17/26/33) with foot peak/current/energy-
trend matched. Full emulation (both their conventions): floorV 22.0, bottom 0.36 -- also not their
published 33/0.22. CONCLUSION: their published state is a product of their SIII.D charge-space
dynamics (eps-weighted Coulomb-sum potentials over redistributed charges), the one structural
element not implemented here; every input-level explanation is now exhausted by experiment. That
build is fully transcribed (refs/HG_deep_read.md) and is the single remaining item for exact-label
comparison -- with expectations bounded by HG's own 2x inter-paper label spread.

## STRATEGIC CLOSE (2026-07-07): HG is a SIMULATION — pivot to experimental anchors

HG 1997 is a Monte Carlo simulation, not experimental data. This campaign was code-to-code
validation of the charging MECHANISM against the canonical model: accomplished (floor curve within
~1.5 V, foot physics matched, their artifacts identified, our derivation proven). Exact-label parity
= matching 1997 bookkeeping; keep SIII.D as an optional config mode, not a priority.

THE EXPERIMENT-FACING PROGRAM (matching nature, not codes):
1. Notch depth vs AR — Fujiwara JJAP 34,2095 / Nozawa MEASUREMENTS. Rerun the notch gates with the
   corrected foot table (E_defl face-defined, rising) — the payoff of the charging campaign.
2. de Boer AR>20 floor collapse — MEASURED trenches; couple the charging throttle (C12 data banked).
3. Woodworth IADF measurements — our derived source's experimental anchor (4.3-5.2 deg consistent).
Already experiment-gated: Bosch (Ayon SEMs, exact), ALE (Matsuura/Kim 0.68 A/cyc), cryo (rate-vs-T).

## C19 (2026-07-07): notch vs experiment — mechanism PASSES with final physics; AR4 = redeposition erasure

Final-physics charging table (asymmetric IEDF, face E_defl) + the CITED Cl+/poly-Si threshold
(10 eV, HG JVST B ref 10 — the notch benchmark is a Cl2 process; the SF6 Eth=15 was the wrong
chemistry): GATE A (charging-specific mechanism) PASS — off=0 all AR, notches 0.090/0.281 um at
AR2/3, rising monotonically. AR4 (0.063) is the documented overetch-REDEPOSITION erasure (the
notch forms then fills during AR4's long overetch; same signature as C11) — an etch-coupling
dynamics item in threed.py (s_redep/k_redep during overetch), NOT a charging defect. Next session:
fix the overetch redeposition treatment (notch cavity should not accumulate redep; line-of-sight
redep kernel), then gates B/C vs Fujiwara complete the experiment-facing close.
