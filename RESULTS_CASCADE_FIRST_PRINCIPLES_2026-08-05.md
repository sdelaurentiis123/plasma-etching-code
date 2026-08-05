# The cascade is correct, and ARDE is structurally impossible in it

Matched-beam audit of the reflection cascade (task: verify the AR-200 claims and
the cascade before anything ships).  Script: `scripts/matched_beam_cascade_grade.py`.
Data: `results/curated/cascade_matched_beam/grade.json`.

## 1. Eq. 2.34 is transcribed correctly

Huang thesis L2336-2341, verbatim:

> `E_s = E_i [(E_i - E_c)/(E_ts - E_c)] [(θ - θ_c)/(90° - θ_c)]`, θ > θ_c and
> E_c < E_i < E_ts ... "Incident particles with Ei > Ets are assumed to retain
> all of their energy.  Incident particles with θ < θc or Ei < Ec are assumed to
> diffusively scatter.  In the studies of etching of high aspect ratio features
> in oxide and ONO stacks in this thesis, Ets = 100 eV, Ec = 10 eV and θc = 70°."

`boundary_transport_3d.py:782-790` implements exactly this: prefactor `E_i`,
angle factor `(θ-θ_c)/(90-θ_c)`, energy factor `(E_i-E_c)/(E_ts-E_c)`, full
retention above `E_ts`, zero below either cutoff.  **No transcription defect.**

## 2. Huang's published fluxes cannot be used as a grading reference

Derivable before running anything.  For any *fixed* angular distribution and
straight walls, direct-ion acceptance of a slot of aspect ratio `A` with entry
uniform across the opening is `F(A) = E[max(0, 1 - A|tan b|)]`, which is
asymptotically `c/A`.  So for two aspect ratios the ratio is **bounded below by
their inverse ratio**.  His oxide AR 0 → 40 sits under a PR of AR 13, i.e.
total AR 13 → 53, giving

    F(53)/F(13) ≥ 13/53 = 0.245

His published decay (L5405-5407) is `0.3e15/2.0e15 = 0.150` — **below the
bound**.  A width sweep confirms it never reaches 0.150 (0.537 at σ=0.83°,
0.260 at 5°, 0.246 at 10°, railing at 0.238 by 30°).

His numbers therefore carry his context, which his own text supplies: the
profile tapers (L5443-5446), the PR erodes (L5586-5596), the feature is a via
not our trench, and every etch-front strike *including re-arrivals* increments
his flux count (L5399-5402).  Grading a straight-wall frozen-geometry scan
against them compares different observables.  Consistent with the funnelling
pass (`3489dbe`), now with a derivation rather than an observation.

## 3. The governing law, derived

In a straight trench with specular walls, a reflection off a vertical wall
flips the lateral velocity and preserves the axial one, so **the polar angle is
conserved along the cascade**.  Eq. 2.34 then returns the incident energy
unchanged for every bounce with `b < 20°` (θ_wall = 90-b > 70° and E > 100 eV),
so **energy is conserved too**.  The only attrition is the leftover-rule
reaction per bounce:

    bounces over aspect ratio AR:  n = AR·tan(b)
    per-bounce loss:               r = p₀·kress(90°-b),  kress → (1+B)·sin b
    ⟹  S(b, AR) = (1-r)ⁿ  →  exp( -p₀(1+B)·AR·b² ) = exp(-9.3·AR·b²)

Validated against the exact bounce product:

| b | AR | exact | exp(-9.3·AR·b²) | rel. err |
|---|---|---|---|---|
| 0.5° | 40 | 0.97099 | 0.97216 | 0.1% |
| 1.0° | 40 | 0.88410 | 0.89319 | 1.0% |
| 2.0° | 40 | 0.57972 | 0.63648 | 9.8% |

(The small-angle form degrades above ~2° where `kress` stops being linear; the
exact product is what the code uses.)

**Attenuation is quadratic in beam width and linear in AR.**  Nothing fitted.

## 4. Consequence: the observed anti-ARDE is correct behaviour, not a bug

Fraction of entering energetic flux still alive at the etch front:

| beam | AR 1 | AR 8 | AR 40 | AR 200 |
|---|---|---|---|---|
| Krüger digitised (σ 0.833°) | 0.998 | 0.983 | 0.922 | 0.734 |
| Kim core only (σ 0.380°) | 1.000 | 0.997 | 0.983 | 0.924 |
| Kim two-component, f=0.50 | 0.997 | 0.974 | 0.900 | 0.735 |
| Kim two-component, f=0.65 | 0.996 | 0.967 | 0.875 | 0.679 |

At the Krüger beam width, `9.3·40·(0.0145)² = 0.078` → **92% of energetic flux
survives to AR 40**.  A specular cascade in a straight feature under a narrow
beam delivers essentially everything that enters, at every depth.  Total
energetic delivery is therefore AR-flat by construction — which is exactly what
the trench scan measured (−3% to AR 16, `074068e`).

## 5. The real finding: ARDE is structurally impossible in this configuration

Two receipted facts combine:

- **Our chemistry is exactly ion-limited.**  `RESULTS_FLOOR_DELIVERY_2026-08-05.md`:
  oxide recession moves **under 1% across a 50× range of neutral delivery** and
  is **exactly proportional to ion flux** (7.54 / 15.10 / 22.78 nm/s at 1× / 2× / 3×).
- **Energetic delivery is AR-flat** (§3-4, derived and measured).

An ion-limited rate law fed by an AR-flat energetic supply **cannot produce
ARDE**.  The +27% rise from AR 0→4 is then geometric concentration on top of a
flat supply, not an anomaly needing its own explanation.

Meanwhile the physical process is in the other regime:

- Krüger's thesis flags his own ion-energy channel as likely overestimated
  (L4884-4888) and describes a neutral-transport-limited process.
- Huang: "thermal neutral species undergo diffusive scattering on the sidewalls
  are conductance limited in reaching the etch front as the AR increases"
  (L5399-5402).
- Our own hole study measured exactly that collapse: thermal delivery **0.656%**
  at AR 200 against 69-95% energetic (`HOLE_STUDY_RESULTS_2026-08-05.md`).

So the missing ARDE and the +29% depth overshoot are **one defect, not two**:
the surface model sits in the ion-limited regime where the real process is
neutral-limited.  Fixing the limiting regime is the single change that should
move both.  This is the same class as the de Boer knee (etchant-starved process
parameters), and it is a *chemistry-side supply* question, not a transport one.

## 6. What this does NOT license

- It does not license adopting a diffusive-scattering term to manufacture ARDE.
  Huang's cascade does carry one — "Reflections from the sidewalls can also
  include a stochastic diffusive component in the reflected velocity" (L5424-5427)
  — and ours does not (`boundary_transport_3d.py:801` is pure specular).  That is
  a real declared difference.  But §5 shows it is not what owns the depth gate,
  and no source tabulates its magnitude, so it stays `[VERIFY]`.
- It does not invalidate the hole study's transport numbers, which are graded
  against exact theory (Clausing 0.656% at AR 200), not against Huang.

## Status of the beam-width claim in the hole study

`HOLE_STUDY_RESULTS_2026-08-05.md` reports the two-component beam giving 5.7×
more ARDE at AR 200 than core-only.  That is an **acceptance/shadowing** effect
(flux that never enters the acceptance cone), which is separate from and much
larger than the **cascade attenuation** tabulated in §4.  Both are real; the
document should not conflate them.  Cross-check: core-only cascade survival at
AR 200 is 0.924 versus two-component 0.679-0.735, i.e. the tail adds ~20-25
points of cascade attenuation on top of its much larger shadowing effect.

## 7. Flag raised by ml21: the class-1 angular form now drives the oxide rows

`4c66df1` applied the Appendix-B angular markers, putting class 1
(`f = (1 + B sin²θ) cos θ`, B = 9.3) on the **SiO₂ bare-sputter row** — a
depth-setting channel.  That transcription is faithful to Krüger's own legend
(∠=1 → Kress 1999), but two facts from the earlier audits collide here:

- Kress 1999 is **Cu/Ar molecular dynamics** — the wrong material system
  (`RESEARCH_VERIFY_HUNT_2026-08-05.md`).
- The only *in-chemistry* angular measurements give peak/normal ≈ **1.3**
  (Cho 2000, CF₄/SiO₂; Schaepkens 1998, SiO₂ V-groove), where B = 9.3 gives
  **4.17**.

| θ | B = 9.3 | B = 1.7 (measurement-bounded) | ratio |
|---|---|---|---|
| 40° | 3.710 | 1.304 | 2.8× |
| 52.6° (peak) | 4.172 | 1.259 | 3.3× |
| 70° | 3.151 | 0.855 | 3.7× |
| 83° | 1.238 | 0.326 | 3.8× |

**Observed consequence.**  ml21 (first run carrying the classes) is integrating
at **dt ≈ 0.06 s against ml19's 0.15 s at matched simulated time — a 2.4×
smaller stable timestep**, with ~65 s of wall per step.  A 60 s endpoint needs
~950 steps ≈ 17 h, against ml19's 305 steps.  The forecast in `4c66df1`
correctly predicted a small effect on the *floor* (+7.3 %, the front sits at
0.92° tilt where both classes are ≈1), but said nothing about **sloped faces**,
where a 3-4× over-peaked yield raises surface velocity and tightens CFL.

**Status.**  Not reverted — it is source-faithful and the gates it shipped with
(normal-incidence invariance, ledger closure, scalar/atom agreement) all hold.
But it is now a *declared risk on the depth channel*, not a neutral
transcription, and it is the most likely cause of the timestep collapse.  The
cheap discriminating test is a frozen-geometry forecast of etch-front and
sidewall velocity under B = 9.3 versus a measurement-bounded B ≈ 1.7 on the
oxide rows only — no run required — before any further 60 s endpoint is bought.
