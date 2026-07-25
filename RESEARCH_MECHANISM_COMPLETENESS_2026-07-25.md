# Mechanism Completeness Audit — petch mixed-layer vs Krüger 2024

Date: 2026-07-25
Scope: `src/petch/mixed_layer.py` + `src/petch/mixed_layer_mechanism.py` (the
element-resolved two-reservoir surface chemistry) audited against the complete
Krüger fluorocarbon SiO2/AC etch mechanism.
Directive: stop discovering missing physics one feature run at a time —
enumerate the COMPLETE published set now, so the next run is executed with the
full mechanism and any remaining miss is attributable to transport/spectrum
resolution, not missing chemistry.

## Sources consulted

- **Krüger, Zhang, Luan, Park, Metz, Kushner, JVST A 42, 043008 (2024)**,
  "Autonomous hybrid optimization of a SiO2 plasma etching mechanism", DOI
  10.1116/6.0003554. Local copy `tmp/pdfs/krueger-2024.pdf` (the 20-pp article).
  Read directly: Table II (complete reaction-class list), Eq. 2 (energy law),
  hot-neutral definition, reflection-off-polymer discussion, crosslink model
  (Fig. 5), O-etch product COFy, mask/AC arm.
- **Krüger PhD thesis (2024)**, "Modeling and Optimization of High Aspect Ratio
  Plasma Etching", DOI 10.7302/23106 — Appendix B (full surface mechanism) and
  Tables I/V. Not held locally as PDF; the exact numeric constants below are
  taken from the repo's own verbatim transcription in
  `src/petch/surface_kinetics.py` (`krueger_2024_reduced_projection`,
  `huang_2019_*` factory), whose provenance strings cite specific Appendix-B /
  Table entries. Constants not directly confirmable in the local article PDF are
  flagged **[VERIFY vs thesis]**.
- **Huang, Nam, Kushner, JVST A 37, 031304 (2019)** — the parent SiO2 mechanism
  Krüger inherits (Eq. 2 energy law; low-energy activation windows; hot neutrals;
  redeposition). Transcribed in `surface_kinetics.py` (`huang_2019` projection).
  Not held locally; **[VERIFY vs paper]** on the numeric activation constants.
- **Belen/Ertl SF6/O2 Si model** (ViennaPS) — `src/petch/belen.py`,
  `src/petch/surface_kinetics.py::SteinbruchelYield`. This is the Si-arm
  reference (sqrt-E threshold yield), not the SiO2/FC arm; used only to confirm
  the Steinbrüchel yield law form.
- Repo context: `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md` (current gate),
  `RESEARCH_MIXED_LAYER_DESIGN_2026-07-23.md`, `AUTONOMOUS_PROGRESS.md`
  (5–70/5–30 eV activation note; hot-neutral production still open).

## The energy law (Krüger Eq. 2 / Huang Eq. 2)

Energetic (ion + hot-neutral) reaction probabilities/yields follow

    p(ε) = p0 · (ε^q − ε_th^q) / (ε0^q − ε_th^q),   p = 0 for ε < ε_th

with per-channel reference probability `p0`, threshold `ε_th`, reference energy
`ε0`, and exponent `q` (called `n` in the repo). Angular dependence multiplies
this by a Kress-1999 form (`(1 + B(1−cos²θ))·cosθ`, B=9.3) or a Chang–Sawin form,
per channel. Thermal (sub-threshold) neutral reactions use a bare probability `p`.

**petch mixed-layer uses a DIFFERENT energy law.** `mixed_layer.py` derives a
deposited-energy ratio `energy_ratio = ε_dep / reference_energy_eV(=1000)` from
the ZBL stopping tables and scales every energetic channel *linearly* by it
(`sputter_total`, `volat_capacity`, `mix_total`). There is **no threshold** and
**no per-channel `q`/`ε_th`/`ε0`**. This is a deliberate "no fitted energy law"
choice, but it means none of Krüger's threshold behavior (the thing that gates
bare-SiO2 sputter off below 70 eV, and gives polymer sputter its sqrt-E shape
above 20 eV) is represented. Flagged in the table wherever it changes an
observable.

## Complete Krüger constant set (from repo transcription of Appendix B / Tables I & V)

SiO2 arm — `krueger_2024_reduced_projection`:

| channel | p0 | ε_th (eV) | q(n) | ε0 (eV) | angular |
|---|---|---|---|---|---|
| Bare SiO2 physical sputter | 0.0909 | 70 | 1 | 140 | Kress B=9.3 |
| SiO2·CxFy complex sputter | 0.1384 | 35 | 1 | 140 | Chang–Sawin |
| Polymer (un-crosslinked) sputter | **0.9** | **20** | **0.5** | **500** | Kress B=9.3 |
| Complex formation (chemisorption) | CF/CF2 0.2729; CF3/C2F3 0.2; C2F4/C3F5/C3F6 0.001 | thermal | — | — | — |
| Polymer dep on polymer | CF/CF2/CF3 0.1; C2F3 0.03 | thermal | — | — | — |
| Polymer dep on bare substrate | CF 0.002; CF2 0.0015; CF3/C2F3 0.001 | thermal | — | — | — |
| O polymer etch (→ COFy) | 0.0628 | thermal | — | — | — |
| C3F4 | inert (in HPEM flux table, absent from surface mechanism) | — | — | — | — |

Huang-2019 activation states (parent mechanism, `huang_2019` projection) —
**[VERIFY vs paper]** on numerics:

| channel | p0 | window / ε_th | notes |
|---|---|---|---|
| Complex low-energy activation | 0.1 | 5–70 eV | LowEnergyActivationYield, Kress 9.3 |
| Polymer low-energy activation | 0.3 | 5–30 eV | LowEnergyActivationYield, Kress 9.3 |
| Energetic FC-ion polymer deposition on complex | 0.1 | 5–70 eV | CF+/CF2+/…/C4F7+ only; Ar+ excluded |
| Activated-site polymer sticking | 10× the unactivated P values | thermal | P* sticking |
| Bare SiO2 sputter (Huang set) | 0.9 | 70 eV, q=0.5, ε0=140 | Krüger re-tuned p0→0.0909, q→1 |

## THE TABLE — mechanism class × Krüger treatment × petch status × observable × priority

Legend for status: **IN** = implemented and correct; **WRONG** = implemented
with an incorrect constant/law; **PARTIAL** = present but incomplete; **ABSENT**.

| # | Mechanism class | Krüger treatment (exact) | petch mixed-layer status | Dominant observable | Priority |
|---|---|---|---|---|---|
| 1 | **Polymer/film physical sputter** | p0=**0.9**, ε_th=**20 eV**, q=**0.5**, ε0=500, Kress B=9.3 | **WRONG.** `film_sputter_yield=0.1384` (mixed_layer.py:56) — that is Krüger's *complex-SiO2* sputter constant, **not** polymer sputter. Also linear-in-ε_dep, no 20 eV threshold, no sqrt-E. Film is sputtered ~**6.5× too weakly** and with the wrong energy shape. | **Mouth equilibrium** (film removal at the lip); sidewall taper | **P0 — CRITICAL** |
| 2 | **Ion→hot-neutral generation (reflection + neutralization)** | Ions that collide with a surface are neutralized into "hot neutrals" that retain energy+angle and continue to sputter/etch; explicitly noted to reflect off polymer curvature and travel deeper (Fig. 7 discussion) | **ABSENT** in chemistry. `mixed_layer_mechanism.advance` compresses the IEAD to one flux-weighted mean (E, cosθ) per face and generates no re-emitted hot-neutral flux. A `GrazingSpecularIonReflection3D` model (p_grazing=0.95, exp=3, retention=0.90) exists in the *transport/charging* layer but is not wired to the etch-front as a reflected hot-neutral source. | **Mouth equilibrium** + **HAR depth saturation** (grazing flux funneled to the floor) | **P0 — CRITICAL** |
| 3 | **Bare-SiO2 physical sputter channel** | p0=0.0909, ε_th=**70 eV**, q=1, ε0=140, Kress | **ABSENT / mis-modeled.** Substrate removal in mixed_layer is `volat_capacity·θ_F` (an F-coverage-gated volatilization, `volatilization_yield=1.0` knob), with no 70 eV threshold and no separate physical-sputter path. The bare (uncomplexed) oxide sputter channel does not exist. | Etch-front rate at low F coverage; **selectivity**; power-sweep | **P1** |
| 4 | **Complex (SiO2·CxFy) sputter as its own state** | p0=0.1384, ε_th=**35 eV**, q=1, ε0=140, Chang–Sawin — a *lower-threshold* channel than bare | **PARTIAL.** mixed_layer has one lumped volatilization capacity; the low-threshold complex channel is effectively what the model represents but without the 35 eV threshold or Chang–Sawin angular form. The bare/complex two-state distinction (central to selectivity) is collapsed. | **Selectivity**; etch-front rate | **P1** |
| 5 | **Low-energy activation states (5–70 eV complex, 5–30 eV polymer)** | LowEnergyActivationYield p0=0.1 (5–70) and 0.3 (5–30); activated sites take 10× polymer sticking | **ABSENT.** No activation-state channel; sticking is not energy-windowed. | Mouth equilibrium; sidewall passivation build-up | **P2** |
| 6 | **Energetic FC-ion direct polymer deposition** | p0=0.1, 5–70 eV, CFx+/CxFy+ family only (not Ar+) | **ABSENT.** Chemisorption channel exists (thermal complex formation) but not the *energetic* FC-ion deposition on complex sites. | Mouth equilibrium (extra deposition at the ion-rich lip) | **P2** |
| 7 | **Crosslink formation** | P(s)+P(s)→PC(s)+PC(s); dose-driven; PC is harder to sputter and lower radical sticking | **PARTIAL.** `n_xl_film` forms via ion-dose (mixed_layer.py:395–406) and lowers *sticking* (crosslinked_deposition p≈0.02). BUT crosslink state does **not** raise the sputter threshold/lower sputter yield — PC and P sputter identically. Krüger's "PC more resistant to sputtering" is not modeled. | Mouth equilibrium (differential fresh-vs-crosslinked at lip); sidewall taper | **P1** |
| 8 | **Crosslink breaking** | PC(s)+PC(s)+M(g)→P(s)+P(s)+M(g) by ions/hot-neutrals/photons | **ABSENT.** `n_xl` only increases (net of proportional removal); no ion/photon-driven de-crosslinking. | Mouth equilibrium dynamics; overetch behavior | **P3** |
| 9 | **SiO2 / etch-product redeposition** | S(s)+SiO2(g)→S(s)+SiO2(s) (Table II); sputtered products redeposit on other surfaces | **ABSENT** (declared omission in both files). | Sidewall taper; **HAR depth saturation**; bottom profile | **P2** |
| 10 | **O polymer etch → COFy** | P(s)+O(g)→COFy(g); only ground-state atomic O; single product COFy; controls necking/clog | **IN (finer).** mixed_layer oxidizes film C (p_ox=0.0628) and layer C, splitting CO/COF2 by local F/C — a *more* resolved product routing than Krüger's single COFy. Constant 0.0628 matches. | O-sweep; mouth equilibrium (clog control knob) | OK |
| 11 | **Complex formation (chemisorption)** | SiO2(s)+CxFy(g)→SiO2CxFy(s); CF/CF2 0.2729, CF3/C2F3 0.2, heavier 0.001 | **IN.** `KRUEGER_2024_CHEMISORPTION_PROBABILITY` matches verbatim (mixed_layer_mechanism.py:467). | Etch-front F delivery on open oxide | OK |
| 12 | **Substrate-dependent polymer deposition split** | on-polymer ~0.1 vs on-bare ~0.001–0.002 (100×) | **IN.** `KRUEGER_2024_DEPOSITION_ON_POLYMER/SUBSTRATE` verbatim. | Mouth narrowing 90→45; sidewall passivation | OK |
| 13 | **Film growth boundary motion** | polymer thickness moves the wall; necking/clog | **IN.** growth velocity = −d(film)/dt (mechanism adapter:291–300). | Mouth narrowing | OK |
| 14 | **Kress angular sputter response** | `(1+B(1−cos²θ))·cosθ`, B=9.3 | **IN** for film sputter (mixed_layer.py:262–264, B=9.3). Not applied to a bare-sputter channel (which is absent). | Sidewall/lip erosion angle | OK (partial coverage) |
| 15 | **Ion-mixing of film into layer** | (Humbird–Graves style; implicit in Krüger's excited-state transfer) | **IN** (`mixing_efficiency=1.0`, mixed_layer.py:275). | Etch-front F supply | OK |
| 16 | **AC mask chemistry (own arm)** | Full parallel mechanism: polymer dep on AC, AC sputter, O etch, crosslink | **PARTIAL.** `substrate="carbon"` path exists (CFx product, no lattice O) and shares the same FC constants, but AC-specific sputter yield/threshold not separately published-anchored. | Mask erosion; selectivity | **P3** |
| 17 | **Electron / photon channels** | Photons break crosslinks (Fig. 5e); no separate thermal electron etch channel in Table II | **ABSENT** (photon de-crosslink). No evidence Krüger's SiO2 arm needs an electron *etch* channel — electrons matter for charging, handled elsewhere. | Minor for neutral etch; relevant to charging | **P4 / out of scope** |
| 18 | **Atomic-F thermal etch** | Krüger's boundary publishes **no** atomic-F wafer flux; F reaches the surface bound in CxFy | **IN by design** (fluorine_species=() for the Krüger set). Correct omission, matches Krüger. | — | OK |

## Deep dive — the current gate (mask mouth equilibrium at 45 nm, landing ~30–40)

The campaign doc's own diagnosis (item 2–3): "Deposition at the lip outruns film
sputter ~50:1 at published fluxes, so nothing stops the closure." This audit
identifies the chemistry causes, in order of impact:

1. **The film sputter constant is the wrong published number (class #1).**
   `film_sputter_yield=0.1384` is Krüger's **complex-SiO2** sputter probability
   (`complex_sio2_yield`, p0=0.1384, ε_th=35 eV, Chang–Sawin), lifted into the
   film-sputter slot. The correct un-crosslinked **polymer** sputter constant is
   **p0=0.9, ε_th=20 eV, q=0.5, ε0=500, Kress B=9.3** (surface_kinetics.py:854–857,
   892–897, from Appendix B). At the ion-rich lip (ε well above 20 eV) the sqrt-E
   law puts the film-sputter yield near-maximal; 0.9 vs 0.1384 alone is ~6.5×
   more removal, and it is applied with the correct threshold shape rather than a
   linear ε/1000 ramp. **This is very likely most of the missing ~50:1.** Fixing
   only this may convert the 45 nm waypoint into an equilibrium.

2. **No hot-neutral / grazing-reflection re-flux (class #2).** Krüger explicitly
   reflects ions off the polymer lip as hot neutrals that (a) add sputter dose to
   the opposite lip and deeper sidewall and (b) funnel flux to the floor. Two
   consequences for the mouth: the lip film sees *more* removal than the compressed
   normal-incidence mean implies, and the equilibrium mouth width is set by a
   grazing-enhanced erosion the current mean-(E,cosθ) compression cannot produce.
   The transport-layer `GrazingSpecularIonReflection3D` (p_grazing=0.95, angular
   exp 3, energy retention 0.90 — Helmer–Graves / microtrench-bounded) is the
   right kernel; it must be delivered to `mixed_layer_mechanism.advance` as a
   reflected hot-neutral population, not collapsed away.

3. **Crosslinked film is not sputter-resistant (class #7, second half).** The
   model already lowers *sticking* on crosslinked film, but Krüger's crosslink
   also *raises* sputter resistance. Right now the lip's crosslinked skin is
   sputtered at the same (wrong, weak) yield as fresh film. Once #1 is fixed, the
   fresh-vs-crosslinked *sputter* differential is what pins the equilibrium: heavy
   ion dose at the lip → crosslinked, harder-to-sputter skin; ion-shadowed deep
   sidewall → fresh, easily-removed film. Without a PC/P sputter split the
   equilibrium mouth width is under-determined by chemistry.

4. **Low-energy activation windows and energetic FC-ion deposition (classes
   #5,#6)** add deposition precisely in the ion-rich mouth, and are the parent
   Huang mechanism's way of shaping the lip. Secondary to #1–#3 but part of
   "complete."

## VERDICT — minimal additions to make the mechanism COMPLETE relative to Krüger

Ranked. Each is a published constant, not a knob. After these, a remaining mouth
or depth miss is honestly attributable to **ion-spectrum compression** (the
declared mean-(E,cosθ) omission) or **transport resolution**, not to missing
surface chemistry.

**P0 — do before the next feature run (likely fixes the mouth by themselves):**

1. **Correct the film sputter constant + law (class #1).** Replace
   `film_sputter_yield=0.1384` / linear-ε_dep with the Appendix-B **polymer**
   sputter law: **p0=0.9, ε_th=20 eV, q=0.5, ε0=500 eV, Kress B=9.3**. Keep 0.1384
   only where a genuine complex-SiO2 sputter channel is added (class #4). This is
   the single highest-leverage change and corrects a confirmed mis-lift.

2. **Wire hot-neutral / grazing specular reflection into the surface chemistry
   (class #2).** Deliver the existing `GrazingSpecularIonReflection3D` output
   (**p_grazing=0.95, angular exponent 3, energy retention 0.90**) to
   `advance` as a reflected hot-neutral flux population, so lip and floor see the
   funneled energetic flux. This removes the mean-(E,cosθ) compression at the one
   place it most distorts the answer.

**P1 — needed for a complete, transferable mechanism (selectivity + mouth pin):**

3. **Two-state oxide sputter (classes #3,#4):** add the **bare-SiO2 physical
   sputter** channel (**p0=0.0909, ε_th=70 eV, q=1, ε0=140, Kress**) and the
   **complex sputter** channel (**p0=0.1384, ε_th=35 eV, q=1, ε0=140,
   Chang–Sawin**) as distinct thresholds, replacing the single knob-scaled
   `volatilization_yield` capacity. This is what makes selectivity emerge from
   thresholds rather than a fitted rate.

4. **Crosslinked-film sputter resistance (class #7):** make PC sputter with a
   higher threshold / lower p0 than fresh P (Krüger states PC is "more resistant
   to sputtering"; the exact PC yield is in Appendix B — **[VERIFY vs thesis]**
   for the numeric factor). Pairs with the already-present sticking reduction to
   pin the 45 nm equilibrium.

**P2 — complete the parent-mechanism inheritance:**

5. **Low-energy activation states (class #5):** complex 5–70 eV (p0=0.1) and
   polymer 5–30 eV (p0=0.3) activation yields; activated sites take 10× polymer
   sticking. **[VERIFY vs Huang 2019]** numerics (transcribed in
   `surface_kinetics.py::huang_2019`).

6. **Energetic FC-ion polymer deposition (class #6):** p0=0.1, 5–70 eV, FC+
   family only (Ar+ excluded).

7. **SiO2 / product redeposition (class #9):** S(s)+SiO2(g)→S(s)+SiO2(s). Matters
   for sidewall taper and HAR floor; currently a declared omission.

**P3 — closure of the model form:**

8. **Crosslink breaking (class #8):** PC+PC+M→P+P by ions/hot-neutrals/photons.

9. **AC-mask-specific sputter constants (class #16):** anchor the carbon-arm
   sputter yield to Krüger's AC values rather than reusing FC-film constants.

**P4 / out of scope for the neutral etch mechanism:**

10. **Photon de-crosslink / electron channels (class #17):** photons only enter
    via crosslink breaking (folds into #8); electrons belong to the charging
    solver, not this surface mechanism.

### Confirmed defect (headline)

`mixed_layer.py:56` `film_sputter_yield: float = 0.1384` **lifted the wrong
published constant** — 0.1384 is Krüger's SiO2·CxFy *complex* sputter probability
(`complex_sio2_yield`), not the polymer/film sputter yield. The correct film
sputter constant is **0.9 @ ε_th=20 eV, q=0.5, ε0=500 eV, Kress B=9.3**
(Appendix B, transcribed at `surface_kinetics.py:892–897`). The film is being
sputtered ~6.5× too weakly with the wrong energy law — the direct chemistry
cause of the ~50:1 deposition-over-sputter mouth clog.
