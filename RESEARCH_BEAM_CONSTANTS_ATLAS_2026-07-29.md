# Beam-constants atlas — which etch chemistries can be decked from measured data alone

Researched 2026-07-29. Question: *for the chemistry-deck framework
(`CHEMISTRY_DECK_DESIGN_2026-07-28.md`, `src/petch/chemistry_deck.py`), which etch
systems have enough **independently measured** surface-kinetics constants that a deck
can be built without adding fitted knobs — and which profile experiments remain free
for blind gates?*

Every citation below was fetched and verified this session unless marked `[VERIFY]`.
Four primary sources were downloaded in full and text-extracted; their numeric tables
are transcribed here (see §7 for local artifact paths).

---

## 0. Headline findings

**F1. "Beam-measured constants" is, in the strict sense, mostly a myth — and the
doctrine needs a finer grain than measured/fitted.** The canonical MIT beam decks are
*not* tables of directly measured rate constants. They are reaction sets whose
coefficients were **regressed against beam-measured etch yields under
mass-spectrometry-measured fluxes**. In Guo's Cl₂/poly-Si mixing-layer deck (Table 2.2,
below) the coefficient provenance column literally reads `Assu` / `Calc.` / `Fitted`,
and *Fitted* dominates (~20 of ~28 rows). The useful distinction is not
measured-vs-fitted but **what the constant was fitted TO**:

| Level | Definition | Blind-gate consequence |
|---|---|---|
| **L0** | Directly measured observable (spontaneous etch rate law, film thickness, product mass spectrum, EPC, surface binding energy) | Free; usable anywhere |
| **L1** | Regressed from **blanket beam yields** with independently measured ion/neutral fluxes (Chang 1998; Kwon TML; Guo 2009) | Free for **profile** gates — the fit never saw a feature |
| **L2** | Computed from an established semi-empirical formula (Yamamura/Bohdansky sputter yield, SRIM/ZBL, Eckstein) | Free; but the formula's own fit corpus must be declared |
| **L3** | **Fitted to feature profiles** (Belen SF₆/O₂; Krüger's five optimizer outputs) | **Contaminates profile gates.** A profile gate on an L3 deck is a self-consistency check, not a blind prediction |

This taxonomy is the single most useful output of this study and should be added to
the deck schema as a per-constant `provenance_level` field.

**F2. petch's "validated Si arm" is L3, and this is stated in the source's own
abstract.** Belen *et al.*, *Feature-scale model of Si etching in SF₆ plasma and
comparison with experiments*, **J. Vac. Sci. Technol. A 23, 99–113 (2005)**, DOI
[10.1116/1.1830495](https://doi.org/10.1116/1.1830495): "experimentally inaccessible
parameters such as the **F sticking coefficient, chemical etch rate constant, and the
ion-enhanced etch yield are determined by matching simulated feature profiles with
those obtained from carefully designed etching experiments." (verified via two
independent renderings of the AIP abstract; AIP blocks direct fetch — `[VERIFY]`
verbatim wording against the PDF). Consequence: the de Boer Figure-9 gate is **not
structurally blind** on those three constants — it is a transfer test across
width/time for a profile-fitted constitutive law. That is still a real and useful
result (it is a *transfer* gate), but it must not be described as a blind
first-principles prediction, and the same three constants are exactly the ones a beam
deck would replace.

**F3. Cl₂/poly-Si is the only system in the world with a complete L1 constant set and
an *independent* profile corpus from the same instrument lineage.** Chang's MIT thesis
(1998) tabulates the regressed sticking (`s`) and ion-enhanced reaction probability
(`β`) at three ion energies for both Ar⁺/Cl and Cl⁺/Cl, plus physical-sputter `Y₀`,
plus a measured product-redeposition suppression experiment (SiCl₂ beam) — and the
Mahorowala/Sawin four-part JVST B 20 (2002) series supplies a *three-variable
fractional-factorial* profile set (Cl₂ **and** HBr, photoresist **and** oxide masks)
that no constant in the deck ever saw. That is the strongest blind-gate architecture
available for any chemistry.

**F4. There is an L1 alternative to the Krüger fluorocarbon deck.** Kwon's
translating-mixed-layer series (JVST A 24, 1906/1914/1920, 2006) + Yin's beam yields
(JVST A 26, 161, 2008) + Guo's 20-reaction SiO₂ table (thesis Table 4.1, 2009) form a
fluorocarbon/SiO₂ deck whose coefficients were fitted to **blanket QCM yields**, never
to a feature. Retiring Krüger's five profile-fitted constants therefore does not
require new measurements — it requires **swapping to the beam-regressed deck and
re-running the Krüger feature as a blind gate**. This is a far cheaper retirement path
than the one currently written in `chemistry_deck.py`'s docstring.

**F5. No curated open database of etch surface-kinetics constants exists.** LXCat is
electron/ion *gas-phase* scattering (Pitchford *et al.*, Plasma Process. Polym. 14,
1600098, 2017); its stated scope mentions plasma-surface interactions but no surface
rate-constant database is published there. The only compilations with database
character are physical-sputtering yield fits (Yamamura & Tawara, At. Data Nucl. Data
Tables **62**, 149–253, 1996, DOI [10.1006/adnd.1996.0005](https://doi.org/10.1006/adnd.1996.0005))
and the fusion-side IAEA/Eckstein sputtering compilations. **The MIT theses are the
database.** They are open on DSpace and machine-readable (post-2005 ones have clean
text layers).

---

## 1. What a petch deck actually requires

Mapping `chemistry_deck.py` sections to physical blocks that must be sourced:

| # | Deck block | Schema key | Physical content |
|---|---|---|---|
| B1 | Species / stoichiometry | `species` | which neutrals+ions reach the wafer, C/F (or Cl/Br/O) counts |
| B2 | Chemisorption on bare substrate | `chemisorption.bare` | reactive sticking of the etchant on clean material |
| B3 | Chemisorption on ion-activated substrate | `chemisorption.activated` | ion-assisted adsorption enhancement |
| B4 | Deposition / passivation sticking | `deposition.*` | polymer or oxide-passivant growth on film / substrate / mask |
| B5 | Substrate sputter law (bare) | `sputter_laws.*_bare` | `Y₀(E,θ)`: p₀, E_th, exponent, angular form |
| B6 | Substrate sputter law (reacted complex) | `sputter_laws.*_complex` | the low-threshold ion-enhanced channel |
| B7 | Mask sputter law | `sputter_laws.*_mask` | selectivity |
| B8 | Film/passivant sputter + state transitions | `sputter_laws.film_*`, `state_transitions` | crosslinking, defluorination, activation probability |
| B9 | Product identity & stoichiometry (F/Cl ledger) | `materials.*.products` | SiF₄ vs SiF₂; SiCl₂ vs SiCl₄; CO/CO₂ |
| B10 | Product redeposition / re-emission | (not yet in schema) | sticking of sputtered products on sidewalls |
| B11 | Spontaneous (ion-free) etch law | (not yet in schema) | Arrhenius `k₀ exp(−Ea/kT)·Γ` |
| B12 | Ion scattering / reflection constants | `scattering` | specular threshold, cutoff angle, cascade depth |
| B13 | Layer scalars | `layer` | displacement energy, mixing efficiency, densities |

B12 and B13 are largely **chemistry-independent physics** (Huang/MCFPM retention rule;
Yamamura/SRIM). Scoring below is over B1–B11.

---

## 2. The Sawin / MIT beam corpus

Herbert Sawin's group ran two generations of apparatus: (i) a **multiple-beam machine**
(independent F-atom, Cl-atom, Cl₂, Ar⁺, Cl⁺, SiCl₂ sources; 1993–1998), and (ii) an
**inductively-coupled-plasma beam** producing a realistic ion+neutral mixture with
mass-spectrometric flux quantification and QCM/ellipsometric yield readout (2001–2009).
Both are the reference sources of `Y(E, θ, N/I)` surfaces.

### 2.1 Papers (all verified via Crossref DOI resolution)

| Ref | Citation | System | What is tabulated |
|---|---|---|---|
| **S-A** | D. C. Gray, I. Tepermeister, H. H. Sawin, *Phenomenological modeling of ion-enhanced surface kinetics in fluorine-based plasma etching*, **J. Vac. Sci. Technol. B 11, 1243–1257 (1993)**, DOI [10.1116/1.586925](https://doi.org/10.1116/1.586925) | F + Ar⁺ on undoped poly-Si and SiO₂ | Multiple-beam, F and Ar⁺ fluxes varied **independently over several orders of magnitude**; XPS of the etching surface; ion-enhancement quantified over a wide energy range (Kaufman + ECR sources). Reported threshold for Ar⁺-enhanced F etching ≈ **4 eV**; physical-sputter threshold of Si by Ar⁺ ≈ **20 eV** (both as cited by Chang 1998 §3.6) |
| **S-B** | J. P. Chang, H. H. Sawin, *Kinetic study of low energy ion-enhanced polysilicon etching using Cl, Cl₂, and Cl⁺ beam scattering*, **J. Vac. Sci. Technol. A 15, 610–615 (1997)**, DOI [10.1116/1.580692](https://doi.org/10.1116/1.580692) | Cl⁺/Cl/Cl₂ on poly-Si | `Y` vs ion energy, ion flux, **N/I ratio**, and impingement angle. `Y ∝ (E^½ − E_th^½)`, **E_th ≈ 10 eV** for Cl⁺; Langmuir-type saturation of `Y` with N/I |
| **S-C** | J. P. Chang, J. C. Arnold, G. C. H. Zau, H.-S. Shin, H. H. Sawin, *Kinetic study of low energy argon ion-enhanced plasma etching of polysilicon with atomic/molecular chlorine*, **J. Vac. Sci. Technol. A 15, 1853–1863 (1997)**, DOI [10.1116/1.580652](https://doi.org/10.1116/1.580652) | Ar⁺ + Cl / Cl₂ on poly-Si | The Ar⁺ companion of S-B; the `s`, `β`, `Y₀` triplet vs energy |
| **S-D** | J. P. Chang, A. P. Mahorowala, H. H. Sawin, *Plasma-surface kinetics and feature profile evolution in chlorine etching of polysilicon*, **J. Vac. Sci. Technol. A 16, 217–224 (1998)**, DOI [10.1116/1.580974](https://doi.org/10.1116/1.580974) | Cl₂/poly-Si | The beam-derived kinetics wired into an MC profile simulator |
| **S-E** | S. A. Vitale, H. Chae, H. H. Sawin, *Silicon etching yields in F₂, Cl₂, Br₂, and HBr high density plasmas*, **J. Vac. Sci. Technol. A 19, 2197–2206 (2001)**, DOI [10.1116/1.1378077](https://doi.org/10.1116/1.1378077) | Si in F₂, Cl₂, Br₂, HBr | `Y` vs **ion energy, ion angle, plasma composition**. Stated explicitly as "a database of experimental values needed for feature profile evolution modeling". Key results: `Y ∝ √E` for all four; **Cl₂ and HBr yields nearly identical** (HBr's lower Si etch rate is a *flux* effect, not a yield effect); Cl₂ `Y(θ)` falls rapidly above **60°**, HBr `Y(θ)` starts falling at small off-normal angles |
| **S-F** | W. Jin, S. A. Vitale, H. H. Sawin, *Plasma–surface kinetics and simulation of feature profile evolution in Cl₂+HBr etching of polysilicon*, **J. Vac. Sci. Technol. A 20, 2106–2114 (2002)**, DOI [10.1116/1.1517993](https://doi.org/10.1116/1.1517993) | Cl₂+HBr/poly-Si | Yields vs composition/energy/angle; **sticking coefficients derived for reactive neutrals in both Cl₂ and HBr**; E_th = **10 eV for both** |
| **S-G** | O. Kwon, H. H. Sawin, *Surface kinetics modeling of silicon and silicon oxide plasma etching. **I.** Effect of neutral and ion fluxes on etching yield of silicon oxide in fluorocarbon plasmas*, **J. Vac. Sci. Technol. A 24, 1906–1913 (2006)**, DOI [10.1116/1.2336225](https://doi.org/10.1116/1.2336225) | SiO₂ in C₄F₈, C₄F₈+80%Ar | Neutral/ion composition by QMS; yield by QCM; **atomic F flux identified as the dominant etch driver, C as the deposition driver** |
| **S-H** | O. Kwon, H. H. Sawin, *… **II.** Plasma etching surface kinetics modeling using translating mixed-layer representation*, **J. Vac. Sci. Technol. A 24, 1914–1919 (2006)**, DOI [10.1116/1.2336226](https://doi.org/10.1116/1.2336226) | (model) | The **TML formalism** — the direct ancestor of petch's `mixed_layer` |
| **S-I** | O. Kwon, B. Bai, H. H. Sawin, *… **III.** Modeling of silicon oxide etching in fluorocarbon chemistry using translating mixed-layer representation*, **J. Vac. Sci. Technol. A 24, 1920–1927 (2006)**, DOI [10.1116/1.2336227](https://doi.org/10.1116/1.2336227) | FC/SiO₂ | The applied deck; also applied to Si/Cl₂ and SiO₂/F |
| **S-J** | Y. Yin, H. H. Sawin, *Impact of etching kinetics on the roughening of thermal SiO₂ and low-k dielectric coral films in fluorocarbon plasmas*, **J. Vac. Sci. Technol. A 25, 802–811 (2007)**, DOI [10.1116/1.2748797](https://doi.org/10.1116/1.2748797) | SiO₂, coral in FC | Yield + roughening vs kinetics |
| **S-K** | Y. Yin, H. H. Sawin, *Surface roughening of silicon, thermal silicon dioxide, and low-k dielectric coral films in argon plasma*, **J. Vac. Sci. Technol. A 26, 151–160 (2008)**, DOI [10.1116/1.2821747](https://doi.org/10.1116/1.2821747) | Ar sputter | Pure-sputter reference arm |
| **S-L** | Y. Yin, H. H. Sawin, *Angular etching yields of polysilicon and dielectric materials in Cl₂/Ar and fluorocarbon plasmas*, **J. Vac. Sci. Technol. A 26, 161–173 (2008)**, DOI [10.1116/1.2821750](https://doi.org/10.1116/1.2821750) | poly-Si (Cl₂/Ar); SiO₂ + coral (FC) | **The angular-yield database.** `Y(θ)` families vs ion energy, feed composition, source pressure. Verified finding: the **neutral-to-ion flux ratio is the primary factor** setting whether `Y(θ)` is sputter-like (peak ~60–70° off-normal, low pressure) or ion-enhanced-like |
| **S-M** | W. Guo, B. Bai, H. H. Sawin, *Mixing-layer kinetics model for plasma etching and the cellular realization in three-dimensional profile simulator*, **J. Vac. Sci. Technol. A 27, 388–403 (2009)**, DOI [10.1116/1.3085722](https://doi.org/10.1116/1.3085722) | Cl₂/poly-Si | The mixing-layer deck + 3D cellular MC; validated against poly-Si-in-Cl₂ beam data across conditions |
| **S-N** | A. P. Mahorowala, H. H. Sawin, R. Jones, A. H. Labun, *Etching of polysilicon in inductively coupled Cl₂ and HBr discharges.* **I** DOI [10.1116/1.1481866](https://doi.org/10.1116/1.1481866) (J. Vac. Sci. Technol. B **20**, 1055–1063, 2002); **II** [10.1116/1.1481867](https://doi.org/10.1116/1.1481867) (1064–1076); **III** [10.1116/1.1481868](https://doi.org/10.1116/1.1481868) (1077–1083); **IV** [10.1116/1.1481869](https://doi.org/10.1116/1.1481869) (1084–1095) | Cl₂ **and** HBr profiles | **The profile-gate corpus.** Part I is a three-variable (ICP power / RF bias / flow) fractional-factorial profile set on photoresist- and oxide-masked poly-Si. Part III: mask faceting, sidewall deposition, microtrenching. Part IV: feature charging |

### 2.2 Theses (open on MIT DSpace; these carry the full parameter tables)

| Ref | Thesis | Handle | Text quality |
|---|---|---|---|
| **T-1** | J. P. Chang, *Study of plasma-surface kinetics and simulation of feature profile evolution in chlorine etching of patterned polysilicon*, PhD MIT ChemE, 1998 (adv. Sawin) | [1721.1/50356](http://hdl.handle.net/1721.1/50356) | **Good OCR** — tables recoverable |
| **T-2** | W. Jin, *Study of plasma-surface kinetics and feature profile simulation of poly-silicon etching in Cl₂/HBr plasma*, PhD MIT ChemE, 2003 (adv. Sawin) | [1721.1/28357](http://hdl.handle.net/1721.1/28357) | **OCR unusable** (scanned microfiche). Abstract recoverable; a library/print copy is needed for the HBr sticking numbers |
| **T-3** | O. Kwon, *Surface kinetics modeling of silicon oxide etching in fluorocarbon plasmas*, PhD MIT MSE, 2004 (adv. Sawin) | [1721.1/28353](http://hdl.handle.net/1721.1/28353) | **OCR unusable** — use the three JVST A 24 papers instead |
| **T-4** | Y. Yin, *Etching kinetics and surface roughening of polysilicon and dielectric materials in inductively coupled plasma beams*, PhD MIT ChemE, 2007 (adv. Sawin) | [1721.1/38973](http://hdl.handle.net/1721.1/38973) | **Clean text.** Contains the mass-spec flux tables for Cl₂/Ar, C₄F₈/Ar, C₄F₈/O₂/Ar, C₂F₆/Ar |
| **T-5** | H. Kawai, *3-dimensional modeling and simulation of surface and sidewall roughening during plasma etching*, PhD MIT, 2008 (adv. Sawin) | [1721.1/43201](http://hdl.handle.net/1721.1/43201) | not pulled |
| **T-6** | W. Guo, *Kinetics modeling and 3-dimensional simulation of surface roughness during plasma etching*, PhD MIT, 2009 (adv. Sawin) | [1721.1/46600](http://hdl.handle.net/1721.1/46600) | **Clean text.** Carries **Table 2.1** (physical-sputter coefficients), **Table 2.2** (full Cl₂/Si reaction deck), **Table 3.1/4.2** (angular forms), **Table 3.2** (Cl₂/Ar composition vs N/I), **Table 4.1** (full C₄F₈/SiO₂ reaction deck), **Table 2.3** (surface-chlorination cross-check vs literature) |

> Note on §2.2 attribution: `RESEARCH_THESIS_MINING_2026-07-23.md` labels handle
> `1721.1/43201` as "Guo, Wei 2008". That is **wrong** — 43201 is Kawai (2008); Guo's
> thesis is 46600 (2009). Correct before citing.

### 2.3 Transcribed constants — Chang 1998 (T-1)

These are the deck-grade L1 numbers for **Cl₂/poly-Si and Cl₂/SiO₂**. Model form
(Table 3.3): overall Cl adsorption with coefficient `s`; ion-enhanced etch with
probability `β` producing SiCl₄; physical sputter `Y₀`; site balance
`θ = sR/(sR + 4β)`; total `Y = c(θ_inc)·[Y₀(1−θ) + βθ]`.

**Table 3.4 — Ar⁺ ion-enhanced poly-Si etching with Cl and Cl₂:**

| Ion / neutral | E_ion (eV) | Y₀ | s | β |
|---|---|---|---|---|
| Ar⁺, Cl | 35 | — | 0.02 | 1.06 |
| Ar⁺, Cl | 60 | 0.007 | 0.12 | 2.15 |
| Ar⁺, Cl | 100 | 0.07 | 0.30 | 3.59 |
| Ar⁺, Cl₂ | 100 | 0.07 | 0.07 | 0.83 |

Both `β` and `s` scale linearly with √E; the ion-enhanced threshold is **E_th = 16 eV**
(vs ~20 eV for physical sputter of Si; vs ~4 eV for Ar⁺-enhanced F etching per Gray).

**Table 5.2 — Cl⁺ ion-enhanced poly-Si etching with Cl:**

| E_ion (eV) | s | β |
|---|---|---|
| 35 | 0.18 | 1.14 |
| 60 | 0.32 | 2.42 |
| 100 | 0.45 | 3.61 |

**Table 4.2 — Ar⁺ ion-enhanced SiO₂ etching with Cl** (products assumed Si, O₂, SiCl₂):

| E_ion (eV) | Y₀ | β | s |
|---|---|---|---|
| 70 | 0.01 | 0.04 | 0.001 |
| 100 | 0.02 | 0.08 | 0.005 |

→ direct measured **selectivity** block: `s(SiO₂) ≈ 0.005` vs `s(poly-Si) ≈ 0.30` at
100 eV — two orders of magnitude, from one instrument.

**Table 3.1 — universal `Y(E) = A(E^½ − E_th^½)` for Ar⁺ on Si:**

| A (eV^−½) | E_th (eV) | Source |
|---|---|---|
| 0.04 | 35 | Chang (this work) |
| 0.04 | 27 | Harper 1981 |
| 0.04 | 31 | Tachi 1986 |
| 0.03 | 64 | Oostra 1988 |

> **Provenance flag for the Krüger deck.** `chemistry_deck.py` labels the SiO₂-complex
> sputter channel `[0.1471, 35, 1.0, 140, "chang_sawin"]` and
> `RESEARCH_MECHANISM_COMPLETENESS_2026-07-25.md` repeats "ε_th = **35 eV**,
> Chang–Sawin". The only 35 eV threshold in the Chang corpus found here is
> **Table 3.1's Ar⁺-on-*Si* physical-sputter threshold**, not an SiO₂-complex
> threshold. Chang's own SiO₂ work (Table 4.2) is measured at 70 and 100 eV and does
> not quote a 35 eV complex threshold. The attribution needs to be traced to Krüger's
> Appendix B primary reference before it can be called measured. `[VERIFY]`

**Table 3.2 — spontaneous (ion-free) etch of Si (block B11), both L0:**

| Neutral | Substrate | Source | k₀ (Å·cm²·s / (#·min)) | Ea (eV) |
|---|---|---|---|---|
| Cl | poly-Si (undoped) | Walker | 2.57 × 10⁻¹⁴ | 0.29 |
| F | Si⟨100⟩ | Flamm | 3.59 × 10⁻¹⁵ | 0.108 |

Chang's evaluation: at N/I = 1000 the spontaneous Cl yield is **0.006**, i.e. 2–3
orders below the ion-enhanced yield → **negligible for undoped poly-Si at 300 K**.
Original F source: D. L. Flamm, V. M. Donnelly, J. A. Mucha, *The reaction of fluorine
atoms with silicon*, **J. Appl. Phys. 52, 3633–3639 (1981)**, rate
`R_F(Si) = (2.91 ± 0.20)×10⁻¹² T^½ n_F e^{−0.108 eV/kT}` Å/min.

**Product redeposition (block B10), measured:** Chang 1998 abstract — a **SiCl₂ beam**
was used to quantify the effect of product deposition: "The presence of SiCl₂
significantly suppressed the etching of the polysilicon by both Cl⁺ and Cl/Cl⁺, by
almost an order of magnitude." This is a direct L0/L1 measurement of the redeposition
block that petch's schema does not yet have a slot for.

**Honest note from the same thesis:** the profile simulator reproduced measured
profiles "with **one adjustable parameter, the chlorine recombination probability on
photoresist surfaces**." That is the single L3 constant in the entire Cl₂ system — and
it is a *mask* constant, not a substrate constant.

### 2.4 Transcribed constants — Guo 2009 (T-6)

**Table 2.1 — physical sputtering (L2, closed-form):**

- `A = 0.0054 (Z_p Z_t)^{1/2} (M_p/M_t)^{...} − 0.0198` (refs 11,12 of the thesis)
- `E_th = 25.2 (M_t/M_p)^{−0.6} + 0.928 (M_t/M_p)` (ref 13)
- `f(θ) = −81.70cos⁵θ + 224.03cos⁴θ − 208.19cos³θ + 67.569cos²θ − 0.711cosθ − 0.0242`
- Validation quoted: analytic yield for Cl⁺ on Si at 1 keV = **0.93 Si/ion** vs SRIM
  **0.97 Si/ion**.

**Table 2.2 — Cl₂/poly-Si mixing-layer deck, provenance column verbatim:**

| Reaction | Coef. class | value | E_th (eV) |
|---|---|---|---|
| Cl⁺(g) → Cl(s) | **Assumed** | 1 | — |
| Cl₂⁺(g) → 2Cl(s) | **Assumed** | 1 | — |
| Cl(g) → Cl(s) (sticking on Si vacancy) | **Fitted** | **0.75** | — |
| Cl₂(g) → 2Cl(s) | **Fitted** | **0.204** | — |
| Si(s)→Si(g) by Ar⁺ / Cl⁺ / Cl₂⁺ | **Calc.** | 0.035 / 0.035 / 0.042 | 33.63 / 31.49 / 46.94 |
| Cl(s)→Cl(g) by Ar⁺ / Cl⁺ / Cl₂⁺ | **Calc.** | 0.045 / 0.045 / 0.055 | 29.44 / 27.62 / 40.86 |
| vacancy generation by Ar⁺ / Cl⁺ / Cl₂⁺ | **Fitted** | 1.8 / 2.1 / 0.598 | 27 / 20.2 / 32.6 |
| vacancy annihilation by Ar⁺ / Cl⁺ / Cl₂⁺ | **Fitted** | 10.0 / 8.459 / 0.001 | 27.0 / 0.0 / 30 |
| Si + 2Cl → SiCl₂ by Ar⁺ / Cl⁺ / Cl₂⁺ | **Fitted** | 8.30 / 7.41 / 3.60 | 26.4 / 25.1 / 33.5 |
| Cl + Cl → Cl₂ by Ar⁺ / Cl⁺ / Cl₂⁺ | **Fitted** | 5.3 / 5.06 / 6.26 | 26.4 / 0.0 / 0 |
| Si-V + Si-V → Si-Si (3 ions); Cl-V + Cl-V (3 ions) | **Fitted** | — | — |

Count: **2 assumed, 6 calculated, ~20 fitted** — all fitted *to beam yields*, i.e. L1.
The dominant ion-induced product is **SiCl₂** (not SiCl₄ as in Chang's simpler 1998
model) — a real, citable disagreement between the two MIT generations on block B9.

**Table 2.3** cross-checks the modelled steady-state surface Si:Cl fractions against
five independent literature sources over 25–150 eV and N/I = 0–400 (e.g. at 50 eV,
N/I = 0: model 0.80/0.20 vs literature 0.82/0.18). This is a ready-made **held-out
composition gate** for petch's mixed-layer state.

**Table 3.2 — measured Cl₂/Ar beam composition vs N/I** (mass spec):

| species fraction | N/I = 3.5 | 20 | 131 |
|---|---|---|---|
| Cl (neutral) | 0.63 | 0.30 | 0.23 |
| Cl₂ (neutral) | 0.37 | 0.70 | 0.77 |
| Cl⁺ | 0.40 | 0.69 | 0.54 |
| Cl₂⁺ | 0.082 | 0.054 | 0.42 |
| Ar⁺ | 0.522 | 0.257 | 0.04 |

**Table 4.1 — C₄F₈/Ar → SiO₂ deck (20 reactions).** Species set (verbatim):
ions `C⁺, O⁺, CO⁺, CF⁺, CF₂⁺, CF₃⁺, C₂F₄⁺, C₂F₅⁺, C₃F₅⁺, C₄F₇⁺`;
neutrals `C, O, CO, CF, CF₂, F, CF₃, C₂F₄, C₂F₅, C₃F₃, C₃F₅, C₄F₃, SiF, SiF₂, SiF₃`.
Coefficients "fitted to the etching yields of oxide in C₄F₈/Ar plasma at different
conditions **measured by Yin**" — i.e. **L1 against S-L/T-4**. Selected values:

| # | Reaction | class | A | E_th (eV) |
|---|---|---|---|---|
| 5–8 | physical sputter of Si / O / C / F | **Calc.** | 0.042 / 0.018 / 0.009 / 0.023 | — |
| 9 | Si(s)+2F(s) → SiF₂(g) ion-enhanced | fitted | 6.75 | — |
| 10 | 2F(s) → F₂(g) | fitted | 0 | 20 |
| 11 | 2O(s) → O₂(g) | fitted | 0.22 | — |
| 12 | Si(s)+O(s) → SiO(g) | fitted | 0.007 | 20 |
| 13 | C(s)+O(s) → CO(g) | fitted | 0.24 | — |
| 14 | C(s)+2O(s) → CO₂(g) | fitted | 0.95 | — |
| 15 | C(s)+2F(s) → CF₂(g) | fitted | 2 | 0 |
| 16 | vacancy creation by ion | fitted | 0.14 | — |
| 17 | densification | fitted | 1.66 | — |
| 18 | dangling-bond annihilation | fitted | 10000 | — |
| 19 | CF₃(s)+F(s) → CF₄(g) recombination | fitted | 0 | — |
| 20 | Si(s)+4F(s) → SiF₄(g) **spontaneous** | fitted | 2.99×10⁻⁵ | — |

(The sticking coefficients `S_I`, `S_F`, `S_N-on-C`, `S_N-on-O` appear as symbols in
rows 1–4; their numeric values were lost in the OCR of the parameter column and must
be re-read from the PDF page 100 or from S-I. `[VERIFY]`)

**Table 4.2 — angular forms**: ion incorporation `f(θ)=cosθ`; physical sputtering
`f(θ) = −141.29cos⁶θ + 641.11cos⁵θ − 1111.3cos⁴θ + 944.63cos³θ − 421.98cos²θ +
95.31cosθ − 5.46`; ion-enhanced etching flat (`f=1`) below ~25°.

---

## 3. Sources outside the Sawin lineage

### 3.1 Coburn–Winters legacy (IBM) — the L0 bedrock

- J. W. Coburn, H. F. Winters, *Ion- and electron-assisted gas-surface chemistry — an
  important effect in plasma etching*, **J. Appl. Phys. 50, 3189–3196 (1979)**, DOI
  [10.1063/1.326355](https://doi.org/10.1063/1.326355). The synergy experiment
  (XeF₂ + 450 eV Ar⁺ on Si, SiO₂, Si₃N₄; also F₂, Cl₂). Etch yields up to ~**25
  Si/ion** reported. The canonical *shape* constraint, not a rate-constant table.
- H. F. Winters, J. W. Coburn, *Surface science aspects of etching reactions*,
  **Surf. Sci. Rep. 14, 161 (1992)** — the review that Chang 1998 cites for "physical
  sputtering rate always drops in the presence of a chemically active gas", the
  physical basis of petch's `Y₀(1−θ)` bare-fraction weighting.
- Winters & Coburn, *Ion-enhanced gas-surface chemistry: the influence of the mass of
  the incident ion*, Surf. Sci. **102**, L1 (1981)-era; mass-scaling of the enhancement.

### 3.2 LLNL molecular-beam group (Levinson, Shaqfeh, Balooch, Hamza)

- *Ion-assisted etching and profile development of silicon in molecular chlorine*,
  **J. Vac. Sci. Technol. A 15, 1902–1912 (1997)**, DOI [10.1116/1.580658](https://doi.org/10.1116/1.580658)
- *Ion-assisted etching and profile development of silicon in molecular and atomic
  chlorine*, **J. Vac. Sci. Technol. B 18, 172–190 (2000)**, DOI [10.1116/1.591170](https://doi.org/10.1116/1.591170)

> **Correction to `CHEMISTRY_EXPANSION_PLAN.md` §W2:** these are attributed there to
> "Chang, Arnold & Sawin, JVST B 18, 172 (2000)". Both papers are by **Levinson,
> Shaqfeh, Balooch, Hamza (LLNL)**, and the 1997 one is JVST **A** 15, 1902. Fix
> before any deck cites them. Their independence from MIT is an *advantage*: they are a
> second, disjoint profile+beam corpus for Cl₂/Si.

Also relevant, independent: T. Kolfschoten, R. Haring, A. Haring, A. de Vries,
*Argon-ion assisted etching of silicon by molecular chlorine*, **J. Appl. Phys. 55,
3813 (1984)**, DOI [10.1063/1.332890](https://doi.org/10.1063/1.332890).

### 3.3 Oehrlein group (IBM → Maryland) — the fluorocarbon film state (block B8)

- T. E. F. M. Standaert, C. Hedlund, E. A. Joseph, G. S. Oehrlein, T. J. Dalton,
  *Role of fluorocarbon film formation in the etching of silicon, silicon dioxide,
  silicon nitride, and amorphous hydrogenated silicon carbide*, **J. Vac. Sci.
  Technol. A 22, 53–60 (2004)**, DOI [10.1116/1.1626642](https://doi.org/10.1116/1.1626642).
  Measured: steady-state FC film thickness (ellipsometry) + composition (XPS) on four
  materials vs self-bias. Verified headline: the correlation between etch rate and FC
  film *thickness* is **weak**; **ion-induced defluorination of the film** is the major
  channel — i.e. the film is a fluorine *source*, not merely an inhibitor.
- T. E. F. M. Standaert *et al.*, J. Vac. Sci. Technol. A **15**, 1881 (1997) — CHF₃
  ICP, the dep→etch transition vs bias (blanket-rate gate).
- D. Metzler, R. L. Bruce, S. Engelmann, E. A. Joseph, G. S. Oehrlein, *Fluorocarbon
  assisted atomic layer etching of SiO₂ using cyclic Ar/C₄F₈ plasma*, J. Vac. Sci.
  Technol. A **32**, 020603 (2014) and the 2016 follow-ups — measured FC thickness
  deposited per cycle, EPC vs ion energy, SiO₂/Si₃N₄ selectivity.
- CFx radical surface loss (block B4, L0): G. Cunge, J. P. Booth, *CF₂ production and
  loss mechanisms in fluorocarbon discharges: fluorine-poor conditions and
  polymerization*, **J. Appl. Phys. 85, 3952–3959 (1999)**, DOI [10.1063/1.370296](https://doi.org/10.1063/1.370296)
  (absolutely-calibrated LIF; CF₂ produced at *all* surfaces via a long-lived
  precursor; heavier CₓF_y radicals stick more than CF₂). Companion measurement of CF₂
  surface interaction during a-C:F deposition from CHF₃: **J. Appl. Phys. 84, 4736
  (1998)**. These are the only *directly measured* CFx sticking numbers in the field
  and they are reactor-averaged, not per-surface-state — a genuine gap.

### 3.4 LTM-Grenoble (Joubert, Cunge, Vallier, Pargon) — products & sidewall layers

- G. Cunge, B. Kogelschatz, N. Sadeghi, *Production and loss mechanisms of SiClₓ etch
  products during silicon etching in a high density HBr/Cl₂/O₂ plasma*, **J. Appl.
  Phys. 96, 4578–4587 (2004)**, DOI [10.1063/1.1786338](https://doi.org/10.1063/1.1786338) — block B9/B10.
- G. Cunge, M. Inglebert, O. Joubert, L. Vallier, N. Sadeghi, *Ion flux composition in
  HBr/Cl₂/O₂ and HBr/Cl₂/O₂/CF₄ chemistries during silicon etching in industrial
  high-density plasmas*, **J. Vac. Sci. Technol. B 20, 2137–2148 (2002)**, DOI
  [10.1116/1.1511219](https://doi.org/10.1116/1.1511219) — block B1 for HBr.
- L. Vallier, J. Foucher, X. Detter, E. Pargon *et al.*, *Chemical topography analyses
  of silicon gates etched in HBr/Cl₂/O₂ and HBr/Cl₂/O₂/CF₄ high density plasmas*,
  **J. Vac. Sci. Technol. B 21, 904–911 (2003)**, DOI [10.1116/1.1563255](https://doi.org/10.1116/1.1563255)
  — measured SiOₓBr_y sidewall passivation composition; block B4 for HBr.
- "Spinning wall" chamber-wall studies (Donnelly/Houston) — HBr ≈ 30% dissociated with
  no bias; SiBr₂/SiBr₄ dominate the neutral product flux while SiBr/SiBr₃ dominate the
  mass-spec ion signal. `[VERIFY]` exact citation before use.

### 3.5 Lam / Gottscho — the ALE-era measured set

- K. J. Kanarik, T. Lill, E. A. Hudson, S. Sriraman, S. Tan, J. Marks, V. Vahedi,
  R. A. Gottscho, *Overview of atomic layer etching in the semiconductor industry*,
  **J. Vac. Sci. Technol. A 33, 020802 (2015)**, DOI [10.1116/1.4913379](https://doi.org/10.1116/1.4913379).
- K. J. Kanarik *et al.*, *Predicting synergy in atomic layer etching*, **J. Vac. Sci.
  Technol. A 35, 05C302 (2017)** — **open PDF**:
  <https://www.osti.gov/servlets/purl/1376399>. Verified measured content: for six
  materials (Si, Ge, C, W, GaN, SiO₂) they report EPC, the per-step spurious rates
  `a` (dose-step etch) and `b` (ion-step sputter), the **ALE energy window edges**, and
  the **synergy** `S`. Quoted values found in the text: Si — window **40–60 eV**,
  `S = 90%` at 50 eV, `E₀ = 4.7 eV`; Ge — window **20–30 eV**, `a = 0.03 nm/cycle`,
  `E₀ = 3.8 eV`; C — window **35–75 eV**, `b = 0.01 nm/cycle` at 50 eV, `E₀ = 7.4 eV`;
  W — 60 eV, `E₀ = 8.9 eV`; GaN — window **50–90 eV**, `b = 0.03 nm/cycle` at 70 eV,
  `E₀ ≈ 8.6 eV`; **SiO₂ — `S = 80%` at 50 eV, EPC = 0.5 nm/cycle, CHF₃/Ar, `a ≈ 0.02
  nm/cycle`, `E₀ ≈ 5 eV`**. Surface binding energies `E₀` are taken from **Yamamura &
  Tawara (1996)**, i.e. independent of the etch measurement — a genuinely orthogonal
  anchor for petch's sputter thresholds.
- Kanarik *et al.*, *Atomic layer etching: rethinking the art of etch*, **J. Phys.
  Chem. Lett. 9, 4814 (2018)**, DOI [10.1021/acs.jpclett.8b00997](https://doi.org/10.1021/acs.jpclett.8b00997).

### 3.6 In-situ ALE surface chemistry (Agarwal/Colorado School of Mines, 2019–2025)

ATR-FTIR of the *actual surface layer* during cyclic etching — the closest thing to a
direct measurement of petch's mixed-layer state variable:
- R. J. Gasvoda *et al.*, *In situ monitoring of surface reactions during atomic layer
  etching of silicon nitride using hydrogen plasma and fluorine radicals*, **ACS Appl.
  Mater. Interfaces 11, 2019**, DOI [10.1021/acsami.9b11489](https://doi.org/10.1021/acsami.9b11489).
- X. Wang, E. A. Hudson, P. Kumar, S. Agarwal, *In situ monitoring surface reactions in
  cryogenic ALE of silicon nitride …*, **Chem. Mater. 36, 11042–11050 (2024)**, DOI
  [10.1021/acs.chemmater.4c01835](https://doi.org/10.1021/acs.chemmater.4c01835).
- Same group, *Effect of O₂ dilution and substrate temperature on etching of SiN in
  C₄F₆/Ar plasma*, **J. Vac. Sci. Technol. A 43(5), 2025** `[VERIFY]` article number.

### 3.7 Physical-sputter compilations (block B5/B7, L2 with a large measured corpus)

- Y. Yamamura, H. Tawara, *Energy dependence of ion-induced sputtering yields from
  monatomic solids at normal incidence*, **At. Data Nucl. Data Tables 62, 149–253
  (1996)**, DOI [10.1006/adnd.1996.0005](https://doi.org/10.1006/adnd.1996.0005).
  Experimental + ACAT points for many ion/target pairs with an empirical fit whose
  parameters are best-fit to the data. This is the closest thing to a *measured
  database* in the whole atlas, and Lam uses it for `E₀`.
- W. Eckstein, R. Preuss, *New fit formulae for the sputtering yield*, J. Nucl. Mater.
  **320**, 209 (2003) `[VERIFY]` — the fusion-community successor.
- Y. Yamamura, Y. Itikawa, N. Itoh, *Angular dependence of sputtering yields of
  monatomic solids*, IPPJ-AM-26 (1983) `[VERIFY]` — the angular counterpart.

### 3.8 Databases — the negative result

There is **no** open, curated database of etch surface rate constants.
- **LXCat** (Pitchford *et al.*, Plasma Process. Polym. **14**, 1600098, 2017; DOI
  [10.1002/ppap.201600098](https://doi.org/10.1002/ppap.201600098)) hosts electron/ion
  scattering cross sections, swarm parameters, ion-neutral potentials, oscillator
  strengths (~33 public databases as of 2022). Plasma-surface interaction is named in
  scope statements, but no surface rate-constant database is published.
- **NIST** offers gas-phase kinetics; no etch-surface set.
- **IAEA AMD / plasma-wall interaction** databases cover sputtering yields and
  reflection coefficients for fusion-relevant materials — usable for B5/B7 only.
- **INPTDAT / Plasma-MDS** (German open plasma data platform, <https://www.plasma-mds.org>)
  is a metadata schema + repository, not a constants corpus. `[VERIFY]` current content.
- Reactor-scale open data does exist (the repo already vendors Zenodo
  [10.5281/zenodo.17122442](https://doi.org/10.5281/zenodo.17122442), Bosch SF₆/C₄F₈
  wafer measurements) but carries no surface constants.

**Practical conclusion: the deck-building corpus is a set of ~15 open PhD theses and
~25 papers, not a database. Building a curated, provenance-levelled machine-readable
version of it is itself a defensible contribution.**

---

## 4. Buildability scores per candidate system

Scoring: for each block B1–B11, the best available provenance level. Score =
`(#blocks at L0/L1/L2) / 11`, with L3-only or absent blocks counted as gaps.

### 4.1 Cl₂ / poly-Si — **9.5 / 11 = 86 %** (BEST)

| Block | Best source | Level |
|---|---|---|
| B1 species/fluxes | Guo T-6 Table 3.2; Yin T-4 Table 3-1 (mass spec, 7 gas mixes) | **L0** |
| B2 chemisorption bare | Chang T-1 Tables 3.4 & 5.2 (`s` at 3 energies, 2 ion types); Guo Table 2.2 (`S_Cl`=0.75, `S_Cl2`=0.204) | **L1** |
| B3 activated chemisorption | implicit — `s ∝ √E` (Chang: high-E ions create adsorption sites) | **L1** |
| B4 passivation/deposition | *only* mask-side (Cl recombination on PR) — Chang's single adjustable knob | **L3 / gap** |
| B5 substrate sputter bare | Chang Table 3.1 (A=0.04, E_th=35 eV) + Guo Table 2.1 (Yamamura closed form) | **L1/L2** |
| B6 reacted-complex channel | Chang `β` at 3 energies, E_th = 16 eV (Ar⁺) / 10 eV (Cl⁺); Guo SiCl₂-formation rows | **L1** |
| B7 mask sputter | Chang ch. 4 (SiO₂ hard mask, Table 4.2) + ch. 4.8 (photoresist w/ Cl) | **L1** |
| B8 film/state transitions | Guo vacancy generation / annihilation / densification rows | **L1** |
| B9 products | Chang: SiCl₄ dominant; Guo: SiCl₂ dominant (a real, citable disagreement) | **L1**, contested |
| B10 redeposition | Chang SiCl₂-beam experiment: ~10× etch suppression | **L0** |
| B11 spontaneous | Walker Arrhenius (k₀ = 2.57e-14, Ea = 0.29 eV); negligible for undoped poly-Si | **L0** |

**Blind profile gates available (none of the constants saw them):**
1. **Mahorowala I** (JVST B 20, 1055): 3-variable fractional-factorial Cl₂ profiles,
   PR and oxide masks — the primary gate.
2. **Mahorowala III** (1077): mask faceting, sidewall deposition, microtrenching —
   secondary/mechanism gates.
3. **Levinson (LLNL)** JVST A 15, 1902 (1997) / JVST B 18, 172 (2000) — a *disjoint
   group's* profile set, the strongest possible independence.
4. Already in-repo: Hwang & Giapis notching (Cl₂ HDP), Nozawa 1995 — these become
   *chemistry-consistent* for the first time.

**Held-out non-profile gate:** Guo Table 2.3 (steady-state surface Si:Cl fraction vs
E and N/I, cross-checked against 5 independent references) → a direct gate on petch's
mixed-layer composition state, which no current gate touches.

**Only real gap:** B4. Cl₂/poly-Si has essentially no polymer/passivant, which is why
this system is *easy* — but the Cl-recombination-on-photoresist constant is genuinely
unmeasured and must be declared L3 in the deck.

### 4.2 C₄F₈ (or C₂F₆) / Ar → SiO₂ — **8 / 11 = 73 %**

| Block | Best source | Level |
|---|---|---|
| B1 | Yin T-4 (mass spec: ion + neutral spectra for C₄F₈/Ar, C₄F₈/O₂/Ar, C₂F₆/Ar); Guo Table 4.1 species list | **L0** |
| B2 | Guo Table 4.1 rows 1–4 sticking symbols; Kwon S-G/S-I | **L1** (values need re-read) |
| B3 | Guo vacancy/ion-incorporation coupling | **L1** |
| B4 deposition | Cunge–Booth CF₂ loss (**L0**, reactor-averaged only); Guo `S_N on C` / `S_N on O` (**L1**); Krüger's per-state map (**L3**) | **L0/L1 partial** |
| B5 | Guo Table 4.1 rows 5–8 (Si 0.042, O 0.018, C 0.009, F 0.023, calc.) | **L2** |
| B6 | Guo row 9 (Si+2F→SiF₂, A = 6.75) | **L1** |
| B7 mask | **gap** — Yin/Guo used blanket wafers; no a-C or PR mask data in this corpus | **gap** |
| B8 film state | Standaert 2004 (measured thickness + XPS composition + defluorination) | **L0** |
| B9 products | Guo rows 9–15: SiF₂, SiF₄, CO, CO₂, CF₂, O₂, SiO — a complete measured-fit branching set | **L1** |
| B10 | **gap** | gap |
| B11 spontaneous | Guo row 20 (Si+4F→SiF₄, A = 2.99e-5) | **L1** |

**Blind profile gates:** Krüger 2024 (already vendored, `data/experimental/krueger_2024/`)
becomes a *genuine blind gate* the moment the deck's constants are L1 instead of L3 —
this is the highest-value experiment in the whole atlas. Also Huard/Kushner JVST A 37,
031304 (2019) HAR Ar/C₄F₈/O₂ (open PDF), and Standaert 1997 blanket dep→etch transition.

**Gaps:** mask sputter (B7) and redeposition (B10). Krüger's a-C mask rows would have
to be carried over as declared-L3, or the gate restricted to the oxide arm.

### 4.3 SF₆–O₂ / Si — **5 / 11 = 45 %** (worst provenance, best existing petch coverage)

| Block | Best source | Level |
|---|---|---|
| B1 | Belen 2005 reactor fluxes; de Boer process conditions | L1-ish (reactor model, not beam) |
| B2 F sticking | **Belen: profile-fitted** | **L3** |
| B3 | not resolved | gap |
| B4 O passivation sticking | **Belen: profile-fitted**; Vallier-style XPS exists only for HBr/Cl₂ | **L3** |
| B5 | Yamamura / Vitale S-E (Si in F₂ plasma, `Y ∝ √E`, angular) | **L1/L2** |
| B6 ion-enhanced yield | **Belen: profile-fitted**; but Gray S-A measured F+Ar⁺ on poly-Si over decades of flux, and Vitale S-E measured F₂-plasma yields vs E and θ | **L3 now, L1 available** |
| B7 mask | gap | gap |
| B8 | SiOₓF_y layer — no measured constants at feature scale | gap |
| B9 products | Flamm 1981 (SiF₂ formation rate-limiting; SiF₄ final) | **L0** |
| B10 | gap | gap |
| B11 spontaneous | **Flamm 1981**: `R_F = 2.91e-12·T^½·n_F·e^{−0.108 eV/kT}` Å/min | **L0** |

**Blind profile gates:** de Boer 2002 (already vendored) — but *only* blind if the
constants are re-sourced from Gray/Vitale rather than Belen. Also Gomez/Belen/Aydil
JVST A 22, 606 (2004) HAR holes; Belen JVST A 24, 350 (2006) SF₆/O₂/HBr.

**Verdict:** the SF₆/O₂ deck is buildable *today* only by importing three L3 constants.
An **L1 re-derivation** is possible (Gray S-A gives the F/Ar⁺ ion-enhancement over
decades of flux ratio; Vitale S-E gives `Y(E,θ)` in a real F₂ plasma) but requires
digitizing two figure sets and re-regressing — a genuine work item, not a citation.

### 4.4 HBr / Si — **4 / 11 = 36 %** (weakest)

| Block | Best source | Level |
|---|---|---|
| B1 | Cunge JVST B 20, 2137 (2002) measured ion flux composition | **L0** |
| B2 | Jin T-2 / S-F derived HBr neutral sticking — **"lower than Cl, probably due to larger Br size and lower reactivity"**; numeric values only in the thesis, whose OCR is unusable | **L1 but not yet extractable** |
| B3 | gap | gap |
| B4 SiOₓBr_y passivation | Vallier JVST B 21, 904 (2003) XPS composition (**L0 composition**, no rate constants) | partial |
| B5 | Vitale S-E (`Y ∝ √E`; angular falls immediately off-normal) | **L1** |
| B6 | Vitale/Jin: `E_th = 10 eV`, HBr yield ≈ Cl₂ yield | **L1** |
| B7 mask | gap | gap |
| B8 | gap | gap |
| B9 products | Cunge JAP 96, 4578 (2004) SiClₓ in HBr/Cl₂/O₂; spinning-wall SiBr₂/SiBr₄ | **L0 qualitative** |
| B10 | gap | gap |
| B11 spontaneous | negligible (no measurable HBr spontaneous etch of Si at 300 K) `[VERIFY]` | — |

**Blind profile gates:** Mahorowala I (the HBr half of the same fractional-factorial
set) is excellent and already needed for Cl₂ — so the marginal cost of adding HBr
*after* Cl₂ is small. But the constants are the blocker.

### 4.5 Summary scoreboard

| System | Score | L3 imports required | Blind profile gate quality | Marginal cost after previous deck |
|---|---|---|---|---|
| **Cl₂ / poly-Si** | **86 %** | 1 (Cl recombination on PR) | **Excellent** (Mahorowala I+III, Levinson ×2, in-repo notch data) | — |
| **FC / SiO₂ (Kwon–Yin–Guo)** | **73 %** | 0 for the oxide arm; mask arm needs Krüger's a-C rows | **Excellent** (Krüger 2024 becomes truly blind; Huard HAR) | medium — reuses the mixed-layer engine as-is |
| **SF₆–O₂ / Si (Belen)** | **45 %** | 3 (F sticking, chemical rate constant, ion-enhanced yield) | de Boer, but **not blind** on those 3 | very low (already implemented) |
| **HBr / Si** | **36 %** | ≥4 | Excellent gate, unavailable constants | low *after* Cl₂ |

---

## 5. Recommendation

### 5.1 Keep Belen SF₆/O₂ as "Deck 2", but retitle its purpose

Moving the existing Belen arm into deck format is worth doing — it proves *engine
generality* (two chemistries, one engine, zero code branches), it is nearly free, and
it is already validated. **But it must not be sold as a beam-measured deck.** Concretely:

```
"declared_fitted": {
  "constants": ["fluorine_sticking", "chemical_etch_rate_constant",
                "ion_enhanced_etch_yield"],
  "provenance": "Belen et al., JVST A 23, 99 (2005), DOI 10.1116/1.1830495 — "
                "abstract: these are 'experimentally inaccessible parameters ... "
                "determined by matching simulated feature profiles with experiments'",
  "provenance_level": "L3_profile_fitted"}
```

and the de Boer claim in `RECONCILIATION.md` / the docs site should be reworded from
"validated" to "**width/time transfer gate on a profile-fitted constitutive law**".
This is a doctrine correction, not a physics regression.

### 5.2 Deck 3 (first *new-science* deck): **Cl₂ / poly-Si, Chang–Sawin**

Rationale: highest buildability (86 %), one L3 constant, and the only system where the
constants and the profile experiments come from *different* papers by construction. It
also makes the existing notching/charging gates chemistry-consistent for the first time
(they currently run SF₆/O₂ surface chemistry under a Cl₂-derived charging table).

**Source list, in order of use:**

*Constants (all L0–L2 except one):*
1. J. P. Chang, PhD MIT 1998, [1721.1/50356](http://hdl.handle.net/1721.1/50356) —
   **Tables 3.1, 3.2, 3.3, 3.4, 4.1, 4.2, 5.1, 5.2** (transcribed in §2.3 above).
2. W. Guo, PhD MIT 2009, [1721.1/46600](http://hdl.handle.net/1721.1/46600) —
   **Tables 2.1, 2.2, 2.3, 3.1, 3.2** (transcribed in §2.4).
3. Chang & Sawin, JVST A **15**, 610 (1997), DOI 10.1116/1.580692 — `E_th = 10 eV`,
   Langmuir N/I saturation.
4. Chang, Arnold, Zau, Shin, Sawin, JVST A **15**, 1853 (1997), DOI 10.1116/1.580652.
5. Yin & Sawin, JVST A **26**, 161 (2008), DOI 10.1116/1.2821750 — the **angular**
   `Y(θ; E, N/I)` family; Yin thesis Table 3-1 for the flux compositions.
6. Guo, Bai & Sawin, JVST A **27**, 388 (2009), DOI 10.1116/1.3085722 — the published
   version of the mixing-layer deck.
7. Flamm/Walker Arrhenius constants via Chang Table 3.2 (block B11).
8. Yamamura & Tawara, ADNDT **62**, 149 (1996) — B5 cross-check + `E₀` for Si.

*Blind profile gates (declare sealed before any run):*
- Mahorowala, Sawin, Jones, Labun, JVST B **20**, 1055 (2002) — fractional-factorial
  Cl₂ profiles (primary).
- Mahorowala & Sawin, JVST B **20**, 1077 (2002) — faceting/microtrenching (mechanism).
- Levinson, Shaqfeh, Balooch, Hamza, JVST A **15**, 1902 (1997) and JVST B **18**, 172
  (2000) — independent-group profiles (the strongest gate).
- Guo Table 2.3 — held-out **surface-composition** gate (novel; no other petch gate
  tests the mixed-layer state directly).
- In-repo: `hwang_giapis_1997`, `nozawa_1995` re-run under the correct chemistry.

*Effort estimate:* the deck is a data file; the engine already carries the mixed-layer
form, threshold-power sputter, Langmuir-type coverage, and an ion-reflection cascade.
The real work is (i) digitizing Yin's `Y(θ)` families, (ii) resolving the SiCl₄
(Chang) vs SiCl₂ (Guo) product disagreement in the F/Cl ledger validator, (iii)
choosing whether `s(E) ∝ √E` is carried as an energy-dependent chemisorption block
(the schema currently assumes energy-independent chemisorption maps — **this is a
schema extension the Cl₂ deck forces**).

### 5.3 Deck 4: **C₄F₈(/C₂F₆)–Ar → SiO₂, Kwon–Yin–Guo — as the Krüger cross-check**

Rationale: this is the *retirement mechanism* for Krüger's five fitted constants that
the current doctrine ("beam measurement > derivation > declared-fitted") implicitly
assumes must be a new measurement. It is not: an L1 deck for the same material system
already exists. Running **two decks (Krüger-L3 and Kwon/Yin/Guo-L1) through one engine
against the same Krüger feature** is a stronger result than either alone — it measures
how much of the feature agreement was carried by the profile-fitted constants.

**Source list:**
1. Kwon & Sawin, JVST A **24**, 1906 (2006) I, DOI 10.1116/1.2336225 — measured
   neutral/ion composition (QMS) and yield (QCM) for C₄F₈ and C₄F₈+80%Ar.
2. Kwon & Sawin, JVST A **24**, 1914 (2006) II, DOI 10.1116/1.2336226 — the TML
   formalism (compare against `petch.mixed_layer` line by line).
3. Kwon, Bai & Sawin, JVST A **24**, 1920 (2006) III, DOI 10.1116/1.2336227.
4. Yin & Sawin, JVST A **26**, 161 (2008) — SiO₂ + coral `Y(θ)`; JVST A **25**, 802
   (2007) — kinetics/roughening; JVST A **26**, 151 (2008) — Ar-only sputter arm.
5. Guo thesis Table 4.1 / 4.2 (20-reaction deck + angular forms) — §2.4 above; and Yin
   thesis ch. 3 for the C₄F₈-percentage and O₂-addition composition series.
6. Standaert *et al.*, JVST A **22**, 53 (2004), DOI 10.1116/1.1626642 — film thickness,
   XPS composition, defluorination (block B8, L0).
7. Cunge & Booth, JAP **85**, 3952 (1999) — CF₂ surface loss (block B4, L0).
8. Kanarik *et al.*, JVST A **35**, 05C302 (2017) (open) — SiO₂ ALE window/synergy as an
   independent energy-threshold anchor.

*Blind gates:* Krüger 2024 base case + transfer observations (in-repo); Huard/Kushner
JVST A 37, 031304 (2019) HAR; Standaert JVST A 15, 1881 (1997) dep→etch transition bias.

### 5.4 Deferred: HBr

Do it **only after Cl₂**, and only if someone obtains a readable copy of Jin's thesis
(MIT 1721.1/28357 — the DSpace OCR is unusable; an ILL/print scan or the S-F paper's
tables are required). The profile gate (Mahorowala I, HBr arm) will already be in hand
from the Cl₂ deck, so the marginal cost is purely constant-extraction.

### 5.5 Schema changes this study forces

1. **`provenance_level` per constant** — `L0_measured | L1_beam_regressed |
   L2_semiempirical | L3_profile_fitted`. Validation should *require* it and should be
   able to report a deck's L3 fraction. A profile gate on a deck with any L3 constant
   in the etched material must be labelled non-blind.
2. **Energy-dependent chemisorption.** Chang measures `s(E) ∝ √E` (0.02 → 0.30 over
   35 → 100 eV). The current `chemisorption.bare` is a scalar map. Either add an
   energy-law form or the Cl₂ deck cannot be expressed.
3. **Per-ion-species rows.** Guo's decks are indexed by *which* ion (Ar⁺ / Cl⁺ / Cl₂⁺)
   with distinct A and E_th. The current schema aggregates ions.
4. **Product-redeposition block (B10).** Chang's SiCl₂-beam result (≈10× suppression)
   has no home in the schema.
5. **Spontaneous-etch block (B11).** Arrhenius `k₀`, `Ea` are L0 for both F/Si and
   Cl/poly-Si and belong in the deck, not the kernel.

---

## 6. Corrections to existing repo documents

| File | Claim | Correction |
|---|---|---|
| `CHEMISTRY_EXPANSION_PLAN.md` §W2 gate 2 + Sources | "Chang, Arnold & Sawin, JVST B 18, 172 (2000) — Ion-assisted etching and profile development of silicon in molecular chlorine" | Authors are **Levinson, Shaqfeh, Balooch, Hamza** (LLNL). Two papers: JVST **A 15**, 1902 (1997), DOI 10.1116/1.580658; JVST **B 18**, 172 (2000), DOI 10.1116/1.591170 (title adds "and atomic") |
| `RESEARCH_THESIS_MINING_2026-07-23.md` §S2 | "Guo, Wei — PhD MIT 2008 — MIT DSpace 1721.1/43201" | 43201 = **Kawai, Hiroyo (2008)**. Guo, Wei (2009) = **1721.1/46600** |
| `RESEARCH_ENERGY_DEPOSITION_ANCHORS_2026-07-23.md` | "Chang & Sawin (beam studies, JVST A ~15, 610 (1997)) `[VERIFY citation]`" | **Verified**: J. Vac. Sci. Technol. A **15**, 610–615 (1997), DOI 10.1116/1.580692 |
| `RESEARCH_MECHANISM_COMPLETENESS_2026-07-25.md` + `chemistry_deck.py` | SiO₂-complex sputter `E_th = 35 eV`, "Chang–Sawin" | The only 35 eV in Chang's corpus is the **Ar⁺-on-Si physical-sputter** threshold (Table 3.1). Chang's SiO₂ measurements (Table 4.2) are at 70/100 eV and quote no complex threshold. Trace to Krüger Appendix B's own reference before calling it Chang–Sawin `[VERIFY]` |
| `PROGRAM_ROADMAP_2026-07-24.md` §1 | "Sawin/Yin 2008 sealed blind campaign (declared #1)" | Yin 2008 is a **blanket beam-yield** dataset — it is a *constitutive* gate, not a profile gate. The profile gates for the same chemistry are Mahorowala JVST B 20 (2002) and Levinson. Both should be named in the sealed declaration |
| docs / `RECONCILIATION.md` | de Boer described as validating the SF₆/O₂ chemistry | Belen's three key constants are **profile-fitted by the source's own abstract**. de Boer is a *transfer* gate, not a blind first-principles prediction |

---

## 7. Local artifacts from this session

Full PDFs + extracted text (scratchpad, not vendored — MIT DSpace open-access, but not
redistributed into the repo):

```
<scratch>/chang1998.pdf   chang1998.txt   # T-1, 1721.1/50356, good OCR, all tables
<scratch>/yin2007.pdf     yin2007.txt     # T-4, 1721.1/38973, clean text
<scratch>/guo2009.pdf     guo2009.txt     # T-6, 1721.1/46600, clean text, Tables 2.1/2.2/2.3/3.2/4.1/4.2
<scratch>/kanarik2017.pdf kanarik2017.txt # open OSTI copy, purl/1376399
<scratch>/jin2003.pdf     jin2003.txt     # T-2, 1721.1/28357 — OCR UNUSABLE
<scratch>/kwon2004.pdf    kwon2004.txt    # T-3, 1721.1/28353 — OCR UNUSABLE
```

Scratchpad root:
`/private/tmp/claude-501/-Users-stanislavdelaurentiis-chip-etch/352207a8-b462-4fe3-a650-c9b37f19523b/scratchpad`

DSpace bitstreams are reachable only via the REST API
(`https://dspace.mit.edu/server/api/core/bitstreams/<uuid>/content`); the `/handle/…`
and `/bitstreams/…/download` routes return 405/empty.

## 8. Open items marked `[VERIFY]`

1. Belen 2005 abstract wording — verify verbatim against the PDF (AIP blocks fetch;
   two independent search renderings agree).
2. Guo Table 4.1 rows 1–4 sticking coefficient *values* (`S_I`, `S_F`, `S_N on C`,
   `S_N on O`) — lost in OCR of the parameter column; re-read thesis p. 100 or S-I.
3. The "35 eV, Chang–Sawin" SiO₂-complex threshold provenance (see §6).
4. Spinning-wall HBr product study — exact citation.
5. Eckstein & Preuss 2003 and Yamamura IPPJ-AM-26 (1983) — page/report numbers.
6. INPTDAT/Plasma-MDS current surface-data content.
7. Whether HBr spontaneously etches Si at 300 K at a measurable rate (assumed no).
8. Agarwal group 2025 JVST A article number for the C₄F₆/Ar SiN paper.
