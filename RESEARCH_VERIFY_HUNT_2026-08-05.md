# Exhaustive `[VERIFY]` source hunt (2026-08-05)

Scope: run down the open `[VERIFY]` items across the research corpus to **tabulated /
verbatim** values. Rule applied throughout: fetch and read, never infer. Every claim below
is either (a) a verbatim quote with a locator, (b) a number digitised from a figure I
rendered myself (method stated), or (c) explicitly marked NOT FOUND with the routes tried.

**Do not commit.** Working artefacts live in
`…/352207a8-b462-4fe3-a650-c9b37f19523b/scratchpad/hunt/`.

---

## 0. Headline result: the primary source for BOTH headline items is the same document

The single highest-yield acquisition this session was:

> **J. P. Chang**, *Study of plasma-surface kinetics and simulation of feature profile
> evolution in chlorine etching of patterned polysilicon*, PhD thesis, MIT Chemical
> Engineering, 1998 (advisor H. H. Sawin). MIT DSpace handle **1721.1/50356**.
> Direct PDF (open, no auth): `https://dspace.mit.edu/bitstreams/5285146e-c906-4f6b-ab72-1b12764d8011/download`
> — 21.8 MB, 1.6 dpi-clean scan, `pdftotext` OCR is good.
> Local: `scratchpad/hunt/chang_thesis.pdf` + `.txt`.

`RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` §T-1 already flagged this thesis as
"Good OCR — tables recoverable" but it had **not been fetched**. It contains, in one
document: the Chang & Sawin 1997 angular curve (Ch. 5), the SiO₂/fluorocarbon angular
curve (Ch. 4), the measured fluorocarbon-**polymer** angular sputter curve (Fig. 4.16),
the beam geometry/flux-normalisation convention (Ch. 2), and the provenance of the
35 eV threshold (Ch. 4).

---

## 1. Crosslink bond multiplicity — **FOUND** (rule + two per-material values; no table exists)

### 1.1 Verdict

There is **no table** of max-bonds-per-material anywhere in the MCFPM lineage. There is a
**stated rule plus two worked per-material values**, and it appears in three independent
renderings of the same text. The value the repo has been carrying as "3, from the example
depicted" is the *figure's example*, **not** a per-material constant — the paper says so
explicitly and then gives the actual per-material numbers, which the thesis omits.

### 1.2 Verbatim — Krüger *et al.* 2024, journal version

F. Krüger, S. Lee, H. Nam(?), M. J. Kushner et al., *Autonomous hybrid optimization of a
SiO₂ plasma etching mechanism*, **J. Vac. Sci. Technol. A 42(4), 043008 (2024)**,
DOI **10.1116/6.0003554**. §IV, p. 043008-6/-7 (`tmp/pdfs/krueger-2024.txt` L389-402):

> "Crosslinking occurs during the deposition of eligible materials [Figs. 5(a) and 5(b)].
> (By material, we refer to incident species, such as CF₃ or C₂F₃, or their counterparts
> in the film.) **Each material has a maximum number of neighbors to which it can
> crosslink based on the number of available bonds (three in the example in Fig. 5).
> For example, CF₂ would have a maximum of two crosslinks and CF₃ would have a maximum
> of a single crosslink.** When a crosslink eligible deposition material is incident onto
> the surface, a randomly ordered search is conducted of nearest neighbors for other
> crosslink eligible materials that have available bonding sites. Based on a prescribed
> crosslinking probability and choice of a random number (0,1) crosslinks are made (or
> not made) with the eligible neighbors."

**Independently confirmed in the DOE accepted manuscript** (a different typesetting, so
this is not an OCR artefact of the two-column layout):
OSTI record **2572026**, full text at `https://www.osti.gov/servlets/purl/2572026`
(`scratchpad/hunt/kr2024_osti.txt` L551-554):

> "…based on the number of available bonds (**3 in the example in Fig. 5**). For example,
> **CF₂ would have a maximum of 2 crosslinks and CF₃ would have a maximum of a single
> crosslink.**"

### 1.3 Verbatim — Krüger thesis §2.2.3 (the source of the current `[VERIFY]`)

`tmp/pdfs/krueger_thesis.txt` L2474-2478, p. 79:

> "Crosslinking occurs during the deposition of eligible materials (Figure 2.2a and b).
> Each material has a maximum number of crosslink partners associated with it, which is
> based on the number of available bonds (**3 in the example depicted in Figure 2.2**).
> During deposition bonds to random eligible cell neighbors can be formed, increasing the
> respective crosslink number, which is tacked [sic] for every cell (Figure 2.2c)."

The thesis genuinely stops at "3 in the example". **The paper does not** — it gives the
rule *and* two instances. This closes the `[VERIFY]` in
`RESULTS_ML19_ENDPOINT_2026-08-05.md` L103 and `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md`
L228/L248 (the "m = 3 → ≈1.4×, m = 8 → ≈1.0×" sensitivity): **m = 8 is excluded by the
source.** The admissible per-material values in this lineage are small integers ≤ 3.

### 1.4 The rule, and what is and is not derivable

Stated rule: `m(material) = number of available (unsatisfied) bonds`.
Stated instances: **CF₂ → 2**, **CF₃ → 1**, and "3" for the generic species drawn in Fig. 5.

Consistent with tetravalent carbon, `m(CF_y) = 4 − y`, which reproduces both stated
instances and makes the Fig.-5 "3" equal to **CF**. That arithmetic is **mine, not
Krüger's** — he never writes `4 − y` and never states `m(CF)` or `m(C₂F₃)`. Treat
`m(CF) = 3` as a strong inference, `m(C₂F₃)` as unknown.

**Corroborating structural evidence from Appendix B (Table B.0.1)**: the crosslinked
species that actually appear in the converged mechanism are exactly
`CF(xs)`, `CF2(xs)`, `CF3(xs)`, `AC(xs)` — **`C2F3(xs)` never appears** in 249 `(xs)`
occurrences (`tmp/pdfs/krueger_thesis.txt`, grep `→ *[A-Za-z0-9]*(xs)`). So in the
delivered mechanism C₂F₃ is deposited (rows `…(s) + C2F3 → …C2F3(s)`, p₀ = 0.03) but is
apparently **not crosslink-eligible**. This is an argument from absence over a table I
read in full; flag it as such if it becomes load-bearing.

### 1.5 What does NOT exist (routes tried, all negative)

| Target | Routes tried | Outcome |
|---|---|---|
| **MCFPM User's Manual** | `cpseg.eecs.umich.edu` root + `projects.html`, `data.html`, `Help.html`, `classes.html`, `short_courses.html`, `preprints.html`, `theses.html` (all fetched and link-extracted); `Projects/MCFPM/MCFPM.htm` fetched and de-tagged in full; two targeted web searches | **NOT FOUND — no such public document.** The MCFPM page is ~200 words of prose. The only downloadables anywhere on the site are `pub/data/e_reactions.pdf`, two short-course decks (`MCSHORT_0502.pdf`, `electronegative_master_class_kushner.pdf`), and 56 theses. No manual, no input deck, no chemistry-file spec. |
| **Public MCFPM input decks** (nanoHUB / GitHub / course material) | web search; `classes.html` → only `eecs517_2012` | **NOT FOUND** |
| **Krüger 2024 supplementary material** ("See supplementary material online for listing of the full reaction mechanism", ref. 86) | `doi.org/10.1116/6.0003554` (403), `pubs.aip.org/avs/jva/article/42/4/043008/3299279` (403), AIP `article-supplement` pattern probe (403), OSTI full text (39 pp — **no supplement appended**) | **NOT FOUND via open routes.** However it is redundant: Krüger thesis **Appendix B / Table B.0.1** is the same mechanism listing and is already in-repo (`RESEARCH_APPENDIX_B_VERBATIM_2026-07-27.md`). Appendix B has **no** crosslink-multiplicity column — the mechanism table simply does not carry this parameter. |
| **Crosslinking probability value** | grep of both Krüger documents for "crosslink(ing) probability" | **NOT PUBLISHED.** Only "Based on a prescribed crosslinking probability…". The number is not in the paper, the thesis, or Appendix B. |
| **Krüger's other papers** | OpenAlex author sweep (25 works); JVST A 40 (2022) DOI 10.1116/6.0002290 downloaded OA from `osti.gov/servlets/purl/2331488` and grepped — **zero** occurrences of "crosslink"; Phys. Plasmas 31, 033508 (2024) DOI 10.1063/5.0189397 (AIP 403; VWT-scope, not mechanism-scope); arXiv API: no hits | Crosslinking module is **new in the 2024 JVST paper + thesis**; no earlier or later restatement. |
| **Later Michigan theses** | `theses.html` full list — newest are Konina, Polito, Krüger, all 2024. Polito thesis downloaded and grepped: 7 "crosslink" hits, **all** polystyrene/APPJ surface chemistry, unrelated to the MCFPM polymer module | **NOT FOUND** |

### 1.6 Adjacent, and useful: the *other* crosslink constant in the lineage

Wang (Mingmei) thesis, `tmp/pdfs/wang_mingmei_phd_thesis.txt` L1247 and L3014 — the same
MCFPM, PR/PMMA rather than fluorocarbon:

> "we included bond-breaking in the PR and conversion to crosslinked sites. **The
> sputtering yield of the cross-linked PR sites is five times smaller** than [that of
> non-crosslinked PR]."

> "…the sputtering yield for cross-linked PMMA is **5 times smaller** than for
> non-crosslinked PMMA."

That is the *property change on crosslinking* (Fig. 5(d) "change in properties"),
tabulated as a factor of 5 for PR. Different material from FC polymer — a scale anchor,
not an import.

---

## 2. Chang & Sawin 1997 angular dependence — **FOUND** (verbatim text + digitised curve)

### 2.1 Citation, confirmed

J. P. Chang and H. H. Sawin, *Kinetic study of low energy ion-enhanced polysilicon etching
using Cl, Cl₂, and Cl⁺ beam scattering*, **J. Vac. Sci. Technol. A 15(3), 610–615 (1997)**,
DOI **10.1116/1.580692**. (Crossref/OpenAlex biblio confirms vol 15, issue 3, pp. 610–615.)

### 2.2 Verbatim: the published abstract carries the two off-normal numbers

Retrieved via OpenAlex (`api.openalex.org/works/doi:10.1116/1.580692`); AIP direct fetch
returns 403. Publisher abstract, final two sentences, verbatim:

> "The angular dependence of ion-enhanced etching yield was also measured. **The etching
> yield was reduced by approximately 30% and 50% when ion impingement angles of 60° and
> 70° off-normal were used, respectively.**"

Also verbatim from the same abstract, and relevant to the deck: *"The ion energy dependence
was a linear function of (E_ion^1/2 − E_th^1/2), where the threshold energy E_th was found
to be approximately 10 eV."*

### 2.3 Verbatim: the thesis prose, which adds the plateau edge and the endpoint

Chang thesis §5.3, **p. 115** (rendered page image `scratchpad/hunt/chang_p115-115.png`;
text `chang_thesis.txt` L5691-5706):

> "The measured etching yields at various ion incident angles at a flux ratio of 200 are
> shown in Figure 5.9. **The etching yield showed no significant change as the ion incident
> angle increased from normal to 40° off normal, but decreased by 30 % and 50 % at 60° and
> 70° off-normal angles. The zero etching yield at glancing angle (90° off-normal) is
> assumed. Maximum etching yield is observed at normal ion incident angle and the etching
> yield starts decreasing at 30-40° off-normal angles**, similar to the angular dependence
> we measured for Cl/Ar⁺ etching of polysilicon. A highly chlorinated layer is produced by
> chlorine ion bombardment and the maximum etching yield is achieved at normal incident
> angle. **The etching yield at off-normal incident angles is expressed in the following
> form: Y(φ) = c(φ)·Y(φ = 0°)** where Y(φ = 0°) is the etching yield measured at normal ion
> incident angle and **c(φ) is a constant extracted from a polynomial fit to the
> experimental data**, as indicated by the line in Figure 5.9."

Fig. 5.9 caption, verbatim:

> "Figure 5.9: Angular dependence of ion enhanced polysilicon etching with Cl and 35 eV Cl⁺
> at a flux ratio of 200. Maximum etching yield is observed at normal ion incident angle and
> the line represents a polynomial fit to the experimental data."

The polynomial **coefficients are not printed** — not in §5.3, not in Ch. 7 (the DSMC
simulator chapter says only "The reaction probability of ions was calculated based on the
measured ion angular dependence", `chang_thesis.txt` ~L6950). The Ch. 4 companion prints
the *form* only: `Y(φ) = Σ_{i=1..6} a_i cos^i(φ) = c(φ)·Y(φ=0°)` (Eq. 4.1, p. 93) with the
`a_i` unlisted. **`a_i` / `c(φ)` closed-form: NOT FOUND** (thesis is the fullest available
rendering; the JVST paper itself is AIP-403 and its abstract quotes only the two
percentages).

### 2.4 The curve, digitised

Method: `pdftoppm -r 400 -gray` of thesis p. 115; plot box located from full-width
runs (`x`: 1106 px = 0°, 2297 px = 90°; `y`: 3403 px = 0, 2280 px = 1.5); **markers**
recovered by 7×7 binary erosion + connected components (kills the thin fit line and axes,
keeps the filled squares); **fit line** recovered by per-column first thin dark run.
Script inline in session; raster `scratchpad/hunt/chang_hi-115.png`.

Measured data points (filled squares), absolute yield in Si/Cl⁺ and normalised:

| θ (deg) | Y (Si/Cl⁺) | Y/Y(0) | thesis prose |
|---|---|---|---|
| 0 | 1.099 | 1.000 | max at normal |
| 20 | 1.101 | 1.001 | "no significant change …to 40°" |
| 40 | 1.100 | 1.001 | " |
| 60 | 0.799 | **0.727** | "decreased by 30 %" |
| 70 | 0.602 | **0.548** | "decreased by 50 %" |
| 90 | 0.007 | 0.006 | "zero … is assumed" |

Digitisation reproduces the prose to within the prose's own rounding (0.727 vs "30 %"
→ 0.70; 0.548 vs "50 %" → 0.50). The *polynomial fit line*, traced at 5° intervals and
normalised to Y(0) = 1.099, against petch's implemented class-2 form:

```
class 2 as implemented (src/petch/mixed_layer.py:227-233)
    f2(theta) = min(1, cos(theta) / cos(45 deg))
```

| θ | Chang fit (digitised) | Chang markers | petch `f2` | petch − Chang |
|---|---|---|---|---|
| 0 | 1.000 | 1.000 | 1.000 | — |
| 10 | 0.980 | — | 1.000 | +0.020 |
| 20 | 1.002 | 1.001 | 1.000 | −0.002 |
| 30 | 1.013 | — | 1.000 | −0.013 |
| 40 | 1.002 | 1.001 | 1.000 | −0.002 |
| 45 | 0.947 | — | 1.000 | **+0.053** |
| 50 | 0.894 | — | 0.909 | +0.015 |
| 55 | 0.825 | — | 0.811 | −0.014 |
| 60 | 0.729 | 0.727 | 0.707 | −0.022 |
| 65 | 0.638 | — | 0.598 | −0.040 |
| 70 | 0.549 | 0.548 | 0.484 | **−0.065** |
| 75 | 0.405 | — | 0.366 | −0.039 |
| 80 | 0.275 | — | 0.246 | −0.029 |
| 85 | 0.142 | — | 0.123 | −0.019 |
| 90 | 0.015 | 0.006 | 0.000 | −0.015 |

**Verdict on the `RESULTS_ANGULAR_CLASSES_2026-08-05.md` `[VERIFY]`:** the chosen
interpolation `min(1, cosθ/cos45°)` is **within 0.065 of the measured curve everywhere**,
max error at 70° where petch is 12 % low, and it reproduces both stated endpoints and the
plateau edge. It slightly over-extends the plateau (Chang's fit has already fallen to 0.947
at 45°, i.e. the real roll-off begins at "30-40°" as the prose says, not at 45°). The
`[VERIFY]` can be **downgraded to a quantified approximation** with the table above as its
receipt. Nothing in the fit motivates changing the form.

### 2.5 CONTEXT WARNING — do not import this curve blindly

Three separate transplant gaps sit between this measurement and petch's use of it:

1. **Chemistry and material.** Chang & Sawin 1997 is **Cl⁺ + Cl on polysilicon at 35 eV,
   neutral/ion flux ratio 200**. Krüger's Appendix B assigns `∠=2` to
   `SiO2CF(s) + Ar+ → SiF + CO2` — **SiO₂·CF_x complex sputter in a fluorocarbon chemistry
   at ~1.5 keV**. Different reactant, different substrate, ~40× the ion energy. Chang
   himself warns against exactly this transplant (thesis p. 93, verbatim): *"The comparison
   shown in Figure 4.7 suggests that **the chemistries used can significantly change the
   angular dependence in etching a certain material**, depending on whether their chemical
   reactivity exceeds the effect of energy and momentum transfer induced by sputtering."*
2. **Flux-ratio conditioning.** The curve is measured at a **fixed N/I = 200**, i.e. on a
   saturated, "highly chlorinated" surface. Chang attributes the normal-incidence maximum
   to that saturated layer. Inside a high-AR feature the local N/I ratio collapses with
   depth, so the shape is not obviously transferable to the deep floor.
3. **Normalisation convention — checked, and it favours petch.** Chang thesis §2, p. 39
   (`chang_thesis.txt` L1897-1901), verbatim: *"When samples are set at off-normal
   positions, the incident laser is sometimes introduced from the top view port… In this
   case, **the measured ionic and neutral flux has to be corrected taking into account the
   angle between the source and the sample surface normal to yield the correct incident
   fluxes.**"* So Y(θ) is **per ion arriving per unit surface area** — the cosθ projection
   is already divided out. That is exactly petch's convention (`atom_flux` carries cosθ,
   `f(θ)` multiplies on top), so the two are commensurate. This was previously assumed;
   it is now sourced.

### 2.6 The in-chemistry corroboration Chang printed on the facing page

**Fig. 4.7, p. 93** (image `scratchpad/hunt/chang_p93-093.png`) plots two SiO₂ angular
curves on one axis. Caption verbatim:

> "Figure 4.7: Etching yield of silicon dioxide as a function of ion incident angle. Solid
> diamonds represent the etching yields measured with Cl and 100 eV Ar⁺ at a flux ratio of
> ~90. **Maximum etching yield is observed at 60° off normal ion incident angle.** The
> hollow triangles represent the etching yields of silicon dioxide measured with **200 eV
> CF_x⁺** and **maximum etching yield is observed at normal incident angle**."

Body text, p. 92-93, verbatim:

> "The maximum etching yield is observed at 60° off-normal ion incident angle, where the
> momentum transfer to the silicon dioxide surface is the most sufficient. **This angular
> dependence indicates the ion enhanced etching of silicon dioxide is more of a physical
> sputtering process than a chemically driven process.**"

> "**However, when SiO₂ is etched by fluorocarbon gases [Mayer, 1981], the angular
> dependence changed drastically: the maximum etching yield is observed at normal ion
> incidence angle**, as shown in Figure 4.7. This result indicates that SiO₂ etching in the
> fluorocarbon chemistry was chemically driven, as CO, CO₂, COF₂, CF₄ and SiF₄ was liberated
> as etching products."

Source of the CF_x⁺ curve: **T. M. Mayer, R. A. Barker, L. J. Whitman**, *Investigation of
plasma etching mechanisms using beams of reactive gas ions*, **J. Vac. Sci. Technol. 18(2),
349 (1981)** (Chang bibliography, `chang_thesis.txt` L7796).

**This is the load-bearing point for petch.** Krüger's `∠=2` assignment to the
SiO₂·CF_x complex row is cited to the *wrong-chemistry* paper (Cl/poly-Si), but the
**right-chemistry** measurement exists and says the same thing: **SiO₂ + CF_x⁺ peaks at
normal incidence and rolls off monotonically.** So the class-2 assignment survives the
context check on physics grounds even though its citation is a chemistry transplant.

---

## 3. NEW — a directly relevant refutation of a standing repo claim (Kress `B = 9.3`)

`RESEARCH_LIP_CERTAINTY_2026-08-04.md` §2.4 states: *"(b) There is **no measured FC-film
angular sputter yield** in the literature to gate the lip law with."*

**That is false.** Chang thesis **Fig. 4.16, p. 104** (image
`scratchpad/hunt/fig416b-104.png`) reprints one, and names it.

Source, from Chang's bibliography (`chang_thesis.txt` L7656):

> **A. M. Barklund and H.-O. Blom**, *Influence of polymer formation on the angular
> dependence of reactive ion beam etching*, **J. Vac. Sci. Technol. A 10(4), 1212–1216
> (1992)**. DOI **10.1116/1.578229** (resolved via Crossref this session).

Barklund & Blom publisher abstract, verbatim (OpenAlex), final sentence:

> "**The angular dependence of physical sputtering with Ar and a more chemical etch with
> CF₄ of such a polymer, deposited from CHF₃, has also been studied.**"

Chang's description, p. 103, verbatim:

> "The silicon samples are coated with approximately **5000 Å thick polymer film formed from
> a CHF₃ plasma**. Argon, CF₄, and CHF₃ with 6 % O₂ are used in the RIBE system to etch the
> polymer film at various ion incident angles; the results are shown in Figure 4.16. **The
> angular dependence measured with argon is typical for a physical sputtering system with a
> maximum at 60° off-normal angle of incidence.** The etching yield measured with CF₄
> exhibits monotonic decrease with increasing angle of incidence…"

Digitised Fig. 4.16 (same method; axes from tick detection: `x` 1296 px = 0°, 2479 px = 100°;
`y` 1589 px = 0, 476 px = 2.0; markers by 5×5 opening + components). Series
identification by x-offset signature and by cross-check against the rendered image.

**Ar⁺ on CHF₃-deposited fluorocarbon polymer** (normalised etch rate, Barklund):

| θ | 0 | 20 | 40 | 50 | 60 | 65 | 70 | 75 | 80 |
|---|---|---|---|---|---|---|---|---|---|
| norm. rate | 1.000 | 1.098 | 1.168 | 1.197 | 1.348 | **1.448** | 1.348 | 1.194 | 0.891 |

**CF₄ on the same polymer** (the chemical channel, monotone):

| θ | 0 | 20 | 40 | 50 | 55 | 60 | 70 | 80 |
|---|---|---|---|---|---|---|---|---|
| norm. rate | 1.000 | 0.902 | 0.751 | 0.697 | 0.647 | 0.599 | 0.498 | 0.302 |

### 3.1 What this does to `B = 9.3`

petch applies `f1(θ) = (1 + 9.3 sin²θ)·cosθ` to the polymer/AC sputter channels
(`src/petch/mixed_layer.py:221-224`, `chemistry_deck.py:178`,
`amorphous_carbon_mask.py:284,290`, `boundary_transport_3d.py:775`), giving
**peak/normal = 4.17 at 52.6°**.

Two independent problems, one soft and one hard:

**(a) Magnitude (soft — convention-dependent).** Barklund's ordinate is labelled
"Normalized Etch Rate" and his abstract says "angular dependence of the etch rate".
Whether the ion flux was tilt-corrected is **not stated in the abstract and I could not
obtain the body** (AIP 403; no OA copy). Two readings:
  - already flux-corrected (Chang calls it "etching yield" in his prose) → measured
    FC-polymer enhancement is **1.45× at 60-65°**; petch is **2.9× too high**.
  - raw rate, `R = Y·cosθ` → implied `Y(60°)/Y(0°) = 1.348/0.5 = 2.70`,
    `Y(65°) = 3.43`; petch is ~1.2-1.5× too high.
  In both readings petch is **above** the only FC-polymer measurement, which is the same
  direction `RESEARCH_LIP_CERTAINTY` reached from the SiO₂ proxies (Cho 2000 ≈ 1.3,
  Schaepkens 1998 ≈ 1.33). `[VERIFY]` the Barklund flux convention against
  JVST A 10, 1212 body before quoting a single number.

**(b) Peak location (hard — convention-independent, and it is a form error).**
For `f(θ) = (1 + B sin²θ)cosθ`, writing `c = cosθ`, `f = (1+B)c − Bc³`, so
`f′ = 0 ⟹ c² = (1+B)/(3B)`. Hence:

| B | 0.6 | 1 | 3 | **9.3** | 50 | →∞ |
|---|---|---|---|---|---|---|
| peak θ | 19.5° | 35.3° | 48.2° | **52.6°** | 54.3° | **54.74°** |
| peak/normal | 1.006 | 1.089 | 1.778 | **4.172** | 19.8 | ∞ |

**The Kress form cannot place its maximum beyond 54.74° for any B > 0.** But the class-1
definition the whole lineage uses is "maximum **at/near 60°**" — Huang thesis L2293
("maximum at 60°"), Huard L2386 ("maximum near a 60° angle of incidence"), Chang on
Barklund's FC polymer ("maximum at 60° off-normal"), Chang on SiO₂/Ar⁺ ("maximum … at 60°
off normal"). The functional form petch uses is **structurally incapable** of matching the
stated class-1 peak position, independent of how B is tuned; raising B to move the peak
right runs the amplitude to infinity. This is a form-selection finding, not a
parameter-fitting one, and it is new.

*(Note the ordinate caveat: if the "60° peak" statements refer to **rate** and petch's
`f1` is a **yield**, the comparison is not apples-to-apples. But Huang/Huard describe
`P(θ)`, a per-particle probability — a yield — and Chang's own yields are flux-corrected
(§2.5 pt 3). So the mismatch stands under the lineage's own convention.)*

**Kress 1999 body / provenance of the literal 9.3: still NOT FOUND.**
Routes: AIP (403), OpenAlex (abstract only, retrieved verbatim — MD of **Cu and Ar on
Cu(111)**, 10-100 eV Cu / 50-250 eV Ar, IPVD context), Semantic Scholar (abstract elided
by publisher), Unpaywall/OA (none). The number 9.3 appears in **no** fetched source.

---

## 4. SiO₂-complex `E_th = 35 eV` — **PROVENANCE FOUND** (and the repo's audit was incomplete)

`RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` L693 records: *"The only 35 eV in Chang's
corpus is the **Ar⁺-on-Si physical-sputter** threshold (Table 3.1). Chang's SiO₂
measurements (Table 4.2) are at 70/100 eV and quote no complex threshold."*

Both halves check out — and there is a **second** 35 eV that the audit missed, and it is
the chemistry-matched one.

**Confirmed, Table 3.1, p. 73** (`chang_thesis.txt` L3109-3131), Ar⁺ on Si physical
sputtering, universal (Steinbrüchel) energy dependence:

| A (eV^−½) | E_th (eV) | Reference |
|---|---|---|
| 0.04 | **35** | Chang (this work) |
| 0.04 | 27 | Harper 1981 |
| 0.04 | 31 | Tachi 1986 |
| 0.03 | 64 | Oostra 1988 |

> "The extrapolated threshold energy from this work is approximately 35 V, which is higher
> than the 20 eV threshold energy reported for physical sputtering silicon with Ar⁺ by Gray
> [Gray, 1993]."

**Confirmed, Table 4.2, p. 96**: SiO₂ model parameters listed only at `E_ion` = 70 and
100 eV (`Y₀` 0.01/0.02, `P` 0.04/0.08, `s` 0.001/0.005) — no threshold column. Chang's own
SiO₂ physical-sputter threshold is **40 eV** (Fig. 4.1, p. 90; and again p. 91 for Ar⁺+Cl).

**THE MISSED ONE — Chang thesis §4.1, p. 90** (`chang_thesis.txt` L4125-4130), verbatim:

> "The extrapolated threshold energy is approximately 40 eV, **which is slightly higher than
> a 35 eV threshold energy reported for etching of silicon dioxide in a CHF₃ plasma
> [Joubert, 1994].**"

Reference resolved (`chang_thesis.txt` L7761 + Crossref):
**O. Joubert, G. S. Oehrlein, M. Surendra**, *Fluorocarbon high density plasma. VI.
Reactive ion etching lag model for contact hole silicon dioxide etching in an electron
cyclotron resonance plasma*, **J. Vac. Sci. Technol. A 12(3), 665–670 (1994)**,
DOI **10.1116/1.578850**.

**Assessment.** A **35 eV threshold for SiO₂ etching in a fluorocarbon plasma**, reported in
a paper that is *specifically about contact-hole RIE lag*, is a far better provenance for
petch's `complex_sio2_yield` (`p0 = 0.1384/0.1471, ε_th = 35 eV, q = 1, ε₀ = 140`) than the
Ar⁺-on-Si coincidence. It is chemistry-matched (fluorocarbon/SiO₂), geometry-matched
(contact hole), and it is reachable from Chang's own text — which plausibly explains how a
"Chang–Sawin" label got attached to it in the Appendix-B lineage.

Remaining gap: I have Chang's **verbatim third-party attribution** of the 35 eV to Joubert,
not Joubert's own printed number. The Joubert abstract (OpenAlex, verbatim, retrieved) does
not contain it — it is a body number. `[VERIFY]` the Joubert 1994 body (AIP paywall; no OA
copy found). Recommend re-labelling the deck row **"35 eV — Joubert/Oehrlein/Surendra 1994
via Chang 1998 p. 90"** rather than "Chang–Sawin", which is demonstrably wrong.

---

## 5. Izawa 2007 body — **NOT FOUND** (all routes exhausted)

M. Izawa, N. Negishi, K. Yokogawa, Y. Momonoi (Hitachi), *Investigation of Bowing Reduction
in SiO₂ Etching Taking into Account Radical Sticking in a Hole*, **Jpn. J. Appl. Phys.
46(12R), 7870 (2007)**, DOI **10.1143/JJAP.46.7870**.

Routes tried this session (on top of the prior session's IOPscience/ResearchGate/J-STAGE):
- **Semantic Scholar API** (`DOI:10.1143/JJAP.46.7870`): record exists
  (`paperId 668d0ac0…`), `openAccessPdf.status = CLOSED`, **abstract field elided by
  publisher**.
- **Unpaywall** (`api.unpaywall.org/v2/...`): `is_oa = false`, `oa_status = closed`,
  `oa_locations = []`.
- **OpenAlex**: abstract available (matches the abstract already recorded verbatim in
  `RESEARCH_LIP_CERTAINTY_2026-08-04.md` §1 word-for-word — independent confirmation that
  the recorded abstract is accurate), no OA URL.
- J-STAGE search: no record.

**Status unchanged: abstract-only.** The `RESEARCH_LIP_CERTAINTY` verdict (do not import
0.004) does not depend on the body, and the three structural objections raised there
(model-inverted fit, net-vs-gross ambiguity, F-rich end of a 125× split) are unaffected.

---

## 6. Corrections and downgrades this session produces

| Location | Current text | Correction |
|---|---|---|
| `RESULTS_ML19_ENDPOINT_2026-08-05.md` L103; `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md` L228, L248 | bond multiplicity `[VERIFY]`, sensitivity spanned "m = 3 … m = 8" | **m = 8 is excluded.** Source rule: m = number of available bonds; **CF₂ → 2, CF₃ → 1**, Fig-5 example = 3. Cite JVST A 42, 043008 (2024) §IV. |
| `RESULTS_ANGULAR_CLASSES_2026-08-05.md` §"The class-2 roll-off is `[VERIFY]`" | roll-off shape unverified; "Chang & Sawin 1997 is paywalled (fetch returned 403)" | **Curve obtained** via Chang MIT thesis Fig. 5.9 (p. 115). `min(1, cosθ/cos45°)` agrees to ≤0.065 absolute; worst point 70° (12 % low). Downgrade to quantified approximation with the §2.4 table as receipt. |
| `RESEARCH_LIP_CERTAINTY_2026-08-04.md` §2.4(b) | "There is **no measured FC-film angular sputter yield** in the literature" | **False.** Barklund & Blom, JVST A 10, 1212 (1992), DOI 10.1116/1.578229 — Ar⁺ and CF₄ RIBE on a 5000 Å CHF₃-plasma polymer film. Reprinted as Chang Fig. 4.16, digitised in §3. |
| `RESEARCH_LIP_CERTAINTY_2026-08-04.md` §2.4 | B = 9.3 flagged only as off-domain (Cu) and ~3× above in-chemistry bounds | Add the **form** finding: `(1+B sin²θ)cosθ` peaks at ≤ **54.74°** for all B, so it cannot reproduce the lineage's own "maximum at 60°" class-1 definition at any B. |
| `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` L693, L726 | "The only 35 eV in Chang's corpus is the Ar⁺-on-Si physical-sputter threshold" | **Incomplete.** Chang p. 90 also cites *"a 35 eV threshold energy reported for etching of silicon dioxide in a CHF₃ plasma [Joubert, 1994]"* = Joubert, Oehrlein & Surendra, JVST A 12, 665 (1994), DOI 10.1116/1.578850. Re-label the deck row. |
| `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` §T-1 | Chang thesis listed as a target | **Acquired.** `1721.1/50356`, direct bitstream URL in §0. |

---

## 7. Still open after this pass

1. **Krüger's crosslinking probability** — not published anywhere (paper, thesis,
   Appendix B, supplementary-by-proxy). Only source would be the authors.
2. **Barklund & Blom 1992 flux-normalisation convention** — decides whether the measured
   FC-polymer peak enhancement is 1.45× or 2.7×. AIP paywalled, no OA copy. Highest-value
   remaining retrieval for the lip law.
3. **Kress 1999 body / origin of the literal 9.3** — still unsourced. (Lower priority now:
   §3.1(b) shows the *form* is wrong for the stated class-1 peak regardless of B.)
4. **Chang's `a_i` polynomial coefficients** (Eq. 4.1 / Fig. 5.9 fit) — never printed; the
   digitised table in §2.4 is the best available rendering.
5. **Joubert 1994 body** — to convert the 35 eV from a verbatim third-party attribution
   into a primary number.
6. **Izawa 2007 body** — abstract-only, all open routes exhausted; would need ILL/paid.
