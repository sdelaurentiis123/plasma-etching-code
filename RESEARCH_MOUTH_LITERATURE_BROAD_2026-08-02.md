# Broad literature survey — what holds the mask-opening / mouth / top-CD in HAR fluorocarbon dielectric etch

Date: 2026-08-02. Author-of-record: autonomous research agent (broad lane).
**Do NOT commit** (per task instruction). Companion to `RESEARCH_MOUTH_EQUILIBRIUM_2026-07-27.md`.

**Deliberate non-overlap:** a parallel agent owns the Krüger 2024 thesis mining. This document
covers *everything else*: other theses, published experiments, mechanism literature, and the
cross-check on whether our specific failure (simulated mouth over-narrows to ~14–25 nm vs 45 nm
experimental at 60 s) has been seen and diagnosed before. Krüger's own numbers are referenced
only where a third-party source cites them.

**Our failure being investigated:** depth reproduced (852 sim vs 825 nm exp, 60 s, Ar/C4F6/O2 CCP,
SiO2 with a-C mask, ~90 nm line opening → 45 nm experimental at 60 s), but simulated mouth
over-narrows to ~14–25 nm. Suspects on entry: (a) mask corner faceting funnelling ions,
(b) in-feature charging, (c) product redeposition, (d) polymer deposition/erosion equilibrium
at the lip.

Citation convention: **verbatim quotes are in quotation marks and were read from a file I fetched
this session** (path given). Publisher-supplied abstracts retrieved via the Crossref/OpenAlex APIs
are labelled *[abstract, publisher metadata]* — these are verbatim but are the abstract, not the
body. `[VERIFY]` marks anything I could not read directly.

---

## 0. Executive answer (read this first)

**The literature is nearly unanimous that the mouth/neck is a POLYMER DEPOSITION–REMOVAL BALANCE
at the lip, and that the mask facet is a SECOND-ORDER effect on the mouth itself — the facet's
first-order role is BOWING (further down the sidewall), not the mouth.** Two independent
industrial groups say this in almost identical words:

- **Lam Research (Kim, Hudson, Cooperberg, Edelberg, Srinivasan, *Thin Solid Films* 515, 4874
  (2007)):** necking ← net polymer deposition rate on the sidewall; bowing ← ion scattering off
  the **secondary** facet; and — the load-bearing negative result for us — **the *primary*
  photoresist facet angle "showed only a small influence on the SiO2 etch profile."**
- **Lam Research 2023 (Shen, Lill *et al.*, JJAP 62, SI0801):** "Polymer deposition on the
  sidewalls of the mask, so-called **necking**, is considered the main root cause of sidewall
  roughness."

The dissenting/complementary experimental claim comes from **Seoul National University + Ajou
(Lee *et al.*, *J. Electrochem. Soc.* 157, D142 (2010))**, who used a **Faraday-cage** experiment
to attribute necking to **redeposition of material sputtered off the mask slope**, and bowing to
**ions reflected off the mask slope** — i.e. mask-facet-mediated *redeposition* is a real,
measured necking channel, distinct from direct plasma polymer flux.

**Ranked answer to "what holds the mouth open", per the literature:**

| Rank | Mechanism | Strength of evidence | Direction on our 14–25 nm |
|---|---|---|---|
| **1** | **Angle-resolved polymer deposition vs angle-resolved etch/sputter at the lip** (the NDR/NER balance). Measured with a Faraday cage; the *whole* narrowing-vs-bowing axis is controlled by how fast deposition falls off with incidence angle relative to how fast the etch yield falls off. | **Strongest — directly measured** (You *et al.*, *Coatings* 13, 1452 (2023); Izawa *et al.*, JJAP 46, 7870 (2007) sticking coefficients). | If our film's deposition is too angle-flat, or our sputter yield falls off too fast off-normal, the mouth over-narrows. **Primary suspect.** |
| **2** | **Angle-resolved *removal* at the grazing lip (sputter/hot-neutral).** NER measured **above** the cosine curve out to 50–60° ⇒ real physical sputtering is *stronger* than cosine at the lip. | **Strong — measured** (You 2023, Fig. 6b). | If we deliver only ∝cosθ removal to the lip, we systematically under-remove there → over-narrow. **Directly matches our ml9b "reflection seals the mouth" symptom.** |
| **3** | **Redeposition of mask-sputtered material onto the lip/upper sidewall** (a *narrowing* channel, not an opening one). | **Strong — Faraday-cage isolated** (Lee 2010). | Wrong sign for us: adding it narrows further. But it means a *self-consistent* model needs it AND a compensating removal. |
| **4** | **Mask facet funnelling / ion reflection off the facet.** | Strong for **bowing**; explicitly **weak for the mouth** (Kim/Hudson: primary facet "only a small influence"). | Adding faceting will mostly widen the *bow* region, not the mouth. **Deprioritize as the mouth fix.** |
| **5** | **Mask erosion relaxing the neutral-conductance limit** (an indirect opener). | Moderate — Huang/Kushner 2019 states it explicitly. | Our 133 nm mask loss should *already* be helping; the fact that it isn't is diagnostic. |
| **6** | **Charging.** | Real for mask faceting/bowing/microtrench (Zhai 2025) and for twisting; **no source found claiming charging sets the mouth width.** | Low priority for the mouth. |

**Has our exact failure been diagnosed before? Partially, twice** — see §4. The closest match is
**Zhang (Michigan, 2015)**: MCFPM-3d reproduced necking + bowing *phenomenology* but **"does not
precisely reproduce the positions of the necking and bowing effect"**, and he attributes it to the
**ion angular distribution**, not to chemistry: *"a slight change in angular distribution may
contribute significantly to different shape evolutions."* No source I found reports a simulated
mouth that seals when the experiment stays open, so our specific magnitude failure appears to be
**undiagnosed in the open literature** — which is consistent with it being an implementation-side
angular-delivery defect rather than a missing chemical channel.

---

## 1. THESES — fetched, mined, and scored

All Kushner-group theses were enumerated from the openly-served directory index
`https://cpseg.eecs.umich.edu/pub/theses/` (fetched with `curl -k`, HTTP 200, 57 entries;
saved parse at `/private/tmp/claude-501/-Users-stanislavdelaurentiis-chip-etch/352207a8-b462-4fe3-a650-c9b37f19523b/scratchpad/theses_index.html`).
**Note:** `huang_shuo_phd_thesis.pdf` on that server is byte-identical (15,584,588 B) to the
repo's existing `tmp/pdfs/huang_thesis.pdf` — the Huang thesis we already have. Skipped per task.

### 1.1 Theses successfully fetched this session

| Thesis | Path | Mask-mouth relevance |
|---|---|---|
| **Wang, Mingmei — *Modeling of complex surface interactions in low and high pressure plasmas*, PhD, Iowa State / Kushner (undated in file; MCFPM + HPEM, Ar/C4F8/O2)** | `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/wang_mingmei_phd_thesis.pdf` (+`.txt`) | **HIGHEST.** Entire thesis chapters on mask (PR) erosion → facet → bowing, and on *preventing* it. |
| **Zhang, Yiting — *Low Temperature Plasma Etching Control through Ion Energy Angular Distribution and 3-Dimensional Profile Simulation*, PhD, U. Michigan, 2015** | `.../tmp/pdfs/zhang_yiting_phd_thesis.pdf` (+`.txt`) | **HIGH.** MCFPM-3d vs Lam-ICP SEM time series; explicit admission of necking/bowing **position** mismatch. |
| **Huard, Chad M. — *Nano-Scale Feature Profile Modeling of Plasma Material Processing*, PhD, U. Michigan, 2018** | `.../tmp/pdfs/huard_chad_phd_thesis.pdf` (+`.txt`) | MEDIUM. Bowing ← off-axis ions + specular reflection off the bowed surface; quantified vs IAD width. Deliberately *suppresses* mask erosion in his ALE studies. |
| **Song, Sang-Heon — *Control of plasma kinetics for microelectronics fabrication*, PhD, U. Michigan, 2014** | `.../tmp/pdfs/song_sangheon_phd_thesis.pdf` (+`.txt`) | MEDIUM-LOW. Contains a printed FC surface-probability set (see §1.5) but deliberately makes PR inert. |
| Tian, Peng — *Controlling Photon and Ion Fluxes in Low Pressure Low Temperature Plasmas*, PhD, U. Michigan, 2018 | `.../tmp/pdfs/tian_peng_phd_thesis.pdf` | LOW — reactor/photon-flux thesis; **zero** facet/mask-opening hits. |
| Qu, Chenhui — *Computational Investigations of Fundamental Plasma Processes in Semiconductor Industrial Applications*, PhD, U. Michigan, 2020 | `.../tmp/pdfs/Qu_Chenhui_PhD_Thesis_2020.pdf` | LOW — HAR mentioned only in framing; **zero** facet hits. |
| Konina, Kseniia — *Atmospheric Pressure Plasma Treatment of Complex Interfaces*, PhD, U. Michigan, 2024 | `.../tmp/pdfs/Konina_Kseniia_PhD_Thesis_2024.pdf` | **NOT RELEVANT** (atmospheric-pressure plasma–liquid). Fetched and ruled out. |
| Lanham, Steven — *Modeling Pulsed Power Plasmas and Applications to In-situ Nanoparticle Growth*, PhD, U. Michigan, 2022 | `.../tmp/pdfs/Lanham_Steven_PhD_Thesis_2022.pdf` | NOT RELEVANT. Fetched and ruled out. |
| Logue, Michael — *Control of Electron and Ion Energy Distributions in ICPs…*, PhD, U. Michigan | `.../tmp/pdfs/logue_michael_phd_thesis.pdf` | NOT RELEVANT. Fetched and ruled out. |

Also fetched (papers, not theses): see §2/§3 for paths.

### 1.2 Wang Mingmei thesis — the facet→bowing lineage, and what it does NOT claim

**Does their simulated mask facet?** Yes, explicitly and centrally. **What holds their mouth?**
*Nothing is claimed to* — Wang's mask (PMMA photoresist) is expected to erode; the thesis is about
*slowing* that erosion, and the failure mode she tracks is **loss of feature definition**, not a
sealed mouth.

Verbatim (`wang_mingmei_phd_thesis.txt`, abstract):
> "Bowing at the top of the feature near the mask-SiO2 interface mainly results from bombardment
> of ions reflecting from eroded mask surface."

Verbatim (§1.3):
> "Bowing is attributed to the change in the acceptance angle of ions into the feature from the
> plasma due to erosion of the photoresist (PR, usually is a C-H based polymer) mask and subsequent
> reflection of ions from the facets of the PR."

Verbatim (§1.3, mechanism + its *timing*):
> "The PR initially usually has a domed shape that does not reflect many energetic ions into the
> trench side walls since the PR is thick compared to the trench size. As etching proceeds, the PR
> is eroded. **When the PR is thin enough, bowing occurs due to reflection off the facets.**"

**Citation chain closed:** Wang's reference for that mechanism is **[59] = J.-K. Lee, H.-Y. Jang,
S.-H. Lee *et al.*, *J. Electrochem. Soc.* 157, D142 (2010)** and her general bowing references are
**[56] = N. Ikegami, A. Yabata, T. Matsui *et al.*, *Jpn. J. Appl. Phys.* Part I **36**, 2470 (1997)**
and **[57] = D. Kim and E. A. Hudson, *Thin Solid Films* **515**, 4874 (2007)** (verbatim from her
bibliography, `wang_mingmei_phd_thesis.txt` L695–700). So the whole Kushner-lineage
facet-reflection-bowing story descends from the two experimental/industrial papers analysed in §3.

**Numbers extractable from Wang (all verbatim):**
- Feature: "The mask opening is **75 nm** wide with a depth of **1500 nm** to a Si stop layer,
  yielding an aspect ratio of **AR = 20**. The photoresist (PR) is initially **dome-shaped** with a
  maximum thickness of **450 nm**." Chemistry Ar/C4F8/O2, dc-augmented dual-frequency CCP.
- "For incident ion energies < **1.3 keV**, polymer deposition effectively stops etching at the Si
  layer." (an *etch-stop* threshold, i.e. the clog boundary in energy).
- Mask hardening: "In the current model, the etching yield for **cross-linked PMMA is 5 times
  smaller** than for non-crosslinked PMMA." Selectivity SiO2/PMMA: **≈10 without cross-linking,
  ≈17 with**; vs ion energy with cross-linking: **S ≈ 3 at 100 eV, ≈13 at 500 eV, ≈16 at 1000 eV**.
- "10-nm-thick PR material (PMMA) is eroded if cross-linking is not considered… There is still
  **5-nm-thick PMMA eroded**" with cross-linking.
- Concept she defines: **LFAR** (limiting feature aspect ratio) = "that depth at which the PR mask
  is fully eroded so that the feature is no long[er] defined"; LFAR rises ~linearly with PMMA
  thickness "up to 25 nm".

**Relevance to petch's PC-sputter-resistance question (our §5b in the July doc):** Wang gives a
concrete, published magnitude for the crosslinked-vs-virgin mask sputter ratio — **5×** — for PMMA
under ion bombardment + VUV. That is an independent order-of-magnitude anchor for the
`[VERIFY vs thesis Appendix B]` PC-vs-P factor in our mixed-layer model, *for a C–H resist*, not
for a-C. `[VERIFY]` before lifting it to a-C.

**Neck/bow-vs-time figures?** Wang has profile-vs-time *sequences* (Fig. 6.5b: "the PMMA continues
to be eroded until feature definition is lost at **96 s**") but **no wm(t) / top-CD(t) trace.**

### 1.3 Zhang Yiting thesis (Michigan 2015) — the only found statement of a necking/bowing model-vs-experiment MISS

**Does their simulated mask facet?** Yes — mask erosion is modelled and shown in SEM+simulation
side by side. **What holds their mouth?** Not addressed; the metric they track is **CDR** (critical
dimension ratio), not the mouth.

Verbatim (`zhang_yiting_phd_thesis.txt`, §7.3 Model Validation, He/Cl2 Si trench, **Lam Research
ICP**, line/pitch 50/100 nm, 60 nm oxide + 60 nm nitride mask, ∆x=∆y=∆z=1.25 nm):
> "The masks show erosion with increasing etch time, which can be seen in the measured SEMs …
> With the thickness of mask continuing to decrease, ions with large horizontal velocities will
> bombard the sidewalls of the feature, and thus causes sidewall etching. After ions strike on the
> surface, there will be high energy neutrals reflecting back to the plasma and bombarding the
> surface again. **This high energy particle reflection brings about the necking and bowing effect**
> as observed in the third column of Fig. 7.2. **There is a difference of necking and bowing
> positions between the experimental measurements and the predicted simulation results.** This is
> mainly due to the absence of reactor scale measurements to validate the ion energy and angular
> distributions and flux ratios. Simulation results in Sec. 7.4 reveal that **a slight change in
> angular distribution may contribute significantly to different shape evolutions.**"

And the conclusion:
> "the MCFPM 3-d **does not precisely reproduce the positions of the necking and bowing effect**,
> which suggests further experimental validation in the prediction of the IEADs and fluxes on the
> wafer or using measured IEADs and fluxes as inputs in the MCFPM 3-d."

Also verbatim, the **mask-erosion → reflection-point-shift** statement, which is the mechanistic
core of "the mouth geometry is set by where reflections land":
> "When the mask is sputtered by high energy particles, **the ion reflection position shifts and
> results in a bowing shape under the mask.**"

His HAR SiO2 metric definition (useful for our scorecard vocabulary): **"CDR is the width at the
center of the feature to the mask opening"**, with **"The desired value of CDR is 1.0 – tapered
profiles have CDR < 1 and bowing profiles have CDR > 1."** Feature: "The width of the mask opening
is **37 nm** and the aspect ratio is **15**. The over-etch was 20%." Ar/CF4/O2 = 75/20/5, 30 mTorr.
He also reports an **etch-stop by polymer** at low power: "With only 300 W at both the HF and LF,
the large flux of low energy ions results in **excessive polymer deposition on sidewalls**. With
this polymer build up, **an etch stop occurs** before reaching the underlying Si."

**Neck/bow-vs-time figures?** Yes for the *validation* case (Fig. 7.2 is a measured-SEM ×
simulated time series with necking and bowing). Not for the fluorocarbon HARC case.

### 1.4 Huard thesis (Michigan 2018) — quantified bowing vs IAD width

Verbatim (`huard_chad_phd_thesis.txt`):
> "Narrower angular distributions produce more tapered side wall profiles than the base case
> (**18% reduction in width at half etch height** for σi = 0.25 compared to the base case). At the
> other extreme, the profile resulting from σi = 1.5 has increased side wall bowing (**30% increase
> in width at half etch height** compared to the base case). The broader ion angular distributions
> enable sites having off normal view-angles to the plasma to intercept a larger fraction of the
> incoming ion flux, the first consequence being **bowing under the mask**. Once this initial bowing
> occurs, **ions specularly reflecting from the bowed surface are more likely to strike the side
> wall deeper in the trench.**"
(Base-case average ion angular spread 2.2°; σi = 0.25 → 0.55°, σi = 1.5 → 3.3°.)

**This is the single most useful *sensitivity* number in the thesis corpus for us:** a ~6× change in
the angular spread moves the mid-depth width by only −18%/+30%. It says **IAD width alone cannot
move a 45 nm mouth to 22 nm** — i.e. our mouth deficit (≈50%) is too large to be an IAD-width
artefact and must come from the **removal law at the lip**, not the source distribution.

Note also Huard deliberately removes the mask from the physics when he wants clean ALE data:
> "Etching of the top surface is prevented by a thin hard mask **to avoid complications caused by
> mask erosion and energetic particles reflecting from the mask.**"
— i.e. the group treats mask erosion + mask reflection as a *coupled pair of complications*, never
one without the other.

### 1.5 Song thesis (Michigan 2014) — a printed FC probability set with a deliberately inert mask

Verbatim (`song_sangheon_phd_thesis.txt`):
> "The probability of reactions with photoresist (PR) was chosen to be **small enough in order to
> eliminate the effect of PR mask erosion on the etch profile**. The sputtering probability of the
> polymer by ion has been assumed to be **20%**. We also considered polymer deposition on top of the
> polymer layer and the sputtering probability for this kind of polymer is assumed to be **25%**,
> which is a little bit larger than that of a normal polymer. The sputtering probability for
> activated SiO2 by ions is assumed to be **90%**. … CF, CF2 and CF3 have **2%, 1%, and 0.3%** polymer
> deposition probabilities on the chamber wall, respectively."

Two things matter here: (i) an *independent* Kushner-lineage polymer-sputter probability
(**0.20–0.25**) that is **2–3× larger** than the Krüger-lineage numbers our mixed-layer uses
(complex sputter 0.1384, bare 0.0909); (ii) a second confirmation that when a Kushner-lineage
worker wants a clean profile study, the standard move is to **turn the mask off**, which is why the
literature is thin on quantitative mask-mouth validation.

### 1.6 Groups asked about where I found NO usable thesis-level mouth/CD source

- **Bogaerts (Antwerp):** no feature-scale HAR fluorocarbon mask-evolution thesis surfaced in any
  search; that group's HAR-relevant output is gas-phase/plasma-chemistry, not profile evolution.
  **Negative result, moderate confidence.** `[VERIFY]` by a direct repository sweep of
  `repository.uantwerpen.be` if this matters.
- **Bochum/RUB (Schulze students):** searched; surfaced only patents and unrelated work. RUB's HAR
  contribution is **reactor-scale (voltage-waveform-tailoring / PIC)**, and the feature-scale
  half of that collaboration is Kushner's (e.g. *Phys. Plasmas* 31, 033508 (2024), local copy at
  `…/scratchpad/PhysPlasmas_31_033508_2024.pdf`). **Negative result.**
- **Hamaguchi (Osaka):** searched; his HAR-relevant modern output that surfaced is
  sputtering/damage MD and *Phys. Plasmas* 33, 033901 (2026) "Computational study of the evolution
  of high-aspect-ratio SiO2 trench and hole features during **physical sputtering**" — i.e.
  *sputtering-only*, no fluorocarbon polymer, so **no mouth equilibrium**. `[VERIFY]` the 2026
  paper's mask treatment if a sputter-only mouth reference is wanted.
- **Kokkoris / Gogolides (NTUA / Demokritos):** their canon is ARDE/RIE-lag/etch-stop and
  the *coverage-balance* surface model — Kokkoris, Gogolides, Boudouvis, *J. Appl. Phys.* **91**,
  2697 (2002), DOI 10.1063/1.1435833; Gogolides *et al.*, *J. Appl. Phys.* **88**, 5570 (2000).
  Their published treatment states charging is handled only "by an increased ion angular spread"
  and **the mask is not evolved**. So: **relevant for the clog boundary, not for the mouth.**
- **TU Wien (Filipovic/Weinbub, ViennaPS):** *is* relevant and *is* fetched — see §3.4.

---

## 2. PUBLISHED EXPERIMENTS — SEM cross-sections, mask evolution, and numbers

### 2.1 The single best quantitative mouth/neck dataset found: You, Yang, Jeon, Chae, Kim, *Coatings* 13, 1452 (2023)

**Full citation:** S. You, H. S. Yang, D. Jeon, H. Chae, C.-K. Kim, "Controlling Bowing and
Narrowing in SiO2 Contact-Hole Etch Profiles Using Heptafluoropropyl Methyl Ether as an Etchant
with Low Global Warming Potential," *Coatings* **13**(8), 1452 (2023). DOI 10.3390/coatings13081452.
Ajou University + SKKU (SAINT). **Open access; fetched:**
`/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/coatings2023_bowing_narrowing.pdf`
(+ `.txt`). ICP, 13.56 MHz source + 13.56 MHz bias, ACL mask, HFE-347mcc3/O2/Ar.

This paper does exactly what we need: it **measures hole diameter vs vertical position** (a
CD-vs-depth trace) across a gas-ratio sweep that walks the profile from **bowing → anisotropic →
narrowing → etch stop**, and it explains the whole axis with **two measured angular functions**.

Verbatim numbers (200 nm-diameter holes, vertical position measured from the **ACL/SiO2
interface = 0**, downward negative):
- **HFE/Ar = 0.40 (8/2/20 sccm):** "the top diameter of the SiO2 hole slightly decreased to
  **197 nm** (from **200 nm** before etching) and **bowing** of the hole was observed."
- **HFE/Ar = 0.47 (9/2/19 sccm):** "**narrowing** of the hole instead of bowing was observed. …
  the hole diameter decreased from **197 nm** at zero vertical position to **178 nm** at a vertical
  position of **−661 nm** and then it increased to **194 nm** at a vertical position of **−1094 nm**."
- **HFE/Ar = 0.56 (10/2/18 sccm):** "the hole diameter shrank to **144 nm** at the vertical position
  of −661 nm. Although the hole diameter increased to **187 nm** at −1094 nm, it decreased again
  rapidly … resulting in a bottom diameter of only **75 nm**."
- **HFE/Ar ≥ 0.65:** "The continuous narrowing inhibits the contact hole from reaching the bottom,
  ultimately resulting in an **etch stop**."
- Best case: "a highly anisotropic and bowing-free 100-nm-diameter contact hole profile"
  with **AR = 24** at 9/2/19 sccm.

**The mechanism, measured with a Faraday cage (this is the load-bearing physics for our mouth):**
> "**Figure 5a** shows the change in the **deposition rates of the fluorocarbon films on SiO2 with
> the ion-incident angles** … The deposition rates of fluorocarbon films determined in this study
> **decreased monotonically with increasing the ion-incident angle** under all conditions …
> Although the deposition rates of the fluorocarbon film decreased with increasing the ion-incident
> angle, **the degree to which the deposition rate was reduced depended on the** HFE-347mcc3/Ar
> **ratio**."

> "The NDR is defined as the deposition rate at a specific angle **normalized with respect to the
> deposition rate on the horizontal surface**. The dotted line in the NDR plot represents a **cosine
> curve**. The NDRs … clearly demonstrated that **the extent of reduction in the deposition rate of
> the fluorocarbon film with ion-incident angle decreased as the** HFE/Ar **ratio increased**. This
> implies the **etch resistance of SiO2 contact holes on the slanted sidewalls increased** with
> increasing … ratio. Therefore, **the contact holes were less etched on the slanted sidewalls at
> higher** … ratios, **leading to the narrowing of the holes.**"

> "In all conditions, **the NERs are above the cosine curve until the ion-incident angle reaches
> 50–60 degrees. This indicates that physical sputtering plays an important role during etching.**
> As the HFE/Ar ratio increased, **the NER decreased more rapidly with increasing the ion-incident
> angle.** The difference in the NERs at low and high ratios … **reached a maximum value at 60°.**
> This implies that changes in the … ratio **primarily impact the etch capability of SiO2 contact
> holes on the slanted sidewalls rather than the bottom plane.**"

**Why this is the most actionable result in this document:** the narrowing/bowing axis is fully
determined by **(NDR(θ) − something) vs (NER(θ))** — two *measured* angular functions with the
**cosine curve as the explicit reference**. If petch's lip removal is delivered ∝ cosθ (or worse,
attenuated at grazing), we sit on the "high-ratio, narrowing" side of this measured axis by
construction. **This is a direct, digitizable, knob-retiring target: NDR(θ) and NER(θ) families
at 5 gas ratios, with the cosine reference drawn on the same axes** (their Figs. 5a/5b, 6a/6b).

### 2.2 Lee, Jang, Lee, Kim, Moon, *J. Electrochem. Soc.* 157, D142 (2010) — mask-slope redeposition ⇒ necking; mask-slope reflection ⇒ bowing

**Full citation:** J.-K. Lee, I.-Y. Jang, S.-H. Lee, C.-K. Kim, S. H. Moon, "Mechanism of Sidewall
Necking and Bowing in the Plasma Etching of High Aspect-Ratio Contact Holes," *J. Electrochem.
Soc.* **157**(3), D142 (2010). DOI 10.1149/1.3276511. C4F6/CH2F2/O2/Ar. **Full text NOT obtained**
(IOP/ECS paywall returned a login page; the stub was deleted). Abstract obtained verbatim from
OpenAlex *[abstract, publisher metadata]*:

> "The mechanism of sidewall necking and bowing during the etching of high aspect-ratio SiO2
> contact holes in a C4F6/CH2F2/O2/Ar plasma was investigated by **monitoring the etch profiles over
> etching time**. As the etching proceeded, **the sidewall necking and bowing became severe, and the
> positions of the necking and bowing moved toward the bottom of the hole**. Although **the mask
> slope, which was responsible for the deviation of the ion direction, was already formed after
> etching for 4 min**, the sidewall bowing appeared **after 9 min** when the necking moved to
> positions below the mask-SiO2 boundary. This suggested that the sidewall bowing resulted from the
> combination of two phenomena: (i) **the protection of the top region of the sidewall by a thick
> CFx film deposited at the necking position** and (ii) the lateral etching of the sidewall below
> the necking position. The effect of the mask slope on the formation of the sidewall necking and
> bowing was examined **using a Faraday cage system**. **The redeposition of particles sputtered from
> the mask slope on the contact-hole sidewall resulted in sidewall necking**, and **the secondary
> etching of the sidewall by ions reflected from the mask slope contributed to the formation of
> sidewall bowing**."

**Four load-bearing facts for us:**
1. **Necking has a *time* and a *place*, and both MOVE DOWN.** The neck is not pinned at the mask
   plane; it migrates below the mask–SiO2 boundary. If our harness measures the mouth at a fixed
   height, or if our neck cannot migrate, we are measuring a different quantity. (This reinforces
   item #5 in the July doc's ranked table.)
2. **The mask slope forms EARLY (4 min) and bowing appears LATE (9 min).** So facet formation is
   *not* rate-limiting for the profile — a **direct experimental corroboration of Kim/Hudson's
   "primary facet has only a small influence."**
3. **Redeposition from the mask slope is a measured, isolated necking channel** — mechanism (d) in
   our suspect list is real, but it is a **narrowing** channel. Adding it makes our mouth worse.
4. **The neck region's thick CFx film *protects* the top of the sidewall.** So in the experiment,
   the mouth region is polymer-covered and *stable*, while the etch action moves below it. Our
   model turns that protective film into a *closure*.

### 2.3 Miyake, Negishi, Izawa, Yokogawa, Oyama, Kanekiyo, *Jpn. J. Appl. Phys.* 48, 08HE01 (2009) — mask TAPER ANGLE controls bowing (Hitachi)

DOI 10.1143/jjap.48.08he01. Institution: Hitachi (Japan). **Full text NOT obtained** (paywall).
Abstract verbatim *[abstract, publisher metadata]*:

> "The evaluation of etching profiles produced with **different taper angle masks** confirmed that
> **the bowing amount and mask selectivity worsened with decreasing mask taper angle**. The
> relationship between mask taper angle and distribution of scattered ion flux on the sidewall of
> a tapered mask was calculated. **The scattered ion flux was heavily concentrated in the upper part
> of the sidewall in the case of a tapered mask, and this was considered to be the main cause of
> the bowing formation.** Direct observation of an etched sidewall by atomic force microscopy (AFM)
> revealed that **the roughness of the necking was strongly related to the roughness of the bottom
> part of the etched sidewall**. … in the case of **nonaxisymmetric necking**, an imbalance of ion
> flux in the bottom of the hole appeared and broke the etching symmetry … causing **twisting**. In
> addition, **the probability of twisting was found to increase with increasing necking growth rate
> irrespective of mask electrification.** Therefore, **mask deformation and nonuniform necking in
> the upper part of the sidewall during HARC etching are considered the main factors causing bottom
> degradation. Accordingly, a vertical and nondeformed mask is very important for a smaller
> critical dimension (CD) and HARC etching.**"

**Two things matter:** (i) this is a **mask-taper-angle DOE against SEM** — the closest thing in the
open literature to a facet-angle experiment; the actual angle values live in the paywalled body
`[VERIFY — get the taper-angle values]`. (ii) The phrase **"irrespective of mask electrification"**
is an explicit *negative* result on charging: they varied/considered mask charging and found the
necking growth rate, not charging, controlled twisting.

### 2.4 a-C / ACL mask evolution — measured behaviour and the numbers that exist

**(a) Kim, Cho, Kim, Jhon, Min, Kim, Yeom, *J. Vac. Sci. Technol. A* 31, 021301 (2013)** —
"Study on the etching characteristics of amorphous carbon layer in oxygen plasma with carbonyl
sulfide." **Samsung Memory Division + SKKU.** DOI 10.1116/1.4780122. **Fetched (author's copy,
open on the SKKU lab server):**
`/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/skku_acl_2013.pdf` (+`.txt`).

This is the best *experimental* description of what an ACL mask hole actually does during opening,
and it is the **mirror image of our problem** (their ACL mouth *widens* uncontrollably):
> "When ACL was etched using O2 only, a **tapered amorphous carbon hole etch profile having a wide
> top hole size and a narrow bottom hole size** was obtained. Also, near the top hole area, a
> **bow-like amorphous carbon hole profile** was also observed. … the oxygen ions **scattered by the
> SiON mask** tended to **enhance the sidewall etching near the top contact area**, resulting in a
> bow-like and tapered amorphous carbon hole etch profile…"
> "as the etch time is increased to (b) 50% overetching and (c) 100% overetching … the tapering and
> bowing … became more serious … Then, the holes finally **merged by the excessive opening of the
> hole top side**."
> "the amorphous carbon mask with a **narrow bottle-neck like profile** due to the bow-like and
> tapered sidewall etching … was etched **faster** than the mask with the improved profile … because
> of the enhanced mask layer removal by the physical bombardment at slant [surfaces]"

Quantified effect of the COS sidewall-passivation additive (50 nm ACL hole, 900 W 60 MHz, 20 mTorr
O2 200 sccm, 300 K):
- **top/bottom opening ratio improved by ~37 %** with O2 + 5% COS vs without;
- **hole distortion reduced by ~6 %**;
- downstream: **"5% improvement in the contact oxide opening and 20% improvement in the mask
  [selectivity]"**;
- ACL etch rate cost of the additive: **−8 %**.

**Direct consequence for our model:** an ACL/a-C mask that is being etched (by O, or by O-rich FC
plasma) develops a **top-widened, bowed, bottle-necked** opening, and the *faster the slant surfaces
are bombarded, the faster the mask goes*. The mouth-relevant lesson is that **the a-C mask edge is
an actively-eroding, angle-sensitive surface, not a static wall** — and it erodes in the
*widening* direction.

**(b) Yeom, Yoon, Choi, Lee, Kim, You, Lee, "Role of oxygen in amorphous carbon hard mask plasma
etching," *ACS Omega* **8**(36), 32450–32457 (2023), DOI 10.1021/acsomega.3c02438.** Open access
(PMC10500572). **Full PDF fetch failed (PMC returned a stub); read via WebFetch of the PMC HTML.**
Mask-opening numbers, verbatim from the fetched page:
- Ar-based plasma: "**mask opening width decreased by 34 nm** compared with the 4% O2 condition";
- Kr mixture at 6.5% O2: "**mask opening width was 4.5 times wider** than that of the reference sample";
- Xe-based plasma: "**mask opening width decreased by 289 nm** compared with the 4% O2 condition".
The same fetch returned an explicit negative: **the paper contains no facet-angle / taper-angle
measurements.** (`[VERIFY]` if the figures contain angles the text does not state.)

Companion result found in search but **not fetched** `[VERIFY]`: a-C:H **etch yields and thresholds**
— "The etch yield in CF4 plasmas is **3.45**, while in O2 plasmas it is **12.3** … threshold energy of
**12 eV** for a-C:H etching … in O2 plasmas, while … **156 eV** … in CF4 plasmas", attributed to
"Ion-Enhanced Etching Characteristics of sp2-Rich Hydrogenated Amorphous Carbons in CF4 Plasmas and
O2 Plasmas" (*Materials*, PMC8198839). **If confirmed, this is a first-principles anchor for
whether our a-C mask should facet at all under our conditions** — a 156 eV threshold in FC plasma
would make the a-C mask nearly inert to FC-ion sputtering while remaining wide-open to O.

**(c) Kwon *et al.*, *Nanomaterials* **14**(2), 209 (2024)** — "Necking Reduction at Low Temperature
in Aspect Ratio Etching of SiO2 at CF4/H2/Ar Plasma," Kwangwoon University, DOI 10.3390/nano14020209.
**Open access; fetched:** `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/nanomaterials2024_necking.pdf`
(+ `.txt`, and a layout-preserved `nanomaterials2024_necking_layout.txt`).
This paper gives the community's **metric definitions**, which we should adopt verbatim:
- **Top CD** = width at the **ACL-mask / SiO2 boundary** after etch;
- **Necking CD** = the **minimum** width after etch; **Bowing CD** = the **maximum** width;
- **Necking ratio = (Top CD − Necking CD) / Top CD**; **Bowing ratio = (Bowing CD − Top CD)/Top CD**;
- **"the effective mask thickness after etching is the mask thickness *without a facet*"** — i.e.
  the community's standard practice is to **subtract the facet** when reporting remaining mask.
- Trend, verbatim: "**The necking CD tended to widen at lower temperatures. However, the necking
  ratio decreased.**" and "The results suggested that the **neutral species reaching the etch front
  of SiO2 had a low sticking coefficient**" *[abstract]*.
(200 nm trench, CF4/H2/Ar, 26 °C → −63 °C; 13.56 MHz + 2 MHz.)

**Direct action for us:** our "opening (nm)" metric should be reported as **Top CD, Necking CD, and
Necking ratio** on this definition. Krüger's wm is a *minimum aperture*; the community's **Necking
CD** is the same object, and the community's **Top CD is a different, larger number**. If any part
of our 22.5-vs-45 gap is definitional, this is where it lives.

### 2.5 Lam Research process reality check: Shen, Lill *et al.*, *Jpn. J. Appl. Phys.* 62, SI0801 (2023)

**Full citation:** M. Shen, T. Lill, J. Hoang, H. Chi, A. Routzahn, J. Church, P. Subramonium,
R. Puthenkovilakam, S. Reddy, S. Bhadauriya, S. Roberts, G. Kamarthy, "Progress report on high
aspect ratio patterning for memory devices," *Jpn. J. Appl. Phys.* **62**(SI), SI0801 (2023).
DOI 10.35848/1347-4065/accbc7. **Open access; fetched:**
`/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/lam_shen_lill_jjap2023.pdf` (+`.txt`).

Verbatim:
> "**Polymer deposition on the sidewalls of the mask, so-called necking, is considered the main root
> cause of sidewall roughness.** … For the conventional etch process, **the polymer necking
> morphology can be irregular**, resulting in uneven species fluxes to the hole that can lead to
> sidewall roughness. This mechanism is suppressed for lean, low-temperature processes as **the mask
> morphology and necking are more consistent from hole to hole.**"

> "One of the main advantages of the semi-empirical feature-scale model is to be able to capture the
> etch profile **as a function of etch time** … In this example, **the time evolution of the hard mask
> shape and the profile (top CD, bow CD, and taper) can be tracked. Parameters can be adjusted to
> match the simulated profile to experimental results at fixed time intervals.**"

> "The aspect ratio is defined as the total height of the hole, including the mask, **to the
> narrowest mask open CD**." … "**Since the carbon hard mask is mainly eroded from the top, the hard
> mask etch rate stays relatively constant.**" … "The flat etch rate at lower aspect ratios can be
> partially attributed to very low IADs (**<0.5° 1σ**)." … "Normalized selectivity is reduced by
> roughly **50%**, going from an aspect ratio of **40 to 140**."
> Reported profile-control gains (BB bias = bow CD − bottom CD): bottom-tier SiN optimization
> **−53%**, top-tier **−18%**, carbon liner **−30%**.

**Two decisive statements for our problem:**
1. **"the carbon hard mask is mainly eroded from the top"** — i.e. in a real 3D-NAND ACL process,
   the a-C mask loses *thickness*, not *corner*, and the industrial model treats it that way.
   That is evidence **against** "mask corner faceting" as our missing mouth mechanism, and
   evidence **for** treating a-C top erosion as a thickness loss with the opening preserved.
2. **The industrial standard practice is exactly the fit we are trying to do:** track top CD /
   bow CD / taper vs time and "adjust parameters to match the simulated profile to experimental
   results at fixed time intervals." **They do not claim it falls out of first principles.**
3. **IAD 1σ < 0.5°** is the operative number in a real HAR tool. Our IADF must be at least this
   narrow or the lip receives spuriously much off-normal flux.

### 2.6 Other experimental anchors (abstract-level, full texts not obtained)

- **Ohiwa, Kojima, Sekine, Sakai, Yonemoto, Watanabe, *Jpn. J. Appl. Phys.* 37, 5060 (1998)**
  (Toshiba), DOI 10.1143/jjap.37.5060. *[abstract]*: "**the redeposition of sputtered species from
  the fluorocarbon polymer on the hole sidewall induces the etch stop** at the bottom of the
  high-aspect hole… the etch stop in a high-aspect-ratio hole is determined by **the balance
  between the effects of high-energy-species bombardment and etch inhibition of carbon species**."
  → the canonical statement that **redeposition of FC polymer is a real transport channel** inside
  the feature, and that the clog is a *balance*, not a threshold.
- **Izawa, Negishi, Yokogawa, Momonoi, *Jpn. J. Appl. Phys.* 46, 7870 (2007)** (Hitachi),
  DOI 10.1143/jjap.46.7870. *[abstract]*, and these are **directly liftable constants**:
  "Sticking coefficients of radicals on the sidewall have been estimated by comparing the observed
  deposition profile with the calculated one. It was found that the coefficients of **C rich
  radicals** and **CFx radicals** were **0.5** and **0.004**, respectively, and that **F radical reaction
  probability to the fluorocarbon polymer is 0.07**. These coefficient values were deduced that the
  **excessive flux of O and F onto the sidewall of a hole causes bowing** during HARC etching. It was
  also indicated that **the bowing can be suppressed by reducing of the flux of oxygen.**"
  → **This is a competing, purely chemical bowing mechanism (O/F flux to the sidewall), with no
  facet in it.** And the CFx sticking coefficient **0.004** is *two orders of magnitude* below the
  0.1–0.0842 class of numbers our mixed layer uses for on-polymer/on-mask deposition. `[VERIFY]`
  the definition mismatch (Izawa's is a sidewall *reactive* sticking on an already-polymerized
  surface) before treating this as a contradiction — but if it survives scrutiny it is a **direct
  candidate explanation for our over-deposition at the lip.**
- **Koehler & Fischer, "Challenges of High Aspect Ratio Oxide Etching," *ECS Transactions* 13,
  47–54 (2008)**, DOI 10.1149/1.3035366 (German DRAM fab lineage). *[abstract]*: "**Profile bowing,
  potato-shaped cross-sections, and increase of length to width ratio were identified to be the
  limiting factors** in patterning of high aspect ratio holes in SiO2." — note the *absence* of
  neck/mouth closure from the limiting-factor list at AR 20:1.
- **Yen, Chang, Chiu, *Microelectronic Engineering* 82, 129–135 (2005)** (Feng Chia Univ. +
  ProMOS), DOI 10.1016/j.mee.2005.07.001. Search-record claim `[VERIFY]`: "etched holes using
  poly-HM masks exhibit **necking and bowing at the middle** of the full depth, whereas those using
  **PR masks have necking and bowing at the top** of the holes." — i.e. **the neck position is set by
  the mask material**, another argument that the mouth is mask-chemistry-coupled.

---

## 3. MECHANISM LITERATURE — who claims what, with what evidence

### 3.1 (b) Polymer deposition–erosion equilibrium at the lip — THE MAJORITY POSITION

**Primary claim: Kim, Hudson, Cooperberg, Edelberg, Srinivasan (Lam Research), *Thin Solid Films*
515(11), 4874–4878 (2007), DOI 10.1016/j.tsf.2006.10.023.** Full text **not obtained** (Elsevier
403). The claim is preserved in two independent, *fetched* secondary sources:

Huang, Huard, Shim, Nam, Song, Lu, Kushner, *J. Vac. Sci. Technol. A* **37**, 031304 (2019),
DOI 10.1116/1.5090606 — read from
`/private/tmp/claude-501/.../scratchpad/JVSTA_37_031304_2019.txt`, verbatim, their ref. 33 = Kim
*et al.*:
> "**A semiempirical profile simulator was employed to investigate the necking and bowing of etching
> of HAR features.³³ The necking resulted from a balance between polymer removal and deposition
> processes, while the bowing was caused by surface scattering of ions from secondary facets.**
> Nonuniform necking was found to cause an imbalance in the ion flux to the bottom of the feature,
> resulting in twisting irrespective of charging.⁶"

And the search-record summary of the Kim/Hudson abstract, consistent across two independent
retrievals `[abstract-level]`:
> "As neutral depositor flux was increased, the resulting profile showed a **monotonic increase in
> necking**; in contrast, the extent of **bowing showed a maximum**, such that minimal bowing was
> obtained at low and high depositor fluxes. **The primary facet angle of photoresist is determined
> by the angular dependence of the erosion rate but has little effect on oxide profile.**"

**This is the most directly relevant single sentence in the entire survey for our task:** the
group that owns the industrial HARC simulator states that (i) necking is monotone in depositor
flux, and (ii) the **primary mask facet angle barely matters for the oxide profile**, while (iii)
the facet's role is via the **secondary** facet reflecting ions into the bow region.

**Supporting/independent:** You *et al.* 2023 (§2.1, measured NDR/NER), Shen/Lill 2023 (§2.5,
"polymer deposition on the sidewalls of the mask, so-called necking"), Ohiwa 1998 (§2.6, balance).

### 3.2 (a) Mask facet funnelling — TRUE FOR BOWING, WEAK FOR THE MOUTH

- **Mahorowala & Sawin, *J. Vac. Sci. Technol. B* 20(3), 1077–1083 (2002)**, DOI 10.1116/1.1481868
  — the canonical facet-formation paper (poly-Si/Cl2-HBr, but the facet physics is material-generic).
  Full text **not obtained** (AIP paywall). *[abstract, verbatim]*:
  > "The simulations suggested that **the top facet angle was controlled by the surface composition
  > at the top of the photoresist lines and the angular dependence for etching of the deposited
  > material; the facet being less steep when there was more deposition of Si-based byproducts.**
  > **The lower facet angle and the polysilicon sidewall profile were governed by the feature aspect
  > ratio, the sticking probabilities, and fluxes of the depositing material and the depositing
  > material etching angular dependence.** Feature bottom microtrenching was strongly linked to
  > sidewall curvature, i.e., bowing… **Scattering of the ions from the curved sidewalls exhibited a
  > focusing effect** on the directional ions concentrating them at a point near the sidewall…"
  → **Facet angle is itself an OUTPUT of the deposition-vs-angular-etch balance**, not an input.
  Same physics as §3.1, applied at the mask corner. (Companion experimental paper: Mahorowala,
  Sawin, Jones, Labun, *JVST B* **20**, 1055 (2002); simulator: *JVST B* **20**, 1064 (2002);
  charging: *JVST B* **20**, 1084 (2002). None obtained.)
- **Zhang, Rauf, Sparks (Motorola), *IEEE Trans. Plasma Sci.* 30(1), 114–115 (2002)**,
  DOI 10.1109/tps.2002.1003950. *[abstract, verbatim]*: "The simulation results show **faceted
  profile evolution for the PR due to preferential ion sputtering at the incident angle to the
  facet**. The faceting also occurs on small defective surface pits, leading to expansion of the
  defect size. **Comparison between simulated and experimental profiles shows good agreement.**"
  → mask faceting *is* reproducible from an angle+energy sputter law + ion activation + F etching.
- **Miyake 2009 (§2.3):** taper angle → scattered-ion concentration in the **upper part of the
  sidewall** → bowing. (Not the mouth.)
- **Huang/Kushner 2019, the controlled facet experiment** (verbatim from the fetched text):
  > "the erosion of the PR results in **increasing the area of the facet at the top of PR**, and
  > scattering of ions at the facet produces hot neutrals into the feature with broad angular
  > distributions. **For reflections from the facets to be a direct source of bowing in SiO2, there
  > should be a line of sight from the facet to the top of SiO2, which is not the case here.**
  > Reflections from the facets here **broaden the angular distribution for subsequent sidewall
  > collisions**."
  and the a0 sweep:
  > "the probabilities of physical sputtering of PR by energetic ions … and thermal etching of PR by
  > O atoms, a0, were both decreased to **50%, 25%, and 0%** … **The selectivity (as indicated by the
  > height of the PR) improves as the PR etch probability decreases while the bowing distortion
  > decreases.** Note that in the case of a0 = 0.5 and 0.25, **there is significant bowing while there
  > is no direct path for an ion reflecting from the PR facet to reach the location of bowing.** The
  > bowing for these conditions results from energetic particles having **multiple reflections**."
  Also, the **counter-intuitive opener**:
  > "There is a **counterintuitive positive contribution to the mask erosion**. With thick masks
  > having a finite AR before reaching the top surface of SiO2, there is already some conductance
  > limit to transport of neutral radicals and sidewall scattering of ions, both of which contribute
  > to ARDE. **As the mask erodes, the conductance limit of neutrals into the SiO2 portion of the
  > feature relaxes and the unimpeded ion flux to the etch front increases. ARDE would be even more
  > severe in the absence of PR erosion.**"
  (Mask numbers: "The PR has been eroded by approximately **450 nm** by the end of the etch when the
  feature reaches the Si stopping layer, yielding a **selectivity of SiO2 over a PR of 10.7**.")

**Net:** in the *only* controlled numerical experiment on facet-vs-profile I could read
(Huang 2019 a0 sweep), turning the mask erosion **down** *reduced* bowing and *improved*
selectivity — and the bowing that remained came from **multiple reflections**, not from a
line-of-sight facet ray. Combined with Kim/Hudson's "primary facet has little effect on the oxide
profile" and Lee 2010's "the mask slope was already formed after 4 min but bowing appeared after
9 min", the literature verdict is: **mask faceting is neither necessary nor sufficient for the
mouth**, and is at best a second-order modifier of it.

### 3.3 (c) Charging — real for facet/bow/microtrench and twist; NOT claimed to set the mouth

- **Zhai, Ge, Hu, Li, Shao, Cheng, Filipovic, Chen, *J. Appl. Phys.* 137(6), 063302 (2025),
  DOI 10.1063/5.0243470.** Full text **not obtained** (AIP 403 on both PDF and HTML).
  *[abstract, verbatim]*:
  > "we propose a novel etching model … incorporating algorithms that **simultaneously account for
  > charging effects and particle reflection mechanisms**. The model is able to **reproduce the
  > real-time profile evolution of both hardmask and substrate layers** in an advanced nanoscale
  > etching process. We calculate the electric field distribution induced by the surface charges
  > accumulated on the hardmask, which affects both the trajectory of individual incident ions and
  > the overall etching profiles. To validate our approach, we perform experiments of **Si etching in
  > Cl2 plasma** and compare the simulated profiles with **scanning electron microscope images**. The
  > model also identifies the impact of the charging effect on profile defects, such as **mask
  > faceting, substrate bowing, and microtrenching**."
  Additional search-record claim `[VERIFY against the body]`: "**The accumulation of negative charges
  on the sidewalls of the hardmask deflects more positive ions to the sides, causing larger
  reflection angles and thus more etching at upper regions of the substrate sidewall.**"
  → **charging widens the upper region**, i.e. its reported sign is *opening*, not closing. If real,
  omitting charging is a plausible (secondary) contributor to an over-narrow simulated mouth — but
  the demonstrated system is **Si/Cl2, not FC/SiO2**, and no mouth-width number is claimed.
- **Huang/Kushner 2019** (fetched, verbatim): "**Positive ions with narrow angular distributions
  typically deposit charge on the bottom of low AR features, producing a maximum in positive
  electric potential on the bottom of the feature. For high AR features, grazing incidence
  collisions of ions on sidewalls depositing charge produce electric potentials with maxima on the
  sidewalls (as opposed to the bottom) of the feature.**"
- **Miyake 2009** (§2.3) explicitly finds twisting probability tracks **necking growth rate
  "irrespective of mask electrification."**
- **Huang/Kushner via ref. 6**, verbatim from the fetched Huang text: "Nonuniform necking was found
  to cause an imbalance in the ion flux to the bottom of the feature, resulting in **twisting
  irrespective of charging**."

**Verdict:** no source in this survey claims charging sets the mouth/neck width in FC HAR oxide
etch. Two sources explicitly de-emphasize charging relative to *necking geometry*. Charging is a
**deprioritized** explanation for our specific miss.

### 3.4 (d) Redeposition — a measured NARROWING channel

- **Lee 2010 (§2.2), Faraday-cage isolated:** "The redeposition of particles sputtered from the mask
  slope on the contact-hole sidewall **resulted in sidewall necking**."
- **Ohiwa 1998 (§2.6):** FC-polymer redeposition induces the etch stop.
→ If we add redeposition without simultaneously adding grazing removal, **our mouth gets worse.**

### 3.5 Modeling papers that EXPLICITLY tune/validate the mouth/CD against experiment

Ranked by directness. (Krüger 2024 is deliberately excluded — other agent's lane.)

1. **Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (Tokyo Electron), *J. Vac. Sci. Technol. A*
   42(4), 6.0003515 (2024), DOI 10.1116/6.0003515.** Full text **not obtained** (AIP paywall).
   *[abstract, verbatim]*: "we created the models for HARC etch with a **cell-based Particle Monte
   Carlo topography simulator** by **fitting both vertical and horizontal cross-sectional profiles
   carefully to the experimental results**. Moreover, we attempted to apply a **model optimization
   algorithm**. By collaboration of human and the algorithm, modeling engineers can minimize a
   try-and-error approach… **the distortion and twisting profiles were reproduced very well**."
   → **The closest published analogue of our optimization-against-a-profile workflow.**
   `[VERIFY]` whether they report a mouth/top-CD residual. **Worth a paid/ILL retrieval.**
2. **Shen/Lill (Lam) 2023 (§2.5)** — "top CD, bow CD, and taper" tracked vs time and matched by
   parameter adjustment "at fixed time intervals".
3. **Rodrigues, Aguinsky, Lenz, Hössinger, Weinbub (TU Wien + Silvaco), *J. Comput. Electron.* 22,
   1558–1563 (2023), DOI 10.1007/s10825-023-02068-y.** **Open access; fetched:**
   `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/tuwien_rodrigues_2023_fc_silica.pdf`
   (+`.txt`). Verbatim: "we are able to **accurately reproduce the etch rates, topography, and
   critical dimensions** of the reported experiments." SiO2 via + **Ru hardmask**, CF4/C4F8 ICP,
   200 nm opening, 100 nm mask, 94 s. Their sticking assumption is explicit and relevant to us:
   "we assume that their sticking values (S_n,p,n/p = **0.1**) are low enough to enable **full reactant
   supply** to the involved CD"; ion source "a sharp **von Mises** source angular distribution with a
   shape parameter of **250**"; J_i = **1.4 × 10¹⁶ cm⁻² s⁻¹**.
   → CDs matched, but with an **inert Ru mask** and a **deliberately low sticking coefficient** —
   i.e. they sidestepped both of our problems.
4. **Manstetten/Toifl/Reiter/Filipovic, "Effect of Mask Geometry Variation on Plasma Etching
   Profiles," *Micromachines* 14(3), 665 (2023), DOI 10.3390/mi14030665.** **Open access; fetched:**
   `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/mask_geometry_micromachines_2023.pdf`
   (+`.txt`). SF6/O2 (not FC) but it is the only **systematic facet-angle parameter sweep** in the
   open literature. Verbatim:
   > "**In mask faceting, ion bombardment sputters away the top mask corner and results in an
   > effective taper angle as well. However, the angles tend to be larger** and the impact is
   > analyzed for thin masks."
   > "**Higher faceting angles resulted in a fully directional vertical etching as no ion[s] were
   > reflected from the mask towards the substrate sidewall, resulting in the same profile as
   > without faceting.**"
   > "the bowing effect … **reached its maximum value between 15° and 20°** while being less
   > pronounced for the thinnest considered mask of 0.1 µm, **since there was a smaller mask sidewall
   > area from which the ions could be reflected**."
   > "**At very small angles, which corresponded to weak faceting, the maximum width occurred just
   > under the mask** as the reflected ions impinged directly at the bottom of the geometry, leaving
   > the top of the substrate most exposed to lateral etching. A further increase in tapering
   > shifted the impact point of the reflected ions from the feature bottom towards the sidewall,
   > resulting in a **steep decrease in the vertical location of the maximum width**."
   > (Also: "The **peak depth** … occurs when the mask is **tapered at about 0.5°**… **The minimum bowing
   > occurs at the peak depth**, and it increases with an increasing taper angle.")
   → **This is the facet-angle response curve.** Two things for us: (i) a **large** facet angle
   behaves like **no facet at all** (nothing reflects into the feature); (ii) **weak faceting puts
   the maximum width immediately under the mask** — i.e. facet-driven widening acts *right at the
   mouth* only in the near-0° limit, and *moves down* as the facet grows. A quantitative sweep we
   can reproduce as a gate.

---

## 4. CROSS-CHECK — has our exact failure been seen and diagnosed before?

**Question:** does anyone report that omitting mask faceting or charging causes an **over-narrowed
simulated mouth**?

**Answer: no exact match found.** The literature contains four near-misses, all instructive:

1. **Zhang (Michigan, 2015) — closest.** A published, honest **model-vs-experiment miss on necking
   and bowing**, with the diagnosis pinned on the **ion angular distribution / IEADF**, not on
   chemistry or the mask: *"There is a difference of necking and bowing **positions** between the
   experimental measurements and the predicted simulation results … **a slight change in angular
   distribution may contribute significantly to different shape evolutions**."* The failure mode is
   **position**, not **closure**. `.../tmp/pdfs/zhang_yiting_phd_thesis.txt` §7.3.
2. **Huang/Kushner 2019 — the inverse experiment.** They *suppressed* mask erosion (a0 → 0.5, 0.25,
   0) and reported **less bowing and better selectivity**, with residual bowing coming from
   multiple reflections. So suppressing the facet did **not** produce a narrowed mouth in their
   model — it produced a **cleaner** profile. This is a mild *counter-indication* that our missing
   ingredient is faceting.
3. **Lee 2010 — the experimental timing argument.** The mask slope was fully formed at 4 min while
   bowing only appeared at 9 min, "when the necking moved to positions below the mask-SiO2
   boundary." → the profile clock is set by the **neck migration**, not by facet formation. Any
   model whose neck cannot migrate downward will mis-time and, if the neck is pinned at the mouth,
   over-close it. **This is the most plausible published analogue of our failure.**
4. **Kim/Hudson 2007 — the explicit null.** "The primary facet angle of photoresist is determined by
   the angular dependence of the erosion rate but **has little effect on oxide profile**"
   *[abstract-level]*; combined with "necking increases **monotonically** with neutral depositor
   flux". → In the industrial simulator, the mouth is a **depositor-flux dial**, and the facet is
   not the dial.

**Nobody publishes a "our mouth sealed and the experiment's didn't" result.** That is expected —
it is a negative result. The nearest structural admissions are Zhang's (position miss) and the
general industry practice (Shen/Lill) of **fitting** top CD to SEMs at fixed time intervals rather
than predicting it.

**Consequence for petch:** the literature does *not* license "add mask faceting and the mouth
opens." It licenses three specific, testable statements, all of which point at the removal law and
the metric, not at the mask geometry:
- the mouth equilibrium is a **deposition-vs-removal balance whose two halves are both
  angle-resolved**, and the *measured* removal is **above cosine out to 50–60°** (You 2023);
- the **neck migrates downward with time** and the top region is *protected* by a thick CFx film
  while the action moves below it (Lee 2010) — our metric must allow the minimum aperture to move;
- **IAD 1σ < 0.5°** in production HAR tools (Shen/Lill 2023), and a **6× change in IAD width moves
  the mid-depth width by only −18%/+30%** (Huard 2018) — so the source angular distribution cannot
  by itself explain a factor-2 mouth error.

---

## 5. WHAT HOLDS THE MOUTH — final ranked answer with the receipt each implies

| Rank | Claim (with owner) | Evidence class | Receipt to run in petch |
|---|---|---|---|
| **1** | **Grazing-angle removal at the lip is ABOVE cosine out to 50–60°** — "the NERs are above the cosine curve until the ion-incident angle reaches 50–60 degrees. This indicates that physical sputtering plays an important role" (You *et al.*, *Coatings* 13, 1452, Fig. 6b). | **Measured**, Faraday cage, on SiO2 with an ACL mask, in an FC/O2/Ar plasma. | Plot petch's *effective* lip removal vs incidence angle (deposition-normalized) against a cosine reference. If it is **below** cosine at 30–70°, that is the bug — and it is exactly what `split_grazing_ion_reflection` (subtractive reflection) does. |
| **2** | **Narrowing vs bowing is set by how fast NDR(θ) falls relative to NER(θ)** — the *same* chemistry that narrows the hole is the one whose deposition stays flat with angle (You *et al.*). | **Measured**, 5-point gas-ratio sweep, with a monotone profile response (bow → anisotropic → narrow → etch stop). | Extract petch's NDR(θ)/NER(θ) at the lip and compare shape (not magnitude) to the published families. **Digitizable held-out target.** |
| **3** | **Necking is monotone in neutral depositor flux; bowing is non-monotone (has a maximum)** (Kim/Hudson, Lam, TSF 515, 4874). | Industrial simulator, abstract-level. | Sweep our depositor flux: mouth must be monotone, bow must show an interior maximum. **A shape gate that needs no absolute calibration.** |
| **4** | **The neck MIGRATES DOWNWARD; the mouth region is protected by a thick CFx film while lateral etching proceeds below it** (Lee 2010, C4F6/CH2F2/O2/Ar, time-resolved SEM). | **Measured**, time series + Faraday cage. | Report **Top CD**, **Necking CD**, **neck depth z_neck(t)**, and **necking ratio** (Kwon 2024 definitions), not a single "opening". If our z_neck is pinned at 0 while the experiment's moves down, part of the gap is definitional/kinematic. |
| **5** | **Mask-slope redeposition is a real necking channel** (Lee 2010) and **FC-polymer redeposition drives etch stop** (Ohiwa 1998). | **Measured** (Faraday cage) / measured. | Do **not** add redeposition alone. If added, it must be paired with #1 or the mouth closes faster. |
| **6** | **Mask erosion RELAXES the neutral conductance limit and helps** — "ARDE would be even more severe in the absence of PR erosion" (Huang/Kushner 2019). | Controlled numerical experiment, fetched verbatim. | Our 133 nm mask loss should be improving lip neutral supply. If lip O-flux is still starved (July doc receipt #2), the deficit is in **transport**, not mask height. |
| **7** | **Facet angle is an OUTPUT of the local deposition/angular-etch balance, not an input** (Mahorowala & Sawin 2002); **primary facet has little effect on the oxide profile** (Kim/Hudson 2007); a **large** facet angle is equivalent to **no** facet (Micromachines 2023). | Simulation + abstract-level + open-access sweep. | Deprioritize "add mask faceting" as the mouth fix. If added, expect it to move the **bow**, and expect weak faceting to put max width *just under the mask* and strong faceting to move it down. |
| **8** | **Charging sign at the mouth is reported as OPENING** ("negative charges on the hardmask sidewalls deflect more positive ions to the sides… more etching at upper regions") (Zhai 2025) — but demonstrated in **Si/Cl2**, and two sources find necking/twisting behaviour "**irrespective of mask electrification / charging**" (Miyake 2009; Huang 2019). | Abstract-level; wrong material system. | Lowest priority for the mouth. |

**One-line synthesis:** *the mouth is held open by grazing-angle physical removal that is stronger
than cosine, balanced against an angle-decaying polymer deposition; the mask facet mostly moves
the bow, redeposition only narrows, and charging is not claimed to set the mouth at all.* Our
over-narrowing is therefore most consistent with **an under-delivered / sub-cosine grazing removal
at the lip** — which is exactly the subtractive-reflection defect already identified in
`RESEARCH_MOUTH_EQUILIBRIUM_2026-07-27.md` §3b — plus a possible **metric mismatch** (Top CD vs
migrating Necking CD).

---

## 6. Fetch ledger

### 6.1 PDFs obtained this session (all under `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/`)

| File | What it is | Bytes |
|---|---|---|
| `wang_mingmei_phd_thesis.pdf` (+`.txt`) | **Wang (Kushner) — PR facet → bowing, cross-linking, Ar/C4F8/O2 HAR SiO2.** Top thesis find. | 7,309,869 |
| `zhang_yiting_phd_thesis.pdf` (+`.txt`) | **Zhang (Michigan 2015) — MCFPM-3d vs SEM; published necking/bowing position miss.** | 15,915,085 |
| `huard_chad_phd_thesis.pdf` (+`.txt`) | Huard (Michigan 2018) — bowing vs IAD width, quantified. | 11,611,696 |
| `song_sangheon_phd_thesis.pdf` (+`.txt`) | Song (Michigan 2014) — printed FC sputter/deposition probabilities. | 13,167,852 |
| `tian_peng_phd_thesis.pdf` (+`.txt`) | Tian (Michigan 2018) — ruled out for mouth. | 16,792,031 |
| `Qu_Chenhui_PhD_Thesis_2020.pdf` (+`.txt`) | Qu (Michigan 2020) — ruled out for mouth. | 21,910,985 |
| `Konina_Kseniia_PhD_Thesis_2024.pdf` | Konina (Michigan 2024) — atmospheric-pressure; ruled out. | 19,655,949 |
| `Lanham_Steven_PhD_Thesis_2022.pdf` (+`.txt`) | Lanham (Michigan 2022) — ruled out. | 14,378,751 |
| `logue_michael_phd_thesis.pdf` (+`.txt`) | Logue (Michigan) — ruled out. | 9,601,217 |
| `coatings2023_bowing_narrowing.pdf` (+`.txt`) | **You *et al.* 2023 — measured NDR(θ)/NER(θ), CD-vs-depth, narrowing→etch-stop sweep.** Best quantitative anchor. | 590,013 |
| `nanomaterials2024_necking.pdf` (+`.txt`, `_layout.txt`) | Kwon *et al.* 2024 — Top/Necking/Bowing CD definitions, necking ratio, facet-corrected mask thickness. | 4,430,658 |
| `lam_shen_lill_jjap2023.pdf` (+`.txt`) | Shen, Lill *et al.* 2023 — Lam HAR progress report; necking = mask sidewall polymer; top-CD-vs-time fitting; IAD <0.5°. | 1,212,265 |
| `skku_acl_2013.pdf` (+`.txt`) | Kim *et al.* 2013 (Samsung + SKKU) — ACL mask hole tapering/bowing/top-opening, COS additive, +37 % top/bottom ratio. | 1,658,960 |
| `mask_geometry_micromachines_2023.pdf` (+`.txt`) | Micromachines 14, 665 — systematic mask taper/facet-angle sweep (ViennaPS). | 3,340,692 |
| `tuwien_rodrigues_2023_fc_silica.pdf` (+`.txt`) | *J. Comput. Electron.* 22, 1558 — FC/SiO2 + Ru mask, CDs matched to experiment. | 1,116,355 |

Kushner thesis directory listing (57 entries) parsed at
`/private/tmp/claude-501/-Users-stanislavdelaurentiis-chip-etch/352207a8-b462-4fe3-a650-c9b37f19523b/scratchpad/theses_index.html`.

### 6.2 Wanted but NOT obtained (paywall/403) — abstract-level only

| Source | Barrier | Why it matters |
|---|---|---|
| **Kim, Hudson, Cooperberg, Edelberg, Srinivasan, *Thin Solid Films* 515, 4874 (2007)**, DOI 10.1016/j.tsf.2006.10.023 | Elsevier 403 | **The** mouth-mechanism paper. Its depositor-flux sweep is the exact experiment we want to reproduce. **Highest-value retrieval.** |
| **Lee *et al.*, *J. Electrochem. Soc.* 157, D142 (2010)**, DOI 10.1149/1.3276511 | IOP/ECS login | The neck-migration timeline and Faraday-cage redeposition numbers. |
| **Miyake *et al.*, *Jpn. J. Appl. Phys.* 48, 08HE01 (2009)**, DOI 10.1143/jjap.48.08he01 | paywall | **The mask-taper-angle DOE — the facet-angle numbers we most want.** |
| **Nishizuka *et al.* (TEL), *JVST A* 42, 6.0003515 (2024)**, DOI 10.1116/6.0003515 | AIP paywall | Closest published analogue of our fit-to-profile workflow; may report a mouth residual. |
| **Zhai *et al.*, *J. Appl. Phys.* 137, 063302 (2025)**, DOI 10.1063/5.0243470 | AIP 403 (despite OA flag) | Charging + reflection → mask faceting/bowing/microtrench, validated vs SEM. |
| **Mahorowala & Sawin, *JVST B* 20, 1055 / 1064 / 1077 / 1084 (2002)** | AIP paywall | Canonical measured facet angles + the four-part experiment/simulator/charging set. |
| **Izawa *et al.*, *JJAP* 46, 7870 (2007)**, DOI 10.1143/jjap.46.7870 | paywall | Sidewall sticking coefficients (C-rich 0.5, CFx 0.004, F→polymer 0.07) — potentially decisive for our lip over-deposition. |
| *Materials* / PMC8198839 a-C:H yields (CF4 3.45, O2 12.3; thresholds 156 eV / 12 eV) | PMC PDF stub | First-principles check on whether our a-C mask should facet under FC ions at all. |

### 6.3 Explicitly checked and found NOT to contain mouth/facet content (negative results, so nobody re-searches them)
Konina 2024, Lanham 2022, Logue, Tian 2018, Qu 2020 (all Kushner theses, fetched and grepped);
Manstetten TU Wien dissertation (flux-calculation methods, no mask evolution); Bogaerts group
(no feature-scale HAR mask thesis surfaced); Bochum/RUB (reactor-scale only); Kokkoris/NTUA
(ARDE/clog, mask not evolved, charging only as "increased ion angular spread").

---

## 7. Next actions this survey licenses (in order)

1. **Adopt the community metric set** — Top CD, Necking CD, neck depth z_neck(t), necking ratio
   (Kwon 2024 definitions) — and re-extract ml9a with them. Free; may recover part of the gap.
   Check whether our neck can migrate below the mask–oxide boundary at all.
2. **Plot petch's lip removal vs incidence angle against a cosine reference** and compare to
   You 2023 Fig. 6b (NER above cosine to 50–60°). This is the single sharpest test of the
   subtractive-reflection hypothesis, and it needs no new physics.
3. **Run the depositor-flux monotonicity gate** (Kim/Hudson): mouth monotone-narrowing in depositor
   flux, bowing non-monotone with an interior maximum. Shape-only, calibration-free.
4. **Retrieve Kim/Hudson TSF 2007 and Miyake JJAP 2009** (ILL / library). These two carry the only
   published facet-angle DOE and the only published depositor-flux→necking curve.
5. **Do NOT** add mask faceting, redeposition, or charging as the mouth fix on current evidence.
   Faceting moves the bow; redeposition narrows; charging is unclaimed for the mouth.
6. Optional constant-check: Izawa's CFx sidewall sticking **0.004** vs our ~0.1/0.0842 class.
   If the definitions are commensurate, our lip is over-depositing by orders of magnitude and the
   diagnosis changes completely. `[VERIFY]` first.
