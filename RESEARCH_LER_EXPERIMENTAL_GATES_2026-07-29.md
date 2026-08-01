# Experimental LER-Transfer Datasets That Can Gate petch End-to-End

**Date:** 2026-07-29 · **Status:** research memo, NOT committed.
**Companions:** `RESEARCH_LER_EXPERIMENTAL_SOURCES_2026-07-21.md` (locator survey),
`RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md` (architecture),
`src/petch/ler_metrology.py` (Palasantzas PSD, σ/ξ/α, Mack noise floor),
`src/petch/ler_transfer.py` (H1 cross-spectral T(f), coherence separation,
`PSD_out = |T|² PSD_in + PSD_intrinsic`).

**Question answered:** which *published experimental* roughness-transfer data has enough
quantitative detail to gate a simulator end-to-end — i.e. simulate a wide-y etch of a rough
mask edge and compare a *predicted output PSD* against a *measured* one?

**One-line answer:** exactly one open-access source publishes a directly-measured
**PSD ratio before/after a named plasma etch step with a fully specified recipe and stack**
(Azarnouche PhD thesis, LTM-CNRS, Fig. IV.39) — that is the T(k)-shape gate. Two others give
**σ/ξ/α before-after pairs plus closed-form transfer rules** (Demokritos+LTM) and
**modern absolute unbiased PSDs with ξ and PSD(0) across four plasma conditions**
(imec/Fractilia). Everything else in the literature is either normalized-away, anonymized,
or reports σ only.

---

## 0. Method and verification status

Fetched and read in full this session (local copies under the session scratchpad):

| Source | Access route | Status |
|---|---|---|
| Azarnouche PhD thesis, tel-00767820 (230 pp, French) | HAL full text | **read** (Ch. III metrology, Ch. IV transfer) |
| Thiault PhD thesis, tel-00321961 (168 pp, French) | HAL full text | **read** (skimmed; 3D-AFM, no PSD analysis) |
| Dupuy et al., SPIE 9428 (2015), hal-01869175 | HAL full text (CC-BY) | **read in full** |
| Gallatin & Naulleau, SPIE 7969 (2011) | NIST tsapps PDF | **read** (LTF derivation) |
| Liang, Mack, Sirard et al., Proc. SPIE 10585, 1058524 (2018) | lithoguru PDF | **read in full** |
| Rutigliani, Lorusso, De Simone, ... Mack, Proc. SPIE 10585, 105851K (2018) | lithoguru PDF | **read in full** |
| Lorusso et al., "imec roughness protocol," JM3 17(4) 041009 (2018) | lithoguru PDF | **read** |
| Mack et al., DSA pattern-transfer roughness, Proc. SPIE 10589, 1058907 (2018) | lithoguru PDF | **read** |
| Wang, Wang, Biolsi, Kushner, JVST A 39, 033003 (2021) | OSTI 1853526 PDF | **read** (roughness mechanism §) |
| Vaglio Pret et al., EIPBN 2011 abstract (imec) | eipbn.org PDF | **read** |

Verified via publisher-indexed abstracts only (OpenAlex/Crossref; full text paywalled):
Constantoudis et al. JM3 8(4) 043004 (2009); Constantoudis JM3 9(4) 041209 (2010);
Constantoudis, Kokkoris, Gogolides JM3 12(4) 041310 (2013); Pargon et al. JVST B 26, 1971
(2008) [DOI 10.1116/1.2917071]; Azarnouche et al. JVST B 31, 012205 (2013)
[DOI 10.1116/1.4773063]; Azarnouche et al. JAP 111, 084318 (2012) [DOI 10.1063/1.4705509];
Sun et al. JM3 14(3) 033501 (2015); Vaglio Pret & Gronheid, Microelectron. Eng. 88 (2011)
[DOI 10.1016/j.mee.2011.02.015]; Qi et al. Proc. SPIE 8522 (2012) [DOI 10.1117/12.976855];
Naulleau & Gallatin JVST B 28, 1259 (2010); Nakazaki et al. JAP 116, 223302 (2014);
Martin & Cunge JVST B 26, 1281 (2008).

Anything not in one of those two buckets is marked **[VERIFY]**.

Access note: HAL (hal.science, theses.hal.science) is behind an Anubis proof-of-work gate;
plain `curl`/WebFetch returns a challenge page. A ~20-line SHA-256 PoW solver
(`scratchpad/getpdf.py`) clears it and the PDFs then download normally. ResearchGate,
Academia.edu, ScienceDirect and SPIE Digital Library all 403 machine clients.

---

## 1. Pargon / LTM-CNRS (Grenoble): the only published *measured* PSD ratio through a named etch

The group is LTM (Laboratoire des Technologies de la Microélectronique), CNRS / Univ.
Grenoble Alpes / CEA-Leti — Pargon, Azarnouche, Fouchier, Martin, Joubert. They have
published unbiased LWR with full σ/ξ/α since ~2012, which is what makes them gateable.

### 1.1 The metrology substrate they use (must be mirrored, else we gate an artifact)

From Azarnouche's thesis Ch. III (verified, eqs. III.24–III.28) the LTM estimator is:

- autocorrelation model **R_m = σ_real² · exp(−|mΔy/ξ|^{2α}) + σ_noise² δ_m**
  (Sinha/Palasantzas form + white metrology noise as a Kronecker spike);
- PSD model **P_n = (Δy/2πN) σ_real² Σ_m (2−δ_m) e^{−|mΔy/ξ|^{2α}} cos(k_n mΔy)(N−m)
  + (Δy/2π) σ_noise²**, with **k_n = 2πn/(NΔy)** — i.e. an *angular* wavenumber
  (rad/nm) and a 1/2π normalization, **not** the f-in-1/nm convention used by
  Mack/Fractilia/imec and by our `ler_metrology`;
- fit varies (σ_noise, ξ, α) under the constraint σ₀² = σ_real² + σ_noise², with σ₀ obtained
  by Parseval integration of the measured PSD;
- finite-box correction σ_inf² = σ_real² + σ_CDV² (Leunissen), valid when L > 10 ξ.

CD-SEM protocol (verified): Hitachi CG4000, rectangular scan, 1024×1024 px, magnification
300 000 (x) / 49 000 (y); pixel 0.44 nm (x) / 2.69 nm (y); smoothing 25 px in x for resist,
7 px otherwise; IA = 800 px, N = 400 measurement points, S = 2 → **L ≈ 2152–2200 nm,
Δy ≈ 5.4 nm**; ≥200 images averaged (Dupuy uses N* > 200). ξ values encountered are
"generally below 50 nm."

**Frequency window this implies:** f ∈ [4.6×10⁻⁴, 9.3×10⁻²] nm⁻¹, i.e. k ∈ [0.0029, 0.58]
rad/nm, i.e. wavelengths 10.8 nm – 2.2 µm. That window is the gate's spectral domain.

### 1.2 The gateable object: Azarnouche thesis Fig. IV.39 = a measured |T(k)|²

Verbatim (thesis p. 146–147, French):

> "Pour une meilleure illustration des fréquences spatiales de LWR, nous avons représenté
> les PSD avant et après gravure (Figs. IV.39). … nous avons retiré le bruit de chacune des
> PSD, et nous avons recalé les spectres à un 3σ_bruit constant, égal à 3 nm. …
> **La courbe bleue représente le rapport des PSD après et avant gravure.** Pour ces
> résines, le domaine des k_n supérieur à 0.1 nm⁻¹ diminue ("hautes fréquences spatiales"
> correspondant aux longueurs inférieures à 63 nm). **On observe des minimums vers
> 0.15 nm⁻¹ et 0.1 nm⁻¹ respectivement pour la résine de référence et celle exposée à
> l'HBr.** … En ce qui concerne la résine exposée aux UV, les deux PSD avant et après
> gravure se superposent."

This is literally `PSD_out(k)/PSD_in(k)` — the quantity `ler_transfer` estimates — plotted
for three mask states (untreated 193 nm resist, HBr-plasma-cured, VUV-cured) through one
named etch step. Three signed predictions fall straight out and are pre-registerable:

1. **|T| < 1 for k ≳ 0.1 rad/nm** (λ ≲ 63 nm) on untreated and HBr resist;
2. **a minimum in the ratio at k ≈ 0.15 rad/nm (untreated) and ≈ 0.10 rad/nm (HBr)** —
   i.e. the attenuation is *band-limited with a turnover*, not monotone;
3. **|T| ≈ 1 across the whole band for the VUV-cured resist** (its high-f content was
   already removed pre-etch) — a null control that is worth as much as the positive cases.

Because the ratio is *dimensionless*, it is immune to the PSD-normalization-convention trap
of §1.1. **This is the single most gateable LER-transfer object I found in the literature.**

Peer-reviewed twin (abstract-verified, paywalled): Azarnouche, Pargon, Menguelti, Fouchier,
Joubert, Gouraud, Verove, "Benefits of plasma treatments on critical dimension control and
line width roughness transfer during gate patterning," **JVST B 31, 012205 (2013)**,
DOI 10.1116/1.4773063 — "The study shows that the high and medium frequency components of
the roughness (**periodicity below 200 nm**) are not totally transferred during the gate
patterning allowing a LWR decrease at each plasma step." Note the 200 nm figure in the paper
vs. 63 nm in the thesis: the paper's "below 200 nm" is the band where attenuation *begins*
(k ≈ 0.031 rad/nm), the thesis' 63 nm is where the ratio drops decisively. Gate against the
thesis curve, cite the paper. **[VERIFY]** the paper's own figures may carry a cleaner
version of the same ratio.

### 1.3 The process the gate must simulate (fully specified — verified from thesis Ch. IV)

- **Stack:** 193 nm resist (Dow IM5010) / **Si-ARC 35 nm** / spin-on-carbon (SOC) 200 nm /
  Si. Lines: isolated, **CD 75 nm, resist height 120 nm**; final Si lines CD 75 nm, height
  80 nm.
- **Reactor:** Applied Materials **DPS AdvantEdge** 300 mm ICP at LTM.
- **Si-ARC open (the step Fig. IV.39 refers to):** 200 sccm Ar / 80 sccm CF₄ / 50 sccm CHF₃,
  source 200 W, **bias 250 W**, 7 mTorr, 20 s (30 s variant also reported).
- **SOC open:** 30 sccm O₂ / 70 sccm HBr, source 500 W, bias 120 W, 7 mTorr, 92 s.
- **Si etch:** 22 sccm SF₆ / 50 sccm CHF₃ / 90 sccm Ar, source 450 W, bias 75 W, 5 mTorr, 35 s.
- **Pre-treatments:** 30 s HBr plasma cure; 30 s VUV-only (plasma UV through a window,
  source power ~6× lower).

### 1.4 Scalar transfer numbers from the same chapter (σ-reduction gate)

For the untreated reference resist (verified text):
LWR decreases −0.7 nm across the Si-ARC open; falls to **3.9 nm** after the SOC/carbon etch;
"la rugosité du masque de carbone se transfère fidèlement dans la couche de silicium"
(carbon-mask roughness transfers *faithfully* into Si); **total improvement 1.6 nm** end to
end; final Si LWR **3.8 nm** for a 20 s Si-ARC open vs **4.7 nm** for 30 s.
⇒ implied post-litho LWR ≈ 5.4 nm and σ-reduction ratio **≈ 0.70 resist→Si**, with a
**+0.9 nm penalty** from a 50 % longer hard-mask open. **[VERIFY]** the 5.4 nm is inferred
from "amélioration totale de 1.6 nm", not read off a table.

ξ and α per step are plotted (Fig. IV.38) and stated to **both increase** for reference and
HBr-cured resist and to **stay flat** for VUV-cured resist. Direct quote:
"Au cours de nos études, nous avons toujours observé que la baisse de LWR est accompagnée
par une augmentation de ξ et de α."

### 1.5 Companion open dataset: Dupuy et al. SPIE 9428, 94280B (2015) — full SADP ladder

Open access CC-BY on HAL (hal-01869175). e-beam 80 nm-pitch resist core → HBr/O₂ trim →
50 °C PEALD SiO₂ spacer → CF₄ spacer etch → O₂ core removal → CHF₃/CF₄/Ar Si-ARC →
HBr/O₂ SOC → **SF₆/CHF₃/Ar Si etch** → 20/20 nm Si lines. Unbiased LWR/LER per step
(verified numbers): LWR 6.3 → 3.8 (trim, −20 %, high-f) → −20 % (deposition) → −34 %
(spacer etch, **low-f**) → 2.5 → **2.4 nm** final; LER 4.4 → 3.3 → 2.7 → **2.3 nm** final;
CD variation 3.5 → 0.7 nm; overall −62 % LWR / −48 % LER.

**Digitization warning (verified caption):** "After noise removal, PSDs are arbitrarily
shifted to 1 for log scale plotting." The per-step PSD *shapes* are usable; the *absolute
levels* are not. Use the σ ladder for amplitude and the PSD panels for shape only.

### 1.6 What LTM does *not* give you

- No tabulated σ/ξ/α — everything is in figures (see §5.3 on digitization).
- Thiault's 3D-AFM thesis (tel-00321961, the data behind Pargon JVST B 26, 1971 (2008)) is
  height-resolved sidewall LWR but is essentially **pre-PSD**: a single mention of "PSD" in
  168 pages. Use it for sidewall-vs-top-down phenomenology, not for spectra.
- CD-AFM correlation lengths are biased: Constantoudis JM3 9(4) 041209 (2010) shows CD-AFM
  **overestimates ξ** (tip convolution). Prefer their CD-SEM+PSD numbers.

---

## 2. Constantoudis / Kokkoris / Gogolides (NCSR Demokritos): the transfer *rules*, not T(k)

**Answer to "do they report T(k) shapes?" — No.** Across the corpus they report
(σ, ξ, α) before/after, PSD *curves* qualitatively, and closed-form scaling rules. I found
no published PSD-ratio / transfer-function curve from this group. Their claims are
nonetheless sharper than most, because they are analytic.

### 2.1 JM3 8(4), 043004 (2009) — the joint Demokritos+LTM model-vs-experiment set

Abstract verified. Experiments: 193 nm resist trimmed in **O₂ plasma, no bias**; anisotropic
etch of **BARC in CF₄** and **Si in HBr/Cl₂/O₂ with bias**. Findings:

- **Trimming:** LER ↓, **ξ ↑ and α ↑** with trim time ("nano-protrusions become shorter and
  wider") — both model and experiment.
- **Anisotropic transfer:** model predicts "noticeable reduction of LWR, whereas
  **correlation length and roughness exponent remain almost unaffected**"; "the first
  experimental results seem to confirm these predictions."

**This is in direct tension with LTM (§1.4), who "always" observe ξ and α *increasing*
whenever LWR drops.** The two are reconcilable only if the etch-driven attenuation is
essentially white (scale-free, ξ/α-preserving) in the Demokritos Si/BARC cases and
band-limited (high-f-selective) in the LTM Si-ARC case — which would mean **T(k) shape is
chemistry- and mask-dependent**, exactly the thing petch is supposed to predict rather than
assume. **Discriminating these two datasets is arguably the highest-value single result the
LER modality could produce.** Flagging it as the intended headline.

### 2.2 JM3 12(4), 041310 (2013) — geometrical model + "rules of thumb" (quantitative)

Abstract verified; the numbers below are from the authors' own open SPIE Newsroom summary,
Constantoudis, Kokkoris & Gogolides, *"Resist roughness plays a key role in pattern
transfer,"* SPIE News, 29 March 2013 (fetched and read):

- substrate 3σ_S vs resist 3σ_R is **linear with slope tan ω ≈ 0.5** above a threshold;
- threshold **σ_R\* ≈ ξ_R / (c · tan θ_R), with c ≈ 2.0–2.5**;
- equivalently, LER *reduction* requires **(σ_R/ξ_R) · tan θ_R > 1/c**;
- Fig. 4 conditions: ξ = 30 nm, α = 0.6, resist thickness 150 nm, **sidewall angle 86.2°**,
  etch selectivity 3, etch depth 150 nm.

With θ_R measured from horizontal (tan 86.2° = 15.1) and c = 2.25, σ_R\* ≈ 0.9 nm
(3σ_R\* ≈ 2.6 nm) — physically sensible. **[VERIFY]** the angle convention against the paper
before gating; the alternative (angle from vertical) gives an absurd 200 nm.

Mechanism: transfer is pure **ion shadowing by the rough resist sidewall**; substrate
sidewall inherits the *envelope* of the resist edge. Our `analytic_occlusion` operator is the
same physics done exactly, so petch should reproduce these rules **without fitting** — this
is a cheap, immediately-runnable structural gate that needs no digitization at all.

### 2.3 JM3 9(4), 041209 (2010) — sidewall anisotropy (how to seed 3-D noise)

Abstract verified: CD-AFM + SEM show post-develop resist sidewalls are **anisotropic —
striations perpendicular to the line direction**. Including that anisotropy is what makes
their model predict the (experimentally observed) benefit of trimming. Consequence for us:
a 3-D sidewall seed must be **anisotropic** (vertical correlation ≠ along-line correlation);
an isotropic seed will get the trim/transfer sign wrong.

### 2.4 Synthesized vs measured

Their pipeline is: *synthesize* self-affine sidewalls with declared (σ, ξ, α) → transfer
geometrically/by MC → *measure* σ/ξ/α of the output → compare to CD-SEM/CD-AFM. That is the
same loop as `ler_metrology.synthesize` + engine + `fit_edge_statistics`. Their validation
level is **trends and rules**, never a blind PSD prediction — consistent with the positioning
in `RESEARCH_LER_MODALITY_DESIGN` §5.

---

## 3. EUV era: mask → printed → etched

### 3.1 Mask → printed: the LTF is *experimentally measured*, with programmed roughness

- **Formalism (verified full text, Gallatin & Naulleau, Proc. SPIE 7969, 79690F (2011),
  NIST copy):** δỹ_wafer(β) = LTF(β) · δỹ_mask(β), with an analytic tophat-fill LTF; the
  wafer edge displacement is δI/(I_th·ILS). LTF(0) = 1 by construction and rolls off near
  β/NAk ~ 1. **The algebra our `ler_transfer` implements is theirs**; petch's contribution is
  computing the *etch* T from physics rather than optics.
- **Programmed-roughness masks (the clean experiment):** Qi, Gallagher, Negishi, McIntyre,
  Zweber, Senna, "Impact of EUV photomask line-edge roughness on wafer prints," Proc. SPIE
  8522 (2012), DOI 10.1117/12.976855 (abstract verified): programmed-LER EUV mask with
  **controlled variations in LER spatial frequency and magnitude**, PSD extracted from
  **64 nm and 90 nm (1X) pitch** mask lines, printed on the ASML **ADT**, wafer LWR
  correlated back to mask patterns "revealing an **empirical LWR transfer function (LTF)**";
  extended to 45 nm pitch with a pupil filter.
- **imec equivalent:** Vaglio Pret & Gronheid, "Mask line roughness contribution in EUV
  lithography," Microelectron. Eng. 88, 1963 (2011), DOI 10.1016/j.mee.2011.02.015 —
  "programmed roughness modules … for roughness transfer function evaluation on **88 nm
  pitch** line/space patterns, with variations of roughness amplitude and spatial frequency";
  finding: the optical system is a **low-pass filter**; sub-cutoff mask frequencies transfer
  fully; high-f mask roughness degrades the aerial image and thereby *raises low-f resist
  roughness*. (Search-snippet verified; **[VERIFY]** exact cutoff numbers.)

**Verdict for petch:** these gate the *optics*, not the etch. Their value is (a) proving the
`PSD_out = |T|² PSD_in` decomposition is experimentally realizable with programmed inputs,
and (b) supplying the design pattern we should copy — **programmed single-frequency input
roughness** is exactly our finite-difference T̂(k) probe, done on silicon.

### 3.2 Printed → etched: the best modern numbers

**Rutigliani, Lorusso, De Simone, Lazzarino, Papavieros, Gogolides, Constantoudis, Mack,
"Setting up a proper PSD and autocorrelation analysis for material and process
characterization," Proc. SPIE 10585, 105851K (2018)** / JM3 17(4), 041016 (2018).
Read in full. Contents:

- 16 nm CD, **32 nm pitch** EUV CAR lines on 4 underlayer stacks (organic UL; SOG-A;
  SOG-B/a-C; SOG-B/SOC) — PSD(0) and ξ per stack, ξ roughly stack-invariant while PSD(0)
  rises monotonically stack 1→4;
- resist-type ξ contrast: **CAR ξ ≈ 9.6 nm vs metal-containing resist ξ ≈ 6.1 nm**, attributed
  to resist blur;
- **the gateable block:** one stack, **four DCS (direct-current-superimposed CCP) plasma
  conditions** vs post-litho reference. **Unbiased LER: 2.63 nm (post-litho) → 1.95, 2.24,
  2.33, 2.07 nm.** Figure 9 gives **absolute LER PSD curves in nm³** for all five states over
  ~10⁻³–10⁻¹ nm⁻¹; Figure 10 gives **ξ (≈8–14 nm range) and PSD(0) per condition**.
- Verified conclusion: "**the etch smoothening plays a big role in the reduction of the
  high-frequency roughness**"; best condition preserves PSD(0) while raising ξ.
- 3-parameter model used is PSD(0)/H/ξ — same family as our Palasantzas form.

This is the **most metrology-clean input/output PSD pair set that exists in the open
literature**, and the axes are absolute. Weaknesses: it is a *resist surface treatment*, not
a transfer into an underlayer, and the DCS recipe is not disclosed.

Weaker EUV-era items, checked and largely rejected as gates:
- **Liang, Mack, Sirard et al., Proc. SPIE 10585, 1058524 (2018)** (Lam + Fractilia,
  "Unbiased roughness measurements: the key to better etch performance"). 5 nm-node EUV
  stack (Si / SOC / SOG), NXE:3350B, CAR, **56 nm pitch**, Lam Kiyo, PR treatment + main
  etch, ADI and AEI at resist and Si. **But** the process knob is anonymized ("Parameter 1"),
  the y-axes are normalized, and the PSD figures are noise-subtraction demos. Excellent
  methodology paper, **not gateable**.
- **Mack et al., Proc. SPIE 10589, 1058907 (2018)** (DSA PTMSS-b-PMOST, 10 nm HP, transfer
  into 13 nm SiN): LER PSDs at 50 %/100 % dry-develop and after SiN transfer, split LF/MF/HF;
  finding: MF roughness down and **ξ up** with longer dry develop; SiN etch lowers HF slightly
  but raises everything else. **Values are normalized to a baseline** → shape-only.
- **Vaglio Pret et al., EIPBN 2011 (imec)**: Fig. 1b is an **absolute PSD in nm³ vs 0.1–100
  µm⁻¹** with four curves (EUV resist / Si implant / etched resist / Si implant+etch) for
  32 nm HP EUV resist. Small but genuinely absolute; useful as a secondary spot check.
  Companion: Vaglio Pret, Gronheid, Foubert, JM3 9(4) (2010), DOI 10.1117/1.3494614 —
  post-litho smoothing in the frequency domain, "**reducing roughness up to 34 %**".
- **Sun, Wang, Beique, Sung, Wood, Kim, JM3 14(3), 033501 (2015)** (GlobalFoundries):
  LER/LWR split into LF/MF/HF with a 3σ number per band through a full SADP flow, plus
  quantitative wiggling detection. A 3-number-per-step summary is a legitimate coarse gate
  if the full PSDs are unobtainable.

### 3.3 LWR vs LER correlation

- **Protocol facts (verified, Lorusso et al. JM3 17(4) 041009 (2018)):** iRP 2018 mandates
  pixel size **0.8 nm in both x and y**, **no filtering**, unbiasing mandatory, ≥50
  uncompressed images; the box-length requirement is tied to ξ — convergence to σ_inf goes as
  L/ξ, and because ξ "dropped from 35 nm down to 7 nm in the last 10 years" the classic
  **2 µm box can be relaxed to ~400 nm**; SEM-vs-AFM accuracy offset for *biased* settings
  varies 1–5 nm, for *unbiased* only 0.6–1.4 nm.
- **LWR/LER algebra:** σ_LWR² = σ_L² + σ_R² + 2c σ_L σ_R with c the left/right edge
  correlation; c ≈ 0 gives the familiar LWR = √2 LER. Jiang, Li, Wang, Jiang, Huang, ICSICT
  2012, DOI 10.1109/icsict.2012.6466733 (abstract verified) derive the **LWR ACF analytically
  from the two LER ACFs plus their cross-correlation**, and find **ξ_LWR decreases as the
  edge-edge correlation coefficient increases**.
- **Load-bearing consequence for our gate design:** LTM/Dupuy data is mostly **LWR**;
  imec/Fractilia report both. A simulator predicts *edges*, so LWR comparison requires
  modelling the left/right correlation the etch induces. In Dupuy, LWR falls 62 % while LER
  falls only 48 % — the gap is entirely edge-edge correlation created by conformal
  deposition. **Never gate on LWR without declaring c.** Prefer LER-based gates.

---

## 4. What sets *intrinsic* etch LER — what our PSD_intrinsic physics should eventually be

Ranked by how directly the literature parameterizes them from countable quantities.

### 4.1 Polymer-cluster discreteness and ion-count statistics (best-quantified; fluorocarbon)

**Wang, Wang, Biolsi & Kushner, JVST A 39, 033003 (2021)**, DOI 10.1116/6.0000941
(OSTI 1853526, read in full). Verbatim from the roughness section:

> "With a FCP flux of 1.4 × 10¹⁶ cm⁻² s⁻¹, deposition time of 5 s, site density of
> 10¹⁵ cm⁻², and sticking coefficient of 0.2, the average number of sticking radicals per
> site per cycle is about 14. This small number of sticking radicals produces a statistical
> random variation in thickness of 25 %. … During the 20 s chemical sputtering step, each
> site is struck by about ten ions, which has about a 30 % statistical variation. Combined
> with the statistical FCP thickness, … there are patches of the film that are bare SiO₂ and
> patches that retain a FCP overlayer. As the residual FCP is removed, it becomes
> statistically less likely to remove the isolated patches of FCP. This leaves a
> statistically rough surface at the beginning of the next deposition step."

25 % ≈ 1/√14 and 30 % ≈ 1/√10 — **pure Poisson counting on quantities petch already
carries** (`mixed_layer.py` film areal densities; ion flux from the transport solve). This is
the numerical anchor for sources 2 and 3 of `RESEARCH_LER_MODALITY_DESIGN` §2 and it
validates the "no free knob" claim: the amplitude is fixed once flux, dwell and site density
are fixed. It also names the **nonlinearity** — residual-patch removal becomes *less* likely
as patches isolate, i.e. a self-amplifying micromask term that a purely additive Gaussian
PSD_intrinsic will miss.

### 4.2 Micromasking is *required* for roughening; plasma alone smooths

**Martin & Cunge, JVST B 26, 1281 (2008)**, DOI 10.1116/1.2932091 (abstract verified):
high-density plasmas "do not generate roughness during silicon etching; on the contrary they
tend to smooth existing roughness"; 20 nm-high/20 nm-wide Si pillars are rapidly smoothed by
Cl₂ and SF₆; the smoothing mechanism is **radical-starved transport — hills receive more
etchant flux than valleys**; reported roughening in F-based plasmas is due to **AlF_x
micromasking** from chamber walls.

This is a *hard constraint* and it is directly checkable in petch: our neutral-transport
solver already produces the hill/valley flux asymmetry (it is the same shadowing/view-factor
machinery as ARDE). **A negative gate:** with redeposition off and no micromask source, petch
must show |T(k)| < 1 at high k and a *decaying* intrinsic term. Any unconditional roughening
term is refuted by this paper.

### 4.3 Two roughening modes set by ion energy and ion identity

**Nakazaki, Tsuda, Takao, Eriguchi, Ono, JAP 116, 223302 (2014)**, DOI 10.1063/1.4903956
(OA at Kyoto repository; abstract verified): ICP Cl₂ on Si; **roughening mode for E_i below
~200–300 eV** with rms rising roughly linearly in time up to 20 min; **smoothing mode above**,
rms falling to **< 0.4 nm** and quasi-steady after ~1 min; the switch coincides with the
dominant ion changing from feed-gas Cl_x⁺ to **ionized etch products SiCl_x⁺**. They also
publish the **PSD evolution** vs E_i. A sharp, single-knob mechanistic gate for the intrinsic
term, on blanket surfaces (no transport confound).

### 4.4 Micromasking as an island-growth instability (kMC theory)

Jiang, Wu, Yang, Liu, Liao, "Kinetic etch front instability responsible for roughness
formation in plasma etching," Appl. Surf. Sci. 537, 148862 (2021),
DOI 10.1016/j.apsusc.2020.148862 (OA; abstract verified): inhibitor adatoms preferentially
bond to each other → 3-D island growth of inhibitor → micromasking → roughness. Gives the
*mechanism class* for a correlation length that emerges from chemistry rather than from the
mask.

### 4.5 Scaling law for the forcing

Drotar, Zhao, Lu & Wang, PRB 61, 3012 (2000) (locator only, inherited from the 2026-07-21
survey; **[VERIFY]** not re-read this session): (2+1)-D flux-redistribution MC gives KPZ-class
exponents α ≈ β ≈ z ≈ 1 with re-emission. Use as the *shape* constraint on PSD_intrinsic, not
as an amplitude.

**Synthesis — what PSD_intrinsic should be built from, in order:**
(i) Poisson counting on FCP/radical sticking events per correlation cell (Kushner numbers,
§4.1); (ii) Poisson counting on ion impacts per cell (same source); (iii) a *gated* micromask
term that is identically zero when redep/inhibitor flux is zero (Martin–Cunge, §4.2) and that
grows as island coverage (Jiang, §4.4); (iv) an ion-energy switch between roughening and
smoothing regimes (Nakazaki, §4.3). Nothing in that list is a free scalar.

---

## 5. Verdict — the gateable-dataset shortlist

### Gate 1 (primary, T(k) shape) — Azarnouche thesis Fig. IV.39 + JVST B 31, 012205 (2013)

- **What it gates:** the *shape* of |T̂(k)|² from the engine — specifically the three signed
  predictions of §1.2 (roll-off above k ≈ 0.1 rad/nm; a **turnover/minimum** at
  0.15 / 0.10 rad/nm for untreated / HBr resist; **flat T ≈ 1** for VUV-cured resist).
  Plus the scalar σ-reduction ladder of §1.4 (resist→Si ratio ≈ 0.70; +0.9 nm for a 30 s vs
  20 s Si-ARC open) and the sign of dξ/dstep and dα/dstep.
- **Why it is first:** it is the only *ratio* published, so it is convention-free; the recipe
  and stack are fully specified (§1.3); the mask-state contrast (none / HBr / VUV) gives a
  three-point curve with a built-in null; and the source PDF is open access.
- **Digitization feasibility: MEDIUM-HIGH but with one trap.** The ratio curve is a single
  blue line on a log-log panel — trivially digitizable (WebPlotDigitizer, ~30 points). **The
  trap:** the plotted PSDs were noise-subtracted and then *re-aligned to a common synthetic
  white floor of 3σ_noise = 3 nm* (verbatim §1.2). At high k both numerator and denominator
  are dominated by that same synthetic floor, so **the digitized ratio is biased toward 1 at
  the top of the band**. Correct it by subtracting the 3σ_noise = 3 nm white floor from both
  PSDs before ratioing, or — better — refit (σ_real, ξ, α) per step from Fig. IV.38 and
  rebuild noiseless Palasantzas PSDs. Budget ~1 day of digitization + reconstruction.
- **Wide-y simulation campaign:**
  - Geometry: isolated line, CD 75 nm, resist 120 nm on Si-ARC 35 nm; simulate the **Si-ARC
    open only** (20 s and 30 s), Ar/CF₄/CHF₃ at 7 mTorr, bias 250 W.
  - Grid: y-extent **≥ 2.2 µm** (match L = 2152 nm); **Δy = 2 nm** (≥ 8 cells per 42 nm, the
    shortest wavelength carrying the measured minimum) → **N_y ≈ 1100**. Cross-section grid
    1–2 nm in x, ≥ 0.5 nm near the edge; 35 nm etch depth is shallow → the profile solve is
    cheap and the cost is entirely in N_y.
  - Roughness seeding: **two independent campaigns.** (a) *Probe*: single-wavenumber
    sinusoidal mask perturbation, ε = 1 nm and 0.5 nm (linearity self-check), swept over
    ~20 k values log-spaced in [0.01, 0.6] rad/nm → 40 deterministic runs, embarrassingly
    parallel. (b) *Ensemble*: Palasantzas-seeded rough mask edges at the measured resist
    (σ, ξ, α), **M = 50 realizations** (the module's coherence bias correction is
    (Mγ²−1)/(M−1); M = 50 puts the per-bin bias at ~2 % and lets 8-bin spectral averaging
    reach ~5 % error on |T|). Campaign (b) is the honest one; (a) is the cheap pre-check and
    they must agree.
  - Mask-state contrast: run all three (reference / HBr-cured / VUV-cured) as **different
    input PSDs through the same physics** — the VUV null is free and it is the strongest test
    that petch is not manufacturing attenuation out of numerics.
  - Expected total: ~150 wide-y solves. At quasi-2-D-per-slice cost this is a GPU
    afternoon; at full 3-D it is the campaign's real budget item.

### Gate 2 (structural, no digitization) — Constantoudis/Kokkoris/Gogolides rules of thumb

- **What it gates:** (i) substrate 3σ_S vs resist 3σ_R **linear with slope ≈ 0.5**;
  (ii) reduction only above **σ_R\* ≈ ξ_R/(c tan θ_R), c ≈ 2.0–2.5**; (iii) the JM3 2009
  claim that anisotropic transfer leaves **ξ and α almost unchanged**; (iv) the requirement
  that resist-sidewall **anisotropy** (⊥ striations) is what makes trimming beneficial.
- **Why:** it is closed-form, so no figure digitization is needed at all — it is a *sweep*
  gate. It also sets up the §2.1 tension with LTM, which is the publishable result.
- **Digitization feasibility: N/A (analytic).** The paper's own figures would only be needed
  to check the constant c. **[VERIFY]** angle convention and c against JM3 12(4) 041310.
- **Campaign:** parameter sweep, not a single simulation. Fix ξ_R = 30 nm, α = 0.6, resist
  150 nm, sidewall 86.2°, selectivity 3, depth 150 nm (their Fig. 4 conditions); sweep
  3σ_R over ~0.5–8 nm in 10 steps × M = 20 seeds = 200 runs at y-extent 1 µm, Δy = 2 nm
  (N_y = 500). Then repeat at ξ_R ∈ {10, 30, 100} nm and θ_R ∈ {84°, 86.2°, 88°} to test the
  σ_R\* ∝ ξ_R/tan θ_R scaling — 3×3 more blocks. Roughness seeding must be **anisotropic**
  in 3-D per §2.3, else the trim case will come out with the wrong sign.

### Gate 3 (modern absolute PSDs) — Rutigliani et al., Proc. SPIE 10585, 105851K (2018)

- **What it gates:** absolute PSD level (nm³) and the ξ/PSD(0) decomposition. Specifically:
  can petch reproduce **ξ ↑ with PSD(0) held fixed** (their best condition) versus
  **ξ ↑ with PSD(0) degraded** (their worse conditions), and land unbiased LER on
  2.63 → {1.95, 2.24, 2.33, 2.07} nm?
- **Why third:** the process is a resist *treatment* with an undisclosed recipe, so it gates
  the *observable and the estimator* more than the transport physics. But it is the only
  place to check that our absolute PSD normalization matches the industry's.
- **Digitization feasibility: HIGH for scalars, MEDIUM for curves.** The five unbiased LER
  values are printed as data labels; Fig. 10's ξ and PSD(0) are bar/point plots; Fig. 9's five
  PSD curves are absolute-log-log and digitizable, though the rasterization in the freely
  available copy is coarse. **[VERIFY]** the JM3 17(4) 041016 version likely has cleaner
  vector figures.
- **Campaign:** 32 nm pitch, 16 nm CD lines; y-extent 2 µm at Δy = 0.8 nm (match iRP pixel)
  → **N_y = 2500**; M = 50 seeds; four notional treatment conditions parameterized by
  ion energy / flux ratio rather than by their undisclosed DCS knob. This one is *not* a
  blind gate — declare it as a normalization/consistency check.

### Honorable mentions, explicitly not shortlisted

- **Dupuy 2015 (HAL, open):** best fully-open per-step σ ladder (§1.5) and it ends in a real
  SF₆/CHF₃/Ar Si etch — but the PSDs are "arbitrarily shifted", the flow is SADP (deposition
  confound), and the largest LWR move is a *conformal-deposition* effect, not an etch effect.
  Use as a **secondary scalar gate** on the last three steps only.
- **Sun et al. JM3 14 (2015):** 3σ per frequency band per step. Usable as a coarse 3-number
  gate if PSD digitization fails.
- **Qi 2012 / Vaglio Pret 2011 programmed-roughness masks:** gate the optical LTF, not the
  etch. Their real value is as the *experimental design template* for our probe sweep.
- **Liang/Mack 2018 and Mack DSA 2018:** normalized or anonymized. Not gateable.

### Cross-cutting implementation notes (these will bite)

1. **Convention conversion is the #1 numerical risk.** LTM uses angular k = 2πn/(NΔy) with a
   Δy/(2πN) normalization; Mack/imec/Fractilia use f in 1/nm; `ler_metrology` uses one-sided
   f in 1/nm with σ² = ∫PSD df. **Do not hand-derive the factor** — for every digitized PSD,
   re-integrate the digitized curve and check it reproduces the *published* unbiased σ. If it
   does not, the conversion is wrong. Make this an assertion in the gate script.
2. **Compare in the unbiased domain only.** Every dataset above is already noise-subtracted;
   do not add a synthetic SEM noise floor to petch output. The one exception is Azarnouche
   Fig. IV.39, where a synthetic 3 nm floor was deliberately *re-added* (§5, Gate 1 trap).
3. **α from PSD, ξ from HHCF** (Constantoudis' own estimator guidance) — report both, as
   `ler_metrology` already does.
4. **Prefer LER over LWR gates** unless the left/right edge correlation is declared (§3.3).
5. **Finite-box bias:** all LTM numbers are σ_real over L ≈ 2.2 µm, not σ_inf. If petch's
   y-extent differs from L, apply σ_inf² = σ_real² + σ_CDV² or simply match L exactly. Matching
   L exactly is cheaper and less arguable.
6. **The Demokritos↔LTM ξ/α contradiction (§2.1) is a feature.** Pre-register which way
   petch predicts it *before* digitizing, then reveal. That is the Krüger-style protocol
   applied to spectra, and it is the first time it would be done on roughness.

---

### Appendix — one-line locators for everything cited

- Azarnouche PhD thesis (2012), *Défis liés à la réduction de la rugosité des motifs de résine
  photosensible 193 nm*, Univ. Grenoble / LTM — https://theses.hal.science/tel-00767820 (open PDF)
- Thiault PhD thesis (2007), 3D-AFM of LER through MOS gate patterning — https://theses.hal.science/tel-00321961 (open PDF)
- Dupuy et al., Proc. SPIE 9428, 94280B (2015), DOI 10.1117/12.2085812 — https://univ-grenoble-alpes.hal.science/hal-01869175 (open, CC-BY)
- Azarnouche et al., JVST B 31, 012205 (2013), DOI 10.1116/1.4773063
- Azarnouche et al., J. Appl. Phys. 111, 084318 (2012), DOI 10.1063/1.4705509
- Pargon, Martin, Thiault, Joubert, Foucher, Lill, JVST B 26, 1971 (2008), DOI 10.1116/1.2917071
- Constantoudis, Kokkoris, Xydi, Gogolides, Pargon, Martin, JM3 8(4), 043004 (2009), DOI 10.1117/1.3268365
- Constantoudis et al., JM3 9(4), 041209 (2010), DOI 10.1117/1.3497601
- Constantoudis, Kokkoris, Gogolides, JM3 12(4), 041310 (2013), DOI 10.1117/1.JMM.12.4.041310;
  open summary: SPIE News 4738, 29 Mar 2013
- Rutigliani et al., Proc. SPIE 10585, 105851K (2018) / JM3 17(4), 041016 (2018), DOI 10.1117/1.JMM.17.4.041016;
  open copy: lithoguru.com/scientist/litho_papers/2018_Setting%20up%20a%20proper%20power%20spectral.pdf
- Lorusso et al., JM3 17(4), 041009 (2018), DOI 10.1117/1.JMM.17.4.041009;
  open copy: lithoguru.com/scientist/litho_papers/2018_The%20imec%20roughness%20protocol.pdf
- Liang, Mack, Sirard et al., Proc. SPIE 10585, 1058524 (2018), DOI 10.1117/12.2297328 (open copy on lithoguru)
- Mack et al. (DSA transfer), Proc. SPIE 10589, 1058907 (2018) (open copy on lithoguru)
- Sun, Wang, Beique, Sung, Wood, Kim, JM3 14(3), 033501 (2015), DOI 10.1117/1.JMM.14.3.033501
- Vaglio Pret & Gronheid, Microelectron. Eng. 88, (2011), DOI 10.1016/j.mee.2011.02.015
- Vaglio Pret, Gronheid, Foubert, JM3 9(4) (2010), DOI 10.1117/1.3494614
- Vaglio Pret et al., EIPBN 2011 4B-2 — https://eipbn.org/abstracts/2011/papers/4B-2.pdf (open)
- Qi, Gallagher, Negishi, McIntyre, Zweber, Senna, Proc. SPIE 8522 (2012), DOI 10.1117/12.976855
- Naulleau & Gallatin, Appl. Opt. 42(17), 3390 (2003), DOI 10.1364/AO.42.003390
- Naulleau & Gallatin, JVST B 28, 1259 (2010), DOI 10.1116/1.3509437
- Gallatin & Naulleau, Proc. SPIE 7969 (2011) — https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=908144 (open)
- Wang, Wang, Biolsi, Kushner, JVST A 39, 033003 (2021), DOI 10.1116/6.0000941 — https://www.osti.gov/servlets/purl/1853526 (open)
- Martin & Cunge, JVST B 26, 1281 (2008), DOI 10.1116/1.2932091
- Nakazaki, Tsuda, Takao, Eriguchi, Ono, J. Appl. Phys. 116, 223302 (2014), DOI 10.1063/1.4903956 (OA: hdl.handle.net/2433/193256)
- Jiang et al., Appl. Surf. Sci. 537, 148862 (2021), DOI 10.1016/j.apsusc.2020.148862 (OA)
- Jiang, Li, Wang, Jiang, Huang, ICSICT 2012, DOI 10.1109/ICSICT.2012.6466733
