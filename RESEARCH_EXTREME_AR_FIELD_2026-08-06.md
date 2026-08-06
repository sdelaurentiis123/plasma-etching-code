# The extreme-AR hole-etch field — a research map (2026-08-06)

Scope: the whole 100–200:1 hole-etch problem, which is the domain of the 200:1
deliverable (`HOLE_STUDY_PLAN_2026-08-05.md`, `HOLE_STUDY_RESULTS_2026-08-05.md`).
Not committed. No code touched.

**Evidence convention used throughout.**
- **[Q]** = verbatim quote extracted by me from the primary text (PDF text layer or
  local repo extract). Reliable to the character.
- **[Q-relay]** = wording relayed through the fetch summariser (a small model). The
  *number* is trustworthy to the extent the source page was rendered; the *phrasing*
  may be paraphrase. Treat as `[VERIFY verbatim]` before quoting in a deliverable.
- **[INF]** = my inference/derivation from cited numbers, not a claim in the source.
- `[VERIFY]` = could not confirm to primary source in this pass.

---

## 0. Headline field map

Five statements that organise everything below.

1. **The commercial problem is 100:1 today and the industry is buying its way out of
   200:1 rather than etching it.** Lam states the ONON channel hole is "patterned to
   about 100 nm hole size, and this pushes the etch aspect ratio (depth/width) above
   100:1" [Q, Shen 2023]. The 1000-layer answer being demonstrated is *not* a 200:1
   hole — Kioxia/SanDisk's VLSI 2026 result stacks **two 218-layer arrays** to make a
   "436-layer-equivalent device" [Q-relay], i.e. two ~100:1 etches, not one 200:1 etch.
   Lam's own roadmap language is "1,000-layer 3D NAND with multiple tiers" [Q,
   Counterpoint/Lam 2024]. **200:1 in a single continuous hole is not a demonstrated
   manufacturing regime — it is a frontier.** That is exactly the whitespace a physics
   engine plus a hardware partner can own, and it is also why nobody has published
   validation data there.

2. **The transport reference the industry quotes is the same Clausing number our
   phase-1 study reproduced.** Lam: "For a cylinder with an aspect ratio of 50:1, only
   2.5% of the incoming flux reaches the other end of the cylinder due to diffuse
   collisions with the sidewall" [Q, Shen 2023]. Our study measured **0.025287 at AR 50**
   (`HOLE_STUDY_RESULTS` §2). Agreement to the digit Lam quotes. This is a strong
   external endorsement of the axisymmetric operator, and it is worth saying explicitly
   in the deliverable: *the number the field uses as its neutral-starvation headline is
   the number our operator computes exactly, and we carry it to AR 200 (0.656 %) where
   the field has no published value at all.*

3. **The field's own ARDE curve has the same shape as ours, and Lam attributes the flat
   part to the same cause we do.** Lam: "The ONON channel hole layer … exhibits a flat
   etch rate up to an aspect ratio of ~60 and decreases at higher aspect ratios… The flat
   etch rate at lower aspect ratios can be partially attributed to very low IADs
   (<0.5° 1σ)" [Q, Shen 2023]. Our phase-1 core-only (narrow) beam declines **0.41 % to
   AR 16 and 5.4 % to AR 200**; our tail-bearing beam declines 30.8 % to AR 200
   (`HOLE_STUDY_RESULTS` §5). Lam is describing our core-only branch and calling it the
   experimentally observed regime up to AR 60. **The disagreement in the field is
   therefore not "is ARDE strong" but "how wide is the beam", which is precisely the one
   parameter we sweep and declare open.**

4. **The dominance of wall-reflected delivery at high AR is settled physics in the
   literature, quantified independently, and our numbers are on the same side of it.**
   Kushner/Samsung: "For an AR greater than 3 … the delivery of power to the etch front
   relies more on the hot neutral flux than on the ion flux… The power density delivered
   by hot neutrals to the etch front is larger than the ions by about 20% for AR > 3"
   [Q, Huang 2019]. Lam: "the etching mechanism for aspect ratios beyond 40 or 50 is
   dominated by fast, chemically active species which carry the mechanical energy…
   These species can be ions or neutrals, which originate from ions that were neutralized
   in collisions with the feature sidewall" [Q, Shen 2023]. Our AR-200 result — 82–90 %
   of etch-front energetic flux arrives by wall reflection — is the extrapolation of the
   published trend into a range no one has published. **Nobody contradicts us; nobody has
   gone there.**

5. **Deep-AR validation data barely exists, and what exists tops out around AR 40–60.**
   The deepest *measured, process-characterised, digitisable rate-or-depth-vs-AR* data
   found this pass is **Nguyen/Jansen (de Boer lineage) 2020, silicon CORE, gap AR ≈ 57
   with a six-point depth-vs-time series to 12.0 µm** (§2.1). The deepest *3D in-feature
   geometry measurement* is **Okawa/Takahashi APL 2026, tender X-ray ptychography-CT of
   ~10 µm HAR holes at 67.5 nm 3D / 38.5 nm 2D resolution** (§2.3) — published five weeks
   ago. The deepest *published coupled profile evolution* is AR 80 (Kushner, simulation
   only, Fig. 9 of Huang 2019). **There is no published validated model at 100:1, and no
   published anything at 200:1.**

---

## 1. State of practice — the real 3D-NAND channel-hole etch

### 1.1 The dimensions, from the vendors

| quantity | value | source |
|---|---|---|
| hole CD | "about 100 nm hole size" [Q] | Shen 2023 (Lam) |
| aspect ratio, current | "pushes the etch aspect ratio (depth/width) above 100:1" [Q] | Shen 2023 |
| depth, cryo production | "etch memory channels as deep as 10 microns" [Q-relay] | Lam Cryo 3.0 (2024-07-31) |
| depth, TEL cryo | "10-µm-deep etch … in just 33 minutes" [Q-relay] | TEL press, VLSI 2023 |
| stack, shipping | "more than 150 layers in 2022" [Q]; 321-layer QLC in MP since Aug 2025 [Q-relay] | Shen 2023; TrendForce |
| stack, roadmap | "500 to 1000 layers are on the horizon with double- or triple-deck architectures" [Q] | Shen 2023 |
| wafer temperature | "temperatures as low as −60 °C" [Q, Lam]; "−70 °C" [Q-relay, TEL] | Counterpoint/Lam 2024; TEL |
| hole count | "forming trillions of perfect channel holes from top to bottom" [Q] | Counterpoint/Lam 2024 |
| etch time per wafer-deck | "roughly an hour to drill the holes through one wafer" at 218 layers [Q-relay] | kantenna/VLSI 2026 [VERIFY] |

**[INF]** 10 µm / 33 min ⇒ **≈ 0.30 µm/min mean rate at AR ≈ 100**. 10 µm in ~60 min at
218 layers ⇒ ≈ 0.17 µm/min. These are the only two public absolute-rate anchors at
AR ≈ 100 I could find. Both are press-level, neither carries a condition table.

### 1.2 The failure modes, named and ranked by the vendors

Lam's own taxonomy [Q, Shen 2023]: "Channel holes at very HARs can exhibit non-ideal
hole behavior such as twisting (off center), roughness (striations), and distortion
(deviation from circle shape)."

Root-cause claims, verbatim:
- **Mask, not the hole, is the first cause of roughness.** "Polymer deposition on the
  sidewalls of the mask, so-called necking, is considered the main root cause of sidewall
  roughness." [Q]
- **Distortion needs ion scattering.** "In the absence of ion scattering with a circular
  hard mask shape (0× ion scatter), the hole is perfectly centered and circular. At 1×
  and 2× ion scatters, hole shape distortion is apparent at etch depths equal to or
  deeper than 4000 nm." [Q]
- **Twisting needs mask evolution *and* scattering.** "In the absence of hard mask
  evolution and for the etch depth studied, no significant twisting of the hole feature
  was observed." [Q] … "Combining mask shape evolution and ion scattering in the hole in
  the model revealed that hole distortion and twisting could occur and were magnified at
  deeper etch depths." [Q]
- **Kushner's competing claim: twisting can be charging-free.** "the feature, resulting
  in twisting irrespective of charging" [Q, Huang 2019 L222]; and separately, with
  charging on, "feature distortion is worse with charging" [Q, Huang 2019 L2104].

**This is a live disagreement in the field** — Lam says mask-morphology + scattering,
Kushner says stochastic flux + charging aggravation. Both are simulation-led. Neither
has an in-feature measurement to settle it. See §5.

### 1.3 The commercial CD specification — and why it is not usable as stated

Lam's headline: "Predictably and repeatedly etch uniform memory channels as deep as 10
microns with less than 0.1% deviation in the channel's critical dimension from the top to
the bottom" [Q, Counterpoint/Lam 2024].

**[INF] This number is ambiguous and must not be adopted as a target without a
definition.** If "0.1 % of CD" then for a 100 nm CD the tolerance is 0.1 nm top-to-bottom,
which is below any inline metrology capability quoted anywhere in this document
(best published inline: sub-nanometre *tilt* precision by SAXS; best 3D imaging:
67.5 nm voxel resolution). If "0.1 % of depth" then it is 10 nm on a 10 µm hole, i.e.
±10 % of CD — plausible and consistent with the rest of the literature. **Flag as
[VERIFY] and quote both readings in any partner-facing document.**

The metric the *technical* literature actually optimises is different and well defined:

> "For HAR profile control, BB bias, defined as the difference between the largest bow CD
> and the bottom CD, is often used as the key metric to define how vertical a hole profile
> is. A near zero BB bias indicates a perfectly vertical profile." [Q, Shen 2023]

Published BB-bias movements (all Lam, experimental, ONON):
- bottom 10 % of SiN tiers replaced with softer film → **BB bias improved by 35 %** [Q]
- bottom 20 % replaced → **BB bias reduction of up to 53 %** [Q]
- top-tier SiN swap → **bow CD reduction of 18 %** [Q]
- carbon liner insertion → **"improvement of up to 30 % in BB bias"** [Q]

**These four numbers are the single best available "commercially meaningful" scorecard
for a profile model at HAR.** They are relative, which is good — they do not require the
absolute-depth channel we know runs +29 % on the Krüger validation case. A model that
reproduces the *sign and rough magnitude* of a bottom-tier-softening BB-bias improvement
is doing something no open-source engine has been shown to do.

### 1.4 Selectivity and the mask, the actual scaling wall

- Lam: "The normalized bulk selectivity is reduced by roughly 50 %, going from an aspect
  ratio of 40 to 140." [Q, Shen 2023, Fig. 7] — **this is the deepest AR axis quantified
  in any published figure I found (AR 140).** It is model output matched to experiment,
  not a direct measurement.
- Kushner: "the PR mask is nearly fully eroded when reaching an AR of 60… This loss of
  selectivity is not because the PR etches faster, but rather because the SiO2 etches
  slower." [Q, Huang 2019]
- Kushner also records a **counterintuitive relief term** we do not model: "As the mask
  erodes, the conductance limit of neutrals into the SiO2 portion of the feature relaxes
  and the unimpeded ion flux to the etch front increases. ARDE would be even more severe
  in the absence of PR erosion." [Q]

**[INF] Consequence for our 200:1 study.** Our phase-1 geometry is a fixed straight
cylinder with no mask erosion. Kushner's statement means a fixed-mask study **overstates
ARDE relative to a real etch** in the direction of the mask term, while a fixed *narrow*
mouth understates the necking term. Both are declared-open in `HOLE_STUDY_RESULTS`; this
is the literature justification for why they are not second-order.

### 1.5 Cryo — what is actually claimed, and what is actually published

Vendor claims (press-level):
- Lam Cryo 3.0: "etch rates two and a half times faster than conventional plasma etch
  technologies" and "two times faster than Lam's previous-generation cryogenic etching
  solutions" [Q-relay]; "over five million wafers" produced on cryo [Q-relay]; "nearly
  1,000 [chambers] in production with cryogenic etch technology" [Q, Counterpoint].
- TEL: 10 µm in 33 min, −70 °C, "84 % reduction of global warming potential" [Q-relay],
  presented at 2023 Symposium on VLSI Technology and Circuits.

Peer-reviewed cryo mechanism (the part that is real physics, not marketing):
- **Shen 2023 (Lam), the only quantitative blanket anchor:** "The blanket silicon oxide
  and nitride etching rates increase by about 20 % for a lean HAR etch process when
  lowering the temperature from 60 °C to −20 °C." [Q] And on ARDE: "The etching rate
  diminishes much more quickly for the high-temperature process than for the
  low-temperature process." [Q, Fig. 3]
- Their two candidate mechanisms, verbatim: "Coburn and Winters predicted a lower ARDE
  when the reactive sticking coefficient is reduced. An additional explanation for the
  improved ARDE performance at lower temperatures is potentially surface diffusion of
  neutrals." [Q]
- **The HF/H₂O claim.** Lill *et al.* (Lam + George/CU Boulder + PPPL), *Low-temperature
  etching of silicon oxide and silicon nitride with hydrogen fluoride*, JVST A **42**,
  063006 (2024) — pure HF etches SiN; SiO₂ needs an added F source (e.g. PF₃).
  *Reaction mechanism for HF based cryogenic plasma etching of SiO₂*, JVST A **44**,
  033006 (2026): "The reaction between HF and SiO₂ generates H₂O, and at cryogenic
  temperatures, the increased physisorption rate enables a several nm thick adsorption
  layer… the adsorbed H₂O layer is thought to be essential for catalyzing HF hydrolysis
  and silanol group formation" [Q-relay]; "At temperatures lower than −100 °C, the
  solidification of H₂O reduces the co-adsorption of HF, decreasing the etch rate"
  [Q-relay]. Operating window quoted as "−10 to −100 °C" [Q-relay].

**Two structural implications for us.**
(a) The cryo chemistry is a *condensed-phase, several-nm adsorbate* mechanism. Our
    fluorocarbon Deck-1 mechanism is a surface-complex mechanism. **Cryo is a new module,
    not a parameter change**, and the 2026 JVST A paper is the mechanism reference to
    build it from.
(b) Cryo's benefit is claimed to come *through the ARDE channel* (reduced sticking →
    deeper neutral penetration). Our engine measures neutral penetration exactly
    (Clausing-gated). **A cryo module bolted onto our exact transport is a clean, testable
    story: predict the ARDE-curve flattening from a sticking-coefficient change, with
    zero fitted transport.**

---

## 2. Measured data at extreme AR that can serve as gates

Ranked by depth reached and by how much of it survives a "can we actually score against
it" test.

### 2.1 The deepest usable rate/depth-vs-AR dataset found: Nguyen & Jansen 2020 (CORE)

Nguyen, Shkondin, Jensen, Hübner, Leussink, Jansen — *Ultrahigh aspect ratio etching of
silicon in SF₆–O₂ plasma: the clear-oxidize-remove-etch (CORE) sequence and chromium
mask* — **J. Vac. Sci. Technol. A 38(5), 053002 (2020)**, DOI 10.1116/6.0000357.
Affiliations DTU Nanolab + MESA+/Twente. Open PDF mirrored at DTU Nanolab labadviser.

Verbatim depth-vs-time series (§III D, extracted from the PDF text layer) [Q]:

> "After 1 h [Fig. 8(a)], the etch pillars are 1.1 μm high and slightly negative tapered.
> After 2 h [Fig. 8(b)], the pillars are 2.1 μm high and almost perfectly straight. 4 h of
> etching [Fig. 8(c)] results in 3.9 μm etch depth with a slightly positive slope. 8 h
> [Fig. 8(d)] gives 7.1 μm etch depth and a more positive slope. 12 h [Fig. 8(e)] results
> in a height of 9.5 μm and an even more positive slope. Finally, 16 h [Fig. 8(f)] results
> in a height of 12.0 μm and the most positive slope."

And the maximum: "an array of 200 nm pillars with 400 nm periodicity is etched to a height
of 11.4 μm. This corresponds to an aspect ratio of 57." [Q]
Abstract summary: "aspect ratios beyond 55 for gaps and up to 200 for pillars" [Q].

**[INF] Derived instantaneous-rate table** (assuming the Fig. 8 series is the 200 nm-dot /
400 nm-pitch DUV array, i.e. 200 nm gap — `[VERIFY]` this pairing against Fig. 8's
caption/scale bar before use):

| interval | Δdepth | mean rate | mid-depth | gap AR at mid-depth | rate / rate(0–1 h) |
|---|---|---|---|---|---|
| 0–1 h | 1.1 µm | 1.10 µm/h | 0.55 µm | 2.8 | 1.000 |
| 1–2 h | 1.0 µm | 1.00 | 1.6 µm | 8.0 | 0.909 |
| 2–4 h | 1.8 µm | 0.90 | 3.0 µm | 15.0 | 0.818 |
| 4–8 h | 3.2 µm | 0.80 | 5.5 µm | 27.5 | 0.727 |
| 8–12 h | 2.4 µm | 0.60 | 8.3 µm | 41.5 | 0.545 |
| 12–16 h | 2.5 µm | 0.625 | 10.75 µm | 53.8 | 0.568 |

**Rate falls ~43 % from AR ≈ 3 to AR ≈ 54 and then flattens.** This is the deepest
*measured* rate-vs-AR curve I located, and its shape — slow decline, plateau at high AR —
is qualitatively the *opposite* of Huang's MCFPM (80 % fall by AR 40) and is closer to our
transport-only decline. That is a genuinely important observation for the honesty appendix
of the hole study: **the one clean literature anchor we currently cite for "our ARDE is too
weak" (Huang's 80 %-by-AR-40) is contradicted by the deepest published measurement, which
declines only ~43 % by AR 54.**

**Caveats that must ship with any use of this data.**
1. It is a **pillar array (open, connected gaps)**, not an enclosed hole. Transport is much
   less confined than a cylinder; do not treat the AR number as equivalent to hole AR.
2. It is **silicon SF₆/O₂ cyclic (CORE)**, not fluorocarbon oxide — matches our Belen/Si
   arm, not Deck 1.
3. Mask retraction is explicit and large: "the erosion rate is closely following the
   expected 2 nm/h. But after that … between 12 and 16 h, it becomes only 0.25 nm/h. The
   authors have no explanation for this." [Q] The profile "becomes more positive" with
   time [Q], so the curve convolves geometry evolution with rate.
4. The repo's existing de Boer gate (`DEBOER_DIVERGENCE.md`: normalised floor rate
   1.0 / 0.43 / 0.29 / 0.20 at AR 0/10/20/40) is from a **different** de Boer/Blauw
   dataset and is *much* steeper than the 2020 CORE series. **These two "de Boer"
   anchors disagree with each other by ~2× at AR 40.** Reconciling them is a concrete,
   cheap, high-value task. (Different chemistry regime, different geometry, and one is a
   trench floor while the other is a pillar gap — but the repo currently leans on one
   without acknowledging the other.)

### 2.2 The deepest *AR axis* in any published figure: Lam Shen 2023

- Etch-rate-vs-AR: flat to AR ≈ 60 then declining, plotted against AR on Fig. 6(a) [Q].
- Normalised ion flux to the etch front vs AR for several IAD widths (1σ Gaussian),
  Fig. 6(b) — **this is directly comparable to our delivery-vs-AR table** and is the
  single most valuable digitisation target found this pass.
- Selectivity vs AR out to **AR 140**, Fig. 7 [Q].

Status: **semi-empirical Monte-Carlo model output, tuned to experiment** ("Parameters can
be adjusted to match the simulated profile to experimental results at fixed time
intervals" [Q]). So it is an **S-grade** reference (simulation), not E-grade. But it is
S-grade *from the tool vendor with the experiment in hand*, at the deepest AR anyone
plots. Digitising Fig. 6(b) gives us a delivery-vs-AR curve **as a function of beam
width**, which is exactly our declared-open axis.

### 2.3 The newest and most valuable in-feature metrology: tender X-ray ptychography-CT

Okawa, … , Takahashi (Tohoku University / SRIS / NanoTerasu BL10U) —
*Nondestructive structural evaluation of high-aspect-ratio-etched holes by tender x-ray
ptychography with computed tomography* — **Appl. Phys. Lett. 129(1), 011109 (2026)**,
published 6 July 2026.

Relayed content [Q-relay, `[VERIFY verbatim]`]:
- "tender x-ray ptychography for the characterization of high-aspect-ratio (HAR)-etched
  hole structures with depths of approximately 10 μm"
- "Reconstructed two-dimensional projection images, with a spatial resolution of 38.5 nm,
  reveal features such as necking at the amorphous carbon mask and tapered hole profiles"
- 3D: "A three-dimensional electron density map was reconstructed using ptychographic
  x-ray computed tomography, achieving a spatial resolution of 67.5 nm, enabling direct
  visualization of etch-rate variations associated with aspect-ratio-dependent etching"
- "provides a powerful and versatile approach for the nondestructive 3D evaluation of
  ARDE-induced structural non-uniformities in HAR hole structures"
- Tender X-ray regime given as 2–5 keV.

**Why this is the most important new item in the whole map.** It is (i) a *hole*, (ii)
~10 µm deep, (iii) with an **amorphous-carbon mask and a resolved neck** — which is our
mouth-equilibrium modality (`RESULTS_MOUTH_EQUILIBRIUM_PROBE`, `RESULTS_NECK_REGRADE`) —
and (iv) 3D, non-destructive, so it can in principle produce *hole-to-hole statistics*
rather than one FIB cross-section. If the paper's supplementary contains a CD-vs-depth
trace, **it is the first publicly available in-feature profile at AR ≈ 100 that our
`regrade_neck_metrics` definitions can be scored against directly.** Chase the full text
and any supplementary; if the trace exists, this becomes the highest-value held-out gate
in the repo.

### 2.4 In-feature *transport* measurement — the closest thing that exists

Cunge, Darnon, Dubois, Bezard, Mourey, Petit-Etienne, Vallier, Despiau-Pujo, Sadeghi —
*Measuring ion velocity distribution functions through high-aspect ratio holes in
inductively coupled plasmas* — **Appl. Phys. Lett. 108, 093109 (2016)**.

Method, verbatim from the Impedans application note VE01 (Gahan, Scullin, Hopkins) [Q]:

> "A series of capillary plates were placed on the RFEA surface, each with a specified AR
> through which the IED was measured. The ion energy distributions were measured at the
> wafer surface in an inductive coupled plasma (ICP) reactor from Applied Materials. A
> helium plasma was used at a pressure of 10 mTorr. ICP power of 750 W and bias power set
> for a DC self-bias voltage of 100 V."

> "A drop in total ion flux (area under the curve) and a change in the shape of the IED as
> the AR is increased are the main observations. An analytical model, in which the ion
> temperature T_i is the main variable, gives excellent agreement with experiment for a
> chosen T_i input."

Impedans' Vertex generalisation replaces the physical capillary with an **electrostatically
synthesised effective AR** (potential difference between grids 2 and 3), giving
AR-resolved IEDs without venting [Q].

**[INF] This is a direct experimental analogue of our cone-acceptance table.** Our study
predicts, from an exact geometric cone plus the two-component beam, the fraction of beam
transmitted at each AR (0.9669 → 0.5736 for core-only from AR 100 → 200; 0.4888 → 0.2421
at tail 0.65). Cunge's method measures exactly that transmitted fraction against AR, and
inverts it for T_i. **A capillary-plate transmission-vs-AR curve is a preregistrable, zero-
chemistry, zero-charging blind gate for our transport operator + beam model, and I found
no evidence anyone has run it above AR ~50.** This is arguably the cheapest experiment a
hardware partner could run that would validate the exact part of our engine that is
strongest.

### 2.5 The beam itself — the measurements our tail fraction hangs on

The Nagoya group (with **KIOXIA Corporation Frontier Technology R&D Institute** as
co-authors on the PSST paper) has published the highest-resolution IAD measurements in the
field, and the PPPL group has just published the theory that explains them.

- Ichikawa *et al.*, *Angular distribution measurement of high-energy argon neutral and
  ion in a 13.56 MHz capacitively-coupled plasma* — **Appl. Phys. Express 14, 126001 (2021)**.
- Kim, Kawamura, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda — *Measurement of
  energy-resolved ion angular distribution in a dual-frequency capacitively coupled argon
  plasma* — **Jpn. J. Appl. Phys. 64, 05SP15 (2025)**.
- Kim, Kawamura, Fujitani, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda — *Influence of
  sheath collisions on the ion angular distributions in a dual-frequency capacitively
  coupled plasma* — **Jpn. J. Appl. Phys. 64, 096002 (2025)**, DOI 10.35848/1347-4065/ae0105.
- Kawamura, Kim, Ichikawa, Suzuki, Iino, Fukumizu, Kurihara, Toyoda — *Angular distribution
  measurement of high-energy ions and neutrals impinging on an rf electrode in a
  dual-frequency capacitively-coupled Ar plasma* — **Plasma Sources Sci. Technol. 34,
  055006 (2025)**, DOI 10.1088/1361-6595/add321. Ion energies "up to approximately 2.5 keV"
  [Q-relay].

**The finding that matters most to us**, from JJAP 64, 096002 [Q-relay,
`[VERIFY verbatim`]:

> "Measured IADs comprise both narrow main-components and broader tail-components…
> Main-to-total component ratios at maximum energy decreased exponentially with pressure…
> collisions between Ar⁺ ions and neutral Ar atoms in the sheath generate the
> tail-component features."

**[INF] This is a route to retire our single declared-open parameter.** `HOLE_STUDY_PLAN`
lists the tail fraction as "the one declared open parameter" and points to an S2 rung of
deriving it from a sheath collision operator. Kim 2025 reports the **measured** main-to-
total ratio as an **exponential function of pressure** — i.e. the tail fraction is a
*measured* quantity with a *measured* scaling law, in the same apparatus whose widths we
already reproduce. Extracting that exponential and quoting the tail fraction as
`f_tail(p)` from a published measurement, rather than as a swept band, converts our
largest honesty caveat into a citation. **This is the single highest-leverage follow-up in
this document.**

And the theory that closes it: Khrabrov & Kaganovich (PPPL), *Ion-neutral and
neutral-neutral scattering in argon at KeV energies and implications for high-aspect-ratio
etching*, **arXiv:2604.04214v2**, 8 April 2026 [Q]:

> "we interpret the 'tail' observed outside of the main peak in the measured angular
> distributions as being due to the finite angular spread in ion-neutral and
> neutral-neutral elastic collisions. The differential scattering in these collisions is
> determined by the interaction potential."

> "In [15], the authors demonstrated that the 'tail' in the ion distribution cannot be
> accounted for within the scattering model of [38]. Such is indeed the case, as the
> elastic scattering predicted by that model at high energies is much smaller than
> predicted by the model based on a realistic potential."

Acknowledgment, verbatim [Q]: "This research was supported in part by **Samsung** and the
U.S. Department of Energy through the PPPL CRADA agreement with **Applied Materials**".

### 2.6 Statistical / metrology datasets (geometry targets, no boundary)

- **CD-SAXS**: *CDSAXS study of 3D NAND channel hole etch pattern edge effects and etched
  hole pattern variance*, Proc. SPIE **12955**, 1295539 (2024), DOI 10.1117/12.3010927.
  Measures "CD, CD profile, tilt and distortion, specifically for holes at pattern edge
  (outer holes) and holes inside pattern (inner holes)" and "hole etch behavior, especially
  versus etch depth, and its impact on the final hole pattern variance" [Q-relay].
  Companion: *Inline metrology of high aspect ratio hole tilt using small-angle x-ray
  scattering*, Proc. SPIE **12053**, 1205312 (2022) — "sub-nanometer precision" on tilt
  without a structural model [Q-relay].
- **FIB-SEM tomography**: Zhang, Lee, Klochkov, Korb, Sorkhabi, Lan, Pichumani,
  Tekleyohannes, Wang, Sallis, Ningen, Kim, Teoh, Polubotko, Pirkle, Foca — *3D
  reconstruction of 3D NAND memory etch profiles using FIB-SEM: identifying variances in
  etched hole patterns*, Proc. SPIE **13426**, 134262J (2025), DOI 10.1117/12.3051156.
  Explicitly frames the cryo failure list: "Cryogenic etch processes can improve throughput
  and reduce costs, but they also introduce issues such as inner/outer hole loading,
  twisting, bowing, and incomplete or partial etching of holes" [Q-relay].
- **TEL experimental 3D HARC profiles**: Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda
  (Tokyo Electron Miyagi) — *Precise and practical 3D topography simulation of high aspect
  ratio contact hole etch by using model optimization algorithm* — **JVST A 42(4), 043003
  (2024)**, DOI 10.1116/6.0003515. "as the aspect ratio increases, novel issues, such as
  'distortion' and 'twisting,' have been highlighted"; models fitted by "matching both
  vertical and horizontal cross-sectional profiles carefully to the experimental results"
  [Q-relay]. **Horizontal (azimuthal) cross-sections are the rare thing here** — almost
  nobody publishes them.

**Twisting magnitudes.** The only numeric twisting figures surfaced were relayed from a
patent-family text ("twisting in the x direction of 2.40 nm in a layout, while a control
that used etching alone was found to have a twisting in the x direction of 6.4 nm"
[Q-relay, `[VERIFY]` — patent source, treat as indicative only]). **There is no
peer-reviewed table of twist-vs-depth statistics in the open literature that I could
find.** That is whitespace (§5).

---

## 3. Modeling state of the art at AR > 50

### 3.1 Who has actually gone deep, and how deep

| group / code | deepest AR shown | what it is | charging? | hot neutrals? | validation claim |
|---|---|---|---|---|---|
| Kushner MCFPM (Michigan + Samsung) | **AR 80** profiles (Huang 2019 Fig. 9); AR 40 base case with full flux/energy diagnostics | 3D voxel MC feature + HPEM reactor | **yes**, Poisson/SOR, per-cell ε and μ± | **yes**, explicit hot-neutral partners | "The predicted etch times here are 36 min … and 48 min … which are in reasonable agreement with experiments" [Q] |
| Lam semi-empirical MC (Shen 2023) | **AR 140** (selectivity curve); ARDE curve past 60 | 3D MC feature, fitted to SEM at fixed time intervals | not stated | ion scattering yes; hot-neutral bookkeeping not described | "the simulated etch profile progression and aspect ratio dependence matched the experimental data" [Q] |
| TEL cell-based particle MC (Nishizuka 2024) | HARC (AR not stated in abstract) | 3D cell MC + model-optimisation algorithm | not stated | not stated | fitted to vertical **and horizontal** SEM cross-sections [Q-relay] |
| K-SPEED (Chungnam Nat. Univ., You group) | HAR trench+hole, physical sputtering | 3D feature-scale, **GPU-parallel, ~200× vs 1 CPU** [Q-relay] | not in the sputtering paper | n/a (pure sputtering) | "the etch front in the trench reaches a greater depth than in the hole" [Q-relay]. Phys. Plasmas **33**, 033901 (2026) |
| Vanraes, Venugopalan, Besemer, Bogaerts (Antwerp) | low AR (≈1 focus) | multiscale + experiment, SiO₂ CHF₃/Ar, CF₄/Ar | explicitly **no** | diffuse reflection yes | "shadowing and diffuse reflection of neutrals" dominant, "without taking into account charging effects and the polymer layer thickness" [Q-relay]. PSST **32**, 065003 (2023) |
| ViennaPS (TU Wien) | published validations are sub-micron trenches | level-set + MC ray tracing (Embree/OptiX) | **no charging module** (see `RESEARCH_SOTA.md`) | no explicit hot-neutral channel | open-source; no HAR-hole validation published |
| Mohanty, Maheshwari, Mohapatra (IIT) | 3D NAND pillar/hole | feature-scale | **yes** ("Impact of Surface Charging on Feature Profile Evolution in 3D NAND Pillar Etching", SISPAD 2026, accepted) | ion reflection (ASMC 2026) | not yet public |

**Bottom line: the published ceiling for coupled profile evolution is AR 80 (simulation,
Kushner 2019). The published ceiling for any AR-resolved curve is AR 140 (Lam, fitted
model). Nothing at 200:1 exists in the open literature.**

### 3.2 What the deep models include and omit

Kushner's AR-40 base case is the most completely instrumented deep-AR calculation in
print. Its diagnostics map almost one-to-one onto ours:

| Kushner AR-40 diagnostic (verbatim) | petch phase-1 counterpart |
|---|---|
| "This shadowing contributes to a decrease in the ion flux to the etch front from 2.0 × 10¹⁵ to 0.3 × 10¹⁵ cm⁻² s⁻¹" [Q] (AR 0→40) | direct-ion delivery falls to 0.0720 at AR 200 (tail 0.65) |
| "The flux of hot neutrals to the etch front increases from 3.1 × 10¹⁵ to 8.0 × 10¹⁵ … as the etch depth increases from 0 to 480 nm (AR = 4)" then "decreases to 1.1 × 10¹⁵ cm⁻² s⁻¹" at AR 40 [Q] | cascaded hot-particle delivery, non-monotone in AR (our funnelling result, `RESULTS_CASCADE_FUNNELLING`) |
| "The average energy of hot neutrals reaching the etch front first decreases with increasing AR and is then maintained at about 400 eV at ARs higher than 3" [Q] | our cascade retains energy via the Kress reaction rule; **we do not currently report an average cascaded energy vs AR — we should, it is directly comparable** |
| "the flux of CFx and CxFy radicals to the etch front decreases from 3.1 × 10¹⁶ to 0.4 × 10¹⁶ cm⁻² s⁻¹" (AR 0→10) and "from 4.4 × 10¹⁵ to 0.9 × 10¹⁵" (AR 10→40) [Q] | thermal delivery 0.656 % of source at AR 200 (Clausing-exact) |
| "the etch rate then decreases by 80 % by the time the AR reaches 40" [Q] | transport-only decline ≈ 10 % at AR 40 (tail 0.65) — the acknowledged gap |
| "the maximum electric potential higher in the feature (AR ≈ 10) is 1100 V, which is about 60 % of the average energy of incident ions" [Q]; "the maximum potentials at ARs of 10–20 are 200–400 V" [Q]; "For an AR of 40, charging produces a decrease of the average energy of ions to the etch front from 1940 to 1050 eV" [Q] | our charging ceiling ("in-feature potential never exceeds the maximum ion energy and is nearly AR-independent above AR ≈ 17", `RESEARCH_CHARGING_DEEP_AR_VALIDATION`) — **consistent in structure and in the ceiling claim** |
| "the etching of HAR vias in silicon oxide having an AR of 40–50 takes as long as 40–50 min" [Q] | the only published absolute-time anchor at AR 40–50 |
| "the ions and hot neutrals have larger fluxes to the etch front than those of CFx and CxFy by a factor of 2–3" at AR > 10 [Q] | our energetic:thermal ratio of **105–144× at AR 200** |

**[INF] The last row is the most consequential comparison in this document.** Kushner
measures energetic:thermal ≈ 2–3× at AR 10–40. We compute 105–144× at AR 200. The
literature offers no intermediate point. **Our claim that the front is "energy-rich and
radical-poor" by two orders of magnitude at 200:1 is a genuinely novel, falsifiable
prediction — it is a 40–70× extrapolation beyond the deepest published ratio.** That
should be stated as the study's boldest scientific claim, with the extrapolation factor
named.

Kushner's charging conclusion also directly supports the study's "not modelled" framing:
"For high AR features, grazing incidence collisions of ions on sidewalls depositing charge
produce electric potentials with maxima on the sidewalls (as opposed to the bottom) of the
feature" [Q, abstract] and "the majority of ions will have collided with the sidewalls
where they deposit charge, then proceeding as hot neutrals which are not affected by the
electric fields" [Q]. **[INF] I.e. at AR ≫ 10 the dominant delivery channel is charge-
immune by construction** — which is the strongest available literature argument that a
charging-free deep-AR energetic-delivery calculation is not a fatal omission, and it is
Kushner's own sentence, not ours. Worth quoting verbatim in the honesty appendix.

### 3.3 Where our engine's combination would stand

Our combination is: **exact axisymmetric transport (Clausing-gated) + measured
two-component beam + E8 thermalised radical return + zero-knob chemistry deck**, with
every claim receipted.

Against the field:

| axis | field SOTA | petch | verdict |
|---|---|---|---|
| transport accuracy | MC ray tracing with statistical noise; Lam quotes a textbook Clausing value | quadrature-exact band algebra, agrees with Santeler to 3.3e-3 at AR 200, flux closure 2e-14 | **best in field, by a clear margin, and provably so** |
| depth of transport validation | AR ~50 (Lam's quoted number), AR 40–80 (Kushner MC) | AR 200 against analytic theory | **deepest gated transport anywhere** |
| beam model | Lam: single Gaussian, 1σ swept. Kushner: HPEM IEAD. Nobody uses a two-component measured IAD in a feature model | measured core+tail reproduced against 0.1°-resolution data | **unique**; and PPPL 2026 now supplies the first-principles justification for the tail |
| hot-neutral / cascade | Kushner: yes, mature, with energy bookkeeping. Lam: ion scattering, less documented | yes, exact, flux-closed to 2e-14, but **bounce cap unconverged at AR 200 (+1.6 % at cap 64)** | **parity on physics, better on receipts, worse on convergence at 200:1** |
| thermalised radical return (E8) | not described in any of the deep papers I read | wired, published sticking, ~2.2× delivery gain at a trench | **plausibly unique; needs the hole geometry to be worth claiming** |
| charging | Kushner: full Poisson, AR 40, quantitative potentials | ceiling result + gate ladder, **not in this pipeline** | **behind at deep AR; the gap is named and bounded** |
| profile evolution at HAR | Kushner AR 80, Lam AR 140 curves, TEL 3D fitted | **not run** — no axisymmetric evolution driver | **the one place we are simply absent** |
| absolute rate | Kushner "reasonable agreement" on 36/48 min vs experiment | +29 % on the validation depth channel; ion-limited surface law | **behind** |
| knob discipline | Lam explicitly fits ("Parameters can be adjusted to match the simulated profile to experimental results at fixed time intervals" [Q]); Kushner runs a physics-informed optimiser | zero fitted knobs, one declared-open swept band | **best in field, and it is a real differentiator, not a slogan** |

**Honest summary sentence for a partner deck.** *"On transport we are the most accurate
and by far the deepest-validated engine in the field, with a receipt at 200:1 that no one
else has attempted. On beam physics we uniquely carry a measured two-component IAD whose
tail was explained from first principles by PPPL three months ago. On profile evolution at
extreme AR we have not yet run, and the field's best (Kushner AR 80, Lam AR 140) is
ahead of us there. On rate we are +29 % and ion-limited, and the field's deep models are
neutral-transport-limited — which is the defect we have already located."*

### 3.4 The one place the field's numbers should change our priors

Lam's Fig. 6(b) — normalised ion flux to the etch front vs AR for several Gaussian IAD
widths — is the field's own version of our sensitivity table, and Lam's conclusion is that
**a sub-0.5° 1σ beam gives a flat ion-flux-vs-AR curve up to AR ≈ 60**, matching their
experiment. Our core-only branch gives 0.41 % decline to AR 16 and 5.4 % to AR 200 — the
same qualitative statement. **Our `HOLE_STUDY_RESULTS` §5 gate ("a configuration whose
delivery curve is flat must not be used to quote depth or profile predictions") is
therefore, read literally, a gate that would disqualify Lam's own published operating
point.** That gate should be re-worded: the flat curve is not unphysical, it is what a
narrow beam genuinely does; what is disqualified is *using a flat-delivery configuration
to attribute depth trends to transport*. Small wording change, removes an unnecessary
conflict with the best experimental reference in the field.

---

## 4. Neutral-beam and alternative-source etch at HAR

This is the comparison set for a neutral-beam tool's pitch.

### 4.1 What the Samukawa lineage actually demonstrated

Primary references (from Khrabrov & Kaganovich's bibliography, verified there):
Samukawa, Sakamoto, Ichiki, JJAP **40**, L779 (2001); Samukawa, JJAP **45**, 2395 (2006);
Samukawa, Appl. Surf. Sci. **253**, 6681 (2007); Kubota, Nukaga, Ueki, Sugiyama, Inamoto,
Ohtake, Samukawa, *200-mm-diameter neutral beam source based on inductively coupled plasma
etcher and silicon etching*, JVST A **28**, 1169 (2010); Miwa, Nishimori, Ueki, Sugiyama,
Kubota, Samukawa, *Low-damage silicon etching using a neutral beam*, JVST B **31**, 051207
(2013); Samukawa, IEEE Open J. Nanotechnol. **3**, 133 (2022).

Performance, as characterised by an independent group (PPPL, 2026) [Q]:

> "Samukawa and co-authors have designed, studied, and successfully applied atomic beam
> sources that couple a surface neutralizer with an inductively coupled discharge… **These
> sources showed good performance at low-to-medium aspect ratio etching of nanometer
> features.** Graphite neutralizer surfaces are more effective in neutralizing negative
> ions… For positive ions, a high rate of neutralization, more than 50%, was also obtained.
> **The angular divergence of the resulting FAB was estimated, apparently based on
> simulations, to be around 5° which is likely near the minimum achievable for
> surface-neutralizing arrays.**"

Quantitative anchors [Q-relay from the JVST abstracts]:
- Kubota 2010: "An Ar neutral beam flux of more than 1 mA/cm² in equivalent current density
  and a neutralization efficiency of more than 99 % were obtained"; F₂-based neutral beam
  Si etch rate "about 47 nm/min"; Cl₂-based "completely no undercut".
- Miwa 2013: "Etch rate decreased with increasing Si trench aspect ratio. This trend was
  minimized by enlarging the aspect ratio of through-holes in the aperture."

**[INF] The verdict a partner needs to hear plainly.** The Samukawa lineage's published
etch demonstrations are **low-to-medium AR** (nanometre features, MEMS trenches), with
**~47 nm/min** rates and a **~5°** beam. A 5° beam at AR 200 has an acceptance cone of
0.286° — **[INF] essentially none of a 5° beam reaches the bottom of a 200:1 hole by line
of sight.** Neutral-beam etching as demonstrated to date is a *damage* technology, not an
*aspect-ratio* technology. Any 100:1+ neutral-beam claim requires a beam an order of
magnitude narrower than the surface-neutraliser state of the art.

### 4.2 The 2026 shift: gas-phase neutralisation, sub-degree beams, and who is funding it

Khrabrov & Kaganovich (PPPL), arXiv:2604.04214v2 (8 Apr 2026), is the pivot paper. Verbatim:

> "Despite their effectiveness, surface-based neutralizers inevitably introduce surface
> sputtering, material contamination, and lifetime limitations, issues that become
> increasingly severe as processing requirements push toward higher FAB energies in the KeV
> range and **extreme aspect-ratio features approaching or exceeding 100:1.**" [Q]

> "For advanced processing applications, the divergence must be quite small, with
> acceptable values typically limited to approximately **1° or less**." [Q]

> "a series of recent experimental works led by authors of Nagoya University demonstrated
> that **angular divergence well below 1° can be achieved** for accelerated ions and for the
> resulting fast neutral beam, and at sufficiently low pressure only limited by the
> transverse thermal motion of the beam ions at the source." [Q]

The AR→tolerance conversion, verbatim and directly usable [Q]:

> "the tolerance parameter for the angular divergence of the beam is set by the value of the
> aspect ratio of the etched features. **For the value of 100, it is roughly 0.5°.** In what
> follows, the tolerance will be set to that angle in the center-of-mass frame of the
> colliding pair, resulting in the angle of **0.25° in the laboratory frame.**"

And the source-design result — the maximum useful conversion efficiency [Q]:

> "For the present case, where θ* = 0.5° implies γ/α = σ*/σcx ≈ 0.3, the optimal value of
> L/λ … is close to unity (≈1.1), and the corresponding flux ratio calculated from Eq. (4)
> is 0.28. Thus, **just over one quarter of the initial ion flux can be transformed into
> fast neutrals that arrive at the target with deflection angles lying within the prescribed
> tolerance.**"

Funding, verbatim [Q]: "supported in part by **Samsung** and the U.S. Department of Energy
through the PPPL CRADA agreement with **Applied Materials**".

Kushner has entered the same space: Cardoso & Kushner, *Controlling energetic neutral beams
produced from inductively coupled plasmas for material processing applications*, 78th Annual
Gaseous Electronics Conference, October 2025 (cited as ref. [4] in the PPPL paper).

**[INF] Field-map reading.** Between April 2025 and April 2026 the neutral-beam HAR problem
went from a niche damage-mitigation topic to an actively funded frontier with **Samsung,
Applied Materials, KIOXIA, PPPL, Nagoya and Michigan all in it**. The physics they are
converging on — beam divergence tolerance set by 1/AR, tail generated by binary elastic
scattering with a realistic repulsive potential, ~25–30 % maximum useful flux conversion —
is **exactly the physics our two-component beam module already implements and gates**.

**[INF] The strategic consequence for a neutral-beam hardware partner.** A charge-free beam
removes the one term our 200:1 study lists as "not modelled" (in-feature charging) — the
study already says "A neutral-beam source is largely exempt from this term". So for a
neutral-beam tool, **our validated subset covers the dominant mechanisms**, and the 200:1
prediction is inside our envelope in a way that it is not for a CCP. That is a defensible,
non-marketing claim and it should be the centrepiece of any neutral-beam pitch.

### 4.3 Other alternative sources

- **GCIB-derived neutral beams.** US 11,199,769 (neutral beam processing from gas-cluster
  ion beam technology) claims processing of features "from 5:1 to 200:1" and that
  "ultra-high aspect ratio etch on devices with features exceeding 100:1 depth to width
  ratios have been demonstrated" [Q-relay]. **This is a patent claim with no accompanying
  peer-reviewed profile data. Do not cite as evidence; cite only as an indication that the
  claim space is occupied.** `[VERIFY]`
- **Metastable-assisted / negative-ion pulsed sources.** Samukawa's pulse-time-modulated
  negative-ion route (Appl. Surf. Sci. 253, 6681, 2007) is the mature academic path to
  high neutralisation efficiency. No HAR profile data at AR > ~20 located.
- **Pulsed-power CCP (the incumbent answer).** Lam: "Lam's pulsed power plasma technology
  utilizes increasing peak power at very short bursts drives ions much deeper with higher
  efficiency alleviating the challenges as the stack gets taller" [Q, Counterpoint 2024].
  Ishikawa 2018 notes required bias power "this year will exceed 25 kW" [Q-relay] — the
  incumbent path is buying AR with brute-force energy, which is what a narrow neutral beam
  would displace.
- **ARDE behaviour of neutral beams.** The only statement located is Miwa 2013's: etch rate
  still decreases with AR, mitigated by widening the aperture through-hole AR [Q-relay].
  **[INF] So neutral beams do not escape ARDE — they escape charging and damage.** A pitch
  that claims ARDE immunity is not supported by the literature; a pitch that claims
  charging immunity plus a sub-degree beam is.

---

## 5. The gaps — what nobody has published

Ordered by (size of the hole in the literature) × (our ability to fill it).

1. **No in-feature measurement of anything at AR > ~50.** The best in-feature transport
   measurement is Cunge 2016's capillary-plate IED, run at modest AR in a He ICP.
   Ishikawa's 2018 review says it outright: "The problem of absolute radical fluxes inside
   HAR holes has not been completely solved yet" [Q-relay]. **Nobody has measured ion
   flux, radical flux, energy spectrum, or potential inside a 100:1 feature.** Every deep
   model in the field, including ours, is unconstrained there.
   *Fillable by us + a hardware partner:* a capillary/AR-plate transmission sweep to
   AR 100–200 is a one-fixture experiment that directly gates our exact cone-acceptance
   table with zero chemistry.

2. **No published validated model at 100:1, and nothing at all at 200:1.** Ceiling is
   Kushner AR 80 (profiles) and Lam AR 140 (curves, fitted). *Fillable by us:* we are
   already the only group with a **receipted** transport calculation at 200:1. The gap is
   the evolution driver, which is exactly what phase 2 is building.

3. **No published energetic:thermal delivery ratio above AR ~40.** Kushner's 2–3× at
   AR 10–40 is the deepest datum. Our 105–144× at AR 200 is a 40–70× extrapolation with no
   competitor. *Fillable:* publish it as a prediction with the extrapolation factor named.

4. **No published statistics of twisting/tilt vs depth in the peer-reviewed literature.**
   CD-SAXS and FIB-SEM papers describe the capability and the qualitative behaviour; the
   only numeric twist values I could surface were from patent text. **A twist-vs-depth
   distribution over N holes is not in the open literature.** *Fillable jointly:* our
   engine's stochastic channels (LER modality, charging, mask morphology) are the natural
   generators of such a distribution; a partner with FIB-SEM or CD-SAXS access supplies
   the measurement.

5. **No published clog/etch-stop *statistics*.** Krüger's thesis gives a clog **boundary**
   vs O₂/C₄F₆ ratio and vs low-frequency power, validated against SEM ("The profiles,
   experiment and simulation, for P_lf = 0 kW produce total clogging of the mask opening"
   [Q]), and notes the honest limitation: "etch depth for the fully clogged feature is not
   reproduced as this depth depends on when the feature was clogged" [Q]. That is a
   boundary, not a rate. **Nobody publishes a not-open rate, a clog probability, or a
   time-to-clog distribution at AR 50–200.** *Fillable:* a clog-probability-vs-condition
   curve from a mouth-equilibrium model is a genuinely new observable and maps onto our
   existing mouth-equilibrium work.

6. **No independent reconciliation of the two conflicting deep ARDE anchors.** Huang's
   MCFPM falls 80 % by AR 40; Nguyen/Jansen's measured CORE series falls ~43 % by AR 54.
   Nobody has reconciled them. *Fillable cheaply and entirely in-repo* — and it partly
   defuses the honesty appendix's current "the magnitude of what is missing" paragraph,
   which leans on Huang alone.

7. **No measured tail fraction in a HARC-relevant reactor.** Kim 2025 measures the
   main-to-total ratio and its exponential pressure scaling in a dual-frequency Ar CCP with
   KIOXIA co-authors — **but nobody has published `f_tail` for a fluorocarbon HARC recipe
   at the wafer.** *Partly fillable:* adopt Kim's measured pressure law and quote our tail
   fraction as a measured function with a stated reactor caveat, rather than a free sweep.

8. **No self-pair-exact general body-of-revolution transport operator.** Our own §2
   limitation ("It plateaus a factor ~1.2 above tolerance and cannot be refined through")
   is, as far as I can tell, an unsolved problem nobody else has even framed, because
   everyone else uses Monte Carlo and never confronts it. *Fillable by us alone:* an
   analytic self-pair kernel would be a small, publishable numerical-methods result and it
   unblocks tapered profiles — i.e. it unblocks the entire commercially meaningful metric
   set (BB bias, taper, bow).

9. **Cryo at extreme AR has no open mechanism-to-profile chain.** The chemistry papers
   (JVST A 42, 063006, 2024; JVST A 44, 033006, 2026) are blanket-film. The profile claims
   (Lam, TEL) are press. **Nobody has published a cryo mechanism propagated through a
   feature model to a 100:1 profile.** *Fillable:* our exact transport + a Langmuir/HF
   condensation module is precisely that chain, and the cryo row in the repo already has
   the Langmuir physics.

10. **Nobody has published a *hole* (as opposed to trench) treatment of thermalised radical
    return.** Our E8 is measured only on a mask-dominated trench (2.2×). Huang's ">95 %
    above AR 10" statement is for a via. **The hole measurement does not exist anywhere.**
    *Fillable immediately* once the axisymmetric evolution driver lands.

---

## 6. Consolidated citation table

| # | citation | role here | access |
|---|---|---|---|
| 1 | Shen, Lill, Hoang, Chi, Routzahn, Church, Subramonium, Puthenkovilakam, Reddy, Bhadauriya, Roberts, Kamarthy (Lam Research) — *Progress report on high aspect ratio patterning for memory devices* — **Jpn. J. Appl. Phys. 62, SI0801 (2023)**, DOI 10.35848/1347-4065/accbc7 | the single best field-map source; 2.5 % @ AR 50, flat-to-AR-60 ARDE, selectivity to AR 140, BB-bias numbers, twisting mechanism | OPEN ACCESS; local extract `research_sources/thesis_extracts/lam_shen_lill_jjap2023.txt` |
| 2 | Huang, Huard, Shim, Nam, Song, Lu, Kushner (Michigan + Samsung) — *Plasma etching of high aspect ratio features in SiO₂ using Ar/C₄F₈/O₂ mixtures* — **J. Vac. Sci. Technol. A 37, 031304 (2019)**, DOI 10.1116/1.5090606 | deepest fully-instrumented deep-AR model (AR 40 base, AR 80 profiles); flux/energy/charging numbers | open PDF `cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf` (cert invalid, use `curl -k`) |
| 3 | Khrabrov & Kaganovich (PPPL) — *Ion-neutral and neutral-neutral scattering in argon at KeV energies and implications for high-aspect-ratio etching* — **arXiv:2604.04214v2** (8 Apr 2026) | first-principles origin of the measured beam tail; AR→tolerance law; neutral-beam source optimum; Samsung/AMAT funded | open |
| 4 | Nguyen, Shkondin, Jensen, Hübner, Leussink, Jansen — *Ultrahigh aspect ratio etching of silicon in SF₆–O₂ plasma: the CORE sequence and chromium mask* — **J. Vac. Sci. Technol. A 38, 053002 (2020)**, DOI 10.1116/6.0000357 | **deepest usable measured depth-vs-time series (gap AR 57)** | open mirror at DTU labadviser |
| 5 | Okawa, …, Takahashi (Tohoku / SRIS / NanoTerasu) — *Nondestructive structural evaluation of high-aspect-ratio-etched holes by tender x-ray ptychography with computed tomography* — **Appl. Phys. Lett. 129, 011109 (2026)** | newest and best in-feature 3D geometry at ~10 µm; 38.5 nm 2D / 67.5 nm 3D | paywalled; press release at nanoterasu.jp |
| 6 | Kim, Kawamura, Fujitani, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda (Nagoya) — *Influence of sheath collisions on the ion angular distributions in a dual-frequency capacitively coupled plasma* — **Jpn. J. Appl. Phys. 64, 096002 (2025)**, DOI 10.35848/1347-4065/ae0105 | **measured main/tail ratio and its pressure law — the route to retiring our tail-fraction knob** | paywalled |
| 7 | Kawamura, Kim, Ichikawa, Suzuki, Iino, Fukumizu, Kurihara, Toyoda (Nagoya + **KIOXIA**) — *Angular distribution measurement of high-energy ions and neutrals impinging on an rf electrode…* — **Plasma Sources Sci. Technol. 34, 055006 (2025)**, DOI 10.1088/1361-6595/add321 | ion *and* fast-neutral angular distributions to ~2.5 keV, memory-maker co-authored | paywalled |
| 8 | Kim, Kawamura, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda — *Measurement of energy-resolved ion angular distribution…* — **Jpn. J. Appl. Phys. 64, 05SP15 (2025)** | the energy-resolved IAD our beam module is gated against | paywalled |
| 9 | Ichikawa, Chu, Moriyama, Nakahara, Suzuki, Iino, Fukumizu, Kurihara, Toyoda — **Appl. Phys. Express 14, 126001 (2021)** | first of the Nagoya series; shows the tail is not explained by Phelps' model | paywalled |
| 10 | Cunge, Darnon, Dubois, Bezard, Mourey, Petit-Etienne, Vallier, Despiau-Pujo, Sadeghi — *Measuring ion velocity distribution functions through high-aspect ratio holes in inductively coupled plasmas* — **Appl. Phys. Lett. 108, 093109 (2016)** | the only in-feature transport measurement method; our transport gate template | paywalled; method described in Impedans app note VE01 (open) |
| 11 | Gahan, Scullin, Hopkins (Impedans) — *Ion energy and ion flux measurements through high-aspect ratio holes using the Vertex system*, App. Note VE01 | AR-resolved IED without venting; commercial availability of the gate experiment | open PDF |
| 12 | Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (Tokyo Electron Miyagi) — *Precise and practical 3D topography simulation of high aspect ratio contact hole etch by using model optimization algorithm* — **J. Vac. Sci. Technol. A 42, 043003 (2024)**, DOI 10.1116/6.0003515 | rare experimental *horizontal* cross-sections of HARC holes; twisting/distortion target | paywalled |
| 13 | Zhang, Lee, Klochkov, Korb, Sorkhabi, Lan, Pichumani, Tekleyohannes, Wang, Sallis, Ningen, Kim, Teoh, Polubotko, Pirkle, Foca — *3D reconstruction of 3D NAND memory etch profiles using FIB-SEM* — **Proc. SPIE 13426, 134262J (2025)**, DOI 10.1117/12.3051156 | cryo failure-mode list; hole-pattern variance | paywalled |
| 14 | *CDSAXS study of 3D NAND channel hole etch pattern edge effects and etched hole pattern variance* — **Proc. SPIE 12955, 1295539 (2024)**, DOI 10.1117/12.3010927; and *Inline metrology of high aspect ratio hole tilt using SAXS* — **Proc. SPIE 12053, 1205312 (2022)** | statistical CD/tilt/distortion vs depth; inner vs outer holes | paywalled |
| 15 | Ishikawa, Karahashi, Ishijima, Cho, Elliott, Hausmann, Mocuta, Wilson, Kinoshita — *Progress in nanoscale dry processes for fabrication of high-aspect-ratio features: how can we control critical dimension uniformity at the bottom?* — **Jpn. J. Appl. Phys. 57, 06JA01 (2018)**, DOI 10.7567/JJAP.57.06JA01 | the field's own statement of the unsolved problems; MaCE AR 160/500 comparison points | open |
| 16 | Lill, Wang, Wu, Oh, Kim, Wilcoxson, Singh, Ghodsi, George, Barsukov, Kaganovich — *Low-temperature etching of silicon oxide and silicon nitride with hydrogen fluoride* — **J. Vac. Sci. Technol. A 42, 063006 (2024)** | cryo HF chemistry, blanket rates vs T | OSTI open (`osti.gov/servlets/purl/2514386`) |
| 17 | *Reaction mechanism for HF based cryogenic plasma etching of SiO₂* — **J. Vac. Sci. Technol. A 44, 033006 (2026)** | the cryo mechanism to build a module from; −10 to −100 °C window; H₂O-catalysed | paywalled |
| 18 | Krüger — *Modeling and Optimization of High Aspect Ratio Plasma Etching* — PhD thesis, University of Michigan (2024); journal instance **J. Vac. Sci. Technol. A 42, 043008 (2024)**, DOI 10.1116/6.0003554 | the clog boundary validated against SEM; our own calibration source; hf = 825 nm | open PDF; local extract `krueger_thesis.txt` |
| 19 | Choi, Kim, Jeong, Cho, Lee, Seong, Choi, You — *Computational study of the evolution of high-aspect-ratio SiO₂ trench and hole features during physical sputtering* — **Phys. Plasmas 33, 033901 (2026)** | K-SPEED; trench-vs-hole difference at HAR; GPU feature-scale | paywalled |
| 20 | Vanraes, Venugopalan, Besemer, Bogaerts — *Assessing neutral transport mechanisms in aspect ratio dependent etching by means of experiments and multiscale plasma modeling* — **Plasma Sources Sci. Technol. 32, 065003 (2023)**, DOI 10.1088/1361-6595/acdc4f | experiment+model ARDE with charging explicitly excluded; inverse RIE lag near AR 1 | paywalled |
| 21 | Huard, Zhang, Sriraman, Paterson, Kanarik, Kushner — *Role of neutral transport in aspect ratio dependent plasma etching of three-dimensional features* — **J. Vac. Sci. Technol. A 35, 05C301 (2017)** | neutral-saturated / ion-starved regime alleviates ARDE — the regime-map reference for our limiting-regime defect | paywalled |
| 22 | Kuboi — *Review and perspective of dry etching and deposition process modeling of Si and Si dielectric films for advanced CMOS device applications* — **Jpn. J. Appl. Phys. 63, 080801 (2024)**, DOI 10.35848/1347-4065/ad5355 | the field's modeling-SOTA review; K-SPEED GPU ~200× claim | open access |
| 23 | Counterpoint Research (commissioned by Lam) — *Scaling to 1,000-Layer 3D NAND in the AI Era*, July 2024 | vendor-level roadmap, cryo numbers, −60 °C, 10 µm, 2.5×, >100:1 with multiple tiers | open PDF (Lam mediaroom) |
| 24 | Lam Research newsroom — *Lam Cryo™ 3.0*, 31 July 2024 | the "<0.1 % CD deviation at 10 µm" headline | open (investor page 403s to automated fetch; newsroom blog works) |
| 25 | Tokyo Electron — *Memory Channel Hole Etch Technology … 10-µm-deep Etching for 3D NAND Flash with Over 400 Layers*, 9 June 2023; underlying paper at **2023 Symp. VLSI Technology & Circuits** | 10 µm / 33 min / −70 °C; the only competing absolute-rate anchor | open press; VLSI paper `[VERIFY]` not retrieved |
| 26 | Cardoso & Kushner — *Controlling energetic neutral beams produced from inductively coupled plasmas for material processing applications* — **78th GEC**, Oct 2025 | Kushner entering neutral beams | abstract only |
| 27 | Kubota, Nukaga, Ueki, Sugiyama, Inamoto, Ohtake, Samukawa — **J. Vac. Sci. Technol. A 28, 1169 (2010)**; Miwa, Nishimori, Ueki, Sugiyama, Kubota, Samukawa — **J. Vac. Sci. Technol. B 31, 051207 (2013)** | the neutral-beam comparison set: 1 mA/cm², >99 % neutralisation, 47 nm/min, ARDE still present | paywalled |
| 28 | Mohanty, Patil, Maheshwari, Mohapatra — *Ionic Reflection-Induced Profile Evolution in 3D NAND Etching* (SEMI ASMC 2026) and *Impact of Surface Charging on Feature Profile Evolution in 3D NAND Pillar Etching* (SISPAD 2026, accepted) | the newest academic entrants on exactly our problem — watch these | not yet published |

---

## 7. Concrete follow-ups this map generates (ranked)

1. **Extract Kim 2025's measured main-to-total ratio vs pressure** (JJAP 64, 096002) and
   re-quote our tail fraction as a measured function. Retires the single declared-open
   parameter of the 200:1 study. *Blocked on: paywalled PDF.*
2. **Chase the Okawa/Takahashi APL 2026 supplementary for a CD-vs-depth trace.** If it
   exists it is the first scoreable in-feature hole profile at AR ≈ 100 in the open
   literature.
3. **Digitise Lam Shen 2023 Fig. 6(b)** (normalised ion flux vs AR for several IAD 1σ).
   Direct, like-for-like comparison against our delivery-vs-AR table at the field's own
   operating point. Open-access source; cheap.
4. **Reconcile Huang's 80 %-by-AR-40 against Nguyen/Jansen's ~43 %-by-AR-54.** Cheap,
   entirely in-repo, and it repairs the weakest paragraph of the current honesty appendix.
5. **Re-word the `HOLE_STUDY_RESULTS` §5 flat-delivery gate** so it does not, read
   literally, disqualify Lam's published operating point (§3.4).
6. **Add "average cascaded-particle energy vs AR" to the phase-1 output.** Kushner reports
   ~400 eV maintained above AR 3; we have the number and do not print it. One line, one
   new literature comparison.
7. **Adopt BB bias (largest bow CD − bottom CD) as a reported metric** alongside top/neck/
   z_neck. It is the industry's key verticality metric and it is relative, so it dodges the
   +29 % absolute-depth defect.
8. **Scope the analytic self-pair kernel** for the general body-of-revolution operator. It
   is the gate on tapered profiles, hence on every commercially meaningful metric, and it
   appears to be a problem nobody else has framed.
9. **For the neutral-beam partner deck:** lead with the PPPL 2026 tolerance law
   (AR 100 ⇒ 0.5° COM / 0.25° lab; ~28 % maximum useful flux conversion), state plainly
   that surface-neutraliser beams are ~5° and therefore not an AR technology, and note that
   a charge-free source is the one configuration where our validated subset covers the
   dominant mechanisms at 200:1.

---

## 8. Addendum — reconciling against phase-2, which landed while this was being written

`HOLE_STUDY_RESULTS_PHASE2_2026-08-06.md` (uncommitted, produced in parallel) now reports a
coupled **rate**-vs-AR curve: rate(200)/rate(1) = 0.3082 at tail 0.65, i.e. **69 % fall by
AR 200**, and 60 % by AR 50. That changes the field comparison materially, so:

**Three-way triangulation of deep-AR rate decline — the field now has exactly three points
and we are between the other two.**

| source | decline | at AR | class |
|---|---|---|---|
| Huang MCFPM (Kushner/Samsung), fluorocarbon SiO₂ via | **80 %** | 40 | S — simulation, HPEM-bounded |
| **petch phase 2**, Deck 1, tail 0.65 | **60 %** | 50 | S — simulation, exact transport, zero knobs |
| Nguyen/Jansen CORE, silicon SF₆/O₂ pillar array | **≈43 %** [INF from Fig. 8 series] | 54 | **E — measured** |
| Lam Shen 2023, ONON cryo | flat to ~60, then declining | 60→140 | S — fitted to experiment |

**[INF] Reading.** petch's coupled decline sits *between* the deepest simulation and the
deepest measurement, and *above* the vendor's own low-temperature curve. That is a much
stronger position than the phase-1 appendix's framing ("our transport-only 10 % vs Huang's
80 % — the gap is everything that converts delivery into rate"). The honest statement now
is: **there are only three deep-AR rate-decline numbers in existence, they span 43–80 %,
and ours is in the middle of them.** Nobody in the field agrees with anybody else to better
than a factor of ~2 at AR 40–60, and no measurement exists past AR 57.

**Two caveats before this triangulation is used externally.**
1. Chemistry and geometry differ across all four rows (fluorocarbon oxide via vs Si pillar
   gap vs ONON cryo hole). This is an order-of-magnitude triangulation, not a like-for-like
   comparison, and must be labelled as such.
2. The Nguyen/Jansen point is my derivation from a depth-vs-time series, not a number the
   authors print, and its AR axis assumes the Fig. 8 array is the 200 nm/400 nm one —
   `[VERIFY]` before quoting.

**Also now unlocked by phase 2:** gap #10 of §5 (nobody has published a *hole* treatment of
thermalised radical return). Phase 2 reports E8 immaterial in the hole (fourth decimal
across the full physical band). **That is a publishable negative result** — it is a
literature-relevant claim about Huang's ">95 % above AR 10" statement, made in the geometry
that statement was made for, and it exists nowhere else.
