# Validating feature charging at deep aspect ratio (AR 20–200)

Literature research, 2026-07-29. Scope question: **what published, quantitative observables can
falsify petch's charging module in the AR 20–200 regime**, and in what order should we gate them.

Companion docs: `CHARGING_PHYSICS_PLAN.md` (W1/W2/W3 mechanism plan, AR≤4 gates),
`DEBOER_DIVERGENCE.md` (charging-throttle hypothesis REFUTED for the de Boer Si floor),
`EXPERIMENTAL_VALIDATION_MATRIX.md` (evidence tiers), `FLOOR_OVERCHARGE_FINDING.md`.

Provenance convention used below:
- **[FULL]** = I read the full text (PDF fetched and parsed); quotes are verbatim.
- **[ABS]** = publisher/Crossref/OSTI abstract fetched; quotes are verbatim from that abstract.
- **[META]** = title/author/venue/DOI confirmed from Crossref or a publisher page; content not read.
- **[VERIFY]** = claim not confirmed from a fetched source; flagged for follow-up.

---

## 0. Headline findings (the five things that change what we do)

1. **There is a hard, published ceiling on the in-feature potential: `V_max < E_ion,max / e`, and it
   is only weakly AR-dependent above AR ≈ 17.** Huang & Kushner (JVST A **44**, 023013, 2026) [FULL]
   report, for preformed SiO₂ trenches at AR = 16.7 / 25 / 50 in an Ar CCP, maximum in-feature
   potentials of 108/112/112 V (max ion energy 125 eV), 211/225/227 V (250 eV) and 315/330/344 V
   (380 eV). Ratio `V_max·e/E_ion,max` = 0.83–0.91, **rising monotonically with AR but never
   reaching 1**. Verbatim: *"The maximum positive potential for all AR is lower than the maximum ion
   energy of 125 eV."* This is the single cheapest and hardest falsifier available for petch's
   deep-AR over-charge failure mode (`FLOOR_OVERCHARGE_FINDING.md`), because it is a self-limiting
   *physical bound*, not a fitted curve. Nine (AR, V₀) points, fully specified reactor conditions,
   **preformed (non-evolving) geometry** ⇒ no chemistry coupling to confound it.

2. **Charging alone does not produce the de Boer-style deep floor collapse — literature agrees with
   our own STEP-2 refutation.** In the best-documented deep dielectric case (Huang *et al.*, JVST A
   **37**, 031304, 2019 [FULL], AR 40, tri-frequency CCP, Ar/C₄F₈/O₂), charging reduces mean ion
   energy at the etch front from 1940 → 1050 eV (−46 %) and lengthens the etch by **only ~30 %**
   (36 min → 48 min, *"which are in reasonable agreement with experiments"* quoting 40–50 min HVM
   times for AR 40–50). Etch *stop* appears only in the low-power arm (2.5 kW, before AR = 40).
   Separately, Ohiwa *et al.* (JJAP **37**, 5060, 1998) [ABS] attribute experimental HARC etch stop
   to **redeposition of sputtered fluorocarbon**, not charging: *"the redeposition of sputtered
   species from the fluorocarbon polymer on the hole sidewall induces the etch stop at the bottom of
   the high-aspect hole"*, and note ions *"are able to bombard the bottom of the hole maintaining
   their high energy"* despite sidewall neutralisation. **Consequence: do not gate charging on the
   de Boer SF₆/O₂ cryo-Si floor** (conductive substrate, thin SiO_xF_y passivation). Gate it on
   dielectric HARC and on notching.

3. **The strongest *statistical* deep-AR observable is twisting incidence, and it has a published
   charging-on/off ablation.** Wang & Kushner (JAP **107**, 023309, 2010) [FULL], AR = 20 SiO₂
   trench (75 nm opening, 1500 nm deep), 41 random seeds per arm: twisting incidence **12 %
   (5/41) without charging → 49 % (20/41) with charging**; polymer conductivity 0.01 Ω⁻¹cm⁻¹ →
   38 %; conductive SiO₂ + insulating polymer → 25 %; both conductive → **12 %, i.e. back to the
   no-charging baseline**; and a monotone `V_dc` sweep 44 % (V_dc = 0) → 10 % (V_dc = −1000 V).
   Reported in-trench maximum potential 100–150 V, *"occurs roughly half-way down the trench"*.
   This is a five-point ablation ladder on an ensemble statistic — far more falsifying than any
   single-profile match.

4. **Direct experimental in-feature potential measurements exist, but only up to AR ≈ 6–7 in
   published detail.** The two instrument families are (a) Samukawa/Ohtake/Jinnai **on-wafer
   monitoring chips** with top and bottom electrodes inside real SiO₂ contact structures — Shimmura
   *et al.* JVST A **22**, 433 (2004) [ABS, AR = 5.7], Ohtake *et al.* JVST A **24**, 2172 (2006)
   [ABS, pulsed vs CW], Ohtake *et al.* JVST B **25**, 400 (2007) [ABS, *IEDF and electron energy at
   the hole bottom*], Jinnai *et al.* JVST B **25**, 1808 (2007) [ABS, *"the charge accumulation
   potential between the top and bottom of the contact-hole structures increased with the aspect
   ratio of the contact holes"*]; and (b) Kamata & Arimoto's **current-through-the-dielectric**
   measurements, JAP **80**, 2637 (1996) and JVST B **14**, 3688 (1996) [ABS], which give
   *"The dc self-bias potential difference reached about 100 V with a substrate rf bias voltage of
   400 V for an aspect ratio of 2"*. A newer route is Park & Chung's **anodic-alumina (AAO) template
   probe**, Rev. Sci. Instrum. **89** (2018) [ABS], which is the only device architecture that can
   in principle reach AR ≫ 20 (self-organised parallel pores) — but the published paper is a device
   demonstration, not an AR sweep. **There is no published direct measurement of V(AR) at AR > ~20.**

5. **Sub-degree coupling is the dominant deep-AR failure channel, and the threshold lateral
   potential is of order 1 V.** Derived below (§4): a lateral potential asymmetry `ΔV_lat` across
   the feature width, acting coherently over the feature depth, throws an ion of axial energy
   `qV_z` out of the AR-set acceptance cone when `ΔV_lat ≳ 2 V_z / AR²`. At AR 50 / 2 kV that is
   1.6 V; at AR 100 / 5 kV, 1.0 V; at AR 200 / 10 kV, **0.5 V** — i.e. a **0.2 % azimuthal asymmetry
   of a 250 V in-feature potential is already lethal.** Meanwhile the *intrinsic* IADF width is
   already at the same scale: Kushner's group (JVST A **43**, 033001, 2025) [FULL] gives the
   collisionless-sheath estimate `θ ≈ tan⁻¹[(k_B T_I / qV_S)^{1/2}]`, *"For an ion approaching the
   sheath with a parallel velocity corresponding to 0.1 eV with a sheath potential of 1000 V, the
   corresponding angle at the wafer is about 0.6°"* — versus acceptance half-angles
   arctan(1/AR) = 2.86° (AR 20), 1.15° (AR 50), 0.573° (AR 100), **0.286° (AR 200)**.

---

## 1. Experimental observables with published quantitative data

### 1.1 Charging potential vs aspect ratio — direct electrical measurement

| Source | Geometry / AR | Measured quantity | Digitizable? |
|---|---|---|---|
| **Kamata & Arimoto**, *Charge build-up in Si-processing plasma caused by electron shading effect*, **J. Appl. Phys. 80, 2637–2642 (1996)**, DOI 10.1063/1.363179 [ABS] | Si substrate + 500 nm SiO₂ line-and-space, varied pattern size; Ar ICP 2–40 mTorr | **Electron and ion currents through the dielectric structure**; floating potential vs pattern size; floating-potential difference vs T_e (2 → 4 eV) | **Yes** — figures are current/potential vs pattern size and vs pressure. This is the *only* published measurement that separates the electron and ion current channels. |
| **Kamata & Arimoto**, *Suppression of electron shading effect by a counter radio frequency bias in plasma etching*, **JVST B 14, 3688–3691 (1996)**, DOI 10.1116/1.588648 [ABS] | Hole patterns, AR ≥ 2, vs substrate rf bias to 400 V | dc self-bias potential *difference* between HAR hole pattern and open area. Verbatim: *"The dc self-bias potential difference reached about 100 V with a substrate rf bias voltage of 400 V for an aspect ratio of 2."* Also: difference *"increased, independent of the hole pattern aspect ratio, at the lower substrate rf bias voltage and they tended to saturate"* at higher bias | **Yes** — V(bias) saturation curves. The *saturation* is itself the ceiling physics of §0.1 measured experimentally. |
| **Shimmura, Suzuki, Soda, Samukawa, Koyanagi, Hane**, *Mitigation of accumulated electric charge by deposited fluorocarbon film during SiO₂ etching*, **JVST A 22, 433–436 (2004)**, DOI 10.1116/1.1649347 [ABS via OSTI] | SiO₂ contact holes, **AR = 5.7**, on-wafer monitoring chip | *"The dc potential of the SiO₂ contact hole top and bottom surfaces were measured during plasma exposure with/without deposited fluorocarbon film in the holes."* | **Yes** — top/bottom dc potential time traces, ± polymer. Directly tests the sidewall-conductivity channel. |
| **Ohtake, Jinnai, Suzuki, Soda, Shimmura, Samukawa**, *Real-time monitoring of charge accumulation during pulse-time-modulated plasma*, **JVST A 24, 2172–2175 (2006)**, DOI 10.1116/1.2362724 [ABS] | SiO₂ contact structure, on-wafer chip | Top-vs-bottom accumulation potential, CW rf-bias vs pulse-TM; time-resolved electron and ion flows. Verbatim: *"the accumulated charge in the pulse-time-modulated operation was drastically decreased"* | **Yes** — CW/pulsed potential traces. This is the pulsed-rescue observable in its cleanest electrical form. |
| **Ohtake, Jinnai, Suzuki, Soda, Shimmura, Samukawa**, *On-wafer monitoring of electron and ion energy distribution at the bottom of contact hole*, **JVST B 25, 400–403 (2007)**, DOI 10.1116/1.2712200 [ABS] | SiO₂ contact hole, Ar UHF plasma | **IEDF at the hole bottom**, plus *"a lower electron density and higher electron temperature at the contact-hole bottom due to the electron-shading effect, as compared with that in the bulk plasma"*; *"The peak energy of IEDF corresponded to the sheath potential."* | **Yes** — this is the single most direct in-feature validation target in the literature: petch predicts exactly these two distributions. |
| **Jinnai, Orita, Konishi, Hashimoto, Ichihashi, Nishitani, Kadomura, Ohtake**, *On-wafer monitoring of charge accumulation and sidewall conductivity in high-aspect-ratio contact holes during SiO₂ etching process*, **JVST B 25, 1808–1813 (2007)**, DOI 10.1116/1.2794050 [ABS] | 8-in wafer, HARC structures *"comparable with the practical interconnect structures of recent DRAM devices"* | Verbatim: *"the charge accumulation potential between the top and bottom of the contact-hole structures increased with the aspect ratio of the contact holes"* — i.e. **the experimental V(AR) curve** | **Yes, and it is the target curve.** Need the figure; AR range not stated in the abstract. **[VERIFY]** the exact AR span. |
| **Ohmori, Goto, Kitajima, Makabe**, APL **83**, 4637–4639 (2003), DOI 10.1063/1.1630163 [ABS]; and **Ohmori & Makabe**, *In situ measurement of plasma charging on SiO₂ hole bottoms and reduction by negative charge injection during etching*, **Appl. Surf. Sci. 254, 3696–3709 (2008)**, DOI 10.1016/j.apsusc.2007.10.070 [META] | SiO₂ contact hole, 2f-CCP CF₄/Ar | *"A reduction in charging voltage is measured in the pulsed operation both of the plasma power source and of the wafer bias"* | **Yes** — the 2008 ASS paper is a 14-page in-situ measurement study; the highest-value un-fetched source in this list. **[VERIFY]** contents. |
| **Park & Chung**, *A monitoring device made of an anodic aluminum oxide template for plasma-induced charging potential measurements in the high-aspect-ratio trench structure*, **Rev. Sci. Instrum. 89 (2018)**, DOI 10.1063/1.5042017 [ABS] | AAO template, self-organised parallel pores, top+bottom floating electrodes | Potential difference top–bottom vs pressure, power, gas species | Device paper; **no published AR sweep**. Architecturally the only path to AR ≫ 20 direct measurement. |

### 1.2 Charging *damage* vs aspect ratio (indirect but the original evidence)

- **Hashimoto**, *Charge Damage Caused by Electron Shading Effect*, **JJAP 33, 6013 (1994)**,
  DOI 10.1143/jjap.33.6013 [ABS] (and JJAP **32**, 6109, 1993). Antenna-MOS capacitors under
  photoresist with HAR openings. Verbatim: *"This damage increased with the pattern's aspect
  ratio, and occurred even when the test wafer was cut into chips about 5 mm square and mounted on
  a wafer with insulation."* And: *"The damaging current from this mechanism increased by a factor
  of more than ten with a decrease in the gate oxide thickness only from 8 nm to 6 nm, implying that
  the degree of shading depends on the gate charging voltage."* Digitizable: breakdown/damage
  fraction vs AR; F–N current vs oxide thickness. **This is the founding electron-shading dataset.**
- **Arita, Akamatsu, Asano**, *Reduction of Charge Build-Up during RIE by Using SOI Structures*,
  **JJAP 36, 1505 (1997)**, DOI 10.1143/jjap.36.1505 [ABS]: MNOS capacitors on SOI; charge build-up
  decreases as buried-oxide thickness increases; etch-rate penalty only 3 %. Digitizable: build-up
  vs BOX thickness — a clean *capacitive-divider* test of the charge-transport boundary condition.
- **Dostalik, Krishnan, Kinoshita, Rangan**, *Electron shading effects in high density plasma
  processing for very high aspect ratio structures*, 3rd Int. Symp. Plasma Process-Induced Damage
  (1998) 160–163, DOI 10.1109/ppid.1998.725599 [META]. **[VERIFY]** — title promises exactly the
  deep-AR damage scaling we want.

### 1.3 Notching depth vs aspect ratio, and the pulsed rescue

- **Fujiwara, Maruyama, Yoneda**, *Profile Control of poly-Si Etching in Electron Cyclotron
  Resonance Plasma*, **JJAP 34, 2095 (1995)**, DOI 10.1143/jjap.34.2095 [ABS]. Local side etch
  (notch) vs outside-space width, vs perpendicular electron temperature `T_ev`, vs ion current
  density, vs rf bias. Verbatim: *"The local side etch with rf bias depends on electron temperature.
  Lower T_ev is effective for r[educing]…"*. Hwang & Giapis (JAP 82, 566) [FULL] cite this as
  *"a monotonic but complex increase in notch depth … with an increase of the AR in the regime
  between 0.7:1 and 2.8:1"* — **that AR range is the digitizable experimental notch-vs-AR curve.**
- **Fujiwara, Maruyama, Yoneda**, *Pulsed Plasma Processing for Reduction of Profile Distortion
  Induced by Charge Buildup in ECR Plasma*, **JJAP 35, 2450 (1996)**, DOI 10.1143/jjap.35.2450
  [ABS]. **The single best rescue dataset.** Verbatim: *"Notch depth reduction was observed in both
  gases with repeated pulsing of the plasma. This reduction strongly depends on the off-time
  length."*; *"In the case of HCl plasma, no aspect-ratio dependence of notching was observed. This
  indicates that accumulated charges on every pattern area were neutralized by exposure to repeated
  pulses of plasma."* Digitizable: notch depth vs off-time, Cl₂ vs HCl; and the **AR-dependence
  collapse** is a structural (sign) test, not a magnitude fit.
- **Maruyama, Fujiwara, Ogino, Yoneda**, JJAP **36**, 2526 (1997) [ABS]; **Maruyama, Fujiwara,
  Ogino, Miyatake**, JJAP **37**, 2306 (1998) [ABS]: *"pulse rf bias in pulse plasma is very
  efficient for reducing the notch and charge build-up. In contrast, continuous rf bias in pulse
  plasma is not effective"* — a two-arm discriminator on *which* modulation matters.
- **Nozawa, Kinoshita, Nishizuka, Narai, Inoue, Nakaue**, *The Electron Charging Effects of Plasma
  on Notch Profile Defects*, **JJAP 34, 2107 (1995)**, DOI 10.1143/jjap.34.2107 [ABS]. Verbatim:
  *"the notch depth increases as the 'perimeter ratio', (i.e. the ratio of the pad perimeter to the
  notch line perimeter), increases"*; and notch occurs *"only outside the line"* for non-connected
  lines but *"at all of the conn[ected lines]"* when connected. **This is a topology/conductor-network
  test** — it validates the conductor-terminal boundary condition in `conductor_terminal_3d.py`
  without touching any energy calibration. Already partially referenced by
  `scripts/nozawa_1995_charge_checkpoint_audit.py`.
- **Nishioka & Fujiwara**, JJAP **34**, 5998 (1995) [ABS]: dc-biased surface potential control;
  *"The profile is improved by setting the applied dc potential equal to floating potential"* —
  a knob-free zero-crossing test.
- **Ohtake & Samukawa** lineage on charge-free etching: Samukawa & Mieno, PSST **5**, 132 (1996)
  [META]; Samukawa, JVST B **12**, 3300 (1994) [META]. Negative-ion afterglow injection.

### 1.4 Twisting statistics

- **Wang & Kushner**, JAP **107**, 023309 (2010) [FULL] — simulation ensemble, but with an explicit
  ablation ladder (see §0.3). Also states the physical framing: *"the statistical variation in
  charged particle fluxes into the feature could charge one side of the feature more than the other
  … produce asymmetrical electric fields that deflect subsequent ions from the vertical"*, and
  *"once a feature begins to twist, it tends to be self-perpetuating"*.
- **Huang, Shim, Nam, Kushner**, *Pattern dependent profile distortion during plasma etching of high
  aspect ratio features in SiO₂*, **JVST A 38, 023001 (2020)**, DOI 10.1116/1.5132800 [FULL]. AR 40
  SiO₂ (AR 53 including mask). Centroid statistics at the stop layer over 10–20 seeds vs pitch:
  *"The centroids of features with patterns below 100 nm vary by as much as 10 nm compared to only a
  few nanometers for the sparsest patterns"*; maximum in-feature potential *"≈1000 V"* at 100 nm
  pitch vs *"≈300 V"* at 200 nm pitch, with a **40 % etch-rate penalty for the dense pattern**;
  and *"The trend of more statistical distortion and random tilting with dense patterns agrees with
  experiments"* (citing Negishi *et al.* 2017). Charge-mobility ablation: raising polymer mobility
  to 50 cm²/V·s (≥10⁹× bulk PTFE) reduces max potential only 585 → 501 V.
- **Negishi, Miyake, Yokogawa, Oyama, Kanekiyo, Izawa**, *Bottom profile degradation mechanism in
  high aspect ratio feature etching based on pattern transfer observation*, **JVST B 35, 051205
  (2017)**, DOI 10.1116/1.4998943 [META] — the experimental counterpart Kushner compares to; also
  the source of the *"200 V at AR = 20"* in-feature potential quoted by Huang 2019. **[VERIFY]**
  whether that 200 V is measured or modelled in Negishi 2017.
- **Miyake, Negishi, Izawa, Yokogawa, Oyama, Kanekiyo**, *Effects of Mask and Necking Deformation on
  Bowing and Twisting in High-Aspect-Ratio Contact Hole Etching*, **JJAP 48, 08HE01 (2009)**,
  DOI 10.1143/jjap.48.08he01 [META] — experimental; competing (non-charging) explanation.
- **Metrology that makes twisting a *quantitative* observable**: Gin, Wormington, Amasay, Grinberg,
  Brady, Reichental, Matney, Zhang, *Inline metrology of high aspect ratio hole tilt and center line
  shift using small-angle x-ray scattering*, **J. Micro/Nanopatterning Mater. Metrol. 22, 031205
  (2023)**, DOI 10.1117/1.jmm.22.3.031205 [META] (and SPIE 12053-58, 2022, DOI 10.1117/12.2614312).
  CD-SAXS gives *ensemble-averaged* tilt and centre-line shift — the natural pairing for a
  simulation ensemble statistic. **[VERIFY]** magnitudes/uncertainties.
- **Shen, Lill, Hoang, Chi, Routzahn, Church, Subramonium, Puthenkovilakam, Reddy, Bhadauriya,
  Roberts, Kamarthy** (Lam Research), *Progress report on high aspect ratio patterning for memory
  devices*, **JJAP 62, SI0801 (2023)**, DOI 10.35848/1347-4065/accbc7 [FULL via publisher summary].
  Industry numbers: channel holes *>100:1* at ~100 nm CD; *"At 50:1, only 2.5 % of incoming neutral
  flux reaches the cylinder terminus"*; *"Narrow IAD (<0.5° 1σ) correlates with flatter etch rates"*;
  selectivity falls ~50 % from AR 40 → 140; twisting *"observed at depths ≥4000 nm"*; digitizable
  Figs. 3, 6, 7 (normalised etch rate / ion flux / selectivity vs AR).

### 1.5 Etch stop vs aspect ratio

- **Ohiwa, Kojima, Sekine, Sakai, Yonemoto, Watanabe**, *Mechanism of Etch Stop in High
  Aspect-Ratio Contact Hole Etching*, **JJAP 37, 5060 (1998)**, DOI 10.1143/jjap.37.5060 [ABS].
  Experimental. *"etch stop occurs at higher aspect ratios for the same hole diameter in oxide films
  with higher boron and phosphorous dopant concentrations"* — i.e. **etch-stop AR is a measurable
  function of film chemistry**, digitizable as AR_stop vs dopant level. Their mechanism is
  redeposition, and ions retain energy at the bottom **despite** sidewall neutralisation. This is
  the strongest published counter-evidence to a charging-dominated etch stop.
- **Matsui, Nakano, Petrović, Makabe**, *The effect of topographical local charging on the etching of
  deep-submicron structures in SiO₂ as a function of aspect ratio*, **Appl. Phys. Lett. 78, 883–885
  (2001)**, DOI 10.1063/1.1347021 [ABS]. Simulation. Verbatim: *"When the aspect ratio is greater
  than seven, the bottom is charged up to a potential sufficient to prevent the influence of all the
  incident ions, with a realistic initial energy of 300 eV for SiO₂ etching within the period
  required for monolayer stripping, resulting in etch stop."* **This AR > 7 etch-stop prediction is
  refuted by manufacturing reality (AR > 100 holes are etched daily), and by Huang & Kushner 2026's
  ceiling result.** It is the canonical example of a charging model that over-charges at depth —
  exactly petch's failure mode. Treat it as a *negative* control: petch must NOT reproduce it.
- **Huang, Huard, Shim, Nam, Song, Lu, Kushner**, JVST A **37**, 031304 (2019) [FULL]: charging-induced
  etch stop predicted *only* in the 2.5 kW low-power arm before AR = 40; at 10 kW the etch completes.
  Verbatim: *"the decrease in ion energy produces an etch stop at 2.5 kW before reaching AR = 40."*
  Digitizable Fig. 20: bottom ion energy vs 5 MHz power, with/without charging (1400→2460 eV without,
  475→1760 eV with; ΔE/E₀ from 66 % to 28 %).

### 1.6 Bosch / cryo DRIE charging evidence

Weak and mostly indirect. What exists:
- SOI **notching / footing** at the buried-oxide interface is the accepted DRIE charging signature;
  Arita *et al.* JJAP **36**, 1505 (1997) [ABS] quantifies the *charge* side (MNOS capacitors, BOX
  thickness dependence) but not the notch geometry.
- Maruyama, Narukage, Onuki, Fujiwara, *High-aspect-ratio deep Si etching in SF₆/O₂ plasma. II*,
  **JVST B 28, 862 (2010)**, DOI 10.1116/1.3466884 [ABS] — explicitly finds the deep-Si lateral etch
  *"seems to depend on time (or etched depth) rather than aspect ratio"* and attributes it to radical
  distribution, **not charging**. Another data point against charging as the deep-Si mechanism.
- Lai, Johnson, Westerman, JVST A **24**, 1283 (2006) [ABS] — Bosch/TDM ARDE lag model; no charging.
- **Conclusion:** the Bosch/cryo literature does **not** supply a quantitative charging observable at
  deep AR. Any petch claim tying the de Boer AR > 20 collapse to charging is unsupported and is
  already refuted internally (`DEBOER_DIVERGENCE.md` STEP 2, and STEP 3 which reassigns the collapse
  to Monte-Carlo deep-floor under-sampling).

---

## 2. How the Kushner-school MCFPM validates charging

### 2.1 The algorithm of record

Krueger, *Modeling and Optimization of High Aspect Ratio Plasma Etching*, PhD thesis, Univ. of
Michigan (2024) [FULL], §2.2.2, describes MCFPM charging verbatim:

> "Electrostatic charging of features results from ions or neutrals which neutralize upon striking
> surfaces and deposit their charge at the impact site. Charge is then retained on the voxel upon
> which it is deposited until neutralized by an opposing charge or transported via conductive charge
> transport."

with `dρ_k/dt = Σ_i q_i w_i / V_k − ∇·(μ_k ρ_k E)`, Poisson solved by red-black SOR on a finite
volume mesh, per-material permittivity and positive/negative mobilities, reflective lateral BCs,
zero-gradient top, grounded bottom, potential updated only *"after a large amount of charged
particles have deposited their charge"*, and *"All positive ions neutralize upon their first
collision with a surface and return as a hot neutral"*. Primary reference: **Wang & Kushner, JAP
107, 023309 (2010)**.

Structural comparison to petch: petch's `charging_poisson_3d.py` / `charging_coupled_3d.py` solve the
same quasi-static Laplace/Poisson problem, but petch's electron supply is analytic (RF-phase weighted)
while MCFPM's electron EADs come from the reactor-scale HPEM/eMCS. **The hot-neutral conversion on
first ion impact is a first-class part of MCFPM's charging model and is what keeps deep floors
etching once the potential saturates** — this is the mechanism that petch must have to reproduce
finding §0.1/§0.2.

### 2.2 What each paper actually claims as agreement

| Paper | Agreement claimed | Strength |
|---|---|---|
| Wang & Kushner, JAP **107**, 023309 (2010) [FULL] | *"These results are consistent with experimental reports of twisting in Si etched in fluorocarbon gas mixtures, particularly if the overlying polymer is thick."* Twisting incidence 12 %→49 % (charging), →12 % (all-conductive). | **Qualitative** vs experiment; **quantitative and reproducible** as an internal ablation ladder. |
| Huang *et al.*, JVST A **37**, 031304 (2019) [FULL] | (i) Total etch time 36 min (no charging) / 48 min (charging) vs *"as long as 40–50 min"* in HVM for AR 40–50 — *"in reasonable agreement with experiments"*. (ii) In-feature potentials 200–400 V at AR 10–20, *"commensurate with"* prior 200 V @ AR 20 (Negishi 2017), 300 V @ AR 10 (Matsui/Maeshige/Makabe, JPD **34**, 2950, 2001), 150 V @ AR 15 (Radjenović, Radmilović-Radjenović, Petrović, IEEE TPS **36**, 874, 2008). (iii) Sidewall polymer 0–20 nm predicted vs 20–100 nm and 5–20 nm measured at AR 10–20. | **Weak** (integral time; order-of-magnitude potential; polymer thickness *under*-predicted). No point-by-point profile fit. |
| Huang *et al.*, JVST A **38**, 023001 (2020) [FULL] | *"The trend of more statistical distortion and random tilting with dense patterns agrees with experiments"* (Negishi 2017). | **Trend only.** |
| Krüger, Lee, Nam, Kushner, JVST A **41**, 013006 (2023), DOI 10.1116/6.0002290 [FULL, downloaded] | VWT controls in-feature charging and defect formation. | Computational; no experimental charging comparison. |
| Huang & Kushner, JVST A **44**, 023013 (2026), DOI 10.1116/6.0005187 [FULL] | *"These results at large ARs agree with past findings"* (refs 5, 17, 18). SEE reduces max potential 16 %/13 %/2 % at AR 16.7/25/50 (V₀ = 450 V), 5 %/3 %/1 % at 300 V, ~0 at 150 V. Concludes *"anisotropic electrons become less able to remediate charging as the potential in the feature is dominated by charge redistribution."* | **Internal-consistency + literature cross-check.** No experiment. |

**Honest read:** the Kushner school has never validated charging against a *measured* in-feature
potential. Its charging validation is (a) order-of-magnitude potential consistency across independent
codes, (b) integral etch-time consistency with HVM, (c) qualitative trend agreement on twisting
statistics. **That is the bar petch has to clear — and the ceiling law (§0.1) plus the twisting
ablation ladder (§0.3) let us clear it more rigorously than the incumbent does.**

Note also for W2 of `CHARGING_PHYSICS_PLAN.md`: the 2026 SEE paper's finding that SEE's effect
**shrinks with AR** (16 % → 2 % from AR 16.7 → 50) is the opposite of what the plan's Memos-based
50 % reduction assumption would extrapolate to. W2's expected effect size should be re-preregistered
as AR-dependent and *small at deep AR*.

### 2.3 The Korean / GPU-accelerated charging Monte Carlo — status and correction

The parent brief assumes a "Korean CUDA-accelerated charging MC". What I can verify:

- **Yook, You, Park, Chang, Kwon, Yoon, Yoon, Shin, Yu, Im**, *Fast and realistic 3D feature profile
  simulation platform for plasma etching process*, **J. Phys. D: Appl. Phys. 55, 255202 (2022)**,
  DOI 10.1088/1361-6463/ac58cf [ABS via publisher]. Jeonbuk National University + KRISS + Samsung-adjacent.
  GPU-parallelised neutral and ion transport, hash-map 3-D level set, adaptive surface meshing,
  two-layer surface reaction model. Verbatim: *"the speedup ratio, as compared to a single central
  processing unit (CPU), is approximately 200"*. Applied to HAR oxide etching and MOSFET SiN spacer.
  Predecessor code is referred to as **K-SPEED**; earlier GPU work presented at GEC 2013
  (*GPU based 3D feature profile simulation of high-aspect ratio contact hole etch process under
  fluorocarbon plasmas*, Bull. APS 2013, GEC MR1.020) and ICOPS 2016 (DOI 10.1109/plasma.2016.7534064).
  **I could NOT verify that this platform contains a charging/Poisson module.** The published
  abstracts describe transport + surface reactions + level set only. **[VERIFY] — this is the single
  most important open citation question in this document.** If it has no charging module, the claim
  "the Korean CUDA charging MC" should be retired from our memory notes.
- The genuine recent **charging** codes with published validation are elsewhere:
  - **Zhai, Ge, Hu, Li, Shao, Cheng, Filipović, Chen**, *Modeling the charging effect of the hardmask
    and silicon substrate during plasma etching in advanced nodes*, **J. Appl. Phys. 137, 063302
    (2025)**, DOI 10.1063/5.0243470 [META, DOI verified via Crossref 2026-08-01; publisher page returns 403]. Chinese
    Academy of Sciences + TU Wien. Validated against **Si etching experiments in Cl₂ plasma with SEM
    profile comparison**; identifies charging contributions to mask faceting, substrate bowing,
    microtrenching. **[VERIFY]** the DOI and the exact validation figures.
  - Same group, *Modeling the Charging Effect on the Twisting Defects During High Aspect Ratio
    Etching of Dielectrics*, IEEE conf. (2025), IEEE Xplore doc 11017954 [META].
  - **Cheong, Lee, Kim, Kim, Whang** (Seoul National Univ.), PSST **23**, 065051 (2014) [META] — this
    is **RIE lag in a magnetized ICP**, *not* charging. Do not cite it as a charging code.
  - **Memos & Kokkoris** lineage (NCSR Demokritos): Micromachines **9**, 415 (2018), DOI
    10.3390/mi9080415 [ABS] — open access, SEEE + ion reflection + charging on PMMA; already the
    W2 anchor in `CHARGING_PHYSICS_PLAN.md`.
  - **Krüger, Wilczek, Mussenbrock, Schulze**, *Voltage waveform tailoring in radio frequency plasmas
    for surface charge neutralization inside etch trenches*, **PSST 28, 075017 (2019)**,
    DOI 10.1088/1361-6595/ab2c72 [ABS] — PIC/MCC + trench charging; the European lineage.
  - **The Implementation of the Surface Charging Effects in Three-Dimensional Simulations of SiO₂
    Etching Profile Evolution**, SCIRP (open access) — already cited in `CHARGING_PHYSICS_PLAN.md`
    as the PIC(XPDC1)-derived boundary-input implementation.

---

## 3. The self-limiting potential law (the load-bearing new physics for petch)

Collecting Huang & Kushner 2026 [FULL], Ar CCP, 10 mTorr, 20 MHz, preformed trenches:

| V₀ (V) | max ion energy (eV) | V_max, AR 16.7 | AR 25 | AR 50 | ratio V_max·e/E_max (16.7 / 25 / 50) |
|---:|---:|---:|---:|---:|---|
| 150 | 125 | 108 | 112 | 112 | 0.86 / 0.90 / 0.90 |
| 300 | 250 | 211 | 225 | 227 | 0.84 / 0.90 / 0.91 |
| 450 | 380 | 315 | 330 | 344 | 0.83 / 0.87 / 0.91 |

Three structural facts, all falsifiable:

1. **Ceiling.** `e·V_max < E_ion,max` at every point. Physically: once the local potential equals the
   ion energy, ions stop arriving there and charging stops. Ceiling is *not* an assumption — it is
   the fixed point of the flux balance.
2. **Saturation with AR.** Going AR 16.7 → 50 (a 3× increase) moves `V_max` by only +4 % (150 V),
   +8 % (300 V), +9 % (450 V). Any model whose `V_max` keeps climbing steeply past AR ≈ 25 is wrong.
   Petch's observed deep-AR over-charge is exactly this failure.
3. **Location migration.** *"The maximum potential for the smaller AR remains at the bottom of the
   feature. The location of the maximum potential for an AR of 50 shifts downward in the feature
   with increasing V₀"*, and at AR 50 with V₀ = 150 V, ions *"first strike and neutralize on the
   sidewalls before reaching the bottom of the feature, producing a maximum in plasma potential of
   112 V at height that corresponds to an AR of approximately 25."* Consistent with Huang 2019
   [FULL]: *"By an AR of about 10, the majority of ions will have collided with the sidewalls where
   they deposit charge, then proceeding as hot neutrals which are not affected by the electric
   fields"*, and Wang & Kushner 2010 [FULL]: max potential *"occurs roughly half-way down the
   trench"*. **The maximum in-feature potential does NOT sit at the etch front at deep AR.** Huang
   2019: *"When the etch front reaches the bottom stop layer, the maximum electric potential higher
   in the feature (AR ≈ 10) is 1100 V, which is about 60 % of the average energy of incident ions.
   Had this maximum been at the bottom of the feature, the ion energy incident onto the etch front
   would have decreased by about 60 %."* — i.e. the *reason* deep features keep etching is that ions
   are decelerated *then re-accelerated* past the mid-feature maximum, plus the hot-neutral channel.

**Direct implication for petch:** any charging solver that pins the potential maximum at the floor
will over-suppress the deep-AR floor flux. `FLOOR_OVERCHARGE_FINDING.md` should be re-read in this
light; the fix is likely the ion-neutralisation/hot-neutral conversion and the sidewall-deposition
bookkeeping, not the electron integrator alone.

---

## 4. Sub-degree coupling: charging deflection vs the acceptance angle

### 4.1 Geometry

Acceptance half-angle for an unobstructed traverse of a feature of width `W`, depth `D`, `AR = D/W`:

| AR | 20 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| `arctan(1/AR)` | 2.86° | 1.15° | 0.573° | 0.382° | **0.286°** |

(centre-launched ions see half of this). The parent's "~0.3° at AR 200" is `arctan(1/200)`. ✔

### 4.2 Intrinsic ion angular width already sits at this scale

- Collisionless sheath (Kushner group, JVST A **43**, 033001, 2025 [FULL]):
  `θ ≈ tan⁻¹[(k_B T_I / q V_S)^{1/2}]`; at `T_I` = 0.1 eV and `V_S` = 1000 V → **0.6°**.
  Reducing to 0.286° requires `V_S ≥ ~4 kV` at the same `T_I`, or `T_I ≤ 0.025 eV` at 1 kV.
- Industry target (Lam, JJAP **62**, SI0801, 2023 [FULL]): *"Narrow IAD (<0.5° 1σ)"*.
- Ar CCP at 300 V amplitude (Huang & Kushner 2026 [FULL]): *"the angular breadth is ±1.5° and the
  average energy is 211 eV"*; at 450 V, *"±1.3°"*.
- dc-augmented high-energy electron beams (Wang & Kushner 2010 [FULL]): *"energies up to a few keV
  with angular spreads of <0.5°"* — note the *electrons* can be made more directional than the ions.

**So the intrinsic IADF alone already exceeds the AR-200 acceptance cone by ~2×, before charging.**
Deep-AR etching therefore cannot proceed by direct line-of-sight ion delivery; it proceeds by
grazing sidewall reflection → hot neutrals (Kushner) plus glancing-angle scattering (Kim *et al.*,
JVST A **40**, 053007, 2022 [META]). **A charging model validated only on axial ion delivery is
validating the wrong channel at AR > 100.**

### 4.3 The charging deflection threshold (derived here — not a literature quote)

Take an ion of axial energy `E_z = qV_z` traversing a depth `L_f` over which a lateral field
`E_x = ΔV_lat / W` acts coherently (`ΔV_lat` = lateral potential difference across the feature width):

    Δv_x / v_z = q E_x L_f / (m v_z²) = (ΔV_lat / W)·L_f / (2 V_z)

so, with `L_f = αD` and `AR = D/W`,

    θ_defl ≈ α · (ΔV_lat / 2V_z) · AR       [radians]

Setting `θ_defl = θ_accept ≈ 1/AR` gives the **critical lateral asymmetry**

    ΔV_lat,crit ≈ 2 V_z / (α · AR²)

| AR | V_z assumed | α = 1 (fully coherent) | α = 0.1 (localised patch) |
|---:|---:|---:|---:|
| 20 | 1 kV | 5.0 V | 50 V |
| 50 | 2 kV | 1.6 V | 16 V |
| 100 | 5 kV | 1.0 V | 10 V |
| 200 | 10 kV | **0.5 V** | 5 V |

Because `θ_defl ∝ AR` while `θ_accept ∝ 1/AR`, the ratio grows as **AR²**. Even in the pessimistic
localised case (α = 0.1), a 5 V azimuthal asymmetry kills a 10 keV ion at AR 200 — against measured/
modelled in-feature potentials of **100–1000 V** (§1, §3). **Deep-AR charging is not a magnitude
problem, it is a symmetry problem: the model must get the *azimuthal residual* of the potential right
to ~0.5 % to predict twisting at AR 200.** This is why every credible deep-AR study (Wang &
Kushner 2010; Huang 2020) treats twisting as a *stochastic ensemble* quantity, not a deterministic
profile — and why petch's `charging_coevolution_3d.py` ensembles are the right instrument.

Corroborating statement (Huang 2020 [FULL]): *"the number of particles/second entering the feature
has decreased to the point that the type (e.g., radical versus ion), energy, and angle of reactants
is subject to statistical noise"*; Wang & Kushner 2010 [FULL]: for a 50 nm via at 10¹⁶ cm⁻²s⁻¹,
*"The time between the arrival of two ions … is 5 μs"*. Discrete-charge shot noise is therefore the
*source* of `ΔV_lat` — see also Petrović/Radjenović, Phys. Plasmas **14**, 103501 (2007),
*Electrostatic potential fluctuation induced by charge discreteness in a nanoscale trench* [META].

### 4.4 Electron shading theory (Hashimoto) in one line

Hashimoto's model (JJAP **32**, 6109, 1993; **33**, 6013, 1994) [ABS]: photoresist with HAR openings
geometrically shades the underlying conductor from *obliquely incident* electrons while transmitting
the *normally incident* ion flux; local flux imbalance charges the bottom positive **without any
wafer-scale potential difference** (proved by the cut-chip-on-insulated-wafer control). The
non-linearity — damage current up >10× for an 8 → 6 nm oxide — implies the shading fraction itself
depends on the accumulated gate voltage, i.e. **the shading is self-consistent, not geometric**.
Kamata & Arimoto (1996) [ABS] then measured the two currents separately and showed the floating
potential difference *grows with T_e* (2 → 4 eV) — the same T_e lever Hwang & Giapis modelled
(PRL **79**, 845, 1997; JAP **81**, 3433, 1997 [ABS]: *"Larger values of T_e cause the potential of
the upper photoresist sidewalls to become more negative; thus, more electrons are repelled back and
the electron current density to the trench bottom decreases"*).

---

## 5. Preregistered validation campaign (5 gates, ordered)

Ordering principle: **cheapest and most falsifying first; geometry-only before chemistry-coupled;
bounds before curves; ensembles before single profiles.** Every gate below is preregistered — the
pass band is fixed *now*, before running.

---

### GATE D1 — Potential ceiling and AR saturation *(preformed geometry, no chemistry)*

**Reference:** Huang & Kushner, JVST A **44**, 023013 (2026), DOI 10.1116/6.0005187, Figs. 10–14 and
the quoted values in §3 above. Conditions fully specified in the paper: Ar, 10 mTorr, 20 MHz CCP,
V₀ = 150/300/450 V, electrodes r = 15 cm, gap 4.0 cm; trenches AR 16.7/25/50; via AR 20, d = 100 nm,
depth 2000 nm.

**Petch configuration:** `charging_coupled_3d` / `charging_poisson_3d` on a *preformed* (frozen)
trench — no level-set evolution, no chemistry. Ion IEAD taken as the published bimodal Ar⁺
distribution (max energy 125/250/380 eV; angular breadth ±1.5° at 300 V, ±1.3° at 450 V), electrons
thermal at T_e from the paper. **SEE off** for the primary arm.

**Pass criteria (all must hold):**
- **D1a (hard bound):** `e·V_max ≤ E_ion,max` for all nine (AR, V₀) points. *A single violation fails
  the gate outright* — this is a physical bound, not a fit.
- **D1b (magnitude):** `V_max` within **±20 %** of the nine published values.
- **D1c (saturation):** `V_max(AR 50) / V_max(AR 16.7)` in **[1.00, 1.20]** at each V₀ (published:
  1.04, 1.08, 1.09). A petch value > 1.5 means the deep-AR over-charge is still live.
- **D1d (location):** at AR 50, the potential maximum lies at a depth corresponding to AR 20–30, not
  at the floor.

**Why first:** no chemistry, no evolution, no digitization, nine points, one hard physical bound, and
it directly targets the known petch failure. Estimated cost: one preformed-geometry sweep.

**Kill criterion:** if D1a fails after the W1 integrator fix, the defect is in the ion boundary
condition (missing neutralise-on-first-impact / hot-neutral conversion), not in the electron model.
Implement hot-neutral conversion before anything else.

---

### GATE D2 — Electron-shading current and floating-potential saturation *(experimental)*

**References:** Kamata & Arimoto, JAP **80**, 2637 (1996), DOI 10.1063/1.363179 (electron/ion
currents through the dielectric structure vs pattern size; floating-potential difference vs T_e
2 → 4 eV); and JVST B **14**, 3688 (1996), DOI 10.1116/1.588648 (dc self-bias potential difference
vs substrate rf bias to 400 V, AR = 2; saturation behaviour; counter-rf-bias suppression).

**Observables:** (i) ratio `I_e/I_i` reaching the substrate through the patterned dielectric, vs
pattern size; (ii) `ΔV_float` vs T_e; (iii) `ΔV_float` vs rf bias amplitude with the ~100 V value at
400 V / AR 2 and the *saturation* at higher bias.

**Pass criteria:**
- **D2a:** predicted `ΔV_float(400 V bias, AR 2)` = **100 V ± 40 %**.
- **D2b (structural):** the bias sweep must **saturate** — d(ΔV)/d(V_bias) at 400 V less than half
  its value at 150 V. Sign/shape test, not a magnitude fit.
- **D2c (structural):** `∂ΔV_float/∂T_e > 0` over 2 → 4 eV, with magnitude within 2× of the digitized
  slope.

**Why second:** it is *measured*, it is cheap (AR 2), and it validates the **source model** (the
electron/ion arrival distributions of W3) which everything at deep AR inherits. D2b is the
experimental shadow of D1a.

**Digitization:** figures must be pulled from the two papers; both are paywalled AIP/AVS.
**[VERIFY]** access route.

---

### GATE D3 — Notch-vs-AR and the pulsed rescue *(experimental, mechanism-level)*

**References:** Fujiwara, Maruyama, Yoneda, JJAP **34**, 2095 (1995), DOI 10.1143/jjap.34.2095
(notch depth vs outside-space width / AR 0.7:1–2.8:1, vs T_ev, vs ion current density);
JJAP **35**, 2450 (1996), DOI 10.1143/jjap.35.2450 (notch depth vs pulse off-time, Cl₂ vs HCl; **AR
dependence of notching disappears in pulsed HCl**); Maruyama *et al.*, JJAP **37**, 2306 (1998)
(pulsed-rf-bias-in-pulsed-plasma effective, CW-bias-in-pulsed-plasma ineffective).
Cross-reference Hwang & Giapis, JAP **82**, 566 (1997), Fig. 8 (predicted notch depth vs AR;
*"in good agreement with experimental trends observed in the work of Fujiwara et al."*).

**Pass criteria:**
- **D3a:** monotone increase of notch depth with AR over 0.7:1 → 2.8:1, with the *rate* of deepening
  increasing above AR ≈ 2 (Hwang & Giapis: *"the rate of notch deepening increases significantly for
  an AR > 2"*). Shape test on normalised depth; ±30 % on relative depths.
- **D3b (rescue):** with the negative-carrier / afterglow channel enabled, notch depth must fall
  monotonically with off-time, **and the AR dependence must flatten to within noise** (the HCl
  result). This is a two-dimensional structural test with no free magnitude.
- **D3c (discriminator):** pulsed rf bias in pulsed plasma reduces notch; CW rf bias in pulsed plasma
  does not. Sign test.
- **D3d (topology, cheap add-on):** Nozawa *et al.*, JJAP **34**, 2107 (1995): notch depth increases
  with pad-perimeter ratio; connected lines notch everywhere, non-connected lines notch only on the
  outer line. Pure conductor-network test of `conductor_terminal_3d.py`; binary pass/fail.

**Why third:** notching is the observable petch already demonstrates (Fujiwara-monotone, HG shape
r = 0.92, commit c07886c). D3 converts that from a single trend into a four-way structural gate and
retires the borrowed HG `E_defl(AR)` table per `CHARGING_PHYSICS_PLAN.md`'s finish line.

---

### GATE D4 — Twisting-incidence ensemble and its ablation ladder

**Reference:** Wang & Kushner, JAP **107**, 023309 (2010), DOI 10.1063/1.3290873 [FULL] (DOI verified via Crossref 2026-08-01; Part I is 10.1063/1.3290870).
Geometry: SiO₂ trench, 75 nm mask opening, 1500 nm to Si stop layer, **AR = 20**, dome-shaped PR
450 nm max thickness; Ar/C₄F₈/O₂ = 80/15/5, 40 mTorr, 300 sccm, 10 MHz rf, 4 kW; **41 random seeds
per arm**; thermal electrons T_e ≈ 2 eV.

**Published ladder (the pass targets):**

| Arm | Twisting incidence |
|---|---|
| No charging (fluorocarbon/SiO₂) | 5/41 = **12 %** |
| Charging, insulating polymer + insulating SiO₂ | 20/41 = **49 %** |
| Charging, polymer σ = 0.01 Ω⁻¹cm⁻¹ | **38 %** |
| Charging, insulating polymer + conductive SiO₂ | **25 %** |
| Charging, both conductive | **≈12 %** (baseline) |
| Charging + HEE flux, V_dc = 0 | 18/41 = **44 %** |
| Charging + HEE flux, V_dc = −750 V | 7/41 = **17 %** |
| V_dc sweep 0 → −1000 V | **44 % → 10 %**, monotone |

**Pass criteria:**
- **D4a (ordering):** petch must reproduce the **strict ordering** of all five conductivity arms.
  Ordering violations fail.
- **D4b (magnitude):** charging-on incidence within **±15 percentage points** of 49 %; all-conductive
  arm within ±10 pp of the no-charging arm.
- **D4c (monotonicity):** incidence decreases monotonically with |V_dc| over 0 → −1000 V.
- **D4d (potential):** in-trench maximum potential 100–150 V at AR 20, located near mid-depth
  (±25 % of depth). Cross-checks D1d at a different AR and chemistry.
- **Ensemble discipline:** ≥ 41 seeds per arm; report a binomial CI, not a point estimate.

**Why fourth:** it is the only published *statistical* charging observable with a clean ablation
ladder, and it is the direct test of the §4.3 azimuthal-symmetry requirement — the actual deep-AR
physics. It is also the most expensive (≈ 200+ 3-D evolving runs).

---

### GATE D5 — Deep-AR ion-energy budget and the *negative* etch-stop control

**References:** Huang *et al.*, JVST A **37**, 031304 (2019), DOI 10.1116/1.5090606, Figs. 19–20
[FULL]. AR 40 SiO₂ (AR 53 from PR top), tri-frequency CCP Ar/C₄F₈/O₂.

**Published targets:**
- Mean ion energy at the etch front, AR 40: **1940 eV (no charging) → 1050 eV (charging)**, −46 %.
- 5 MHz power sweep 2.5 → 10 kW: bottom ion energy 1400 → 2460 eV (no charging), 475 → 1760 eV
  (charging); `ΔE/E₀` = 66 % → 28 %.
- Total etch time: 36 min (no charging) → 48 min (charging), **+33 %**, vs 40–50 min HVM.
- Hot-neutral flux to the etch front is **nearly independent of charging**; hot-neutral mean energy
  falls 50–150 eV; hot-neutral power to the front falls ~15 %.
- Etch stop occurs **only** at 2.5 kW, before AR 40.

**Pass criteria:**
- **D5a:** `ΔE/E₀` at AR 40 = **46 % ± 15 pp**.
- **D5b (channel separation):** hot-neutral flux to the etch front changes by **< 10 %** between
  charging-on and charging-off. If petch's charging suppresses the neutral/hot-neutral channel, the
  neutralisation bookkeeping is wrong.
- **D5c (etch-time budget):** charging lengthens the etch by **20–50 %**, not by 2×+.
- **D5d (NEGATIVE control — must FAIL to reproduce):** petch must **not** reproduce Matsui/Nakano/
  Petrović/Makabe's AR > 7, 300 eV etch stop (APL **78**, 883, 2001). At 300 eV incident and AR ≥ 7
  in an oxide trench, petch must still deliver non-zero power to the etch front. Reproducing that
  etch stop is a **failure**, because the manufacturing record (AR > 100 holes) refutes it.
- **D5e (scope discipline):** petch must **not** claim charging as the mechanism for the de Boer
  SF₆/O₂ cryo-Si AR > 20 floor collapse. `DEBOER_DIVERGENCE.md` STEP 2 refuted it internally;
  Maruyama *et al.* JVST B **28**, 862 (2010) and Ohiwa *et al.* JJAP **37**, 5060 (1998) refute it
  externally. Any charging-throttle arm run against de Boer must be labelled a negative result.

---

### Ordering rationale and stopping rules

```
D1 (bound, preformed, 9 pts)  ──► D2 (measured source model, AR 2)
        │                                    │
        └──────────► D3 (notch structure + rescue, AR 0.7–2.8) ──► D4 (twisting ensemble, AR 20)
                                                                        │
                                                                        └──► D5 (deep budget, AR 40)
```

- **If D1a fails:** stop. Fix ion neutralisation / hot-neutral conversion. Everything downstream is
  confounded.
- **If D1 passes but D2b fails:** the source (arrival-distribution) model is the defect ⇒ execute W3
  of `CHARGING_PHYSICS_PLAN.md` (1-D RF sheath MC) before D3.
- **If D3a/D3d pass but D3b fails:** the negative-carrier channel is missing, not the electrostatics.
- **If D4a ordering fails but D4b magnitude passes:** the conductivity/charge-transport model is
  wrong (mobilities, polymer, substrate), not the field solver.
- **If D5b fails:** the neutral channel is being incorrectly coupled to the field — the most likely
  root cause of a spurious deep-AR etch stop in petch.

### What none of these gates can do

No published dataset gives **measured in-feature potential above AR ≈ 20**. The campaign above
validates the *mechanism* at AR 2–50 and validates *bounds and structure* at AR 40–50; extrapolation
to AR 100–200 remains model-based. The honest paper statement is: petch's charging is validated
against measurement to AR ≈ 6 (on-wafer monitors), against measured structure to AR ≈ 3 (notching)
and AR ≈ 2 (shading currents), and against the best available simulation consensus + physical bounds
to AR 50. **Closing the AR > 20 measurement gap would require an AAO-template probe (Park & Chung,
RSI 89, 2018) AR sweep, which nobody has published — that is a genuine experimental opportunity for
a design partner with a fab.**

---

## 6. Complete source list

**[FULL] — full text fetched and read**
1. G. S. Hwang and K. P. Giapis, *Aspect-ratio-dependent charging in high-density plasmas*,
   J. Appl. Phys. **82**, 566–571 (1997). DOI 10.1063/1.365616.
   Open PDF: https://authors.library.caltech.edu/records/je8bd-j6v68/files/HWAjap97b.pdf
2. M. Wang and M. J. Kushner, *High energy electron fluxes in dc-augmented capacitively coupled
   plasmas. II. Effects on twisting in high aspect ratio etching of dielectrics*,
   J. Appl. Phys. **107**, 023309 (2010). PDF: https://cpseg.eecs.umich.edu/pub/articles/jap_107_023309_2010.pdf
3. S. Huang, C. Huard, S. Shim, S. K. Nam, I.-C. Song, S. Lu, M. J. Kushner, *Plasma etching of high
   aspect ratio features in SiO₂ using Ar/C₄F₈/O₂ mixtures: A computational investigation*,
   J. Vac. Sci. Technol. A **37**, 031304 (2019). PDF: https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
4. S. Huang, S. Shim, S. K. Nam, M. J. Kushner, *Pattern dependent profile distortion during plasma
   etching of high aspect ratio features in SiO₂*, J. Vac. Sci. Technol. A **38**, 023001 (2020).
   DOI 10.1116/1.5132800. PDF: https://cpseg.eecs.umich.edu/pub/articles/JVSTA_38_023001_2020.pdf
5. C. Huang and M. J. Kushner, *Consequences of secondary electron emission on charging of SiO₂
   features in capacitively coupled plasmas having sinusoidal and tailored bias waveforms*,
   J. Vac. Sci. Technol. A **44**, 023013 (2026). DOI 10.1116/6.0005187.
   PDF: https://cpseg.eecs.umich.edu/pub/articles/JVSTA_44_023013_2026.pdf
6. M. J. Kushner group, *Consequences of low bias frequencies in inductively coupled plasmas on ion
   angular distributions for high aspect ratio plasma etching*, J. Vac. Sci. Technol. A **43**,
   033001 (2025). DOI 10.1116/6.0004250. PDF: https://cpseg.eecs.umich.edu/pub/articles/JVSTA_43_033001_2025.pdf
7. F. Krüger, *Modeling and Optimization of High Aspect Ratio Plasma Etching*, PhD thesis, University
   of Michigan (2024). https://cpseg.eecs.umich.edu/pub/theses/Krueger_Florian_PhD_Thesis_2024.pdf
8. F. Krüger, H. Lee, S. K. Nam, M. J. Kushner, *Voltage waveform tailoring for high aspect ratio
   plasma etching of SiO₂ using Ar/CF₄/O₂ mixtures: Consequences of ion and electron distributions
   on etch profiles*, J. Vac. Sci. Technol. A **41**, 013006 (2023). DOI 10.1116/6.0002290.
9. M. Shen *et al.* (Lam Research), *Progress report on high aspect ratio patterning for memory
   devices*, Jpn. J. Appl. Phys. **62**, SI0801 (2023). DOI 10.35848/1347-4065/accbc7.
10. M. J. Kushner group, *Future of plasma etching for microelectronics: Challenges and
    opportunities*, J. Vac. Sci. Technol. B **42**, 041501 (2024).
11. M. J. Kushner group, *Voltage waveform tailoring … Consequences of low fundamental frequency
    biases*, Phys. Plasmas **31**, 033508 (2024).

**[ABS] — abstract fetched verbatim (Crossref JATS or OSTI)**
12. K. Hashimoto, *Charge Damage Caused by Electron Shading Effect*, Jpn. J. Appl. Phys. **33**,
    6013 (1994). DOI 10.1143/jjap.33.6013. (Also JJAP **32**, 6109, 1993.)
13. T. Kamata and H. Arimoto, *Charge build-up in Si-processing plasma caused by electron shading
    effect*, J. Appl. Phys. **80**, 2637–2642 (1996). DOI 10.1063/1.363179.
14. T. Kamata and H. Arimoto, *Suppression of electron shading effect by a counter radio frequency
    bias in plasma etching*, J. Vac. Sci. Technol. B **14**, 3688–3691 (1996). DOI 10.1116/1.588648.
15. N. Fujiwara, T. Maruyama, M. Yoneda, *Profile Control of poly-Si Etching in Electron Cyclotron
    Resonance Plasma*, Jpn. J. Appl. Phys. **34**, 2095 (1995). DOI 10.1143/jjap.34.2095.
16. N. Fujiwara, T. Maruyama, M. Yoneda, *Pulsed Plasma Processing for Reduction of Profile
    Distortion Induced by Charge Buildup in ECR Plasma*, Jpn. J. Appl. Phys. **35**, 2450 (1996).
    DOI 10.1143/jjap.35.2450.
17. T. Maruyama, N. Fujiwara, S. Ogino, M. Yoneda, *Reduction of Charge Build-up with High-Power
    Pulsed ECR Plasma*, Jpn. J. Appl. Phys. **36**, 2526 (1997). DOI 10.1143/jjap.36.2526.
18. T. Maruyama, N. Fujiwara, S. Ogino, H. Miyatake, *Reduction of Charge Build-Up with
    Pulse-Modulated Bias in Pulsed ECR Plasma*, Jpn. J. Appl. Phys. **37**, 2306 (1998).
    DOI 10.1143/jjap.37.2306.
19. K. Nishioka and N. Fujiwara, *Effect of Electric Field on ECR Plasma Etching*, Jpn. J. Appl.
    Phys. **34**, 5998 (1995). DOI 10.1143/jjap.34.5998.
20. T. Nozawa, T. Kinoshita, T. Nishizuka, A. Narai, T. Inoue, A. Nakaue, *The Electron Charging
    Effects of Plasma on Notch Profile Defects*, Jpn. J. Appl. Phys. **34**, 2107 (1995).
    DOI 10.1143/jjap.34.2107.
21. T. Shimmura, Y. Suzuki, S. Soda, S. Samukawa, M. Koyanagi, K. Hane, *Mitigation of accumulated
    electric charge by deposited fluorocarbon film during SiO₂ etching*, J. Vac. Sci. Technol. A
    **22**, 433–436 (2004). DOI 10.1116/1.1649347.
22. H. Ohtake, B. Jinnai, Y. Suzuki, S. Soda, T. Shimmura, S. Samukawa, *Real-time monitoring of
    charge accumulation during pulse-time-modulated plasma*, J. Vac. Sci. Technol. A **24**,
    2172–2175 (2006). DOI 10.1116/1.2362724.
23. H. Ohtake *et al.*, *On-wafer monitoring of electron and ion energy distribution at the bottom of
    contact hole*, J. Vac. Sci. Technol. B **25**, 400–403 (2007). DOI 10.1116/1.2712200.
24. B. Jinnai, T. Orita, M. Konishi, J. Hashimoto, Y. Ichihashi, A. Nishitani, S. Kadomura,
    H. Ohtake, *On-wafer monitoring of charge accumulation and sidewall conductivity in
    high-aspect-ratio contact holes during SiO₂ etching process*, J. Vac. Sci. Technol. B **25**,
    1808–1813 (2007). DOI 10.1116/1.2794050.
25. T. Ohmori, T. K. Goto, T. Kitajima, T. Makabe, *Negative charge injection to a positively charged
    SiO₂ hole exposed to plasma etching in a pulsed two-frequency capacitively coupled plasma in
    CF₄/Ar*, Appl. Phys. Lett. **83**, 4637–4639 (2003). DOI 10.1063/1.1630163.
26. J. Matsui, N. Nakano, Z. Lj. Petrović, T. Makabe, *The effect of topographical local charging on
    the etching of deep-submicron structures in SiO₂ as a function of aspect ratio*, Appl. Phys.
    Lett. **78**, 883–885 (2001). DOI 10.1063/1.1347021.
27. T. Ohiwa, A. Kojima, M. Sekine, I. Sakai, S. Yonemoto, Y. Watanabe, *Mechanism of Etch Stop in
    High Aspect-Ratio Contact Hole Etching*, Jpn. J. Appl. Phys. **37**, 5060 (1998).
    DOI 10.1143/jjap.37.5060.
28. J. C. Arnold and H. H. Sawin, *Charging of pattern features during plasma etching*, J. Appl.
    Phys. **70**, 5314–5317 (1991). DOI 10.1063/1.350241.
29. K. Arita, M. Akamatsu, T. Asano, *Reduction of Charge Build-Up during Reactive Ion Etching by
    Using Silicon-On-Insulator Structures*, Jpn. J. Appl. Phys. **36**, 1505 (1997).
    DOI 10.1143/jjap.36.1505.
30. J.-H. Park and C.-W. Chung, *A monitoring device made of an anodic aluminum oxide template for
    plasma-induced charging potential measurements in the high-aspect-ratio trench structure*,
    Rev. Sci. Instrum. **89**, 113503 (2018). DOI 10.1063/1.5042017.
31. G. S. Hwang and K. P. Giapis, *The influence of electron temperature on pattern-dependent
    charging during etching in high-density plasmas*, J. Appl. Phys. **81**, 3433–3439 (1997).
    DOI 10.1063/1.365039.
32. G. S. Hwang and K. P. Giapis, *Electron irradiance of conductive sidewalls: A determining factor
    for pattern-dependent charging*, J. Vac. Sci. Technol. B **15**, 1741–1746 (1997).
    DOI 10.1116/1.589364.
33. G. S. Hwang and K. P. Giapis, *On the link between electron shadowing and charging damage*,
    J. Vac. Sci. Technol. B **15**, 1839–1842 (1997). DOI 10.1116/1.589336.
34. F. Krüger, S. Wilczek, T. Mussenbrock, J. Schulze, *Voltage waveform tailoring in radio frequency
    plasmas for surface charge neutralization inside etch trenches*, Plasma Sources Sci. Technol.
    **28**, 075017 (2019). DOI 10.1088/1361-6595/ab2c72.
35. G. Memos, E. Lidorikis, G. Kokkoris, *Roughness Evolution and Charging in Plasma-Based Surface
    Engineering of Polymeric Substrates: The Effects of Ion Reflection and Secondary Electron
    Emission*, Micromachines **9**, 415 (2018). DOI 10.3390/mi9080415. (Open access.)
36. T. Maruyama, T. Narukage, R. Onuki, N. Fujiwara, *High-aspect-ratio deep Si etching in SF₆/O₂
    plasma. II.*, J. Vac. Sci. Technol. B **28**, 862–868 (2010). DOI 10.1116/1.3466884.
37. Y. G. Yook *et al.*, *Fast and realistic 3D feature profile simulation platform for plasma
    etching process*, J. Phys. D: Appl. Phys. **55**, 255202 (2022). DOI 10.1088/1361-6463/ac58cf.
38. T. Nishizuka, R. Igosawa, T. Yokoyama, K. Sako, H. Moki, M. Honda (TEL), *Precise and practical
    3D topography simulation of high aspect ratio contact hole etch by using model optimization
    algorithm*, J. Vac. Sci. Technol. A **42**, 043003 (2024). DOI 10.1116/6.0003515.

**[META] — citation confirmed only**
39. N. Negishi, M. Miyake, K. Yokogawa, M. Oyama, T. Kanekiyo, M. Izawa, *Bottom profile degradation
    mechanism in high aspect ratio feature etching based on pattern transfer observation*,
    J. Vac. Sci. Technol. B **35**, 051205 (2017). DOI 10.1116/1.4998943.
40. M. Miyake *et al.*, *Effects of Mask and Necking Deformation on Bowing and Twisting in
    High-Aspect-Ratio Contact Hole Etching*, Jpn. J. Appl. Phys. **48**, 08HE01 (2009).
    DOI 10.1143/jjap.48.08he01.
41. J. Matsui, K. Maeshige, T. Makabe, *Effect of aspect ratio on topographic dependent charging in
    oxide etching*, J. Phys. D: Appl. Phys. **34**, 2950–2955 (2001). DOI 10.1088/0022-3727/34/19/304.
42. B. M. Radjenović, M. D. Radmilović-Radjenović, Z. Lj. Petrović, IEEE Trans. Plasma Sci. **36**,
    874 (2008).
43. T. Ohmori and T. Makabe, *In situ measurement of plasma charging on SiO₂ hole bottoms and
    reduction by negative charge injection during etching*, Appl. Surf. Sci. **254**, 3696–3709
    (2008). DOI 10.1016/j.apsusc.2007.10.070.
44. P. Gin *et al.*, *Inline metrology of high aspect ratio hole tilt and center line shift using
    small-angle x-ray scattering*, J. Micro/Nanopattern. Mater. Metrol. **22**, 031205 (2023).
    DOI 10.1117/1.jmm.22.3.031205.
45. W. W. Dostalik, S. Krishnan, T. Kinoshita, S. Rangan, *Electron shading effects in high density
    plasma processing for very high aspect ratio structures*, 3rd Int. Symp. Plasma Process-Induced
    Damage (1998), 160–163. DOI 10.1109/ppid.1998.725599.
46. Y. Zhai, R. Ge, Z. Hu, J. Li, H. Shao, J. Cheng, L. Filipović, R. Chen, *Modeling the charging
    effect of the hardmask and silicon substrate during plasma etching in advanced nodes*,
    J. Appl. Phys. **137**, 063302 (2025).
47. Y. Zhai *et al.*, *Modeling the Charging Effect on the Twisting Defects During High Aspect Ratio
    Etching of Dielectrics*, IEEE (2025), IEEE Xplore document 11017954.
48. M. Converse, J. Booske, A. Wendt, S. Gearhart, *In-situ measurements of polycrystalline silicon
    feature voltages for investigation of the notching effect during plasma etching*, ICOPS 1998,
    p. 223. DOI 10.1109/plasma.1998.677743. **[VERIFY]** whether a full journal paper exists.
49. K. Ishikawa *et al.*, *Progress in nanoscale dry processes for fabrication of high-aspect-ratio
    features: How can we control critical dimension uniformity at the bottom?*, Jpn. J. Appl. Phys.
    **57**, 06JA01 (2018). DOI 10.7567/jjap.57.06ja01.
50. Z. Lj. Petrović *et al.*, *Electrostatic potential fluctuation induced by charge discreteness in
    a nanoscale trench*, Phys. Plasmas **14**, 103501 (2007).

## 7. Open [VERIFY] items, ranked

1. **Does the Yook/Im (K-SPEED) GPU platform contain a charging module?** If not, retire "Korean
   CUDA charging MC" from project memory and re-point Route A at the Kushner MCFPM algorithm
   (Wang & Kushner 2010) as the reference implementation.
2. **Jinnai 2007 (JVST B 25, 1808) AR span and V(AR) figure values** — this is the only measured
   V(AR) curve; obtain the figure.
3. **Ohmori & Makabe, Appl. Surf. Sci. 254, 3696 (2008)** — a 14-page in-situ measurement study;
   highest-value unread source.
4. **Negishi 2017 (JVST B 35, 051205)** — is the "200 V at AR 20" measured or modelled?
5. **Kamata & Arimoto figures** (both 1996 papers) — needed for D2 digitization.
6. **Gin 2023 CD-SAXS tilt magnitudes and repeatability** — needed to turn D4 into an
   experiment-anchored gate rather than a simulation-anchored one.
7. **Zhai 2025 JAP validation figures** (publisher blocks WebFetch). DOI now RESOLVED:
   10.1063/5.0243470, *J. Appl. Phys.* **137**, 063302 (2025), Zhai, Ge, Hu, Li, Shao, Cheng,
   Filipovic, Chen — abstract confirms validation against Si/Cl2 SEM cross-sections.
8. **Dostalik 1998 P2ID** — potentially the only deep-AR shading-damage scaling dataset.
