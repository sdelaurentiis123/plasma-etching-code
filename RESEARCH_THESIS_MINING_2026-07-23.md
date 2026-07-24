# PhD-thesis & measured-constants mining for petch — knob-retirement + blind-validation catalog

Research pass 2026-07-23. Author-of-record: autonomous research agent.
Status: survey document (reviewed and committed 2026-07-24).

## Purpose & lane

This catalog mines **PhD theses and their closely-attached papers** from the five major
plasma-etch simulation / surface-science groups, for two petch programs:

- **(a) Blind held-out validation** — datasets with *fully specified boundary conditions*
  (fluxes, ion energies, gas ratios), the completeness that made the Krüger 2024 blind test
  possible.
- **(b) Knob retirement** — *measured* chemical/physical constants (sticking coefficients,
  mixed-layer thickness vs ion energy, F/C/O/Si composition, G-values, yield thresholds,
  angular yield curves, SEE/charging parameters) that replace fitted petch parameters with
  derived/measured physics.

**Deliberate non-overlap with `RESEARCH_INDUSTRY_DATASETS_2026-07-23.md`.** That file owns the
*industry-affiliated reactor+feature datasets* (Krüger 2024 thesis A1, Huang/Huard 2019 A2,
Krüger waveform A4, TEL/Lam/IMEC). This file owns the **student theses and the measured-constants
literature attached to them** — the beam-yield databases, MD mixed-layer studies, and
surface-science film measurements that the industry papers *cite* as their kinetic inputs. Where a
thesis is already in the industry file (Krüger, Huard's 2019 paper), I cross-reference and do not
re-score it; I add the *other* theses from the same group.

Citation convention: `[VERIFY]` marks an exact locator (DOI / volume / page / year) reconstructed
from the search record and not confirmed against the publisher's page this pass. Author, group,
venue, and the physics claim are checked against the search record; only the precise locator
carries the flag. Most AIP/AVS bodies (`pubs.aip.org`) 403 automated fetch; Kushner theses are
openly mirrored at `cpseg.eecs.umich.edu/pub/theses/`.

petch knob vocabulary (targets to retire), from `RESEARCH_MIXED_LAYER_DESIGN_2026-07-23.md` and
memory: **yield energy law** (`Y=A(√E−√E_th)` / `(E−E_th)^n`, fitted A & E_th), **mixed-layer
ledger** (element-resolved C/F/O/Si of top 1–3 nm), **crosslink E_x** (ion-processed-skin dose
competition), **oxygen saturation**, **complex_formation_probability** (Krüger 0.2729 etc.),
**sticking coefficients** (CFx radicals), **IADF/IEADs**, **charging/SEE**, **LER/PSD**, **ARDE**.

---

## RANKED CATALOG

Ranked within each group by (measured-data richness × boundary completeness × directness of the
knob it retires). Group order follows the task priority (Kushner → Graves → Sawin → Oehrlein →
Demokritos), but the **highest knob-retirement leverage is in Graves/Sawin** (beam + MD give the
constants), while the **highest blind-validation readiness is in Kushner** (complete boundaries).

---

### 1. KUSHNER group (Michigan) — complete-boundary feature-scale theses

The single best *boundary-characterized* source: HPEM-solved reactor → MCFPM feature model with
printed IADFs and flux ratios. Theses live openly at `cpseg.eecs.umich.edu/pub/theses/`.
(Krüger 2024 thesis is already catalogued as **industry-file A1** — extend its unused sweeps;
not re-scored here.)

#### K1. Huard, C. M. — *Nano-Scale Feature Profile Modeling of Plasma Material Processing* — PhD, U. Michigan (EE), 2018
- PDF: `cpseg.eecs.umich.edu/pub/theses/huard_chad_phd_thesis.pdf`.
- MEASURED/simulated: the **MCFPM methodology thesis** — voxel feature model, ALE of Si (Cl2+Ar)
  and SiO2 (fluorocarbon), wafer-scale ALE uniformity, mask/charging interplay. Companion paper
  Huard, Lanham, Kushner, *J. Phys. D* 51 (2018) (ALE uniformity in ICP). Also the source of the
  **complex-formation / polymer-inventory mechanism structure** petch's `surface_kinetics.py`
  descends from.
- Retires / addresses: **complex_formation_probability**, **mixed-layer/polymer-inventory** model
  form, ALE synergy structure. Digitizable: EPC-vs-energy and ALE synergy curves; profile-vs-cycle.
- Boundary: full (HPEM IADFs). Blind readiness: **medium-high** for an ALE-synergy campaign
  (0-D/1-D), lower for absolute profiles (methodology, not a single DOE).

#### K2. Zhang, Yiting — 3-D feature profile / ion-tilting thesis — PhD, U. Michigan, ~2017
- PDF: `cpseg.eecs.umich.edu/pub/theses/zhang_yiting_phd_thesis.pdf`. Core paper: Zhang, Huard,
  Sriraman, Belen, Paterson, Kushner, *JVST A* 35(2), 021303 (2017), "Investigation of feature
  orientation and consequences of ion tilting during plasma etching with a 3-D feature profile
  simulator" [VERIFY DOI].
- MEASURED/simulated: **ion-tilt → sidewall-angle / microtrench / feature-orientation** mapping;
  Cl2/HBr poly-Si and dielectric. Quantifies how IEADF tilt propagates to profile asymmetry.
- Retires / addresses: **IADF/IEADs** boundary consumption, tilt→profile (petch's tilt/twist row).
  Digitizable: sidewall-angle vs ion-tilt-angle curves.
- Boundary: full. Blind readiness: **medium** (clean tilt axis to hold out).

#### K3. Huang, Shuo — HAR SiO2 dielectric etch thesis — PhD, U. Michigan, ~2018
- CV: `cpseg.eecs.umich.edu/pub/vita/CV_Shuo_Huang_Aug_2018.pdf`. Core papers are **industry-file
  A2** (Huang/Huard 2019, Ar/C4F8/O2, AR≤80) and **A3** (Huang 2020, pattern-dependent distortion
  + ONO). Listed here only to mark the thesis as the umbrella artifact; use A2/A3 scoring.
- Retires / addresses: fluorocarbon SiO2 HAR mechanism, ARDE, pattern-dependent charging.
- Cross-ref industry-file A2/A3; **do not double-count.**

#### K4. Kushner group — *Modeling the charging effect of the hardmask and silicon substrate during plasma etching in advanced nodes* — *J. Appl. Phys.* 137(6), 063302 (2025)
- `pubs.aip.org/aip/jap/article/137/6/063302`. Recent Kushner-group charging model (not a thesis,
  but the group's current charging-constants paper — the attached-literature companion to the
  thesis line).
- MEASURED/simulated: charging-driven **mask faceting, bowing, microtrenching**; simultaneous
  charging + particle-reflection algorithm; positive-ion redirection to sidewalls. Prior
  companion: dc-augmented CCP high-energy-electron flux → twisting in HAR dielectric.
- Retires / addresses: **charging/SEE** module structure (petch's charging frontier, "Open").
  Digitizable: microtrench-depth / bow vs charging-on/off; field maps.
- Boundary: full. Blind readiness: **medium** (mechanism/sign target, not absolute-rate).

---

### 2. GRAVES group (Berkeley / Princeton) — MD mixed-layer + beam theses (the knob-retirement mother lode)

Graves' MD and beam work is where petch's **mixed-layer ledger and crosslink E_x** get their
first-principles anchors. These measure/compute the *very quantities* the mixed-layer design doc
wants: a-Si:C(:F) mixed-layer depth vs ion energy, F transport (G-value) through FC films, and
steady-etch composition windows.

#### G1. Humbird, D. — *Computational Studies of Plasma–Surface Interactions* — PhD, UC Berkeley (ChemE), 2004
- MEASURED/computed (MD): the **mixed amorphous a-Si:C top layer** that forms under ion impact;
  its **depth as a function of ion energy and incident radical C/F composition** — this is the
  physical object petch's mixed-layer ledger models. Attached papers:
  - Humbird & Graves, *Fluorocarbon plasma etching of silicon: factors controlling etch rate*,
    *J. Appl. Phys.* 96(1), 65 (2004) — etch-rate control by mixed a-Si:C layer.
  - Humbird & Graves, *Mechanism of Si etching in the presence of CF2, F, and Ar+*,
    *J. Appl. Phys.* 96(5), 2466 (2004).
  - Humbird & Graves, *MD simulations of Ar+-induced transport of fluorine through fluorocarbon
    films*, *Appl. Phys. Lett.* 84(7), 1073 (2004) — **F-transport / effective G-value** for
    ion-driven F delivery through the FC overlayer.
- Retires / addresses: **mixed-layer ledger** (a-Si:C thickness vs E_ion, C/F), **crosslink E_x**
  (ion-processed-skin), **sticking + F-transport constants**. Digitizable: layer-depth-vs-energy
  and F-flux-through-film curves.
- Boundary: MD boundary conditions fully specified (ion species, energy, flux ratio, C/F) →
  **directly usable as a first-principles calibration target for the stopping-based mixed-layer
  module**, not a reactor blind test. Readiness for *knob retirement*: **HIGHEST in this file.**

#### G2. Végh, J. J. — MD of Si etch in fluorocarbon + Ar+ (steady-etch-with-film window) — Graves group, Berkeley, ~2005–08
- Core papers: Végh, Humbird, Graves, *Silicon etch by fluorocarbon and argon plasmas in the
  presence of fluorocarbon films*, *JVST A* 23(6), 1598 (2005); *Silicon etch in the presence of
  a fluorocarbon overlayer: the role of fluorocarbon cluster ejection*, *JVST A* 26(1), 52 (2008).
- MEASURED/computed (MD): the **input window (C/F ratio, neutral/ion flux ratio, ion energy)** in
  which *steady Si etch coexists with a steady FC film* — outside it, either film with no etch or
  etch with no film. Plus **large-FC-cluster ejection** as a removal channel (a distinct volatile
  pathway from monomer SiFx).
- Retires / addresses: **complex_formation / polymer-inventory steady-state closure**, the
  etch↔clog boundary from first principles; **cluster-ejection removal term** (missing volatile
  channel in petch's ledger). Digitizable: steady-etch phase boundary in (C/F, N/I, E) space.
- Boundary: MD, specified. Readiness for knob retirement: **high** (defines the ledger's
  steady-state existence region).

#### G3. Gray / Coburn / Graves — *Vacuum beam studies of fluorocarbon radicals and argon ions on Si and SiO2 surfaces* — Graves group beam apparatus, ~2004 [VERIFY exact author list & venue]
- MEASURED (experiment, not simulation): Si and SiO2 exposed to c-C4F8 ± Ar+, and to a
  *characterized* CxFy radical+ion mixture (CF, CF2, CF3, heavy CxFy), with **neutral-FC/Ar+ flux
  ratio and Ar+ energy varied**. Reports: **CF2 sticking coefficient on FC layers vs on SiO2**;
  the **fluorine content of the steady-state FC layer drops as Ar+ is added** (ion-driven F
  depletion / mixing); FC-layer thickness (ellipsometry) + composition (XPS).
- Retires / addresses: **sticking coefficients** (CFx, the biggest un-measured petch input),
  **mixed-layer F-composition vs ion flux/energy**, **oxygen/F saturation** of the film.
  Digitizable: sticking-coefficient values, F-fraction-vs-Ar+ curves, thickness-vs-flux-ratio.
- Boundary: beam BCs fully specified. Readiness for knob retirement: **HIGH — this is the direct
  experimental sticking-coefficient source** (pairs with G1/G2 MD to cross-check).

---

### 3. SAWIN group (MIT) — the canonical angular-yield rate-table theses

Sawin's ICP-beam apparatus is the **canonical source of angle-resolved etch-yield databases** —
exactly the `Y(E,θ)` surfaces petch's yield law approximates. These theses ARE the rate tables.

#### S1. Yin, Yunpeng — angular-yield-database thesis — PhD, MIT (ChemE), ~2007–08
- Core paper: Yin & Sawin, *Angular etching yields of polysilicon and dielectric materials in
  Cl2/Ar and fluorocarbon plasmas*, *JVST A* 26(1) (Jan 2008), OSTI biblio 21020898; companion
  Yin & Sawin, *Impact of etching kinetics on roughening of thermal SiO2 and low-k coral in FC
  plasmas*, *JVST A* 25(4), 802 (2007).
- MEASURED (ICP-beam experiment): **angular etch-yield curves Y(θ)** for **poly-Si (Cl2/Ar)** and
  **thermal SiO2 + low-k coral (C4F8/Ar)** as functions of **ion energy, incidence angle, and
  effective neutral-to-ion flux ratio**. Key result: N/I ratio is the primary control — low N/I →
  sputter-like curve peaking ~60–70° off-normal; high N/I → ion-enhanced-etch-like curve peaking
  ~65°. Grazing-angle roughening via ion channeling.
- Retires / addresses: **yield energy law + angular factor** directly (the fitted `A(√E−√E_th)`
  and angular `f(θ)` become *measured* curves); **N/I-ratio dependence** of the yield shape.
  Digitizable: Y-vs-angle families at several energies and N/I ratios — a full `Y(E,θ,N/I)` grid.
- Boundary: beam BCs fully specified. Readiness: **HIGHEST for the yield-law knob** — this is the
  reference rate table for both the Si (Cl2/Ar) and SiO2 (FC) arms.

#### S2. Guo, Wei — *3-Dimensional Modeling and Simulation of Surface and Sidewall Roughening During Plasma Etching* — PhD, MIT, 2008
- MIT DSpace `1721.1/43201`. Companion: Guo & Sawin, *Etching of SiO2 in C4F8/Ar plasmas. II.
  Simulation of surface roughening and local polymerization*, *JVST A* 28(2), 259 (2010).
- MEASURED/simulated: 3-D MC profile simulator vs **reactive-ion-beam-etch (RIBE) roughening
  experiments** — surface morphology vs ion incidence angle (smooth 0–45°, striations ⟂ beam at
  60–75°); local-polymerization roughening mechanism.
- Retires / addresses: **LER/PSD** (sidewall-roughening angular dependence), local-polymerization
  coupling to roughness. Digitizable: RMS-roughness / striation-onset vs angle.
- Boundary: RIBE BCs specified. Readiness: **medium-high** for the LER/roughening arm
  (complements Demokritos D-group below with a beam-calibrated angular dataset).

---

### 4. OEHRLEIN group (Maryland / prior IBM) — steady-state FC-film surface-science theses

Oehrlein's XPS/ellipsometry work measures the **steady-state fluorocarbon film thickness and C/F
composition vs bias/ion energy** — the reservoir sitting above petch's mixed layer, and the
SiO2/SiN selectivity data.

#### O1. Metzler, Dominik — fluorocarbon-ALE thesis — PhD, U. Maryland (MSE), ~2016
- Core papers: Metzler, Li, Engelmann, Bruce, Joseph, Oehrlein, *Fluorocarbon assisted ALE of
  SiO2 and Si using cyclic Ar/C4F8 and Ar/CHF3 plasma*, *JVST A* 34(1), 01B101 (2016), OSTI
  1225188; Li, Metzler, Lai, Hudson, Oehrlein, *FC-based ALE of Si3N4 and SiO2/Si3N4 selectivity*,
  *JVST A* 34(4), 041307 (2016); *Scaling of ALE of SiO2: transient etching and surface
  roughness*, *JVST A* 39(3), 033003 (2021).
- MEASURED (experiment): **FC film thickness deposited per cycle**, **etch-per-cycle vs ion
  energy and step length**, **SiO2/Si3N4 selectivity**, XPS surface chemistry. Key quantitative
  anchor: CHF3 deposits a *thicker* FC layer on Si than on SiO2 (selectivity origin).
- Retires / addresses: **oxygen saturation / FC-reservoir thickness**, **complex_formation** (FC
  film per cycle), **selectivity** (SiO2 vs SiN vs Si). Digitizable: film-thickness-per-cycle and
  EPC-vs-energy tables, selectivity ratios.
- Boundary: ALE BCs (Ar+ energy, C4F8 dose) specified. Readiness: **high for a FC-ALE blind
  campaign** (0-D/1-D EPC + film-thickness targets).

#### O2. Oehrlein-group steady-state-FC-film corpus (Standaert / Schaepkens / Hua / Bruce / Engelmann theses) — Maryland/IBM, 1997–2010
- Anchor papers: *Role of fluorocarbon film formation in the etching of Si, SiO2, Si3N4, and
  a-SiC:H*, *JVST A* 22(1), 53 (2004); *Role of steady-state fluorocarbon films in the etching of
  SiO2 using CHF3 in an ICP reactor*, *JVST A* 15(4), 1881 (1997) [Standaert et al.].
- MEASURED (XPS + ellipsometry): **steady-state FC film thickness vs self-bias / ion energy** and
  **C/F ratio** across Si / SiO2 / SiN / SiC; the etch-suppression-by-thick-film mechanism
  underlying selectivity.
- Retires / addresses: **FC-reservoir thickness vs bias**, **C/F composition** input to the mixed
  layer, **selectivity ordering**. Digitizable: film-thickness-vs-bias curves, C/F-vs-bias, etch-
  rate-vs-film-thickness. [VERIFY exact per-thesis authorship — Standaert & Schaepkens are the
  1997–2001 steady-state-film theses; Hua/Bruce/Engelmann/Weilnboeck are the 2004–2010 roughening
  + composition theses.]
- Boundary: reactor-recipe BCs (bias, pressure, flow) — **not** feature-plane distributions → for
  *chemistry-constant extraction* (film thickness/composition), not a feature blind test.
  Readiness for knob retirement: **high** (film-reservoir + C/F constants).

---

### 5. GOGOLIDES / KOKKORIS group (NCSR Demokritos) — LER/PSD + charging-constants theses

The reference phenomenological FC surface model and the charging/SEE constants live here; also the
canonical LER-transfer + PSD experimental datasets.

#### D1. Kokkoris, George — detailed-FC-surface-model + combined-simulator thesis — NCSR Demokritos / NTU Athens, ~2005
- Anchor papers: Gogolides, Vauvert, Kokkoris, Turban, Boudouvis, *Etching of SiO2 and Si in
  fluorocarbon plasmas: a detailed surface model accounting for etching and deposition*, *J. Appl.
  Phys.* 88, 5570 (2000), DOI 10.1063/1.1311808; Kokkoris et al., *A combined simulator coupling
  surface etching, local flux calculation, and profile evolution* (SiO2/Si feature etching).
- MEASURED/modeled: the **canonical phenomenological FC coverage model** (θ_poly, θ_CFx, θ_F
  balance) reproducing the etch→deposition (clog) transition and RIE lag/ARDE. This is the model
  petch's `fluorocarbon_lamagna.py` / `surface_kinetics.py` lineage descends from — its rate
  constants and site-balance structure are a lift-able reference.
- Retires / addresses: **complex_formation_probability / coverage-balance structure**, **ARDE**.
  Digitizable: clog-boundary curves, ARDE-vs-AR. (Already discussed in
  `RESEARCH_MIXED_LAYER_DESIGN_2026-07-23.md` §1.1 — cross-ref, use for constants.)
- Boundary: modeled. Readiness: reference model / constants, not a blind test.

#### D2. Memos, Georgios — charging + SEE + ion-reflection thesis — NCSR Demokritos, ~2018
- Anchor papers: Memos & Kokkoris, *Roughness evolution and charging in plasma-based surface
  engineering of polymeric substrates: the effects of ion reflection and secondary electron
  emission*, *Micromachines* 9(8), 415 (2018); Memos, Gerardis, Kokkoris, *Charging effect in
  basic and complex mask patterns during plasma etching*, *Plasma Chem. Plasma Process.* (2022),
  DOI 10.1007/s11090-022-10277-9; Memos & Kokkoris, *Modeling of charging on unconventional
  surface morphologies of PMMA during Ar plasma etching* (2016).
- MEASURED/modeled: **charging Monte Carlo with explicit SEE yield + ion-reflection** parameters;
  surface-potential build-up in trench/mask patterns; the constants (SEE yield curves, reflection
  coefficients) that petch's charging plan (W2 "SEE Memos-Kokkoris") explicitly names.
- Retires / addresses: **charging/SEE** (directly — this is the plan's cited source), ion
  reflection. Digitizable: SEE-yield-vs-energy, reflection-coefficient tables, potential maps.
- Boundary: specified for the charging sub-model. Readiness: **high for the charging/SEE knob**
  (constants + validation morphologies).

#### D3. Constantoudis, Vassilios — LER/PSD metrology + transfer thesis/corpus — NCSR Demokritos
- Anchor papers: Constantoudis, Gogolides et al., *LER transfer during plasma etching: modeling
  approaches vs experiment* (SPIE 7273, 2009); *3-D geometrical modeling of plasma transfer effects
  on LER: comparison with experiments and rules of thumb*, *J. Micro/Nanolith. MEMS MOEMS* 12(4),
  041310 (2013); *2D and 3D photoresist line roughness characterization* (2013).
- MEASURED (experiment): **PSD-resolved LER/LWR of photoresist and its transfer to underlayers**;
  3-D CD-AFM sidewall reconstruction; correlation-length + roughness-exponent extraction.
- Retires / addresses: **LER/PSD** (the experimental transfer dataset for petch's LER modality —
  pairs with `RESEARCH_LER_EXPERIMENTAL_SOURCES_2026-07-21.md`). Digitizable: PSD curves,
  LER-in vs LER-out, correlation length vs etch time.
- Boundary: recipe-level (not feature-plane distributions) → LER-transfer target + PSD constants.
  Readiness: **high for the LER arm** (statistical held-out targets).

---

## TOP 3 NEXT BLIND CAMPAIGNS (ranked by blind-validation readiness)

Ranked by (fully-specified BCs × chemistry match to a live petch mechanism × held-out richness).
"Blind" here means: calibrate on a subset, predict a held-out subset the model never saw.

1. **Sawin/Yin angular-yield grid (S1) as a held-out `Y(E,θ,N/I)` prediction.** petch's yield law
   is currently a fitted `A(√E−√E_th)·f(θ)`. Calibrate A & E_th on *one* (energy, N/I) slice of
   Yin's poly-Si (Cl2/Ar) or SiO2 (C4F8/Ar) beam data, then **blind-predict the other energies,
   angles, and N/I ratios**. BCs are beam-complete (ion energy, angle, N/I all specified) — this
   is the cleanest yield-law blind test available and it retires the yield knob simultaneously.
   Zero new chemistry. **Do this first.**

2. **Oehrlein/Metzler FC-ALE (O1) EPC + film-per-cycle held-out.** Calibrate the FC-reservoir /
   complex-formation on one ion-energy + one C4F8-dose point, blind-predict EPC-vs-energy and
   film-thickness-per-cycle across the rest, plus **SiO2/Si3N4 selectivity** as an independent
   held-out. Directly exercises oxygen-saturation + FC-reservoir + selectivity together. BCs (Ar+
   energy, dose) are specified. Small new module (cyclic FC-ALE driver) but same chemistry family.

3. **Kushner Krüger-2024-thesis unused sweeps (industry-file A1) — still the top *reactor* blind
   test.** Already flagged in the industry file; repeated here because among *theses* it remains
   the only artifact with a publication-complete HPEM feature-plane boundary AND many swept
   conditions beyond the base case petch calibrated on. Hold out the O2-ratio and low-frequency-
   power sweeps. **Highest reactor-scale realism; overlaps industry file (coordinate, don't
   double-run).** If a fresh reactor thesis is wanted instead, K2 (Zhang ion-tilt) gives a clean
   held-out tilt axis.

Runner-up (constants-first, not "blind" but highest knob-retirement value): **Graves G1 (Humbird
MD mixed-layer) + G3 (beam sticking coefficients)** — calibrate the ZBL-stopping mixed-layer
module against MD a-Si:C depth-vs-energy, then check the *measured* CF2 sticking coefficient falls
out. This is the derive-then-measure loop the knob-elimination doctrine wants.

---

## CONSTANTS WE CAN LIFT DIRECTLY

Measured/computed constants extractable from the above, mapped to the petch knob each replaces.
`[digitize]` = value lives in a figure needing extraction; `[table]` = printed numeric value.

| Constant (measured) | Source | petch knob retired | Form | Readiness |
|---|---|---|---|---|
| Angular etch-yield curves `Y(θ)` for poly-Si (Cl2/Ar) & SiO2 (C4F8/Ar) vs ion energy & N/I ratio | Yin & Sawin, *JVST A* 26(1) 2008 (S1) | yield energy law + angular `f(θ)` | `[digitize]` Y-vs-θ families | **direct** |
| Sputter-vs-ion-enhanced peak angle (~60–70° off-normal), N/I threshold for regime switch | Yin & Sawin (S1) | angular factor regime switch | `[digitize]` | direct |
| a-Si:C(:F) **mixed-layer depth vs ion energy & incident C/F** | Humbird thesis 2004 / *JAP* 96(1) 65 (G1) | mixed-layer ledger thickness; crosslink E_x | `[digitize]` depth-vs-E | **direct** |
| **F-transport (effective G-value)** through FC film under Ar+ | Humbird & Graves, *APL* 84(7) 1073 (G1) | ion-driven F delivery term in ledger | `[digitize]` F-flux-vs-E | direct |
| **CF2 sticking coefficient** on FC layer vs on SiO2 | Gray/Coburn/Graves vacuum-beam ~2004 (G3) [VERIFY] | sticking coefficients (CFx) | `[table]`/`[digitize]` | direct |
| **F-fraction of steady FC layer decreases with Ar+ addition** (ion-driven F depletion) | G3 beam study | mixed-layer F composition vs ion flux; O/F saturation | `[digitize]` F%-vs-Ar+ | direct |
| Steady-etch **existence window** in (C/F, N/I, E_ion) for Si-with-FC-film | Végh/Humbird/Graves *JVST A* 23(6) 1598 (G2) | complex_formation steady-state closure | phase boundary | high |
| **FC-cluster-ejection** removal channel (non-monomer volatile) | Végh & Graves *JVST A* 26(1) 52 (G2) | missing volatile term in ledger | mechanism + rate | medium |
| **Steady-state FC film thickness vs self-bias / ion energy** (Si, SiO2, SiN, SiC) | Oehrlein *JVST A* 22(1) 53 (2004); 15(4) 1881 (1997) (O2) | FC-reservoir thickness; selectivity | `[digitize]` thickness-vs-bias | direct |
| **FC film thickness per cycle & EPC vs ion energy**; SiO2/Si3N4 selectivity | Metzler/Oehrlein *JVST A* 34(1) 01B101 (2016) (O1) | complex_formation; oxygen saturation; selectivity | `[table]`/`[digitize]` | direct |
| **C/F ratio of FC film vs bias** (XPS) | Oehrlein corpus (O2) | mixed-layer C/F input | `[digitize]` | direct |
| **SEE yield curves + ion-reflection coefficients** (charging MC) | Memos & Kokkoris *Micromachines* 9(8) 415 (2018); PCPP 2022 (D2) | charging/SEE (plan W2) | `[digitize]`/`[table]` | direct |
| Coverage-balance rate constants (θ_poly/θ_CFx/θ_F), clog-boundary, ARDE-vs-AR | Gogolides/Kokkoris *JAP* 88 5570 (2000) (D1) | complex_formation structure; ARDE | model constants | reference |
| **PSD / correlation-length / roughness-exponent** for LER transfer (PR → underlayer) | Constantoudis/Gogolides SPIE 7273 (2009); JMM 12(4) 041310 (2013) (D3) | LER/PSD | `[digitize]` PSD curves | high |
| ALE EPC + synergy % (Si Cl2/Ar; SiO2 FC) from MCFPM methodology | Huard thesis 2018 (K1) | complex_formation; ALE synergy | `[digitize]` EPC-vs-E | medium |

---

## Verification / next-pass to-dos
- Confirm G3 exact author list + venue/year of "Vacuum beam studies of fluorocarbon radicals and
  argon ions on Si and SiO2 surfaces" (search attributed to Graves group ~2004; author order
  unconfirmed — likely Gray/Coburn/Graves lineage). Pull the CF2 sticking-coefficient numeric.
- Confirm per-thesis authorship of the O2 steady-state-FC-film corpus (Standaert & Schaepkens =
  1997–2001; Hua/Bruce/Engelmann/Weilnboeck = 2004–2010).
- Verify DOIs flagged: K2 (021303/2017), and pull exact page for Yin & Sawin 2008 (OSTI 21020898).
- Fetch Kushner theses directory `cpseg.eecs.umich.edu/pub/theses/` (cert/verify issue this pass)
  to enumerate any additional recent charging/HAR theses beyond Huard/Zhang/Huang/Krüger.
- Digitize priority: S1 (yield grid) → G1 (mixed-layer depth) → O1 (FC-ALE EPC) — these three feed
  the top-2 blind campaigns and the constants-first loop.
