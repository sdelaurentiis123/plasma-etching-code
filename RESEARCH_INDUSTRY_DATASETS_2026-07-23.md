# Industry-affiliated quantitative etch datasets — research catalog

Author: research pass 2026-07-23. Purpose: expand petch's stock of *industry-affiliated*
quantitative experimental etch datasets for (a) preregistered blind validation campaigns
(calibrate on a subset, predict held-out conditions) and (b) chemistry-parameter extraction.
petch has already blind-validated against Krüger 2024 (Kushner group / Samsung / TEL America
conditions). This document catalogs the next tranche.

Status: reviewed and committed 2026-07-24. Contains no proprietary material; all entries are published literature.

## Conventions & caveats

- Citation format: authors — title — venue vol(issue), article/page (year). DOI.
- `[VERIFY]` marks a fact I could not confirm to primary source in this pass (affiliation,
  exact DOI, exact figure number). Most journal bodies (AIP `pubs.aip.org`) 403 to automated
  fetch; where a Kushner-group PDF is openly mirrored at `cpseg.eecs.umich.edu/pub`, I note it.
- "Boundary characterizability" = can we build at least a Tier-1 global-model boundary
  (source plane IADF/IADF-ratio, neutral-to-ion flux ratio) from the paper? Kushner-group
  papers are strongest here because the reactor is solved with HPEM and the feature boundary
  (energy/angle distributions, flux ratios) is printed or reconstructable. Pure fab/vendor
  experimental papers usually give recipe knobs (power, pressure, flows) but not the
  feature-plane distributions — those need a Tier-1 CCP/ICP global surrogate to bridge.
- Evidence-class language follows `EXPERIMENTAL_VALIDATION_MATRIX.md` (E4 held-out, E3
  calibrated, E2 trend, S simulation-reference). A simulation paper with industry conditions
  is an **S**-reference for boundary/mechanism; its embedded experimental anchors are the
  E-grade data.

---

## Tier A — Kushner-group integrated reactor+feature (industry co-authored/funded, boundary fully characterizable)

These are computational papers, but they are the single best *boundary-characterized* source
we have: HPEM-solved reactor → MCFPM feature model, with printed IADFs, flux ratios, and
process sweeps, and Samsung/TEL co-authors or funding. They are simulation references (S) for
transport/boundary, and several carry experimental anchors. Krüger's own PhD thesis is the
richest single artifact.

### A1. Krüger 2024 — the paper petch already validated (extend it)
- **Krüger, F.** — *Modeling and Optimization of High Aspect Ratio Plasma Etching* — PhD
  thesis, University of Michigan (2024). Open PDF: `cpseg.eecs.umich.edu/pub/theses/Krueger_Florian_PhD_Thesis_2024.pdf`.
  Journal instance: **JVST A 42, 043008 (2024)**, DOI 10.1116/6.0003554 (C4F6/Ar/O2 SiO2;
  Table IV base case hf=825 nm, wt=90°, wm=45°; clog↔etch boundary vs ion energy incl. Fig 16).
- Measured/reported: SiO2 HAR trench/hole depth & profile vs low-frequency bias power, O2/C4F6
  ratio (rate rises-then-saturates), ion energy (clog boundary). Base-case metrics already in
  repo (`data/krueger_2024/base_case_metrics.csv`, `transfer_observations.csv`).
- Chemistry: fluorocarbon SiO2 HAR CCP — **exact match to petch's Krüger-calibrated mechanism.**
- Blind suitability: **HIGHEST.** The thesis contains many more swept conditions than the single
  base case petch calibrated on; hold out the O2-ratio and low-freq-power sweeps as E2/E4
  transfer. Boundary: fully characterized (HPEM IADFs printed).
- Action: this is not "new data" — it is *unused data we already partly own*. Digitize the
  thesis sweep figures NOT used in base calibration; preregister them as held-out.

### A2. Huang, Huard, Shim, Nam, Song, Lu, Kushner 2019 — Ar/C4F8/O2 SiO2 HAR
- **JVST A 37(3), 031304 (2019)**, DOI 10.1116/1.5090606. Open PDF:
  `cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf`. Co-authors Shim/Nam/Song =
  Samsung Electronics.
- Measured: SiO2 etch of HAR features to **AR up to 80**, tri-frequency CCP; depth, CD,
  twisting, contact-edge roughening, ARDE vs frequency/phase and flux composition. Multiple
  figures of profile vs AR and vs process knob — digitizable.
- Chemistry: fluorocarbon (C4F8) SiO2 — **match** (C4F8 vs C4F6 is a small mechanism delta).
- Blind suitability: **very high.** Several conditions → hold some out. Boundary characterized
  (tri-freq CCP solved; IADFs given). Complements A1 by extending AR range and to C4F8.

### A3. Huang, Shim, Nam, Kushner 2020 — pattern-dependent profile distortion
- **JVST A 38(2), 023001 (2020)**, DOI 10.1116/1.5125568 [VERIFY exact DOI]. Shim/Nam = Samsung.
- Measured: twisting, tilting, surface roughening of HAR SiO2 and **SiO2–Si3N4–SiO2 (ONO)
  stacks** vs pattern layout (isolated vs dense, edge vs interior), attributing distortion to
  flux randomness + charging + pattern dependency. Profile-distortion metrics vs AR and pattern.
- Chemistry: fluorocarbon SiO2 + ONO — **match + first ONO stack** (relevant to 3D NAND).
- Blind suitability: high for a *pattern-dependence* campaign (our pattern-dependent-charging/
  tilt row is currently "Open"). Boundary characterized. Would exercise the charging module.

### A4. Krüger, Lee, Nam, Kushner 2023 + 2024 — voltage-waveform-tailored SiO2 HAR
- **JVST A 41(1), 013006 (2023)**, DOI 10.1116/6.0002290 (ion/electron distributions →
  profiles). Companion: **Phys. Plasmas 31(3), 033508 (2024)**, DOI 10.1063/5.0189675 [VERIFY]
  (low fundamental-frequency biases). Lee/Nam = Samsung Mechatronics R&D; funded Samsung + TEL
  America. Ar/CF4/O2, 40 mTorr, geometrically asymmetric CCP.
- Measured: DC self-bias, IEDF, and resulting SiO2 HAR trench profile vs tailored-waveform
  harmonic content and fundamental frequency. Rich IEDF→profile mapping.
- Chemistry: fluorocarbon (CF4) SiO2 — **match**. Boundary: **best in class** (waveform → IEDF
  → profile all printed; ideal Tier-1 boundary).
- Blind suitability: high — waveform/frequency is a clean swept axis to hold out.

### A5. (supporting) Physical-sputtering evolution of HAR SiO2 trench/hole
- *Computational study of the evolution of high-aspect-ratio SiO2 trench and hole features
  during physical sputtering* — **Phys. Plasmas 33(3), 033901 (2026)** [VERIFY authors/Kushner].
  Use as a mechanism reference for the sputter channel of our fluorocarbon law; not a primary
  blind dataset.

---

## Tier B — Tokyo Electron (TEL) — real experimental 3D profiles

### B1. Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (TEL) 2024 — 3D HARC topography
- **JVST A 42(4), 043003 (2024)**, DOI 10.1116/6.0003515. All authors Tokyo Electron.
- Measured: **experimental 3D HARC hole profiles** (vertical *and* horizontal cross-sections,
  i.e. 3D-imaged), reproducing "distortion" and "twisting" (non-axisymmetric features along the
  hole axis). A cell-based Particle Monte Carlo simulator is fitted to the SEM/3D data via a
  model-optimization algorithm. Data lives in the profile figures (digitize vertical + azimuthal
  cross-sections) — this is genuine *experimental profile* data from a tool vendor.
- Chemistry: fluorocarbon HARC (dielectric) — **match**.
- Blind suitability: **high value but boundary-limited.** Real experimental profiles (E-grade),
  but the TEL reactor boundary is not fully published → needs a Tier-1 CCP surrogate. Best used
  as a held-out *profile-shape* target (twist/distortion) rather than absolute-rate calibration.

### B2. TEL cryogenic HAR dielectric etch — the "10 µm / 33 min" claims
- Primary technical venue: **2023 Symposium on VLSI Technology & Circuits** (TEL Miyagi);
  press/blog corroboration: `tel.com/news/product/2023/20230609_001.html`,
  `tel.com/blog/all/20241021_001.html`. TELAVES cryo etcher: >10 µm depth, ~2.5–3× faster,
  −40% power, −84% GWP, 400+ layer 3D NAND channel hole.
- Measured (published): depth (10 µm), time (33 min) → mean rate ~0.3 µm/min; limited swept
  conditions in the *public* material.
- Chemistry: cryogenic HF-based dielectric HAR — **needs a new HF/cryo condensation module**
  (petch's cryo row is a narrow calibrated anchor only).
- Blind suitability: **low as published** (headline numbers, not a DOE). Chase the VLSI 2023
  paper and any JVST/JJAP follow-up for the actual condition table.

---

## Tier C — Lam Research

### C1. Kanarik, Tan, Yang, Kim, Lill, Kabansky, Hudson, Ohba, Nojiri, Yu, Wise, Berry, Pan, Marks, Gottscho (Lam) 2017 — predicting ALE synergy
- **JVST A 35(5), 05C302 (2017)**, DOI 10.1116/1.4979019. Open PDF: `osti.gov/servlets/purl/1376399`.
- Measured: directional plasma ALE **etch-per-cycle (EPC) and synergy % for Si, Ge, C, W, GaN,
  SiO2**; energy dependence of EPC/synergy; e.g. Si (Cl2 + Ar+) EPC ≈ 0.5 nm/cycle; **W ALE at
  60 eV: synergy 95%, EPC 0.21 nm/cycle, sputter 0.01 nm/cycle**. Synergy scales with surface
  binding energy. Table + EPC-vs-energy curves are directly digitizable.
- Chemistry: ALE across materials incl. Si/Cl2/Ar+ (petch has a released 100 eV Si-Cl2-Ar+
  slice) and metals (W) — Si arm is a **match**; W/GaN/Ge/C need new ALE tables.
- Blind suitability: **excellent for ALE-module calibration** (0-D EPC/synergy). Hold out
  energies/materials. This is *the* quantitative ALE-synergy reference from a tool vendor.
- Companion: Kanarik et al., *Overview of atomic layer etching* — JVST A 33(2), 020802 (2015),
  DOI 10.1116/1.4913379 (context + more EPC anchors).

### C2. Lill, Wang, Wu, Oh, Kim, Wilcoxson, Singh, Ghodsi (Lam) + George (CU Boulder) + Barsukov, Kaganovich (PPPL) 2024 — low-T HF etching of SiO2/SiN
- **JVST A 42(6), 063006 (2024)**, DOI 10.1116/6.0003813 [VERIFY]. Open: `osti.gov/servlets/purl/2514386`.
- Measured: blanket-film **etch rate of SiO2 and SiN vs substrate temperature** with HF (and the
  chemistry split — pure HF etches SiN; SiO2 needs an added F source e.g. PF3). Rate-vs-T curves
  digitizable; the enabling mechanism (H2O catalysis) is discussed.
- Chemistry: HF cryogenic dielectric — **needs a new HF/low-T condensation-etch module.**
- Blind suitability: medium — blanket rate vs T (not profiles). Good for *chemistry extraction*
  of the HF/cryo module we don't yet have; pairs with the IMEC/TEL cryo direction.

---

## Tier D — Fab / IMEC / metrology experimental datasets

### D1. Deep-trench charging in 3D NAND (fab-authored)
- *Influence of accumulated charges on deep trench etch process in 3D NAND memory* —
  **Semicond. Sci. Technol. 35, 055012 (2020)**, DOI 10.1088/1361-6641/ab73e7. Affiliation
  [VERIFY — fab, likely YMTC/SK Hynix]. Channel-hole (CH) etch traps charge → common-source-line
  (CSL) slit **tilting**; charging model validated against experiment; mitigation via backside
  poly-Si.
- Measured: CSL tilt vs charge accumulation; AR > 50:1. E2/E3-grade experimental charging effect.
- Chemistry: fluorocarbon HAR dielectric + **charging** — directly relevant to petch's
  notching/charging frontier (currently "Open"/"partial"). Boundary not characterized → Tier-1
  surrogate needed; best as a qualitative charging-sign/held-out tilt-direction target.

### D2. IMEC HF-cryo mechanism paper
- *Reaction mechanism for HF based cryogenic plasma etching of SiO2* — **JVST A 44(3), 033006
  (2026)** [VERIFY affiliation — IMEC/collaborators]. Mechanistic (H2O-catalyzed HF etch of
  SiO2, PF3 as F source). Supports the HF/cryo module (D-tier chemistry). Pair with C2.

### D3. CD-SAXS 3D NAND channel-hole profile-variance metrology
- *CDSAXS study of 3D NAND channel hole etch pattern edge effects and etched hole pattern
  variance* — **Proc. SPIE 12955, 1295539 (2024)**, DOI 10.1117/12.3010927 [VERIFY affiliation —
  metrology/fab]. Companion: *3D reconstruction of 3D NAND etch profiles using FIB-SEM*, Proc.
  SPIE 13426 (2025).
- Measured: **CD, CD-profile, tilt, distortion, center-line-shift, pitch/hole-hole variance vs
  depth and vs pattern position (edge vs interior)** across process conditions — quantitative,
  statistical, held-out-rich profile metrology.
- Chemistry: fluorocarbon HAR dielectric — **match** (profile metrics). Boundary: none (metrology
  only) → Tier-1 surrogate needed. Excellent *held-out geometric-target* dataset (CD-vs-depth
  and tilt statistics), weak for absolute-rate calibration.

---

## Tier E — Hitachi & modeling references

### E1. Hitachi charging/notch simulation-validation lineage
- *Charging effect simulation model used in simulations of plasma etching of silicon* — **J.
  Appl. Phys. 112(8), 084308 (2012)**, DOI 10.1063/1.4759313 [VERIFY]. 2D MC of poly-Si
  overetch; validated against the **Nozawa et al. notch measurements** (line-connectivity and
  open-area-width dependence of notch depth). Directly relevant to petch's charging-driven-notch
  row (already E2/partial). Digitize notch-depth vs open-area and vs connectivity.
- Review (chase for a dataset index): *Review and perspective of dry etching and deposition
  process modeling of Si and Si dielectric films for advanced CMOS* — **Jpn. J. Appl. Phys.
  (2024)**, DOI 10.35848/1347-4065/ad5355 [VERIFY authors — Hitachi]. Use to find more
  vendor-validated condition tables.

### E2. Academic DOE with decoupled ion energy/flux (chemistry extraction)
- *Contribution of Ion Energy and Flux on High-Aspect-Ratio SiO2 Etching in a Dual-Frequency
  CCP Ar/C4F8: individually controlled ion energy and flux* — Nanomaterials (2023),
  PMC10222222 [VERIFY exact cite]. **Independently swept ion energy vs ion flux** → ideal for
  disentangling the ion-enhanced vs physical-sputter terms of our fluorocarbon SiO2 law. SKKU
  (academic, not vendor) but the decoupled-DOE structure is rare and high-value.

### E3. de Boer / Blauw SF6-O2 silicon ARDE (already partly in repo — keep for the Si arm)
- de Boer et al., *Ultrahigh AR etching of Si in SF6-O2 (CORE sequence, Cr mask)* — **JVST A
  38(5), 053002 (2020)**, DOI 10.1116/6.0000254 [VERIFY]. Blauw et al., *Kinetics and crystal
  orientation dependence in HAR Si dry etching* — **JVST B 18(6), 3453 (2000)** (Knudsen-fit
  ARDE). Si SF6/O2 — **matches petch's Belen mechanism.** Note: the repo's prior E4 claim on the
  fitted Fig-9 curve was withdrawn; use for the *silicon* ARDE arm only, never re-use the exposed
  Fig-9 points as held-out.

---

## Priority ranking — top 8 for the next blind campaign

Ranked by (data richness × boundary characterizability × chemistry match to current mechanisms).
"New module?" flags chemistry we do not yet have.

| # | Dataset | Richness | Boundary | Chem match | New module? |
|---|---|---|---|---|---|
| 1 | **A1 Krüger 2024 thesis** (unused sweeps beyond base case) | high | full (HPEM) | exact (FC SiO2) | none — extend existing calibration |
| 2 | **A2 Huang/Huard 2019** Ar/C4F8/O2 AR≤80 | high | full | match (C4F8) | small FC delta |
| 3 | **A4 Krüger 2023/24** voltage-waveform SiO2 | high | best (IEDF→profile) | match (CF4) | none |
| 4 | **B1 Nishizuka/Honda TEL 2024** 3D HARC profiles | high (real 3D) | partial (needs Tier-1) | match | none (profile-shape target) |
| 5 | **A3 Huang 2020** pattern-dependent distortion + ONO | med-high | full | match + ONO | exercises charging module |
| 6 | **C1 Kanarik 2017 ALE synergy** (Lam) | high (EPC/synergy table) | 0-D (n/a) | Si arm match | W/GaN/Ge ALE tables |
| 7 | **D3 CD-SAXS 3D NAND** profile-variance metrology | high (statistical) | none (needs Tier-1) | match (metrics) | none (geometric target) |
| 8 | **E3 de Boer/Blauw SF6-O2** Si ARDE | med | partial | match (Si) | none |

**Immediately actionable, zero new chemistry:** #1, #2, #3, #4, #8 — all fluorocarbon-SiO2 or
Si-SF6/O2, matching petch's live mechanisms. #1 is the highest-leverage because it is *data we
already partially hold* (Krüger thesis) whose non-base-case sweeps have never been scored as
held-out.

**Needs new module before use:** C1 metals arm (W/GaN/Ge ALE tables), C2/D2 (HF low-T /
condensation module), B2 TEL cryo (HF/cryo + condition table not yet public).

**Boundary caveat:** Tier-A (Kushner/Samsung/TEL-funded) papers are the only ones that give a
publication-complete feature-plane boundary. Every fab/metrology dataset (B1, D1, D3) needs a
Tier-1 CCP/ICP global-model surrogate to turn recipe knobs (power/pressure/flow) into the
IADF + flux-ratio boundary petch consumes — build that surrogate once and it unlocks B1/D1/D3
together.

## Leads still to chase (next pass)
- Krüger PhD thesis (Michigan, open PDF): extract every non-base sweep figure → preregister.
- TEL VLSI 2023 cryo paper (and any JJAP/JVST follow-up) for the cryo condition table.
- Confirm affiliations flagged [VERIFY]: D1 (fab), D2/D3 (IMEC/metrology), E1/E2 review authors.
- Applied Materials Sym3-era etch: no clean vendor-authored quantitative JVST hit this pass —
  worth a dedicated search (AMAT authors publish less openly than Lam/TEL/Kushner).
- IMEC-authored (not just IMEC-adjacent) HAR/patterning DOE with printed CD-vs-depth tables.
