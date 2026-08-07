# petch internal literature library

Every source this project has drawn on, cited to source, with the claims we
actually consumed. Built so no future pass has to re-hunt a paper we already read.

## Two conventions (binding)

1. **Fetch ⇒ extract + entry in the same commit.** If you fetch a PDF, extract it to
   `research_sources/thesis_extracts/<name>.txt` *and* create/update
   `research_sources/library/<bibkey>.md` in the same commit. A fetched-but-unarchived
   source does not exist for the next agent.
2. **Provenance strings use bibkeys.** Deck and code provenance should name the bibkey
   (e.g. `gray-1993-thesis`), so any constant is one grep from its claim row, its
   verbatim quote, and its retrieval route.

## How to use this (future agents)

- **Before hunting anything:** `grep -ril "<topic>" research_sources/library/` and read the
  matching entries. The claim table records what we already extracted and which doc consumed it.
- **Before landing a constant:** find it in the reverse index below. If it is not there,
  it is unsourced — add its entry first.
- **Full text beats abstracts.** Entries marked `ABSTRACT-ONLY`, `RELAY` or `QUARANTINED`
  cannot support a landed constant (standing rule from the campaign). Upgrade them by
  fetching, then follow convention 1.
- **Claim rows marked `[unquoted — verify on next use]`** are mentions our docs recorded
  without a verbatim quote. Quote them from the archived full text before use.
- Grep-ready full texts live in `research_sources/thesis_extracts/` (30+ files);
  digitized figure data in `research_sources/digitized/`.

## Index — 111 sources by topic

### Fluorocarbon/SiO₂ surface mechanism (the Krüger–Kushner lineage)

| bibkey | source | status |
|---|---|---|
| [`benck-2003-c4f6`](library/benck-2003-c4f6.md) | Benck, Goyette & Wang, absolute mass-resolved C4F6/Ar ion flux and IEDs | PRIMARY FULL TEXT + FIGURES 9/10 QUANTITATIVELY PIL-AUDITED |
| [`bruce-graves`](library/bruce-graves.md) | Bruce / Graves, ion-dose crosslinking | VIA RESEARCH DOCS |
| [`hiwasa-2022-apex`](library/hiwasa-2022-apex.md) | Hiwasa et al. (KIOXIA), APEX 15, 106002 (2022) | FETCHED |
| [`huang-2019-jvsta`](library/huang-2019-jvsta.md) | Huang, Huard, ... Kushner, JVST A 37, 031304 (2019) | FULL TEXT (via thesis lineage) |
| [`huang-thesis`](library/huang-thesis.md) | Huang, PhD thesis, Univ. Michigan | FULL TEXT: research_sources/thesis_extracts/huang_thesis.txt |
| [`humbird-2004-apl`](library/humbird-2004-apl.md) | Humbird, Graves, Hua & Oehrlein, APL 84, 1073 (2004) | PRIMARY FULL TEXT ONLINE + VERIFIED EXCERPT |
| [`huard-thesis`](library/huard-thesis.md) | Huard, PhD thesis, Univ. Michigan | FULL TEXT: research_sources/thesis_extracts/huard_chad_phd_thesis.txt |
| [`izawa-2007-jjap`](library/izawa-2007-jjap.md) | Izawa et al., JJAP 46, 7870 (2007) | ABSTRACT-ONLY — NOT IMPORTABLE (model-inverted) |
| [`kim-2021-coatings`](library/kim-2021-coatings.md) | Kim et al., measured C4F6/Ar neutral + ion mass/energy spectra | FULL TEXT + FIGURES VISUALLY AUDITED |
| [`krueger-2024-jvsta`](library/krueger-2024-jvsta.md) | Krüger et al., JVST A 42, 043008 (2024) | FULL TEXT: research_sources/thesis_extracts/krueger-2024.txt, kr2024_osti.txt |
| [`krueger-2024-thesis`](library/krueger-2024-thesis.md) | Krüger, PhD thesis, Univ. Michigan (2024) | FULL TEXT: research_sources/thesis_extracts/krueger_thesis_2024.txt (469k), krueger_thesis.txt (older OCR) |
| [`li-2002-c4f6-c4f8`](library/li-2002-c4f6-c4f8.md) | Li et al., matched C4F6/Ar versus C4F8/Ar reactor/surface board | PRIMARY FULL TEXT ONLINE; FIGURES NOT YET DIGITIZED |
| [`omura-2019-jjap`](library/omura-2019-jjap.md) | Omura et al., JJAP (2019) | ABSTRACT/relay |
| [`metzler-2016-thesis`](library/metzler-2016-thesis.md) | Metzler, UMD thesis (2016), cyclic C4F8/Ar-ion ALE | FULL TEXT + PIL-AUDITED FIGURES |
| [`standaert-oehrlein`](library/standaert-oehrlein.md) | Standaert/Oehrlein, mixed-layer selectivity | OFFICIAL OSTI FULL REPORTS + VERIFIED EXCERPTS |
| [`wang-thesis`](library/wang-thesis.md) | Wang (Mingmei), PhD thesis, Univ. Michigan | FULL TEXT: research_sources/thesis_extracts/wang_mingmei_phd_thesis.txt |
| [`woo-2024-c4f6-thesis`](library/woo-2024-c4f6-thesis.md) | Woo, CF4/C4F6/He ICP diagnostics, patterned rates, and SEM profiles | PRIMARY FULL TEXT + PIL-DIGITIZED FIGURE 4.1 |

### Beam-measured yields, thresholds, sticking (the provenance floor)

| bibkey | source | status |
|---|---|---|
| [`chae-2003-jvsta`](library/chae-2003-jvsta.md) | Chae, Vitale & Sawin (2003) | PRIMARY FULL TEXT ONLINE + VERIFIED EXCERPT |
| [`chang-1997-jvsta`](library/chang-1997-jvsta.md) | Chang & Sawin, JVST A 15, 610 (1997) | ABSTRACT-ONLY (curve digitized from thesis p.115) |
| [`chang-thesis`](library/chang-thesis.md) | Chang, PhD thesis, MIT (1721.1/50356) | FULL TEXT: research_sources/thesis_extracts/chang_thesis.txt |
| [`gray-1993-thesis`](library/gray-1993-thesis.md) | Gray, PhD thesis, MIT (1993) | FULL TEXT (OCR sections): research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt |
| [`guo-thesis`](library/guo-thesis.md) | Guo, PhD thesis, MIT (1721.1/46600) | FETCHED + PIL-AUDITED; PDF not redistributed |
| [`joubert-1994-jvsta`](library/joubert-1994-jvsta.md) | Joubert, Oehrlein & Surendra, JVST A 12, 665 (1994) | VIA CHANG THESIS (body [VERIFY]) |
| [`karahashi-2007`](library/karahashi-2007.md) | Karahashi, ion-beam SiO2/CFx | FULL TEXT: research_sources/thesis_extracts/karahashi_2007_sio2_cfx_ionbeam.txt |
| [`kwon-2006-jvsta`](library/kwon-2006-jvsta.md) | Kwon et al., JVST A 24, 1906/1914/1920 (2006) | PARTIAL (Fig 3.4 replot of Gray 1993) |
| [`levinson-llnl`](library/levinson-llnl.md) | Levinson, Shaqfeh, Balooch & Hamza controlled chlorine-beam corpus | PRIMARY FULL TEXT ONLINE; FEATURE TIME/FLUENCE AND ORIGINAL PIXELS MISSING |
| [`mahorowala-1998-thesis`](library/mahorowala-1998-thesis.md) | Mahorowala, fixed-time Cl2/HBr poly-Si feature corpus and simulator | PRIMARY FULL THESIS READ + TABLE 2.2/FIGURE 2.4 PIL-AUDITED |
| [`mahorowala-2002`](library/mahorowala-2002.md) | Mahorowala & Sawin, JVST B 20, 1055/1077 (2002) | ABSTRACT/relay |
| [`steinbruchel`](library/steinbruchel.md) | Steinbrüchel, sqrt-E yield form | VIA ARM + THESES |
| [`takada-2005-tmrsj`](library/takada-2005-tmrsj.md) | Takada, Toyoda & Sugai, stable C5F8 molecule / Ar+ co-incidence on SiO2 | FULL TEXT + PIL-AUDITED FIGURE 3 |
| [`tachi-1982-jjaps`](library/tachi-1982-jjaps.md) | Tachi et al., mass-selected fluorocarbon ions on elemental Si through 3 keV | ABSTRACT-ONLY — TARGET-MISMATCHED LEAD, NOT IMPORTABLE |
| [`yin-2008-jvsta`](library/yin-2008-jvsta.md) | Yin & Sawin, JVST A 26, 161 (2008) | not-fetched |

### Atomistic surface physics and event kernels

| bibkey | source | status |
|---|---|---|
| [`an-2026-apsusc`](library/an-2026-apsusc.md) | An et al., DFT-trained NNP/MD for HFC ions on SiO2 and Si3N4 | FULL TEXT READ + RELEASED DATA PINNED |
| [`cagomoc-2023-thesis`](library/cagomoc-2023-thesis.md) | Cagomoc, Osaka thesis, CF3+/radical SiO2 MD and nanohole product escape | FULL TEXT + PIL-AUDITED FIGURES |
| [`tinacba-2021-jvstb`](library/tinacba-2021-jvstb.md) | Tinacba et al., DFT-informed SF5+ MD versus mass-selected Si/SiO2 beam depths | FULL TEXT READ + PIL-AUDITED FIGURE 8 |

### Angular yield laws (class-1 physical, class-2 chemical)

| bibkey | source | status |
|---|---|---|
| [`arts-2021-apr`](library/arts-2021-apr.md) | Arts et al., atomic-level plasma-processing review / angular provenance audit | FULL TEXT READ + VERIFIED EXCERPT |
| [`barklund-1992-jvsta`](library/barklund-1992-jvsta.md) | Barklund & Blom, JVST A 10, 1212 (1992) | VIA CHANG THESIS FIG 4.16 |
| [`cho-2000-jvsta`](library/cho-2000-jvsta.md) | Cho et al. (2000) | ABSTRACT-ONLY |
| [`kress-1999-jvsta`](library/kress-1999-jvsta.md) | Kress et al., JVST A 17, 2819 (1999) | ABSTRACT-ONLY (title verified: Cu/Ar MD, wrong system) |
| [`lee-2002-jvsta`](library/lee-2002-jvsta.md) | Lee et al. (2002), grazing-angle etch rates | ABSTRACT/relay |
| [`mayer-1981`](library/mayer-1981.md) | Mayer (1981), SiO2 + CFx+ angular yield | VIA CHANG THESIS FIG 4.7 |
| [`schaepkens-1998-jvsta`](library/schaepkens-1998-jvsta.md) | Schaepkens et al., JVST A 16, 3281 (1998) | ABSTRACT-ONLY |
| [`you-2023-coatings`](library/you-2023-coatings.md) | You et al., Coatings 13, 1452 (2023) | FULL TEXT: research_sources/thesis_extracts/coatings2023_bowing_narrowing.txt |

### Ion angular/energy distributions and sheath collisions

| bibkey | source | status |
|---|---|---|
| [`cunge-2016-apl`](library/cunge-2016-apl.md) | Cunge et al., APL 108, 093109 (2016) | ABSTRACT-ONLY |
| [`kawamura-2025-psst`](library/kawamura-2025-psst.md) | Kawamura et al. (Nagoya+KIOXIA), PSST 34, 055006 (2025) | ABSTRACT-ONLY |
| [`khrabrov-2026-arxiv`](library/khrabrov-2026-arxiv.md) | Khrabrov & Kaganovich, arXiv:2604.04214 (2026) | FETCHED (arXiv) |
| [`kim-2025-jjap-05sp15`](library/kim-2025-jjap-05sp15.md) | Kim et al., JJAP 64, 05SP15 (2025) | ABSTRACT/relay [VERIFY verbatim] |
| [`kim-2025-jjap-096002`](library/kim-2025-jjap-096002.md) | Kim et al., JJAP 64, 096002 (2025) | ABSTRACT-ONLY |

### Reactor-scale and sheath closure models

| bibkey | source | status |
|---|---|---|
| [`edelberg-1999`](library/edelberg-1999.md) | Edelberg & Aydil (1999) | not-fetched |
| [`kawamura-1999-psst`](library/kawamura-1999-psst.md) | Kawamura et al., PSST 8, R45 (1999) | not-fetched (citation corrected from p.313) |
| [`benyoucef-yousfi-2014-ion-transport`](library/benyoucef-yousfi-2014-ion-transport.md) | Benyoucef & Yousfi, semiclassical Ar+/Ar, O2+/O2, N2+/N2 transport validation | AUTHOR-PROVIDED FULL TEXT READ VIA HTML; FIGURES NOT LOCALLY ARCHIVED |
| [`lee-lieberman-1994-global`](library/lee-lieberman-1994-global.md) | Lee & Lieberman, conserved Ar/O2/Cl2 global plasma model and argon rate deck | PRIMARY FULL TEXT + EQUATIONS/TABLE 3 VISUALLY AUDITED |
| [`miller-1997`](library/miller-1997.md) | Miller & Riley (1997) sheath model | not-fetched |
| [`nist-asd-argon`](library/nist-asd-argon.md) | NIST ASD neutral-argon ionization and 4s metastable energies | PRIMARY NIST DATABASE QUERIES |
| [`nist-tn-2279-gas-diffusion`](library/nist-tn-2279-gas-diffusion.md) | NIST evaluated gas self-/binary-diffusion correlations | AUTHORITATIVE FULL PDF READ |
| [`phelps-1994-ar-ion-scattering`](library/phelps-1994-ar-ion-scattering.md) | Phelps, consistent Ar+-Ar momentum-transfer/scattering model | PRIMARY PUBLISHER RECORD; EQUATION CROSS-CHECKED |
| [`raja-linne`](library/raja-linne.md) | Raja & Linne | not-fetched (DOI corrected from 1.1519941) |

### Transport references and analytic benchmarks

| bibkey | source | status |
|---|---|---|
| [`clausing`](library/clausing.md) | Clausing, transmission probability | ANALYTIC REFERENCE |
| [`coburn-winters`](library/coburn-winters.md) | Coburn & Winters, conductance/ARDE | VIA THESES |
| [`santeler`](library/santeler.md) | Santeler, transmission probability accuracy | ANALYTIC REFERENCE |
| [`zbl`](library/zbl.md) | Ziegler-Biersack-Littmark stopping | IMPLEMENTED FORM |

### Feature charging, notching, electron shading

| bibkey | source | status |
|---|---|---|
| [`fujiwara-notching`](library/fujiwara-notching.md) | Fujiwara et al., JJAP 34, 2095 / 35, 2450 | ABSTRACT-ONLY |
| [`hashimoto-shading`](library/hashimoto-shading.md) | Hashimoto, electron shading theory | not-fetched |
| [`huang-2026-jvsta-charging`](library/huang-2026-jvsta-charging.md) | Huang & Kushner, JVST A 44, 023013 (2026) | FULL TEXT READ |
| [`jinnai-2007`](library/jinnai-2007.md) | Jinnai et al., JVST B 25, 1808 (2007) | ABSTRACT-ONLY |
| [`kamata-arimoto`](library/kamata-arimoto.md) | Kamata & Arimoto, JAP 80, 2637 / JVST B 14, 3688 | ABSTRACT-ONLY |
| [`krueger-2019-psst`](library/krueger-2019-psst.md) | Krüger & Schulze, PSST 28, 075017 (2019) | not-fetched |
| [`matsui-makabe`](library/matsui-makabe.md) | Matsui & Makabe, APL 78, 883 | ABSTRACT-ONLY (negative control) |
| [`nozawa-1995`](library/nozawa-1995.md) | Nozawa et al., JJAP 34, 2107 (1995) | ABSTRACT-ONLY |
| [`ohtake-2007`](library/ohtake-2007.md) | Ohtake et al., JVST B 25, 400 (2007) | ABSTRACT-ONLY |
| [`shimmura-2004`](library/shimmura-2004.md) | Shimmura et al., JVST A 22, 433 (2004) | ABSTRACT-ONLY |
| [`wang-2010-kushner`](library/wang-2010-kushner.md) | Wang & Kushner (2010) twisting | READ (via research pass) |

### HARC / extreme-AR field practice and ARDE measurements

| bibkey | source | status |
|---|---|---|
| [`gottscho-arde`](library/gottscho-arde.md) | Gottscho, Jurgensen & Vitkavage, ARDE review | ABSTRACT/relay |
| [`ishikawa-2018-jjap`](library/ishikawa-2018-jjap.md) | Ishikawa et al., JJAP 57, 06JA01 (2018) | FETCHED |
| [`kim-2007-tsf-lam`](library/kim-2007-tsf-lam.md) | Kim, Hudson, Cooperberg, Edelberg, Srinivasan (Lam), Thin Solid Films 515, 4874 (2007) | ABSTRACT-ONLY [top follow-up] |
| [`lam-shen-2023-jjap`](library/lam-shen-2023-jjap.md) | Shen, Lill et al. (Lam), JJAP 62, SI0801 (2023) | FULL TEXT: research_sources/thesis_extracts/lam_shen_lill_jjap2023.txt |
| [`lee-2010-jes`](library/lee-2010-jes.md) | Lee et al., JES 157, D142 (2010) | ABSTRACT/relay |
| [`maruyama-2010`](library/maruyama-2010.md) | Maruyama et al., JVST B 28, 862 (2010) | ABSTRACT-ONLY |
| [`nguyen-2020-jvsta`](library/nguyen-2020-jvsta.md) | Nguyen, ... Jansen, JVST A 38, 053002 (2020) | FETCHED (rate table derived) |
| [`nishizuka-2024-jvsta`](library/nishizuka-2024-jvsta.md) | Nishizuka et al. (TEL), JVST A 42, 043003 (2024) | ABSTRACT-ONLY |
| [`ohiwa-1998`](library/ohiwa-1998.md) | Ohiwa et al., JJAP 37, 5060 (1998) | ABSTRACT-ONLY |
| [`samukawa-nbe`](library/samukawa-nbe.md) | Samukawa lineage, neutral-beam etching | ABSTRACT/relay |
| [`vanraes-2023-psst`](library/vanraes-2023-psst.md) | Vanraes et al., PSST 32, 065003 (2023) | ABSTRACT-ONLY |

### Modeling state of the art and competitor codes

| bibkey | source | status |
|---|---|---|
| [`hoekstra-2002-jvstb`](library/hoekstra-2002-jvstb.md) | Hoekstra & Kushner, JVST B 20, 1077 (2002) | ABSTRACT/relay |
| [`kokkoris`](library/kokkoris.md) | Kokkoris et al., ARDE/profile modeling | ABSTRACT/relay |
| [`kuboi-2024-jjap`](library/kuboi-2024-jjap.md) | Kuboi, JJAP 63, 080801 (2024) | FETCHED |
| [`rodrigues-2023`](library/rodrigues-2023.md) | Rodrigues et al. (TU Wien), FC/silica (2023) | FULL TEXT: research_sources/thesis_extracts/tuwien_rodrigues_2023_fc_silica.txt |
| [`yook-2022-jphysd`](library/yook-2022-jphysd.md) | Yook/Im (K-SPEED), J. Phys. D 55, 255202 (2022) | ABSTRACT-ONLY (K-SPEED: no charging module — prior-art claim refuted) |
| [`zhai-2025-jap`](library/zhai-2025-jap.md) | Zhai, Filipović, Chen, JAP 137, 063302 (2025) | not-fetched |
| [`zhang-thesis`](library/zhang-thesis.md) | Zhang (Yiting), PhD thesis, Univ. Michigan (2015) | FULL TEXT: research_sources/thesis_extracts/zhang_yiting_phd_thesis.txt |

### Atomic-F density / flux measurements (the supply band)

| bibkey | source | status |
|---|---|---|
| [`chun-2015-tsf`](library/chun-2015-tsf.md) | Chun, Efremov, Yeom & Kwon, Thin Solid Films 579, 136 (2015) | FULL TEXT READ |
| [`cunge-2001-jap`](library/cunge-2001-jap.md) | Cunge, Chabert & Booth, JAP 89, 7750 (2001) | QUARANTINED — abstract only, no numeric [F] |
| [`jenq-1994-psst`](library/jenq-1994-psst.md) | Jenq, Ding, Taylor & Hershkowitz, PSST 3, 154 (1994) | RELAY [Q-relay][VERIFY] |
| [`kawai-1997-jjap`](library/kawai-1997-jjap.md) | Kawai, Sasaki & Kadota, JJAP 36, L1261 (1997) | RELAY [Q-relay][VERIFY] |
| [`sankaran-2005-jap`](library/sankaran-2005-jap.md) | Sankaran & Kushner, JAP 97, 023307 (2005) | FULL TEXT (Table II) |
| [`sasaki-1997-jap`](library/sasaki-1997-jap.md) | Sasaki et al., JAP 82, 5938 (1997) | ABSTRACT-ONLY |

### Line-edge roughness: metrology, transfer, experimental gates

| bibkey | source | status |
|---|---|---|
| [`azarnouche-thesis`](library/azarnouche-thesis.md) | Azarnouche, PhD thesis LTM-CNRS (tel-00767820) | FETCHED (Fig IV.39 = measured |T|^2) |
| [`constantoudis`](library/constantoudis.md) | Constantoudis, Kokkoris & Gogolides, LER transfer | ABSTRACT/relay |
| [`dupuy-2015`](library/dupuy-2015.md) | Dupuy et al. (2015) SADP LWR | FETCHED (demoted: PSDs shifted to 1) |
| [`kushner-2021-jvsta-ler`](library/kushner-2021-jvsta-ler.md) | Kushner, JVST A 39, 033003 (2021) | FETCHED (Poisson anchor) |
| [`liang-2018`](library/liang-2018.md) | Liang, Mack (Lam) 2018 LER | REJECTED as gate (knob anonymized, axes normalized) |
| [`mack-ler`](library/mack-ler.md) | Mack, LER metrology / noise floor | VIA LER DOCS |
| [`martin-cunge`](library/martin-cunge.md) | Martin & Cunge, plasma smoothing | ABSTRACT/relay |
| [`palasantzas`](library/palasantzas.md) | Palasantzas, self-affine PSD form | VIA LER DOCS |
| [`pargon-2013-jvstb`](library/pargon-2013-jvstb.md) | Pargon et al. (LTM), JVST B 31, 012205 (2013) | ABSTRACT/relay |
| [`rutigliani-2018-spie`](library/rutigliani-2018-spie.md) | Rutigliani, Lorusso, De Simone & Mack, Proc. SPIE 10585 (2018) | FETCHED (absolute nm^3 PSDs) |

### SF₆/C₄F₈ and SF₆/O₂ on silicon (the partner-relevant arm)

| bibkey | source | status |
|---|---|---|
| [`belen-2005-jvsta`](library/belen-2005-jvsta.md) | Belen et al., JVST A 23, 99 (2005) | ABSTRACT (self-declares L3 profile-fitted) |
| [`deboer-2002`](library/deboer-2002.md) | de Boer et al. (2002), cryo SF6/O2 Si | FULL TEXT: research_sources/thesis_extracts/deboer-2002.txt |
| [`micromachines-2023`](library/micromachines-2023.md) | TU Wien/ViennaPS, Micromachines (2023) | FULL TEXT: research_sources/thesis_extracts/mask_geometry_micromachines_2023.txt |
| [`yoshie-2023-apsusc`](library/yoshie-2023-apsusc.md) | Yoshie et al., cyclic C4F8/SF6 Si bias timing and ARDE | FULL TEXT ONLINE + PIL-AUDITED FIGURES |

### In-feature and profile metrology

| bibkey | source | status |
|---|---|---|
| [`okawa-2026-apl`](library/okawa-2026-apl.md) | Okawa/Takahashi, APL 129, 011109 (2026) | ABSTRACT-ONLY |

## Reverse index — petch constant / law / decision → source

| petch constant, law or decision | source | what it fixed / why it is there |
|---|---|---|
| Conserved global-model power/particle balance and cylindrical effective-loss area | [`lee-lieberman-1994-global`](library/lee-lieberman-1994-global.md) | chemistry-agnostic 0-D reactor skeleton; electrons explicit so plasma shorthand cannot hide charge error |
| Argon Table-3 electron-impact rates and Bohm/diffusive wall-loss forms | [`lee-lieberman-1994-global`](library/lee-lieberman-1994-global.md) | first no-fit reactor verification deck; empirical/regressed evidence retained, not called first-principles |
| Argon ionization/metastable event energies `15.7596119 / 11.54835442 / 4.21125748 eV` | [`nist-asd-argon`](library/nist-asd-argon.md) | physical power ledger kept separate from `18.68 / 15.06 / 4.95 eV` rate-fit exponents |
| Ar+-Ar `Qm(E)` with explicit center-of-mass energy and Maxwellian rate averaging | [`phelps-1994-ar-ion-scattering`](library/phelps-1994-ar-ion-scattering.md) | replaces the fixed ion mean-free-path sensitivity with an energy- and pressure-dependent collision frequency |
| Ar-in-Ar `D(T,p)` and 298 K reference `0.182 cm2/s` | [`nist-tn-2279-gas-diffusion`](library/nist-tn-2279-gas-diffusion.md) | supplies neutral bulk diffusion before Lee's Knudsen harmonic; 600 K use is explicitly flagged as extrapolated |
| Complex-channel sputter yield `Y = 0.0139(√E − √18)` | [`gray-1993-thesis`](library/gray-1993-thesis.md) | measured absolute law; replaced Krüger's 2-row linear/Sigmund anomaly (RESULTS_GRAY_ANCHORING) |
| Chemical branching `β = 0.053(√E − √4)` | [`gray-1993-thesis`](library/gray-1993-thesis.md) | Table 5-10 parenthetical B₀ column; co-regressed partner of s₀ |
| Oxide F sticking `s₀ = 0.02` | [`gray-1993-thesis`](library/gray-1993-thesis.md) | printed p.246; landed only as the (s₀,B₀) pair |
| Silicon F sticking `s₀ = 0.2` + `B₀ = 0.687(√E − √4)` | [`gray-1993-thesis`](library/gray-1993-thesis.md) | Si-side pair for the SF₆ arm upgrade (RESEARCH_SF6_RELEVANCE) |
| Cyclic C4F8/Ar-ion Si/SiO2 Å-per-cycle + ion-normalized yield + XPS state board | [`metzler-2016-thesis`](library/metzler-2016-thesis.md) | 42 depth + 25 cycle-normalized yield + 32 XPS markers; finite transfer-depth constraint, retrospective not blind |
| Measured mixed-layer depth: C/F/O/Si to ~15 Å, F/O tail to ~25 Å | [`metzler-2016-thesis`](library/metzler-2016-thesis.md) | direct SIMS constraint on the silicon-under-film state |
| Ar+-induced C-F destruction → Si-F formation, nearly 1:1 | [`humbird-2004-apl`](library/humbird-2004-apl.md) | atom-balanced film-to-substrate transfer channel |
| Energy-dependent Si-C-F depth (~15 Å at 20 eV; ~30 Å at 200 eV) | [`humbird-2004-apl`](library/humbird-2004-apl.md) | atomistic constraint on mixed-layer capacity/depth |
| Neutral radical sticking: CF 0.4; CF2 0.01-0.07 | [`standaert-oehrlein`](library/standaert-oehrlein.md) | species-resolved film deposition; prevents one precursor-sticking constant |
| CF2 afterglow wall sticking `(8±6)e-3` at 298 K | [`standaert-oehrlein`](library/standaert-oehrlein.md) | unbiased-wall cross-check only, not ion-activated wafer sticking |
| Exponential Si-yield attenuation vs measured FC-film thickness | [`standaert-oehrlein`](library/standaert-oehrlein.md) | measured film-shielding topology; no decay length invented from the report text |
| C2HF5 plasma SiO2 yield/deposition topology | [`chae-2003-jvsta`](library/chae-2003-jvsta.md) | QCM evidence that reactive CFx+ identity and ion-neutral balance matter; fitted site yield not importable |
| Cyclic C4F8/SF6 blanket/feature depth board | [`yoshie-2023-apsusc`](library/yoshie-2023-apsusc.md) | 7 same-reactor poly-Si blanket rates + 49 held-out bulk-Si feature rates; material/history transfer warning |
| Beam gate N1 (floor 0.25 / half-rise / plateau) | [`gray-1993-thesis`](library/gray-1993-thesis.md) | the measured yield-vs-F/Ar⁺ curve petch is graded against |
| Species-resolved pure-reactive-ion yield ladder | [`karahashi-2007`](library/karahashi-2007.md) | direct F+/CF+/CF2+/CF3+ constraint; forbids species-agnostic validation and extrapolation above measured support |
| DFT-trained HFC-ion event physics + product-escape defect | [`an-2026-apsusc`](library/an-2026-apsusc.md) | no-yield-fit atomistic transfer candidate; forces explicit mixed-layer/film and product-escape closure rather than a scalar depth fit |
| DFT-informed SF5+ atomistic depth-per-dose provider | [`tinacba-2021-jvstb`](library/tinacba-2021-jvstb.md) | four independent mass-selected MD/beam overlaps on Si and SiO2; 5.88% mean, 15.04% maximum error; sulfur-surface chemistry omission retained |
| Condition-unknown CFx+/Ar+ angular species board | [`karahashi-2007`](library/karahashi-2007.md), [`arts-2021-apr`](library/arts-2021-apr.md) | 20 PIL-audited markers; normal points strongly imply 1000 eV, but the source/review do not report energy, so production use is forbidden |
| Stable-parent molecule / ion co-incidence envelope | [`takada-2005-tmrsj`](library/takada-2005-tmrsj.md) | C5F8 analog proves molecule-assisted yields can exceed a pure-CF3+ value; explicitly not a C4F6 law |
| C4F6 parent-signal + ion-mixture existence | [`kim-2021-coatings`](library/kim-2021-coatings.md) | wafer-facing mass/energy spectra directly show C4F6 and multiple CFx+/CxFy+ species in another CCP; qualitative existence constraint only, not Krüger flux calibration |
| Absolute C4F6/Ar total and mass-resolved ion flux/IED board | [`benck-2003-c4f6`](library/benck-2003-c4f6.md) | Faraday-cup-normalized GEC-ICP mixture and pressure boards; quantitative reactor-model validation, not a Krüger boundary transplant |
| Matched C4F6/Ar versus C4F8/Ar ion/radical/film/yield board | [`li-2002-c4f6-c4f8`](library/li-2002-c4f6-c4f8.md) | same 600 W ICP and diagnostics expose feed-specific fragmentation and surface-film differences; figures require original-pixel digitization before constants land |
| CF4/C4F6/He patterned SiO2/ACL rate + reactor-diagnostic board | [`woo-2024-c4f6-thesis`](library/woo-2024-c4f6-thesis.md) | 10 original-pixel rate points plus Te/current/self-bias/OES/profile constraints; two source-reporting conflicts quarantined and equal-depth SEM timing excluded from blind depth validation |
| Class-1 angular form `(1+B sin²θ)cosθ` | [`kress-1999-jvsta`](library/kress-1999-jvsta.md) | Krüger's cited source — **wrong system (Cu/Ar MD)**; peak 4.17 vs measured ~1.3 |
| Class-1 bound `B = 1.7` on oxide/mask rows | [`cho-2000-jvsta`](library/cho-2000-jvsta.md) | in-chemistry peak/normal 1.30 |
| Class-1 bound cross-check | [`schaepkens-1998-jvsta`](library/schaepkens-1998-jvsta.md) | peak/normal 1.33, 54.7° V-groove |
| FC-film angular curve (peak 1.448 @ 65°) | [`barklund-1992-jvsta`](library/barklund-1992-jvsta.md) | the only measured FC-*film* angular yield; selects the yield reading (eliminated 48/64 models) |
| Class-2 chemical roll-off `min(1, cosθ/cos45°)` | [`chang-1997-jvsta`](library/chang-1997-jvsta.md) | unity to 45° then monotone roll-off; digitized curve within 0.065 absolute |
| Class-2 digitized curve + thresholds | [`chang-thesis`](library/chang-thesis.md) | p.115 rendered at 400 dpi; Table 3.1 sputter thresholds |
| Ar+/Cl and Ar+/Cl2 poly-Si site balance at 100 eV | [`chang-thesis`](library/chang-thesis.md) | Tables 3.3--3.4 and Eqs. 3.9--3.11; beam-regressed `s`, `Y0`, and `beta`, with explicit SiCl4 simplification |
| SiO₂-complex threshold 35 eV | [`joubert-1994-jvsta`](library/joubert-1994-jvsta.md) | provenance corrected from 'Chang–Sawin' (via Chang p.90) |
| Reflection cascade Eq. 2.34 (E_ts=100, E_c=10, θ_c=70°) | [`huang-thesis`](library/huang-thesis.md) | verbatim retention rule + leftover-probability selection |
| Crosslink formation = deposition-driven (ion-broken) | [`krueger-2024-thesis`](library/krueger-2024-thesis.md) | §2.2.3 + Table 6.2 — the lip inversion fix (ml18: 5.07×→1.83×) |
| Crosslink bond multiplicity (CF:3, CF₂:2, CF₃:1) | [`krueger-2024-jvsta`](library/krueger-2024-jvsta.md) | sec. III verbatim — closed the campaign's last [VERIFY] integer |
| All Appendix-B reaction probabilities + Table 6.5 converged set | [`krueger-2024-thesis`](library/krueger-2024-thesis.md) | the deck's chemistry rows |
| Site-turnover asymmetry (~140× re-passivation drop) | [`huang-thesis`](library/huang-thesis.md) | L10214-10222; independently in Krüger L6556-6564 |
| Finite-site-turnover regime statement | [`huard-thesis`](library/huard-thesis.md) | Gottscho quote L3763-3770 — why validated models are neutral-limited |
| Two-component IADF (core + collisional tail) | [`kim-2025-jjap-05sp15`](library/kim-2025-jjap-05sp15.md) | measured core 0.044 eV / tail 0.57 eV at 0.1° resolution |
| Tail-fraction pressure law (route to retiring the sweep) | [`kim-2025-jjap-096002`](library/kim-2025-jjap-096002.md) | main/tail ratio falls exponentially with pressure |
| Sheath-collision cross sections (anisotropic Born–Mayer) | [`khrabrov-2026-arxiv`](library/khrabrov-2026-arxiv.md) | first-principles origin of the tail; AR→tolerance law (0.25° lab frame @ AR 100) |
| Transport benchmark at AR 200 (0.656%, err ≤0.7%) | [`clausing`](library/clausing.md) | analytic transmission probability |
| Transport accuracy band | [`santeler`](library/santeler.md) | ≤0.7% error envelope petch sits inside at every AR |
| Deposited-energy factor (ZBL shape) | [`zbl`](library/zbl.md) | the DEKNOB-derived energy law that retired a fitted knob |
| AR-50 delivery cross-check (2.5%) | [`lam-shen-2023-jjap`](library/lam-shen-2023-jjap.md) | field's published number; petch measures 0.025287 |
| ARDE band 43–80% (our 60% @ AR 50 sits inside) | [`nguyen-2020-jvsta`](library/nguyen-2020-jvsta.md) | deepest published *measurement*, ~43% by AR 54 |
| ARDE upper anchor (80% by AR 40) | [`huang-2019-jvsta`](library/huang-2019-jvsta.md) | simulated; paired with Nguyen as the band |
| Charging potential ceiling gate D1 (e·V_max ≤ E_ion,max) | [`huang-2026-jvsta-charging`](library/huang-2026-jvsta-charging.md) | 9 published points, near AR-independent above AR 17 |
| Charging twisting ensemble ladder D4 | [`wang-2010-kushner`](library/wang-2010-kushner.md) | 12%→49%→38%→25%→12% ablation ordering |
| Electron-shading gate D2 (100 V @ 400 V bias, AR 2) | [`kamata-arimoto`](library/kamata-arimoto.md) | measured current-through-dielectric |
| Notch-vs-AR + pulsed rescue gate D3 | [`fujiwara-notching`](library/fujiwara-notching.md) | monotone notch(AR), knee at AR≈2 |
| Negative control: must NOT reproduce AR>7 etch stop | [`matsui-makabe`](library/matsui-makabe.md) | refuted by manufacturing reality |
| de Boer charging claim retired (redeposition owns etch stop) | [`ohiwa-1998`](library/ohiwa-1998.md) | measured HARC etch stop attributed to redeposition |
| Atomic-F supply band Γ_F = 2×10²⁰–1×10²¹ m⁻²s⁻¹ | [`sankaran-2005-jap`](library/sankaran-2005-jap.md) | HPEM-published F flux anchor (Table II) |
| F-band ceiling (measured densities) | [`chun-2015-tsf`](library/chun-2015-tsf.md) | C₄F₈/O₂ inverts the O₂-raises-F premise (−9.3×) |
| F-band relay bounds | [`jenq-1994-psst`](library/jenq-1994-psst.md) | [Q-relay][VERIFY] RIE CF₄ densities |
| F wall-loss probability ~1e-3 | [`sasaki-1997-jap`](library/sasaki-1997-jap.md) | brackets Krüger's 0.01 polymer rows |
| Sidewall CFx sticking 0.004 — **NOT importable** | [`izawa-2007-jjap`](library/izawa-2007-jjap.md) | model-inverted effective quantity spanning a 125× F-rich/C-rich axis |
| Grazing NER above cosine (30–60°) — does NOT apply at our 86–89° lip | [`you-2023-coatings`](library/you-2023-coatings.md) | Faraday-cage measured; net deposition >80° |
| LER PSD synthesis / σ,ξ,α | [`palasantzas`](library/palasantzas.md) | self-affine form in ler_metrology |
| LER noise floor subtraction | [`mack-ler`](library/mack-ler.md) | Mack floor in ler_metrology |
| LER Gate 2 (σ-transfer slope 0.5, threshold form) | [`constantoudis`](library/constantoudis.md) | static shadowing gives 1.0 — bounded negative, erosion is the missing term |
| LER Gate 1 (measured |T(k)|² + null control) | [`azarnouche-thesis`](library/azarnouche-thesis.md) | Fig. IV.39 PSD ratio; VUV-cured resist = free null |
| LER intrinsic-noise Poisson anchor | [`kushner-2021-jvsta-ler`](library/kushner-2021-jvsta-ler.md) | ~14 radicals/site → 25%, ~10 ions/site → 30% = 1/√N |
| LER absolute normalization check | [`rutigliani-2018-spie`](library/rutigliani-2018-spie.md) | 16 nm CD EUV, absolute nm³ PSDs |
| SF₆/O₂ arm constants (declared L3 profile-fitted) | [`belen-2005-jvsta`](library/belen-2005-jvsta.md) | paper's own abstract: parameters fitted to feature profiles |
| SF₆/O₂ cryo validation case | [`deboer-2002`](library/deboer-2002.md) | transfer test, not blind; cryo≠Gray room-temp caveat |
| √E yield form in the SF₆ arm | [`steinbruchel`](library/steinbruchel.md) | already-correct structural choice in belen.py |
| γ_F = 0.7 default (same ViennaPS default, not corroboration) | [`micromachines-2023`](library/micromachines-2023.md) | citation corrected in RESEARCH_SF6_RELEVANCE |
| Ion-dose crosslinking (zero-knob dose law) | [`bruce-graves`](library/bruce-graves.md) | ion-processed-skin dose competition |
| Selectivity from lattice oxygen | [`standaert-oehrlein`](library/standaert-oehrlein.md) | mixed-layer selectivity mechanism |
| Neutral conductance / ARDE framing | [`coburn-winters`](library/coburn-winters.md) | Coburn–Winters conductance used in the F-delivery estimate |

## Quarantine / do-not-use-as-evidence

| source | why |
|---|---|
| [`izawa-2007-jjap`](library/izawa-2007-jjap.md) | 0.004 sticking is a *model-inverted* effective coefficient spanning a 125× F-rich/C-rich axis — ordering constraint only |
| [`kress-1999-jvsta`](library/kress-1999-jvsta.md) | Krüger's cited class-1 source is Cu/Ar MD — wrong system; retained only where Krüger's own fit depends on it |
| [`liang-2018`](library/liang-2018.md) | knob anonymized, axes normalized — rejected as an LER gate |
| [`cunge-2001-jap`](library/cunge-2001-jap.md) | right method and reactor class, but body unobtainable (403 + Anubis) — no numeric [F] read |
| [`yook-2022-jphysd`](library/yook-2022-jphysd.md) | K-SPEED prior-art claim **refuted** — no charging/Poisson module in the publication |
