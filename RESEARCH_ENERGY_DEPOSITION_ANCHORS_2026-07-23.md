# Energy-Deposition Physics: Literature Anchors for petch
Date: 2026-07-23. Purpose: replace fitted yield-saturation knobs with derived
ZBL/Lindhard stopping + Sigmund deposition-in-layer yields + radiation-chemistry
crosslinking. All numbers below are citable; items marked [VERIFY] need a
page-level check before being quoted in a paper.

---

## 1. Stopping / range: formula constants and validation anchors

### 1.1 ZBL universal nuclear stopping (implementation-verification constants)

Primary source: J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and
Range of Ions in Solids*, Vol. 1, Pergamon Press, New York (1985) ("ZBL85").
Modern summary: J. F. Ziegler, M. D. Ziegler, J. P. Biersack, "SRIM — The
stopping and range of ions in matter (2010)", *Nucl. Instrum. Methods B* 268,
1818–1823 (2010), DOI: 10.1016/j.nimb.2010.02.091.

**Universal screening function** (x = r/a_U):

    phi(x) = 0.1818 e^(-3.2 x) + 0.5099 e^(-0.9423 x)
           + 0.2802 e^(-0.4029 x) + 0.02817 e^(-0.2016 x)

(4-digit rounding as quoted on Wikipedia "Stopping power (particle radiation)";
ZBL85 full precision: 0.18175/3.19980, 0.50986/0.94229, 0.28022/0.40290,
0.02817/0.20162.)

**Universal screening length** (a0 = Bohr radius = 0.529 Å):

    a_U = 0.8854 a0 / (Z1^0.23 + Z2^0.23)

**Reduced energy** (E in keV, M in amu):

    eps = 32.53 M2 E / [ Z1 Z2 (M1 + M2) (Z1^0.23 + Z2^0.23) ]

**Reduced nuclear stopping** (ZBL85 universal fit; confirmed in arXiv:1203.4620
and the SR-NIEL long write-up, sr-niel.org):

    Sn(eps) = ln(1 + 1.1383 eps) / [ 2 (eps + 0.01321 eps^0.21226 + 0.19593 eps^0.5) ]   for eps <= 30
    Sn(eps) = ln(eps) / (2 eps)                                                          for eps > 30

Fit constants to check in our code: **1.1383, 0.01321, 0.21226, 0.19593, 32.53,
0.8854, 0.23**.

**Conversion to physical units** [VERIFY against ZBL85 Eq. 2-88 before paper use]:

    Sn(E) = 8.462e-15 * Z1 Z2 M1 * Sn(eps) / [ (M1+M2)(Z1^0.23+Z2^0.23) ]   eV cm^2/atom

### 1.2 Lindhard–Scharff (LSS) electronic stopping

Primary sources: J. Lindhard, M. Scharff, "Energy dissipation by ions in the
keV region", *Phys. Rev.* 124, 128 (1961), DOI: 10.1103/PhysRev.124.128;
J. Lindhard, M. Scharff, H. E. Schiøtt, *Mat. Fys. Medd. Dan. Vid. Selsk.* 33,
No. 14 (1963).

Velocity-proportional cross section (per target atom), as implemented e.g. in
LAMMPS (A. Sand et al., "Inclusion and validation of electronic stopping in the
open source LAMMPS code", arXiv:2005.11940, Eq. 3):

    Se(E) = 8 pi e^2 a0 * Z1^(1/6) * Z1 Z2 / (Z1^(2/3) + Z2^(2/3))^(3/2) * (v / v0)

with xi_e ≈ Z1^(1/6) (Lindhard's electronic factor), v0 the Bohr velocity in
the original LS formulation (the LAMMPS paper substitutes a Fermi-velocity
variant — note the difference when cross-checking). Valid for v << v0 Z1^(2/3),
i.e. the entire 0.1–10 keV Ar regime.

Reduced (LSS) form used in range theory:

    Se(eps) = k eps^(1/2),
    k = xi_e * 0.0793 * Z1^(1/2) Z2^(1/2) (A1 + A2)^(3/2)
        / [ (Z1^(2/3) + Z2^(2/3))^(3/4) A1^(3/2) A2^(1/2) ],   xi_e ≈ Z1^(1/6)

Constant to check: **0.0793** (k typically 0.1–0.2 for heavy ions).
[VERIFY exact k expression against LSS 1963 or Nastasi/Mayer/Hirvonen,
*Ion–Solid Interactions*, Cambridge (1996), Ch. 5 — two equivalent groupings of
the Z1 exponents circulate (Z1^(1/6)·Z1^(1/2) = Z1^(2/3)).]

### 1.3 Accuracy expectations, analytic ZBL vs full SRIM at sub-10-keV

- The universal ZBL potential is a fit across pairs: standard deviation of the
  fit is **18% above 2 eV** (ZBL85, quoted on Wikipedia "Stopping power").
  Pair-specific error for Ar on Si/O/C is typically better but O(10%).
- SRIM's overall stopping accuracy claim is ~4% (Ziegler 2010, NIMB 268, 1818),
  but that statistic is dominated by >25 keV/u data; sub-keV heavy-ion data are
  sparse and the electronic partition relies on LS-type velocity scaling.
- For **shallow implants, SRIM systematically overestimates depth by 2–6 nm**
  vs SIMS/RBS, growing with energy: "Comparing SRIM simulations and
  experimental results for shallow implants", OSTI report/paper,
  https://www.osti.gov/servlets/purl/1408290. See also "Evaluation of the
  Accuracy of Stopping and Range of Ions in Matter",
  https://www.osti.gov/servlets/purl/1830531.
- Practical target for petch: analytic ZBL+LS magnetic-formula ranges should
  match SRIM full-cascade Rp within ~10–20% at 0.25–5 keV; discrepancies beyond
  that in compound targets (SiO2, CFx) usually trace to Bragg-rule additivity
  and the assumed target density, not the integrator.

### 1.4 Published Rp / damage-depth anchors, Ar into Si / SiO2 / carbon

A clean published table at exactly 0.25/0.5/1/2/3/5 keV does not exist in the
open literature — SRIM is treated as the source of truth and papers quote spot
values. Citable spot anchors found:

| System | Energy | Quantity | Value | Source |
|---|---|---|---|---|
| Ar -> Si | 0.25 keV | amorphous/damage layer | ~1 nm | MD study, *Beilstein J. Nanotechnol.* 13 (2022), art. 86 (water-contamination sputtering study), beilstein-journals.org/bjnano/articles/13/86 [VERIFY exact table] |
| Ar -> Si | 0.5 keV | SRIM mean range | ~2.3 nm (displacement depth ~2 nm) | arXiv:2412.11470 (UHV annealing of sputtered Si(111)), SRIM-cited |
| Ar -> Si (XPS profiling context) | 0.5 keV | projected range | 1.6 ± 0.5 nm | "Evaluation methods for XPS depth profiling: A review" (2025), ScienceDirect S2666523925001825 |
| Ar -> Si | 3 keV | amorphous layer | ~2.1 nm (MD); RBS/channeling damage depths larger (several nm) | Beilstein 2022 (MD); classic: "Low-energy (2–5 keV) argon damage in silicon" (RBS/channeling) [VERIFY authors/values] |
| Ar -> SiO2 | 17 keV | SRIM damage profile | full profile published | "Engineering Silicon Oxide by Argon Ion Implantation...", *Front. Mater.* 9 (2022), DOI: 10.3389/fmats.2022.813407 |
| Ar -> SiO2 | 160 keV | concentration peak depth | ~190 nm (peak 1e19 cm^-3 at 2e14 cm^-2); nuclear deposited energy ~flat from surface to 160 nm then falls | R. Charavel, J.-P. Raskin, *Electrochem. Solid-State Lett.* 9(7), G245–G247 (2006), DOI: 10.1149/1.2200307 |
| Ar -> a-C (DLC) | 4–16 keV | SRIM energy/angular distributions (sputter geometry, 45°) | published as heatmaps, no Rp table | arXiv:2406.07144 (ECR-IBD DLC, SRIM-2013) |

Gaps: no published sub-keV Ar Rp for SiO2 specifically, and nothing for
a-C:F/PTFE-like targets at these energies (fluoropolymer SRIM runs exist only
at >>10 keV in the ion-track literature).

**Action for petch**: generate our own SRIM-2013 tables for Ar -> SiO2
(rho = 2.20 g/cm3) and Ar -> a-C:F (use PTFE, rho = 2.0–2.2 g/cm3, and a
CF1.2 "plasma polymer" at rho ≈ 1.6–1.9 g/cm3) at 0.25/0.5/1/2/3/5 keV, commit
the raw outputs as the validation fixture, and check them against (a) the spot
anchors above and (b) our analytic ZBL+LS integrator (expect ≤10–20% and the
known SRIM shallow-depth positive bias of 2–6 nm at the top end). Rule-of-thumb
consistency check: damage depth stays 1–3 nm below ~1 keV and ~1 nm per keV
slope thereafter for Ar in Si-based targets (consistent with all anchors above
and with Standaert's "ion penetration ~1 nm" at 100s of eV, §2).

---

## 2. Reactive mixed-layer thickness in fluorocarbon etching of SiO2 / Si

| Source | Conditions | Layer | Thickness | Energy dependence |
|---|---|---|---|---|
| N. R. Rueger, J. J. Beulens, M. Schaepkens, M. F. Doemling, J. M. Mirza, T. E. F. M. Standaert, G. S. Oehrlein, *J. Vac. Sci. Technol. A* 15, 1881–1889 (1997) — "Role of steady state fluorocarbon films in the etching of silicon dioxide using CHF3 in an inductively coupled plasma reactor" | CHF3 ICP, SiO2 | steady-state CFx film | order 1 nm | "slight variations in the film thickness, on the order of 1 nm, can result in ... ~400 nm/min" etch-rate change; thickness falls with self-bias |
| T. E. F. M. Standaert, M. Schaepkens, N. R. Rueger, P. G. M. Sebel, G. S. Oehrlein, J. M. Cook, *J. Vac. Sci. Technol. A* 16, 239–249 (1998), DOI: 10.1116/1.580978 — "High density fluorocarbon etching of silicon in an inductively coupled plasma: Mechanism of etching through a thick steady state fluorocarbon layer" | CHF3/C2F6/C3F6(/H2) ICP, 1400 W, self-bias to −150 V | CFx on **Si**: 2–7 nm; CFx on **SiO2**: <1.5 nm (typically <1 nm); SiFy reacted layer under CFx measured by XPS+ellipsometry (Fig. 4), ~1–2 nm scale; F/C gradient over top 2–3 nm | CFx thins with increasing self-bias (ion energy); paper states "the penetration depth of ions in the energy range investigated is about 1 nm" | direct quotes verified from full text |
| M. Schaepkens, T. E. F. M. Standaert, N. R. Rueger, P. G. M. Sebel, G. S. Oehrlein, J. M. Cook, *J. Vac. Sci. Technol. A* 17, 26–37 (1999) — "Study of the SiO2-to-Si3N4 etch selectivity mechanism..." | CHF3, C2F6/C3F6, C3F6/H2 ICP | FC film on all substrates | substrate etch rate **inversely proportional to FC film thickness**; SiO2 thinnest, Si thickest | selectivity = differential FC thickness |
| T. E. F. M. Standaert, C. Hedlund, E. A. Joseph, G. S. Oehrlein, T. J. Dalton, *J. Vac. Sci. Technol. A* 22, 53–60 (2004), DOI: 10.1116/1.1626642 | Si, SiO2, Si3N4, SiCH; multiple FC gases | steady-state FC film | thickness NOT the sole rate-controlling parameter; **ion-induced defluorination** of the FC film plays a major role | key caveat for a pure thickness-based model |
| M. E. Barone, D. B. Graves, *J. Appl. Phys.* 77, 1263–1274 (1995) — "Chemical and physical sputtering of fluorinated silicon" (MD) | Ar+ 20/50/200 eV on SiFx layers | fluorinated Si mixed layer | ~1–2 nm at these energies; weakly bound SiFx (x=1–3) formed in layer | total yield ∝ sqrt(E); weak-bond formation threshold extrapolates to ≤4 eV |
| D. Humbird, D. B. Graves, *J. Appl. Phys.* 96, 65 (2004) — "Fluorocarbon plasma etching of silicon: Factors controlling etch rate" (MD) | FC radicals + Ar+ on Si | FC/SiFx mixed layer | few nm, grows with ion energy | MD counterpart of Standaert 1998 picture |
| R. L. Bruce, F. Weilnboeck, ... G. S. Oehrlein, D. B. Graves, *J. Appl. Phys.* 107, 084310 (2010) | 100 eV Ar+ on polystyrene | densified, crosslinked ion-damaged layer | "a few nanometers" | layer thickness tracks ion damage range; VUV modifies ~100 nm |

Consensus for the model: reacted/mixed layer (SiOxFyCz on SiO2, SiFx/CFx on
Si) is **1–3 nm at 50–500 eV**, scaling with ion range (≈ sqrt(E) in this
regime), with the FC overlayer on Si up to 2–7 nm at low bias. Any
deposition-in-layer yield should use d ≈ 1–2 nm for SiO2 and let d grow with
ion range.

---

## 3. Radiation-chemistry G-values (events per 100 eV absorbed)

### 3.1 Measured values

| Material | Condition | G(scission) | G(crosslink) | Other | Source |
|---|---|---|---|---|---|
| PTFE | gamma, vacuum, 20 °C | ≈ 2 | ~0 (scission-dominant) | G(radicals) 0.15–0.19; G(total gas) ≈ 0.3 (SiF4* 0.12–0.16, CO2 0.06–0.12, CF4 0.004–0.009) | R. E. Florin, L. A. Wall, *J. Res. NBS A* 65A, 375–387 (1961), PMC5287141. (*SiF4 from glass-ampoule F attack) |
| PTFE | gamma, air, RT | ≈ 10 | ~0 | oxidative amplification | Florin & Wall 1961 |
| PTFE | EB/gamma, **molten 613 K**, O2-free | drastically reduced | **0.35 (lower limit)** | G(S) rises with T up to 600 K | A. Oshima, S. Ikeda, H. Kudoh, T. Seguchi, Y. Tabata, *Radiat. Phys. Chem.* 50, 611–615 (1997), DOI: 10.1016/S0969-806X(97)00103-5 |
| PTFE (heat-treated, 2 MGy molten-state) | 19F MAS NMR | — | up to **1.85** | branching Y-units observed | *Radiat. Phys. Chem.* (1999), DOI: 10.1016/S0969-806X(98)00250-3 [VERIFY authors] |
| PCTFE | gamma, vacuum or air | 0.67 (constant) | none detected | G(R) ≈ 1.0 | Florin & Wall 1961 |
| TFE-PMVE copolymer | radiolysis, NMR quantified | **1.4** | **0.9** | crosslink/scission ratio 0.64 | "NMR Study of the Radiation-Induced Cross-Linking of Poly(TFE-co-PFMVE)", *Macromolecules* (1997), DOI: 10.1021/ma970656l |
| TFE-PMVE (PFA-like) EB | volatile G-values | — | — | G(CF4)=0.93, G(COF2)=0.31, G(CO2)=0.055, G(CF3OCF3)=0.14 | *Macromolecules*, DOI: 10.1021/ma960668r |
| H-containing fluoropolymers (VDF-HFP etc.) | gamma | secondary | **dominant initially** | G(H2) 0.11–0.27 | Florin & Wall 1961 |

Spread to encode in the model: fluoropolymer G(S) ≈ 0.7–2 (vacuum, RT, CF2-rich
backbone), G(X) ≈ 0–0.9 rising to ≈ 0.35–1.85 when the matrix is mobile and
F-depleted; any C–H content flips the balance toward crosslinking (classic
Charlesby rule: PTFE degrades, PE crosslinks).

### 3.2 Plasma-deposited a-C:F films and ion irradiation

- **No direct G(X)/G(S) measurements for plasma-polymerized fluorocarbon films
  were found** — this is a genuine literature gap. The proxy chain is:
  (a) plasma polymers are already crosslinked and F-deficient (F/C ≈ 1–1.2 vs
  2 for PTFE), which by the H/branching rule shifts them crosslink-dominant;
  (b) under 100 eV Ar+ a densified, dehydrogenated/defluorinated crosslinked
  skin of a few nm forms and **greatly reduces sputter yield** (Bruce et al.,
  JAP 107, 084310 (2010) — polystyrene; MD by Graves group), i.e. the etch-side
  evidence of ion-induced crosslinking during processing;
  (c) Standaert 2004 (10.1116/1.1626642): ion-induced **defluorination** of the
  steady-state FC film controls etching — same physics, XPS-observed.
- Ion vs gamma/electron: crosslinking G rises with LET, and Charlesby–Pinner
  underestimates G(X) for ions because intra-track crosslinks don't contribute
  to gelation. Track-overlap / dose-competition framework: L. Calcagno,
  G. Compagnini, G. Foti, "Structural modification of polymer films by ion
  irradiation", *Nucl. Instrum. Methods B* 65, 413–422 (1992) [VERIFY pages];
  see also L. Calcagno, G. Foti, "Ion irradiation of polymers" (100 keV He /
  200 keV Ne / 400 keV Ar on polystyrene, gel-fraction analysis).
- Closest published "model" of ion-induced crosslinking during etch-like
  exposure: the Bruce/Oehrlein/Graves 2010 picture (crosslinked layer thickness
  = ion damage range; competition between crosslinking dose and sputter
  removal) — a DPA-style dose-competition model, not G-value-parameterized.
  petch formalizing this with G-values + ZBL deposited dose would be new.

---

## 4. Prior art: etch yield as energy-deposited-in-surface-layer

| Who | What | Status |
|---|---|---|
| P. Sigmund, *Phys. Rev.* 184, 383 (1969), DOI: 10.1103/PhysRev.184.383 | Y = Λ · F_D(0): yield ∝ nuclear energy deposited at the surface; linear cascade theory; Λ ≈ 0.042/(N·U0) Å/eV | the root formula; valid keV+, breaks near threshold |
| C. Steinbrüchel, *Appl. Phys. Lett.* 55, 1960 (1989), DOI: 10.1063/1.102336 | showed Sn(E) ∝ sqrt(E) at low E ⇒ universal Y = A(√E − √Eth) for **both** physical sputtering and ion-enhanced chemical etching (Si, SiO2, metals; noble + reactive ions) | the standard low-energy reduction; used by essentially every plasma-etch surface model since (incl. Kushner-group HPEM/MCFPM rate coefficients) |
| M. E. Barone, D. B. Graves, *J. Appl. Phys.* 77, 1263 (1995) | MD on fluorinated Si: total (physical+chemical) yield follows sqrt(E); chemical channel = weakly bound SiFx created **in the mixed layer** by deposited energy; thresholds ≤4 eV (chem) vs ~20 eV (phys) | microscopic justification that "energy deposited in the reacted layer" drives chemical yield |
| D. C. Gray, I. Tepermeister, H. H. Sawin, *J. Vac. Sci. Technol. B* 11, 1243 (1993) | phenomenological ion-neutral synergy model of ion-enhanced F/Si etching (beam data) | the kinetics wrapper the deposition-in-layer model must reproduce |
| J. P. Chang, H. H. Sawin (beam studies, JVST A ~15, 610 (1997)) [VERIFY citation] | measured Y(E) = A(√E − √Eth) for Cl/poly-Si ion-enhanced etching | quantitative sqrt-law confirmation |
| Y. Yamamura, H. Tawara, *At. Data Nucl. Data Tables* 62, 149 (1996) | empirical low-energy sputter-yield formula (threshold-corrected Sigmund) | validation dataset for the physical-sputter limit |
| Wittmaack / ICRU 49 (1993) | ZBL nuclear stopping adopted for standard-reference stopping | pedigree for using ZBL Sn as F_D input |

What worked historically: Sigmund's Λ·F_D(0) with ZBL Sn is quantitative above
a few hundred eV; below that, the working recipe collapsed to Y = A(√E − √Eth)
(Steinbrüchel) with A **fitted** per surface chemistry — A is exactly the knob
petch is deriving: A ∝ (energy fraction deposited within the reacted layer of
thickness d, per §2) / (effective binding of the fluorinated/crosslinked
matrix, modifiable by the G-value chemistry of §3). No published model derives
A from ZBL deposition-in-layer + radiation-chemistry crosslinking — the
combination appears to be open territory.

---

## Sources (primary)

1. Ziegler, Biersack, Littmark, *The Stopping and Range of Ions in Solids*, Pergamon (1985); Ziegler et al., NIMB 268, 1818 (2010), DOI: 10.1016/j.nimb.2010.02.091.
2. Lindhard & Scharff, Phys. Rev. 124, 128 (1961), DOI: 10.1103/PhysRev.124.128; LSS, Mat. Fys. Medd. 33 no. 14 (1963).
3. ZBL constants cross-checks: arXiv:1203.4620; sr-niel.org long write-up; Univ. Helsinki radiation-damage course notes (mv.helsinki.fi/home/knordlun/rad_dam_course/str_skador4.pdf); Wikipedia "Stopping power (particle radiation)" (screening function, a_U, 18% fit spread).
4. SRIM shallow-implant bias: OSTI 1408290; OSTI 1830531.
5. Rueger et al., JVST A 15, 1881 (1997). 6. Standaert et al., JVST A 16, 239 (1998), DOI: 10.1116/1.580978 (full text verified). 7. Schaepkens et al., JVST A 17, 26 (1999). 8. Standaert et al., JVST A 22, 53 (2004), DOI: 10.1116/1.1626642.
9. Barone & Graves, JAP 77, 1263 (1995). 10. Humbird & Graves, JAP 96, 65 (2004). 11. Bruce et al., JAP 107, 084310 (2010).
12. Florin & Wall, J. Res. NBS 65A, 375 (1961) (PMC5287141, full text verified). 13. Oshima et al., Radiat. Phys. Chem. 50, 611 (1997), DOI: 10.1016/S0969-806X(97)00103-5. 14. Macromolecules DOI: 10.1021/ma970656l (TFE-PMVE G(S)=1.4, G(X)=0.9); DOI: 10.1021/ma960668r. 15. Calcagno/Compagnini/Foti, NIMB 65, 413 (1992).
16. Sigmund, Phys. Rev. 184, 383 (1969). 17. Steinbrüchel, APL 55, 1960 (1989), DOI: 10.1063/1.102336. 18. Yamamura & Tawara, ADNDT 62, 149 (1996).
19. Charavel & Raskin, Electrochem. Solid-State Lett. 9, G245 (2006), DOI: 10.1149/1.2200307 (full text verified — nuclear-deposited-energy vs etch-rate correlation in SiO2).
20. Kaler, Lou, Donnelly, Economou, J. Phys. D 50, 234001 (2017), DOI: 10.1088/1361-6463/aa6f40 (FC-ALE of SiO2 context).
