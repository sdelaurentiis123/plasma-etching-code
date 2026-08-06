# Atomic-fluorine supply at Krüger base-case conditions — sourced band, budget impact, and an information request

**Status: research note. NOT committed. No constant landed. No gate moved.**
Written 2026-08-06 against `VALIDATION_DOSSIER_KRUEGER_2026-08-05.md`,
`BENCHMARK_CERTIFICATION_2026-08-06.md`, `RESULTS_E8_COUPLED_2026-08-05.md`,
`RESULTS_ML23_GRAY_LAWS_2026-08-05.md`.

**Evidence tags used throughout.**

| tag | meaning |
|---|---|
| `[Q]` | verbatim quote read directly from the primary document (local PDF/text extract, or a full text I downloaded and converted myself) |
| `[Q-relay]` | value reached me through an abstract page, a search-tool summary, or a secondary citation — **not** read in the primary body. Treated as weaker than `[Q]`. |
| `[INF]` | my inference/arithmetic, with the chain shown |
| `[VERIFY]` | not confirmed; must not be used to move a gate |

---

## 0. The hole this note fills

`chemistry_deck.py` (Krüger 2024 deck) carries, verbatim in the source:

```python
"fluorine": [],          # the boundary publishes no atomic-F flux
```

and `mixed_layer_mechanism.py`:

```python
# Krueger 2024 species channels: fluorocarbon radicals with explicit C/F
# stoichiometry (the boundary publishes NO atomic-F flux — fluorine reaches
# the layer through the film AND by direct complex-formation chemisorption).
```

That is accurate transcription. Krüger's Table I is `[Q]`
(`research_sources/thesis_extracts/krueger-2024.txt` L289-297, and thesis
L4342-4349):

> TABLE I. Base case fluxes to wafer.
> C3F4 9.5 × 10^16 / C2F3 6.8 × 10^16 / CF 4.4 × 10^16 / CF2 9.4 × 10^16 /
> CF3 8.4 × 10^15 / O 7.7 × 10^16 / Ions 1.2 × 10^16   (all cm⁻² s⁻¹)

and the sentence that *selects* that list is `[Q]` (same file, L273-277):

> "These fluxes include the fluorocarbon radicals most **responsible for
> polymer deposition** as well as atomic oxygen, which etch and remove the
> polymer."

**So the table is a deposition-balance selection, not a species inventory.**
Absence of F from Table I is not a statement that F is absent from the HPEM.

**The mechanism he published requires an F flux the boundary does not carry.**
Counted directly in `research_sources/thesis_extracts/krueger_thesis.txt` with
`grep -cE "^\s*\S+\(s\)\s+\+\s+F\s+→"` → **13 rows** with gas-phase atomic F as
a reactant. `[Q]`, exact lines:

| line | row | probability |
|---|---|---|
| 5895 | `AC(s) + F → CF(s)` | 0.01 |
| 5896 | `AC(xs) + F → CF(s)` | 0.001 |
| 5897 | `CF(s) + F → CF2(s)` | 0.01 |
| 5898 | `CF2(s) + F → CF3(s)` | 0.01 |
| 5899 | `CF3(s) + F → CF4` | 0.01 |
| 6021 | `SiO2*(s) + F → SiO2(s) + F` | **1** (deactivation — a *loss*) |
| 6548-6555 | `SiO2CF(s)+F→SiO2CF2(s)`, `SiO2CF2(s)+F→SiO2CF3(s)`, `SiO2C2F3(s)+F→SiO2C2F4(s)`, `SiO2C3F5(s)+F→SiO2C3F6(s)` and their `*` (activated) twins | **0.1** each |

Krüger's mechanism inherits from Huang, whose Appendix E species table is `[Q]`
(`huang_thesis.txt` L10140-10143) and lists `F+, F(h), F` and `F2+, F2(h), F2`
as gas-phase species, with a "Fluorination of passivated surface" block at
p₀ = 0.1 (L10214-10222).

Consequence: **the published boundary and the published mechanism are
inconsistent with each other as printed.** Running the mechanism with Γ_F = 0
(what petch does today, faithfully) silently zeroes 13 rows including all four
complex-fluorination channels and the SiO₂\* deactivation channel.

---

## 1. Model-published (HPEM) F fluxes — the strongest anchors

These are the same code family (HPEM/MCFPM, Kushner group) that produced
Krüger's boundary. All three were read directly from the primary document.

### 1.1 ANCHOR A — Huang, Ar/C₄F₈/O₂ tri-frequency CCP, 25 mTorr

**Conditions (verbatim, `[Q]`, `huang2019.txt` L749 / thesis L6122):**
> "Operating conditions: Ar/C4F8/O2 = 75/15/10, 25 mTorr, 500 sccm,
> 80/10/5 MHz power = 0.4/2.5/5 kW."
Fluxes sampled at the wafer; parametric IEAD study at r = 7.5 cm — the same
radial position Krüger reports.

**Value (verbatim, `[Q]`):**
> "CF2 and CF3 are the dominant reactive fluorocarbon species incident onto the
> wafer with average fluxes of 1.4 × 10¹⁷ and 0.9 × 10¹⁷ cm⁻² s⁻¹, respectively.
> O and F atoms are mainly produced through electron impact dissociation of O2
> and CFx species and diffuse to the wafer with average fluxes of 1.2 × 10¹⁷ and
> **0.8 × 10¹⁷ cm⁻² s⁻¹**, respectively."

**Independently cross-checked in two documents.**
- S. Huang PhD thesis, U. Michigan, `huang_thesis.txt` **L5272-5275**.
- S. Huang, C. Huard, S. Shim, S. K. Nam, I.-C. Song, S. Lu, M. J. Kushner,
  "Plasma etching of high aspect ratio features in SiO₂ using Ar/C₄F₈/O₂
  mixtures: A computational investigation", *J. Vac. Sci. Technol. A* **37**,
  031304 (2019), DOI [10.1116/1.5090606](https://doi.org/10.1116/1.5090606) —
  downloaded from the group's open archive
  (`https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf`),
  **L726-729**. Numbers identical to the thesis to the last digit. **AGREE.**

**Total fluorocarbon radical flux at the same plane (verbatim, `[Q]`,
huang2019 L~950 / thesis L5451):**
> "The flux of CFx and CxFy to the top of the PR is 3.0 × 10¹⁷ cm⁻² s⁻¹"

**Ions at the same plane (verbatim, `[Q]`):** Ar⁺ 3.9 × 10¹⁵, C₂F₄⁺
2.0 × 10¹⁵ cm⁻² s⁻¹; "Fluxes of ions to the wafer surface are lower than those
of radicals by 1–2 orders of magnitudes".

**Derived ratios `[INF]`** (arithmetic shown, all in cm⁻² s⁻¹):
- F/O = 0.8e17 / 1.2e17 = **0.6667**
- F/CF₂ = 0.8e17 / 1.4e17 = **0.5714**
- F/ΣFC(reactive) = 0.8e17 / 3.0e17 = **0.2667**

**Internal consistency check `[INF]`:** CF₂ + CF₃ = 1.4e17 + 0.9e17 = 2.3e17,
which is 77 % of the stated 3.0e17 total, leaving 0.7e17 for CF/C₂F₃/etc.
Self-consistent. Note C₂F₄ is explicitly *excluded* from the "reactive"
count — `[Q]` "Although C₂F₄ is not particularly reactive in oxide etching".

### 1.2 ANCHOR B — Wang, Ar/C₄F₈/O₂ dc-augmented CCP, 40 mTorr

**Conditions (verbatim, `[Q]`, `wang_mingmei_phd_thesis.txt` L2543-2546):**
> "The base case operating conditions are Ar/C4F8/O2 = 80/15/5 at 40 mTorr with
> a flow rate of 300 sccm (standard cubic centimeter per minute at STP). The
> substrate is biased at 10 MHz delivering a power of 4 kW and the dc electrode
> delivers 200 W."
Reactor source cited by that chapter as ref. [3] = M. Wang and M. J. Kushner,
*J. Appl. Phys.* **107**, 023309 (2010).

**Value (`[Q]`, thesis "Table I. Fluxes of radicals and ions incident onto the
wafer", L2849-2887, units 10¹⁵ cm⁻² s⁻¹):**

| ions | | neutrals | |
|---|---|---|---|
| Ar⁺ | 4.23 | CF | 0.32 |
| CF₃⁺ | 0.26 | CF₂ | 1.86 |
| CF₂⁺ | 0.15 | C₂F₄ | 29.6 |
| F⁺ | 0.08 | C₃F₅ | 0.32 |
| O⁺ | 0.04 | **F** | **3.2** |
| C₂F₄⁺ | 4.28 | F₂ | 0.04 |
| C₃F₅⁺ | 1.23 | | |
| C₃F₆⁺ | 0.002 | | |
| C₄F₇⁺ | 0.006 | | |

**Derived `[INF]`:** Σions = 10.278e15; ΣCₓF_y(neutral) = 32.10e15.
F/CF₂ = 3.2/1.86 = **1.720**; F/ΣCₓF_y = 3.2/32.10 = **0.0997**;
F/Σions = 3.2/10.278 = **0.3113**; ΣCₓF_y/Σions = **3.12**.

**Not cross-checkable against a second document** (the thesis table is the only
copy I read). Recorded here because it is the only *published atomic-F wafer
flux for a CCP with an O₂-bearing C₄F₈ feed*, but it is **excluded from the
band derivation** for a stated reason (§3.3): its neutral/ion ratio of 3.12 is
8× below Krüger's 25.8, i.e. it is a far less dissociated plasma and its ratios
do not transfer.

### 1.3 ANCHOR C — Sankaran & Kushner, pure c-C₄F₈ ICP, 20 mTorr (zero-O₂ floor)

**Citation:** A. Sankaran and M. J. Kushner, "Etching of porous and solid SiO₂
in Ar/c-C₄F₈, O₂/c-C₄F₈ and Ar/O₂/c-C₄F₈ plasmas", *J. Appl. Phys.* **97**,
023307 (2005), DOI [10.1063/1.1834979](https://doi.org/10.1063/1.1834979).
Downloaded from `https://cpseg.eecs.umich.edu/pub/articles/jap_97_023307_2005.pdf`
and converted locally; read directly. `[Q]`

**Conditions (verbatim table caption):**
> "TABLE II. Fluxes of selected species to the center of the wafer in a C4F8
> plasma (600-W ICP, 20 mTorr, 40 sccm, −120-V dc bias)."

**Values (`[Q]`, cm⁻² s⁻¹):** CF₃⁺ 3.91e14, CF₂⁺ 1.35e15, CF⁺ 6.17e14,
C₂F₄⁺ 3.75e15, CF₂ 4.08e17, CF 5.59e16, **F 6.73e15**, C₂F₃ 1.46e16,
C₂F₄ 3.55e17.

**Derived `[INF]`:** ΣCₓF_y = 4.08e17 + 5.59e16 + 1.46e16 + 3.55e17 = 8.335e17.
F/CF₂ = **0.01650**; F/ΣCₓF_y = **0.00807**.

**Role:** this is the **zero-oxygen limit** of the same code. It bounds the
band from below *only for a no-O₂ feed*; Krüger's feed is 30.4 % O₂.

### 1.4 The O₂ ordering inside the HPEM family `[INF]`

Normalising the feed by **O₂ molecules per feed carbon atom** (the quantity
that controls carbon scavenging, hence F liberation):

| anchor | feed | O₂/C in feed | F/ΣFC at wafer | p (mTorr) |
|---|---|---|---|---|
| C — Sankaran 2005 | c-C₄F₈ (no O₂) | 0.000 | 0.0081 | 20 |
| B — Wang thesis | Ar/C₄F₈/O₂ 80/15/5 | 0.0833 | 0.0997 | 40 |
| A — Huang 2019/thesis | Ar/C₄F₈/O₂ 75/15/10 | 0.1667 | 0.2667 | 25 |
| **Krüger base case** | **C₄F₆/Ar/O₂ 140/100/105** | **0.1875** | **?** | **10** |

Arithmetic for Krüger: 105 / (140 × 4 C per C₄F₆) = 0.1875. For Huang:
10 / (15 × 4) = 0.1667. For Wang: 5 / (15 × 4) = 0.0833.

**Krüger's feed sits 12 % above Huang's on O₂/C, and 25 % below on F/C**
(C₄F₆ carries 6 F per 4 C = 1.5; C₄F₈ carries 8/4 = 2.0; ratio 0.75).
That places Krüger *immediately adjacent to Anchor A on both axes* — which is
why Anchor A carries the band and B/C only bound it.

### 1.5 Krüger's own published O₂ sweep — what it does and does not say

`[Q]`, `krueger-2024.txt` L946-951 (peer-reviewed text, cross-checked against
thesis L4816-4820, identical):

> "The variation in the O2/C4F6 ratio was achieved, experimentally and in the
> simulation, by adjusting the rate of oxygen inflow only. As a consequence,
> the fluxes of the major contributing fluorocarbons remain nearly constant.
> The O-atom flux increases from 4.1 × 10¹⁶ to 1.5 × 10¹⁷ cm⁻² s⁻¹, a factor of
> 3.6, with an increase of the oxygen inflow by factor of 5."

This is his own model saying the FC fluxes are ~invariant across the transfer
sweep while O rises 3.6×. **He does not report what F did.** Figure 16(a) of the
paper (digitised in `scripts/digitize_krueger_2024_transfer_boundary.py`) has
exactly **seven** series — CF₂, C₃F₄, C₂F₃, CF, O, ions, CF₃ — and no F trace.

**Also relevant `[Q]`** (`krueger-2024.txt` L266-270, thesis L4324-4329): the
HPEM resolves F⁻ with "a maximum density of 2.2 × 10¹⁰ cm⁻³", the largest
negative-ion density in the discharge. F⁻ is produced overwhelmingly by
dissociative attachment to F-bearing neutrals, so the model demonstrably tracks
a substantial gas-phase fluorine budget it does not publish at the wafer. `[INF]`

---

## 2. Measured F densities and fluxes in fluorocarbon plasmas

Every row states the technique, the reactor, and whether the number is measured
or model-derived. **No number in this section is used to set the band centre**;
they are used to bound it and to sanity-check the magnitude.

### 2.1 Unit conversion, shown once, explicitly

Thermal (one-sided) flux from an isotropic density:

    Γ = n · v̄ / 4 ,   v̄ = sqrt(8 k_B T / (π m))

with k_B = 1.380649e-23 J/K, m(F) = 18.998403 u × 1.66053907e-27 kg/u
= 3.15477e-26 kg. Evaluated:

| T_gas (K) | v̄(F) (m/s) | v̄/4 (m/s) | Γ for n = 1e12 cm⁻³ (= 1e18 m⁻³) |
|---|---|---|---|
| 300 | 578.2 | 144.6 | 1.446e20 m⁻² s⁻¹ |
| 400 | 667.7 | 166.9 | 1.669e20 |
| **500** | **746.5** | **186.6** | **1.866e20** |
| 600 | 817.7 | 204.4 | 2.044e20 |
| 800 | 944.2 | 236.1 | 2.361e20 |

Also note: 1 cm⁻² s⁻¹ = 1e4 m⁻² s⁻¹, and 1 cm⁻³ = 1e6 m⁻³. Krüger's Table I is
in cm⁻² s⁻¹; petch's deck is in m⁻² s⁻¹; the factor is **1e4**, e.g. CF₂
9.4e16 cm⁻² s⁻¹ = 9.4e20 m⁻² s⁻¹ ✓ (matches `chemistry_deck` / `floor_delivery_scan.py`).

**Conversion validated against the source's own paired numbers `[INF]`:** Huang
gives a volume-averaged C₂F₄ density of 3.2 × 10¹³ cm⁻³ `[Q]` and a total FC
radical flux of 3.0 × 10¹⁷ cm⁻² s⁻¹ `[Q]`. Applying the formula with
m(C₂F₄) = 100.02 u at 500 K: v̄ = 325.3 m/s, Γ = 3.2e19 m⁻³ × 81.3 m/s
= 2.60e21 m⁻² s⁻¹ = 2.60e17 cm⁻² s⁻¹ — the same order as his stated *total*,
consistent with C₂F₄ being the dominant dissociation product and with the
"3.0e17" figure counting only the reactive subset. The formula is therefore not
misapplied by an order of magnitude. (At 300 K the same arithmetic gives
2.02e17 cm⁻² s⁻¹.)

### 2.2 Measured absolute F densities — RIE / CCP class

| source | technique | reactor + conditions | measured [F] | tag |
|---|---|---|---|---|
| J.-S. Jenq, J. Ding, J. W. Taylor, N. Hershkowitz, "Absolute fluorine atom concentrations in RIE and ECR CF₄ plasmas measured by actinometry", *Plasma Sources Sci. Technol.* **3**, 154 (1994), DOI [10.1088/0963-0252/3/2/005](https://doi.org/10.1088/0963-0252/3/2/005) | Ar actinometry | **RIE**, CF₄, 20 sccm, **13–82 mTorr**, 250–1500 W | **(0.8 – 4.2) × 10¹³ cm⁻³** | `[Q-relay]` (abstract page only) |
| same | same | **ECR**, 0.5–3.5 mTorr, 500–900 W µwave, −50 V bias, at wafer stage | (0.4 – 4) × 10¹² cm⁻³ | `[Q-relay]` |

Reported reproducibility ±10 % with a stated potential factor-of-two systematic
from cross-section uncertainty. `[Q-relay]`

**Flux equivalent `[INF]`:** the RIE band 0.8–4.2e13 cm⁻³ at 500 K →
Γ = (0.8–4.2)e19 m⁻³ × 186.6 = **1.5e21 – 7.8e21 m⁻² s⁻¹**. This is the
*most F-rich* feed (pure CF₄, F/C = 4) and therefore an upper envelope, not a
target.

### 2.3 Measured/model F densities in C₄F₈-based mixtures, and the O₂ direction

**Source:** I. Chun, A. Efremov, G. Y. Yeom, K.-H. Kwon, "A comparative study
of CF₄/O₂/Ar and C₄F₈/O₂/Ar plasmas for dry etching applications", *Thin Solid
Films* **579**, 136–143 (2015), DOI
[10.1016/j.tsf.2015.02.060](https://doi.org/10.1016/j.tsf.2015.02.060).
PDF downloaded and converted locally; read directly. `[Q]`

**Conditions (verbatim):** planar ICP, "gas flow rate (q = 40 sccm), gas
pressure (p = 6 mTorr or 0.8 Pa)"; 50 % fluorocarbon held constant while Ar is
substituted by O₂ from 0 to 50 %. **Method is a 0-D global model constrained by
double-Langmuir-probe measurements of T_e and J₊ — the F densities are
model outputs, not direct F measurements.** Stated verbatim by the authors:
> "the changes in O and F atom densities are quite close to those obtained
> experimentally by Takashi et al. [22]."

**Values (verbatim `[Q]`):**
- CF₄/Ar/O₂: "an increase in F atom density (nF = 5.8 × 10¹² – 2.6 × 10¹³ cm⁻³
  for 0–12 % O₂)" — a **4.5× rise**.
- C₄F₈/Ar/O₂: "the F atom density decreases monotonically toward O₂-rich plasmas
  (nF = 8.0 × 10¹² – 8.6 × 10¹¹ cm⁻³ for 0–50 % O₂)" — a **9.3× fall**.
- Mechanism given verbatim: "the low densities of the O and O(1D) atoms limit
  the formation rates for FO and CFO species… even in the presence of oxygen,
  the F atom formation rate is composed only of R11–R14 and has no support from
  both atom-molecular reactions and electron-impact dissociation of FO (R21)
  and CFO (R22)."

**This is a direct sign-conflict with the HPEM ordering in §1.4 and must not be
papered over.** In this 6 mTorr ICP with 50 % C₄F₈, adding O₂ *lowers* F; in the
Kushner HPEM CCP family, higher O₂/C correlates with higher F/ΣFC. Both can be
true — the Chun case substitutes O₂ for **Ar** at fixed FC fraction in a
low-density discharge where O₂ is destroyed by CF/C before it can dissociate,
whereas the HPEM cases differ in reactor, pressure, power density and total
dissociation. **Operational consequence, stated as a limit on the claim:** the
band below is derived at Krüger's *base-case* O₂/C only. It must **not** be
extrapolated across his O₂/C₄F₆ = 0.5 → 2.5 transfer sweep; on that axis the
sign of dΓ_F/d(O₂) is **unresolved** `[VERIFY]`.

### 2.4 Actinometry reliability for F in C₄F₈ — a caution that constrains the whole measured literature

**Source:** Y. Kawai, K. Sasaki, K. Kadota, "Comparison of the Fluorine Atom
Density Measured by Actinometry and Vacuum Ultraviolet Absorption
Spectroscopy", *Jpn. J. Appl. Phys.* **36**, L1261 (1997), DOI
[10.1143/JJAP.36.L1261](https://doi.org/10.1143/JJAP.36.L1261). `[Q-relay]`
(abstract page).

Reported: high-density CF₄ and C₄F₈ helicon plasmas; relative F density from
actinometry accurate "within a factor of ~2 for the same working gas"; the
actinometric constant K "varied significantly between CF₄ and C₄F₈ plasmas,
making absolute density determination difficult except for rough estimates
using K ≈ 2."

**Consequence:** actinometric absolute F densities in *fluorocarbon* plasmas
carry a factor-≈2 systematic, and cross-gas comparisons (CF₄ vs C₄F₆/C₄F₈)
carry more. This is why §3 builds the band on same-code HPEM ratios and uses
the measured densities only as an envelope.

### 2.5 Measured F surface loss probability — needed for the transport step

**Source:** K. Sasaki, Y. Kawai, C. Suzuki, K. Kadota, "Kinetics of fluorine
atoms in high-density carbon–tetrafluoride plasmas", *J. Appl. Phys.* **82**,
5938–5943 (1997), DOI [10.1063/1.366495](https://doi.org/10.1063/1.366495).
Abstract verbatim via OpenAlex record for that DOI `[Q-relay]`:

> "Reaction processes of fluorine (F) atoms in high-density carbon–tetrafluoride
> (CF4) plasmas were investigated using vacuum ultraviolet absorption
> spectroscopy… The surface loss probability of F atoms on the chamber wall was
> evaluated from the decay time constant in the late afterglow, and was **on the
> order of 10⁻³**."

**Independent bracket from Krüger's own mechanism `[Q]`:** his AC-mask and
polymer fluorination rows all carry p = 0.01 (§0 table), and the complexed-oxide
rows carry 0.1. So the per-strike loss of F on a *polymer-coated* feature
sidewall is bracketed **10⁻³ (measured, chamber wall) to 10⁻² (Krüger's own
polymer rows)**, and on the *complexed oxide floor* is 0.1 (Krüger's own rows).
Two independent sources agree that F is a **low-sticking, deeply-penetrating**
radical — the property that makes it matter at high AR. `[INF]`

### 2.6 Beam-side measurement that anchors what F does to the yield

`BENCHMARK_CERTIFICATION_2026-08-06.md` already carries Gray's measured
constants: SiO₂ ion-enhanced yield `β_e = 0.053(√E − √4)` "the number of SiF₄
molecules removed from **fluorine saturated** surface regions per incoming ion"
and F adsorption `s0 = 0.02` (Gray MIT thesis 1993, Table 5-10 + p.246,
co-regressed). The certification also records: "Gray's own printed parameters
predict [a half-rise at F/Ar⁺ of] **96–127** across the plausible energy
assignment". Those are the numbers §4 grades the band against.

**Independent measured corroboration that F, not CFₓ, is the oxide etchant:**
O. Kwon and H. H. Sawin, *J. Vac. Sci. Technol. A* **24**, 1906–1913 (2006),
DOI [10.1116/1.2336225](https://doi.org/10.1116/1.2336225) — neutral and ion
composition by QMS, yield by QCM, in C₄F₈ and C₄F₈+80 % Ar; already logged in
`RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` as "**atomic F flux identified as
the dominant etch driver, C as the deposition driver**". `[Q-relay]` (atlas
entry + abstract; the paper body was not read for this note).

---

## 3. The sourced band for Γ_F at Krüger base-case conditions

### 3.1 Inference chain, stated plainly

1. Krüger's boundary was produced by the HPEM. **`[Q]`**
2. Huang's Ar/C₄F₈/O₂ tri-frequency CCP at 25 mTorr is the closest published
   HPEM case to Krüger's reactor: same code, same group, same class (high-bias
   multi-frequency CCP for HARC SiO₂), and the *nearest neighbour on the
   O₂-per-feed-carbon axis* (0.1667 vs Krüger's 0.1875). **`[INF]`, §1.4**
3. Absolute fluxes do not transfer between reactors (different pressure, power,
   dissociation). **Ratios to species Krüger himself publishes** do transfer far
   better, because they divide out the overall dissociation level. Krüger
   publishes O, CF₂, and the FC sum — three independent normalisers. **`[INF]`**
4. Apply each of Huang's three ratios to Krüger's corresponding published flux.
5. Apply one chemistry correction: C₄F₆ carries 0.75× the fluorine per carbon
   that C₄F₈ does. **`[INF]`**
6. Bound the result against (a) the zero-O₂ HPEM point, (b) the measured density
   envelope, (c) the unresolved O₂ sign question.

### 3.2 The arithmetic, every step shown

All in cm⁻² s⁻¹ unless marked. Krüger Table I values as printed.

```
Krüger ΣFC(all 5 FC rows) = 9.5e16 + 6.8e16 + 4.4e16 + 9.4e16 + 8.4e15 = 3.0940e17
Krüger ΣFC(reactive; C3F4 declared inert in petch) = 6.8e16+4.4e16+9.4e16+8.4e15
                                                   = 2.1440e17

Estimate 1, via O:            (0.8e17 / 1.2e17) × 7.7e16 = 0.66667 × 7.7e16 = 5.133e16
Estimate 2, via CF2:          (0.8e17 / 1.4e17) × 9.4e16 = 0.57143 × 9.4e16 = 5.371e16
Estimate 3, via reactive-FC:  (0.8e17 / 3.0e17) × 2.1440e17 = 0.26667 × 2.1440e17 = 5.717e16
```

**Three independent normalisers agree to ±5.4 % about 5.41e16 cm⁻² s⁻¹.**
That mutual agreement is the single strongest piece of evidence in this note:
had the ratio method been unsound, O-, CF₂- and ΣFC-normalisation would have
scattered.

Convert to petch units (×1e4):

```
Estimate 1 = 5.133e20 m^-2 s^-1
Estimate 2 = 5.371e20
Estimate 3 = 5.717e20
mean       = 5.41e20 m^-2 s^-1     (uncorrected, i.e. "as if the feed were C4F8")
```

Apply the C₄F₆ fluorine-per-carbon correction (×0.75):

```
Estimate 1 → 3.850e20 ; Estimate 2 → 4.029e20 ; Estimate 3 → 4.288e20
mean       → 4.06e20 m^-2 s^-1
```

Krüger's O₂/C is 1.12× Huang's, which pushes the other way by an unquantified
amount. The two corrections partially cancel; I therefore take the **band
centre as 5.4 × 10²⁰ m⁻² s⁻¹ (uncorrected mean)** and let the C₄F₆-corrected
mean 4.1 × 10²⁰ sit inside the lower half of the band rather than picking one.

### 3.3 Why Anchor B is excluded from the band (and the disagreement it exposes)

The F/ion normaliser was tested and **rejected**. `[INF]`
```
Wang  F/Σions  = 3.2 / 10.278 = 0.311
Huang F/Σions ≈ 0.8e17 / ~6e15–1e16 = 8 – 13
```
These differ by **25–40×**, so F/ion does not transfer between HPEM cases and
must not be used. The reason is visible in the neutral/ion ratio:
```
Wang    ΣCxFy / Σions = 32.10 / 10.278 = 3.12
Krüger  ΣFC(all) / ions = 3.0940e17 / 1.2e16 = 25.8
Huang   ΣFC / Σions   = 3.0e17 / (6e15–1e16) = 30 – 50
```
**Krüger's discharge is 8× more radical-rich than Wang's and within 1.2–1.9× of
Huang's.** That is the quantitative justification for weighting Anchor A and
discarding Anchor B's ratios. Recorded as a **disagreement**, not smoothed.

### 3.4 The band

> **Γ_F (atomic F, thermal, at the wafer plane, Krüger base case
> C₄F₆/Ar/O₂ = 140/100/105 sccm, 10 mTorr, P_lf/P_hf = 8.0/2.5 kW, r = 7.5 cm):**
>
> ### 2 × 10²⁰ — 1 × 10²¹ m⁻² s⁻¹, centre 5.4 × 10²⁰ m⁻² s⁻¹
>
> equivalently **2 × 10¹⁶ — 1 × 10¹⁷ cm⁻² s⁻¹, centre 5.4 × 10¹⁶ cm⁻² s⁻¹**
>
> i.e. **F ≈ 0.5–0.7 × the published O flux ≈ 0.4–0.6 × the published CF₂ flux**,
> and **≈ 0.25 × the reactive fluorocarbon total**.

**What bounds it from BELOW, and why.**
- Anchor C (Sankaran & Kushner 2005, zero-O₂ c-C₄F₈ ICP, F/ΣCₓF_y = 0.0081)
  shows that with *no* oxygen the same code produces F ≈ 1 % of the FC flux;
  applied to Krüger's ΣFC that would be 1.7e18 m⁻² s⁻¹. **That is not the lower
  bound**, because Krüger's feed is 30.4 % O₂ — it is the "oxygen removed"
  limit, and it is stated here so the band cannot be read as insensitive to O₂.
- The operative lower bound **2 × 10²⁰** is set by compounding the C₄F₆
  fluorine-per-carbon penalty (×0.75 → 4.06e20) with a further ×0.5 allowance
  for (i) Krüger's 10 mTorr vs Huang's 25 mTorr and (ii) the Chun 2015 finding
  that in C₄F₈-based discharges F can *fall* with O₂ (§2.3). Rounding
  4.06e20 × 0.5 = 2.03e20 → **2 × 10²⁰**.

**What bounds it from ABOVE, and why.**
- The uncorrected ΣFC-normalised estimate is 5.72e20; allowing ×1.7 for
  Krüger's higher O₂/C and for the possibility that C₃F₄ (9.5e16, declared
  inert in petch but *not* inert in his mechanism) belongs in the denominator
  gives ≈ 1.0e21.
- **Cross-checked against the measured envelope:** 1.0e21 m⁻² s⁻¹ at 500 K
  corresponds to n_F = 1.0e21 / 186.6 = 5.36e18 m⁻³ = **5.4 × 10¹² cm⁻³**.
  That sits *below* the highest C₄F₈-based value in Chun 2015 (8.0 × 10¹² cm⁻³
  at 0 % O₂) and ~1.5–8× below the measured CF₄ RIE range of Jenq 1994
  (0.8–4.2 × 10¹³ cm⁻³). A C₄F₆ discharge is the *least* fluorine-rich of these
  feeds, so the top of the band is consistent with, and comfortably under, the
  measured ceiling. **AGREE.**

**Density equivalents of the band (at 500 K, Γ/186.6):**

| band point | Γ (m⁻² s⁻¹) | n_F (cm⁻³) @300 K | @500 K | @800 K |
|---|---|---|---|---|
| low | 2.0e20 | 1.38e12 | 1.07e12 | 0.85e12 |
| centre | 5.4e20 | 3.74e12 | 2.89e12 | 2.29e12 |
| high | 1.0e21 | 6.92e12 | 5.36e12 | 4.24e12 |

All inside the measured/model literature range for fluorocarbon discharges
(§2.2, §2.3). **The band is not exotic; it is unremarkable.**

**Confidence.** Band centre `[INF]` grade: strong (three normalisers, ±5 %,
peer-reviewed same-code anchor cross-checked in two documents). Band *width*
`[INF]` grade: moderate — the ×5 span is driven by the unresolved O₂ sign
(§2.3) and by the absence of any published C₄F₆ HPEM F flux. **This band is a
research input, not a deck constant. It is not landed.**

---

## 4. What the band does to the floor fluorine budget

All inputs re-derived from the committed docs, not from memory. Sources named
per line.

### 4.1 Demand side

```
SiO2 formula-unit density n_SiO2 = 2.2e28 m^-3
   [CHEMISTRY_DECK_DESIGN_2026-07-28.md L19: "densities": {"formula_m3": 2.2e28}]

Krüger target rate = 825 nm / 60 s = 13.75 nm/s
   [VALIDATION_DOSSIER §2 "etch depth 825 ± 5 %"; RESULTS_EARLY_TRANSIENT L44]
R_Krüger = 13.75e-9 m/s × 2.2e28 m^-3 = 3.0250e20 formula units m^-2 s^-1
   [reproduces RESULTS_LIMITING_REGIME_2026-08-05 L32-33 "3.03e20" exactly]

petch ml23 (final mechanism) rate = 6.97 nm/s, 60 s projection 418 nm
   [RESULTS_ML23_GRAY_LAWS §3: "Floor rate 6.97 nm/s"; "Projected 60 s depth 418 nm"]
   [BENCHMARK_CERTIFICATION (f): "depth, feature keV (final) ~418 projected (-49%)"]
R_petch = 6.97e-9 × 2.2e28 = 1.5334e20 formula units m^-2 s^-1

THROUGHPUT DEFICIT = 3.0250e20 - 1.5334e20 = 1.4916e20 units m^-2 s^-1
                   = 6.78 nm/s = 407 nm over 60 s
```

Fluorine demanded per removed formula unit — **two sourced values, both carried,
because the sources disagree and it matters by 2×**:

- **4 F/unit** — Kwon ScD MIT 2004 p.76 stoichiometric anchor
  `SiO2(s) + 4 F(s) → SiF4(g) + O2(g)`
  [`RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md` L784, `[Q]` in that doc].
- **≈2 F/unit** — Krüger's *actual* removal rows produce SiF, SiF₂, SiF₃
  (`SiO2CF(s)+Ar+ → SiF + CO2 + Ar#`, `SiO2CF2(s) → SiF2`, `SiO2CF3(s) → SiF3`,
  `krueger_thesis.txt` L6023-6026, `[Q]`), i.e. 1–3 F/unit, mean ≈2.

```
F demand at Krüger's rate : 2F → 6.050e20 ; 4F → 1.210e21 m^-2 s^-1
F demand at petch's rate  : 2F → 3.067e20 ; 4F → 6.134e20
F DEFICIT to be supplied  : 2F → 2.983e20 ; 4F → 5.966e20 m^-2 s^-1
```

**Scale check `[INF]`:** the band centre 5.4e20 m⁻² s⁻¹ is *the same order as
the entire fluorine deficit*. This is not a rounding term.

### 4.2 Transport of F to the floor — with the estimator validated first

Coburn–Winters / Gottscho conductance form, already in the repo as element E7
(Huard Eq. 4.1, `RESEARCH_NEUTRAL_LIMITED_REGIME` §5.2):

    Γ_floor / Γ_wafer = K / (K + S(1 − S_geom-correction))  →  K / (K + S(1−K))

with K the geometric (Clausing) transmission and S the per-strike loss.

**K calibrated on petch's own certified transmission gates**
[`BENCHMARK_CERTIFICATION` §(c): AR 50 → 0.025287; AR 100 → τ = 0.013]:
```
1.27 / 50  = 0.0254   vs petch 0.0253   (0.4 % )
1.27 / 100 = 0.0127   vs petch 0.013    (2.3 % )
=> K(AR) ≈ 1.27/AR ;  K(21.4) = 0.0593
```
AR 21.4 is petch's own Krüger-cell total aspect ratio
[`RESULTS_E8_COUPLED_2026-08-05` §4 table: "Krüger cell | oxide AR 12 | mask
0.85 µm | total AR 21.4"].

**Estimator validated against a petch measurement `[INF]`:** same table gives
plasma CF₂ delivered at that geometry = 1.485e20 from a source of 9.4e20, i.e.
ratio 0.1580. Inverting the formula: S_eff(CF₂) = 0.336 — physically sensible
(bare-oxide chemisorption 0.278 plus polymer deposition). Forward-checking with
the literal chemisorption value S = 0.278 gives 0.185 vs measured 0.158, i.e.
the estimator is accurate to **1.17×** here. Every F number below therefore
carries a stated ±20 % transport uncertainty.

**F delivery**, using the sidewall loss bracket established in §2.5
(10⁻³ measured chamber wall — 10⁻² Krüger's own polymer rows; I use 0.02 as the
central value and show the sensitivity):

| S_wall | Γ_floor/Γ_wafer | ratio to CF₂'s 0.158 |
|---|---|---|
| 0.01 | 0.863 | **5.5×** |
| **0.02** | **0.759** | **4.8×** |
| 0.03 | 0.678 | 4.3× |
| 0.10 | 0.387 | 2.4× |

**This is the physically important structural point.** Atomic F reaches the
floor of an AR-21 feature **2.4–5.5× more efficiently than CF₂ does**, and the
advantage *grows* with AR because K falls linearly while the low-S denominator
does not. Any missing-F error is therefore an **AR-amplified** error — small on
a blanket, largest exactly where petch's depth gate fails.

### 4.3 Reacted F at the floor and the depth impact

Reaction at the floor uses **Krüger's own published probability for the
complex-fluorination rows, p = 0.1** (§0, `krueger_thesis.txt` L6548-6555).
No new constant.

```
Γ_floor       = Γ_wafer × 0.759          (S_wall = 0.02)
Γ_reacted     = Γ_floor × 0.1
Δv            = Γ_reacted / (N_F × 2.2e28)      [m/s]
```

| band point | Γ_wafer (m⁻²s⁻¹) | Γ_floor | Γ_reacted | Δv @2F | Δdepth/60 s | % of the 407 nm deficit | Δv @4F | % @4F |
|---|---|---|---|---|---|---|---|---|
| low | 2.0e20 | 1.52e20 | 1.52e19 | 0.345 nm/s | 21 nm | **5.1 %** | 0.173 | 2.5 % |
| centre | 5.4e20 | 4.10e20 | 4.10e19 | 0.932 nm/s | 56 nm | **13.7 %** | 0.466 | 6.9 % |
| high | 1.0e21 | 7.59e20 | 7.59e19 | 1.726 nm/s | 104 nm | **25.4 %** | 0.863 | 12.7 % |

> **Result: the sourced atomic-F band recovers 2.5 – 25 % of the missing depth
> (central estimate ≈ 7 – 14 %), entering through Krüger's own published
> fluorination rows at his own probability.**

**It does not close −49 %.** That is a genuine negative result and it is the
most useful sentence in this note: it removes "the missing F flux" from the
list of candidate single-cause explanations for the depth gate, exactly as E8
was removed by `RESULTS_E8_COUPLED`. It also **beats E8 by two orders of
magnitude** — E8 moved the floor rate −0.42 % at its physical upper bound
[`RESULTS_E8_COUPLED` §5], where atomic F moves it +5 to +25 %.

### 4.4 The coverage route, and why it also fails to close the gap

The second way F could act is not stoichiometric but through coverage: Gray's
β_e = 2.99 SiF₄/ion applies "from fluorine **saturated** surface regions", so the
question is whether the band saturates the floor.

```
Floor energetic flux = (2.2 – 2.9) × 9.6e19 = 2.11e20 – 2.78e20 m^-2 s^-1
   [RESULTS_FLOOR_DELIVERY_2026-08-05: funnelled delivery 2.2–2.9× source;
    ION_SRC = 9.6e19 m^-2 s^-1 per RESEARCH_NEUTRAL_LIMITED_REGIME L163]
```

| band point | Γ_F,floor | F / energetic at the floor |
|---|---|---|
| low | 1.52e20 | 0.55 – 0.72 |
| centre | 4.10e20 | 1.47 – 1.94 |
| high | 7.59e20 | 2.73 – 3.60 |

Against Gray's own printed half-rise of **96–127** in F/Ar⁺
[`BENCHMARK_CERTIFICATION` §(e) half-rise note], a floor ratio of 0.55–3.6 is
**27–230× below half-rise**. **Atomic F at any point in the sourced band cannot
fluorine-saturate this floor.** The floor's saturation, if it happens, must come
from chemisorbed CFₓ (bare 0.278 / activated 0.85–0.9), which is 3–9× more
sticky per strike and carries 1–3 bound F per arrival. Consistent with
`RESULTS_ML23_GRAY_LAWS` §5: "the chemical channel alone carries the observed
etch **if the floor is fluorine-saturated**".

### 4.5 Cross-channel effects that would have to be graded, not assumed

Landing an F flux is **not** a one-sided depth gain. Krüger's mechanism also
routes F into:
- `CF3(s) + F → CF4` at 0.01 — **polymer removal**, which would thin the lip and
  widen the mouth, moving the `w_m` gate (currently +13 %) and the closure/etch
  ratio (currently 0.1478 vs 0.0310).
- `AC(s) + F → CF(s)` at 0.01 — **mask fluorination**, touching the 850.2 nm
  mask-remaining PASS.
- `SiO2*(s) + F → SiO2(s) + F` at probability **1** — a *deactivation* channel
  that destroys ion-activated sites and therefore **reduces** the activated
  chemisorption path (0.85–0.9), i.e. a term of the **opposite sign** to §4.3.

The §4.3 numbers are therefore an **upper** estimate of the net depth gain
`[INF]`. The deactivation row at p = 1 is potentially large: at the centre band
it destroys activated sites at 4.10e20 m⁻² s⁻¹ against an energetic activation
supply of 2.11–2.78e20 m⁻² s⁻¹ — **F would consume activated sites faster than
ions create them**, which could *reduce* the complex channel by more than the
fluorination rows add. **This is the single most important thing to check before
any implementation**, and it cannot be settled on paper. `[VERIFY]`

---

## 5. What I recommend, in order

1. **Do not land a number from this note.** The band's centre rests on a
   ratio transfer between two different fluorocarbon feeds; that is
   research-grade, not deck-grade. It fails the knob-elimination doctrine's
   "derive, measure, or declare as a fab-measurable constant" test — it is an
   *inferred* constant.
2. **Ask the source** (§6). One number from Krüger/Schulze converts the whole
   of §3 from `[INF]` to `[Q]` and costs an email.
3. **Meanwhile, run the falsifiable 0-D probe that costs nothing**: with the
   frozen-geometry forecaster, sweep an injected Γ_F over the band
   [2e20, 1e21] m⁻² s⁻¹ with **all thirteen** of Krüger's F rows enabled
   (including the p = 1 deactivation row), and report the net floor rate. The
   preregistered prediction from §4.3/§4.5 is: **+5 % to +25 % from the
   fluorination rows, partially or wholly cancelled by the SiO₂\* deactivation
   row; net sign unknown.** If the net is negative, the missing F flux is not
   merely insufficient to close depth — it makes it worse, and that is a
   publishable statement about the published mechanism.
4. **Record, either way, the AR-amplification result** (§4.2): F is delivered
   2.4–5.5× better than CF₂ at AR 21 and the gap widens with AR. That is
   independent of the band's absolute value and is directly relevant to the
   AR-80/AR-200 HARC regime petch is aiming at.

---

## 6. Draft email to Krüger / Schulze — **NOT SENT**

> **To:** Dr. Peter Krüger; Prof. Dr. Sebastian Schulze *(confirm current
> affiliations and addresses before sending — `[VERIFY]`)*
> **Cc:** Prof. Mark J. Kushner *(optional; the HPEM boundary is his code)*
> **Subject:** Atomic-F wafer flux for the Ar/C₄F₆/O₂ base case in JVST A 42, 043008 (2024)
>
> Dear Dr. Krüger, dear Prof. Schulze,
>
> I am writing about your paper "Multiscale modelling of high aspect ratio
> etching of SiO₂ in Ar/C₄F₆/O₂ capacitively coupled plasmas" (*J. Vac. Sci.
> Technol. A* **42**, 043008, 2024) and the accompanying thesis. We have
> implemented your Appendix-B surface mechanism faithfully and without fitting
> any constant to your metrics — the reflection cascade, the deposition-driven
> crosslinking with ion breaking and the per-material bond multiplicities from
> §IV, the substrate-dependent deposition probabilities, the complex-formation
> and complex-removal rows, and the two-component IEAD reconstructed from your
> Fig. 4 — and it reproduces your feature closely: 850.2 nm of remaining mask
> against your 850 nm, an aperture of 41.8 nm at 270 nm depth against your
> 41.1 nm, a sealed feature at O₂/C₄F₆ = 0.5 and an open one at 2.5, and the
> low-frequency power ratios r(4/6) = 0.90 and r(8/6) = 1.02 inside their
> published bands. Our remaining discrepancy is etch depth, and our analysis
> localises it to the fluorine budget and coverage at the etch front rather than
> to ion yield or to transport.
>
> That brings me to the request. Table I of the paper (Table 6.1 of the thesis)
> reports C₃F₄, C₂F₃, CF, CF₂, CF₃, O and the total ion flux, and states that
> these are the species "most responsible for polymer deposition"; atomic
> fluorine is not listed. The Appendix-B mechanism, however, consumes gas-phase
> atomic F in thirteen rows — the fluorination of the passivated complex
> (SiO₂CF/SiO₂CF₂/SiO₂C₂F₃/SiO₂C₃F₅ + F, p = 0.1, and their activated twins),
> the fluorination of the amorphous-carbon mask and of the polymer
> (AC(s)/CF(s)/CF₂(s)/CF₃(s) + F, p = 0.01), and the deactivation of ion-excited
> SiO₂\*(s) by F at probability 1. Without a published Γ_F these rows are
> inactive, which we suspect materially changes the fluorine supply at the etch
> front and, through the deactivation row, the activated-complex population as
> well. **Would you be willing to share the HPEM-predicted atomic-fluorine flux
> to the wafer at r = 7.5 cm for the base case (C₄F₆/Ar/O₂ = 140/100/105 sccm,
> 10 mTorr, P_lf/P_hf = 8.0/2.5 kW)?** If it is easy to include, the F₂ and F⁺
> fluxes and the value of Γ_F across the O₂/C₄F₆ = 0.5–2.5 transfer sweep would
> be equally valuable, since the sign of dΓ_F/d(O₂) in C₄F₆ chemistry is not
> settled in the literature we have found. A single number would be enough for
> our purposes; we currently infer 2 × 10¹⁶ – 1 × 10¹⁷ cm⁻² s⁻¹ by ratio
> transfer from Huang *et al.*, JVST A **37**, 031304 (2019), and we would much
> rather cite you than infer.
>
> We are happy to share what we have found in return: a per-channel comparison
> of the Appendix-B SiO₂ ion rows against Gray's absolutely-calibrated beam
> yields (MIT thesis, 1993), which indicates that the two SiO₂ rows carrying
> n = 1 behave differently from his measured √E laws at multi-keV front
> energies; a quantified delivery advantage of atomic F over CF₂ at AR ≈ 20
> arising from its much lower surface loss probability; and our reconstruction
> of the Fig. 4 IEAD, including an axisymmetric lift correction. We would of
> course cite the paper and thesis for every constant used, and we are glad to
> acknowledge any data you provide in whatever form you prefer. Thank you for
> publishing the mechanism in full — it is rare, and it is what made this
> comparison possible at all.
>
> With best regards,
> *(signature block)*

**Notes before sending:** (i) verify the paper's exact title and both authors'
current affiliations `[VERIFY]`; (ii) decide whether to cc Prof. Kushner — the
HPEM is his group's code and the request is for a reactor-model output;
(iii) the offer in the last paragraph should only be made if we are prepared to
follow through.

---

## 7. Single-source, unverified — quarantined from the band

Nothing in this table is used in §3 or §4. Listed so the note is auditable.

| item | source | why quarantined |
|---|---|---|
| Wang thesis Table I (F = 3.2e15 cm⁻² s⁻¹, Ar/C₄F₈/O₂ 80/15/5, 40 mTorr) | `wang_mingmei_phd_thesis.txt` L2849-2887 | single document; its F/ion ratio disagrees with Huang's by 25–40× (§3.3); its discharge is 8× less radical-rich than Krüger's |
| Jenq 1994 RIE CF₄ [F] = (0.8–4.2)e13 cm⁻³ | abstract page only, `[Q-relay]` | body not read; different chemistry (CF₄, F/C = 4); used only as an upper envelope |
| Kawai/Sasaki/Kadota 1997 actinometry-vs-VUV factor-2 caveat | abstract page only, `[Q-relay]` | body not read; used only as a caution |
| Sasaki 1997 F wall loss "order of 10⁻³" | OpenAlex abstract for DOI 10.1063/1.366495, `[Q-relay]` | body not read; chamber wall, not a feature sidewall; used only as the lower bracket on S, alongside Krüger's own 0.01 rows |
| Kwon & Sawin 2006 "atomic F flux identified as the dominant etch driver" | `RESEARCH_BEAM_CONSTANTS_ATLAS` entry + abstract, `[Q-relay]` | body not read for this note; qualitative use only |
| C₆F₆/Ar/O₂ CCP OES intensity ratios CF₂/F = 0.8–1.1 at 30 mTorr (Materials **15**, 1300, 2022) | WebFetch summary of the PMC article, `[Q-relay]` | OES *intensity* ratios, not densities; no cross-section correction; different feed gas |
| Baek/Efremov 2023 (*Materials* **16**, 5043, DOI 10.3390/ma16145043) F trends in CF₄/CHF₃/C₄F₈ + O₂ | full text fetched via Europe PMC; numeric F values live in figures, not text | no numeric value extractable without digitising figures |
| Cunge, Chabert & Booth, *J. Appl. Phys.* **89**, 7750 (2001) — absolute [F] in a CF₄ CCP from CF₂ loss kinetics + LIF | abstract verbatim via OpenAlex `[Q]`; **body inaccessible** (publisher 403, HAL mirror behind Anubis) | the method is exactly right (CCP, absolute, LIF-anchored) but no numeric [F] was read. **Highest-value single follow-up if a library copy is available.** |
| Krüger Fig. 16(a) transfer-sweep fluxes | digitiser `scripts/digitize_krueger_2024_transfer_boundary.py` colour masks | confirms only that **no F series exists** in that figure (7 series: CF₂, C₃F₄, C₂F₃, CF, O, ions, CF₃) |

---

## 8. Reproduction

Local primary sources read for this note:
- `research_sources/thesis_extracts/krueger-2024.txt` (L266-297, L945-951)
- `research_sources/thesis_extracts/krueger_thesis.txt` (L4290-4360, L4810-4825, L5895-5899, L6021, L6548-6555)
- `research_sources/thesis_extracts/huang_thesis.txt` (L5257-5290, L5440-5460, L10140-10240)
- `research_sources/thesis_extracts/wang_mingmei_phd_thesis.txt` (L2540-2560, L2849-2905)

Downloaded and converted for this note (scratchpad, not committed):
- `JVSTA_37_031304_2019.pdf` — Huang et al. 2019, open archive copy
- `jap_97_023307_2005.pdf` — Sankaran & Kushner 2005, open archive copy
- Chun/Efremov/Yeom/Kwon, *Thin Solid Films* **579**, 136 (2015), open author copy
- Europe PMC full text for *Materials* **16**, 5043 (2023)

All arithmetic in §2.1, §3.2, §4.1–4.4 is a single short Python evaluation; the
key intermediate results are printed inline in each section so they can be
checked by hand.
