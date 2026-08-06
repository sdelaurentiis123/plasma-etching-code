# SF6/O2 relevance pass — what the archived corpus gives the silicon arm (2026-08-06)

Audit and spec only; no physics changed this pass. Every number below is either
quoted verbatim from an archived source with a line reference (`[Q]`), verified by
my own arithmetic on a quoted table (`[Q+calc]`), or my inference (`[INF]`).

Motivating question: the fluorocarbon/SiO2 campaign closed with depth bounded-open,
owned by a floor fluorine-supply term with no published constant. **Does that class
of problem exist for SF6/O2 — the partner-relevant chemistry?** Answer in §5: no,
and the arithmetic says why.

---

## 1. Gray's silicon-side pair, extracted

`research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt`.

### 1.1 The co-regression statement (L1460-1466) `[Q]`

> "where B, is regressed from the intercept of a plot of 1/C versus 1/R, and s, is
> derived from the slope. Parameters regressed from available yield data are
> presented in Tables 5-9 and 5-10. We found that to a good approximation, the
> available etching yield data could be well represented by allowing B, to vary,
> while setting **s,=0.2 and 0.02 for the cases of silicon and SiO, etching
> respectively**."

So the silicon half of the co-regressed pair is `s0_Si = 0.2` — the same statement
that supplied `0.02` to the oxide arm (landed as `_THERMAL_F_STICKING`).

### 1.2 Table 5-9, Ar+/F-Si model parameters (L1490-1497) `[Q]`

| E (eV) | b | s0 (fitted) | s0 (fixed) | B0 (fitted) | B0 (at s0=0.2) |
|---|---|---|---|---|---|
| 20 | 0.040 | 0.216 | (0.2) | 2.46 | (2.5) |
| 150 | 0.110 | 0.203 | (0.2) | 6.95 | (7.0) |
| 250 | 0.142 | 0.308 | (0.2) | 7.71 | (8.5) |
| 350 | 0.168 | 0.231 | (0.2) | 9.58 | (11.0) |
| 500 | 0.201 | 0.126 | (0.2) | 15.20 | (14.0) |
| 1000 | 0.284 | 0.177 | (0.2) | 21.19 | (21.0) |

Footnote `[Q]`: "values in parentheses were adjusted to fit based on the constant
sticking coefficient of 0.2". The parenthetical column is therefore the one paired
with `s0 = 0.2`, exactly as the oxide column was paired with `0.02`.

### 1.3 The silicon energy law `[Q]` + `[INF]`

Gray states the scaling verbatim (L1521-1523):

> "In the case of silicon etching, B, is found to scale well with root ion energy
> from 1 keV down to **a threshold value near 4 eV**. It is significant that the
> Ar'/F etching threshold values predicted are on the order of Si-F bond energies,
> versus physical sputtering thresholds"

My free regression of the parenthetical B0 column against sqrt(E) returns slope
**0.6872** and threshold **4.06 eV** `[INF]` — i.e. it recovers Gray's stated ~4 eV
without being told it. The silicon analogue of the landed oxide law is therefore

```
B0_Si(E) = 0.687 * (sqrt(E) - sqrt(4))          [INF, regressed from Table 5-9]
B0_SiO2(E) = 0.053 * (sqrt(E) - sqrt(4))        [already landed, verified below]
```

Silicon's prefactor is **12.9x** the oxide's — the quantitative form of "silicon
etches much faster than oxide at the same F/ion ratio".

### 1.4 Verification of the already-landed oxide law `[Q+calc]`

The oxide law in the codebase reproduces Gray's own Table 5-10 parenthetical column:

| E (eV) | Table 5-10 | 0.053(sqrt(E)-2) |
|---|---|---|
| 20 | 0.13 | 0.13 |
| 150 | 0.55 | 0.54 |
| 250 | 0.60 | 0.73 |
| 350 | 0.85 | 0.89 |
| 500 | 1.10 | 1.08 |
| 2000 | 2.25 | 2.26 |

Five of six points inside 4%; the 250 eV point is the outlier in Gray's own data
(his fitted s0 there, 0.022, is his largest deviation from 0.02). An independent
free regression of that column gives slope 0.0532 — the landed 0.053 to three
figures. **The oxide law now has a second, independent confirmation.**

### 1.5 Branching ratios — OCR repaired by arithmetic `[Q+calc]`

Eqs. (5-27)/(5-28) are OCR-mangled in the extract ("bs, = 9x10° (E,)"?"). Dividing
the tabulated b columns by sqrt(E) recovers them unambiguously:

- Si: b/sqrt(E) = 0.00894, 0.00898, 0.00898, 0.00898, 0.00899, 0.00898 → **b_Si = 9e-3 sqrt(E)**
- SiO2: b/sqrt(E) = 0.00693, 0.00702, 0.00702, 0.00700, 0.00702, 0.00700 → **b_SiO2 = 7e-3 sqrt(E)**

Both match the OCR'd leading digits (9x10^-3, 7x10^-3). Constant to <1% across a
50x energy range, so the functional form is confirmed, not assumed.

### 1.6 The site-balance and half-rise forms `[Q]`

Eq. (5-30), the fluorinated-site fraction:

```
theta_F = s0*R / (s0*R + 2*B0*(1+b))        R = Gamma_F / Gamma_ion
```

giving a half-rise at `R_1/2 = 2*B0*(1+b)/s0` — the closed form used in §5.

### 1.7 What Gray's Si model does NOT give us

- **Angular dependence.** Verbatim (L1517-1519): B0 "may also be expected to vary
  within a factor of two or so over ion incidence angles of 0-60°" — an estimate,
  not a measurement. No Si angular curve.
- **Reactive-ion extension.** Verbatim (L1520): "Extension of this concept to
  reactive ion systems is not straight-forward and is beyond the present scope".
  Gray's ions are Ar+; SF6 plasmas deliver SF_x+ .
- **Temperature.** Gray is room-temperature beam work (see §4 for why this matters).

---

## 2. SF6 and SiO_xF_y passivation in the archived corpus — a negative result

Counted across all 24 archived extracts:

| file | SF6 mentions | passivation mentions |
|---|---|---|
| `deboer-2002.txt` | 5 | many (the experiment itself) |
| `mask_geometry_micromachines_2023.txt` | 18 | — |
| huang / huard / zhang / Qu / song / others | 0-3 | Cl2- or FC-context |

Regex for `siof|sioxfy|SiO_xF|oxyfluoride` across all extracts: **zero hits.**
Huard has 254 "passivation" mentions but one SF6 mention — his passivation chemistry
is Cl2/Si.

**Conclusion:** no Kushner-lineage thesis in the archive carries an SF6/O2 SiO_xF_y
passivation mechanism. The passivation half of the SF6/O2 model has no beam-measured
source in our corpus and remains Belen-lineage (profile-fitted). This is the single
largest sourcing gap for a provenance-graded Deck 2.

De Boer names the layer explicitly (`deboer-2002.txt` L721) `[Q]`: "sivation layer,
probably composed of SiO F ." — note *probably*; even the validation experiment does
not pin the composition.

---

## 3. The SF6/O2 arm as it stands

### 3.1 Two independent code paths (verified)

```
legacy PoC path :  scripts/* -> src/petch/belen.py -> src/petch/chemistry.py -> src/petch/params.py (PAR)
common engine   :  scripts/deboer_feature3d.py -> src/petch/silicon_sf6o2.py (BelenSiliconParameters)
```

`belen.py` and `chemistry.py` import nothing from `mixed_layer`; `grep` for
`mixed_layer` in `belen.py`, `chemistry.py`, `params.py` returns **zero hits**.

### 3.2 Current constants and their provenance

| quantity | value | where | provenance |
|---|---|---|---|
| F sticking on Si | **0.7** | `params.py:betaE` | ViennaPS `psSF6O2Etching.hpp` default (code, not measurement) |
| F sticking on Si | **0.20** | `params.py:s_F` | "PoC values" — undocumented |
| F sticking on Si | **0.06** | `deboer_feature3d.py:49` default | **fitted** — selected by RMSE against the de Boer knee (L359-370) |
| F sticking on Si | **0.5** | `deboer_feature3d.py:82` builder default | undocumented |
| F sticking sweep | 0.08 / 0.14 / 0.22 | `deboer_calibrate_predict` | fit grid |
| O sticking | 1.0 (belen) / 0.30 (PoC) | `params.py` | ViennaPS default |
| ion-enhanced yield | `A_ie=7.0`, `Eth_ie=15.0 eV` | `params.py` | ViennaPS default |
| physical sputter | `A_sp=0.0337`, `Eth_sp=20.0 eV` | `params.py` | ViennaPS default |
| passivation-sputter | `A_p=3.0`, `Eth_p=10.0 eV` | `params.py` | ViennaPS default |
| chemical etch const | `k_sigma=300.0` | `params.py` | ViennaPS default |
| passivation const | `beta_sigma=0.04` | `params.py` | ViennaPS default |
| angular class (sputter) | `kress_1999`, B=9.3 | `deboer_feature3d.py:157` | declared |
| angular class (ion-enh.) | `chang_sawin_1997` | `deboer_feature3d.py:159` | declared |
| F per removed Si | 4.0 (bounds fixed 4.0-4.0) | `deboer_feature3d.py` | stoichiometric |
| Knudsen wall loss | `2.9` | `params.py` | **fitted knob**, self-documented: "LEGACY replay-only calibration... the former de Boer wafer/held-out AR40 label is **withdrawn**; this is not evidence" |
| Knudsen front-loss trio | 8.0 / 0.5 / 15.0 | `params.py` | "tuned on the de Boer static-AR floor probe" — **fitted** |

**Five different values (0.06, 0.20, 0.5, 0.7, plus a 0.08-0.22 sweep) are in the
tree for one physical quantity.** They are reachable from different entry points and
none of them is currently traceable to a measurement.

### 3.3 The atlas's L3 verdict, re-verified `[Q]`

Belen et al., JVST A 23, 99 (2005), abstract:

> "experimentally inaccessible parameters such as the **F sticking coefficient,
> chemical etch rate constant, and the ion-enhanced etch yield** are determined by
> matching simulated feature profiles with those obtained from carefully designed
> etching experiments."

The three parameters Belen calls *experimentally inaccessible* are precisely the
three Gray measured by beam experiment twelve years earlier: s0 (§1.1), the thermal
channel k_p, and B0 (§1.2). That is the whole upgrade thesis in one sentence.

### 3.4 The Micromachines 2023 citation, re-verified `[Q]`

Bobinac, Reiter, Piso, Klemenschits, Baumgartner, Stanojevic, Strof, Karner,
Filipovic, "Effect of Mask Geometry Variation on Plasma Etching Profiles",
*Micromachines* **14**, 665 (2023), DOI 10.3390/mi14030665. Parameter table
(L477-489): "F sticking on Si γF = **0.7**" across all four cases; "O sticking on Si
γO = 1.0"; coverage ODEs (L265-278) are the Belen form.

**Correction to how this was cited previously:** these authors are the TU Wien
ViennaPS group. Their 0.7 is *the same ViennaPS default already in `params.py`*, not
an independent measurement of it. It is a second appearance of one number, not a
corroboration. Evidence class: declared simulation input, not L0/L1.

Useful for us anyway: their flux table gives Γ_F = 3.0-5.5e18 cm^-2 s^-1 against
Γ_i = 1e16 cm^-2 s^-1, i.e. **R = 300-550** — an independent statement of how
neutral-rich this chemistry is (§5).

### 3.5 Structural point in the arm's favour

The arm already uses `Y = A*(sqrt(E) - sqrt(Eth))` (`chemistry.py:_ied_yield`) — the
same Steinbruechel form Gray's B0 obeys — and already assigns **Kress class 1 to
physical sputter and Chang-Sawin class 2 to ion-enhanced** (`deboer_feature3d.py`
157-161). That is the per-material/per-channel angular split the fluorocarbon side
spent this week discovering. The SF6 arm had it right first.

---

## 4. Where the newly-sourced constants slot in

| target | current | measured replacement | class change | caveat |
|---|---|---|---|---|
| F sticking on Si | 0.06 / 0.2 / 0.5 / 0.7 | **0.2** (Gray L1463) | L3-fitted / code-default → **L1 beam-regressed** | must land **with** its B0 partner |
| ion-enhanced yield | `A_ie=7.0`, `Eth=15 eV` | **B0 = 0.687(sqrt(E)-sqrt(4))**, Eth ≈ **4 eV** (Gray L1521 + Table 5-9) | ViennaPS default → **L1** | different normalisation; needs a unit bridge before comparison |
| branching ratio | absent | **b = 9e-3 sqrt(E)** (Gray Eq. 5-27, §1.5) | — → **L1** | new channel, not a replacement |
| F wall loss | absent (Knudsen knob 2.9 stands in) | ~1e-3 (Sasaki, JAP 82, 5938) | fitted knob → sourced | `[Q-relay]`, abstract only |
| passivation (SiO_xF_y) | Belen-fitted | **none available** | — | §2: nothing in the corpus |

**The pair rule applies here exactly as it did on the oxide side.** Gray's s0 and B0
are two halves of one regression (§1.1); landing 0.2 without the B0 partner would
repeat the error the oxide pass caught and reverted. `Eth_ie: 15 eV → ~4 eV` is the
larger physical change of the two and must be forecast before it is landed.

**Temperature caveat, and it is not small.** De Boer 2002 is *cryogenic* — title
(L7): "Using Fluorine High-Density Plasmas at **Cryogenic** Temperatures"; the
mechanism (L122) `[Q]`: "atures enhances passivation by reducing the chemical
reactivity". Gray's beam data is room-temperature. Cryo suppresses exactly the
spontaneous/thermal-F channel Gray's s0 governs. **Gray's Si pair therefore transfers
cleanly to room-temperature SF6/O2 (the partner-relevant case) but must not be landed
against the de Boer cryo validation without a declared temperature treatment.** Two
different operating points; the arm currently serves both from one constant set.

---

## 5. Does the fluorocarbon depth-class problem exist for SF6/O2? **No.**

The oxide failure is a supply failure: the floor sits far below the F/ion ratio at
which the chemistry saturates, so removal is ion-limited and cannot respond to
radical delivery. Gray's own half-rise form (§1.6) makes this a one-line test.

**Silicon, at the arm's own delivered ratio** (`params.py`: `Fflux=1800`,
`ionFlux=12` → R = 150) `[Q+calc]`:

| E (eV) | B0_Si | b | R_1/2 = 2B0(1+b)/s0 |
|---|---|---|---|
| 50 | 3.47 | 0.063 | 36.9 |
| 100 | 5.49 | 0.090 | **59.8** |
| 150 | 7.03 | 0.110 | 78.0 |
| 300 | 10.51 | 0.156 | 121.5 |
| 500 | 13.98 | 0.201 | 167.8 |

At the arm's `Emean = 100 eV`: **R = 150 sits 2.5x ABOVE half-rise**, implying
`theta_F = 0.715` — the saturated, neutral-rich branch. Micromachines' independent
flux table (§3.4) gives R = 300-550, further above.

**Oxide, at the Krueger floor** (Gray SiO2, s0 = 0.02):

| E (eV) | R_1/2 |
|---|---|
| 350 | 100.2 |
| 1500 | 247.4 |
| 3406 (measured front) | 420.7 |

against a measured floor ratio of 0.55-3.6 → **117-765x BELOW half-rise.**

**The two chemistries sit on opposite sides of the same measured curve**, and the
separation is ~3 orders of magnitude in R/R_1/2. Two independent reasons, both
physical: SF6 dissociates directly to atomic F (the oxide recipe delivers fluorine
bundled in CF_x, with atomic F absent from Krueger's published boundary entirely),
and s0 is **10x larger** on Si than on SiO2 (0.2 vs 0.02, Gray L1463), which lowers
the half-rise by the same factor.

**Consequences.**
1. The bounded-open depth item does **not** propagate to the partner-relevant
   chemistry. It is specific to fluorocarbon/SiO2 at keV.
2. The SF6 arm should be *more* sensitive to radical delivery, not less — ARDE from
   neutral shadowing should appear naturally, which is consistent with de Boer being
   an ARDE experiment in the first place.
3. **Watch item:** at 500 eV, R_1/2 = 168 > R = 150. High-bias SF6 crosses into the
   starved branch. Any high-bias silicon recipe needs this checked, not assumed.

---

## 6. Deck 2 (SF6/O2, provenance-graded) — what it needs beyond relabelling

**Available now, sourced:**
1. F sticking on Si = 0.2 `[L1, Gray L1463]` — with its B0 partner.
2. Ion-enhanced yield law B0 = 0.687(sqrt(E)-sqrt(4)) `[L1, Table 5-9 + L1521]`.
3. Branching ratio b = 9e-3 sqrt(E) `[L1, Eq. 5-27, arithmetic-verified §1.5]`.
4. Stoichiometry: 4 F per removed Si — already fixed by bounds.
5. Angular classes — already declared correctly (§3.5).
6. F wall loss ~1e-3 `[L2, Sasaki, relay]` — to retire part of the Knudsen knob.
7. Physical sputter threshold ~20 eV — consistent with Gray's remark that sputter
   thresholds exceed the ~4 eV ion-enhanced threshold `[Q, L1523]`.

**Missing, no source in corpus:**
1. **SiO_xF_y passivation kinetics** — the whole O-side (§2). Belen-fitted only. This
   is the blocker; a deck cannot claim provenance grading with a fitted O channel.
2. **SF_x+ vs Ar+ ion chemistry** — Gray explicitly declines the extension `[Q, L1520]`.
3. **Angular dependence for Si** — Gray offers "a factor of two or so", not a curve.
4. **Cryogenic temperature dependence** — required for de Boer specifically (§4).
5. **O sticking on Si** — 1.0 is a ViennaPS default; no measurement located.
6. **Absolute F flux for de Boer's reactor** — same class of gap as the oxide side.

**Recommended sequence** (no physics landed this pass):
1. Unify the five F-sticking entry points onto one declared constant with a
   provenance field — a pure hygiene fix, no behaviour change if each caller keeps
   its current value explicitly.
2. Land Gray's Si pair (s0 = 0.2 **with** B0 and Eth ≈ 4 eV) behind a flag, gated on
   the room-temperature branch only; forecast before any run.
3. Re-grade de Boer under an explicit cryo declaration rather than silently reusing
   room-temperature constants.
4. Deck 2 ships as *partially graded*: Si channels L1, O/passivation channels
   declared L3-fitted with the Belen quote attached. Honest, and better than the
   current uniform L3.

---

## 7. Files inspected

Sources: `research_sources/thesis_extracts/{gray_thesis_1993_ocr_sections,
deboer-2002, mask_geometry_micromachines_2023, huang_thesis, huard_chad_phd_thesis,
wang_mingmei_phd_thesis, zhang_yiting_phd_thesis, Qu_Chenhui_PhD_Thesis_2020,
Lanham_Steven_PhD_Thesis_2022, logue_michael_phd_thesis, song_sangheon_phd_thesis,
tian_peng_phd_thesis, Konina_Kseniia_PhD_Thesis_2024}.txt`

Code: `src/petch/{belen,chemistry,params,silicon_sf6o2,mixed_layer}.py`,
`scripts/deboer_feature3d.py`, `tests/{test_silicon_sf6o2,
test_deboer_direct_validation,test_experimental_data}.py` (40 passed at audit time).
