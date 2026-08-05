# The measured energy dependence of the two SiO2 ion channels

Question posed: what is the MEASURED `Y(E)` for (a) physical sputtering and (b)
chemically-enhanced (F-saturated) etching of SiO2 under Ar+/fluorocarbon, over
~100 eV to >1 keV — and does it follow ZBL deposited energy, the published
linear `(E-Eth)/(Er-Eth)` with n=1, Steinbruchel `sqrt(E)-sqrt(Eth)`, or
something else?

**Answer: it is measured, it is `sqrt(E) - sqrt(Eth)`, and the measurement is in
the same source that gave us the absolute magnitudes.** Gray's own thesis
tabulates the F-saturated SiO2 yield at **six ion energies from 20 eV to
2000 eV** and publishes the fitted law. Nothing here is inferred from a review
or a modelling paper — it is the primary beam data.

---

## 0. New primary source acquired

`research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt` (2063 lines,
committed to the archive, NOT to git yet).

> Gray, D.C., *"Beam simulation studies of plasma-surface interactions in
> fluorocarbon etching of Si and SiO2"*, PhD thesis, MIT Dept. of Chemical
> Engineering (1993). MIT DSpace handle **1721.1/13187**, bitstream
> `26680838-MIT.pdf`, 32.3 MB, 618 pages, image-only (no text layer).

Refetch route (open bitstream, no auth):

```
https://dspace.mit.edu/server/api/core/bitstreams/d1b2c1c0-6709-495b-b725-df3470cf5b2f/content
```

OCR'd with `pdftoppm -r 250 -gray` + `tesseract --psm 6`. **Printed page number
== PDF page number in the body.** Sections captured: 5.1 (155-163), 5.3
(163-190), 5.6 (229-256), 6.4 (300-312). Line refs below are into the archive
file; page refs are the thesis's own.

This is the thesis behind Gray, Tepermeister & Sawin, *JVST B* **11**, 1243
(1993) — the source of the 0.28 / 1.10 absolute yields in
`RESULTS_ABSOLUTE_YIELD_2026-08-05.md`.

---

## 1. The published law, verbatim, from three independent levels

**(a) The journal paper's own abstract** (Crossref, DOI `10.1116/1.586925`,
retrieved verbatim):

> "The etching ''enhancement'' effect of normally incident Ar+ ions has been
> quantified over a wide range of ion energy through the use of Kaufman and
> electron cyclotron resonance-type ion sources. **The increase in per ion
> etching yield of fluorine saturated silicon and silicon dioxide surfaces with
> increasing ion energy (Eion) was found to scale as (Eion^1/2 - Eth^1/2)**,
> where Eth is the etching threshold energy for the process."

**(b) The thesis body** (Gray 1993 p. 252; archive L1612-1618):

> "The E_i^1/2 scaling is also observed in the case of SiO2 etching **over the
> range of 20 eV to 2000 eV**. Best fit lines to this data give the following
> energy dependencies,
>
>     beta_e,Si   = 0.687 (E_i^1/2 - 4^1/2)      (5-34)
>     beta_e,SiO2 = 0.053 (E_i^1/2 - 4^1/2)      (5-35)"

`beta_e` is defined on p. 248 (archive L1519-1520):

> "In the simplest interpretation, **beta_e is the number of SiF4 molecules
> removed from fluorine saturated surface regions per incoming ion** and is a
> strong function of ion energy."

That is exactly our chemically-enhanced channel: yield per ion at saturated
fluorine coverage.

**(c) The canonical form**, Steinbruchel, *Appl. Phys. Lett.* **55**, 1960-1962
(1989), DOI `10.1063/1.102336` — abstract verbatim via Crossref:

> "**Physical and ion-enhanced chemical etch yields are shown to be a linear
> function of the square root of the ion energy down to the etching threshold**
> for self-sputtering of metals, sputtering of metals by noble gas ions,
> **sputtering of Si and SiO2 by noble gas and reactive ions**, and ion beam
> enhanced chemical etching of Si. The threshold energy must be taken into
> account for a quantitative description of etch yields **even at intermediate
> ion energies**. The relationship between the dependence of etch yields on ion
> energy and incident angle is also discussed."

Chang (MIT ScD 1999, `research_sources/chang_thesis.pdf`) writes the same form
explicitly (his Eq. 3.2-3.4, thesis text extract L2266-2280):

> "The universal energy dependence of ion bombardment induced etching processes
> proposed by Steinbruchel [Steinbruchel, 1989] is therefore used to model the
> experimental measurements. ... where the modified nuclear stopping function
> was proposed by Wilson [Wilson, 1977] and f(Eth/Eion) was a function proposed
> by Matsunami [Matsunami, 1980], and the etching yield can be expressed as:
>
>     Y(Eion) = A * (Eion^1/2 - Eth^1/2)        (3.4)"

---

## 2. The measured points

### 2a. Chemically-enhanced (F-saturated) SiO2 — six points, 20 eV to 2000 eV

Gray 1993 **Table 5-10, "Ar+/F-SiO2 Model Parameters"** (p. 247; archive
L1492-1500). Verbatim:

```
Table 5-10. Ar*/F-SiO, Model Parameters
E,        b       S,             B,
20      0.031   0.013 (0.02)   0.143 (0.13)
150     0.086   0.017 (0.02)   0.587 (0.55)
250     0.111   0.022 (0.02)   0.536 (0.6)
350     0.131   0.033 (0.02)   0.740 (0.85)
500     0.157   0.045 (0.02)   0.892 (1.10)
2000    0.313   0.013 (0.02)   2.845 (2,25)
```

Reading: `E_i` (eV) | `b` | `s_0` regressed (value in parens = refit at the
constant `s_0 = 0.02`) | `beta_e` regressed (parens = refit at constant `s_0`).
Gray states the convention on p. 246 (archive L1461-1466):

> "We found that to a good approximation, the available etching yield data could
> be well represented by allowing beta_e to vary, while setting s_0 = 0.2 and
> **0.02** for the cases of silicon and SiO2 etching respectively. ''Best fit''
> values of beta_e using this constant initial sticking coefficient
> approximation are also provided in Tables 5-9 and 5-10."

Two side-notes worth banking:

- `s_0 = 0.02` on SiO2 is the **site-limited F adsorption coefficient** that
  `RESULTS_ION_CHANNEL_SOLVE` §3 inverted from Gray's half-rise as
  `s ~ 0.06` and recorded `[VERIFY]`. Gray's own value is **0.02**, and the
  regressed column shows it rising with ion energy (0.013 -> 0.045 over 20-500 eV).
  The inverted 0.06 is the right order; the source number is **0.02-0.045**.
  This closes the C1b `[VERIFY]` on the low side and it is *smaller* than
  inverted, i.e. the knee should sit even further right.
- `b` is exactly `0.007 * sqrt(E)` at all six energies (0.031/sqrt20 = 0.00693;
  0.313/sqrt2000 = 0.00700). The auxiliary parameter is itself sqrt-scaled.

### 2b. Physical sputtering of SiO2

Gray 1993 **Table 5-1, "Si and SiO2 Low Energy Sputtering Model Parameters"**
(p. 159; archive L132-135), verbatim:

```
Table 5.1 Si and SiO, Low Energy Sputtering Model Parameters
        A          E_th
Si      0.0337     ~20 eV
SiO,    0.0139     ~18 eV [ref. Chapman, 1980]
```

used in his Eq. (5-7) (p. 159, archive L127):

> `R_sp = Y_sp = A (E_i^0.5 - E_th^0.5) I+`

So: **`Y_sputter,SiO2(E) = 0.0139 (sqrt(E) - sqrt(18))` SiO2 per Ar+.**

The supporting data set is Figure 5-2 (p. 162), plotted twice — panel (a) yield
vs `E` over **0-1200 eV**, panel (b) yield vs `E^0.5` over `0-35` (i.e. to
1225 eV) — with four sources: Gray's own 1992 data, Chapman 1980, Oostra 1986,
Steinbruchel 1989. Its caption:

> "Figure 5-2. Comparison of the a) Sigmund [1969] and b) Steinbruchel [1989]
> models for the low energy Ar+ sputtering of SiO2."

**Independent corroboration that `A = 0.0139`, `Eths = 18 eV` is the live
in-chemistry constant:** the TU Wien / ViennaPS fluorocarbon-silica model
carries exactly these numbers. Rodrigues et al., *J. Comput. Electron.* **22**,
1558-1563 (2023) — `research_sources/thesis_extracts/tuwien_rodrigues_2023_fc_silica.txt`
L212-270 — Eqs. (9)-(10) and Table 1:

```
Y_{n,n/p}(E,theta) = A_{n,n/p} (sqrt(E) - sqrt(Eth)) cos(theta)          (9)
Y_s(E,theta)       = A_s (sqrt(E) - sqrt(Eths)) (1 + B sin^2(theta)) cos(theta)  (10)

Eth   4 eV        An    0.0361 eV^-1/2
Eths  18 eV       An/p  0.1444 eV^-1/2
                  As    0.0139 eV^-1/2
                  B     9.3
```

`As = 0.0139`, `Eths = 18` is Gray's Table 5-1 row transported verbatim; `Eth =
4 eV` is Gray's Eq. (5-35) threshold. **Both of ViennaPS's SiO2-in-fluorocarbon
channels are sqrt-form.** (Aside, free: their `B = 9.3` in `(1 + B sin^2 th)
cos th` gives peak/normal = **4.172** — bit-identical to the "Kress B=9.3" shape
in `RESULTS_ION_CHANNEL_SOLVE` §4. Same functional form, two different
attributions; recorded, not chased here.)

### 2c. The fluorocarbon-radical channel — a real caveat, cutting the other way

Gray 1993 **Table 6-1, "Ar+/CF2 SiO2 Etching Model Parameters"** (p. 305;
archive L1873-1878):

```
        s_o     beta_e
150 eV  0.19    0.41
250 eV  0.66    0.45
350 eV  1.50    0.52
```

and p. 306-307 (archive L1909, L1919-1921), verbatim:

> "The values of B_e for the Ar+/CF2 system are only **weakly dependent on the
> ion energy**, much less so than in the case of Ar+/F etching. Due to the size
> of the CF2 radical as compared to F, it has a much lower mobility in the
> ''mixing zone'' ... In fact, the physical sputtering yield, p_s, (Ar+ only)
> shows a **stronger dependence on ion energy than does the parameter B_e**.
> Therefore, as the ion energy increases, **the etching process becomes less
> chemical and more physical**."

Over 150 -> 350 eV that is `0.41 -> 0.52`, a factor **1.27**, against sqrt's
1.63 and linear's 2.74. Fitted as a power law it is `E^0.28` — **sub-sqrt**.
Three points over 2.3x in energy is thin leverage and it is a *different*
mechanism (CF2-transport-limited, not F-saturated), but it is the closest
in-chemistry analogue to Krueger's `SiO2CFxy` complex row and it says the
complex channel's energy dependence is **at most** sqrt, never linear.

Also recorded from the same page (relevant to the class-1/class-2 angular
question, not to this one): Gray citing **Mayer et al. [1981]** on SiO2 etched
by CF3+ — "the angular dependence of the sputtering yield was peaked at normal
incidence for ion energies less than 400 eV, but became peaked at an increasing
angle off normal for energies greater than 400 eV". The angular class is itself
energy-dependent, with the switch at ~400 eV. `[VERIFY]` against Mayer's
original.

---

## 3. Discrimination: which form fits the measured points

Fitted to Gray's six SiO2 `beta_e` points (script run inline; both columns of
Table 5-10 tested because the raw regressed column carries the `s_0` scatter):

**Constant-`s_0` column (0.13, 0.55, 0.60, 0.85, 1.10, 2.25) — Gray's own
preferred fit:**

| form | free-fit parameters | R^2 | RMSE |
|---|---|---|---|
| **`A(sqrt(E) - sqrt(Eth))`** | A = 0.0532, Eth = 6.6 eV | **0.9943** | **0.051** |
| Gray Eq. (5-35) `0.053(sqrt(E)-2)` (no free params) | — | 0.9927 | 0.057 |
| `A(E - Eth)` best possible linear | A = 0.00125, Eth -> 0 | 0.7440 | 0.338 |
| Krueger's published `(E-35)/(350-35)`, anchored at 350 eV | none | **-2.532** | 1.254 |
| petch's ZBL-deposited-in-layer, anchored at 350 eV | none | **-1.843** | — |

**Raw regressed column (0.143, 0.587, 0.536, 0.740, 0.892, 2.845):**

| form | R^2 |
|---|---|
| `A(sqrt(E)-sqrt(Eth))` free | 0.9543 |
| `A(E-Eth)` free | 0.9436 |
| Gray Eq. (5-35) | 0.9055 |
| Krueger `(E-35)/(350-35)` anchored at 350 | 0.2792 |
| petch ZBL-in-layer anchored at 350 | 0.4813 |

Honest reading: on the raw column a free-parameter linear fit is nearly as good
as sqrt (0.944 vs 0.954), because the 2 keV raw point (2.845) sits above the
sqrt curve. On Gray's own preferred column the sqrt form wins decisively
(0.994 vs 0.744). **In both columns the two forms petch actually runs — the
Appendix-B linear anchored at its own threshold, and the ZBL-deposited-in-layer
shape — score negative or near-zero R^2**, because both predict roughly 5x at
2 keV where the measurement is 2.25-2.85.

Concretely, at the one point that carries the extrapolation:

| | measured 2000 eV | predicted from the 350 eV value |
|---|---|---|
| Gray, `beta_e` | **2.25** (adj) / 2.845 (raw) | — |
| sqrt, Eth = 4 | | 2.18 |
| Krueger linear, Eth = 35, n = 1 | | **5.30** |
| petch ZBL-in-layer | | **4.99** |

**The linear and ZBL forms over-predict the measured 2 keV yield by 1.8-2.4x.**
That is a direct, in-material, in-chemistry refutation, at an energy only 1.7x
below our feature front.

### Why "ZBL vs published-linear" was never the real question

petch's `_complex_energy_factor` uses `eps(E)/eps(140)` where `eps` is
`nuclear_energy_in_layer_eV` evaluated over the **projected range** — which
grows with E. Measured directly:

| E (eV) | Rp (nm) | eps_dep (eV) | eps_dep / E |
|---|---|---|---|
| 70 | 0.677 | 51.7 | 0.738 |
| 140 | 1.107 | 114.0 | 0.814 |
| 350 | 1.969 | 302.7 | 0.865 |
| 1500 | 4.747 | 1333.8 | 0.889 |
| 3406 | 7.988 | 3013.0 | 0.885 |

Because the integration depth tracks the range, this recovers ~0.885 x E and is
**linear in E to within 2% above ~700 eV**. petch's "ZBL" factor and Krueger's
`n = 1` are the *same family*: 9.95x vs 10.70x over 350 -> 3406 eV, a 7%
difference. The Sigmund sputter driver is not the range-integrated energy, it is
the deposited-energy **density** in the escape layer. Evaluating the same ZBL
module over a *fixed* near-surface layer:

| fixed layer | eps(350) | eps(3406) | ratio |
|---|---|---|---|
| 0.5 nm | 108.2 | 251.8 | **2.33** |
| 1.0 nm | 196.2 | 498.5 | **2.54** |

which matches an analytic ZBL reduced-stopping evaluation
(`sn(eps)` at `eps = 1.407e-5 * E`; ratio `sn(0.0479)/sn(0.00492) = 2.23`).
So the genuine ZBL/Sigmund answer is **~2.2-2.5x**, not 9.95x. The naming in
`mixed_layer.py:274` ("ZBL deposited-energy shape") is accurate about the module
but misleading about the physics: **it is a linear-in-E factor wearing a ZBL
label.** That is a finding in its own right.

### The form the literature says is wrong, and why it looked right

Gray states the linear-in-E origin explicitly (p. 157-158, archive L80-83):

> "When the incoming particle energies are <1 KeV, the binary particle
> interactions are well described by a Born-Mayer potential, and **sputtering
> yield is found to be linear in the incoming particle energy, E_i**."

and then (p. 158, archive L101-102):

> "what is most important to note is that **the linear cascade theory does not
> accurately reflect measured physical sputtering yields when the incoming ion
> energies fall below 200 eV or so**. This shortcoming is due to the lack of
> treatment of sputtering ''threshold'' effects..."

and his verdict for Si (p. 161, archive L140) and for SiO2 (p. 161, archive
L199-202):

> "It is readily apparent that the available data is more accurately
> represented by equation (5-7) [the sqrt form]..."
>
> "As demonstrated in Figure 5-2, **equation (5-7) was again found to give a
> much more accurate representation of SiO2 sputtering yields as a function of
> incoming Ar+ ion energy. Agreement of our sputtering data with that found in
> the literature is excellent.**"

**Gray ran exactly the experiment we are arguing about — linear-in-E (Sigmund)
against sqrt (Steinbruchel), on SiO2, over 0-1200 eV — and rejected linear.**

### An internal-consistency finding in Krueger's own mechanism

`research_sources/thesis_extracts/krueger_thesis.txt`. The Appendix B column
order is stated at L5329-5345: `p0 | eps_th | n | eps_0 | angular-class`.
Counting every energy-dependent row in Appendix B (L5345-6100):

| `n` | rows | which |
|---|---|---|
| **0.5** | **416** | everything else (all polymer / CFx / crosslink / mask rows) |
| 1 | 144 | **only** `0.0852 70 1 140 1` (68 rows, bare SiO2) and `0.1471 35 1 140 2` (76 rows, SiO2CFxy complex) |

**The two rows this investigation is about are the only two in the entire
mechanism that depart from n = 0.5.** And three Kushner-group theses state the
group default in words:

- Huang L2287-2289: *"where Eth is the threshold energy, Er is a reference
  energy, p0 is the yield at the reference energy, **n is the energy dependent
  exponent (typically 0.5)**"*
- Qu L2698-2699: *"The **exponential term n is typically 0.5**."*
- Huard L2381-2382: *"p0 is the sputtering yield at the reference energy and **n
  is the energy dependency exponent (typically 1/2)**."*

Huang's own SiO2 mechanism table (Table E.2, L10194-10200 and L10362-10363)
lists `Eth = 70 / Er = 140` for SiO2 sputtering and `Eth = 35 / Er = 140` for
the complex — the same thresholds Krueger uses — but **prints no `n` column at
all**, i.e. the group default 0.5. The `n = 1` appears to enter at Krueger.

And Krueger flags the consequence himself (thesis L4884-4888):

> *"the effect of ion energy (for example in sputter yield or related processes)
> might be overestimated in the mechanism"*

Zhang (Kushner group) states the scaling and cites it — L5154-5156, and the
reference resolves at L5519:

> "The etch rate of SiO2 generally scales as `(eps_ion^1/2 - eps_th^1/2)`, where
> eps_th is a threshold energy that depends on the details of the chemical
> system.[28]"
>
> "28. D. C. Gray, I. Tepermeister and H. H. Sawin, J. Vac. Sci. Technol. B 11,
> 1243 (1993)."

---

## 4. The 350 -> 3406 eV transfer factor, per form

| form | Eth | transfer 350 -> 3406 eV |
|---|---|---|
| Krueger linear n=1, bare | 70 | **11.91** |
| Krueger linear n=1, complex | 35 | **10.70** |
| petch ZBL-deposited-in-layer (range-integrated) | — | **9.95** |
| Steinbruchel sqrt, Krueger's bare threshold | 70 | 4.83 |
| Steinbruchel sqrt, Krueger's complex threshold | 35 | 4.10 |
| **Steinbruchel sqrt, Gray's own sputter threshold** | **18** | **3.74** |
| **Steinbruchel sqrt, Gray's own chem threshold** | **4** | **3.37** |
| sqrt to 1 keV then true ZBL surface-density above | 18 / 4 | 2.78-2.86 / 2.60-2.68 |
| pure ZBL/Sigmund surface deposited-energy density | — | 2.33-2.54 |

The **form** moves the answer by ~3x; the **threshold** moves it by ~30%. Any
argument about Eth is second-order to the argument about n.

### Carrying Gray's measured magnitudes to the feature front

Gray's laws, in absolute units (SiO2 per incident Ar+):

```
Y_sputter(E)  = 0.0139 (sqrt(E) - sqrt(18))                 Table 5-1
beta_e(E)     = 0.053  (sqrt(E) - sqrt(4))                  Eq. (5-35)
b(E)          = 0.007 sqrt(E)                               Table 5-10, exact
Y_ion-enh(E)  = beta_e (1 + b)                              Eq. (5-31) high-R limit
```

| | petch | Gray (measured law) | petch / Gray |
|---|---|---|---|
| **at 350 eV** | | | |
| bare / physical sputter | 0.341 | **0.201** | **1.69x too strong** |
| complex / ion-enhanced | 0.391 | **1.002** | **0.39x — 2.6x too weak** |
| (complex vs total saturated yield 1.203) | 0.391 | 1.203 | 0.32x |
| bare / complex | 0.872 | 0.201 | 4.3x inverted |
| **at 3406 eV (the front)** | | | |
| bare / physical sputter | 4.060 | **0.752** | **5.40x too strong** |
| complex / ion-enhanced | 3.890 | **4.207** | **0.92x — essentially right** |
| (complex vs total saturated yield 4.960) | 3.890 | 4.960 | 0.78x |
| bare / complex | 1.044 | **0.179** | **5.8x inverted** |

Two notes on the anchors:

1. `RESULTS_ABSOLUTE_YIELD` reads the Gray/Kwon floor as **0.28** at 350 eV.
   Gray's own sputtering model gives **0.201** at 350 eV, so the bare row is
   **1.69x** too strong against the primary source, not 1.22x against the
   replot. The 0.28 reading corresponds to `Y_sputter` at ~600 eV under Gray's
   own constants. Recommend re-reading Kwon Fig. 3.4's stated energy — the Kwon
   thesis DSpace text layer is unusable OCR, so this stays `[VERIFY]`.
2. The 1.10 plateau anchor is **confirmed and located**: Gray's total saturated
   SiO2 yield at 350 eV is `beta_e(1+b) + Y_sputter = 0.885*1.131 + 0.201 =
   1.20`, and the constant-`s_0` table entry at 500 eV is literally `1.10`.
   The plateau is real; the 350 vs 500 eV attribution should be tightened.

**Under Gray's own measured energy laws, the bare/complex balance at the feature
front inverts from 1.04:1 to 0.18:1 — the chemical channel becomes ~5.6x the
physical one.** The fork's forecast (0.30:1 under sqrt with Krueger's
thresholds) was directionally right and conservative; with Gray's own constants
the inversion is stronger.

---

## 5. Above 1 keV: is there saturation or roll-over?

This is the honest extrapolation limit, and it is smaller than feared.

**Measured coverage.**
- Chemically-enhanced SiO2: Gray's own data runs **to 2000 eV** and the sqrt law
  holds there (p. 252, quoted above). 3406 eV is only **1.70x** beyond his top
  measured point.
- Physical sputter of SiO2: Gray's Figure 5-2 spans **0-1225 eV** (Gray 1992,
  Chapman 1980, Oostra 1986, Steinbruchel 1989). 3406 eV is 2.8x beyond.
- Gray's beam itself ran 20-500 eV (p. 155, archive L9; and his Table 5-3 header, archive L1081: "E_i = Ion energy (20-500 eV)"); the
  1000 eV and 2000 eV points come from Gerlach-Meyer 1981 and Tu et al. 1981
  respectively (Figs. 5-5, 5-13).

**Roll-over.** No source shows the SiO2 yield saturating or falling below
~10 keV. Literature compilations put the Ar+ -> SiO2 yield **maximum at ~10-20
keV** (~1.8 atoms/ion), and Ar+ -> Si at ~10 keV; measured yields agree with
Sigmund linear-cascade theory through the intermediate range. So **3406 eV is on
the rising, flattening part of the curve** — the sqrt form does not break down
qualitatively, it just becomes optimistic. `[VERIFY]` on the exact 10-20 keV
maximum: it comes from search summaries of Seah et al. 2010 / Radiation Effects
21 (1974), not from a document I read end-to-end.

**How optimistic.** The ZBL reduced nuclear stopping is `∝ sqrt(eps)` only while
`eps << 0.01`. For Ar -> SiO2, `eps = 1.407e-5 * E`, so `eps = 0.0049` at 350 eV
(sqrt regime, good) and `eps = 0.0479` at 3406 eV (already flattening). The
surface deposited-energy-density ratio is 2.33-2.54 against the sqrt form's
3.37-3.74, i.e. **the sqrt law over-predicts by ~1.4-1.6x at 3406 eV.**

Two independent statements of the ceiling:

- Gray p. 157 (archive L77): *"Sigmund [1969] has developed analytical
  solutions ... in the medium **(0.1-1 KeV)** and high **(>1 KeV)** energy
  ranges"* — 1 keV is the regime boundary.
- Chang thesis L2249-2251: *"In the low ion energy regime **(< 1 kV)**, the
  binary particle interactions can be characterized by a Born-Mayer-type cross
  section, and the sputtering yield is linear to the square root of ion incident
  energy"*.

**Bracket for 3406 eV, ordered by defensibility:**

| | bare | complex |
|---|---|---|
| Steinbruchel sqrt (Gray's own thresholds), extrapolated | 3.74 | 3.37 |
| sqrt to 1 keV, ZBL surface density above | 2.78-2.86 | 2.60-2.68 |
| pure ZBL surface density | 2.33 | 2.54 |

The sqrt value is the **upper** end of the defensible band. Whatever is chosen
inside 2.3-3.7, it is 3-5x below what petch currently applies.

---

## 6. Verdict

**(a) Physical sputter channel of SiO2.**
Doctrine-compliant form: `Y = A (sqrt(E) - sqrt(Eth))`, with
**`A = 0.0139 eV^-1/2`, `Eth ~ 18 eV`** — Gray, MIT PhD thesis (1993) Table 5-1
p. 159, attributed to Chapman 1980, supported by his Fig. 5-2 against Chapman
1980 / Oostra 1986 / Steinbruchel 1989 to 1225 eV; independently carried by
ViennaPS (Rodrigues et al. 2023, Table 1 `As = 0.0139`, `Eths = 18`).
Canonical form: Steinbruchel, *APL* **55**, 1960 (1989), DOI 10.1063/1.102336.
**Transfer factor 350 -> 3406 eV = 3.74** (upper bound; 2.3-2.9 if the ZBL
flattening above 1 keV is honoured).

Threshold caveat, recorded not smoothed: `Eth` for SiO2 physical sputter is
reported anywhere from 18 eV (Gray/Chapman, a sqrt-fit extrapolation) to
70 +/- 5 eV (a 2024 low-pressure-ICP measurement over 20-200 eV) with Chang's
beam fit at ~40 eV, Joubert 1994 (CHF3 plasma) at 35 eV and Thomas III 1990 at
72 +/- 5 eV on native oxide. Krueger's 70 eV is inside the plasma-measured band.
Choosing 70 with the sqrt form gives 4.83 instead of 3.74 — the threshold choice
is a 30% question, the form is a 3x question. **Do not spend a run on Eth.**

**(b) Chemically-enhanced (F-saturated) channel of SiO2.**
Doctrine-compliant form: `beta_e = 0.053 (sqrt(E) - sqrt(4))`, measured over
**20-2000 eV** — Gray thesis Eq. (5-35) p. 252, Table 5-10 p. 247; published as
Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), DOI 10.1116/1.586925,
whose abstract states the law. Fits its own data at **R^2 = 0.993**. The full
saturated yield is `beta_e (1 + 0.007 sqrt(E)) + Y_sputter`.
**Transfer factor 350 -> 3406 eV = 3.37.**

Caveat: for the **CF2-mediated** (as opposed to atomic-F-mediated) channel Gray
measures a *weaker* dependence still, `~E^0.28` over 150-350 eV (Table 6-1
p. 305), because CF2 mobility in the mixing zone does not scale with the mixing
depth. Our chemistry is C4F8, so the complex row is CF2-mediated. **The sqrt
form is therefore an upper bound on the complex channel too**; nothing in the
measured record supports anything steeper.

**(c) Can the Gray magnitudes now be carried to keV?**
**Yes.** The blocker named in `RESULTS_ABSOLUTE_YIELD` — "these are at ONE energy
(350 eV) ... they cannot be applied until the energy scaling is settled" — is
lifted. Gray does not give one energy; he gives six from 20 to 2000 eV, plus the
fitted law, plus the absolute constants, all in the same document. The magnitudes
and the scaling come from the same measurement, so carrying them to 3406 eV is a
single sourced act, not a fit-and-hope. The honest residual is the 3406/2000 =
1.7x extrapolation beyond the top measured point, bounded above by sqrt and below
by ZBL surface density — a 1.4x band, against the 3x error the current form
carries.

**(d) What this predicts, unspent.**
Replacing both rows with Gray's laws:
- bare drops **5.4x** at the front — which is the ion-energy overestimate that
  three prior routes (cascade audit ~1.4x, Krueger's own statement, the joint
  solve's 1.24-1.35x depth gap) each saw a fragment of, now measured whole;
- complex lands within **8%** of its current value at the front, so the
  chemically-enhanced magnitude was never the problem *at keV* — it only looked
  2.8x weak because it was being read at 350 eV under a linear law;
- the bare/complex balance inverts **5.8x**, making the coverage-dependent
  channel dominant at the front, which is the structural precondition for
  neutral-limited behaviour that `RESULTS_LIMITING_REGIME` §2 found unreachable
  by re-weighting.

Forecast before spend: this is a two-line change (`_complex_energy_factor` and
the bare row's `n`) into hooks that already exist, and it should be graded
against C1a (dynamic range), C5 (depth factor) and C7 (ARDE sign) in coupled
mode — the beam-mode forecast failure of `RESULTS_ANGULAR_CONVENTION` §4 applies.

---

## 7. NOT FOUND / open, with routes

- **Kwon Fig. 3.4 p. 76 stated measurement condition.** The MIT DSpace text
  bitstream for Kwon's thesis (`15b8c719-4edf-483b-9d83-11699d97003f`, bundle
  `be0b8937-...`) is unusable OCR garbage. Route: OCR the ORIGINAL bundle the
  same way this doc did Gray's, or read Kwon & Sawin, *JVST A* **24**, 1906
  (2006) DOI 10.1116/1.2336225 and *JVST A* **24**, 1914 (2006) DOI
  10.1116/1.2336226 (the published thesis). Needed to settle whether the 0.28
  floor is at 350 eV (contradicting Gray's own 0.201) or at a higher energy.
- **Steinbruchel 1989 full text** — AIP returns 403. Abstract secured verbatim
  via Crossref; the functional form is independently transcribed in Chang's
  thesis Eq. (3.4). Route: interlibrary or Chang's Table 3.1 (A, Eth for Ar+/Si:
  Chang 0.04/35, Harper 0.04/27, Tachi 0.04/31, Oostra 0.03/64).
- **Mayer et al. 1981** (SiO2 by CFx+, angular peak moving off-normal above
  400 eV) — quoted only through Gray p. 307. Route: JVST/JAP 1981, chase via
  Gray's reference list (thesis pp. 376-384, not OCR'd).
- **The 10-20 keV yield maximum for Ar+ -> SiO2** — from search summaries of
  Seah et al. 2010 and *Radiation Effects* **21** (1974) 5-30 keV data, not read
  directly. Immaterial to the verdict (3406 eV is far below any maximum under
  every source), but marked `[VERIFY]`.
- **Gray's Figure 5-2 digitized points.** The OCR recovers the axes and the
  legend but not the marker coordinates. If the sputter magnitude ever needs to
  be defended point-by-point rather than through Gray's own fitted A/Eth,
  digitize pp. 160 and 162 of `gray_thesis.pdf` at 600 dpi (same route as
  `research_sources/digitized/krueger_fig7*`).
