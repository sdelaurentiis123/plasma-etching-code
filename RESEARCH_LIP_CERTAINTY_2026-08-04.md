# Lip certainty: what actually closes the top-of-mask 15% (2026-08-04)

Scope: settle, before implementing, (1) whether Izawa's 0.004 is importable, (2) whether
You-2023's above-cosine removal applies to FC film on a mask lip at 30–70°, (3) what removes
polymer at the mask **top** specifically, (4) the single highest-certainty zero-knob change.

Convention: **"quoted"** = read verbatim this session from a file I fetched (path given).
*[abstract]* = publisher/OpenAlex abstract, verbatim but abstract-only. `[VERIFY]` = not read.
**[INFERENCE]** = my derivation from things I did read, not a quotable claim.

**DO NOT COMMIT** (per task instruction).

---

## 0. Headline (read this and nothing else)

The 15% is **not a chemistry-constant error and not a missing mechanism**. It is a
**wall-slope** error in the probe's prescribed geometry, and it is quantitatively closed by
**0.8–1.0° of local wall tilt**.

Both petch's removal channels at a mask sidewall carry `cos θ` twice — once in the *areal* ion
flux delivered to the face (`∝ cos θ`) and once in the Kress angular yield
(`kress = (1+9.3 sin²θ)·cos θ`) — while thermal deposition carries **neither**. Writing the wall
tilt from vertical as α (θ = 90°−α), petch's ion removal at a mask wall goes as

```
R_ion(α)  ∝  sin²α · (1 + 9.3 cos²α)          [INFERENCE, derived from src/petch/mixed_layer.py:256-301 + areal-flux convention]
```

which is **0.0053 at α = 1.3° and 0.356 at α = 10.9° — a 67× range over 10 degrees of slope.**

`scripts/mouth_equilibrium_probe.py` builds the aperture as a *Gaussian* bump
(`half_width = 0.5(W − (W−w)·exp(−(Δ/σ)²))`, σ = 0.10 µm, Δ measured from the neck plane).
A Gaussian has **zero slope at its own minimum and zero slope in the far field** — so the
prescribed geometry pins α ≈ 0 both at the neck apex and across the entire top 100 nm, i.e. at
exactly the two places the probe reports as unbalanced. Predicted vs measured wall angle:

| depth below mask top | α from the Gaussian (predicted) | incidence predicted | incidence **measured** (`depth_profile_45nm.json`) | mean net (nm/s) |
|---|---|---|---|---|
| 25 nm | 0.37° | 89.6° | **88.7°** | −0.0469 |
| 75 nm | 2.11° | 87.9° | **87.6°** | −0.0122 |
| 125 nm | 6.72° | 83.3° | **83.9°** | −0.00145 |
| 175 nm | **10.89°** | 79.1° | **80.0°** | **−0.00018** |
| 225 nm | 6.03° | 84.0° | **83.0°** | −0.00063 |

The predicted and measured incidence agree to ≤1° in every band. **The net velocity is a
monotone function of the local wall slope, not of depth and not of aperture.** The band the
probe called "the most stable part of our mask wall" (150–200 nm, −0.00018 nm/s ≈ 99.5% balanced)
is simply the band where the Gaussian happens to pass through α ≈ 11°. The band that "closes
30–50× faster" is the band the Gaussian forces to α ≈ 0.4°.

Two consequences that reframe the whole result:

1. **`equilibrium_aperture_nm: None` is a property of the shape family, not of petch's
   chemistry.** Any *smooth* minimum has α = 0 at the minimum; at α = 0 both removal factors
   vanish and deposition does not; therefore *no* aperture of a smooth neck can balance **at its
   own apex**. The probe searched the wrong axis.
2. **The equilibrium of this mechanism is an equilibrium in ANGLE, not in aperture.** With
   the published Table-I fluxes and Appendix-B probabilities the balance point sits at
   **α ≈ 9–11° from vertical** — which is precisely where the probe already measures balance.
   The "85.5% removal / 15% deficit" at the 39 nm neck corresponds to the neck band's
   flux-weighted α ≈ 8°; getting from 0.196 → 0.239 on the `sin²α(1+9.3cos²α)` curve is
   **α: 8.0° → 8.85°.**

Verdict, stated plainly: **do not import Izawa, do not change the grazing removal law, and do
not add a mechanism, on this evidence.** §4 gives the ordered change list and gates.

---

## 1. Izawa et al., *JJAP* **46**, 7870 (2007) — **NOT importable**

**Citation.** M. Izawa, N. Negishi, K. Yokogawa, Y. Momonoi (Hitachi), "Investigation of Bowing
Reduction in SiO₂ Etching Taking into Account Radical Sticking in a Hole," *Jpn. J. Appl. Phys.*
**46**(12), 7870–7874 (2007), DOI 10.1143/JJAP.46.7870.
**Full text NOT obtained** — IOPscience serves the abstract HTML but returns HTML (not PDF) for
`/pdf`; OpenAlex/Unpaywall: `oa_status: closed`, no repository copy; ResearchGate 403; J-STAGE
404. Saved abstract HTML: `…/scratchpad/lip/izawa_iop.html`.

**Abstract, verbatim** (read from `izawa_iop.html` this session):

> "The bowing mechanism in high-aspect-ratio contact hole (HARC) etching was investigated by
> taking into account reactive sticking on the sidewall of the hole. **Sticking coefficients of
> radicals on the sidewall have been estimated by comparing the observed deposition profile with
> the calculated one.** It was found that the coefficients of **C rich radicals** and **CF x
> radicals** were **0.5** and **0.004**, respectively, and that **F radical reaction probability
> to the fluorocarbon polymer is 0.07**. These coefficient values were deduced that the excessive
> flux of O and F onto the sidewall of a hole causes bowing during HARC etching. It was also
> indicated that the bowing can be suppressed by reducing of the flux of oxygen. These findings
> were confirmed by the results of experiments using an **ultra-high frequency-electron cyclotron
> resonance (UHF-ECR) plasma**."

### 1.1 What is actually measured

Nothing is measured directly. The sentence **"estimated by comparing the observed deposition
profile with the calculated one"** makes 0.004 an **inverse-fit parameter of Izawa's own
in-hole transport model**, extracted from a deposition-thickness-vs-depth profile. Its numerical
value is conditioned on that model's (unstated, unread) assumptions: the re-emission angular law
for non-sticking radicals, whether multiple wall bounces are tracked, how the incident radical
flux was normalised, the assumed hole geometry, and which measured species were lumped into
"CF x". Change any of those and the fitted number moves. It is the same class of object as
Krüger's `pd,poly-AC = 0.094`, which is *also* a fitted number (§3.2) — not a beam-measured
constant.

### 1.2 The 0.5 / 0.004 split is C-rich vs F-rich, not "C-rich" vs "the CFx class"

Third-party reading of Izawa by a group that used it, verbatim from
`…/scratchpad/lip/10_35848_1882-0786_ac8d46.txt` (Hiwasa, Kataoka, Sasao, Kuboi, Iino, Kurihara,
Fukumizu, **KIOXIA**, *Appl. Phys. Express* **15**, 106002 (2022), DOI 10.35848/1882-0786/ac8d46;
their ref. 12 **is** Izawa 2007):

> "This can be attributed to the difference in the sticking coefficients between the C-rich CFx
> radicals and F-rich CFx radicals. **The sticking coefficient of F-rich CFx radicals such as CF2
> is reported to be lower than that of C-rich radicals;¹²)** hence, F-rich CFx radicals can reach
> high-AR regions."

> "F-rich CFx radicals are less reactive to nucleophilic attack,³²⁾ **leading to a low sticking
> coefficient** because most molecules are covered with electron-rich F… the unpaired electron on
> C of the C-rich CFx radicals was less hindered, **leading to a high sticking coefficient**."

So Izawa's two numbers are the *ends of an F/C-ratio axis*: **0.004 belongs to CF₂/CF₃-class
(F-rich) radicals and 0.5 to C-rich radicals** (C₂F₃, C₃F₃, C₆F₅…). They are 125× apart and
**there is no single number to import.**

**This inverts petch's deck.** `src/petch/chemistry_deck.py:113` has
`on_polymer = {CF 0.1, CF2 0.1, CF3 0.1, C2F3 0.03}` — the **C-rich** species (C₂F₃, F/C = 1.5)
gets the **lowest** probability, 3.3× below CF₃ (F/C = 3). Izawa/KIOXIA say the ordering should
be the other way round. That is a real, checkable inconsistency between the two provenances; it
is **not** resolved by adopting either number, and it is **not** the mouth bug (§0).

### 1.3 Commensurability with petch's `deposition_probability_on_polymer`

| axis | Krüger Appendix B row (petch's source) | Izawa 0.004 |
|---|---|---|
| what the number multiplies | one incident radical pseudoparticle onto one surface **cell** | incident radical flux in Izawa's transport model |
| surface state | explicit — `CF(s)`, `CF₂(s)`, `CF₃(s)` (fresh polymer), separate rows for `(xs)` crosslinked (0.02) and `AC(s)` mask (0.2) | "the sidewall", state unspecified |
| removal | **separate rows** — `CF(s)+ion→EP` 0.9@20 eV, `CF(s)+O→EP` 0.0423 | not separated; a *net* deposition profile was fitted |
| ion bombardment during measurement | n/a (mechanism, not measurement) | **unknown** — HARC sidewall during etch is under grazing ion + hot-neutral flux |
| species resolution | 4 species, individually rowed | 2 lumped classes |

Verbatim Appendix B rows (`tmp/pdfs/krueger_thesis.txt` L5350-5390, read this session):

```
CF(s)   + CF  ->  CF(s)  + CF(s)     0.1
C2F3(s) + CF  ->  C2F3(s)+ CF(s)     0.03
AC(s)   + CF  ->  AC(s)  + CF(s)     0.2      <- appendix-converged (paper Table 6.5: 0.094)
CF(xs)  + CF  ->  CF(xs) + CF(s)     0.02
CF(s)   + O   ->  EP                 0.0423
CF(s)   + Ar+ ->  EP  + Ar#          0.9  20  0.5  500  1
AC(s)   + Ar+ ->  C   + Ar#          0.001 200 0.4 250  1
AC(s)   + O   ->  CO                 1.00E-05
```

(The trailing `1` is the angular-form column. Appendix B.1 defines **∠=1 = Kress et al., *JVST A*
**17**, 2819 (1999)** and **∠=2 = Chang & Sawin, *JVST A* **15**, 610 (1997)** — verbatim from
`krueger_thesis.txt` L6611-6616. See §2.3.)

**VERDICT — NOT COMMENSURATE. Do not import 0.004.** Reasons, in order of force:
1. It is a **model-inverted effective coefficient**, not a beam/surface measurement; Krüger's
   0.1/0.094 are **gross per-event site probabilities** with removal carried by separate rows.
   Substituting one for the other double-counts (or un-counts) removal.
2. It is not one number: **0.004 is the F-rich end, 0.5 the C-rich end.** Importing 0.004 across
   `{CF, CF2, CF3, C2F3}` would apply the F-rich value to C-rich species that Izawa puts 125×
   higher.
3. Whether Izawa's sidewall was under concurrent ion bombardment is **unstated in the abstract
   and I could not read the body** `[VERIFY — body not obtained]`. If it was, 0.004 is a *net*
   quantity and importing it into a mechanism that *also* has a sputter row is a double subtraction.
4. **Partial importability (the only defensible use):** as a *shape* constraint —
   `p(C-rich) / p(F-rich) ≈ 125` — i.e. a monotone-in-F/C ordering of the deposition
   probabilities. petch currently has the opposite ordering for C₂F₃. That is a legitimate
   held-out ordering gate; the magnitudes are not liftable.

**Corroborating negative on the "mask vs polymer" question** — Omura *et al.* (Toshiba Memory +
Nagoya), *JJAP* **58**, SEEB02 (2019), DOI 10.7567/1347-4065/ab163c, read verbatim from
`…/scratchpad/lip/10_7567_1347-4065_ab163c.txt`:

> "this inference is based on the hypothesis that **the sticking probabilities of fluorocarbon
> radicals to carbon mask and dielectric films are comparable because the surface will be
> activated by ion bombardment and covered with similar fluorocarbon films.**"

i.e. the industrial default is exactly petch's `on_mask ≈ on_polymer` (0.094 vs 0.1). No change
warranted there either.

---

## 2. You et al., *Coatings* **13**, 1452 (2023) — the above-cosine window is **30–70°**, and our lip is at **80–89°**

### 2.1 Exact conditions (read from `tmp/pdfs/coatings2023_bowing_narrowing.txt`)

| quantity | value | verbatim source line |
|---|---|---|
| reactor | ICP, 13.56 MHz source + 13.56 MHz bias | §2.1 |
| chemistry | **HFE-347mcc3 / O₂ / Ar**, total 30 sccm, O₂ fixed 2 sccm, HFE 8→12 sccm | "The flow rate of O2 was fixed at 2 sccm while the flow rate of HFE-347mcc3 was varied from 8 to 12 sccm." |
| pressure / bias | **1.33 Pa**, source **250 W**, **DC bias −1200 V** | "The source power and the direct current (DC) bias voltage were 250 W and −1200 V" |
| electrode T | 15 °C | ibid. |
| angle control | **Faraday cage**, closed stainless box, grid pitch 0.229 mm, sample tilted inside | §2.2 |
| **NDR surface** | **fluorocarbon film deposited on SiO₂**, **bias = 0 V** | "The process conditions for fluorocarbon film deposition were the same as those for SiO2 contact hole etching, **except that no DC bias voltage was applied**… (source power = 250 W, DC bias voltage = 0 V…)" |
| **NER surface** | **SiO₂ substrate**, bias −1200 V | Fig. 6 caption: "Change in (a) the etch rate and (b) the normalized etch rate of **SiO2** with the ion-incident angle" |

### 2.2 Does the above-cosine result apply to FC film on a mask sidewall under keV Ar⁺? **NO.**

Three verbatim facts, each disqualifying:

1. **The measured quantity is SiO₂ etch rate, not FC-film sputter rate.**
   > "In all conditions, **the NERs are above the cosine curve until the ion-incident angle reaches
   > 50–60 degrees.** This indicates that physical sputtering plays an important role during etching."

   Nothing in You 2023 measures the removal of a fluorocarbon film. It measures oxide removal
   *through* one.

2. **Above 80° the measurement says NET DEPOSITION, not enhanced removal.**
   > "When the ion-incident angles were greater than 80°, **the etch rates exhibited negative
   > values, indicating a net deposition at these angles.** A net fluorocarbon-film deposition
   > instead of substrate etching occurred at ion-incident angles greater than 80° **because the
   > flux of ions at high incident angles is negligible.**"

   Our top lip sits at **86–89°** and our neck band at **79–84°** (§0 table). You-2023's own data
   at those angles says the surface **deposits**. Importing "above-cosine removal" from the
   30–70° window into an 80–89° lip is an extrapolation *against* the measurement's own trend.

3. **The Faraday-cage coupon is not a deep-feature lip.** A tilted flat coupon in an open cage
   sees the **full 2π thermal radical flux** at every tilt angle; a mask-wall face at depth sees a
   small view factor. The NDR/NER pair therefore describes the balance at unshadowed radical
   supply, which is the *worst* case for deposition. **[INFERENCE]**

### 2.3 Other measured angular-yield data — the corpus, and what it bounds

Three independent Faraday-cage / V-groove measurements *in fluorocarbon chemistry*, abstracts
verbatim from OpenAlex:

- **Cho, Hwang, Lee, Moon, *JVST A* 18, 2791 (2000)**, DOI 10.1116/1.1318193, CF₄, 5 mTorr, bias
  −100 to −800 V, Faraday cage: *[abstract]* "The normalized etch-yield curves showed virtually
  the same angular dependence regardless of the ion incident energy. The curve shape was similar
  to that of physical sputtering **except that the ratio of the maximum yield to that at 0° was as
  low as about 1.3.** … partly attributed to the fluorocarbon polymer film, which existed **as a
  few monolayers-thick film on the substrate surface at low angles near 0° but as a submonolayer
  at high angles between 45° and 75°.**"
- **Lee, Hwang, Min, Moon, *JVST A* 20, 1808 (2002)**, DOI 10.1116/1.1503786, CHF₃, −20 to −600 V:
  *[abstract]* "**When the absolute value of the bias voltage was smaller than 200 V, the
  normalized etch rate … changed following a cosine curve** … When the magnitude of the bias
  voltage was larger than 200 V, **the NER was deviated to higher values from those given by a
  cosine curve at ion angles between 30° and 70°, and then drastically decreased at angles higher
  than 70° until a net deposition was observed at angles near 90°.**"
- **Schaepkens, Oehrlein, Hedlund, Jonsson, Blom, *JVST A* 16, 3281 (1998)**, DOI 10.1116/1.581534,
  V-groove, high-density FC plasma: *[abstract]* "The SiO₂ etch rate on 54.7° inclined surfaces is
  lower than on flat surfaces, while **the SiO₂ etch yield (atoms/ion) is a factor of 1.33
  higher.** … **The fluorocarbon deposition is decreased at 54.7° whereas the fluorocarbon etching
  rate is increased at 54.7°.** This produces a thinner steady-state fluorocarbon film on the
  inclined … surface."
- **Chae, Vitale, Sawin, *JVST A* 21, 381 (2003)**, DOI 10.1116/1.1539085, QCM, direct yields vs
  energy / ion-to-neutral ratio / **angle**: *[abstract]* "**Angular yield measurement shows that
  when fluorocarbon deposition is relatively severe, etching yield decreases significantly as the
  incident angle increases and deposit fluorocarbon at a high incident angle above 60°.**" Also:
  "Two fluorocarbon deposition mechanisms are identified: neutral deposition and **ion-enhanced
  deposition**. The low-energy ions … enhance the deposition rates by creating active sites."

**Negative result, stated explicitly so nobody re-searches it: I found no published measurement
of the angular sputter yield of a fluorocarbon FILM by itself.** Every angle-resolved dataset in
this literature is a *substrate* (SiO₂ / Si₃N₄ / photoresist) etch rate measured through an FC
film. Schaepkens is the only source that states the FC-film removal rate is angle-enhanced, and
it gives no number.

### 2.4 Where petch's B = 9.3 actually comes from, and how it compares

`kress_1999`, `f(θ) = (1 + B(1−cos²θ))·cos θ`, `B = 9.3`
(`src/petch/surface_kinetics.py:65-67`, `mixed_layer.py:256-257,297-298`,
`chemistry_deck.py:166`). In petch this factor is applied to the **polymer sputter channel only** —
`kernel_complex`, `kernel_bare` and `kernel_ac` carry **no** angular factor
(`mixed_layer.py:270-310`).

- **Provenance of the *form*: solid.** Krüger Appendix B marks the polymer and AC sputter rows
  ∠=1, and B.1 defines ∠=1 as **Kress et al., *JVST A* 17, 2819 (1999)**. Huang's thesis (the
  parent MCFPM) describes it verbatim (`tmp/pdfs/huang_thesis.txt` L2290-2296):
  > "For physical sputtering, f(θ) is an empirical function with a **maximum at 60°, reduced
  > probability at normal incidence and zero probability at grazing incidence.** For chemically
  > enhanced etching, f(θ) is unity for normal incidence and angles up to 45°, with a monotonic
  > roll-off to zero probability at grazing incidence."

  petch's B=9.3 curve peaks at **52.6°** and is zero at 90° — the right shape.
- **Provenance of the *value 9.3*: none found, and the source paper is off-domain.** Kress,
  Hanson, Voter, C. L. Liu, X.-Y. Liu, Coronell, *JVST A* **17**, 2819 (1999) is titled
  *[abstract, OpenAlex]* **"Molecular dynamics simulation of Cu and Ar ion sputtering of Cu(111)
  surfaces"**, for "ionized physical vapor deposition, used in **Cu interconnect technology**",
  over "**10–100 eV for Cu ions and 50–250 eV for Ar ions**". It is an MD study of a **metal**
  at **≤250 eV**. petch applies its angular parameter to a **fluorocarbon polymer** at
  **~1.5 keV**. `[VERIFY]` — I did not obtain the Kress body and found no source printing 9.3.
- **Magnitude check against the in-chemistry measurements.** petch's `kress` is the yield per
  incident ion (petch's `atom_flux` is the *areal* flux onto the face and already carries cos θ —
  confirmed numerically: the probe's flat mask-top face reads 9.68e19 m⁻²s⁻¹ against a Table-I
  incident 1.2e20, and an 88.7° face reads 1.25e18 ≈ cos(88.7°)×9.68e19). So petch's
  yield-enhancement is **max/normal = 4.17 at 52.6°** and **3.99 at 60°**, against
  **≈1.3 measured** (Cho 2000, CF₄ plasma, all bias voltages) and **1.33 measured** (Schaepkens
  1998, 54.7° V-groove). Those are SiO₂ yields, not FC-film yields, so this is a **bound, not a
  refutation** — but it is the only in-chemistry bound that exists, and petch sits **~3× above it.**

**VERDICT on Q2.** (a) The You-2023 above-cosine result is **SiO₂ etch rate at 30–60°, bias
−1200 V, on an unshadowed Faraday-cage coupon** — it does **not** license extra grazing removal
of FC film on a mask lip at 80–89°; at those angles You's own data reads *net deposition*.
(b) There is **no measured FC-film angular sputter yield** in the literature to gate the lip law
with. (c) petch's B = 9.3 is already **above** every in-chemistry angular measurement by ~3×, so
the direction of any evidence-driven change to the lip removal law is **down, not up** — which
would make the mouth worse, not better. This is the strongest reason to reject "add 15% more
grazing removal" as the fix.

---

## 3. Top-of-mask: what actually removes polymer there

### 3.1 The ion channel is *structurally* dead at the top lip — and that is correct physics

`depth_profile_45nm.json`, per-face: the 0–50 nm band has **median incidence 88.7°** (the
"3.14e19 mean ion flux" quoted in the probe report is the *mean*, dragged up by a handful of
flat mask-top faces at 1.1° incidence; the **median face flux in that band is 1.25e18**, 1.3% of
the mask-top value). At 88.7°, `kress = 0.234` **and** the areal flux is 2.3% of incident. The
product `cos²θ(1+B sin²θ)` is **0.0053** — a **190× collapse** vs a flat top surface.

So the lip is not "the highest ion flux region"; on a near-vertical wall it is one of the
lowest. Nothing needs fixing about `kress(1) = 1` — the flat mask top is fine (its net velocity
is −4.7e−9 nm/s, i.e. exactly balanced). The near-vertical wall immediately below it is where
the ion channel disappears.

### 3.2 Krüger names the top-of-mask remover explicitly, and it is **O radicals**

Verbatim, `tmp/pdfs/krueger-2024.txt` (Krüger, Huang, Kushner *et al.*, *JVST A* **42**, 043008
(2024), DOI 10.1116/6.0003554), the paragraph immediately preceding Fig. 6:

> "Since the control of necking and clogging by polymer deposition is of special interest, we
> discuss an important control parameter. The amount of polymer growth is determined by the
> relative contributions of deposition and removal. **A steady state polymer thickness occurs when
> these contributions balance. Polymer removal occurs by sputtering and O-radical based etching.
> In oxygen rich gas mixtures, polymer removal can be dominated by O-radical based etching, mostly
> ground state atomic oxygen.** In the mechanism discussed here, **necking and clogging (the amount
> of polymer deposition in the mask region) can ultimately be controlled by the reaction
> probability of the O based polymer etching.** … **An etch probability of 0.005 results in net
> polymer growth ultimately leads to a complete clog at the top of the feature.** Higher removal
> probabilities, while still producing significant necking, do not fully clog the feature…"

and (`krueger-2024.txt` L891, and thesis L4740):

> "The decreased **polymer deposition at the top of the mask** is a consequence of the increasing
> **polymer etch by oxygen**."

Fig. 6 caption, verbatim: "Etch features for different O based polymer etch probabilities of
(a) **0.005** and (b) **0.02**." Converged value, Table 6.5, verbatim: **`pe,poly = 0.0423`**.

**This is the answer to "what removes polymer at the mask TOP specifically": ground-state atomic
O, at reaction probability 0.0423 per site, non-angular.** Krüger's Fig. 6(a) failure mode —
"complete clog **at the top of the feature**" — is *literally* petch's symptom, and his diagnosis
is insufficient O-based polymer etch. petch already carries 0.0423 (`chemistry_deck.py:134`), so
if the O channel is the problem it is in **O flux delivery or O-term scaling**, not the constant.

Also relevant, from the wider Kushner lineage — the reason the *normal-incidence* removal at a
polymer-covered surface is not suppressed (`…/scratchpad/lip/kushner_jap97.txt`, Kushner group,
*J. Appl. Phys.* **97**, 023307 (2005)):

> "In this work, f(θ) is a semiempirical function, typical of chemically enhanced sputtering with
> a maximum near θ = 60°. This angular dependence is appropriate for the fundamental chemically
> enhanced sputtering which occurs at the surface or interface between, for example, polymer and
> underlying SiO2. **Etch yields which peak at normal incidence results from the penetration of
> activation energy from ions through the overlying polymer.**"

### 3.3 Is the experiment's mask top polymer-covered? — **Yes, and open**

Krüger's Fig. 7 (paper page 8; rendered from `tmp/pdfs/krueger-opening/page-08.png`, crops in
`…/scratchpad/lip/fig7_top.png`, `fig7_sem_top.png`, inspected this session):

- **Fig. 7(a)** (his simulated feature): the polymer band (blue/green) is **continuous from the
  mask top surface downward**. It is *thin* at the very top, thickens to the `wm` minimum roughly
  25–30% down the mask, then thins again. So **his simulation also puts polymer on the mask top —
  and it is thinnest there.** The metric definition is verbatim: "**The width of the mask opening
  wm including deposition** stands as a measure of the necking and clogging."
- **Fig. 7(b)** (the base-case SEM): the a-C mask top is **rounded/dished**, curving smoothly into
  the opening, with a bright thin rim running continuously over the shoulder and down both
  sidewalls. The narrowest aperture is clearly **below** the shoulder. **[INFERENCE from the image
  — SEM contrast cannot certify "polymer" vs "topographic edge", so the *thickness* is not
  readable; the *shape* is.]**
- **The mask does not recede.** Table IV target metrics, verbatim: `hm  Remaining mask thickness
  **850 nm**`, against "A SiO₂ substrate is covered by a **850 nm thick** AC film patterned to
  ideally yield a straight walled opening with an **initial width of 90 nm**." Initial = final.
  Consistent with Appendix B (`AC + ion → C`, p₀ = 0.001, ε_th = 200 eV; `AC + O → CO`, 1e-5).
  So **"the retreating mask carries the polymer away" is refuted** — there is no retreat.

So: the experiment's mask top **is** polymer-covered and **is** open. petch's closure at that
band is not "right but slow" — a 45 nm neck at 200–270 nm is the experiment's answer, and the top
100 nm stays wide.

### 3.4 The per-site budget at the top lip, from published numbers only

Using Krüger Table-I base-case boundary fluxes
(`data/experimental/krueger_2024/base_case_boundary_fluxes.csv`, cm⁻²s⁻¹: CF 4.4e16, CF₂ 9.4e16,
CF₃ 8.4e15, C₂F₃ 6.8e16, O 7.7e16, ions 1.2e16) with Appendix-B probabilities, **per site**:

```
deposition  = 0.1(4.4e16 + 9.4e16 + 8.4e15) + 0.03(6.8e16) = 1.668e16 sites cm^-2 s^-1
O-etch      = 0.0423 (7.7e16)                              = 3.257e15 sites cm^-2 s^-1
O-etch / deposition                                        = 0.195
```

Both O and CF*x* are thermal, so they share the same view factor at any face and the ratio is
**geometry-free**: **the O channel alone can only ever supply ~20% of what deposition brings.**
The remaining 80% must come from ions — which is exactly why the balance is so violently
sensitive to wall slope (§0) and why the near-vertical band cannot hold. **[INFERENCE, arithmetic
above.]**

Fitting the probe's own 39 nm neck receipt (removal/deposition = 0.855) to this budget gives the
face view factor V ≈ 0.15 and the split **ion 77% / O 23%** of removal at the neck. **[INFERENCE]**

---

## 4. VERDICT — the single highest-certainty change, and its gate

### 4.1 THE CHANGE (rank 1, zero knobs, published provenance, no new physics)

**Stop searching for an equilibrium *aperture* and search for the equilibrium *wall angle*.
Re-run the probe with the wall slope as the swept variable, and with the top-of-mask shoulder
taken from the SEM instead of from a Gaussian tail.**

Concretely:

1. `scripts/mouth_equilibrium_probe.py:114-115` prescribes
   `half_width = 0.5(W − (W−w)·exp(−(Δ/σ)²))`, σ = 0.10 µm. This forces α → 0 at the neck apex
   **and** across the whole top 100 nm. The probe's neck-face selector
   (`neck_face_mask`, ±0.08 µm) then averages over α ∈ [0°, 10.7°], landing at a flux-weighted
   α ≈ 8°. **The reported "15% deposition surplus" is the balance at α = 8°, and the balance
   point of the published mechanism is α ≈ 9–11°.**
2. Replace the sweep variable: hold a fixed depth band, prescribe the **local wall tilt α**, and
   bisect on α for `net = 0`. This is the same cost (one transport gather per geometry) and it
   asks the question the mechanism actually answers.
3. Replace the top-of-mask shoulder with the **digitised Krüger Fig. 7(b) SEM shoulder** (rounded,
   finite slope at depth 0) instead of the Gaussian far-field tail (slope ≈ 0.37° at 25 nm).

**Why this is the highest-certainty item:** every input is already in the repo and already
verified. The predicted-vs-measured incidence table in §0 agrees to ≤1° in all five bands, the
net velocity is monotone in α over a 250× range, and the required correction (**+0.8–1.0° of
wall slope**) is smaller than the difference between the Gaussian tail and the SEM shoulder.
It changes **no constant**, imports **no** literature number, and adds **no** mechanism.

**Gate data (all held-out, all already digitisable):**
- `net(α)` must cross zero at **α = 9–11°** at the neck depth, with the published Table-I fluxes.
- The digitised Krüger Fig. 7(b) shoulder angle in the top 100 nm must be **≥ the equilibrium α**
  (if the SEM shoulder is steeper than the balance angle, the top lip erodes — the experiment's
  answer — and petch reproduces it with no new physics).
- `z_neck` must come out at **200–271 nm** (SEM / his MCFPM) rather than 0–100 nm, and `wm` at
  **45 nm** (Table IV target), once α is the free variable.

### 4.2 Rank 2 — the one open units question, and its receipt (do **not** implement blind)

**[INFERENCE, high leverage, needs a receipt first.]** Krüger's Appendix B is written **per
surface cell / per monomer**: `CF(s)+CF → CF(s)+CF(s)` adds one cell; `CF(s)+Ar⁺ → EP` removes one
cell; `CF(s)+O → EP` removes one cell. petch's mixed layer is written **per atom** on the removal
side but **per monomer** on the deposition side:

```
src/petch/mixed_layer.py:337-350   dep_c = p_dep * Gamma_CFx          # 1 C  ... plus dep_f = dep_c * F/C  -> 2.5 atoms/event
src/petch/mixed_layer.py:357-359   sput_c/f = kernel_sputter*theta*x_c/x_f   # 0.9*f*kress atoms/ion total
src/petch/mixed_layer.py:362-364   ox_c = p_ox * Gamma_O * theta * x_c ; ox_f = ox_c*(F/C)   # 1 atom/event total
```

With F/C = 1.5, x_c = 0.4: **one deposition event adds 2.5 atoms; one sputter or O event removes
1 atom.** Under Krüger's per-cell convention all three are one cell. If the per-cell reading is
right, petch's removal/deposition ratio is **2.5× low** — which would take the 39 nm neck from
0.855 to 2.1 (over-correcting badly, i.e. it cannot be applied naively to both channels).

**Do not change this on argument.** The receipt that settles it, with no feature run:
- Predict petch's **blanket** FC-film growth rate and steady-state thickness at the Table-I
  fluxes with bias off and on, and compare against a measured blanket FC deposition/etch rate at
  known flux (Schaepkens/Standaert/Oehrlein blanket data, or Chae-Vitale-Sawin QCM yields per ion,
  DOI 10.1116/1.1539085 — the only direct FC deposition-and-etch **yield-per-ion** dataset found).
- Independently: check that petch's polymer sputter yield at normal incidence and 1.5 keV
  (`0.9·((1500−20)/480)^0.5 = 1.58` in petch's atoms/ion, or 1.58 cells ≈ 4 atoms/ion under the
  per-cell reading) brackets the measured yield. Only one of the two will.

### 4.3 Rank 3 — small, published, and currently missing

- **`kernel_ac` carries no angular factor** (`mixed_layer.py:309`), although Appendix B marks
  **every** `AC(s)/AC(xs) + ion/# → C + #` row **∠=1 (Kress)**. Zero-knob fix (apply the same
  `atom_kress` already computed two lines above). Effect is small in magnitude (a-C p₀ = 0.001,
  ε_th = 200 eV ⇒ ≲0.6 nm of corner erosion in 60 s **[INFERENCE, arithmetic]**), which is *itself*
  a useful result: it confirms the a-C mask cannot facet during this 60 s etch, so the SEM's
  rounded shoulder is **inherited from the mask-open step**, not created by the oxide etch — which
  is the physical justification for item 4.1.3.
- **Deck ordering check** (§1.2): `C2F3 = 0.03` vs `CF3 = 0.1` on polymer inverts the
  Izawa/KIOXIA F/C ordering. Flag as a declared divergence from the Hitachi/KIOXIA evidence, keep
  Krüger's values (they are the deck's declared provenance), and note it in the deck provenance
  string rather than changing anything.

### 4.4 Explicitly REJECTED on this evidence

| candidate | why rejected |
|---|---|
| Import Izawa 0.004 as the lip deposition probability | not commensurate (§1.3): model-inverted, net-vs-gross ambiguous, and it is the **F-rich** end of a 125× split whose C-rich end is 0.5 |
| Add above-cosine grazing removal per You 2023 | the above-cosine window is **30–70°**; our lip is at **80–89°**, where You/Lee/Cho all measure **net deposition** (§2.2) |
| Raise the Kress B to strengthen the lip | petch is **already ~3× above** the only in-chemistry angular measurements (1.3 / 1.33) — the evidence points **down** (§2.4) |
| Mask-top recession carrying polymer away | refuted by Krüger Table IV: `hm` target **850 nm** = initial thickness (§3.3) |
| Mask faceting during the etch as the mouth mechanism | a-C sputter 0.001 @ 200 eV gives ≲0.6 nm in 60 s (§4.3); also refuted independently in `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md` §3.2 |
| Charging / redeposition | unchanged from the two prior passes: wrong sign or unclaimed for the mouth |

---

## 5. Fetch ledger (this session)

Obtained, converted, and read (`…/scratchpad/lip/`):

| file | what | route |
|---|---|---|
| `izawa_iop.html` | Izawa 2007 **abstract only** (body paywalled) | iopscience HTML, 200 |
| `10_7567_jjap_57_06ja01.pdf/.txt` | Ishikawa *et al.*, *JJAP* **57**, 06JA01 (2018) HAR progress review (OA) | iopscience `/pdf`, 1.88 MB |
| `10_35848_1882-0786_ac8d46.pdf/.txt` | Hiwasa *et al.* (KIOXIA), *APEX* **15**, 106002 (2022) — a-CF*x* vs AR, C-rich/F-rich sticking, cites Izawa | iopscience `/pdf`, 0.92 MB |
| `10_7567_1347-4065_ab163c.pdf/.txt` | Omura *et al.* (Toshiba Memory), *JJAP* **58**, SEEB02 (2019) — striation; mask-vs-dielectric sticking equality | iopscience `/pdf`, 1.40 MB |
| `kushner_jap97.pdf/.txt` | Kushner group, *J. Appl. Phys.* **97**, 023307 (2005) — f(θ) definition + refs 22/24/25 | cpseg.eecs.umich.edu, 0.81 MB |
| `fig7_top.png`, `fig7_sem_top.png` | Krüger Fig. 7(a)/(b) crops from `tmp/pdfs/krueger-opening/page-08.png` | local render |

Abstracts obtained verbatim via OpenAlex (bodies **not** obtained):
Cho 2000 (10.1116/1.1318193), Lee 2002 (10.1116/1.1503786), Schaepkens 1998 (10.1116/1.581534),
Chae-Vitale-Sawin 2003 (10.1116/1.1539085), Kress 1999 (10.1116/1.581948).

Wanted, **not** obtained: Izawa 2007 body (IOP paywall — the single remaining `[VERIFY]` on Q1);
Kress 1999 body (AIP paywall — the provenance of B = 9.3 remains unsourced).
