# What opens and holds the 45 nm mouth? — Krüger's own evidence (2026-08-02)

Research pass on the standing ml13 residual: petch reproduces the Krüger fig-7
depth (852 vs 825 nm) but the mask opening stalls at 20–25 nm against the
experimental 45 nm. This document goes to Krüger's own text, tables, appendix
and figures — plus the Kushner-group MCFPM literature the mechanism came from —
and asks what *he* says sets that 45 nm.

**Not committed. Every quantitative claim below carries a source line/page.**

## Sources (all local unless noted)

| tag | artifact | notes |
|---|---|---|
| **[T]** | `tmp/pdfs/krueger_thesis.txt` (extraction of `krueger_thesis_2024.pdf`, 6614 lines / 477 kB) | F. Krüger, PhD thesis, Univ. of Michigan, 2024. Ch. 6 = the fig-7 study; App. B = converged mechanism |
| **[P]** | `tmp/pdfs/krueger-2024.pdf` (text at `<scratch>/krueger_paper.txt`) | F. Krüger, D. Zhang, P. Luan, M. Park, A. Metz, M. J. Kushner, "Autonomous hybrid optimization of a SiO2 plasma etching mechanism", *J. Vac. Sci. Technol. A* **42**, 043008 (2024), doi:10.1116/6.0003554 (submitted 19 Feb 2024, published 28 Jun 2024; UMich + TEL Technology Center America). Thesis Ch. 6 == this article; article **Fig. 7** == thesis **Fig. 6.7** |
| **[H]** | `tmp/pdfs/huang_thesis.txt` | S. Huang, PhD thesis (Kushner group) — origin of the SiO₂ mechanism |
| **[H19]** | `<scratch>/JVSTA_37_031304_2019.txt` | Huang *et al.*, *JVST A* **37**, 031304 (2019), "Plasma etching of high aspect ratio features in SiO₂ using Ar/C₄F₈/O₂" |
| **[HK02]** | Hoekstra & Kushner, *JVST B* **20**, 1077 (2002) | "…III. Photoresist mask faceting, sidewall deposition, and microtrenching" — abstract only (paywalled; content via search snippet, marked [VERIFY-quote]) |
| **[D]** | `tmp/mouth_profiles/` | **new**: aperture-vs-depth digitized from published Fig. 7(a) MCFPM output and Fig. 7(b) SEM (script + CSVs + crops) |

---

## Executive verdict

**None of the four ranked suspects is what holds 45 nm in Krüger's work.** The
answer his own paper gives is blunter and more useful:

> **w_m = 45 nm is not a prediction. It is one of six *fitted target metrics*, and
> two of his five tuned parameters are precisely the two levers that set it.**

The physics he names as *controlling* the neck is a **deposition-vs-removal
balance of the fluorocarbon film inside the mask**, where removal has two terms —
**O-radical etching of the polymer** (the parameter he calls "ultimately"
controlling) and **ion/hot-neutral physical sputtering of the polymer** (proved by
his own P_lf = 0 kW clog). Faceting, charging and product redeposition are, in
his own words and tables, either absent, explicitly disclaimed, or acting in the
closing direction.

Suspect ranking after this pass:

| suspect | verdict from Krüger's own work |
|---|---|
| (1) mask corner faceting | **REFUTED as his mechanism.** The word "facet" occurs **0 times** in the thesis and **0 times** in the paper. His converged mask is straight-walled with a target of *zero* erosion (h_m = 850 of 850). He explicitly lists **mask taper as an unreproduced feature** of the experiment. New measurement below shows the experimental constriction *re-opens to 90 nm below it*, which no mask facet can do |
| (2) in-feature charging | **REFUTED for this pipeline.** "charg" occurs **0 times** in thesis Ch. 6 and **0 times** in the entire article. MCFPM *has* a charging module (§2.2.2) and he uses it in Ch. 4 — where he states mask charging is **"nominal"** for HAR and the potential maximum sits **0.77–1.2 µm deep**, i.e. it is a *depth* effect, not a mouth effect |
| (3) sputtered-product redeposition | **PRESENT BUT MINOR AND CLOSING.** Exactly three redeposition rows survive into the converged Appendix B, all `+ C → AC(s)` at **p = 0.01**. Redeposition in this literature is invoked to *lower* etch rates, never to open a mouth |
| (4) mouth dynamics in time | **NO EVIDENCE EITHER WAY — he never measured it.** "metrics are only derived from the final feature after a fixed etch time" [P p. 043008-14]. Clogging is treated as terminal and diagnosed only indirectly through depth |

**The one thing his work *does* say is missing, and that petch can act on:** the
polymer **curvature** at the neck, and its consequence for the incidence angle of
ions on the neck — flagged verbatim by Krüger as an acknowledged, un-optimized
gap. Combined with his own base-case flux table, the neck removal budget is
**~4:1 ion-sputter-dominated over O-radical etch**, which makes the mouth a
*transport/geometry* observable, not a chemistry one. That is exactly what the
P0/P1a angular-closure work concluded from the other direction.

---

## 1. Faceting, bevel, taper, corner

### 1.1 Word-level search

```
$ grep -ci facet krueger_thesis.txt   -> 0
$ grep -ci facet krueger_paper.txt    -> 0
$ grep -ci "bevel" krueger_thesis.txt -> 0
$ grep -ci facet huang_thesis.txt     -> 12
```

"taper"/"tapering" occurs 7 times in [T] (lines 3483, 4049, 4203, 4763, 4767,
4771, 5314). **Every occurrence in the fig-7 chapter is a statement that his
model does *not* reproduce the taper.**

### 1.2 The verbatim admission ([T] 4757–4774; [P] p. 043008-13)

> "Overall, the final simulated etch feature at epoch 200 is in good agreement
> with the experimental counterpart in terms of the metrics used in the
> optimization process. There are, however, some key differences between the
> predicted profiles and the experimental SEM. For example, **there are
> differences in the vertical position of the minimum in necking and the taper of
> the mask. The mismatch in these properties is due to the fact that they were not
> part of the optimization metrics.** Assuming that the optimization process
> produced good model parameters, these differences indicate might that a) **the
> mechanism lacks processes that would otherwise determine necking location or
> taper**, b) the non-optimized physics parameters are not accurate and/or c) the
> solution is not unique and a second solution might better capture these
> phenomena."

and, immediately after:

> "For example, **taper of the feature is known to be sensitive to the chemical
> sputtering probability as a function of angle of incidence of energetic
> particles. This physics parameter was included in our mechanism but was not part
> of the optimization process.**"

So: mask taper is a *known miss* in his own run, and he names the responsible
physics as the **angular dependence of the chemical sputtering probability**.

### 1.3 Why his mask cannot facet: the constants

Appendix B [T] 5397 ff.:

```
AC(s) + Ar+  →  C + Ar#      0.001  200  0.4  250  1     (∠=1 = Kress 1999)
AC(s) + O    →  CO           1.00E-05
```

and the target metric [T] 4583:

```
h_m   Remaining mask height   850 nm      (initial thickness 850 nm)
```

i.e. **the optimization target is literally zero net mask erosion.** With
p₀ = 0.001 at E_th = 200 eV and an O-inertness of 1e-5, the a-C mask is
sputter-armoured by construction and no facet can develop within 60 s. This is
the same armour petch's ml14 run applied — and ml14 narrowed the mouth
(19–20 nm), consistent with the mask playing no opening role.

### 1.4 The MCFPM *does* facet when the mask is soft — and it causes bowing, not mouth

Huang, same code family, photoresist mask ([H] 5380–5388; [H19] lines 874–880):

> "the erosion of the PR results in **increasing the area of facet at the top of
> PR** and scattering of ions at the facet produces hot neutrals into the feature
> with broad angular distributions. **For reflections from the facets to be a
> direct source of bowing in the SiO₂, there should be line-of-sight from the
> facet to the top of the SiO₂, which is not the case here.** Reflections from the
> facets here broaden the angular distribution for subsequent sidewall
> collisions. The PR has been eroded by approximately 450 nm by the end of etch …
> yielding a selectivity of SiO₂ over PR of 10.7."

and [H] 5145–5148:

> "**The necking resulted from a balance between polymer removal and deposition
> processes** while the bowing was caused by surface scattering of ions from
> secondary facets. Non-uniform necking was found to cause an imbalance in the ion
> flux to the bottom of the feature, resulting in twisting **irrespective of
> charging**."

That is the cleanest single sentence in this whole pass: in this modelling
lineage, **facets → bowing; deposition/removal balance → necking**.

### 1.5 New measurement — the experimental constriction is a polymer neck, not a facet

`tmp/mouth_profiles/extract_mouth_profiles.py` digitizes aperture-vs-depth from
the two panels of [P] Fig. 7. Vertical calibration: h_m = 850 nm. Horizontal
calibration for the simulated panel: aperture at the mask/oxide interface = the
w_t = 90 nm target (met to <1 %). SEM aperture is a valley-FWHM between the two
bright sidewall lines, so it systematically under-reads a few nm.

| depth into 850 nm mask (nm) | **Fig. 7(a) MCFPM** aperture (nm) | **Fig. 7(b) SEM** aperture (nm) |
|---|---|---|
| 60 | 73.8 | 82.3 |
| 120 | 58.0 | 59.5 |
| 160 | 48.8 | 47.5 |
| **200** | 43.4 | **39.0 (minimum, 0.24 h_m)** |
| 240 | 40.8 | 69.7 |
| **271** | **38.8 (minimum, 0.32 h_m)** | ~83 |
| 340 | 46.1 | 91.3 |
| 400 | 60.4 | 89.5 |
| 600 | 79.8 | 85.3 |
| 800 | 86.8 | 75.1 |

Three consequences.

1. **The narrowest point is ~200–270 nm *inside* the mask, in both.** The
   aperture at the very top of the mask is ~90 nm in the simulation and ≥82 nm
   in the SEM. A mask corner facet would put the *widest* point at the top and
   narrow monotonically; the top region does taper, but —
2. **the experiment re-opens to the full 90 nm by 280–340 nm depth and stays
   there.** Eroded mask cannot regrow. The constriction is therefore *deposited
   material*, i.e. the fluorocarbon neck, in the experiment as well as in the
   model.
3. **The mismatch is neck *shape*, not neck *depth of aperture*.** Both land at
   39–43 nm minimum. But the experimental neck is short (full width at 1.5× the
   minimum ≈ 85 nm axially) while the simulated neck is long and blunt (≈ 240 nm
   axially) and sits ~70 nm deeper. That is exactly the "vertical position of the
   minimum in necking" mismatch he names in §1.2, now quantified.

**Verdict on suspect (1): mask faceting is not the mechanism holding the mouth
at 45 nm — in his model (no faceting at all) or in the experiment (constriction
re-opens below).** Petch's own observation that a *softer* mask yields a wider
mouth (24.8) than an armoured one (19–20) is a real effect but is *not* the
route Krüger took, and pursuing it would be reproducing his known error class,
not his result.

---

## 2. Charging

### 2.1 The module exists ([T] §2.2.2, lines 2421–2458, verbatim)

> "Electrostatic charging of features results from ions or neutrals which
> neutralize upon striking surfaces and deposit their charge at the impact site.
> Charge is then retained on the voxel upon which it is deposited until
> neutralized by an opposing charge or transported via conductive charge
> transport. … The electric potential is obtained by implicitly solving Poisson's
> equation, −∇·ε∇φ = ρ, using finite volume techniques. Computationally, this is
> performed using the Successive-Over-Relaxation algorithm employing a parallel
> red-black technique. Each material included in the simulation is assigned a
> dielectric constant and mobilities for positive and negative charge transport …
> Compared to the other aspects of the simulation, directly solving Poisson's
> equation is extremely computationally expensive."

Reference given: Wang *et al.* [18]. Boundary conditions: reflective in x–y,
zero-gradient on top, grounded bottom.

### 2.2 It is not used, or at least never mentioned, for fig-7

```
$ awk 'NR>=4142 && NR<=5180' krueger_thesis.txt | grep -ci charg   -> 0   (all of Ch. 6)
$ grep -ci charg krueger_paper.txt                                 -> 0   (whole article)
$ grep -ci "poisson\|electrostatic" krueger_paper.txt              -> only the reactor-scale FKPM
```

The article's MCFPM description ([P] p. 043008-5) mentions only material
identity, reaction probabilities, the Eq. (2) energy law and the fact that ions
neutralize on impact. There is **no feature-scale field, potential or Poisson
solve anywhere in the fig-7 methods**. [VERIFY] whether charging was silently
enabled — the text does not say either way; the absence is uniform across thesis
chapter and article.

### 2.3 And where he *does* use it, he says the mask/mouth is not where it acts

Ch. 4 ([T] 3320–3326, verbatim):

> "When using the VWT generated and thermal EEADs, **the location of maximum
> potential is not at the bottom of the feature but at depths of 1.2 µm for the
> VWT EEAD and 0.77 µm for the thermal EEAD. The charging of the mask is nominal
> compared to the energies of ions and EFR accelerated electrons. The mask
> charging may be more important for low aspect ratio (AR) features where the
> voltages and powers are lower. However, for HAR features, this is typically not
> an issue.**"

and the magnitude of the effect he *does* attribute to charging ([T] 3292–3296):

> "The maximum potential with the VWT produced EEAD is 227 V and the maximum
> potential produced with thermal electrons is 346 V. **The larger positive
> potential with thermal electrons reduces the final etch depth by 30 %** by
> slowing the positive ions incident into the feature."

**Verdict on suspect (2): charging is a depth/twisting mechanism in Krüger's own
framing, explicitly *not* a mask-region mechanism, and it is absent from the
fig-7 pipeline entirely.** Our depth already matches without it; adding it would
be expected to *reduce* depth by up to ~30 % while leaving the mouth alone.

---

## 3. Redeposition of sputtered products

### 3.1 What the high-level table promises ([T] 4385, 4447; [P] p. 043008-5, Table II)

> "SiO₂ can be removed through physical sputtering by energetic ions and hot
> neutrals. **The sputtered products can be redeposited on other surfaces.**"

Table II row (article wording): `S(s) + SiO₂(g) → S(s) + SiO₂(s)   SiO₂ redeposition`
(thesis Table 6.2 prints the same row as `SiO2(g) + S(s) → SiO2(s) + S(s)`).

### 3.2 What the converged mechanism actually contains

Appendix B is "the surface reaction mechanism used in this work **after
convergence**" ([T] 5330). Scanning all 1229 reaction rows for a gas-phase etch
product appearing as a *reactant*, exactly three survive ([T] 5900–5902):

```
AC(s)   + C  →  AC(s)   + AC(s)     0.01
SiO2(s) + C  →  SiO2(s) + AC(s)     0.01
CF2(s)  + C  →  CF2(s)  + AC(s)     0.01
```

There are **no** `+ SiO2 →`, `+ SiF →`, `+ SiF2 →`, `+ SiF3 →`, `+ CO →`,
`+ CO2 →` or `+ CF4 →` rows. The named sputter/etch products are
`SiO2`, `SiF`, `SiF2`, `SiF3`, `CO`, `CO2`, `CF4`, plus the catch-all `EP`
which the appendix defines as "**a generic etch product that is assumed to be
inert and not tracked further**" ([T] 5340). [VERIFY] whether MCFPM applies a
default sticking to named-but-unrowed products, or whether the Table II SiO₂
redeposition row was dropped during convergence; the appendix header warns only
that *surface neutralization* reactions are omitted from the listing.

**So the only surviving redeposition channel is sputtered mask carbon
redepositing as amorphous carbon at p = 0.01** — one part in a hundred, onto a
material (AC) that is then O-inert (1e-5) and sputter-armoured (0.001@200 eV).
Note the sibling carbonization channel, which produces the same armour *without*
transport ([T] 5464): `CF(s) + Ar+ → AC(s) + F + Ar#  0.01 20 0.5 500 1`.

### 3.3 Direction of the effect in this literature

[H] 1040–1043:

> "These large number of collision may result in **redeposition of the etch
> products**. Conductance limits result in decreased neutral etchant delivery from
> the top to the bottom of the feature, and **redeposition of etch products in
> transport from the bottom to the top of the feature, both of which can produce
> lower etch rates.**"

**Verdict on suspect (3): redeposition is present, tiny, and acts to add material
to the upper feature — a closing, not an opening, mechanism.** It cannot be what
holds 45 nm. (It *is* a plausible contributor to why the experimental neck sits
higher and sharper than his simulated one — [VERIFY], speculative.)

---

## 4. His own mouth dynamics in time

**He has none.** [P] p. 043008-14 / [T] 4778–4788, verbatim:

> "In the event of clogging, the etch process is stopped which prevents evaluation
> and optimization of the etch process as a whole. **Since in this work metrics
> are only derived from the final feature** [article: "after a fixed etch time"],
> **there is no mechanism to differentiate between two fully clogged features
> based on the mask metrics alone. The width of mask opening is w_m = 0 in both
> cases even if the rate of the deposition were different and the clogging
> occurred at different times during the etch.** However, differences in rates of,
> for example, deposition can affect the time during which the feature is not
> clogged, which translates to differences in total etch depth. **Through this
> mechanism, the etch depth acts as a secondary metric for the clogging
> mechanism**, without which prior attempts to optimize the isolated mask
> mechanism failed."

And the caption disclaimer on the only feature *sequence* he publishes ([T]
4741–4744):

> "(These features are predictions at the end of the full etch period for a
> particular set of model parameters corresponding to an epoch. **The sequence of
> features is not the temporal evolution of the feature itself**.)"

**Verdict on suspect (4): monotone-narrow vs open-then-close vs close-then-open
is unanswerable from Krüger.** He measures a single 60 s endpoint. His only
dynamical statement is that clogging is *terminal and irreversible* — the
feature never re-opens once sealed, which is why he needed etch depth as a
proxy. This matters directly for petch: our campaign-1 result (mouth passes
through 45 at ~6 s and seals at ~33 s) is not comparable to anything he
published, and his 45 nm is a 60 s endpoint only.

---

## 5. What his own work says *does* set 45 nm

### 5.1 It is a fitted target, and the fit had exactly the right two knobs

Table 6.4 / Table IV, target metrics ([T] 4578–4584):

| symbol | description | target |
|---|---|---|
| w_m | Width of mask opening | **45 nm** |
| w_t | Width at the top of the feature | 90 nm |
| w_f | Maximum width of the feature | 90 nm |
| h_f | Etch depth | 825 nm |
| h_m | Remaining mask height | 850 nm |
| a_h | Asymmetry | 0 |

Table 6.3 / Table III, tuning parameters ([T] 4571–4577): `ps,SiO2`,
`ps,SiO2CFXY`, `pp,SiO2`, **`pe,poly` (O based polymer etch)**, **`pd,poly-AC`
(polymer deposition probability on mask)**.

Converged values ([T] 4729–4735 and [P] Table V):

| parameter | single feature (the fig-7 run) | 4-feature array |
|---|---|---|
| ps,SiO2 | 0.0852 | 0.0909 |
| ps,SiO2CFXY | 0.1471 | 0.1384 |
| pp,SiO2 | 0.278 | 0.2729 |
| **pe,poly** | **0.0423** | 0.0628 |
| **pd,poly-AC** | **0.094** | 0.0842 |

Gradient descent + Nelder-Mead, 200 epochs, ~20 CPU-hours per optimization,
final loss "on the order of a single cell or less" at 1 nm voxels ([T] 4726;
[P] p. 043008-14).

He describes the resulting w_m as a *consequence of the fit*, not a prediction
([T] 4747–4752):

> "During epochs 0 – 20, belonging to stage I (initial descent), the monotonic
> behavior is reflected in the increasing etch depth and **widening of the
> necking**. … **The decreased polymer deposition at the top of the mask is a
> consequence of the increasing polymer etch by oxygen.**"

### 5.2 The mechanism he names: deposition/removal balance, O-etch as the master knob

[T] 4415–4433, verbatim:

> "Since the control of necking and clogging by polymer deposition is of special
> interest, we discuss it as an important control parameter. **The amount of
> polymer growth is determined by the ratio of deposition and removal. A steady
> state polymer thickness occurs when these contributions balance.** Polymer
> removal occurs by sputtering and O-radical based etching. In oxygen rich gas
> mixtures, polymer removal can be dominated by O-radical based etching, mostly
> ground state atomic oxygen. …
> **In the mechanism discussed here, necking and clogging (the amount of polymer
> deposition in the mask region) can ultimately be controlled by the reaction
> probability of the O based polymer etch.** The final etch profiles for otherwise
> identical process conditions are shown in Figure 6.6 while varying the
> probability of polymer etching by O-atoms. **An etch probability of 0.5 %
> results in net polymer growth which ultimately leads to a complete clog at the
> top of the feature. Higher removal probabilities, while still producing
> significant necking, does not fully clog the feature** and allows for continued
> etching throughout the entire process."

(0.5 % vs 2.0 % in Fig. 6.6 / [P] Fig. 6; converged value 4.23 %.)

### 5.3 The second removal term: ion energy, proven by his own transfer case

[T] 4879–4882, verbatim:

> "**The profiles, experiment and simulation, for P_lf = 0 kW produce total
> clogging of the mask opening, indicating that ion energy plays an important role
> in removing excess polymer.** The remaining cases (P_lf = 4, 6 and 8 kW) have
> unclogged features and full etching with unexpected little variation as a
> function of LF power."

and — note this for the petch scorecard's power-ratio misses —

> "These trends indicate that above a certain threshold energy the etch
> progression and the mask removal process are not ion starved, but rather limited
> by neutral gas transport. … These outcomes indicate that **the effect of ion
> energy (for example in sputter yield or related processes) might be
> overestimated in the mechanism.**"

### 5.4 How big is the ion term at the neck? (derived, not quoted)

From base-case Table I ([P] p. 043008-4):

| species | flux (cm⁻² s⁻¹) |
|---|---|
| C₃F₄ | 9.5 × 10¹⁶ |
| C₂F₃ | 6.8 × 10¹⁶ |
| CF | 4.4 × 10¹⁶ |
| CF₂ | 9.4 × 10¹⁶ |
| CF₃ | 8.4 × 10¹⁵ |
| O | 7.7 × 10¹⁶ |
| Ions | 1.2 × 10¹⁶ |

**Derived (mine, not his):** on fresh polymer, deposition
= 0.1·(CF+CF₂+CF₃) + 0.03·C₂F₃ ≈ **1.67 × 10¹⁶** site-events cm⁻² s⁻¹.
O-radical removal = 0.0423 · 7.7e16 ≈ **3.3 × 10¹⁵**. Ion/hot-neutral polymer
sputter has p₀ = 0.9 at E_th = 20 eV, E_r = 500 eV, q = 0.5 ([T] 5463), so at the
several-keV energies of Fig. 6.4 the energy factor saturates the probability;
the ceiling on ion removal is the ion flux itself, **1.2 × 10¹⁶ × f_Kress(θ)**.

⇒ **At normal-ish incidence the removal budget at the neck is ~78 % ions /
~22 % O radicals**, and deposition and removal are within ~10 % of each other.
Two consequences:
- the neck equilibrium is *acutely* sensitive to how much ion flux reaches the
  neck and at what angle — a transport observable;
- an engine that under-delivers ion flux (or mis-orients the local normal at the
  neck) closes the mouth even with perfect chemistry constants. This is the
  quantitative bridge between the mouth residual and the P0/P1a √2 azimuthal
  closure finding.

### 5.5 The gap he himself flags: polymer curvature → ion incidence angle

[T] 4470–4478, verbatim (the paragraph that defines w_m):

> "The width of the mask opening w_m **including deposition** stands as a measure
> of the necking and clogging, a process which is dependent on the ration [sic] of
> the fluxes of polymer depositing to removing species. The narrowed opening can
> impede neutral gas transport into the feature, trap etch products inside the
> feature and shadow the trajectories ions and photons. **Although the shape
> (curvature) of the polymer deposition is not a shape parameter in this
> investigation, the curvature of the polymer affects the angle of the trajectory
> with which ions (hot neutrals) reflect from its surface deeper into the surface.
> An improved approach might include shape (curvature) in the optimization
> process.**"

and, on the mask ([T] 4486–4490):

> "**The thickness of the mask determines the degree to which ions having broad
> angular spread are shadowed by and reflect off the mask prior to entering the
> feature.**"

Together with §5.1's admission that the *angular* dependence of chemical
sputtering "was included in our mechanism but was not part of the optimization
process", this is Krüger's own shortlist for the neck: **local surface angle at
the neck × angular yield × ion delivery.**

### 5.6 Crosslinking gives the neck its anisotropy

[T] 2459–2470 and 4405–4412; Appendix B rows [T] 5695/5697:

```
CF(xs) + Ar+ → EP     + Ar#   0.6  50  0.5  500  1   (crosslinked sputter — harder)
CF(xs) + Ar+ → CF(s)  + Ar#   0.3   8  0.5  500  1   (ion de-crosslinking)
CF(s)  + Ar+ → EP     + Ar#   0.9  20  0.5  500  1   (fresh sputter — softer)
```
Deposition on crosslinked film is 0.02 vs 0.1 fresh ([T] 5376 vs 5368-neighbours).
Krüger: "Since these bonds can be broken by exposure to highly energetic particles
… **this spatially discriminate activation can result in anisotropic shapes if
deposits in the neck area of the feature**" ([T] 2464–2466).

---

## 6. Cross-checks against petch's current configuration

These fell out of the pass and are directly actionable.

### 6.1 Grid resolution — the single largest un-listed suspect

[P] p. 043008-5, verbatim: "**MCFPM resolves the bulk and surface properties on
3D cubic mesh. Here, the mesh cells—voxels—have equal side lengths of 1 nm.**"
(Thesis Ch. 4 and Ch. 5 use 5 nm for their 100 nm-CD cases: [T] 3103, 3665.)

The preregistered ml16 pilot runs `--dx-um 0.01` = **10 nm**. His 45 nm neck is
**45 voxels** wide with **22 nm ≈ 22 voxels** of film per side; ours is
**4.5 cells** wide with 2–3 cells of film per side. One cell of film per side is
a 20 nm swing in the aperture — **larger than the entire 45→24.8 discrepancy.**
His own convergence criterion was that the residual be "on the order of a single
cell or less" *at 1 nm*.

**Recommended gate (cheap, decisive): rerun ml13 at dx = 2.5 nm and 5 nm on the
mask band only (or globally on a shortened etch) and plot w_m(dx).** If w_m rises
monotonically toward 45 as dx falls, the mouth residual is discretization, and no
mechanism work is warranted. This is a pure numerics ablation with no new
physics, and it is the one control Krüger's paper makes explicit.

### 6.2 Mixed constant columns in `src/petch/mixed_layer.py`

`MixedLayerParams` defaults (`src/petch/mixed_layer.py:51,53,58`):

```python
sticking_probability: float = 0.0842   # 4-feature column pd,poly-AC
oxidation_probability: float = 0.0423  # SINGLE-feature column pe,poly
film_sputter_yield: float = 0.1384     # = 4-feature ps,SiO2CFXY (complex sputter)
```

The neck is set by the *ratio* deposition/removal. Self-consistent published
pairs are **0.094 / 0.0423 = 2.22** (single feature — the column that produced
Fig. 7) and **0.0842 / 0.0628 = 1.34** (4-feature array). The dataclass default
is **0.0842 / 0.0423 = 1.99**, i.e. ~48 % more deposition-dominant than the
4-feature pair it half-belongs to. `amorphous_carbon_mask.py` (the router the
pilot actually calls, `scripts/krueger_2024_trench_pilot.py:781`) *is*
self-consistently 4-feature (0.0842/0.0628), so the pilot may be unaffected —
but the dataclass default is a live footgun and the mixing is not receipted as
such. **Check which pair the ml13 run of record actually applied, and declare
one column.** (Both columns target the *same* w_m = 45 / h_f = 825 SEM, so
either is defensible; mixing is not.)

### 6.3 Appendix B contradicts Table 6.5 on the mask deposition constant

Appendix B row [T] 5368: `AC(s) + CF → AC(s) + CF(s)   0.2` (and CF₂, CF₃, C₂F₃
identically 0.2, plus the same on AC(xs)). Table 6.5 / Table V converged
`pd,poly-AC` = **0.094**. Appendix B reproduces the single-feature column exactly
for the other four tuned parameters (0.0852, 0.1471, 0.278, 0.0423 — see [T]
5388, 5903-neighbours, 6006) — so the AC rows appear to carry the
**pre-optimization default of 0.2 and were not updated**, a 2.1× over-deposition
on the mask relative to the converged value. This is very likely the reason ml15
("COMPLETE verbatim converged set") collapsed to opening 7.2 while ml13 (paper
Table V) held 24.8. **The literal Appendix-B set should not be treated as the
fig-7 set.**

### 6.4 Already correct in petch, recorded for completeness

- **C₃F₄ is inert.** It is the *largest* neutral flux in Table I (9.5e16) and
  appears **zero times** in the 1229-row Appendix B. petch already declares it
  inert (`src/petch/amorphous_carbon_mask.py:295`,
  `src/petch/surface_kinetics.py:899`). Good — treating it as a depositing
  radical would have inflated the depositing flux by ~44 %.
- Kress/Chang-Sawin attribution confirmed verbatim: Appendix B "∠=1
  corresponding to the results obtained by [1] and ∠=2 … by [2]" with
  [1] = Kress *et al.*, *JVST A* **17**, 2819 (1999) and [2] = Chang & Sawin,
  *JVST A* **15**, 610 (1997) ([T] 5334–5336 and B.1 references at end of file).
  Polymer and AC sputter rows are ∠=1 (Kress); SiO₂-complex sputter rows are ∠=2.

### 6.5 Declared omission that bears on the neck

`amorphous_carbon_mask.py` known_omissions includes "**crosslinked-film sputter
resistance is not represented**". Krüger's neck sits exactly where the ion dose
crosslinks the film, and his crosslinked polymer is *harder* (0.6@50 eV) than
fresh (0.9@20 eV) while also accepting 5× less deposition (0.02 vs 0.1). Absent
the resistance term, our neck film is uniformly the *soft* species — which
should, if anything, open the mouth, so this omission does not explain the
residual and is listed as a direction-checked non-suspect.

---

## 7. Broader literature — a-C / mask faceting under energetic ion bombardment

Quick pass; none of it overturns §1.

- **Hoekstra & Kushner, *JVST B* 20, 1077 (2002)**, "Etching of polysilicon in
  inductively coupled Cl₂ and HBr discharges. III. Photoresist mask faceting,
  sidewall deposition, and microtrenching." [VERIFY-quote — paywalled, abstract
  via search index]: "**The top facet angle was controlled by the surface
  composition at the top of the photoresist lines and the angular dependence for
  etching of the deposited material; the facet being less steep when there was
  more deposition of Si-based byproducts.**" and "the lower facet angle and the
  polysilicon sidewall profile were governed by the feature aspect ratio, the
  sticking probabilities, and fluxes of the depositing material and the depositing
  material etching angular dependence." Note the framing: even in the canonical
  Kushner-group faceting paper, faceting is entangled with *deposition*, and the
  consequences drawn are microtrenching and bowing.
- **Huang *et al.*, *JVST A* 37, 031304 (2019)** [H19]: facets → hot neutrals with
  broad angular distributions → **bowing**; direct line-of-sight from facet into
  the feature is required for direct bowing, otherwise the effect is via
  broadened sidewall-collision angles.
- **Yeom *et al.*, "Role of Oxygen in Amorphous Carbon Hard Mask Plasma Etching",
  *ACS Omega* 8 (2023), doi:10.1021/acsomega.3c02438** (open access,
  PMC10500572): a-C etch rate 25→76 nm/min as O₂ goes 3.5 %→6.5 %; "**the mask
  opening width increased 1.57 times** for the Kr mixture in the 6.5 % O₂
  condition". This is a *mask-open* process, not HARC, but it makes the point
  that a-C aperture widening is an O-chemistry effect on the mask, and Krüger's
  a-C is O-inert at 1e-5 by construction — consistent with his zero-erosion
  target. It matches his own transfer observation ([T] 4835–4838): "**The lack of
  polymer film on the mask leads to increased mask erosion with increasing O₂
  inflow. The amount of mask erosion in the simulated features is not in
  quantitative agreement with the experimental results** … For example, direct
  oxidation of the AC mask by O₂ was not included in the mechanism but may become
  important at larger O₂ flow rates."
- No source found reporting a **measured a-C facet angle under keV Ar⁺ in HARC
  dielectric etch**. [VERIFY] — if this matters later, the search should target
  Lam/TEL/imec HARC papers and a-C sputter-yield-vs-angle measurements rather
  than the simulation literature.

---

## 8. Recommended next actions, in cost order

1. **w_m(dx) convergence ablation** (§6.1). 10 → 5 → 2.5 nm on the ml13 config.
   Krüger ran 1 nm; our 10 nm quantizes the film to ±10 nm of aperture, which
   alone brackets the entire residual. Pre-register: if w_m rises monotonically
   with refinement, close the residual as numerics.
2. **Declare one constant column** (§6.2) and confirm which pair ml13 used.
   Cheap, and the deposition/removal ratio differs by 1.34 vs 2.22 across the
   two published columns.
3. **Run ml16 with the P1a-corrected lift** as already preregistered. §5.4 gives
   the independent physical argument for why it should move the mouth: ~78 % of
   the neck's removal budget is ion-delivered, and the corrected lift raises
   AR-9 sidewall flux by ≈×1.40.
4. **Report the aperture profile, not the scalar.** Both Krüger's model and the
   SEM have a neck minimum 200–270 nm inside the mask with full recovery below.
   If petch's minimum is at the mask top instead, the failure mode is
   "lip over-deposits / lip under-sputters", which is diagnosable per-face and is
   a different bug from "film too thick everywhere". Digitized reference curves
   are in `tmp/mouth_profiles/*.csv`.
5. **Do not pursue mask faceting or feature charging for the mouth.** Neither is
   in his fig-7 model; the first is his admitted un-reproduced feature and the
   second he characterizes as nominal at the mask for HAR.

## 9. Open [VERIFY] items

- Whether MCFPM applies a default sticking to named gas-phase products
  (`SiO2`, `SiF_x`, `CO`, `CO2`) that carry no Appendix-B surface row (§3.2).
- Whether feature charging was silently enabled in the Ch. 6 runs — the absence
  is uniform across thesis chapter and article, but never stated (§2.2).
- The Hoekstra & Kushner 2002 facet-angle quotes are from a search index, not
  from the article text (§7).
- The Fig. 7(a) horizontal scale is assumed uniform and calibrated on w_t = 90 nm;
  the panel is labelled "Not to Scale", which is read here as applying to the
  vertical truncation of the oxide. The *shape* conclusions in §1.5 do not depend
  on this; the absolute nm values do.
- SEM aperture is a valley-FWHM and under-reads the true void by a few nm; the
  39 nm minimum should be read as consistent with, not a correction to, the
  published w_m = 45 nm target.

## 10. Artifacts produced

- `tmp/mouth_profiles/extract_mouth_profiles.py` — digitizer, reproducible from
  `tmp/pdfs/krueger-2024.pdf` p. 8.
- `tmp/mouth_profiles/krueger_fig7a_simulated_aperture.csv` (1070 rows) —
  MCFPM aperture vs depth into the mask.
- `tmp/mouth_profiles/krueger_fig7b_experimental_aperture.csv` (1390 rows) —
  SEM aperture vs depth into the mask.
- `tmp/mouth_profiles/krueger_fig7{a,b}_*_crop.png` — the 600 dpi crops used.

(`tmp/` is gitignored; nothing here is staged for commit.)
