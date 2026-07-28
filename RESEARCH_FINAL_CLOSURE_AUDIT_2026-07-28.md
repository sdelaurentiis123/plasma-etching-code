# Final Closure Audit — Krüger Appendix B vs petch mixed-layer (2026-07-28)

Authority: Krüger PhD thesis 2024 (`tmp/pdfs/krueger_thesis.txt`), Appendix B
Table B.0.1 (pp. 232–260, full table parsed) + Chapter 6 (optimization,
pp. 171–211). petch implementation: `src/petch/mixed_layer.py`,
`src/petch/mixed_layer_mechanism.py`, reflection cascade
`boundary_transport_3d.split_grazing_ion_reflection`.

Columns in the table: `p0  Eth(eV)  n(=q)  E0(eV)  ∠`. `∠=1` = Kress-type
angular law (ref [1]); `∠=2` = second angular law (ref [2]). `(s)`=fresh
surface species, `(xs)`=crosslinked, `X+`=ion, `X#`=hot neutral, `X*`=excited
(activated). "Not listed: surface neutralization for every charge species. EP
= inert generic etch product."

---

## PART 1 — Complete class list from Appendix B (every unique reaction class)

| # | Class | Reaction template | p0 | Eth | n | E0 | ∠ |
|---|-------|-------------------|----|----|---|----|---|
| D1 | Deposition on fresh polymer (s) | CFx(s) + {CF,CF2,CF3} → +radical(s) | 0.1 | — | — | — | — |
| D1b | Deposition on fresh polymer, C2F3 | CFx(s) + C2F3 → +C2F3(s); C2F3(s)+radical | 0.03 | — | — | — | — |
| D2 | Deposition on crosslinked (xs) | CFx(xs) + radical → +radical(s) | 0.02 | — | — | — | — |
| D3 | Deposition on AC mask | AC(s/xs) + {CF,CF2,CF3,C2F3} → AC + radical(s) | **0.2** | — | — | — | — |
| O1 | O-etch of polymer | CFx(s), CFx(xs) + O → EP | **0.0423** | — | — | — | — |
| O2 | O-etch of AC mask | AC(s/xs) + O → CO | 1.0e-5 | — | — | — | — |
| S1 | Polymer sputter, fresh (s) | CFx(s) + ion/# → EP + # | 0.9 | 20 | 0.5 | 500 | 1 |
| S1c | Carbonization branch (CF only) | CF(s) + ion/# → **AC(s)** + F + # | 0.01 | 20 | 0.5 | 500 | 1 |
| S1d | Sputter-defluorination branch | CF2(s)→CF(s)+F ; CF3(s)→CF2(s)+F (+ion/#) | 0.01 | 20 | 0.5 | 500 | 1 |
| S2 | C2F3(s) sputter (heavy, resistant) | C2F3(s) + ion/# → EP + # | 0.8 | 50 | 0.5 | 500 | 1 |
| S3 | Crosslinked polymer sputter (xs) | CFx(xs) + ion/# → EP + # | 0.6 | 50 | 0.5 | 500 | 1 |
| S3d | **De-crosslink** branch (xs→s) | CFx(xs) + ion/# → CFx(**s**) + # | 0.3 | 8 | 0.5 | 500 | 1 |
| S4 | Radical release on O# strike | CFx(s/xs) + O# → CF + O# (non-energetic) | 0.1 | — | — | — | — |
| M1 | AC mask sputter | AC(s/xs) + ion/# → C + # | 0.001 | 200 | 0.4 | 250 | 1 |
| F1 | Fluorination ladder | AC(s)+F→CF(s) 0.01; AC(xs)+F→CF(s) 0.001; CF(s)+F→CF2(s) 0.01; CF2(s)+F→CF3(s) 0.01; CF3(s)+F→CF4 0.01 | — | — | — | — | — |
| C1 | Carbon (sputtered C) redeposition | AC(s)+C→2AC(s) 0.01; SiO2(s)+C→SiO2(s)+AC(s) 0.01; CF2(s)+C→CF2(s)+AC(s) 0.01 | 0.01 | — | — | — | — |
| A1 | **SiO2 ion-activation** | SiO2(s) + ion/# → **SiO2\*(s)** + # | **0.9** | — | — | — | — |
| SP1 | Bare SiO2 physical sputter | SiO2(s) + ion/# → # + SiO2(volatile) | **0.0852** | 70 | 1 | 140 | 1 |
| SP1b | Activated-SiO2\* physical sputter | SiO2\*(s) + ion/# → # + SiO2 | 0.0852 | 70 | 1 | 140 | 1 |
| CH1 | Chemisorption on **bare** SiO2 | SiO2(s)+CF→SiO2CF(s) 0.278; +CF2 0.278; +CF3 0.2; +C2F3 0.2; +C2F4/C3F5/C3F6 0.001 | **0.278**… | — | — | — | — |
| CH2 | Chemisorption on **activated** SiO2\* (two-state, w/ reflection) | SiO2\*+CF→SiO2CF(s) **0.8** (refl 0.2); +CF2 **0.85** (0.15); +CF3 **0.9** (0.1); +C2F3 **0.9** (0.1) | 0.8–0.9 | — | — | — | — |
| DA1 | SiO2\* deactivation (quench) | SiO2\*(s)+F→SiO2(s)+F ; SiO2\*(s)+O→SiO2(s)+O | 1.0 | — | — | — | — |
| EX1 | Complex ion-assisted etch (the oxide-removal channel) | SiO2CFx(s), SiO2CFx\*(s) + ion/# → SiFx + {CO2,CO,COF} + # (+ SiOCF3(s), SiO2CF(s) intermediates) | **0.1471** | 35 | 1 | 140 | **2** |
| EX2 | Final SiOCF3 etch | SiOCF3(s) + ion/# → # + CO + SiF3 | 0.1471 | 35 | 1 | 140 | 2 |

Oxide etch is **polymer-mediated two-path**: (a) direct physical sputter of
bare/activated SiO2 at 0.0852@70eV (no F cost); (b) fluorocarbon radical forms
a surface complex — weakly on bare oxide (0.2–0.278), **strongly on
ion-activated SiO2\*** (0.8–0.9) — then ions sputter the complex at
0.1471@35eV to volatilize SiFx + COx. Activation (A1, 0.9) + two-state
chemisorption (CH2) is the dominant complex-formation route.

---

## PART 2 — Diff against petch (mixed_layer.py + mixed_layer_mechanism.py)

Legend: **IMPL** implemented / **PART** partial or coarse-grained / **ABS** absent.

| Class | Status | Where / constant in petch | Delta vs thesis |
|-------|--------|---------------------------|-----------------|
| D1 on-polymer | IMPL | `KRUEGER_2024_DEPOSITION_ON_POLYMER` CF/CF2/CF3 0.1, C2F3 0.03 | exact ✓ |
| D2 on-crosslinked | IMPL | `KRUEGER_2024_DEPOSITION_ON_CROSSLINKED` 0.02 | exact ✓ |
| D3 on-mask | PART | `KRUEGER_2024_DEPOSITION_ON_MASK` = **0.0842** | thesis Appendix=0.2, **Table 6.5 converged=0.094**; petch value (0.0842) matches neither |
| (deposition on bare oxide) | EXTRA | `KRUEGER_2024_DEPOSITION_ON_SUBSTRATE` CF 0.002…0.001 | **not in thesis** — thesis has NO polymer deposition on bare SiO2 (radicals chemisorb, CH1). petch invention; small |
| O1 polymer O-etch | PART | `oxidation_probability=0.0628` | thesis = **0.0423** (wrong set — see verdict) |
| O2 AC O-etch (1e-5) | PART | mask carbon substrate is O-inert to *substrate* removal ✓, but film-on-mask O-etch uses 0.0628 | mask *material* inertness OK; film constant wrong |
| S1 polymer sputter | IMPL | `kernel_sputter` p0=0.9, Eth=20, q=0.5, E0=500, Kress B=9.3 | exact ✓ (∠=1) |
| S1c carbonization CF→AC | **ABS** | — | no in-situ polymer→AC conversion on oxide |
| S1d sputter-defluorination | PART | film is lumped (n_c_film/n_f_film); sputter removes both, no stepwise CFx→CF(x-1) ladder | coarse-grained |
| S2 C2F3 heavy sputter (0.8@50) | PART | single film sputter law (0.9@20) for all film | C2F3 distinct resistance not modeled |
| S3 xs sputter (0.6@50) | PART | n_xl tracked but sputter law is single 0.9@20; xl only lowers deposition attachment | crosslinked sputter-resistance not differentiated in removal |
| S3d **de-crosslink (0.3@8)** | **ABS** | crosslinking `xl_rate` is one-way; no ion reversion xs→s | missing reverse channel |
| S4 radical release on O# | ABS | — | negligible |
| M1 AC mask sputter | IMPL | `kernel_ac` p0=0.001, Eth=200, q=0.4, E0=250 | exact ✓ |
| F1 fluorination ladder | PART/NA | `absorb_f` = fluorine_film_sticking·F·θ; but Krüger config sets `fluorine_species=()` (no atomic-F flux) → absorb_f=0. F reaches layer only bound to C via CH1 | ladder effectively absent in Krüger config (by construction — no free F species) |
| C1 sputtered-C redeposition | **ABS** | no gas-phase C species tracked | sidewall/mask carbon redep missing |
| A1 **SiO2 activation (0.9)** | **ABS** | no SiO2\* state in mixed layer | missing |
| SP1 bare SiO2 sputter | PART | `kernel_bare` p0=**0.0909**, Eth=70, E0=140, n=1 | threshold/E0/n ✓; magnitude thesis=**0.0852** (wrong set) |
| SP1b activated-SiO2\* sputter | N/A | folded into SP1 (no SiO2\* state) | same yield anyway (0.0852) |
| CH1 bare chemisorption | PART | `KRUEGER_2024_CHEMISORPTION_PROBABILITY` CF/CF2 **0.2729**, CF3/C2F3 0.2, heavies 0.001 | thesis CF/CF2=**0.278**; petch flat-0.2729 (wrong set) |
| CH2 **activated chemisorption (0.8–0.9)** | **ABS** | only bare CH1 applied everywhere | the DOMINANT complex channel is missing |
| DA1 SiO2\* deactivation by F/O | ABS | (no SiO2\* state) | missing with A1 |
| EX1 complex ion-etch | PART | `kernel_complex` = 0.1384·flux·εdep/εdep(140) | thesis p0=**0.1471** Eth=35 n=1 ∠=**2**; petch uses ZBL deposited-energy shape (no Eth=35 threshold, no ∠=2 law) and magnitude 0.1384 |
| EX2 SiOCF3 final etch | PART | folded into single complex removal → SiF4 | intermediate speciation (SiOCF3, COF) not resolved (declared) |

Reflection cascade (`split_grazing_ion_reflection`): **IMPL** and faithful —
Huang Eq. 2.34 (E_ts=100, E_c=10, θ_c=70°), leftover-rule continue weight
`1 − clip(0.9·kress, 0,1)` (interpretation B'), multi-bounce (≤8),
diffusive-drop below cutoffs, y-symmetrize. Hot-neutral rows generated. The
`react=0.9·kress` uses the polymer-sputter p0 as a single surface-agnostic
consumption prob (modeling choice, not surface-resolved).

---

## PART 3 — SET QUESTION VERDICT

**The thesis fig-7 / base-case validation (45 nm mouth, 825 nm depth, 850 nm
remaining mask) was generated with the APPENDIX-B = converged constants, NOT
the JVST-paper "optimized" set petch currently uses.**

Evidence chain (thesis Chapter 6):

1. Appendix B header (l. 5330): *"Table B.0.1 contains the surface reaction
   mechanism used in this work **after convergence**."* — "convergence" = the
   Chapter 6 hybrid gradient-descent → Nelder-Mead optimization.

2. **Table 6.4 (target metrics)** = the validation target itself:
   `wm=45 nm, wt=90 nm, wf=90 nm, hf=825 nm, hm=850 nm, ah=0`. This IS the
   45 nm/825 nm base case.

3. **Table 6.3 (the 5 tuned parameters)** are *exactly* the disputed constants:
   `ps,SiO2` (bare sputter), `ps,SiO2CFXY` (complex sputter), `pp,SiO2`
   (complex formation), `pe,poly` (O polymer etch), `pd,poly-AC` (deposition on
   mask).

4. **Table 6.5 (final converged values)** — the set that produced the epoch-200
   validation feature (Fig. 6.16, matched to the SEM):

   | Parameter | Table 6.5 (converged, = fig-7) | petch current | thesis Appendix B row |
   |-----------|-------------------------------|---------------|-----------------------|
   | ps,SiO2 (bare sputter) | **0.0852** | 0.0909 ✗ | 0.0852 |
   | ps,SiO2CFXY (complex sputter) | **0.1471** | 0.1384 ✗ | 0.1471 |
   | pp,SiO2 (complex formation) | **0.278** | 0.2729 ✗ | 0.278 (bare) |
   | pe,poly (O polymer etch) | **0.0423** | 0.0628 ✗ | 0.0423 |
   | pd,poly-AC (dep on mask) | **0.094** | 0.0842 ✗ | 0.2 |

So on the three cleanest levers (bare sputter, complex sputter, O-etch) the
converged/fig-7 values are the **Appendix** numbers (0.0852 / 0.1471 / 0.0423),
and petch is currently on the *other* set (0.0909 / 0.1384 / 0.0628). petch's
own code comment (`mixed_layer_mechanism.py` l. 523-525: "the optimized set is
what produced fig-7; appendix-converged 0.2 is the alternative") has the
attribution **backwards** — Table 6.5 IS the converged set and it matches
Appendix, not the JVST-paper optimized values.

Note the one genuine internal thesis inconsistency: **pd,poly-AC** = 0.094
(Table 6.5) vs 0.2 (Appendix B AC-deposition rows). Table 6.5 is the
optimizer output that produced the validation figure, so **0.094** is the
fig-7 value; the Appendix 0.2 is likely a stale/pre-optimization row. petch's
0.0842 matches neither — retarget to 0.094.

The activated chemisorption question resolves cleanly too: the fig-7 mechanism
uses the FULL two-state chemisorption — bare `pp,SiO2 = 0.278` (the optimized
CH1) **plus** the fixed activated CH2 (0.8/0.85/0.9/0.9) and A1 activation
(0.9). CH2/A1 were NOT free parameters (not in Table 6.3), so they sit at their
Appendix values in every converged run. petch implements only a single flat
chemisorption (0.2729) and no activation — i.e. petch is missing the dominant
half of the fig-7 complex-formation channel.

---

## PART 4 — ABSENT/PARTIAL ranked by expected mouth impact, with direction

Residual targets: **mouth 25→45 nm** (petch too NARROW, needs more lateral
oxide removal near the top), **mask-remaining →850 nm** (hm), **O-shape dip**
(necking depth/shape).

### 1. CH2 + A1 — ion-activated SiO2\* two-state chemisorption  [ABSENT — TOP LEVER]
- Physics: ion/hot-neutral bombardment activates oxide (A1, 0.9); activated
  sites chemisorb CF/CF2/CF3/C2F3 at **0.8–0.9** vs bare **0.2–0.278** (≈3–4×).
- Mouth direction: **WIDENS (helps 25→45).** The upper sidewall receives the
  highest grazing-ion + reflected-hot-neutral flux (the cascade petch already
  builds), so it activates most, chemisorbs ~3× more complex, and the complex
  is sputtered off → strong *lateral* oxide removal exactly at the mouth. This
  is the single largest missing lever and the most likely cause of the 25→45
  gap. Magnitude guess: **+10–20 nm** mouth (closes most of the residual).
- Mask/O-shape: modest; concentrates removal near ion-rich regions, may
  slightly deepen necking.

### 2. Constant-set correction to the fig-7/Table-6.5 values  [PARTIAL — LARGE]
- O-etch 0.0628→**0.0423** (lower): LESS polymer removed by O → thicker
  sidewall/mask polymer. Mouth direction: **NARROWS** (works against 25→45 on
  its own) but **protects mask (helps hm→850)** and deepens necking (O-shape).
  This is why it can only be lowered *together with* CH2/A1, which over-supply
  the lateral etch — matching how the thesis converged to 0.0423 AND 45 nm.
- Complex sputter 0.1384→**0.1471** (+6%) and bare sputter 0.0909→**0.0852**
  (−6%): small net; complex up slightly raises vertical + complex-mediated
  lateral etch (mild mouth widen), bare down slightly reduces pure physical
  sputter. Magnitude: few nm.
- Chemisorption 0.2729→**0.278** flat + fixing CF3/C2F3: negligible alone (the
  activated CH2 dwarfs it).
- Mask deposition 0.0842→**0.094**: slightly more mask polymer → better mask
  survival (hm), slightly narrower mouth. Few nm.

### 3. S3d — de-crosslinking (xs→s, 0.3@8 eV)  [ABSENT — MEDIUM]
- Ions revert crosslinked skin to fresh polymer, which re-etches by O and
  re-sputters more readily (fresh 0.9@20 vs xs 0.6@50). Net: keeps the
  ion-exposed upper sidewall polymer *thinner/more mobile*.
- Mouth direction: **WIDENS slightly** (less protective crosslinked skin at the
  mouth). Magnitude: **+2–5 nm**. Also modulates O-shape (necking turnover).

### 4. S1c — carbonization CF(s)+ion→AC(s)+F (0.01@20)  [ABSENT — MEDIUM, opposes]
- Converts ion-struck sidewall polymer to graphitic AC that is near
  sputter-immune (0.001@200) and O-inert (1e-5) → passivating carbon skin.
- Mouth direction: **NARROWS / protects** (opposes widening) but builds a
  robust cap that helps **mask survival (hm→850)** and can pin the necking
  location (O-shape). Magnitude on mouth: **−2–5 nm**; on hm meaningful. Include
  for correct mask/necking even though it costs a little mouth.

### 5. S2/S3 differentiated film sputter resistance (C2F3 0.8@50, xs 0.6@50)  [PARTIAL — SMALL]
- Heavy/crosslinked polymer harder to sputter than fresh (petch uses one
  0.9@20 law). Thickens protective film where ion energy is marginal (deep
  sidewall, mouth shoulders). Mouth: small **narrowing**; mainly O-shape/mask.
  Magnitude: 1–3 nm.

### 6. C1 — sputtered-carbon redeposition as AC (0.01)  [ABSENT — SMALL]
- Mask-sputtered C redeposits as AC on sidewalls/oxide → extra passivation and
  apparent mask-material recycling. Mouth: small **narrowing**; helps hm.
  Magnitude: 1–3 nm.

### 7. EX1 angular law ∠=2 + Eth=35 threshold form  [PARTIAL — SMALL/SHAPE]
- petch's complex removal uses a ZBL deposited-energy shape rather than the
  thesis threshold-power law (Eth=35, n=1) with the ∠=2 angular dependence.
  Mostly affects the *taper/angle* response (the thesis itself flags taper as
  sensitive to the chemical-sputter angular law, l. 4790s) → O-shape, minor
  mouth.

**One-commit batch, in impact order:** (1) add SiO2\* state + A1 activation
(0.9) + CH2 activated chemisorption (0.8/0.85/0.9/0.9 with reflection
remainder) + DA1 F/O quench; (2) retarget the five constants to Table 6.5
(0.0852 / 0.1471 / 0.278 / 0.0423 / 0.094) and drop the JVST-paper set; (3) add
S3d de-crosslink (0.3@8) and (4) S1c carbonization (0.01@20). Items 5–7 are
polish for O-shape/taper and mask, lower priority.
