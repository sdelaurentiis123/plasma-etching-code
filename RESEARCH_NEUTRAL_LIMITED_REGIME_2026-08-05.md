# How validated fluorocarbon-HARC models produce a NEUTRAL-LIMITED floor

**Date:** 2026-08-05
**Scope:** source-side only. A parallel code fork audits petch's F-budget receipts.
**Question:** our mixed layer is measured exactly ion-limited (<1% response over a 50x radical
range) at the floor, where Krueger and Huang say the real process is neutral-transport-limited.
What mechanism makes *their* models neutral-responsive, what N/I does the published literature
measure as the ion-limited / neutral-limited crossover, and what is the doctrine-compliant fix shape?

**Convention used throughout.** `N` = thermal reactive neutral (radical) flux at the surface in
question. `I` = energetic-particle flux (ions **plus** hot neutrals; in the MCFPM lineage they carry
identical surface reactions — Huang App. E note (a), L10352-10355: *"All ions neutralize on
surfaces, returning to gas phase as their hot neutral partner. Ions and their hot neutral partners
have the same surface reactions with the same probability."*). Where a source reports ions only,
it is stated. All local line numbers refer to files in `tmp/pdfs/`.

---

## 0. Executive answer

1. **Mechanism.** Neutral limitation in Krueger/Huang (MCFPM) and in Kwon/Sawin (TML) is *not* a
   flux threshold or a starvation switch. It is **finite site turnover**: a radical arrival
   *consumes* a pristine substrate site to create a complex; the ion removes the complex and
   **regenerates a pristine site**, which cannot be removed again by the chemical channel until
   another radical arrives. Steady-state coverage is
   `theta = s*Gamma_n / (s*Gamma_n + Y*Gamma_i)`, and the etch rate is coverage-weighted between a
   fast chemical channel and a slow physical-sputter channel. The removal law is **linear in layer
   halogen concentration** (Kwon TML `r4 = beta1 * x_F * ...`) with **no cap** — the cap is one
   non-reflective event *per incident particle*, not per site.
2. **Measured crossover.** For SiO2 in an F/Ar+ beam at 350 eV (Gray 1993, replotted as Kwon
   Fig. 3.4): sputter floor `Y = 0.28` at `F/Ar+ -> 0`, half-rise at `F/Ar+ ~ 25-30`, ~90% of
   plateau at `F/Ar+ ~ 100`, plateau `Y ~ 1.1`. **Full neutral-off-to-saturated dynamic range =
   3.9x.** For Si in a Cl/Ar+ beam (Chang 1997, Huard Fig. 1.10) the yield "depends strongly on
   `Gamma_n/Gamma_i` for low and moderate values (< 100)" and saturates above.
3. **Floor target.** Huang's HARC base case delivers **N/I (radicals / (ions+hot neutrals)) = 0.33-0.5
   at AR > 10** and ~0.64 at AR 40, against **30 at the wafer** — a ~75x differential collapse. The
   floor is deep in the neutral-limited branch of the measured curve.
4. **Verdict on petch's <1% response.** No published SiO2 surface model or beam dataset is that
   flat *anywhere*. A 50x radical sweep must move the yield by roughly **-55 to -60%** on the Gray
   1993 curve regardless of which end of the sweep you anchor. `<1%` is a model defect, not physical
   saturation. The doctrine-compliant repair is Section 5.

---

## 1. KRUEGER — verbatim

Source: `tmp/pdfs/krueger_thesis.txt` (P. Krueger, PhD thesis, U. Michigan, Kushner group;
HPEM reactor + MCFPM feature scale; HAR SiO2 etch with amorphous-carbon mask,
90 nm opening, 850 nm AC mask, 60 s etch — L4371-4374).

### 1.1 The anchor passage (L4879-4889) — the process is neutral-transport-limited

> "The profiles, experiment and simulation, for Plf = 0 kW produce total clogging of the mask
> opening, indicating that ion energy plays an important role in removing excess polymer. The
> remaining cases (Plf = 4, 6 and 8 kW) have unclogged features and full etching with unexpected
> little variation as a function of LF power. In the experimental data a doubling of the Plf (4 kW to
> 8 kW) produces few differences in the final features. **These trends indicate that above a certain
> threshold energy the etch progression and the mask removal process are not ion starved, but rather
> limited by neutral gas transport.** To some degree this trend is reproduced by the simulations where
> etch depth does increase with increasing low frequency power however the rate of increase is
> substantially sublinear. **These outcomes indicate that the effect of ion energy (for example in
> sputter yield or related processes) might be overestimated in the mechanism.**"

Two load-bearing facts here, and the second is the one that matters most for petch:

- The *experiment* is flat in ion energy across a 2x LF power doubling — the signature of a
  neutral-limited floor.
- **Krueger's own ML-optimized mechanism is diagnosed by its author as over-weighting the ion-energy
  channel.** petch inherits Krueger Appendix B row-for-row. petch's ion-limited floor is therefore
  an *amplification of a defect the source author already flagged*, not a departure from the source.

### 1.2 Krueger's regime taxonomy (L1547-1571)

> "In HAR etch process, as the etch front propagates downwards, a drastic reduction in etch rate is
> often observed [95-97]. Generally, this is a consequence of transport limitations of neutrals as
> well as ions. Neutral transport to the bottom of the feature is limited by a host of factors. In
> scenarios where the mean free path of a neutral thermal particle is on the same scale or larger
> than geometric restrictions, as is the case in narrow high aspect ratio features, transport is
> governed by Knudsen diffusion which severely limits neutral gas transport deep into the feature.
> This is coupled with the loss of radicals during sidewall collisions, the number of which increases
> with increases etch depth and aspect ratio [98,99]. ... The reduction of etch rate due to lack of
> neutral and ion transport to the bottom of the feature are respectively known as the **neutral and
> ion starved regime**."

### 1.3 The mechanism requires radicals as reactants (L4015-4019)

> "The SiO2 etch mechanism contains a SiO2-polymer complex **which requires fluorocarbon radical
> fluxes as reactants**. Since these fluxes have a small increase with phi (Figure 5.15a), the
> increased availability of reactants could also play a role in the increased etch rate **if the etch
> progression is flux limited as is often the case in HAR features**."

And (L4384-4392):

> "SiO2 can be removed through physical sputtering by energetic ions and hot neutrals. The sputtered
> products can be redeposited on other surfaces. **Unsaturated fluorocarbons can chemisorb on the
> SiO2 to form an oxide-fluorocarbon complex. This complex is in turn easier to sputter based on a
> modified threshold (reduced total binding energy) and an overall higher reaction probability as the
> site has a lower binding energy.**"

### 1.4 Which rows saturate — Krueger Appendix B (the rows petch replays)

Thermal-radical passivation of pristine oxide (L6006-6012), probabilities as optimized:

| row | p0 |
|---|---|
| `SiO2(s) + CF -> SiO2CF(s)` | 0.278 |
| `SiO2(s) + CF2 -> SiO2CF2(s)` | 0.278 |
| `SiO2(s) + CF3 -> SiO2CF3(s)` | 0.2 |
| `SiO2(s) + C2F3 -> SiO2C2F3(s)` | 0.2 |
| `SiO2(s) + C2F4 / C3F5 / C3F6` | 0.001 |

Further passivation of an *already-complexed* site (L6556-6564) is **0.0002-0.002**, i.e. ~140x
smaller than first passivation. Fluorination of the complex (L6548-6554) is 0.1.

**This is the saturation.** Once a site is complexed it is essentially inert to further carbon
arrivals. Extra radical flux buys nothing; only *turnover* (ion removal exposing a fresh pristine
site) creates new demand. That is the structural origin of neutral limitation.

Removal rows (L6180+, all identical across every ion and hot-neutral partner):

```
SiO2CF(s)    + I+/I#  ->  SiF  + CO2       0.1471  Eth=35  n=1  Er=140  angclass=2
SiO2CF2(s)   + I+/I#  ->  SiF2 + CO2       0.1471  35  1  140  2
SiO2CF3(s)   + I+/I#  ->  SiF3 + CO2       0.1471  35  1  140  2
SiO2C2F3(s)  + I+/I#  ->  SiOCF3(s) + CO   0.1471  35  1  140  2
SiOCF3(s)    + I+/I#  ->  SiF3 + CO        0.1471  35  1  140  2
```

versus bare-oxide physical sputtering (L5904-5941): `p0 = 0.0852, Eth = 70, Er = 140, angclass = 1`.

**Stoichiometry carried:** one carbon-bearing radical arrival per Si removed. `SiO2CF -> SiF + CO2`
consumes 1 C + 1 F and removes one SiO2 formula unit; `SiO2CF3 -> SiF3 + CO2` consumes 1 C + 3 F.
The `C2F3` branch takes **two** ion strikes per Si (`SiO2C2F3 -> SiOCF3(s) + CO`, then
`SiOCF3 -> SiF3 + CO`). So: **1 radical arrival per removal, 1-3 F per removal, and for the
2-carbon branch 2 ion strikes per removal.** No row anywhere requires *multiple* radical arrivals
to make one removal — the multiplicity that produces neutral limitation is **arrival-per-turnover**,
not arrivals-per-event.

**Energy scaling at the floor.** Krueger's law is `p0 * (E^q - Eth^q)/(Er^q - Eth^q)` with q = 1
here (col 3 = 1). At Krueger's wafer ion energies (up to 4800 eV, L4364) the complex row evaluates
to `0.1471 * (1500-35)/(140-35) = 2.05` and the bare row to `0.0852 * (1500-70)/(140-70) = 1.74`,
both far above unity. **They do not scale freely — they are capped**, by the MCFPM
reaction-selection normalisation (Huang L2311-2320, quoted in §2.5). One incident particle causes at
most one non-reflective event. After capping, the complex:bare contrast at the floor is set by the
*ratio* and by the angular class, not by energy.

### 1.5 Krueger base-case wafer fluxes (Table 6.1, L4342-4349)

| species | flux (cm^-2 s^-1) |
|---|---|
| C3F4 | 9.5e16 |
| C2F3 | 6.8e16 |
| CF   | 4.4e16 |
| CF2  | 9.4e16 |
| CF3  | 8.4e15 |
| O    | 7.7e16 |

Total FC neutrals = 3.29e17 cm^-2 s^-1. Ion flux is not tabulated but Krueger states the regime
(L2519-2521): *"These conditions produce ion fluxes to the wafer of 10^15-10^16 cm^-2 s^-1."*
Hence **wafer-plane N/I = 33 to 330**.

petch's replay (`scripts/floor_delivery_scan.py` L19-20) uses CF 4.4e20, CF2 9.4e20, C2F3 6.8e20,
CF3 8.4e19, C3F4 9.5e20 (declared inert), O 7.7e20 m^-2 s^-1 and `ION_SRC = 9.6e19 m^-2 s^-1`
(= 9.6e15 cm^-2 s^-1, inside Krueger's stated band). Reactive-FC N/I at the wafer = **22.1**;
including C3F4, **32.0**. That matches Huang's base case (30) exactly. **The wafer plane is right.
Everything at issue is what happens between the wafer and the floor.**

---

## 2. HUANG — verbatim

Source: `tmp/pdfs/huang_thesis.txt` (C. Huang, PhD thesis, U. Michigan, Kushner group; Ar/C4F8/O2
CCP, HAR SiO2 via, photoresist mask of AR 13, oxide etched to AR 40+).

### 2.1 The anchor passage (L5391-5403) — why ions and neutrals attenuate differently

> "By definition, ions striking the etch front are those not having had collisions with the sidewalls
> as, in this mechanism, ions neutralize upon striking a surface and become hot neutrals. A given ion
> can only strike the etch front once, so the flux of ions to the etch front to some degree
> represents the decreasing view angle of the etch front subtending the ion angular distribution at
> that AR. For this reason alone, the ion flux to the etch front decreases with AR. **On the other
> hand, thermal neutral species undergo diffusive scattering on the sidewalls are conductance limited
> in reaching the etch front as the AR increases.** If a neutral has a non-unity reaction probability,
> once at the bottom of the feature that particle can strike the etch front multiple times. **Each
> strike of the etch front by a neutral particle (either hot or thermal) increments the flux-count.**"

(That last sentence is a flux-accounting caveat: Huang's "flux" is a strike count, not a first-hit
count. Any grade against his absolute numbers must match the observable —
`scripts/matched_beam_cascade_grade.py` already records this.)

### 2.2 The floor N/I numbers (L5430-5478) — the quantitative target

Verbatim, with the numbers extracted:

> "As the etch depth increases from 0 to 1,200 nm (AR = 10), the flux of CFx and CxFy radicals to the
> etch front decreases from 3.1 x 10^16 to 0.4 x 10^16 cm^-2 s^-1. This decrease is due to
> consumption by deposition as polymer on the sidewalls and diffusive scattering which reflects the
> neutrals out of the feature before reaching the etch front (**neutral conduction limit**)."
> (L5456-5462)

> "As the etch depth further increases from 1,200 to 4,800 nm (AR from 10 to 40), the flux of CFx and
> CxFy radicals to the etch front decreases from 4.4 x 10^15 to 0.9 x 10^15 cm^-2 s^-1, resulting in
> a surface complex layer with little overlying polymer at the etch front. Due to the limited
> availability of CFx and CxFy to thicken the polymer, the polymer is rapidly removed by ions and hot
> neutrals. **Deep in the feature (AR > 10), the ions and hot neutrals have larger fluxes to the etch
> front than CFx and CxFy by a factor of 2 - 3.**" (L5470-5478)

> "This shadowing contributes to a decrease in ion flux to the etch front from 2.0 x 10^15 to 0.3 x
> 10^15 cm^-2 s^-1. ... The flux of hot neutrals to the etch front increases from 3.1 x 10^15 to 8.0
> x 10^15 cm^-2 s^-1 as the etch depth increases from 0 to 480 nm (AR = 4). As the feature further
> deepens to 960 nm (AR = 8), the flux of hot neutrals surpasses the fluxes of thermal CFx and CxFy
> radicals whose fluxes are conductance limited. ... As the etch depth increases from 480 to 4,800 nm
> (AR = 40), the flux of hot neutrals to the etch front decreases to 1.1 x 10^15 cm^-2 s^-1."
> (L5433-5450)

> "The flux of CFx and CxFy to the top of the PR is 3.0 x 10^17 cm^-2 s^-1." (L5451-5453)

**Assembled N/I ladder (radicals / (ions + hot neutrals)):**

| location | radicals (cm^-2 s^-1) | ions | hot neutrals | **N/I** |
|---|---|---|---|---|
| wafer / top of PR | 3.0e17 | 4.2e15 (PR height 0, L5413-5415) | 0 | **~71**; base-case (CFx+CxFy)/ion is stated as **30** (L502-504, L5696) |
| top of oxide (AR 0, under PR AR 13) | 3.1e16 | 2.0e15 | 3.1e15 | **6.1** |
| **oxide AR 10** | **4.4e15** | ~1e15 | ~8e15 | **0.33 - 0.5** (verbatim "factor of 2 - 3") |
| oxide AR 40 | 0.9e15 | 0.3e15 | 1.1e15 | **0.64** |

**Differential collapse, wafer -> AR 10: N/I falls from 30 to ~0.4 = ~75x.**
**Absolute floor target at AR 10-40: energetic:radical = 1.5 : 1 to 3 : 1.**

This is the single most useful number in this document. If petch's floor sits at energetic:radical
= 100-144 : 1, it is **30-90x too energetic-rich (equivalently too radical-poor) relative to Huang's
published HARC base case**. If instead petch's floor N/I is 100-144 (radical-rich), it is
**200-400x too radical-rich**. Either reading is a two-order-of-magnitude miss against the same
target.

### 2.3 Where the floor radicals actually come from (L5716-5748) — the likely petch gap

> "These energetic ions and their hot neutral partners can remove oxide directly through physical
> sputtering or remove the complex through chemically enhanced sputtering. **After losing energy
> through several collisions with the sidewalls and etch front, these energetic species become
> thermal CFx and CxFy radicals, which can passivate the oxide surface or deposit as polymer.** In the
> base case, the majority of the CFx and CxFy radicals at low AR (< 5) ... originate from the thermal
> neutrals incident into the feature from the plasma. **As the AR increases to greater than 10, the
> neutralized and thermalized CFx+ and CxFy+ ions become the main source (> 95%) of radicals reaching
> the etch front.**"

> "**As the majority of the CFx and CxFy radicals reaching deep into the feature (AR > 10) originates
> from CFx+ and CxFy+, the flux of CFx and CxFy to the etch front at high AR increases by
> approximately 3 times when increasing fractional CFx+ and CxFy+** [15% -> 60% of total ion flux].
> This increase in fluorocarbon radicals to the etch front results in more surface passivation and
> polymer deposition which requires additional ions to remove the oxide. The neutral and thermalized
> partners of other ions are non-reactive species and diffuse out of the feature with no surface
> reactions (only scattering at the surface). ... **Thus, the etch rate of oxide decreases by about
> 30% when increasing the fraction of CFx+ and CxFy+**."

Three consequences, all directly actionable:

- **The HARC floor's radical budget is a product of the ion cascade, not of conducted thermal
  radicals.** Above AR 10 it is >95% ion-derived. A model that does not thermalize spent energetic
  FC species back into the *radical* ledger will show an unphysically radical-poor floor. A model
  that returns them without any conductance loss will show an unphysically radical-rich floor.
- Base case is 30% CFx+/CxFy+ of total ion flux; that fraction is the control knob on the floor
  radical budget and it is a *plasma* quantity, not a surface one.
- Note the sign: **more floor radicals -> ~30% LOWER etch rate**, because polymer thickening starts
  to cost ions. The neutral-limited branch is not monotone-forever; the HARC optimum sits near the
  knee. Any fix that makes petch neutral-responsive must reproduce *both* signs.

### 2.4 The mechanism, as reactions (L3235-3272, App. E Table E.2 L10194-10350)

```
(3.28)  SiO2(s)          + I+   -> SiO2 + I(h)                    physical sputter, p0=0.9, Eth=70, Er=140
(3.29)  SiO2(s)          + CFx  -> SiO2CFx(s)                     passivation  (CF 0.4, CF2 0.3, CF3 0.2, C2F3 0.2)
        SiO2*(s)         + CFx  -> SiO2CFx(s)                     on ACTIVATED oxide: 0.9
(3.30a) SiO2CFx(s)       + I+   -> SiFx + CO2 + I(h)              chem. sputter, p0=0.75, Eth=35, Er=140
(3.30b) SiO2C2F3(s)      + I+   -> SiOCF3(s) + CO + I(h)
(3.30c) SiOCF3(s)        + I+   -> SiF3 + CO + I(h)
(3.31)  SiO2(s)/complex  + I+   -> activated site + I(h)          activation, p0=0.9 (oxide) / 0.1 (complex, Eth=5, Em=70)
(3.32)  SiO2CmFn*(s)     + CFx  -> SiO2CmFn(s) + P(s)             polymer nucleation on activated complex, 0.001-0.002
(3.33)  P(s)             + CFx  -> P(s) + P(s)                    polymer growth, 0.001-0.002 (activated P*: 0.01-0.02)
(3.34)  P(s) + F -> CF2 ;  P(s) + O -> COF                        polymer erosion (F: 0.001/0.03, O: 0.5/0.9)
(3.35)  P(s) + I+ -> I(h) + CF2                                   polymer sputter, p0=0.3, Eth=30
```

"Further passivation of complex" rows (L10214-10222) are all **1e-4** — the same saturation Krueger
carries at 2e-4. Fluorination of the complex (L10224-10227) is 0.1.

**Does complex formation require multiple radical arrivals per removal?** No — one carbon-bearing
arrival creates one removable complex. But the *pristine* sticking probabilities are 0.2-0.4, so
**2.5 to 5 radical arrivals per successful passivation** at an un-activated site (0.9 -> ~1.1
arrivals at an ion-activated site). And the site is single-use: after `SiO2CFx + I+ -> SiFx + CO2`
the site returns to pristine `SiO2(s)`. **Arrivals per removal is therefore ~1.1 to 5 depending on
activation state, and the site must be re-passivated after every single removal event.**

### 2.5 Yield capping at the floor (L2277-2320) — the answer to "does yield scale with energy?"

> "the reaction yield for a particle incident onto a surface with an incident energy of Ei and an
> incident angle of theta with respect to the local surface normal is determined by
> `p(Ei, theta) = p0 * ((Ei^n - Eth^n)/(Er^n - Eth^n)) * f(theta)`, where Eth is the threshold
> energy, Er is a reference energy, p0 is the yield at the reference energy, n is the energy
> dependent exponent (typically 0.5), and f(theta) is the relative probability at angle of
> incidence theta."

> "The angular dependence is often different between direct physical sputtering and chemically
> enhanced etching. ... **For physical sputtering, f(theta) is an empirical function with a maximum
> at 60 deg, reduced probability at normal incidence and zero probability at grazing incidence.
> For chemically enhanced etching, f(theta) is unity for normal incidence and angles up to 45 deg,
> with a monotonic roll-off to zero probability at grazing incidence.**" (L2290-2296)

> "**Since the probability of a particle striking the surface upon arrival is, by definition, unity,
> the following procedure is followed to normalize selection of reaction probabilities. The
> cumulative yield of all allowed processes for the energy and angle of incidence is computed.** For
> all such interactions, there is an elastic collision -- meaning a reflection from the surface
> without changing the state of the surface. If the cumulative yield of non-reflective processes is
> less than unity, then the elastic scattering yield is increased so that the cumulative yield is
> unity. If the cumulative yield is greater than unity, the elastic yield is reduced so that the
> cumulative yield is unity. **If after scaling the elastic yield to zero, the cumulative yield is
> still greater than unity, then the yields of all processes are scaled to provide a unity cumulative
> yield.**" (L2311-2319)

So at HARC energies **the yield is capped, not scaled** — one non-reflective event per incident
particle. petch already implements this (`src/petch/mixed_layer.py::_unit_cumulative_scale`, which
quotes this passage verbatim at L186-192).

**Consequence that matters for petch's flat floor.** With energy capped, the *only* surviving
neutral sensitivity in the Kushner mechanism at keV is the ratio of the complex row to the bare row
and the difference in their angular classes:

- Krueger: complex 0.1471 vs bare 0.0852 = **1.73x**.
- Huang: complex 0.75 vs bare 0.9 = **0.83x** (!) — Huang's bare-oxide sputter p0 is *higher*; his
  neutral sensitivity at high energy comes entirely from the lower threshold (35 vs 70 eV), which
  is irrelevant at 1-4 keV, and from the angular class.
- **Angular class is therefore the dominant remaining lever.** At the HARC floor the ion angular
  distribution is <4 deg (Huang L5485) i.e. essentially normal incidence, where the physical-sputter
  class is at its *minimum* and the chemical class is at unity. In a Kress-type class-1 form
  normalised to its 60-deg peak, `f_phys(0)/f_phys(60) ~ 1/4`. petch currently normalises **both**
  classes to `f(0) = 1` and flags this explicitly as a declared risk
  (`src/petch/mixed_layer.py` L236-242, commit `0a1c809`). Under the peak-normalised convention the
  complex:bare contrast at the HARC floor goes from ~1.7x to ~**7x** — a factor-of-4 recovery of
  neutral sensitivity that is *entirely inside the published table*, no new constant.

### 2.6 Huang's own ARDE attribution — honest caveat

Huang does **not** attribute his HARC ARDE primarily to neutral supply. L5573-5586:

> "The instantaneous etch rate of SiO2 as a function of AR shown in Fig. 6.8 shows typical ARDE...
> **The etch rate then decreases by 80% by the time the AR reaches 40**, following a similar trend as
> the fluxes and power densities of energetic species to the etch front... With the delivery of both
> energetic species and power to the etch front decreasing with increasing AR, the etch rate
> decreases resulting in ARDE. There is also a contribution to ARDE from redeposition of etch
> products, **however the more dominant source of ARDE is the reduction in power to the etch front.**"

and the chapter summary (L5998-5999):

> "**ARDE results from both lack of neutral radicals by conductance limits and by a decrease in the
> power delivered to the etch front as AR increases.**"

and (L5990-5997):

> "For sufficiently large AR, the neutral radicals reaching the bottom of the feature originate from
> neutralized ions which can overcome conductance limits by virtue of their initially anisotropic
> trajectories. **The dominant oxide removal process transitions from chemical to physical sputtering
> as the fluxes of energetic species to the etch front surpass that of radicals.** As the AR further
> increases, even the energetic species with narrow angular distributions are scattered by the
> sidewalls, losing energy and resulting in an etch stop."

Read carefully this is *consistent* with, not contrary to, the neutral-limited thesis: the etch
transitions **into** the physical-sputter branch precisely because the surface can no longer be kept
passivated. Both mechanisms must be present for the profile to be right. Also L5662-5664:

> "**With the flux of thermal neutrals being conduction limited to the bottom of the feature and not
> significantly increasing with 5 MHz power, the etching is neutral-starved at the higher powers.**"

---

## 3. THE BEAM LITERATURE — measured yield vs N/I, and the crossover

### 3.1 SiO2 + F + Ar+ : the primary number (Gray 1993, replotted by Kwon)

- **Data:** D. C. Gray, I. Tepermeister, H. H. Sawin, "Phenomenological modeling of ion-enhanced
  surface kinetics in fluorine-based plasma etching," *J. Vac. Sci. Technol. B* **11**, 1243-1257
  (1993). Multiple-beam apparatus; F-atom and Ar+ fluxes varied independently over several orders of
  magnitude on undoped polysilicon and SiO2.
- **Replot + model:** O. Kwon, "Surface kinetics modeling of silicon oxide etching in fluorocarbon
  plasmas," ScD thesis, MIT Dept. of Materials Science and Engineering, 2004 (advisor H. H. Sawin),
  https://hdl.handle.net/1721.1/28353 — **Fig. 3.4, p. 76**: "A translating mixed-layer model
  calculation result for silicon oxide etching with fluorine chemistry compared with experimental
  data. **Ion bombardment energy is 350 eV.**" Axes: `Etching Yield [SiO2/Ar+]` vs `Flux Ratio
  [F/Ar+]`, x from 0 to 600.

**Digitized from the figure** (data markers, "Gray (1993)"; my read of the plotted points —
mark as [VERIFY] against the original JVST B figure if a decimal point ends up load-bearing):

| F/Ar+ | Yield (SiO2/Ar+) |
|---|---|
| ~0 | 0.28 |
| 5 | 0.51 |
| 10 | 0.54 |
| 20 | 0.63 |
| 40 | 0.78 |
| 50 | 0.84 |
| 70 | 0.89 |
| 100 | 1.05 |
| 200 | 1.07 |
| plateau (model) | ~1.10 |

**The crossover, stated three ways:**
- **Half-rise** (Y = (0.28+1.10)/2 = 0.69) at **F/Ar+ ~= 25-30**.
- **90% of plateau** at **F/Ar+ ~= 100**.
- **Full dynamic range, zero-neutral to saturated: 0.28 -> 1.10 = 3.9x.**

Below `F/Ar+ ~ 25` the surface is strongly neutral-limited; above ~100 it is ion-limited (each ion
sees an equivalently passivated surface).

### 3.2 SiO2 + CF2 + F + Ar+ : the three-beam surface plot (Butterbaugh 1991)

J. W. Butterbaugh, D. C. Gray, H. H. Sawin, "Plasma-surface interactions in fluorocarbon etching of
silicon dioxide," *J. Vac. Sci. Technol. B* **9**, 1461-1470 (1991), DOI 10.1116/1.585451.
Reproduced as **Kwon Fig. 2.6, p. 36**: "Oxide etching yield measured with three beams (CF2, F and
Ar+)". Axes: `Etching Yield (SiO2/Ar+)` 0-0.7; `Flux Ratio (F/Ar+)` 0-160; `Flux Ratio (CF2/Ar+)`
0-20. The `F/Ar+ = 0, CF2/Ar+ = 0` corner is labelled "**Ar+ physical sputtering**" at Y ~ 0.1;
the high-`F/Ar+` edge plateaus at Y ~ 0.6-0.65; a ridge labelled "**competitive reaction effect**"
runs along the CF2 axis.

Kwon's verbatim description of the three regimes (**p. 35**):

> "Butterbaugh et al.^2 suggested a surface plot representing oxide etching yield as a function of
> incoming fluxes, where F, CF2 neutrals and Ar+ ion beams are used to mimic a real oxide etching
> process. (Figure 2.6)
> In this plot, there are three regimes. **In the regime where F flux is low, the surface kinetics is
> dominated by CF2 flux because the surface is covered with CF2. In the regime where F flux is
> moderate, the surface is covered with both CF2 and F. F is considered as the primary etchant and
> the etching yield is relatively independent of CF2 flux, which is due to the CF2 reduction by F
> flux. In the regime where F flux is high, the CF2 reduction by F is very rapid and surface is
> primarily covered with F.**
> This picture, however, is a case when the effect of deposition is not significant. It has been
> reported by many researchers that when carbon fluxes (CFx) increase, deposition component becomes
> greater^{3,4}."

and **p. 36**:

> "It is a modification of Butterbaugh's surface plot, where **there is a reduction in etching yield
> as the CFx flux increases due to the deposition on the surface.**"

**Sign structure to reproduce:** yield rises with F/Ar+ and *falls* with excess CFx/Ar+. Two
opposite-signed neutral responses. A model that is flat in radicals is wrong in both directions at
once.

### 3.3 Si + Cl + Ar+ : the independent crossover (Chang 1997)

J. P. Chang, J. Arnold, G. Zau, H.-S. Shin, H. H. Sawin, "Kinetic study of low energy argon
ion-enhanced plasma etching of polysilicon with atomic/molecular chlorine," *J. Vac. Sci. Technol.
A* **15**, 1853 (1997). Reproduced as **Huard Fig. 1.10, p. 36** (`huard_chad_phd_thesis.txt`
L177-179, L1449-1451; PDF page 54):

> "Fig. 1.10 Etch yield per Ar+ ion as a function of `Gamma_n/Gamma_i` (Cl/Ar+). The etching yield
> increases rapidly for low `Gamma_n/Gamma_i` and **saturates when the reaction becomes ion starved**.
> Reproduced from Chang et al.[22]"

Read from the figure: `Etching Yield (Si/Ar+)` vs `Flux Ratio (Cl/Ar+)`, x 0-800, three ion
energies.

| Ar+ energy | plateau yield | approx. half-rise Cl/Ar+ | saturated by |
|---|---|---|---|
| 100 eV | ~3.4 | ~50-70 | ~400 |
| 60 eV | ~2.0 | ~100 | ~600 |
| 35 eV | ~0.8+ (still rising at 800) | ~250 | not reached in range |

Huard's text statement (L1063-1074):

> "One significant consequence of relying on passivation reactions to increase sputtering yield is
> the synergistic coupling of ion and radical fluxes. **Chang et al. explored this coupling using
> molecular beam experiments and found that the etching yield depends strongly on the neutral to ion
> flux ratio, `Gamma_n/Gamma_i`, for low and moderate values of `Gamma_n/Gamma_i` (< 100)**, as shown
> in Fig. 1.10.[22] **The dependence saturates at larger values of `Gamma_n/Gamma_i` when the etching
> reaction becomes ion starved. In this ion starved regime the surface can be considered to be
> completely passivated, presenting each impinging ion with a similarly passivated surface and
> creating a constant etching yield. For lower values of `Gamma_n/Gamma_i` the surface is not fully
> passivated, and the coverage of passivated sites depends directly on ion flux, ion energy and
> radical flux.**"

**Inference (not verbatim, flagged):** the half-rise ratio moves *right* as ion energy falls in the
Chang data, and the Gray SiO2/F half-rise (~27 at 350 eV) is far left of the Chang Si/Cl half-rise
(~50-70 at 100 eV). The Langmuir crossover is `Gamma_n/Gamma_i = Y_removal/s_sticking` — it moves
with the per-ion neutral *consumption*. At HARC keV energies the per-ion consumption is larger
(uncapped `sqrt(E)` scaling would give ~3x from 350 eV to 3.4 keV), so the crossover ratio for
HARC oxide should be **larger than 100**, not smaller. Marked **[VERIFY]** — no source states this
directly, and the MCFPM per-particle cap complicates the scaling.

### 3.4 Reactor-side N/I, for scale (Kwon Figs. 2.10 and 2.20)

Kwon Fig. 2.10 (p. 40), C2F6 ICP, 5 mTorr, 300 eV DC bias: overall neutral-to-ion flux ratio
falls from **~230 at 200 W to ~75 at 500 W** as ion current rises 11 -> 33 microA.
Kwon Fig. 2.20 (p. 52), C4F8 + 80% Ar, 5 mTorr, 300 eV: **~165 at 300 W to ~70 at 500 W**.

**Real oxide-etch reactors deliver N/I ~ 70-230 at the wafer.** Consistent with Huang's 30, Krueger's
33-330, and petch's 22-32. That is the *open field*, and it sits at or above the measured saturation
knee — which is exactly why blanket etch-rate data is ion-limited and tells you nothing about the
floor.

### 3.5 The trap, stated by a source

Huard L4306-4310:

> "**The available experimental data indicating that large values of neutral-to-ion flux ratio
> produces a saturation in etch rate were predominantly obtained from measurements in the open field
> (non-patterned wafers).** Although these data have provided extremely important insights, **the data
> do not address the possible coupling of neutral and ion fluxes that may occur within features.**"

---

## 4. ARDE MECHANISM ATTRIBUTION IN MEASURED / MODELLED HARC

### 4.1 Coburn & Winters 1989 — the conductance baseline

J. W. Coburn, H. F. Winters, "Conductance considerations in the reactive ion etching of high aspect
ratio features," *Appl. Phys. Lett.* **55**, 2730-2732 (1989), DOI 10.1063/1.101937. Abstract
(publisher): vacuum conductance is adequate for *products* to leave the feature without gas-phase
collisions becoming important, **but "the conductance can be expected to limit the flow of the
reactive species to the bottom of the feature where the etching is taking place, thus creating the
possibility of an etch rate dependence on the aspect ratio of the etched feature."**

Their result as restated by Huard (L3731-3759):

> "Coburn and Winters, for example, determined that if the etch front consumes incident radicals with
> a reactive sticking coefficient, `Sn`, well known vacuum conductance concepts can be used to provide
> insights to ARDE. They proposed that the ratio of neutral flux arriving at the etch front,
> `Gamma_f`, compared to the incoming flux from the bulk plasma, `Gamma_in`, could be approximated as
> **`Gamma_f / Gamma_in = K / (K + Sn - K*Sn)`** (4.1)
> where `K` is the aspect ratio dependent probability that a randomly directed neutral particle
> incident on the feature opening will reach the etch front. The parameter `K` is analogous to
> Clausing's transmission probability... The assumptions required for this expression to be valid are
> that the side walls do not consume neutrals, and that the reactive sticking coefficient on the etch
> front, `Sn`, is constant in time and aspect ratio."

petch's `ARDE_PHYSICS_REFERENCE.md` already validates the transport side of this (`K`) three ways
against an independent ray-trace and an independent particle MC. **The transport is not the open
question; `Sn` is.**

### 4.2 Gottscho, Jurgensen & Vitkavage 1992 — neutral-ion synergy (THE fix mechanism)

R. A. Gottscho, C. W. Jurgensen, D. J. Vitkavage, "Microscopic uniformity in plasma etching,"
*J. Vac. Sci. Technol. B* **10**, 2133-2147 (1992), DOI 10.1116/1.586180. Their abstract states
that although ARDE and pattern-dependent etching have been observed across many material systems,
"the fundamental causes underlying these effects are poorly understood."

Their model as restated by Huard (L3760-3778) — **this is the load-bearing paragraph of the whole
document**:

> "Gottscho et al. expanded on the model proposed by Coburn and Winters by introducing the concept of
> **neutral-ion synergy**.[104] The basic premise of this model is that `Sn` is not constant, but
> rather simultaneously depends on the neutral and ion fluxes at the etch front. **In this model, each
> surface site can only be passivated once, and will then require an ion impact event before the
> underlying site can accept a new passivating radical. Therefore, high neutral fluxes will
> progressively produce smaller values of `Sn` since the available sites may already be occupied by
> reaction with earlier fluxes. Neutral starved regimes (that is, a low neutral radical flux relative
> to the ion flux) will have larger values of `Sn` since there is higher likelihood that sites will be
> empty. Since the neutral flux reaching the etch front depends on AR, then `Sn` also depends on AR.**
> Gottscho et al. developed analytical expressions for two extreme cases -- **perfectly diffusive
> walls (molecular flow), and perfectly absorbing walls (neutral shadowing)**. Most actual etching
> processes fall between these two extremes."

And the coupling it creates (Huard L4336-4341):

> "The synergy model of Gottscho et al. predicts that there is a coupling between the ion and neutral
> fluxes such that **changing the incoming ion flux can change the neutral flux to surfaces deep in the
> feature with no change in the neutral flux entering the feature**.[104] This process occurs due to
> the ion flux changing the steady state surface chlorine coverage on the etch front, which changes
> `Sn`. The change in `Sn` then impacts the neutral flux through Eq. 4.1."

### 4.3 Huard 2019 (MCFPM, Kushner group) — quantified attribution and an N/I sweep

C. M. Huard, PhD thesis, U. Michigan (`tmp/pdfs/huard_chad_phd_thesis.txt`); the underlying paper is
C. M. Huard, Y. Zhang, S. Sriraman, A. Paterson, M. J. Kushner, "Role of neutral transport in aspect
ratio dependent plasma etching of three-dimensional features," *J. Vac. Sci. Technol. A* **35**,
05C301 (2017). Ar/Cl2 = 80/20 ICP, 20 mTorr, 150 W, dc bias -113 V, IEAD bimodal at 125/325 eV with
2.2 deg HWHM; wafer fluxes Cl = 2.2e17, Cl+ = 1.0e16, Cl2+ = 8.3e15, Ar+ = 6.1e14 cm^-2 s^-1
(L3842-3845) — **base-case N/I = 11.6**, deliberately reduced from HPEM values.

**Attribution (L3813-3821, L4301-4305):**

> "We found that for the conditions investigated in this chapter, **the dominant cause of ARDE is the
> depletion of neutral species reaching the etch front relative to ions as the AR increases due to
> neutral transport issues.** We also found that increasing the neutral flux reaching the etch front
> relative to the ion flux can alleviate ARDE for small to moderate AR ( < 8 ), but at the expense of
> producing more tapered features and sometimes increasing ARDE at larger AR."

> "The results of this computational investigation suggest that **a dominant cause of ARDE is the
> decrease in neutral radical flux reaching the etch front with increasing AR, provided the process is
> not already in a neutral saturated regime.**"

**Ion vs neutral flux behaviour with AR (L3863-3880):**

> "The power density delivered to the horizontal surfaces by energetic ions ... is nearly uniform
> along the length of the trenches, as well as in the open field. This lack of sensitivity to AR
> results from the ion angular distribution being sufficiently anisotropic ... On the other hand, the
> flux of chlorine radicals incident onto surfaces ... decreases from the top of the feature to the
> bottom... **The relative insensitivity of power density to the etch front as a function of AR,
> compared to the strong dependence of neutral flux, indicates that neutral transport likely dominates
> the ARDE process for this mechanism.**"

**The N/I sweep — the decisive quantitative result (L4014-4057):**

> "the magnitude of the incident ion flux (`Gamma_i`) was varied while keeping the shape of the IEAD
> and the incident neutral flux (`Gamma_n`) constant. ... the initial etch rates for the first 6 s in
> the 50 nm trench are 1.66, 1.04, 0.58 and 0.41 times that of the base case for
> `Gamma_n/Gamma_i` = 5, 10, 20 and 30. **This trend indicates that the etch rate in the open field is
> dominantly ion starved -- there is sufficient Cl radical flux to nearly fully passivate sites and so
> the etch rate increases nearly linearly with increasing ion flux.** As discussed previously, for a
> given ion flux ARDE depends on the arrival of neutral flux at the etch site. **As a result, ARDE
> behaves as though the process is neutral starved**, higher ratios of `Gamma_n/Gamma_i` reduce the
> dependence of etch rate on etch depth."

> "**For AR=10, the normalized etch rate monotonically increases by 46%, from 0.50 to 0.73, for
> `Gamma_n/Gamma_i` increasing from 5 to 30. At AR = 6 ... The normalized etch rate increases by 66%,
> from 0.70 to 1.16 times the etch rate in the open field, for `Gamma_n/Gamma_i` increasing from 5 to
> 30.**"

**The same model is ion-limited in the open field and neutral-limited at AR 10.** That is the exact
duality petch must reproduce and currently does not.

### 4.4 Vanraes / ASML 2023 — experiment + multiscale model, charging excluded

P. Vanraes, S. Parayil Venugopalan, M. Besemer, A. Bogaerts, "Assessing neutral transport mechanisms
in aspect ratio dependent etching by means of experiments and multiscale plasma modeling," *Plasma
Sources Sci. Technol.* **32**, 064004 (2023), DOI 10.1088/1361-6595/acdc4f. Verbatim abstract
excerpt:

> "we accordingly assessed the neutral transport mechanisms in ARDE by means of experiments and
> multiscale modeling for SiO2 etching with CHF3/Ar and CF4/Ar plasmas. The experiments revealed a
> **local maximum in the etch rate for an aspect ratio around unity**, i.e. the simultaneous occurrence
> of regular and inverse reactive ion etching lag for a given etch condition. **We were able to
> reproduce this ARDE trend in the simulations without taking into account charging effects and the
> polymer layer thickness, suggesting shadowing and diffuse reflection of neutrals as the primary
> underlying mechanisms.** ... this work supports the insight that **physisorption may be more important
> in plasma etching at room temperature than originally thought**."

**Attribution: neutral shadowing + diffuse reflection; charging explicitly NOT needed.**

### 4.5 Lam Research 2023 — industry HARC data and attribution

M. Shen, T. Lill, et al. (Lam Research), *Jpn. J. Appl. Phys.* **62**, SI0801 (2023) — progress
review; local copy `tmp/pdfs/lam_shen_lill_jjap2023.txt`.

Mechanism statement (L179-199):

> "This etching mechanism encounters challenges for HARs because **the flux of neutral precursors for
> polymer formation attenuates rapidly as the aspect ratio increases.**^{13,14} ... The transport of
> neutrals inside the etch feature can be described by Knudsen transport. **For a cylinder with an
> aspect ratio of 50:1, only 2.5% of the incoming flux reaches the other end of the cylinder due to
> diffuse collisions with the sidewall.**^{18} The reason for this loss of neutrals is that they don't
> have a preferred direction of movement. Ions are accelerated in the normal direction to the wafer
> surface when they pass the plasma sheath. **Therefore, the etching mechanism for aspect ratios
> beyond 40 or 50 is dominated by fast, chemically active species which carry the mechanical energy
> and possess the chemical composition to initiate etching. These species can be ions or neutrals,
> which originate from ions that were neutralized in collisions with the feature sidewall. The
> mechanisms can be described as more akin to chemical sputtering rather than synergistic neutral-ion
> etching.**"

> "Attenuation of ion and neutral fluxes as a function of aspect ratio reduces the etching rate. This
> undesired effect is called aspect ratio-dependent etching (ARDE). **To reduce ARDE, the transport of
> neutrals and ions must be enhanced.** For ions, this can be achieved by increasing the ion energy."
> (L200-204)

**Measured decay curve** (Fig. 3, L218-226): normalized etch rate vs depth for 100 nm circular holes,
conventional vs lean low-temperature process:

> "**The etching rate diminishes much more quickly for the high-temperature process than for the
> low-temperature process.** One possible explanation is that chemical reactions are suppressed at
> lower temperatures. **Coburn and Winters predicted a lower ARDE when the reactive sticking
> coefficient is reduced.**^{19} An additional explanation for the improved ARDE performance at lower
> temperatures is potentially surface diffusion of neutrals.^{20}"

Ref 19 is Coburn & Winters APL 55, 2730 (1989). **Lam attributes the measured ARDE difference between
two production HARC recipes to the reactive sticking coefficient of the neutrals** — i.e. to `Sn` in
Eq. 4.1, exactly the Gottscho variable.

Simulated ARDE (Fig. 6a, L291-309):

> "The ONON channel hole layer ... **exhibits a flat etch rate up to an aspect ratio of ~60 and
> decreases at higher aspect ratios.** ... The flat etch rate at lower aspect ratios can be partially
> attributed to very low IADs (<0.5 deg 1 sigma). **Similar curves have been reported by Huang et
> al.**^{13} In Huang's SiOx feature-scale model, the SiOx etch rate at relatively lower aspect ratios
> reaches a maximum before decreasing."

Also relevant to petch's mixed layer at the floor: Fig. 2 XPS shows the lean low-T process leaves
**"the absence of a measurable polymer film on the surface of SiO2 for low-temperature HAR etching
processes versus conventional processes"** (L209-212). In modern HARC the etch front is a bare/mixed
layer, not a polymer-blocked one.

Companion transport paper: T. Panagopoulos, T. Lill, "Neutral transport during etching of high
aspect ratio features," *J. Vac. Sci. Technol. A* **41**, 033006 (2023), DOI 10.1116/6.0002468 —
**"For an aspect ratio of depth to diameter of 100:1, the flux at the bottom of the feature is only
1.3% of the incoming flux."** [VERIFY — from abstract via search; full text paywalled.]

### 4.6 Joubert & Oehrlein 1994 — RIE lag correlates with the *chemical* term

O. Joubert, G. S. Oehrlein, Y. Zhang, "Fluorocarbon high density plasma. V. Influence of aspect ratio
on the etch rate of silicon dioxide in an electron cyclotron resonance plasma," *J. Vac. Sci.
Technol. A* **12**, 658-664 (1994). Publisher summary: **"The magnitude of the RIE lag is correlated
with the deposition rate of fluorocarbon film on an unbiased sample, showing that chemical effects
are important to understand the mechanisms of RIE lag in high density plasmas."** [VERIFY — abstract
only; full text paywalled.]

### 4.7 Attribution scoreboard

| source | system | measured / modelled decay | stated dominant mechanism |
|---|---|---|---|
| Coburn & Winters 1989 (APL 55, 2730) | general | analytic `K/(K+Sn-K*Sn)` | neutral conductance; products are NOT limiting |
| Gottscho 1992 (JVST B 10, 2133) | review | analytic, two wall limits | neutral-ion synergy; `Sn` is flux-ratio dependent |
| Joubert/Oehrlein 1994 (JVST A 12, 658) | SiO2, ECR FC | RIE lag vs feature size | **chemical** (FC deposition rate), not pure geometry |
| Huard/Kushner 2017-19 (JVST A 35, 05C301) | Si, Ar/Cl2 | AR10 rate 0.50-0.73 of open field over N/I 5-30 | **neutral depletion relative to ions**; ion power density ~AR-flat |
| Huang (thesis ch. 6) | SiO2, Ar/C4F8/O2 CCP | **-80% by AR 40** | **both**: power to the etch front (dominant) + neutral conductance |
| Vanraes/ASML 2023 (PSST 32, 064004) | SiO2, CHF3/Ar, CF4/Ar | local max at AR~1, regular+inverse lag | **shadowing + diffuse reflection of neutrals**; charging & polymer NOT needed |
| Lam 2023 (JJAP 62, SI0801) | ONON HARC | flat to AR ~60 then decay; conventional decays much faster than lean low-T | **neutral sticking coefficient / Knudsen transport** (cites Coburn-Winters); IAD sets the flat region |
| Panagopoulos & Lill 2023 (JVST A 41, 033006) | cylinders | AR 100:1 -> 1.3% transmission | Knudsen conductance |

**Charging is attributed a role in *profile distortion* (notching, twisting, bowing) across the
literature but is NOT the attributed cause of the rate decay** in any of the above; Vanraes
explicitly reproduces the ARDE trend with charging switched off.

---

## 5. VERDICT — the doctrine-compliant fix shape

The published, citable form that makes a mixed layer neutral-responsive **without a fitted knob** is
the **translating mixed-layer (TML) / Gottscho-synergy** formulation. Every element below has a named
source; none introduces a free constant.

### 5.1 The reference implementation, verbatim (Kwon TML, SiO2 + F)

Kwon thesis §3.3 (pp. 72-75), the model that reproduces Gray's beam data in Fig. 3.4:

**Assumptions (p. 65):**
> "1. Within the translating mixed-layer, the composition is homogeneous; i.e., well mixed by ion
> bombardment so that there is no concentration gradient.
> 2. **The total number of atoms in the translating mixed-layer is conserved.** The difference between
> the adsorption flux ... and the removal flux ... equals the movement flux (the flux from the
> substrate volume to the surface translating mixed-layer)...
> 3. All the reactions, which include ion bombardment, adsorption reaction, ion enhanced chemical
> etching surface reaction, physical sputtering, and deposition, occur in the translating mixed-layer
> volume."

Layer thickness (p. 64): *"The thickness of the reaction volume can be as much as about 10 monolayers
or 30-40 Angstrom, depending on experimental conditions."*

**Adsorption — site-limited with explicit stoichiometric capacity (p. 72):**
```
F(g) -> F(s)     r1 = s1 * R_F * ( x_Si - (1/2) * x_F * 2*x_Si/(2*x_Si + x_O) )
```
> "**assuming there are two sites available per silicon for fluorine adsorption.** It is also assumed
> that the chemical affinities of fluorine to oxygen and to silicon are equal. The term in the
> parenthesis `x_F * 2*x_Si/(2*x_Si + x_O)` is the concentration of fluorine in the translating
> mixed-layer incorporated with silicon."
```
F(g) -> F(s)     r2 = s2 * R_F * ( x_O - x_F * x_O/(2*x_Si + x_O) )
```
> "because **number of available adsorption sites per oxygen is assumed to be one.**"
```
F+(g) -> F(s)    r3 = s3 * R_F+     "where the sticking coefficient is unity"   (implantation-like)
```

**Removal — LINEAR in layer halogen concentration (p. 73):**
```
Si(s) + 2F(s) -> SiF2(g)   r4 = beta1 * x_F * 2*x_Si/(2*x_Si + x_O)
O(s)  +  F(s) -> OF(g)     r5 = beta2 * x_F *   x_O/(2*x_Si + x_O)
```
> "which means **the production of SiF2 is proportional to the concentration of fluorine incorporated
> with silicon in the translating mixed-layer.**"

**Physical sputtering — mole-fraction weighted, per species (p. 74):**
```
Si(s) -> Si(g)   r6 = Y_Si * x_Si
O(s)  -> O(g)    r7 = Y_O  * x_O
F(s)  -> F(g)    r8 = Y_F  * x_F
```

**Closure — the etch rate IS the conservation residual (p. 74):**
```
Si(sub) + 2 O(sub) -> Si(s) + 2 O(s)
r9 = (1/3) * ( -r1 - r2 - r3 + 3*r4 + 2*r5 + r6 + r7 + r8 )
```
> "and **the movement flux is essentially the etching yield because it is the net flux removed.**"

**Stoichiometric anchor stated on p. 76:** `SiO2(s) + 4 F(s) -> SiF4(g) + O2(g)`, though the TML
splits Si and O removal into separate reactions (r4/r5).

**Result:** the Fig. 3.4 curve — sputter floor 0.28, half-rise at F/Ar+ ~ 27, plateau 1.1.

### 5.2 The seven required elements, each with its source

| # | element | why it produces neutral limitation | source |
|---|---|---|---|
| **E1** | **Site-limited adsorption with a stated stoichiometric capacity**, `r_ads = s * Gamma_n * (available sites)` | at high `Gamma_n` the bracket -> 0, so `Sn` self-limits; at low `Gamma_n` sites are empty and every arrival counts | Kwon r1/r2 (2 F per Si, 1 F per O); Huang "further passivation of complex" p0 = 1e-4 (L10214-10222); Krueger 2e-4 (L6556-6564) |
| **E2** | **Removal linear in layer halogen/complex concentration**, `r_etch = beta * x_F`, with `beta ∝ Gamma_i * Y(E, theta)` | rate follows coverage, which follows `s*Gamma_n/(s*Gamma_n + Y*Gamma_i)` | Kwon r4/r5, p. 73 |
| **E3** | **Finite site turnover: removal regenerates a pristine site** | each removal creates a fresh radical demand; the demand rate is set by the ion flux, so the *ratio* is what matters | Gottscho synergy via Huard L3763-3770: *"each surface site can only be passivated once, and will then require an ion impact event before the underlying site can accept a new passivating radical"*; Huang Eq. 3.30 `SiO2CFx(s) + I+ -> SiFx + CO2` leaves `SiO2(s)` |
| **E4** | **Stoichiometric enforcement of F (or C) per removal** — 4 F per SiF4, or Huang's 1 C + 1-3 F per SiO2 formula unit | closes the ledger; makes the F budget, not a rate constant, the limiter at depth | Kwon p. 76 `SiO2 + 4F -> SiF4 + O2`; Huang App. E rows (`SiO2CF -> SiF + CO2` etc., L10229-10236) |
| **E5** | **One non-reflective event per incident particle** (cap, not scale) | prevents keV energy from papering over an unpassivated surface | Huang L2311-2319 (already in petch, `_unit_cumulative_scale`) |
| **E6** | **Physical-sputter channel weighted by the UNREACTED fraction, with the off-normal-peaked angular class** | supplies the `(1 - theta)` floor and makes it *small at normal incidence*, which is the HARC floor condition | Kwon r6 = `Y_Si * x_Si`; Huang L2290-2296 (phys. sputter peaks at 60 deg, reduced at normal; chem. sputter unity to 45 deg); Krueger Appendix B angular-class column (1 = Kress 1999, 2 = Chang & Sawin 1997) |
| **E7** | **Sticking coefficient fed back into transport**, `Gamma_f/Gamma_in = K/(K + Sn - K*Sn)` with `Sn = Sn(Gamma_n/Gamma_i)` | closes the loop Gottscho identified: changing the *ion* flux changes the *neutral* flux at depth | Coburn & Winters APL 55, 2730 via Huard Eq. 4.1 (L3734-3745); Gottscho JVST B 10, 2133 via Huard L3760-3778 |

**Plus one delivery-side element that is not surface kinetics at all but sets the floor's N/I:**

| **E8** | **Spent energetic FC species must thermalize back into the RADICAL ledger, subject to feature conductance** | above AR 10, **>95%** of the floor's radicals are ion-derived, not conducted | Huang L5716-5727, and the 15%->60% CFx+ sweep giving ~3x floor radical flux and **-30% etch rate** (L5742-5754) |

### 5.3 What is NOT allowed under the knob-elimination doctrine

- No `neutral_response_exponent`, no `arde_gain`, no `depth_scaled_sticking`. Every constant above is
  either (a) a stoichiometric integer (2 sites/Si, 1 site/O, 4 F/SiF4, 1 C/SiO2), (b) a published
  probability from Krueger Appendix B / Huang Appendix E, (c) a measured layer thickness
  (10 monolayers / 30-40 Angstrom, Kwon p. 64), or (d) a conservation identity (r9).
- The one genuinely open convention is the **angular-class normalisation** (§2.5). That is a
  *reading* of a published table, resolvable by reproducing the beam data, not a free parameter.

### 5.4 Falsifiable gate to preregister (this is the actionable output)

A 0-D beam test with **published data on both axes**, no feature, no transport:

> **GATE N1 (SiO2 / F / Ar+, 350 eV).** Drive the mixed layer with an Ar+ beam at 350 eV and an F
> beam, sweeping `F/Ar+` over `0 -> 500`. Required, from Gray 1993 as plotted in Kwon Fig. 3.4:
> - `Y(F/Ar+ -> 0) / Y(saturated) = 0.28/1.10 = 0.25 +/- 0.05`
> - half-rise at `F/Ar+ = 27 +/- 8`
> - `Y(F/Ar+ = 100) >= 0.90 * Y(plateau)`
> - monotone, no overshoot

> **GATE N2 (deposition branch).** With a CF2 beam added at fixed `F/Ar+`, the yield must *decrease*
> with `CF2/Ar+` over `0 -> 20` (Butterbaugh three-beam plot via Kwon Figs. 2.6/2.7; Kwon p. 36:
> *"there is a reduction in etching yield as the CFx flux increases due to the deposition on the
> surface"*).

> **GATE N3 (HARC floor ratio).** Instrument the AR 10-20 floor and report `radical / (ion + hot
> neutral)`. Published target from Huang L5470-5478: **0.33 - 0.5 at AR > 10**, ~0.64 at AR 40,
> against 30 at the wafer. Any value within 30x of the wafer-plane ratio means the delivery, not the
> kinetics, is the first defect.

### 5.5 Direct answer to "is <1% over 50x physically possible?"

**No.** Take the parent's floor at either reading:
- If the floor sits at `N/I ~ 144`, a 50x radical reduction lands at `N/I = 2.9`. On the Gray 1993
  curve that is `Y: 1.06 -> ~0.44`, i.e. **-58%**.
- If the floor sits at `N/I ~ 100`, 50x down lands at `N/I = 2.0`, `Y: 1.05 -> ~0.42`, **-60%**.
- If the floor is *already* radical-starved at `energetic:radical = 100-144` (i.e. `N/I ~ 0.007-0.01`),
  then the layer should be sitting essentially **at the sputter floor**, `Y ~ 0.28/1.10 = 25%` of the
  wafer-plane yield, and the etch rate should already show a ~4x ARDE from surface kinetics alone
  before any transport effect is counted.

There is no window in the measured data in which SiO2 is flat to <1% over a 50x radical sweep. The
flatness is a model defect, and E1/E2/E3/E6 above are where it lives.

---

## 6. Source list

**Primary, local (`tmp/pdfs/`):**
1. P. Krueger, PhD thesis, U. Michigan (Kushner group), `krueger_thesis.txt`. Key: L1547-1571
   (regime taxonomy), L2519-2521 (ion flux band), L4015-4019 (mechanism requires radicals),
   L4342-4349 (Table 6.1 wafer fluxes), L4364 (IEAD to 4800 eV), L4384-4392 (mechanism overview),
   **L4879-4889 (neutral-transport-limited + "ion energy might be overestimated")**, L5904-5941
   (bare sputter rows), L6006-6012 (passivation rows), L6180+ (complex removal rows),
   L6548-6564 (fluorination / further passivation).
2. C. Huang, PhD thesis, U. Michigan (Kushner group), `huang_thesis.txt`. Key: L2277-2320 (yield law,
   angular classes, unity-cumulative normalisation), L3235-3295 (mechanism), **L5391-5403 (ion vs
   neutral transport asymmetry)**, **L5430-5478 (floor flux ladder, "factor of 2-3")**, L5573-5605
   (ARDE attribution), L5662-5664 (neutral-starved at high power), **L5716-5754 (>95% ion-derived
   floor radicals; +3x radicals -> -30% rate)**, L5980-6012 (chapter summary), L10140-10385
   (Appendix E mechanism table).
3. C. M. Huard, PhD thesis, U. Michigan (Kushner group), `huard_chad_phd_thesis.txt`. Key:
   **L1063-1074 (Chang crossover statement)**, L177-179/L1449-1451 + PDF p. 54 (Fig. 1.10 Chang
   data), L1100-1128 (Standaert polymer-thickness mediation), **L3731-3778 (Coburn-Winters Eq. 4.1 +
   Gottscho synergy)**, L3804-3898 (attribution), **L4006-4118 (N/I sweep, AR10 0.50->0.73)**,
   L4300-4349 (parameters affecting ARDE).
4. M. Shen, T. Lill et al. (Lam Research), *Jpn. J. Appl. Phys.* 62, SI0801 (2023),
   `lam_shen_lill_jjap2023.txt`. Key: L179-226 (Knudsen 50:1 -> 2.5%; sticking-coefficient
   attribution; Fig. 3 decay), L291-344 (Fig. 6a flat to AR 60; Fig. 7 selectivity -50% AR 40->140).

**Primary, fetched:**
5. O. Kwon, "Surface kinetics modeling of silicon oxide etching in fluorocarbon plasmas," ScD thesis,
   MIT DMSE, 2004, advisor H. H. Sawin. https://hdl.handle.net/1721.1/28353
   (local copy: scratchpad `kwon_thesis.pdf`, PDF pages = printed page + 1). Key: **p. 76 Fig. 3.4
   (SiO2 yield vs F/Ar+, Gray 1993 data, 350 eV)**, p. 36 Fig. 2.6 (Butterbaugh three-beam),
   p. 35 (three-regime text), pp. 64-65 (TML assumptions), **pp. 72-75 (SiO2/F TML rate laws)**,
   p. 40 Fig. 2.10 and p. 52 Fig. 2.20 (reactor N/I 70-230).

**Cited, not fetched in full (marked where a number depends on them):**
6. J. W. Butterbaugh, D. C. Gray, H. H. Sawin, *J. Vac. Sci. Technol. B* **9**, 1461 (1991).
   DOI 10.1116/1.585451.
7. D. C. Gray, I. Tepermeister, H. H. Sawin, *J. Vac. Sci. Technol. B* **11**, 1243 (1993).
8. J. P. Chang, J. Arnold, G. Zau, H.-S. Shin, H. H. Sawin, *J. Vac. Sci. Technol. A* **15**, 1853
   (1997); and J. P. Chang, H. H. Sawin, *J. Vac. Sci. Technol. A* **15**, 610 (1997).
9. O. Kwon, H. H. Sawin, "Surface kinetics modeling of silicon and silicon oxide plasma etching.
   I / III," *J. Vac. Sci. Technol. A* **24** (2006); Part III at 24(5), 1920-1927.
10. J. W. Coburn, H. F. Winters, *Appl. Phys. Lett.* **55**, 2730 (1989). DOI 10.1063/1.101937.
11. R. A. Gottscho, C. W. Jurgensen, D. J. Vitkavage, *J. Vac. Sci. Technol. B* **10**, 2133 (1992).
    DOI 10.1116/1.586180.
12. O. Joubert, G. S. Oehrlein, Y. Zhang, *J. Vac. Sci. Technol. A* **12**, 658 (1994).
13. C. M. Huard, Y. Zhang, S. Sriraman, A. Paterson, M. J. Kushner, *J. Vac. Sci. Technol. A* **35**,
    05C301 (2017).
14. P. Vanraes, S. Parayil Venugopalan, M. Besemer, A. Bogaerts, *Plasma Sources Sci. Technol.*
    **32**, 064004 (2023). DOI 10.1088/1361-6595/acdc4f.
15. T. Panagopoulos, T. Lill, *J. Vac. Sci. Technol. A* **41**, 033006 (2023).
    DOI 10.1116/6.0002468.
16. T. E. F. M. Standaert, C. Hedlund, E. A. Joseph, G. S. Oehrlein, T. J. Dalton, *J. Vac. Sci.
    Technol. A* **22**, 53 (2004) — steady-state FC film thickness mediates radical and ion-energy
    transport to the etch front.

**[VERIFY] items:** (a) the digitized Gray 1993 point values in §3.1 — read off Kwon Fig. 3.4, not
from the original JVST B figure; (b) the Chang half-rise ratios in §3.3 — read off Huard Fig. 1.10;
(c) the energy scaling of the crossover ratio in §3.3 — my inference, no source states it;
(d) Panagopoulos & Lill's 1.3% at AR 100 — from abstract only; (e) Joubert/Oehrlein's RIE-lag /
deposition-rate correlation — from abstract only.
