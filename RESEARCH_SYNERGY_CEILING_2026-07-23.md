# Ion–Neutral Synergy and the Supply Ceiling: Literature Verification for Retiring the Yield-Scale Knob

Date: 2026-07-23
Author: literature-research pass (petch / Krüger-2024 blind validation follow-up)
Status: **not committed** — decision-support memo.

## 0. The question on the table

Our blind Krüger-2024 transfer missed one held-out behavior: measured etch depth
**saturates sharply between 6 and 8 kW LF power** at fixed neutral fluxes while the
digitized fig-16b mean ion energy keeps rising (~+20%/step: 2551 → 2998 → 3593 eV
for 4/6/8 kW). petch predicted **614 / 723 / 875 nm**; MCFPM reference **635 / 715 / 720 nm**;
experiment shows "few differences" 4→8 kW. The 8 kW point is the miss.

Desk analysis: with the **published, unscaled** yields (Krüger, complex removal
p0 = 0.1384, Eth = 35 eV, Er = 140 eV, n = 1), ion **removal capacity**
(ion flux ~1.2e16 cm⁻² s⁻¹ × complex yield at per-power energies) crosses the neutral
**supply ceiling** (complex-formation rate from the fluorocarbon fluxes, ~5.3e16 cm⁻² s⁻¹)
right between 6 and 8 kW → the etch becomes neutral-limited and depth flattens. Our fitted
`oxide_etch_yield_scale = 0.5586` divided capacity down below the ceiling everywhere,
disabling the crossover, so the law ran unbounded and overshot at 8 kW.

This memo verifies each leg of that mechanism against literature before we retire the knob.

---

## 1. Ion–neutral synergy / flux-ratio-limited etching (the ceiling is canonical)

**1.1 The founding synergy experiment.**
J. W. Coburn & H. F. Winters, "Ion- and electron-assisted gas-surface chemistry — an
important effect in plasma etching," *J. Appl. Phys.* **50**, 3189–3196 (1979),
DOI: 10.1063/1.326355. The canonical Si + XeF₂ + Ar⁺ toggling experiment: XeF₂ (neutral F)
alone and 450 eV Ar⁺ alone each give small etch rates; **combined the rate jumps ~an order
of magnitude above their sum** — the defining evidence that etching is a *product* of a
neutral-supply channel and an ion-activation channel, so **whichever channel is scarce caps
the rate.** (Exact per-curve rate values from the original fig. are [VERIFY] against the
paper; the "≈ order-of-magnitude synergy / sum-exceeding" statement is firmly established and
reproduced in every review since.)

**1.2 The quantitative flux-ratio model.**
D. C. Gray, I. Tepermeister & H. H. Sawin, "Phenomenological modeling of ion-enhanced surface
kinetics in fluorine-based plasma etching," *J. Vac. Sci. Technol. B* **11**, 1243–1257 (1993),
DOI: 10.1116/1.586925. Independently varied F-atom and Ar⁺ fluxes over several orders of
magnitude on poly-Si and SiO₂. The etch rate is modeled as **ion-driven removal of a
reactive-layer coverage θ that is fed by neutral adsorption** — a Langmuir-type site balance.
Two regimes: **neutral(reactant)-limited** at low neutral/ion ratio (rate ∝ neutral flux) and
**ion-limited / saturated** at high neutral/ion ratio (rate ∝ ion flux, independent of extra
neutrals). The saturation of one regime *is* the ceiling. (Exact crossover flux-ratio values
and the θ(Γn/Γi) fit constants are [VERIFY] at the page/figure level; the two-regime structure
and Langmuir-coverage form are the paper's central result.)

**1.3 Beam confirmation at fixed coverage.**
J. P. Chang & H. H. Sawin, "Kinetic study of low-energy ion-enhanced polysilicon etching using
Cl, Cl₂, and Cl⁺ beam scattering," *J. Vac. Sci. Technol. A* **15**, 610–615 (1997),
DOI: 10.1116/1.580692. Directly measured: **the ion-enhanced yield rises with the Cl/Cl⁺
flux ratio then saturates once the surface is chlorine-saturated** — i.e. at fixed (saturated)
coverage the per-ion yield stops responding to more neutrals; the rate is then set by ion flux.
This is the microscopic statement of §1.2 and is exactly the coverage-capped behavior our
supply ceiling encodes. (This is also the source of petch's `chang_sawin_1997` angular model.)

**1.4 SiO₂-specific: rate flattening / high threshold vs bias.**
T. E. F. M. Standaert, G. S. Oehrlein and co-workers on fluorocarbon SiO₂ etching (e.g.
"SiO₂ and Si etching in fluorocarbon plasmas: a detailed surface model coupled with a complete
plasma and profile simulator," and the CHF₃/ICP selectivity studies) establish that oxide
etching proceeds through a **steady-state fluorocarbon/complex layer**, and that etch rate vs
average ion energy shows a **threshold then a compressive (sub-linear) response** — SiO₂ has a
higher energy threshold than the polymer, and the oxide rate is gated by how fast the FC layer
can be supplied and converted, not by ion energy alone. Representative modern quantitative
confirmation: dual-frequency CCP Ar/C₄F₈ HAR-SiO₂ study, *Materials* (MDPI) PMC10222222,
which separates ion-energy vs ion-flux contributions and shows the oxide rate saturating with
ion energy at fixed chemistry. (Specific Standaert JVST volume/page and the ellipsometry rate
numbers are [VERIFY]; the FC-layer-mediated, threshold-then-saturating picture is standard.)

**Takeaway (1):** "etch rate saturates when ion capacity exceeds neutral supply" is not a
modeling convenience — it is the *defining* structure of ion-enhanced etching from Coburn-Winters
onward, quantified by Gray-Tepermeister-Sawin as a Langmuir coverage balance and confirmed at
the beam level by Chang-Sawin. A feature-scale model **without** a supply ceiling is missing the
first-order physics of the process.

---

## 2. Does the Krüger/Huang-Kushner MCFPM carry an explicit supply ceiling? (Yes.)

Primary sources read directly:
- F. Krüger, *Modeling and Optimization of High Aspect Ratio Plasma Etching*, PhD thesis,
  Univ. of Michigan (2024), DOI: 10.7302/23106
  (cpseg.eecs.umich.edu/pub/theses/Krueger_Florian_PhD_Thesis_2024.pdf). §2.1.3, §2.2.1, §6.9.2.
- S. Huang, S. Shim, S. K. Nam, M. J. Kushner et al., "Plasma etching of high aspect ratio
  features in SiO₂ using Ar/C₄F₈/O₂ mixtures: A computational investigation," *J. Vac. Sci.
  Technol. A* **37**, 031304 (2019), DOI: 10.1116/1.5090904
  (cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf). The SiO₂ mechanism Krüger uses
  is "based on previous work by Huang et al." (thesis §6.4).

**2.1 Explicit site balance in the surface kinetics.** Krüger §2.1.3 (Surface Kinetics
Module / **Surface Site Balance Model, SSBM**): reaction rate

    R_im = α_i · φ_Am · θ_Bm            (Eq. 2.39)

where φ is the incident gas-phase flux and **θ_Bm is the fractional coverage of surface site B**,
integrated from all reaction rates (3rd-order Runge-Kutta) — an explicit, bounded (θ ≤ 1)
site-occupancy ceiling implemented as a *variable sticking coefficient*. When a site is consumed
faster than it is replenished, θ falls and the rate self-limits. This is a hard supply ceiling.

**2.2 Two-step complex chain at the feature scale (MCFPM).** Thesis §6.4 / Table 6.2: SiO₂ is
removed via a **two-step supply-gated chain** —
(i) unsaturated fluorocarbons chemisorb on SiO₂ to form a **SiO₂-CₓFy complex** ("significantly
lowering the binding energy," thesis §3, line ~3805), then
(ii) an ion / hot neutral chemically sputters that complex via the energetic yield law. The
enhanced-removal channel **cannot fire faster than the complex is formed** from the FC flux →
neutral-supply-limited when ion capacity is high. This is precisely the desk-analysis mechanism.

**2.3 The energetic yield law petch reproduces.** Thesis Eq. 2.40:

    p(E,θ) = p0 · [ (E − Eth)/(Er − Eth) ]^n · f(θ)

matching petch's `EnergeticYield`. With **n = 1** it is *unbounded* in E — nothing in this factor
saturates, so saturation in the model must come from the coverage/site balance (§2.1–2.2),
**not** from the yield factor.

**2.4 Krüger states the high-power saturation IS a neutral ceiling — and that his own model
under-captures it.** Thesis §6.9.2 (Variation of Low Frequency Power, Plf = 0/4/6/8 kW), verbatim:

> "The remaining cases (Plf = 4, 6 and 8 kW) have unclogged features and full etching with
> unexpected little variation as a function of LF power. In the experimental data a doubling of
> the Plf (4 kW to 8 kW) produces few differences in the final features. **These trends indicate
> that above a certain threshold energy the etch progression and the mask removal process are not
> ion starved, but rather limited by neutral gas transport.** To some degree this trend is
> reproduced by the simulations where etch depth does increase with increasing low frequency
> power however **the rate of increase is substantially sublinear. These outcomes indicate that
> the effect of ion energy (for example in sputter yield or related processes) might be
> overestimated in the mechanism.**"

Two decisive points:
1. Krüger attributes the 6→8 kW plateau to a **neutral (supply/transport) limit**, exactly the
   ceiling in the desk analysis.
2. He concedes his **own MCFPM only partially reproduces the plateau** (sublinear, not flat) and
   diagnoses the residual as **ion-energy/sputter-yield effect being overestimated**. petch's
   875 nm-at-8 kW overshoot is the *same* residual, magnified because our fitted yield-scale knob
   removed the ceiling that would otherwise have clamped it.

**2.5 The fig-16 reference depths.** Krüger's own MCFPM: 635 / 715 / 720 nm at 4 / 6 / 8 kW —
a clear 6→8 kW plateau (+0.7%) after a 4→6 kW rise (+12.6%). petch: 614 / 723 / 875 nm — matches
through 6 kW (within ~1%) then diverges upward at 8 kW. The divergence appears exactly where the
capacity/supply crossover is predicted.

**Takeaway (2):** the reference model has an explicit site-balance supply ceiling (SSBM Eq. 2.39)
and a two-step complex chain (Table 6.2); the author explicitly names the high-power saturation as
neutral-limited and explicitly says the ion-energy effect is over-weighted. Our residual is the
Krüger residual, amplified by the knob.

---

## 3. Yield linearity (n = 1) — defensible, and not where the fix belongs

**3.1 Low-energy law is sqrt(E), not linear.** C. Steinbrüchel, "Universal energy dependence of
physical and ion-enhanced chemical etch yields at low ion energy," *Appl. Phys. Lett.* **55**,
1960–1962 (1989), DOI: 10.1063/1.102336: physical and ion-enhanced chemical yields are
**linear in √E** — Y = A(√E − √Eth) — for metals, Si and SiO₂ by noble-gas and reactive ions,
across the sub-keV regime. petch already carries this as the separate `SteinbruchelYield`
(A·(√E − √Eth)·f(θ)).

**3.2 Krüger's n = 1 is a fitted network exponent, not a universal law.** In Eq. 2.40 n is a
**per-reaction tuning exponent** (physics-informed gradient descent, thesis §6.8), calibrated to a
single base case, not asserted as the true E-dependence. At the keV energies here (2.5–3.6 keV
mean, up to 4.8 keV) **neither √E nor linear is physically reliable** — both are low-energy fits.
The real high-energy behavior is a *compression*: keV ions deposit an increasing fraction of their
energy **below** the thin reactive/mixing layer where the complex lives, so the *useful* yield per
incremental eV falls (see companion memo `RESEARCH_ENERGY_DEPOSITION_ANCHORS_2026-07-23.md`:
ZBL/LSS stopping + Sigmund deposited-energy-in-layer). This is the "energy-deposition softening"
term.

**3.3 Beam evidence for per-ion yield saturation at fixed coverage.** Chang-Sawin 1997 (§1.3):
at chlorine saturation the yield stops responding — per-ion yield **saturates at fixed neutral
coverage**. Gray-Sawin 1993 (§1.2): in the ion-limited regime, adding neutrals does nothing.
Together these say the observed saturation is **coverage/energy-deposition**, *not* a steeper vs
shallower choice between √E and E in the bare yield prefactor.

**Takeaway (3):** keeping n = 1 is defensible *as a low-energy prefactor*, but it is the wrong knob
to carry the saturation. The saturation belongs in (a) the coverage/supply ceiling and (b) an
energy-deposition softening above a keV-scale knee — **not** in a global multiplicative yield
scale, and not in switching the exponent.

---

## 4. Grazing-ion reflection

**4.1 Standard feature-scale reflection treatment (what Krüger/MCFPM uses).** Thesis §2.2.1,
Eq. 2.41: post-collision energy

    Es(θ) = Ei · [ (Ei − Ec)/(Ets − Ec) ] · [ (θ − θc)/(90° − θc) ]

with a **specular branch** for Ei > Ets (energy preserved, exit angle = incidence — i.e. grazing
ions glance off and continue down the feature, redepositing energy deeper) and a **diffusive
branch** for Ei < Ec or θ < θc (fractional energy loss + randomized angle). This threshold
specular/diffuse split is the standard treatment in profile simulators (MCFPM, and the same idea
in the empirical Kress/Chang-Sawin closures petch carries).

**4.2 Angular yield = deposition-rise minus reflection-loss.** Empirical angular etch-yield curves
decompose into (i) a rise toward ~60–70° off-normal from increased energy deposition near the
surface and (ii) a collapse toward grazing from **ion reflection** carrying energy away:
- Chang-Sawin 1997 (DOI 10.1116/1.580692): measured poly-Si ion-enhanced **yield reduced ≈30% at
  60° and ≈50% at 70°** off-normal — the reflection-loss shoulder. (petch `chang_sawin_1997`.)
- J. D. Kress, D. E. Hanson, A. F. Voter et al. (MD of energetic halogen/Ar reflection off
  Si/SiO₂, ~1999), cited in petch as DOI 10.1116/1.581948, supplies the reflected-energy /
  reflection-probability angular dependence used for the `kress_1999` yield model.
  (Exact JVST volume/page/year for that DOI is [VERIFY].)

**4.3 Binary-collision / measurement picture at grazing.** BCA and sputter-yield data: yield
rises to a maximum near ~70–80° then falls toward zero at grazing as the **reflection coefficient
approaches unity**; near-100% reflection is predicted only for very glancing angles and is *not*
observed below ~300 eV (so reflection is energy-dependent, weaker at low E). Representative:
Y. Yamamura's angular-yield formalism and the sputter-yield-vs-angle reviews
(e.g. *Nucl. Instrum. Methods B* angular-dependence literature). (Specific citations [VERIFY];
the qualitative rise-then-collapse-with-reflection is universal.)

**Takeaway (4):** the specular/diffuse energy-threshold model (Eq. 2.41) is the accepted feature-
scale reflection treatment, and empirical angular yields (Chang-Sawin, Kress) already fold the
reflection loss into f(θ). Grazing reflection is a *sidewall/HAR-transport* effect; it is not the
mechanism of the 6→8 kW **depth** saturation (that is the vertical supply ceiling), but it must
stay engaged so the extra 8 kW energy that reflects off sidewalls is not spuriously converted to
floor removal.

---

## 5. VERDICT

**Retire the `oxide_etch_yield_scale = 0.5586` knob in favor of unscaled published yields + an
engaged supply ceiling + a keV-knee energy-deposition softening. Literature supports this. Strong.**

Rationale, point by point:

1. **A supply ceiling is mandatory physics, not a fit.** Coburn-Winters (synergy = product of two
   channels), Gray-Tepermeister-Sawin (Langmuir coverage balance, ion-limited saturation regime),
   and Chang-Sawin (yield saturates at fixed coverage) make the ceiling the first-order structure of
   ion-enhanced SiO₂ etching. The reference model itself implements it (SSBM Eq. 2.39; two-step
   complex chain, Table 6.2).

2. **The knob suppresses exactly this physics.** A global 0.559× on the yield amplitude pushes ion
   removal capacity below the neutral supply ceiling *at every power*, so the crossover that should
   occur between 6 and 8 kW never happens and the (unbounded, n = 1) law runs free and overshoots
   (875 nm). Removing the knob and using the published p0 = 0.1384 lets capacity cross the
   ~5.3e16 cm⁻² s⁻¹ complex-formation supply between 6 and 8 kW — reproducing the plateau.

3. **The reference author confirms both the mechanism and the failure mode.** Krüger §6.9.2:
   high-power saturation is "not ion starved, but rather limited by neutral gas transport," and his
   own model's residual is that "the effect of ion energy … might be overestimated." petch's
   overshoot is the same residual, amplified by the knob. Fixing it via the ceiling addresses the
   named cause; the knob only masks it and mis-calibrates absolute rate.

4. **Do not fix it in the yield exponent.** n = 1 vs √E is a low-energy (<~1 keV) distinction
   (Steinbrüchel); at 2.5–3.6 keV neither is trustworthy and the physical compression is
   energy-deposited-below-the-reactive-layer, i.e. an **energy-deposition softening / keV knee**
   (ZBL/LSS, companion memo), not a change of prefactor. Keep n = 1 as the low-energy prefactor;
   add the knee.

5. **Keep grazing reflection engaged** (Eq. 2.41 specular/diffuse; Chang-Sawin / Kress f(θ)). It is
   the accepted treatment and prevents reflected high-power energy from being spuriously counted as
   floor etch, but it is a secondary/lateral effect, not the primary saturation driver.

**Caveats before rerun:**
- The ~5.3e16 cm⁻² s⁻¹ "complex-formation ceiling" is a chemistry-weighted subset of the Table 6.1
  fluorocarbon fluxes (Σ unsaturated FC ≈ 3.1e17 cm⁻² s⁻¹ raw; ceiling = Σ φ_FC × sticking/complex
  conversion). Its absolute value must be **derived from the mechanism's complex-formation
  probabilities, not fit to the 8 kW point** — otherwise we replace one knob with another. The
  held-out 8 kW datum may verify but must not *select* the ceiling.
- With the knob gone, the **base-case 825 nm depth re-calibrates**: absolute rate was the knob's
  declared job (protocol R1.3). A knob-free rerun needs the ceiling + knee to hit base depth *and*
  the O2-ratio and power sweeps simultaneously — that is the actual test of whether this is physics
  or a two-parameter re-fit. Preregister before revealing.
- Exact numbers flagged **[VERIFY]** above (Coburn-Winters per-curve rates; Gray-Sawin crossover
  flux ratio; Standaert JVST page; Kress JVST volume/page) should be page-checked before any of
  them is quoted in a paper. The *structural* claims are all sourced.

---

## 6. Source list

- J. W. Coburn, H. F. Winters, *J. Appl. Phys.* **50**, 3189 (1979). DOI: 10.1063/1.326355.
- D. C. Gray, I. Tepermeister, H. H. Sawin, *J. Vac. Sci. Technol. B* **11**, 1243 (1993).
  DOI: 10.1116/1.586925.
- J. P. Chang, H. H. Sawin, *J. Vac. Sci. Technol. A* **15**, 610 (1997). DOI: 10.1116/1.580692.
- C. Steinbrüchel, *Appl. Phys. Lett.* **55**, 1960 (1989). DOI: 10.1063/1.102336.
- J. D. Kress et al., MD reflection off Si/SiO₂ (~1999), petch DOI: 10.1116/1.581948 [venue VERIFY].
- S. Huang, M. J. Kushner et al., *J. Vac. Sci. Technol. A* **37**, 031304 (2019).
  DOI: 10.1116/1.5090904.
- F. Krüger, PhD thesis, Univ. Michigan (2024). DOI: 10.7302/23106.
  §2.1.3 (SSBM Eq. 2.39), §2.2.1 (Eq. 2.40 yield, Eq. 2.41 reflection), §6.4 / Tables 6.1–6.2
  (fluxes + mechanism), §6.9.2 (LF-power saturation statement).
- T. E. F. M. Standaert, G. S. Oehrlein et al., fluorocarbon SiO₂ surface-model work
  (JVST; exact vol/page [VERIFY]).
- Dual-frequency Ar/C₄F₈ HAR-SiO₂ ion-energy/flux study, *Materials* (MDPI), PMC10222222.
- Companion: `RESEARCH_ENERGY_DEPOSITION_ANCHORS_2026-07-23.md` (ZBL/LSS + Sigmund
  deposited-energy-in-layer, for the keV knee).
