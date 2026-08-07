# Tier-1 reactor-scale model — research + design (2026-07-24)

Author: research/design pass 2026-07-24 (Stage E of `PROGRAM_ROADMAP_2026-07-24.md`).
Status: **design doc, not built.** No code committed. `[VERIFY]` marks a locator (DOI, volume,
figure, affiliation) I could not confirm to primary source in this pass — most AIP/AVS bodies
(`pubs.aip.org`) 403 to automated fetch, so several DOIs are reconstructed from abstract pages and
must be checked against the paper before they enter a receipt.

## 0. Purpose and the one contract this must satisfy

petch consumes a **feature-plane boundary**: for each species, a wall flux (m⁻² s⁻¹) plus a
phase-space law — a thermal half-Maxwellian for neutrals, a joint energy–angle distribution
(IEAD/IEADF) for ions, a Maxwellian for the electron closure. The authoritative object is
`PlasmaBoundaryState` (a tuple of `SpeciesBoundaryState`), built today by
`src/petch/reactor_boundary.py`. The Tier-1 model's **only job** is to turn a *recipe* — RF
powers/frequencies (W), gas flows (sccm), pressure (mTorr), gap (mm), wall/wafer temperature (K) —
into exactly that boundary object, with the same provenance and evidence-grade discipline the file
already enforces.

Why this unlocks the roadmap: every industry dataset that publishes a recipe but **not** a
feature-plane boundary (TEL 3D-HARC profiles B1, 3D-NAND charging D1, CD-SAXS D3 in
`RESEARCH_INDUSTRY_DATASETS_2026-07-23.md`) is currently unusable because the boundary is missing.
Build the recipe→boundary surrogate once and B1/D1/D3 all open together. The Kushner/Krüger
Tier-A papers already publish their HPEM boundary, so they need **no** surrogate — which is
precisely why Krüger is the perfect *closed-loop test* (§4): we independently know both the recipe
and the boundary HPEM produced from it, so we can score the surrogate against a known answer before
ever pointing it at a dataset where the boundary is unknown.

### The contract, concretely (what the model must emit)

Reading `reactor_boundary.py`, the boundary is assembled from these existing primitives — the
Tier-1 model must produce their **inputs**, nothing more, nothing less:

| Boundary element | Existing constructor | Tier-1 must supply |
|---|---|---|
| Neutral species flux + thermal law | `_thermal_flux_species` / `_direction_marginalized_thermal_flux_species` | per-neutral wall flux Γ (m⁻²s⁻¹), gas temperature T_g (K) |
| Positive-ion flux + IEAD | `Krueger2024DigitizedIEAD.development_species` (2-D signed-angle grid, optional axisymmetric azimuth) | total ion flux Γ_i, joint p(E,θ) on the wafer plane, effective/aggregate ion mass + mixture closure |
| Electron closure (charging current) | `append_global_current_balance_maxwellian_electrons` | T_e (eV) — flux is fixed by global charge neutrality |
| Sheath-resolved ion/electron boundary | `build_diagnostic_virtual_sheath_boundary` (+ `sheath.py` `CollisionlessWaveformSheath`) | n_e, T_e, ion mass, and a **sheath-voltage waveform** V_sh(t) |
| Deck of raw fluxes (provenance-bound) | `TabulatedReactorFluxDeck` / `build_tabulated_reactor_boundary` | the full species→flux vector + a `source_sha256` + a validation reference |

The cleanest integration point is: **the Tier-1 model emits a `TabulatedReactorFluxDeck` (the flux
vector) plus a sheath-voltage waveform**, and the existing `build_diagnostic_virtual_sheath_boundary`
+ `build_tabulated_reactor_boundary` machinery converts those into the kinetic boundary. That means
the Tier-1 model does **not** invent a new IEAD sampler — it feeds the sheath module the plasma-side
state (n_e, T_e, ion mass) and V_sh(t), and the *already-audited* `CollisionlessWaveformSheath`
produces the finite-transit IEAD. The Tier-1 model owns two things: **(a) the global chemistry that
gives fluxes, densities, T_e; (b) the sheath-voltage waveform from the bias recipe.** Everything
downstream is reuse.

---

## 1. Literature — state of the art

### 1.1 Global (0-D / volume-averaged) model formalism — Lieberman & Lichtenberg

The canonical formalism is Lieberman & Lichtenberg, *Principles of Plasma Discharges and Materials
Processing*, 2nd ed. (Wiley, 2005), ch. 10 ("Particle and Energy Balance in Discharges") — the
reference every global-model paper builds on. The founding worked example for etch gases is
**Lee, Graves, Lieberman, Hess — "Global model of Ar, O₂, Cl₂, and Ar/O₂ high-density plasma
discharges" — JVST A 13(2), 368 (1995)** ([NASA/ADS](https://ui.adsabs.harvard.edu/abs/1995JVSTA..13..368L/abstract)),
which fixes the template still used today:

- **Particle balance** per species s: `d n_s/dt = Σ (generation − loss)` from volume reactions
  (electron-impact ionization/dissociation/attachment, ion–ion/ion–neutral, metastable pooling)
  plus **wall losses** for charged species (ambipolar/Bohm flux to area A_eff) and for radicals
  (surface recombination/sticking, coefficient γ). Flow in (sccm→particles/s) and pumping out close
  the neutral balance.
- **Electron energy (power) balance**: `P_abs = e Σ_s Γ_s A [ε_c(T_e) + ε_i + ε_e]` — absorbed power
  equals the collisional energy cost ε_c(T_e) per electron–ion pair created (summed over inelastic
  channels), plus the ion kinetic energy lost to the sheath ε_i and the electron energy ε_e ≈ 2T_e
  carried to walls. This is the equation that **closes T_e**: given P_abs, geometry, and cross
  sections, T_e is the root of the power balance; densities then follow from the particle balances.
- **Edge-to-center factors** h_l, h_R (Godyak/Lieberman) convert center density to the wall flux
  through the effective area A_eff — the single most important geometric closure, and pressure
  dependent.

This is a **factor-of-2 tool by construction** (0-D averages, assumed EEDF, γ's uncertain) — stated
honestly by the authors and central to our accuracy declaration in §4.

### 1.2 Gudmundsson's global-model series — the liftable chemistry

Gudmundsson's group has published the most complete, *documented* reaction sets with rate
coefficients for exactly our gases:

- **Ar/O₂** — Gudmundsson et al., "Oxygen discharges diluted with argon: dissociation processes,"
  *Plasma Sources Sci. Technol.* **16**, 399 (2007), DOI [10.1088/0963-0252/16/2/025](https://iopscience.iop.org/article/10.1088/0963-0252/16/2/025).
  Full O₂ set with O(³P), O(¹D), O₂(a¹Δ), O⁻, O₂⁺, O⁺, plus the Ar dilution channels.
- **CF₄** — "A global model study of low pressure high density CF₄ discharge," *Plasma Sources Sci.
  Technol.* **28**, 025005 (2019) [VERIFY vol/page], DOI [10.1088/1361-6595/aaf412](https://iopscience.iop.org/article/10.1088/1361-6595/aaf412).
  Revised CF₄ set: neutrals CF₄/CF₃/CF₂/CF/F₂/F/C, metastables, ions CF₃⁺/CF₂⁺/CF⁺/F⁺/F⁻, with
  updated rate coefficients — the modern replacement for older CF₄ sets.
- **O₂/F₂** (SF₆ decomposition proxy) — "Global model of plasma chemistry in a low-pressure O₂/F₂
  discharge," *J. Phys. D: Appl. Phys.* **35**, 328 (2002) [VERIFY page], DOI
  [10.1088/0022-3727/35/4/308](https://iopscience.iop.org/article/10.1088/0022-3727/35/4/308).
- **SF₆/O₂(/Ar)** — the directly liftable set for the Si/SF₆-O₂ arm is the ICP global model
  "Modeling of inductively coupled plasma SF₆/O₂/Ar plasma discharge: effect of O₂ on the plasma
  kinetic properties," *JVST A* **32**(2), 021303 (2014) [VERIFY DOI 10.1116/1.4859376]
  ([abstract](https://pubs.aip.org/avs/jva/article-abstract/32/2/021303/985312/)); complemented by
  "Numerical Study of SF₆/O₂ Plasma Discharge for Etching Applications," *Plasma Chem. Plasma
  Process.* (2021), DOI [10.1007/s11090-021-10170-x](https://link.springer.com/article/10.1007/s11090-021-10170-x)
  (78 species reaction set). Also the older Kokkoris/Turner-lineage "A global model for SF₆ plasmas
  coupling reaction kinetics in the gas phase and on the surface of the reactor walls" (ResearchGate
  mirror) for the wall-loss closure.

These are the **published rate tables ready to lift** for the SF₆/O₂ (silicon, Belen arm) and CF₄
(fluorocarbon proxy) directions.

### 1.3 Kushner HPEM — the reference standard (target, not dependency)

The Hybrid Plasma Equipment Model (Kushner group, U. Michigan) is the field's reference
2-D/3-D reactor simulator: a module hierarchy coupling fluid drift-diffusion for bulk ions/electrons
with a Monte-Carlo treatment of the hot sheath-accelerated electrons and a Boltzmann/EEDF solver,
producing self-consistent densities, fluxes, and **printed wafer IEADs** — see Kushner, "Hybrid
modelling of low temperature plasmas for fundamental investigations and equipment design," *J. Phys.
D: Appl. Phys.* **42**, 194013 (2009), DOI [10.1088/0022-3727/42/19/194013](https://www.researchgate.net/publication/230979059)
and the group tutorial ([cpseg intro PDF](https://cpseg.eecs.umich.edu/pub/odp_files/hpem_intro.pdf)).
HPEM is the **standard we validate against**, not a dependency: Krüger's thesis fluxes/IEAD (already
in `data/experimental/krueger_2024/`) *are* HPEM output, so reproducing them (§4 Gate 1) is a direct
Tier-1-vs-HPEM benchmark on an identical recipe. We are **not** reimplementing HPEM — we are building
the 0-D surrogate that, on a factor-of-2 band, lands where HPEM lands.

### 1.4 Volume-averaged CCP with multi-frequency bias — Turner, Chabert

For CCP (Krüger's reactor is a dual-frequency CCP, 1/40 MHz), the extra physics beyond an ICP
global model is **capacitive sheath power coupling** and the **self-bias / voltage-division** that
sets ion energy. Canonical references:

- Chabert & Braithwaite, *Physics of Radio-Frequency Plasmas* (Cambridge, 2011) — chapters on the
  global model, RF sheaths, and single-/multi-frequency CCP; the homogeneous-discharge model that
  gives sheath voltage, stochastic + ohmic heating split, and the current-driven vs voltage-driven
  regimes. This is the primary CCP reference.
- Turner & Chabert, "A radio-frequency sheath model for complex waveforms,"
  [arXiv:1212.2612](https://arxiv.org/pdf/1212.2612) — analytic sheath charge–voltage relation for
  arbitrary (multi-harmonic, tailored) waveforms; the basis for turning a *bias recipe* (amplitudes,
  frequencies, phases) into V_sh(t).
- The **dual-frequency functional separation** (a high frequency controls plasma density/flux, a low
  frequency controls ion energy) is the design principle behind Krüger's 40 MHz-source/1 MHz-bias
  split; it justifies our "flux from source power, energy from bias waveform" architecture (§3).

### 1.5 Sheath IEDF/IEAD from recipe without a full HPEM — the state of the art

This is the crux of the whole model and the part petch already partly owns. The lineage:

- **Benoit-Cattin & Bernard, J. Appl. Phys. 39, 5723 (1968)** — the collisionless bimodal IED:
  constant sheath width, uniform field, sinusoidal sheath voltage, ions entering at Bohm speed give
  the classic **two-horn (bimodal) distribution** whose splitting ΔE ∝ (1/ωτ_ion), narrowing with
  frequency and ion mass. Its high-frequency limit is the **arcsine distribution** (density ∝
  1/√(1−((E−Ē)/ΔE)²), horns at ±ΔE/2) — **this is exactly the form petch already implements** in
  `chemistry.py::_ied_yield` (`'bimodal'` branch, cited to Kawamura 1999).
- **Kawamura, Vahedi, Lieberman, Birdsall — "Ion energy distributions in rf sheaths: review,
  analysis and simulation" — Plasma Sources Sci. Technol. 8, 313 (1999)**, DOI
  [10.1088/0963-0252/8/3/202](https://iopscience.iop.org/article/10.1088/0963-0252/8/3/202) — the
  definitive review spanning collisionless→collisional and low→high ωτ_ion; the analytic bimodal
  form and its corrections. **The reference petch's IED form is already keyed to.**
- **Collisional corrections**: the "sheath model and IED for all radio frequencies" analytic model
  (Economou group; [UH PDF](https://www.chee.uh.edu/sites/chbe/files/faculty/economou/jap_sheath.pdf))
  and the charge-exchange-collision low-energy tail (secondary peaks from ions born inside the
  sheath). At 10–40 mTorr HARC pressures the sheath is **mildly collisional** — the collisionless
  arcsine is a good leading model but the low-energy CX tail should be a declared correction.
- **Angular distribution (IADF)**: Manenschijn/Goedheer lineage and specifically
  **"Analytical model for ion angular distribution functions at rf biased surfaces with collisionless
  plasma sheaths," J. Appl. Phys. 92, 7032 (2002)**, DOI [10.1063/1.1519941](https://pubs.aip.org/aip/jap/article-abstract/92/12/7032/)
  [VERIFY] — the angle spread comes from **conserved transverse (thermal) velocity** at the sheath
  edge divided by the normal velocity gained across the sheath: θ ≈ atan(v⊥/v‖), with v⊥ set by ion
  temperature T_i and v‖ by the sheath voltage. Narrow at high energy, wider at low — the physics
  petch's sheath module already encodes via `ion_tangential_temperature_eV`.
- **Measured-IEDF anchor / ground truth**: **Sobolewski (NIST)** — "Ion energy distributions at
  rf-biased wafer surfaces" and the noninvasive rf current/voltage IED method (validated in Ar and
  CF₄ ICP, 10 mTorr, 0.1–20 MHz bias) — [NIST Plasma Process Metrology](https://www.nist.gov/programs-projects/plasma-process-metrology);
  and **Edelberg & Aydil** measured IEDFs at rf-biased surfaces. These are the experimental datasets
  that tell us how good the analytic sheath IED actually is (Gate-2-adjacent).

**State of the art answer**: you can predict a usable IEAD from recipe *without* HPEM by (i) getting
n_e, T_e, ion flux from a global model, (ii) getting the sheath-voltage waveform V_sh(t) from the
bias recipe via an analytic RF-sheath charge–voltage model (Turner/Chabert), and (iii) pushing ions
through a finite-transit collisionless sheath (Benoit-Cattin/Kawamura form for energy; conserved-T_i
transverse velocity for angle). Accuracy is **factor-of-2 on absolute fluxes** and **~10–20% on
mean ion energy / IED horn positions** when the sheath voltage is right — good enough for
calibrated-transfer feature prediction, not blind absolute rates. This is **already the exact
architecture of petch's `build_diagnostic_virtual_sheath_boundary`** — the Tier-1 model just feeds it.

---

## 2. Chemistry sets with published rate coefficients ready to lift

Ranked by immediacy to petch's live mechanisms (fluorocarbon-SiO₂ first, Si-SF₆/O₂ second).

| Mixture | Primary liftable source (reaction set + rate table) | DOI | Notes / match |
|---|---|---|---|
| **Ar/CF₄/O₂ (+CHF₃)** — Krüger's mix | Base CF₄ set: Gudmundsson, "Global model study of low-pressure high-density CF₄," *Plasma Sources Sci. Technol.* 28 (2019). Ar/CF₄/O₂/CHF₃ ICP densities + Si etch: **Efremov et al.**, "A Comparison of CF₄, CHF₃ and C₄F₈ + Ar/O₂ ICPs for Dry Etching," *Plasma Chem. Plasma Process.* (2021) | [10.1088/1361-6595/aaf412](https://iopscience.iop.org/article/10.1088/1361-6595/aaf412) ; [10.1007/s11090-021-10198-z](https://link.springer.com/article/10.1007/s11090-021-10198-z) | **Exact match** to Krüger's Ar/CF₄/O₂ waveform-tailoring papers (A4). Efremov papers hand you Langmuir-probe-validated densities to gate against. |
| **Ar/C₄F₈/O₂** | **Kokkoris et al.**, "A global model for C₄F₈ plasmas coupling gas phase and wall surface reaction kinetics," *J. Phys. D: Appl. Phys.* 41, 195211 (2008) — full set + wall-loss/sticking closure, validated on F/CF₂ densities. Reaction mechanism II: "Properties of c-C₄F₈ ICPs II: plasma chemistry and reaction mechanism for Ar/c-C₄F₈/O₂" (Vasenkov/Kushner lineage) | [10.1088/0022-3727/41/19/195211](https://iopscience.iop.org/article/10.1088/0022-3727/41/19/195211) | **Match** to Huang/Huard 2019 (A2, C₄F₈ HARC AR≤80). **Not a C₄F₆ closure:** Benck 2003 and Kim 2021 measure C₄F₆-specific ion composition, IEDs, electronegativity, and parent/fragment signals. Porting the C₄F₈ network is a prior that must be revalidated reaction-by-reaction, not a “small delta.” |
| **Ar/C₄F₆/O₂** — Krüger base case (A1) | No standalone published C₄F₆ *global* set found this pass — build C₄F₆ as a **C₄F₈-set delta** (heavier unsaturation → more polymerizing CF/CF₂, less F). The self-consistent boundary already exists as HPEM output in `data/experimental/krueger_2024/` (Table I), so for Gate 1 we **do not need** the C₄F₆ gas-phase set — we validate the *sheath+flux-partition* layer against Krüger's HPEM fluxes directly. | — [VERIFY: search dedicated C₄F₆ set] | The gas-phase C₄F₆ set is the **hardest gap**; deferred until after Gate 1 because Gate 1 tests the recipe→boundary layer with HPEM fluxes as the reference. |
| **SF₆/O₂(/Ar)** — Si arm (Belen/de Boer) | ICP SF₆/O₂/Ar global model, *JVST A* 32, 021303 (2014) [VERIFY 10.1116/1.4859376]; + *Plasma Chem. Plasma Process.* 2021 (78-species set) 10.1007/s11090-021-10170-x; + Gudmundsson O₂/F₂ J. Phys. D 35 (2002) for the F₂ channels | [10.1007/s11090-021-10170-x](https://link.springer.com/article/10.1007/s11090-021-10170-x) | **Match** to petch's Belen Si mechanism and de Boer/Blauw SF₆-O₂ ARDE (E3). Highly electronegative — exercises the negative-ion balance. |

**Immediately liftable, zero new physics gap**: CF₄ (Gudmundsson 2019) and C₄F₈ (Kokkoris 2008) — both
have complete printed rate tables and independent density validation. **SF₆/O₂** liftable but
electronegative-heavy. **C₄F₆** is the one genuine gas-phase gap, and Gate 1's design specifically
routes around it (validate against HPEM fluxes, not a from-scratch C₄F₆ chemistry).

---

## 3. Architecture

```
recipe (W_src, W_bias, f_src, f_bias, φ, sccm_i, p, gap, T_wall)
        │
        ▼
┌─────────────────────────────────────────────┐
│ (A) GLOBAL CHEMISTRY CORE  (0-D ODEs)         │
│   • species particle balances dn_s/dt         │
│   • electron power balance  → T_e             │
│   • electronegative closure (n_-, α)          │
│   • edge factors h_l,h_R → wall fluxes Γ_s    │
│   outputs: {Γ_s}, n_e, T_e, ion flux Γ_i,     │
│            ion-species mix, gas temp T_g      │
└───────────────┬─────────────────────────────┘
                │ n_e, T_e, Γ_i, ion mass
                ▼
┌─────────────────────────────────────────────┐        ┌──────────────────────────┐
│ (B) POWER/SHEATH-VOLTAGE MODULE               │        │ (C) EEDF assumption knob │
│   • ICP: P_abs → source; sheath = floating    │        │  Maxwellian (default)    │
│   • CCP: voltage division → self-bias Vdc      │◄───────┤  two-temperature         │
│   • multi-freq: Σ harmonics → V_sh(t) waveform │        │  measured EEDF (import)  │
│   outputs: PeriodicSheathVoltage V_sh(t)      │        └──────────────────────────┘
└───────────────┬─────────────────────────────┘
                │ V_sh(t), n_e, T_e, ion mass, Γ_i
                ▼
┌─────────────────────────────────────────────┐
│ (D) REUSE petch sheath+boundary machinery     │
│   PlasmaDiagnosticState + PeriodicSheathVoltage│
│   → build_diagnostic_virtual_sheath_boundary  │
│   → CollisionlessWaveformSheath → IEAD         │
│   + TabulatedReactorFluxDeck (neutrals)        │
│   + append_global_current_balance_electrons    │
└───────────────┬─────────────────────────────┘
                ▼
         PlasmaBoundaryState   ── petch feature engine
```

### (A) Global chemistry core — the ODE system

State vector: number densities {n_s} for every neutral, ion, metastable, and electron; plus T_e as
an algebraic unknown closed by power balance (or a coupled ODE for electron energy density).

- **Particle balance** (per species): `dn_s/dt = Σ_r ν_r,s k_r(T_e) Π n_reactants − L_wall,s`.
  Rate coefficients k_r(T_e) come from the lifted sets (§2) — either Arrhenius/`A·T_e^b·exp(−E/T_e)`
  fits or tabulated k(T_e) from cross-section integration. Electron-impact k's assume the chosen EEDF
  (module C).
- **Wall loss**: charged species leave at the **Bohm flux** `Γ = h·n·u_B` through effective area
  A_eff = h_l·(2πR²) + h_R·(2πRL) with Godyak edge factors h_l(p,L), h_R(p,R); radicals leave with
  **surface loss probability γ_s** (recombination/deposition on chamber walls — the dominant declared
  uncertainty, §5). For fluorocarbons, the wall is a *reactive polymer film* (Kokkoris' surface
  submodel) — γ is coverage-dependent; v1 uses a constant declared γ per radical.
- **Electron power balance**: `P_abs = e·A_eff·Σ_ions Γ_ion·[ε_c(T_e) + ε_i] + Γ_e·ε_e`, solved for
  T_e. ε_c(T_e) = Σ (inelastic energy cost per channel)·(k_channel/k_iz) is the collisional energy
  per ion-pair; ε_i is the ion energy hitting the sheath (≈ V_sheath); ε_e ≈ 2T_e.
- **Electronegative discharge handling** (O₂, SF₆/O₂, CF₄/O₂): add negative-ion species (O⁻, F⁻) with
  attachment (source) and detachment/recombination (loss); electronegativity α = n₋/n_e reshapes the
  ambipolar field and the edge factors (Lieberman & Lichtenberg ch. 10.4; Gudmundsson O₂ set). At
  high α the discharge can go to an ion–ion core with a **flat-topped** density profile → different
  h factors. v1: implement the two-region (electronegative core + electropositive edge) h-factor
  correction from L&L; declare it as the electronegative closure.
- **Gas temperature T_g**: v1 declares T_g from the recipe (wall/wafer T + a fixed rise), used only
  to set the neutral thermal law and gas density N = p/(k T_g). A T_g energy balance is a later
  refinement.

Numerics: stiff ODE integrate to steady state (`scipy.integrate.solve_ivp` BDF), or Newton on the
steady-state residual. ~20–60 species → trivial CPU cost (ms), so this layer is **not** the compute
frontier and can be made differentiable later (Stage F) via the same dual-number path.

### (B) Power deposition split & sheath-voltage waveform

- **ICP source**: absorbed power P_abs from the coil sets the plasma density (transformer/collisional
  heating); the wafer sheath is a **floating (thin, low-voltage) sheath** unless separately biased.
  Ion energy ≈ plasma potential + any bias.
- **CCP**: power couples through the sheaths. The **self-bias V_dc** on the powered electrode comes
  from the capacitive voltage division between the two sheaths (area-ratio asymmetry) — Chabert &
  Braithwaite ch. on single-frequency CCP. The **sheath voltage waveform V_sh(t)** is the object we
  need for the IEAD.
- **Multi-frequency / tailored waveform**: V_sh(t) = Σ_h V_h·cos(2π f_h t + φ_h). Dual-frequency
  functional separation: high f sets flux (feeds back into A's power balance), low f sets energy (the
  horns of the IED). Turner–Chabert complex-waveform sheath model maps (amplitudes, frequencies,
  phases) → the charge–voltage relation → V_sh(t). petch's `PeriodicSheathVoltage` **already accepts a
  harmonic list** (`waveform.harmonic_number`) — the Tier-1 module just fills its coefficients.
- v1 pragmatic closure: for a CCP bias recipe, take the sheath-voltage *amplitude* from the DC
  self-bias + RF amplitude (analytic voltage division), assemble the harmonic waveform, and feed it
  straight to `CollisionlessWaveformSheath`. Where a paper prints the sheath voltage or self-bias
  (Krüger does), use that as the ground-truth check on module B.

### (C) EEDF assumption — an explicit, swappable knob

Three declared options, recorded in provenance:
1. **Maxwellian** (default) — rate coefficients pre-integrated over a Maxwellian at T_e; simplest,
   standard for global models, adequate at higher pressure/collisionality.
2. **Two-temperature / bi-Maxwellian** — captures the depleted high-energy tail (Druyvesteyn-like) at
   low pressure; changes ionization/dissociation balance. Optional.
3. **Measured / imported EEDF** — if a paper supplies an EEDF or T_e, override the assumption and
   integrate k(T_e) against it. This is the honest path when validating against a specific reactor.

The EEDF choice **is** a declared uncertainty (§5); v1 ships Maxwellian and records it.

### (D) Mapping onto petch's boundary — exact

This is the load-bearing section. The Tier-1 model produces:

1. **Neutral fluxes** → build a `TabulatedReactorFluxDeck` with one `ReactorSpeciesFlux(role="neutral",
   charge_number=0, mass_amu=…, flux_m2_s=Γ_s, evidence_kind="published_reactor_model_output")` per
   radical, `source=` recipe id, `source_sha256=` hash of the recipe+model version. Then
   `build_tabulated_reactor_boundary(..., neutral_temperature_K=T_g)` emits the half-Maxwellian neutral
   species via the existing `_thermal_flux_species` / `_direction_marginalized_thermal_flux_species`.
   **No new neutral code.**
2. **Ion flux + IEAD** → construct a `PlasmaDiagnosticState(n_e, T_e, ion_name, ion_mass_amu,
   ion_flux_m2_s=Γ_i, source=…)` and a `PeriodicSheathVoltage` from module B, then call
   `build_diagnostic_virtual_sheath_boundary(diagnostic, waveform, reference_plane_m=…,
   collisionless_justification=…, claim_mode="development")`. The existing `CollisionlessWaveformSheath`
   produces the finite-transit ion IEAD **and** the Maxwellian electron closure. **No new sheath code.**
3. **Aggregate vs resolved ions**: like Krüger's Table I, v1 emits an **aggregate positive-ion mixture**
   with an effective mass + declared `mixture_closure` (reusing the `positive_ion_mixture` role and the
   `development_species` effective-mass path). Species-resolved ion fluxes (CF₃⁺ vs CF₂⁺ …) are a later
   upgrade that unlocks `claim_mode="predictive"`.
4. **Electron charging current**: `append_global_current_balance_maxwellian_electrons(boundary,
   electron_temperature_eV=T_e, …)` — flux fixed by global neutrality, exactly as the file already does.

The Tier-1 model therefore adds **two new small modules** (global chemistry core A, sheath-voltage
builder B) and a thin **adapter** that fills the existing constructors. It touches `reactor_boundary.py`
only to add a `build_global_model_boundary(recipe, ...)` entry point that internally runs A+B and calls
the existing builders. The reference shape to mirror is `build_krueger_2024_development_boundary` — same
provenance keys, same `claim_mode`/`supports_prediction` discipline, same evidence-kind gating.

---

## 4. Validation ladder

Global models are **factor-of-2 tools**. Stated plainly: absolute wall fluxes and densities are
expected within ~2× of HPEM/experiment; **ratios** (flux compositions, O/F balance, electronegativity
trends) and **mean ion energy / IED horn positions** are better (~10–30%) when the sheath voltage is
pinned. **Downstream consequence for feature claims**: at first this supports **calibrated-transfer**
prediction (calibrate one condition's absolute scale, predict the *shape* of a sweep), **not blind
absolute rate prediction**. That is the honest ceiling and it must be declared in every receipt.

### Gate 1 — reproduce Krüger's HPEM boundary from his recipe (closed loop)
- **Recipe** (from Table I provenance in `load_krueger_2024_reactor_flux_deck`): 10 mTorr C₄F₆/Ar/O₂,
  140/100/105 sccm, 1/40 MHz CCP; base case Table IV geometry.
- **Reference** (already in repo): Table I wall fluxes (`base_case_boundary_fluxes.csv`), Figure-4
  combined IEAD (`digitized_figure4_iead.csv`), Figure-16 flux/IEAD sweeps.
- **What we score, and why it dodges the C₄F₆ gap**: Gate 1 has **two sub-gates**.
  - **1a (sheath layer, do first)**: take Krüger's *published* aggregate ion flux + n_e/T_e and drive
    module B+D to reproduce the **Figure-4 IEAD** (mean energy, horn positions/split, angular width)
    within a declared band (target: mean energy within 15%, horn split within 25%). This validates the
    recipe→sheath→IEAD path using HPEM's own plasma state — **no C₄F₆ chemistry needed**.
  - **1b (chemistry layer)**: run module A on the C₄F₆/Ar/O₂ recipe and compare the **neutral flux
    vector** (C₃F₄, C₂F₃, CF, CF₂, CF₃, O) and total ion flux to Table I within a factor of 2, and the
    **flux ratios** (O/CFx, CF₂/CF₃) within ~30%. This is where the C₄F₆-as-C₄F₈-delta set is exercised;
    a miss here localizes the gas-phase gap rather than contaminating the sheath validation.
- **Declared band**: preregister the bands above **before** running (receipt machinery, §5). Gate 1
  passes if 1a passes and 1b lands within the factor-of-2 flux / 30% ratio band.
- This is the **perfect closed-loop test**: recipe *and* boundary both known, from the same paper petch
  already validated the feature engine against.

### Gate 2 — independent global-model measured densities
- Reproduce an **independent** paper's *measured* densities/T_e, not HPEM output. Best targets:
  Gudmundsson CF₄ (2019) or O₂/Ar (2007) density-vs-pressure curves (the paper validates its own set),
  and **Efremov** Ar/CF₄/O₂ / CHF₃ Langmuir-probe densities (2021) — these give F, CF₂, CF₃, ion-density
  vs mixing ratio. SF₆/O₂ arm: the ICP SF₆/O₂/Ar 2014 model's density-vs-%O₂. Score density trends
  (E2/E3) within factor-of-2 absolute, correct monotonic trends.
- Also fold in a **measured-IEDF** check against Sobolewski/Edelberg-Aydil CF₄ IEDs (10 mTorr) to grade
  the sheath module against experiment, not just HPEM.

### Gate 3 — one TEL-recipe blind attempt
- Take a **TEL** published recipe (B1 Nishizuka 2024 3D-HARC, or the B2 cryo condition table if it
  surfaces) with **no published boundary**, generate the Tier-1 boundary, run the feature engine, and
  compare **profile shape** (twist/distortion/CD-vs-depth) to the SEM/3D data. Declared **blind** and
  **calibrated-transfer** (calibrate absolute scale on one condition; predict the sweep). This is the
  first true recipe→profile prediction on a dataset whose boundary we never had — the payoff of the
  whole Tier-1 build. Expected outcome: correct *trends and shapes*, not absolute depths, at first.

**Accuracy declaration (standing)**: Tier-1 outputs are labelled `claim_mode="development"` and
`supports_prediction=False` until (a) species-resolved ion fluxes replace the aggregate closure and
(b) Gate 2 passes on measured (not HPEM) data. Even then, feature-scale claims built on a Tier-1
boundary are **calibrated-transfer** grade (E3), never E4 blind-absolute, unless a measured wafer IEAD
is supplied.

---

## 5. Knob discipline

Per the standing doctrine (no fitted knobs — derive, measure, or declare as fab-measurable), classify
every input:

**Physics (derived / from cited sources — not knobs):**
- Electron-impact cross sections & rate coefficients k_r(T_e): from the lifted sets — Gudmundsson CF₄
  (2019), Kokkoris C₄F₈ (2008), SF₆/O₂ (2014/2021), O₂/Ar (2007). **Cited, versioned, hashed.**
- Ion mass, Bohm speed, sheath charge–voltage relation, finite-transit ion dynamics: first-principles
  (already in `sheath.py`, receipted).
- Edge-to-center factors h_l, h_R: Godyak/Lieberman analytic forms (derived, pressure-dependent).
- Global charge neutrality electron flux: exact constraint (already enforced).

**Declared uncertain (measurable, carried in provenance, swept — not tuned to fit):**
- **Wall/surface loss coefficients γ_s** (radical recombination/deposition). The single largest
  uncertainty. Declare a literature range per radical; **sweep it**, don't fit it; report sensitivity.
  For fluorocarbons, γ is film-coverage-dependent (Kokkoris) — v1 uses constant declared γ and flags it.
- **EEDF shape** (module C): declared Maxwellian in v1; two-temperature/measured are labelled
  alternatives. The choice is recorded, its effect on T_e/densities is reported.
- **Power coupling efficiency** (fraction of nominal RF power actually absorbed): declared, ideally
  pinned by a measured density or self-bias where a paper gives one (Krüger prints self-bias).
- **Gas temperature T_g**: declared from recipe in v1; a measurable quantity, not a fit parameter.
- **Ion mixture effective mass + closure**: explicit declared closure (exactly as
  `build_krueger_2024_development_boundary` already does), until species-resolved.

**Preregistration / receipt machinery at reactor scale** (reuse what `reactor_boundary.py` already
enforces):
- Every input source (rate set, recipe) carries a `source_sha256` and a citable `source_location`;
  the `TabulatedReactorFluxDeck` / boundary provenance already demand this.
- `claim_mode ∈ {development, predictive}` and `supports_prediction` gate exactly as today — a Tier-1
  boundary cannot claim prediction while any input is `assumed`/aggregate.
- **Preregister the Gate-1 bands before running** (write the numeric bands + declared γ range to a
  committed doc, then execute) — the same freeze→reveal discipline used for the Krüger blind campaign.
- Sensitivity to every declared-uncertain input is a **required output** of a Tier-1 run (a small
  ensemble over the γ range and EEDF choice), so the factor-of-2 band is *reported*, never hidden.

---

## 6. First build slice (v1, ~1–2 sessions → Gate 1)

Smallest thing that reaches Gate 1a and starts 1b. Order matters — do the sheath layer first because
it needs no new chemistry and uses Krüger's own plasma state.

**Session 1 — sheath/boundary adapter (Gate 1a), zero new chemistry:**
1. `src/petch/reactor_tier1/sheath_voltage.py` — module B v1: map a CCP bias recipe (V_lf, V_hf, f_lf,
   f_hf, phases, area asymmetry) → DC self-bias + harmonic `PeriodicSheathVoltage`. For Gate 1a, seed
   it with Krüger's published self-bias/sheath-voltage amplitude (ground-truth check).
2. `build_global_model_boundary(...)` adapter in a new `reactor_tier1/boundary_adapter.py` (mirrors
   `build_krueger_2024_development_boundary`): given {Γ_s}, n_e, T_e, Γ_i, ion mass, V_sh(t) → assemble
   `TabulatedReactorFluxDeck` (neutrals) + `PlasmaDiagnosticState`+waveform →
   `build_diagnostic_virtual_sheath_boundary` → `append_global_current_balance_maxwellian_electrons`.
   All reuse; the adapter only fills existing constructors.
3. **Gate 1a test** (`tests/test_reactor_tier1_sheath.py`): feed Krüger's Table-I ion flux + published
   n_e/T_e + sheath voltage; assert the produced IEAD's mean energy within 15% and horn split within
   25% of the Figure-4 digitization already in the repo. Preregister these bands in a committed stub.

**Session 2 — minimal global chemistry core (start Gate 1b):**
4. `src/petch/reactor_tier1/global_model.py` — module A v1 for **CF₄/O₂/Ar** (Gudmundsson 2019 set,
   the cleanest complete table), electropositive→mildly-electronegative: particle balances + electron
   power balance (Maxwellian EEDF) + Godyak h-factors + Bohm wall flux + constant declared radical γ.
   Steady-state solve. Output {Γ_s}, n_e, T_e, Γ_i.
5. **Gate 2-lite test**: reproduce a Gudmundsson/Efremov CF₄(/O₂/Ar) density-vs-pressure or
   density-vs-mixing point within factor-of-2 + correct trend (validates the ODE core independently
   before trusting it on C₄F₆).
6. **Gate 1b (stretch)**: add the C₄F₆-as-C₄F₈-delta neutral set; run Krüger's recipe; compare neutral
   flux vector + ratios to Table I within the preregistered factor-of-2 / 30% band. If it misses, the
   residual localizes the C₄F₆ gas-phase gap — documented, not fitted around.

**Deliverable of the slice**: a committed `RESEARCH_REACTOR_TIER1_GATE1_*.md` scorecard (preregistered
bands + results), the two new modules + adapter, and passing Gate-1a / Gate-2-lite tests. Gate 1b
either passes or produces a documented, localized chemistry gap — both are acceptable slice outcomes
under the doctrine.

**Explicitly deferred** (not in v1): species-resolved ion fluxes (predictive-mode unlock),
two-temperature/measured EEDF, collisional CX low-energy IED tail, fluorocarbon coverage-dependent γ,
T_g energy balance, C₄F₆ from-scratch set, and any GPU/differentiable port (Stage F).

---

## Appendix — citation ledger (verify before receipting)

- Lee, Graves, Lieberman, Hess — *Global model of Ar, O₂, Cl₂, Ar/O₂ HD discharges* — JVST A 13(2), 368 (1995).
- Lieberman & Lichtenberg — *Principles of Plasma Discharges and Materials Processing*, 2nd ed., Wiley (2005), ch. 10.
- Gudmundsson et al. — *Oxygen discharges diluted with argon: dissociation processes* — Plasma Sources Sci. Technol. 16, 399 (2007). DOI 10.1088/0963-0252/16/2/025.
- Gudmundsson et al. — *Global model study of low-pressure high-density CF₄ discharge* — Plasma Sources Sci. Technol. 28 (2019). DOI 10.1088/1361-6595/aaf412. [VERIFY vol/page]
- Gudmundsson — *Global model of plasma chemistry in a low-pressure O₂/F₂ discharge* — J. Phys. D 35, 328 (2002). DOI 10.1088/0022-3727/35/4/308. [VERIFY page]
- (SF₆/O₂/Ar ICP) — *Modeling of ICP SF₆/O₂/Ar discharge* — JVST A 32(2), 021303 (2014). [VERIFY DOI 10.1116/1.4859376]
- *Numerical Study of SF₆/O₂ Plasma Discharge for Etching Applications* — Plasma Chem. Plasma Process. (2021). DOI 10.1007/s11090-021-10170-x.
- Kokkoris et al. — *A global model for C₄F₈ plasmas coupling gas phase and wall surface reaction kinetics* — J. Phys. D 41, 195211 (2008). DOI 10.1088/0022-3727/41/19/195211.
- Efremov et al. — *A Comparison of CF₄, CHF₃ and C₄F₈ + Ar/O₂ ICPs for Dry Etching* — Plasma Chem. Plasma Process. (2021). DOI 10.1007/s11090-021-10198-z.
- Kushner — *Hybrid modelling of low temperature plasmas…* — J. Phys. D 42, 194013 (2009). DOI 10.1088/0022-3727/42/19/194013.
- Krüger, Zhang, Luan, Park, Metz, Kushner — *Autonomous hybrid optimization of a SiO₂ plasma etching mechanism* — JVST A 42(4), 043008 (2024). DOI 10.1116/6.0003554.
- Krüger et al. — *Voltage waveform tailoring for HAR etching of SiO₂ using Ar/CF₄/O₂: low fundamental-frequency biases* — Phys. Plasmas 31(3), 033508 (2024). DOI 10.1063/5.0189675. [VERIFY]
- Chabert & Braithwaite — *Physics of Radio-Frequency Plasmas* — Cambridge (2011).
- Turner & Chabert — *A radio-frequency sheath model for complex waveforms* — arXiv:1212.2612.
- Benoit-Cattin & Bernard — *Anomalies in the energy of positive ions extracted from HF discharges* — J. Appl. Phys. 39, 5723 (1968).
- Kawamura, Vahedi, Lieberman, Birdsall — *Ion energy distributions in rf sheaths: review, analysis and simulation* — Plasma Sources Sci. Technol. 8, 313 (1999). DOI 10.1088/0963-0252/8/3/202.
- (Economou group) — *Plasma sheath model and ion energy distribution for all radio frequencies* — J. Appl. Phys. [VERIFY cite]. UH mirror.
- *Analytical model for ion angular distribution functions at rf biased surfaces with collisionless plasma sheaths* — J. Appl. Phys. 92, 7032 (2002). DOI 10.1063/1.1519941. [VERIFY]
- Sobolewski (NIST) — *Ion energy distributions at rf-biased wafer surfaces*; noninvasive rf I/V IED method. NIST Plasma Process Metrology.
