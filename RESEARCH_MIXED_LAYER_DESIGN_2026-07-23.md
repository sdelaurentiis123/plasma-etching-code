# Resolved thin mixed-layer model for fluorocarbon etching in petch — literature review + v1 design

Research + design synthesis, 2026-07-23. Author-of-record: autonomous research agent.
Status: **design document, not code.** Do not commit unless the maintainer decides to.

Scope: replace the lumped surface-site-balance parameters in petch's fluorocarbon chemistry
(`src/petch/surface_kinetics.py` Krüger reduced mechanism; `src/petch/fluorocarbon_lamagna.py`
La Magna–Garozzo three-coverage parity model) with a **resolved thin mixed-layer** per surface
patch: element-resolved composition (C, F, O, Si) of the top ~1–3 nm, ion-induced mixing driven
by the new ZBL stopping module `src/petch/ion_energy_deposition.py`, volatilization, and a
fluorocarbon film reservoir above — all under petch's deterministic transport engine with exact
conservation ledgers.

Citation convention: `[VERIFY]` marks a DOI/volume/page reconstructed from memory and not
re-confirmed against the publisher record during this session. Titles, authors, venues, and the
physics claims are checked against the search record; only the exact locator carries the flag.

---

## 0. What petch has today (the thing being replaced)

- **Belen SF6/O2 on Si** (`chemistry.py`): competitive-Langmuir θ_F/θ_O fixed point, yield
  `Y = A·max(√E − √E_th, 0)`, IED-integrated. This is the *base calibration* chemistry
  (de Boer cryo, `data/experimental/deboer_2002/`). One adsorbing-species-pair Langmuir.
- **La Magna–Garozzo three-coverage** (`fluorocarbon_lamagna.py`): θ_pe (polymer-on-etchant),
  θ_p (polymer), θ_e (etchant) algebraic fixed point + a saturated-polymer finite film
  inventory; ViennaPS 4.6.1 parity. DOI 10.1149/1.1602084.
- **Krüger reduced MCFPM mechanism** (`surface_kinetics.py`): oxide–fluorocarbon *complex*
  coverage + finite polymer inventory; energetic yield `Y = Y_ref·max((E−E_th)/(E_ref−E_th),0)^n·f(θ)`
  (Eq. 2.40 of Krüger 2024, DOI 10.1116/6.0003554); calibrated probabilities (complex-formation
  0.2729, bare-SiO2 sputter 0.0909, complex sputter 0.1384, O-polymer etch 0.0628, polymer dep
  0.0842). These are **experiment-calibrated lumped coefficients**, not first-principles.

All three are **coverage** models: fractional site occupancies, not element inventories. None of
them carries a fluorine *budget* that must close, and none derives the ion channel from a stopping
integral — the energy law is a fitted `(√E−√E_th)` or `(E−E_th)^n` with a fitted threshold and a
fitted magnitude. That is exactly the lumping the mixed-layer model removes.

---

## 1. Phenomenological surface / mixed-layer models (state variables, equations, what they predicted)

### 1.1 Gogolides–Kokkoris detailed surface model (the reference phenomenological FC model)
E. Gogolides, P. Vauvert, G. Kokkoris, G. Turban, A. G. Boudouvis, "Etching of SiO2 and Si in
fluorocarbon plasmas: A detailed surface model accounting for etching and deposition,"
*J. Appl. Phys.* **88**, 5570 (2000), DOI 10.1063/1.1311808. Earlier profile-coupled version:
"SiO2 and Si etching in fluorocarbon plasmas: A detailed surface model coupled with a complete
plasma and profile simulator," *Microelectron. Eng.* (1999/2001) — same group/lineage.

- **State variables (surface coverages):** polymer coverage θ_poly, CFx coverage θ_CFx, atomic-F
  coverage θ_F on the reactive surface; a polymer film that grows/erodes. Fluxes in: F atoms,
  CFx radicals (deposition precursor), ions (flux + energy).
- **Balance form:** competitive site balance — deposition from CFx sticking, removal of polymer by
  F etching + ion sputtering, F consumption into volatile SiFx, ion-enhanced etch of exposed oxide.
  Yields carry an ion-energy `(√E−√E_th)`-type law and an angular factor.
- **What it predicted:** the **etch→deposition transition** as a function of F and CFx radical
  concentration vs ion flux and energy (the "clog/no-clog" boundary in modern language), and,
  coupled to a profile simulator, **RIE lag** (aspect-ratio-dependent rate) for features of
  different AR. This is the canonical phenomenological demonstration that a coverage/film model
  reproduces both the deposition cliff and ARDE.

### 1.2 Kushner-lineage surface-site-balance module (SKM in HPEM/MCFPM)
D. Zhang & M. J. Kushner, "Investigations of surface reactions during C2F6 plasma etching of SiO2
with equipment and feature scale models," and the SURFACE REACTION MECHANISMS thesis
(D. Zhang, Univ. Michigan). Modern umbrella: multiscale review, *Appl. Phys. Rev.* **8**, 041305
(2021) [VERIFY DOI 10.1063/5.0058196].

- **State variables:** a steady-state fluorocarbon **layer thickness** obtained numerically from a
  *transient multistep site balance* of FC deposition and consumption, plus site fractions of
  reactive surface complexes. Polymerization = neutral sticking + low-energy-ion-assisted
  deposition; consumption = ion sputtering + F-atom etching + O removal.
- **Role assigned to the polymer:** "regulates the etch by limiting the availability of activation
  energy and reactants, and providing the fuel for removal of oxygen" — i.e. the film both
  *attenuates ion energy to the interface* and *supplies carbon to scavenge lattice O as CO/COF2.*
- **What it predicted:** abnormal HAR profile behavior (bowing, twisting) once coupled to MCFPM;
  the C2F6/SiO2 rate vs bias/pressure trends.

### 1.3 Standaert–Oehrlein "defluorination-controlled" picture (the mixed-layer quantitative anchor)
- T. E. F. M. Standaert, M. Schaepkens, N. R. Rueger, P. G. M. Sebel, G. S. Oehrlein, J. M. Cook,
  "High density fluorocarbon etching of silicon in an inductively coupled plasma: Mechanism of
  etching through a thick steady state fluorocarbon layer," *J. Vac. Sci. Technol. A* **16**, 239
  (1998), DOI 10.1116/1.580978.
- T. E. F. M. Standaert et al., "Role of fluorocarbon film formation in the etching of silicon,
  silicon dioxide, silicon nitride, and amorphous hydrogenated silicon carbide," *J. Vac. Sci.
  Technol. A* **22**, 53 (2004), DOI 10.1116/1.1614270 [VERIFY].
- M. Schaepkens & G. S. Oehrlein, "Selective etching of SiO2 over polycrystalline silicon using
  CHF3 in an inductively coupled plasma reactor," *J. Vac. Sci. Technol. A* **17**, 2492 (1999),
  DOI 10.1116/1.581991 [VERIFY]; angular companion, *J. Vac. Sci. Technol. A* **16**, 3281 (1998)
  [VERIFY DOI].

Key **quantitative closures** (these are the ones petch should reproduce, not fit):
1. A **thin steady-state FC film** (few Å to a few nm) is *always present* during etching; it is a
   **fluorine source**, not merely an inhibitor.
2. **Etch rate is controlled by ion-induced defluorination of that film and by the energy that
   reaches the film/substrate interface**, not by film thickness per se. Empirically the
   **substrate etch rate is inversely related to the steady-state FC layer thickness**
   (Oehrlein's inverse-thickness / diffusion picture: ER ∝ 1/d, or ER ∝ exp(−d/λ) with the ion
   energy attenuated through the film) — this is the single most transferable phenomenology.
3. **Selectivity (SiO2 over Si/Si3N4)** arises because the film is *thicker on Si* (no lattice O to
   consume carbon) → less energy to the interface → lower Si rate, while SiO2's lattice O keeps the
   film thin by forming CO/COF2. Selectivity is thus an *emergent* film-thickness difference.

### 1.4 Coburn–Winters reviews (the mechanism ground truth)
- J. W. Coburn & H. F. Winters, "Ion- and electron-assisted gas-surface chemistry—An important
  effect in plasma etching," *J. Appl. Phys.* **50**, 3189 (1979), DOI 10.1063/1.326355.
- J. W. Coburn & H. F. Winters, "Plasma etching—A discussion of mechanisms," *J. Vac. Sci.
  Technol.* **16**, 391 (1979), DOI 10.1116/1.569958 [VERIFY].
- H. F. Winters & J. W. Coburn, "Surface science aspects of etching reactions," *Surf. Sci. Rep.*
  **14**, 161 (1992), DOI 10.1016/0167-5729(92)90009-Z.
Establish: ion energy deposited + mixing *is* the rate multiplier (the Ar+/XeF2 synergy); F atoms
etch the CxFy film as a bulk process; O removes carbon. These are the qualitative constraints every
term below must respect.

---

## 2. MD-informed layer models (Graves lineage) — extractable closures

- **Barone & Graves** (yield-vs-deposited-energy, mixed-layer thickness):
  M. E. Barone & D. B. Graves, "Molecular-dynamics simulations of direct reactive ion etching of
  silicon by fluorine and chlorine," *J. Appl. Phys.* **78**, 6604 (1995), DOI 10.1063/1.360681
  [VERIFY]; and "Chemical and physical sputtering of fluorinated silicon," *J. Appl. Phys.* **77**,
  1263 (1995), DOI 10.1063/1.358928 [VERIFY].
  - Extractable: F+/Cl+ 10–50 eV normal incidence → **quasi-steady halogenated Si layer** forms;
    at 10 eV ≈ **2 equivalent monolayers**, penetration ≈ **15 Å**; roughness ~1–2 nm at 25–50 eV.
    Etch **yield increases with ion energy** following the threshold-`√E` shape. → sets `L_mix`
    (mixed-layer depth) scaling and confirms the `√E`-threshold yield the stopping module derives.
- **Humbird & Graves** (Ar+-induced F transport through FC — the *mixing* term):
  D. Humbird & D. B. Graves, "Improved interatomic potentials for silicon-fluorine and
  silicon-chlorine chemistry and the simulation of silicon etching by fluorine," *J. Chem. Phys.*
  **120**, 2405 (2004), DOI 10.1063/1.1637035 [VERIFY]; and MD of Ar+-induced transport of F
  through fluorocarbon films (APL, 2004) [VERIFY].
  - Extractable: **F is transported from the FC film into the SiFx/Si layer at a rate proportional
    to the ion flux and set by the ion's energy deposition** — a direct, quantifiable *mixing flux*
    J_mix ∝ Φ_ion × (nuclear energy deposited in layer). This is the physical term the ZBL module
    was built to supply.
- **Vegh, Humbird & Graves** (FC-film-mediated steady Si etch):
  J. Vegh, D. Humbird, D. B. Graves, "Silicon etch by fluorocarbon and argon plasmas in the
  presence of fluorocarbon films," *J. Vac. Sci. Technol. A* **23**, 1598 (2005),
  DOI 10.1116/1.2049303 [VERIFY]; cluster-ejection follow-up, *J. Vac. Sci. Technol. A* **26**, 52
  (2008) [VERIFY].
  - Extractable: **steady Si etch requires the FC film thickness to fluctuate**; F must be driven
    through the film by ions; removal is partly as **CFx/SiFx cluster ejection**, not only SiF4.
    Confirms the two-reservoir (film + mixed layer) topology and the ion-gated coupling between them.

**Net MD closures for v1:** (a) yield ∝ nuclear energy deposited in the top ~1.5 nm with a `√E`
threshold near 10–35 eV; (b) F-mixing flux from film→layer ∝ ion flux × deposited energy;
(c) mixed-layer depth `L_mix` ≈ ion projected range in the (fluorinated) substrate, ~1–2 nm at
100s eV, growing with E; (d) volatile products SiF4 (Si) and CO/COF2 (C+O). Every one of (a)–(c) is
an integral of the ZBL curves already in `ion_energy_deposition.py` — no new fitted energy law.

---

## 3. Element-resolved / element-conserving reduced surface models (prior art for petch's ledgers)

Petch's differentiator is **exact element ledgers**. Two published lines already track *elements*,
not just coverages, and are the direct antecedents:

- **Wang & Kushner implantation-and-mixing model:**
  M. Wang & M. J. Kushner, "Modeling of implantation and mixing damage during etching of SiO2 over
  Si in fluorocarbon plasmas," *J. Vac. Sci. Technol. A* **29**, 051306 (2011),
  DOI 10.1116/1.3231450 [VERIFY]. Implantation + mixing built into an MCFPM: energetic ions carry
  F/C into the lattice and mix Si/O/C/F across a near-surface region during HAR SiO2 etch. This is
  the *published* precedent for "ion-induced mixing feeds a resolved near-surface composition."
- **Kuboi voxel-slab reaction-layer model (Sony) — the closest structural analog:**
  N. Kuboi et al., "Insights into different etching properties of continuous wave and atomic layer
  etching processes for SiO2 and Si3N4 films using voxel-slab model," *J. Vac. Sci. Technol. A*
  **37**, 051004 (2019), DOI 10.1116/1.5090606 [VERIFY]; review, IOP *Jpn. J. Appl. Phys.* (2024),
  DOI 10.35848/1347-4065/ad5355 [VERIFY].
  - **Surface = two stacked layers: (i) a C–F polymer layer, (ii) a reactive layer subdivided into
    ~0.25 nm slabs.** In the reactive layer, **ion bombardment generates Si dangling bonds that
    react with F to make SiF2/SiF4 by-products, branching on incident ion energy.** This is almost
    exactly the petch v1 topology (film reservoir + depth-resolved mixed layer) — but Kuboi's slab
    fractions are updated by tabulated MD probabilities, *not* by a closed element ledger.
- **Chang et al. unified semi-global surface reaction model:**
  W. S. Chang et al., "A unified semi-global surface reaction model of polymer deposition and SiO2
  etching in fluorocarbon plasma," *Appl. Surf. Sci.* **515**, 145975 (2020),
  DOI 10.1016/j.apsusc.2020.145975 [VERIFY]. Describes **deposition→etch transition self-
  consistently without ad-hoc assumptions**, and explains HAR bowing/twisting — a phenomenological
  target for "clog/no-clog emergent."

**Answer to Q3:** No published model enforces a *closed, exact* per-element budget (separate F and C
budgets, F entering as CxFy and leaving as SiF4/COF2) as a conservation invariant. Kuboi and
Wang–Kushner resolve composition but *advance* it with probabilistic MD rates; conservation is
approximate/statistical. **An element-conserving reduced surface model with exact F/C/O/Si ledgers
is a genuine gap, and it is the natural extension of petch's engine** — the ledger becomes the
model's backbone rather than a diagnostic.

---

## 4. Oxygen pathways (how O consumes the film / defluorinates; model treatments; saturation)

Mechanism (Coburn–Winters, Standaert, Kushner SKM):
- O (atomic O and O2-dissociation products) at the surface **oxidizes the carbon of the FC film to
  volatile CO / CO2 / COF2**, thinning the film. Thinner film → more ion energy to the interface →
  higher SiO2 etch rate. This is the additive that keeps the film thin and prevents clog.
- **Branching:** C leaves as **CO / CO2** (carbon + lattice/gas O) or as **COF2** (carbon + O +
  2 F) — the COF2 channel *also removes fluorine*, so at high O the film is both defluorinated and
  decarbonized. Weak-oxidant chemistry (CO/CO2 additions form COF2, scavenging F) is the gas-phase
  mirror of the same C–O–F coupling.
- **Saturation (Langmuir-type):** O adsorption/coverage saturates; once the FC film is fully
  consumed, additional O (a) over-oxidizes the substrate and (b) increasingly consumes F as COF2,
  so etch rate **rises then plateaus/rolls over** with O2 fraction. Published treatments carry an
  O-site coverage θ_O with Langmuir isotherm and an O-driven polymer-removal term (Krüger's
  O-based polymer-etch probability 0.0628 is the lumped version of exactly this term).

For petch, oxygen is **not** a separate empirical knob: it is the O element in the ledger. The
"oxygen saturation curve" must fall out of (O-supply vs C-available-to-oxidize) competition, not a
fitted saturation constant. This is validation target Rung 2 (Krüger O2 sweep).

---

## 5. Design verdict — petch mixed-layer v1

### 5.1 Minimal state vector (per surface patch)
Two stacked, element-resolved reservoirs (Kuboi topology, but as *closed budgets*):

**(A) Fluorocarbon film reservoir** — areal densities [atoms·m⁻²]:
  `n_C^film`, `n_F^film`.  Film thickness `d_FC = (n_C^film + n_F^film)/ρ_FC` with
  `ρ_FC = 7.5e28 m⁻³` (already the `FLUOROCARBON_FILM` density in `ion_energy_deposition.py`).

**(B) Mixed reaction layer** — top `L_mix` of the substrate, areal densities [atoms·m⁻²]:
  `n_Si`, `n_O`, `n_C`, `n_F`.  `L_mix = R_p(E_ion; SIO2)` from
  `projected_range_nm(E, z1, m1, SIO2)` (E-dependent, ~1–2 nm; **not fitted**).

That's **6 scalars per patch** (`n_C^film, n_F^film, n_Si, n_O, n_C, n_F`). The SF6/O2 base
chemistry is the degenerate case `n_C^film=n_F^film=n_C=0` → a pure Si/F/(O) layer, which must
reproduce Belen/de Boer exactly (Rung 0 regression).

### 5.2 Governing equations (all rates areal, per timestep; every term traceable)
Incoming fluxes: `J_CFx` (deposition precursor, F/C ratio `y`≈1.5 for C4F6), `J_F` (atomic F),
`J_O` (atomic O/O2 products), `J_ion(E,θ)`.

Ion energy actually reaching the film/substrate interface (Standaert defluorination law):
```
E_iface = E_ion · exp(−d_FC / λ_FC),   λ_FC = R_p(E_ion; FLUOROCARBON_FILM)   # from stopping module
```
Nuclear energy deposited in the mixed layer (the yield/mixing driver — no free knob):
```
ε_dep = nuclear_energy_in_layer_eV(E_iface, cosθ, L_mix, z1, m1, SIO2)        # ion_energy_deposition.py
```

Film reservoir:
```
dn_C^film/dt = s_p·J_CFx           − Y_sput(ε_dep_top)·J_ion·x_C^film  − k_O·n_O_surf·n_C^film
dn_F^film/dt = s_p·y·J_CFx + J_F·(1−θ_F) − Y_sput(ε_dep_top)·J_ion·x_F^film − 2·k_O·n_O_surf·n_C^film − J_mix
```
- `s_p` sticking (Gogolides/Kushner; Krüger 0.0842). `Y_sput ∝ ε_dep_top` (ion sputter of film,
  from stopping module, reference-anchored). `k_O·n_O·n_C^film` = O-oxidation of film carbon
  (→CO/COF2); the `2·k_O` on F is the COF2 F-loss branch. `J_mix` = F driven into the mixed layer.

Ion-induced F-mixing flux (Humbird–Graves):
```
J_mix = η_mix · J_ion · (ε_dep / E_ref) · n_F^film/(n_C^film+n_F^film)      # η_mix from MD, O(1)
```

Mixed layer (per element; SiF4 and COF2/CO stoichiometry fixed → ledger closes by construction):
```
R_SiF4 = k_v · g(n_F, n_Si) · (ε_dep/E_ref)          # Si volatilization, √E-threshold via ε_dep
dn_F/dt  = J_mix + J_F·θ_penetrate − 4·R_SiF4 − 2·R_COF2
dn_Si/dt = ρ_SiO2·v_recession/3 (Si liberated by recession) − R_SiF4        # SiO2: 1 Si : 2 O
dn_O/dt  = 2·ρ_SiO2·v_recession/3 − R_CO − R_CO2 − R_COF2 + J_O·θ_O
dn_C/dt  = (C mixed from film) − R_CO − R_CO2 − R_COF2
```
Recession velocity (front motion) from the Si leaving as SiF4:
```
v_recession = R_SiF4 / n_SiO2^formula        # n_SiO2^formula = 2.2e28 m⁻³ (SiO2 formula units)
```

**Exact ledgers enforced each step (the invariant, not a diagnostic):**
```
F:  s_p·y·J_CFx + J_F         = 4·R_SiF4 + 2·R_COF2 + (dn_F^film+dn_F)/dt      (+HF if enabled)
C:  s_p·J_CFx                 = R_CO + R_CO2 + R_COF2 + C_sputtered + dn_C^film/dt + dn_C/dt
O:  J_O + O_from_lattice      = R_CO + 2·R_CO2 + R_COF2 + dn_O/dt
Si: Si_from_lattice(recession)= R_SiF4 + Si_redeposited
```
Any residual is a bug, not a knob — this is petch's exact-conservation contract applied to chemistry.

### 5.3 Why the three target behaviors are *emergent*, not fitted
- **(a) Ion-capacity-vs-neutral-supply ceiling cliff:** `v_recession = min(ion-driven SiF4 capacity,
  F actually supplied)`. Raising `J_CFx` thickens `d_FC` → `E_iface = E_ion·exp(−d_FC/λ)` collapses
  → recession → 0 (deposition/clog). Lowering neutral F starves `R_SiF4`. The cliff is the crossing
  of these two limbs; nothing sets its location by hand — it moves with `J_CFx/J_ion` and `E_ion`.
- **(b) Oxygen saturation emergent:** O enters only as the `k_O·n_O·n_C^film` and `J_O·θ_O` terms
  and the COF2/CO product branches. Rate rises while O thins the film (more `E_iface`), then rolls
  over once C is exhausted and added O consumes F as COF2. The saturation constant is *derived* from
  C availability, not fitted.
- **(c) Clog/no-clog emergent:** steady state `dn^film/dt = 0` has a solution `d_FC*<∞` (etch) or no
  finite solution (`d_FC→∞`, clog) depending on whether deposition exceeds the ion+O removal that
  can reach the film top. The boundary in `(J_CFx/J_ion, E_ion, J_O)` space is the model output.

### 5.4 Parameter provenance table (every parameter → literature or the stopping module)
| Parameter | Meaning | Source (traceable) |
|---|---|---|
| `λ_FC`, `L_mix` | ion range in FC / SiO2 | `ion_energy_deposition.projected_range_nm` (ZBL, **no free param**) |
| `ε_dep`, `Y_sput`, `R_SiF4(E)` shape | energy deposited → yield | `nuclear_energy_in_layer_eV` (ZBL/Sigmund); Barone–Graves `√E`-threshold |
| `η_mix` | F mixing efficiency | Humbird–Graves MD (Ar+-induced F transport), O(1) |
| `s_p`, `y` | FC sticking, precursor F/C | Gogolides/Kokkoris, Kushner SKM; Krüger 0.0842; gas chemistry |
| `k_O` + COF2/CO branch | O oxidation of film C | Standaert 2004; Kushner SKM; Krüger O-etch 0.0628 |
| `k_v`, `g(n_F,n_Si)` | SiF4 volatilization | Barone–Graves MD; anchored to base calibration magnitude |
| SiF4 / COF2 / CO stoichiometry | product elements | **fixed by chemistry — zero free params** (this is the ledger) |
| `ρ_FC`, `n_SiO2^formula` | densities | `ion_energy_deposition.py` (`7.5e28`, `6.6e28` atomic) |

Free-parameter count vs today: the *energy law* (threshold, exponent, magnitude shape) becomes
**derived** (was 3 fitted numbers per channel in Krüger/Belen); the surviving fitted numbers are the
sticking `s_p`, the oxidation rate `k_O`, the mixing efficiency `η_mix`, and the volatilization
prefactor `k_v` — each with a named literature anchor and a stated uncertainty, and all magnitudes
reference-anchored to the existing base calibration so v1 cannot silently move a validated number.

### 5.5 Validation ladder (mapped to datasets already in the repo)
| Rung | Dataset (in repo) | Target | Gate |
|---|---|---|---|
| 0 (regression) | `deboer_2002/` + Belen base | degenerate F-only layer == current SF6/O2 | de Boer Fig 9 evolving RMSE ≤ 0.05; notch gates byte-identical |
| 1 (calibration) | `krueger_2024/base_case_metrics.csv` | C4F6/Ar/O2 SiO2 base (Table IV: hf=825 nm, wt=90, wm=45) | reproduce within existing Krüger replay tolerance |
| 2 (O saturation, held-out) | `krueger_2024/transfer_observations.csv` (O2/C4F6 sweep) | rate-vs-O2 rises then saturates — **emergent** | qualitative monotone-then-rollover; `reference_only`, never scored as experiment |
| 3 (ion-energy / clog, held-out) | `krueger_2024` low-freq power sweep + `digitized_figure16*` | clog↔etch boundary vs ion energy | boundary exists and moves correct direction; `reference_only` |
| 4 (ARDE + ledger) | `deboer_2002/digitized_figure9.csv` | ARDE knee under AR flux attenuation; **F/C/O/Si ledger closes at every depth** | keep existing ARDE agreement; ledger residual < 1e-9 (petch conservation contract) |
| 5 (selectivity, stretch) | Standaert/Schaepkens (external) | SiO2/Si selectivity from `d_FC` difference | qualitative: thicker film on Si, lower Si rate |

Rungs 2–3 respect the repo's held-out contract (`krueger_2024/README.md`): transfer trends are
**calibration-excluded**; the mixed-layer model earns a predictive claim only by making the O2 and
power sweeps *emergent* after being calibrated solely on Rung 0–1. Ledger closure (Rung 4) is the
one gate no coverage model can offer and petch's engine makes free.

---

## 6. Recommendation
Build v1 as the **two-reservoir, element-resolved (C,F,O,Si) mixed layer** of §5.1–5.2, driven by
`ion_energy_deposition.py` for every energy/range/mixing quantity, with the four conservation
ledgers enforced as invariants. It subsumes La Magna–Garozzo and the Krüger reduced mechanism as
lumped limits, turns their fitted energy law into a stopping integral, and makes the ceiling cliff,
oxygen saturation, and clog/no-clog *outputs* rather than inputs. The Kuboi voxel-slab model and the
Wang–Kushner implantation-mixing model are the published structural precedents; the **exact element
ledger is the novel, defensible contribution** — no prior FC surface model closes F and C budgets
exactly.

---

### Source list (markdown)
- Gogolides et al., *J. Appl. Phys.* 88, 5570 (2000) — https://doi.org/10.1063/1.1311808
- Standaert et al., *J. Vac. Sci. Technol. A* 16, 239 (1998) — https://doi.org/10.1116/1.580978
- Standaert et al., *J. Vac. Sci. Technol. A* 22, 53 (2004) — https://pubs.aip.org/avs/jva/article-abstract/22/1/53
- Schaepkens & Oehrlein, *J. Vac. Sci. Technol. A* 17, 2492 (1999) — https://pubs.aip.org/avs/jva/article-abstract/17/5/2492
- Coburn & Winters, *J. Appl. Phys.* 50, 3189 (1979) — https://doi.org/10.1063/1.326355
- Winters & Coburn, *Surf. Sci. Rep.* 14, 161 (1992) — https://doi.org/10.1016/0167-5729(92)90009-Z
- Barone & Graves, *J. Appl. Phys.* 77, 1263 (1995); 78, 6604 (1995) — [VERIFY DOIs]
- Humbird & Graves, *J. Chem. Phys.* 120, 2405 (2004) — [VERIFY DOI 10.1063/1.1637035]
- Vegh, Humbird & Graves, *J. Vac. Sci. Technol. A* 23, 1598 (2005) — https://pubs.aip.org/avs/jva/article-abstract/23/6/1598
- Wang & Kushner, *J. Vac. Sci. Technol. A* 29, 051306 (2011) — https://pubs.aip.org/avs/jva/article-abstract/29/5/051306
- Kuboi et al., *J. Vac. Sci. Technol. A* 37, 051004 (2019) — https://pubs.aip.org/avs/jva/article-abstract/37/5/051004
- Chang et al., *Appl. Surf. Sci.* 515, 145975 (2020) — https://www.sciencedirect.com/science/article/abs/pii/S0169433220307315
- Krüger, *Modeling and Optimization of High Aspect Ratio Plasma Etching* / JVSTA 42, 043008 (2024) — https://doi.org/10.1116/6.0003554
- de Boer et al., *J. Microelectromech. Syst.* 11, 385 (2002) — https://doi.org/10.1109/JMEMS.2002.800928
- La Magna & Garozzo, *J. Electrochem. Soc.* 150, F178 (2003) — https://doi.org/10.1149/1.1602084
