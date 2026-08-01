# Sub-degree IADF physics + Tier-1 reactor completeness — research pass (2026-07-29)

Status: **research document, not built, not committed.** Companion/critique of
`RESEARCH_REACTOR_TIER1_DESIGN_2026-07-24.md`.
`[VERIFY]` marks a locator or number I could not confirm against the primary full text in this pass
(most AIP/AVS/IOP bodies paywall or 403 to automated fetch; abstracts and Crossref metadata were
fetched successfully and every DOI below was resolved through `api.crossref.org` unless flagged).

Every petch number in §A.7 was **measured by running the repo in this session**, not estimated.
Reproduction commands are inline.

---

## 0. Executive summary (read this if nothing else)

1. **The thermal ion temperature is not the problem.** petch's default
   `ion_tangential_temperature_eV = 0.026` is only ~1.3× narrower than the measured value
   (0.044 eV, Kim *et al.* 2025). The real deficit is **structural**: the collisionless sheath model
   produces a *single-Gaussian, purely thermal* IADF, while the measured and HPEM-computed IADFs at
   HAR conditions are **bi-Gaussian — a thermal core plus a much wider non-thermal tail produced by
   ion–neutral collisions inside the sheath**, and at AR ≥ 100 the tail carries essentially all of
   the sidewall flux.
2. **Measured, at HAR-relevant energy, with better-than-0.1° resolution, this decomposition is:**
   core T⊥ = 0.044 eV, tail T⊥ = 0.57 eV (Kim *et al.*, JJAP **64**, 05SP15, 2025), Ar, 2.4 Pa,
   V_pp = 2.7 kV, V_dc = 950 V, E_ion = 1.4–2.0 keV. Sub-degree structure **is real, is measured,
   and is predictable** — but only from a *collisional* sheath model.
3. **Measured petch gap (this session):** driving `build_diagnostic_virtual_sheath_boundary` at
   Krüger-like conditions gives σ_θ(1-D) = **0.148°**; the repo's own digitized Krüger Fig-4 HPEM
   IEAD gives σ_θ(1-D) = **0.833°** (repo's own digitization band 0.822–0.860°). That is **5.6×
   too narrow**. Translated to AR 200: petch's virtual sheath puts **99%** of the ion flux inside
   the acceptance cone; the Krüger reference puts **28%** there. The sidewall flux is under-delivered
   by **~70×** on that path (~30× against the conservative measured bi-Gaussian).
   This composes multiplicatively with the known ~1.5× azimuth-compression deficit.
4. **The whole angular chain in petch is quantized at 0.25°**, which *equals* the AR-200 acceptance
   half-angle (0.286°). Boundary source (digitized Fig-4 bins = 0.25°), production transport
   quadrature (`joint ion-IEAD bins of 250 eV by 0.25 degrees`,
   `KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md`), and the 3-node Gauss–Hermite sheath
   representation all sit at the same scale. **AR 200 needs ≈0.05°.**
5. **Biggest Tier-1 design-doc gap:** the design doc explicitly *defers* the gas-temperature energy
   balance and *assumes* a collisionless sheath. Both are load-bearing for sub-degree IADF:
   T_gas **is** the core width (θ_core = √(kT_gas/E)), and sheath collisions **are** the tail.
   The doc also has **no fast-neutral (NEAD) species** in the boundary contract, although Krüger's
   own MCFPM treats hot neutrals as a first-class energetic flux and the measured neutral angular
   width is *wider* than the ion width at every condition.

---

# A. Sub-degree IADF physics

## A.1 The governing formula, and its two exact forms

Ions cross a (planar, collisionless) sheath under a purely normal field. Transverse velocity is
conserved; normal velocity is set by the sheath drop. Therefore

> tan θ = v⊥ / v∥ , with v∥ = √(2eV_sh/M) and v⊥ drawn from the sheath-edge transverse distribution.

For a Maxwellian transverse population at temperature T⊥ and a normal energy E:

- **2-D radial form** (what Khrabrov & Kaganovich write): f(θ) ∝ θ·exp(−θ² E/T⊥); the characteristic
  thermal spread is
  > **θ_th = √(T⊥ / E_b)  rad**
  — stated verbatim in Khrabrov & Kaganovich, arXiv:2604.04214v2 (2026), §6.2: *"The differential
  cross-section only needs to be accurately known down to the angle on the order of the thermal
  spread θ_th = √(T⊥/E_b) rad or approximately 0.4° in the case at hand"* (their case: T⊥ = 0.044 eV,
  E_b = 1 keV → 0.38°). **Verified quote.**
- **1-D projected form** (what a signed-angle IEAD plot such as Krüger Fig. 4 shows) is a Gaussian
  with
  > **σ_θ(1-D) = √( T⊥ / (2E) )  rad = θ_th/√2**

Both are mass-independent (M cancels). This is confirmed experimentally: Ar⁺ and Kr⁺ show
*"similar main-component angular widths"* attributable to thermal motion (Kim *et al.* 2025b).

**Consequences for HAR.** The geometric acceptance half-angle of a straight feature is
arctan(1/AR):

| AR | acceptance half-angle | σ needed for 5 bins across the half-angle |
|---|---|---|
| 30 | 1.909° | 0.382° |
| 50 | 1.146° | 0.229° |
| 100 | 0.573° | 0.115° |
| **200** | **0.286°** | **0.057°** |

At 3465 eV (Krüger base-case mean ion energy) the *thermal core alone* is σ = 0.144° for
T⊥ = 0.044 eV. **The thermal core is already ~half the AR-200 acceptance half-angle.** There is no
regime at AR 200 where the IADF can be treated as a delta function.

## A.2 What sets T⊥ — and why it is not a knob

The measured core temperature, 0.044–0.045 eV, corresponds to **510–520 K**
(kT = 0.044 eV ⇒ T = 0.044/8.6173e-5 = 511 K). That is the **neutral gas temperature** of an
industrial HAR CCP. Three independent lines say T⊥,core *is* T_gas, not a fit parameter:

- **Kinematic argument.** The sheath field is normal; it does no work on v⊥. Absent collisions, the
  transverse distribution arriving at the wafer is the transverse distribution at the sheath edge,
  which is thermalised to the neutral gas by charge exchange in the bulk/presheath.
- **Measurement.** Kim *et al.* (JJAP **64**, 05SP15, 2025) report *effective ion temperature (main
  component) = 0.044 eV*, and Khrabrov & Kaganovich record that *"for T⊥ = 0.045 eV reported by the
  authors as the value consistently seen for Ar⁺ ions in all cases"* — i.e. **energy- and
  condition-independent**, the signature of a thermodynamic property rather than a fitted width.
- **Transport-coefficient route (knob-elimination path).** In a drift field, the transverse
  temperature is fixed by the measured ratio of transverse diffusion to mobility,
  kT⊥ = e·D_T/μ (generalised Einstein relation), which for Ar⁺ in Ar is *measured* over the full
  E/N range: Stefánsson & Skullerud, *J. Phys. B* **32**, 1057 (1999),
  DOI [10.1088/0953-4075/32/5/001](https://doi.org/10.1088/0953-4075/32/5/001). Wannier-type
  field heating is therefore an **evaluable** correction, not a free parameter.

**Numerical sensitivity — why T_gas is load-bearing.** At E = 3465 eV, thermal core only:

| T_gas | kT (eV) | σ_core | P(\|θ\| < 0.286°) at AR 200 |
|---|---|---|---|
| 300 K | 0.0259 | 0.111° | 0.990 |
| 500 K | 0.0431 | 0.143° | 0.955 |
| 800 K | 0.0689 | 0.181° | 0.887 |
| 1000 K | 0.0862 | 0.202° | 0.844 |

A 300 K ↔ 1000 K gas-temperature ambiguity moves the AR-200 bottom-arrival fraction from 0.99 to
0.84, i.e. it changes the *sidewall* flux by **~16×** (0.010 → 0.156). The Tier-1 design doc lists
"T_g energy balance" under **Explicitly deferred**. For sub-degree IADF that deferral is not
tolerable — see §B.1.

## A.3 The dominant term: sheath collisions and the non-thermal tail

**This is the finding that changes the build order.** The measured IADF at HAR conditions is not
one Gaussian. Kim *et al.* fit a **bi-Gaussian**: a narrow main (thermal) component plus a weak,
much broader tail.

- **Measured decomposition** (Kim, Kawamura, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda,
  *"Measurement of energy-resolved ion angular distribution in a dual-frequency capacitively coupled
  argon plasma"*, **Jpn. J. Appl. Phys. 64, 05SP15 (2025)**,
  DOI [10.35848/1347-4065/adce84](https://doi.org/10.35848/1347-4065/adce84)):
  - effective ion temperature, **main component: 0.044 eV**
  - **tail component: 0.57 eV** (13× hotter ⇒ 3.6× wider)
  - Ar, **2.4 Pa (18 mTorr)**, V_pp = 2.7 kV, **V_dc = 950 V**, ion energies **1.4–2.0 keV**
  - orifice→MCP drift 29 cm, drift chamber < 5×10⁻⁵ Pa, ICCD gate 100 ns (~150 eV energy resolution
    at 2 keV), bi-Gaussian parameters reconstructed by gradient descent, normalised error < 0.015
  - main-component angular width decreases monotonically (≈logarithmically) with incident energy —
    the √(1/E) law of §A.1
- **Mechanism identified** (Kim, Kawamura, Fujitani, Naito, Iino, Fukumizu, Kurihara, Suzuki
  [+ Toyoda], *"Influence of sheath collisions on the ion angular distributions in a dual-frequency
  capacitively coupled plasma"*, **Jpn. J. Appl. Phys. 64, 096002 (2025)**,
  DOI [10.35848/1347-4065/ae0105](https://doi.org/10.35848/1347-4065/ae0105)):
  tail components originate from *"Ar⁺ collisions with Ar atoms in the sheath"*; the
  **main-to-total component ratio at the maximum energy decreases exponentially with pressure** —
  a Beer–Lambert survival law exp(−s/λ) in the sheath thickness *s* over the ion–neutral mean free
  path λ. Ar and Kr have the *same* main-component width (thermal) but pressure-dependent tails
  (collisional). **[VERIFY] the printed tail-fraction values and the fitted λ — full text paywalled.**
- **Cross-section physics** (Khrabrov & Kaganovich, arXiv:2604.04214v2, PPPL, Apr 2026; supported by
  Samsung + DOE/PPPL CRADA with Applied Materials): the tail is **anisotropic small-angle elastic
  scattering** governed by the Born–Mayer repulsive pair potential. Two model-quality findings that
  matter for us:
  - The widely used Phelps (1994) sheath ion-flux cross-section model
    (*J. Appl. Phys.* **76**, 747, DOI [10.1063/1.357820](https://doi.org/10.1063/1.357820))
    is *"essentially an identity-switch charge exchange"* at 1 keV and **"fails to predict appreciable
    elastic scattering of ions"** — i.e. it produces *no* ion angular tail. Nanbu–Kitatani likewise.
    **A petch sheath collision operator built on Phelps-1994 alone would reproduce the CX energy
    tail but none of the angular tail.**
  - The Phelps–Greene–Burke (2000) Ar⁰–Ar⁰ model
    (*J. Phys. B* **33**, 2965, DOI [10.1088/0953-4075/33/16/303](https://doi.org/10.1088/0953-4075/33/16/303))
    over-broadens: *"a half-width of about 1°… the number of particles scattered with 1° < θ < 5° is
    approximately twice the number for the Born-Mayer model"*, because a screened-Coulomb Born DSC
    fitted to total/momentum/viscosity cross-sections gets an artificially high screening angle.
  - Verdict quoted verbatim: *"the scattering of both ion and neutral species is anisotropic at all
    energies down to just above the thermal range and needs to be treated carefully in order to
    obtain meaningful results."*

**Is the Krüger base case actually collisional? Yes.** Order-of-magnitude with explicit inputs
(σ_CX(Ar⁺/Ar) ≈ 2–3 × 10⁻¹⁹ m² at keV **[VERIFY exact value against Phelps 1994 Fig.]**, n = p/kT_gas):

| p | T_gas | λ | sheath s | s/λ | P(no collision) |
|---|---|---|---|---|---|
| 10 mTorr | 400 K | 13.8–20.7 mm | 8.7 mm (petch Child law, this session) | 0.42–0.63 | 0.53–0.66 |
| 10 mTorr | 400 K | 13.8–20.7 mm | 14 mm | 0.68–1.01 | 0.36–0.51 |
| 10 mTorr | 600 K | 20.7–31.1 mm | 8.7–14 mm | 0.28–0.68 | 0.51–0.76 |

Krüger's own sheath widths are centimetre-scale: thesis §3, *"In pure argon, the sheath thickness
during the cathodic phase is 9.5 mm… For Ar/O₂ = 50/50, the sheath thickness is 15 mm"*
(`tmp/pdfs/krueger_thesis.txt` L2798–2800). Independent literature value quoted in search:
*"ion mean free path at 10, 20, and 40 mTorr is 10, 5, and 2.5 mm"* **[VERIFY source]**.
**Conclusion: 30–65% of ions charge-exchange or scatter inside the sheath in Krüger's base case.
"Collisionless" is a ~50% error on the population that sets sidewall flux.**

**Independent consistency check performed this session.** Fitting a bi-Gaussian to the marginal
angular distribution of the repo's `digitized_figure4_iead.csv` returns a **tail fraction of 0.65**.
The collided fraction 1−exp(−s/λ) for s = 14 mm, λ = 13.8 mm is **0.64**. The tail weight in
Krüger's own HPEM output is quantitatively what sheath collisionality predicts. (The fitted widths,
σ_core = 0.600°, σ_tail = 0.946°, are *upper bounds* — the digitization grid is 0.25° so the core
cannot be resolved. σ_tail = 0.946° is ~1.8× the 0.520° implied by the measured T_tail = 0.57 eV at
3465 eV; the residual is plausibly log-colorbar digitization over-weighting the far tail, or genuine
extra HPEM broadening. **Open.**)

## A.4 RF phase modulation, and the energy–angle correlation

Two effects, both already partly present in petch and both second-order compared with §A.3:

- **Finite ion transit (ωτ_ion).** Krüger's own thesis defines
  S = (d_s/2)·√(2m_i/(qV_S))·2π f_RF (Eq. 5.1) — S > 1 ⇒ ions do not respond to intra-cycle field
  change and arrive near the mean sheath potential; S < 1 ⇒ dynamic response and the bimodal
  (Benoit-Cattin) split. He finds **S > 1 at f₀ = 10 MHz and S < 1 at f₀ = 1 MHz** for his geometry.
  petch's `CollisionlessWaveformSheath` integrates this exactly (velocity-Verlet on a Child profile)
  — this part is **sound and does not need work**.
- **Energy–angle correlation.** Because θ ∝ 1/√E and E is phase-dependent, the IADF is *not*
  separable: the low-energy horn is systematically wider. petch's tensor construction in
  `collisionless_sheath_boundary_state` (`src/petch/boundary_state.py:853`) *does* get this right —
  each phase node's energy is paired with the same transverse node set, so θ = atan(v⊥/√E) varies
  per phase. **Keep this structure when replacing the quadrature.**
- **Tangential drift → a shifted, azimuthally-asymmetric IADF.** Raja & Linne,
  *"Analytical model for ion angular distribution functions at rf biased surfaces with collisionless
  plasma sheaths"*, **J. Appl. Phys. 92, 7032–7040 (2002)**,
  DOI [10.1063/1.1524020](https://doi.org/10.1063/1.1524020) (note: the design doc's ledger gives
  10.1063/1.1519941 — **that DOI is wrong**; correct is 10.1063/1.1524020, authors Raja & Linne,
  Lund/Colorado). Their result: *"Tangential ion drift velocities introduce azimuthal angle
  dependence on the IADF and a shift in the peak IADF to off-normal polar angles… the shift in peak
  IADF in the polar angle depends on both the drift velocity as well as the bias frequency"*, and
  under DC bias the IADF shape depends strongly on bias voltage and ion temperature but is *"relatively
  independent of the plasma electron temperature, ion density, and the ion mass"*. **This is the
  physical bridge to the azimuth-compression problem: a real IADF is not azimuthally isotropic.**

## A.5 Measured IADFs — the anchor table (with angular resolution)

| Source | Reactor / gas | E_ion | Angular resolution | Reported width | T⊥ |
|---|---|---|---|---|---|
| Woodworth, Riley, Meister, Aragon, Le, Sawin, **J. Appl. Phys. 80, 1304–1311 (1996)**, DOI [10.1063/1.362977](https://doi.org/10.1063/1.362977) | ICP Ar, GEC cell, 2.4–50 mTorr | 12–29 eV (plasma potential; unbiased) | not sub-degree | half-widths **5°–12.5°** | **0.2–0.5 eV** (broadens with p and RF power) |
| Woodworth, Riley, Miller, Hebner, Hamilton, **J. Appl. Phys. 81, 5950–5959 (1997)**, DOI [10.1063/1.364383](https://doi.org/10.1063/1.364383) | ICP Cl₂, 20–60 mTorr, grounded electrode | 9–13 eV | not sub-degree | half-widths **6°–7.5°** | **0.13–0.21 eV** |
| Sharma, Gahan, Scullin, Daniels, Hopkins, **Rev. Sci. Instrum. 86, 113501 (2015)**, DOI [10.1063/1.4934808](https://doi.org/10.1063/1.4934808) | planar RFA, variable-aspect-ratio aperture | — | **~3° (best)** | — | — |
| Ichikawa, Chu, Moriyama, Nakahara, Suzuki, Iino, Fukumizu, Kurihara [+ Toyoda], **Appl. Phys. Express 14, 126001 (2021)**, DOI [10.35848/1882-0786/ac33c4](https://doi.org/10.35848/1882-0786/ac33c4) | 13.56 MHz CCP Ar, MCP 2-D beam imaging | self-bias > 300 V | sub-degree (MCP) | **< 1°**, monotonically narrowing with V_dc; **neutral width > ion width at every condition** | — |
| **Kim, Kawamura, Naito, Iino, Fukumizu, Kurihara, Suzuki, Toyoda, Jpn. J. Appl. Phys. 64, 05SP15 (2025)**, DOI [10.35848/1347-4065/adce84](https://doi.org/10.35848/1347-4065/adce84) | dual-freq CCP Ar, 2.4 Pa, V_pp 2.7 kV, V_dc 950 V | **1.4–2.0 keV** | **better than 0.1°** | bi-Gaussian | **core 0.044 eV / tail 0.57 eV** |
| Kim *et al.*, **Jpn. J. Appl. Phys. 64, 096002 (2025)**, DOI [10.35848/1347-4065/ae0105](https://doi.org/10.35848/1347-4065/ae0105) | as above, Ar and Ar/Kr, pressure scan | keV | better than 0.1° | tail = sheath collisions; **main/total ratio falls exponentially with p** | — |
| Kawamura, Kim, Ichikawa, Suzuki, Iino, Fukumizu, Kurihara, Toyoda, **Plasma Sources Sci. Technol. 34, 055006 (2025)**, DOI [10.1088/1361-6595/add321](https://doi.org/10.1088/1361-6595/add321) | 40.7/2 MHz dual-freq CCP Ar | up to **~2.5 keV** | better than 0.1° | ion **and neutral** widths ∝ √(T_i/E); both broaden with pressure, strongest at high V_pp | — |

**Answer to "is sub-degree structure real and predictable?"** — *Real*: yes, resolved better than
0.1° by MCP 2-D beam imaging at industrially relevant keV energies (Nagoya/Toyoda group, 2021–2025).
*Predictable from sheath models*: **the thermal core, yes** (θ = √(T_gas/E), zero free parameters);
**the tail, only from a collisional sheath with an anisotropic differential cross-section** — and
the two most-used cross-section models bracket the truth badly (Phelps-1994 gives no ion tail;
Phelps–Greene–Burke-2000 gives ~2× too much 1°–5° scattering). Khrabrov & Kaganovich's Born–Mayer
model still *"shows a steeper fall-off than the experimentally measured one"* even at L/λ_cx = 4.
**The wide-angle tail is a genuine open frontier — a place where petch can be first, not last.**

Beware Woodworth's T⊥ = 0.2–0.5 eV: those are **single-component fits at 12–29 eV** where core and
tail overlap; they are not comparable to the 0.044 eV *core* of a bi-Gaussian fit at 2 keV. Also
relevant: presheath LIF shows ion heating from ~room temperature in the bulk to **0.13 eV 0.5 cm
from the plate** **[VERIFY primary source; from search summary of LIF presheath literature]**, and
ion-acoustic instability heating in presheaths is documented. Do **not** propagate Woodworth's
0.2–0.5 eV into a keV HAR IADF; do use it as the upper bound for low-energy ICP arms.

## A.6 A larger angular error than the width: mean-angle tilt

For AR 200 the *centroid* of the IADF matters more than its width. Sheath bending at the wafer edge
(focus-ring erosion, height step) produces a **systematic tilt of the mean incidence angle**:

- Seong, Lee, Kim, Lee, Cho, Lee, Jeong, You *et al.*, *"Characterization of an Etch Profile at a
  Wafer Edge in Capacitively Coupled Plasma"*, **Nanomaterials 12, 3963 (2022)**,
  DOI [10.3390/nano12223963](https://doi.org/10.3390/nano12223963): at a wafer/electrode height step
  ≤ 0.6 mm the profile tilt is **≤ 2.5°**, growing with the step; and explicitly *"a 2-degree tilt
  can be a defect"*. Tilt is **independent of contact-hole diameter** — i.e. it is set purely by the
  ion trajectory, not by feature-scale transport.

**2.5° is ~9× the AR-200 acceptance half-angle.** A 0-D global model cannot produce this: it is a
2-D sheath-geometry quantity. Any Tier-1 claim of "recipe → IEAD" is a *center-of-wafer* claim only.
This must be stated in the boundary provenance.

## A.7 petch audit — measured in this session

Reproduce with (from `plasma-etching-code/`):

```
python3 -c "
import sys;sys.path.insert(0,'src')
from petch.sheath import PeriodicSheathVoltage
from petch.reactor_boundary import PlasmaDiagnosticState, build_diagnostic_virtual_sheath_boundary
import numpy as np
d=PlasmaDiagnosticState(electron_density_m3=3e16, electron_temperature_eV=3.6, ion_name='ion',
    ion_mass_amu=40.0, source='probe', ion_flux_m2_s=1e20, ion_flux_evidence_kind='assumed')
wf=PeriodicSheathVoltage.sinusoidal(dc_v=2500.0, amplitude_v=2300.0, frequency_hz=1e6, source='probe')
for nt in (3,5,9,17):
    b=build_diagnostic_virtual_sheath_boundary(d,wf,reference_plane_m=0.0,
        collisionless_justification='probe', n_transverse_ion=nt)
    v=b.get('ion').velocity_sqrt_eV; w=b.get('ion').weight
    thx=np.degrees(np.arctan2(v[:,0],v[:,2])); th=np.degrees(np.arctan2(np.hypot(v[:,0],v[:,1]),v[:,2]))
    print(nt, np.sqrt(np.average(thx**2,weights=w)), np.abs(th).max(), w[np.abs(thx)<0.2865].sum()/w.sum())
"
```

| n_transverse_ion | sheath s | ⟨E⟩ | **σ_θ(1-D)** | **max\|θ\| represented** | w(\|θ_x\| < 0.286°) |
|---|---|---|---|---|---|
| 3 (default) | 8.73 mm | 2958 eV | **0.1478°** | 0.771° | 0.938 |
| 5 | 8.73 mm | 2958 eV | 0.1478° | 1.272° | 0.934 |
| 9 | 8.73 mm | 2958 eV | 0.1478° | 2.009° | 0.942 |
| 17 | 8.73 mm | 2958 eV | 0.1478° | 3.065° | 0.947 |
| **Krüger Fig-4 digitized (repo)** | — | 3465 eV | **0.833°** (repo band 0.822–0.860) | 2.86° | **0.249** |

Three pathologies, in order of severity:

1. **Missing tail (5.6×).** The IADF is thermal-only. Everything the measurement community has
   found about HAR IADFs since 2021 lives in the component petch does not have.
2. **Quadrature-order-dependent angular support.** σ_θ is *exactly invariant* under Gauss–Hermite
   order (GH preserves the second moment) while `max|θ|` runs 0.77° → 3.07°. For AR ≥ 100 the
   sidewall-flux *profile* is determined by exactly the part that is quadrature-determined, not
   physics-determined. **A Gauss–Hermite tensor product over transverse velocity is the wrong
   discretisation for HAR** — it is designed to converge moments, and HAR needs tails.
3. **9 discrete directions.** At n=3 the entire 2-D IADF is 3×3 = 9 transverse nodes. The single
   off-axis node sits at 0.36°, within 25% of the AR-200 acceptance angle — so the bottom/wall split
   at AR 200 is a step function of where one quadrature node lands.

**Acceptance-fraction table** (analytic, σ in 1-D projection at ⟨E⟩ = 3465 eV):

| IADF model | σ_core | σ_tail | AR 30 | AR 50 | AR 100 | **AR 200** |
|---|---|---|---|---|---|---|
| petch virtual sheath, T⊥ = 0.026 eV, thermal only | 0.111° | — | 1.000 | 1.000 | 1.000 | **0.990** |
| thermal only, T⊥ = 0.044 eV (measured core) | 0.144° | — | 1.000 | 1.000 | 1.000 | **0.953** |
| measured bi-Gaussian, tail frac 0.50 | 0.144° | 0.520° | 1.000 | 0.986 | 0.865 | **0.686** |
| measured bi-Gaussian, tail frac 0.65 (= Krüger s/λ) | 0.144° | 0.520° | 1.000 | 0.982 | 0.824 | **0.606** |
| Krüger Fig-4 digitized bi-G fit | 0.600° | 0.946° | 0.971 | 0.834 | 0.527 | **0.283** |

Sidewall flux at AR 200 = 1 − column. petch: 0.010. Measured bi-Gaussian: 0.31–0.39. Krüger HPEM:
0.72. **Ratio 30–70×.** Note this is the *virtual-sheath* path (the Tier-1 target path), **not** the
path used for the existing Krüger feature validation, which consumes the digitized IEAD directly.

**The 0.25° wall.** Independent of the sheath model, three places in petch quantise angle at 0.25°:
the digitized Krüger IEAD grid (verified: unique-angle clusters at exactly −2.75, −2.50, … +2.75);
the production transport quadrature (`joint ion-IEAD bins of 250 eV by 0.25 degrees`, per
`KRUEGER_2024_VALIDATION_PROTOCOL_2026-07-16.md`); and the effective angular support of the 3-node
GH set. AR 200 needs 0.057°. **The angular chain must be refined ~5× end-to-end, and the boundary
IEAD must come from something with better than 0.25° resolution — the digitized Krüger figure
physically cannot supply it** (25.0% of its weight lies in the single central bin spanning ±0.125°,
and the fraction inside |θ| < 0.143° and inside |θ| < 0.286° are *identical*, which is the
signature of an unresolved core).

## A.8 Hot/fast neutrals are a first-class angular species

Every measurement above reports **neutrals as well as ions**, and every one finds the **neutral
angular width larger than the ion width** (Ichikawa 2021, explicitly at all voltage conditions;
Kawamura 2025, explicitly). Physically: a charge-exchange event both creates the fast neutral and
(via momentum conservation with elastic scattering) sets its angle; downstream neutrals are then
scattered further with no field to re-collimate them.

Krüger's models treat this as standard: HPEM's PCMCM records *"the Ion/Neutral energy angular
distribution (IEAD and NEAD respectively)"* and the MCFPM launches *"all incoming ions or hot
neutrals resulting from neutralization of ions when striking surfaces"*
(`tmp/pdfs/krueger_thesis.txt` L2353–2356, L3081). Older literature: at 10–100 mTorr in CF₄,
*"the hot atom and ion fluxes to the substrate are comparable"* **[VERIFY primary source]**.

**petch's `PlasmaBoundaryState` has no fast-neutral species and the Tier-1 design doc's contract
table (§0) does not list one.** At AR 200 with a ~50% sheath-collision probability, roughly half of
the energetic flux is being modelled with the wrong charge state (charging!) and the wrong angular
law. This is the single largest *missing species* in the boundary contract.

---

# B. Tier-1 completeness — gap list against `RESEARCH_REACTOR_TIER1_DESIGN_2026-07-24.md`

The design doc's architecture (0-D global model → sheath-voltage waveform → reuse petch's sheath →
`PlasmaBoundaryState`) is **the right architecture and is confirmed by the literature**: the modern
precedent is exactly "volume-averaged global model + fluid/analytic sheath + Monte-Carlo collision
model in the sheath to get IEDF *and* IADF" (see §B.2). What follows are the load-bearing omissions.

## B.1 Gaps that block AR-200 capability (must fix)

**G1 — Gas temperature is deferred but is the IADF core width.**
Doc §6 "Explicitly deferred: … T_g energy balance"; doc §3(A) "v1 declares T_g from the recipe
(wall/wafer T + a fixed rise)". But θ_core = √(kT_gas/E). A declared T_g is a **declared IADF
width**, and §A.2 shows 300 K ↔ 1000 K is a 16× swing in AR-200 sidewall flux. *Fix:* promote T_g
from "used only to set the neutral thermal law" to a **primary output with a sensitivity band**, and
make `ion_tangential_temperature_eV` a **derived** quantity `= k·T_gas` (plus an evaluable
Wannier/D_T-μ correction), never a default constant. Today it is a keyword default of 0.026 with no
source — a knob, in violation of the standing doctrine.

**G2 — The collisionless sheath is assumed, but Krüger's own base case is 30–65% collisional.**
Doc §1.5 concedes "At 10–40 mTorr HARC pressures the sheath is *mildly collisional* — the
collisionless arcsine is a good leading model but the low-energy CX tail should be a declared
correction," and §6 defers "collisional CX low-energy IED tail". That framing is **energy-only**.
The literature says the important consequence is **angular**, not energetic, and it is not mild
(§A.3). `build_diagnostic_virtual_sheath_boundary` even takes a mandatory
`collisionless_justification` string — at HAR conditions there is no true justification to write.
*Fix:* the sheath module must carry an ion–neutral collision operator with an **anisotropic
differential cross section**, and must emit both the ion tail and the fast-neutral species.

**G3 — No fast-neutral (NEAD) species in the boundary contract.** §A.8. Blocks any honest AR-200
charging claim as well, since the neutralised fraction deposits no charge.

**G4 — The IADF representation cannot express a tail at any cost.** The 3-node Gauss–Hermite tensor
product is moment-converging by construction. Doc §0 says "the Tier-1 model does **not** invent a
new IEAD sampler … the *already-audited* `CollisionlessWaveformSheath` produces the finite-transit
IEAD." That reuse decision is the reason the gap exists. *Fix:* the reuse boundary should be drawn
at the **energy** law (finite-transit Child integration — genuinely good, keep) and **not** at the
angular law (replace).

**G5 — Gate 1a cannot validate what it claims to validate.** Doc §4 Gate 1a scores "mean energy
within 15%, horn split within 25%, angular width" against the digitized Fig-4. But the Fig-4
digitization is on a **0.25° grid** = the AR-200 acceptance half-angle (§A.7), so "angular width"
can only be validated to ±0.25° — the entire quantity of interest at AR 200. Worse, the reference is
**HPEM output, not measurement**, and Krüger's thesis says so explicitly: *"the fluxes and EADs from
the HPEM are treated as the ground truth"* (L4283) and *"These trends also rely on the fluxes
produced by the HPEM being truth and accurately representing experimental fluxes, which adds
additional variability"* (L4852). **There is no reactor-scale experimental validation anywhere in
Krüger's thesis.** *Fix:* add a **Gate 1a′** against the Nagoya measured bi-Gaussian (core/tail
temperatures, width-vs-energy law, tail-fraction-vs-pressure exponential) — real measurements, at
keV, at better than 0.1° resolution, which is the only dataset in the world that can grade a
sub-degree IADF model.

**G6 — Mean-angle tilt / radial non-uniformity is out of scope but not declared.** §A.6: measured
tilt ≤ 2.5° at the edge, "a 2-degree tilt can be a defect", vs 0.286° acceptance. A 0-D Tier-1
boundary is a **center-of-wafer** object. *Fix:* make this an explicit provenance field
(`radial_position_m`, `mean_incidence_angle_deg = 0 (declared, 0-D)`), refusing edge-die claims.

## B.2 Gaps in the literature survey (add these)

- **`CollisionlessWaveformSheath` has a named academic sibling worth benchmarking against:**
  Miller & Riley, *"Dynamics of collisionless rf plasma sheaths"*, **J. Appl. Phys. 82, 3689–3709
  (1997)**, DOI [10.1063/1.365732](https://doi.org/10.1063/1.365732) — semianalytic RF sheath from an
  approximate first integral of Poisson's equation, embeddable in an external-circuit model,
  predicting IEDs; **and it was tested against an experiment**, which petch's version has not been.
  Missing entirely from the design doc.
- **The sheath↔circuit closure**: Edelberg & Aydil, **J. Appl. Phys. 86, 4799–4812 (1999)**,
  DOI [10.1063/1.371446](https://doi.org/10.1063/1.371446) — self-consistent dynamic sheath coupled to
  an equivalent-circuit model, compared to measurement. The design doc mentions Edelberg & Aydil only
  as "measured IEDFs", missing that this is the canonical **module-B** (recipe → V_sh(t)) model.
- **The exact Tier-1 architecture, already published and validated:** the "global model + fluid ion
  sheath + Monte-Carlo collision model to obtain IEDF *and* IADF" pattern for biased ICPs
  (Plasma Process. Polym. 2016, DOI [10.1002/ppap.201600100](https://doi.org/10.1002/ppap.201600100))
  **[VERIFY authors/volume — Wiley 403'd]**; and the multi-dimensional-sheath + Monte-Carlo IEADF
  route, Denpoh, **Jpn. J. Appl. Phys. 53, 080304 (2014)**,
  DOI [10.7567/jjap.53.080304](https://doi.org/10.7567/jjap.53.080304). Both include the MC
  collision step the design doc omits.
- **A modern hybrid comparator on our exact chemistry class:** Ducluzaux, Ristoiu, Cunge,
  Despiau-Pujo, *"Impact of plasma operating conditions on the ion energy and angular distributions
  in dual-frequency capacitively coupled plasma reactors using CF₄ chemistry"*,
  **J. Vac. Sci. Technol. A 42, 013002 (2024)**,
  DOI [10.1116/6.0003291](https://doi.org/10.1116/6.0003291).
- **Kaganovich/PPPL is actively working the same tail** (arXiv:2604.04214, Apr 2026, funded by
  Samsung + DOE/PPPL CRADA with Applied Materials). This is a *competitive* signal: the wide-angle
  IADF tail is currently a named open problem at a national lab with industrial funding, and their
  own model still under-predicts the measured tail. Petch has a differentiable-transport advantage
  here that PPPL's Monte-Carlo does not.

## B.3 Citation-ledger corrections to the design doc

| Design-doc entry | Correction |
|---|---|
| "Analytical model for ion angular distribution functions…" DOI 10.1063/1.1519941 [VERIFY] | **Wrong DOI.** Correct: **10.1063/1.1524020**, Raja & Linne, *J. Appl. Phys.* **92**, 7032–7040 (2002). |
| Kawamura, Vahedi, Lieberman, Birdsall, PSST **8**, 313 (1999) | **Wrong pages.** Correct: PSST **8**, **R45–R64** (1999); title is *"Ion energy distributions in rf sheaths; review, analysis and simulation"*. DOI 10.1088/0963-0252/8/3/202 is right. |
| Krüger *et al.*, Phys. Plasmas 31(3), 033508 (2024), DOI 10.1063/5.0189675 [VERIFY] | **DOI does not resolve.** Correct: **10.1063/5.0189397**; authors Krüger, Lee, Nam, Kushner. Volume/article 31, 033508 confirmed. |
| Krüger, Zhang, Luan, Park, Metz, Kushner, JVST A 42(4), 043008 (2024), DOI 10.1116/6.0003554 | **Confirmed** via Crossref (title and authors exact). |
| Gate-1 recipe "10 mTorr C₄F₆/Ar/O₂, 140/100/105 sccm, 1/40 MHz CCP" | **Incomplete.** Thesis §6.3: electrodes **r = 15 cm**, gap **4 cm**, **P_lf = 8.0 kW**, **P_hf = 2.5 kW**, and **a −500 V DC bias applied to the *top* electrode delivering 650 W**. Reference plasma state: **n_e,max = 9.5×10¹⁰ cm⁻³**, **T_e = 3.4–3.8 eV over the wafer**, **[F⁻]_max = 2.2×10¹⁰ cm⁻³**, fluxes tabulated **at r = 7.5 cm**, ion energies to **4800 eV**. Using n_e,max (bulk centre) in the Child-law thickness would be a systematic error — the sheath-edge density is the correct input. |

## B.4 What the design doc gets right (keep, unchanged)

- The `TabulatedReactorFluxDeck` + `PlasmaDiagnosticState` + waveform adapter decomposition.
- Finite-transit Child-profile integration for the **energy** law (matches Krüger's own S-parameter
  framing; genuinely first-principles).
- Factor-of-2 honesty on absolute global-model fluxes, and the calibrated-transfer claim ceiling.
- Gate 1b's routing around the missing C₄F₆ gas-phase set.
- Preregistered bands + `claim_mode` gating.
- The dual-frequency "flux from source power, energy from bias waveform" split — Krüger's own
  reactor uses exactly this (500 W @ 80 MHz top / VWT 1 MHz bottom in Ch. 4; 2.5 kW @ 40 MHz /
  8.0 kW @ 1 MHz in Ch. 6).

---

# C. Combined staged plan — how sub-degree IADF + Tier-1 + azimuth transport compose into AR-200

**Design principle.** All three deficits are the *same* deficit seen from three sides: petch resolves
angle at 0.25° and represents it with moment-converging quadratures, while AR 200 is a
**tail-and-0.05°** problem. Fix the representation once, and the three fixes compose; fix them
separately and they will keep multiplying (today: 5.6× IADF width × 1.5× azimuth ⇒ ~30–70× sidewall
flux at AR 200).

### S0 — Build the ruler first (½–1 day, no physics)
An angular-convergence harness, before anything else. For a fixed straight feature at AR = 30 / 100 /
200, report sidewall flux vs depth as a function of (a) IADF polar-bin width, (b) ion-azimuth node
count, (c) Gauss–Hermite order, (d) compressed-azimuth depth cells. Outputs: the *measured*
AR-dependence of the 1.5× azimuth deficit, and the bin width at which each observable converges to
1%.
**Gate S0:** publish the convergence tables. Nothing downstream is believable without them.
*Why first:* it is the only way to distinguish "the tail is missing" from "the transport loses the
tail", and it directly tests the standing hypothesis that the ml13 mouth-opening gap
(24.8 vs 45 nm, `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md`) is a wide-angle-tail transport loss —
wide-angle ions land near the mouth (a 5° ion strikes the wall at depth ≈ 11×W).

### S1 — Replace the angular law of the boundary (2–3 days)
New `IonAngularEnergyDistribution` object replacing the transverse Gauss–Hermite tensor in
`collisionless_sheath_boundary_state` (`src/petch/boundary_state.py:853`). Requirements:
- explicit (E, θ, φ) representation with **arbitrary polar resolution** and support out to ≥ 5°;
- **two-component** structure (thermal core + non-thermal tail) so each can be integrated against the
  geometric acceptance **analytically** — an erf per component, exactly as in the §A.7 table — giving
  exact bottom/wall splits at any AR with no quadrature noise;
- keeps the existing phase↔energy pairing (§A.4) so the energy–angle correlation survives;
- `T_perp` **derived** from `T_gas` (G1), never a keyword default;
- azimuthal dependence allowed (Raja & Linne drift term) so S3 has something to consume.
**Gate S1 (measured, not HPEM):** reproduce Kim 2025 — core 0.044 eV, tail 0.57 eV, main width
∝ E^(−1/2) monotone over 1.4–2.0 keV; and reproduce the repo's Krüger Fig-4 marginal σ within its
own digitization band 0.822–0.860° when driven with the fitted bi-Gaussian.

### S2 — Derive the tail instead of fitting it (1–2 weeks) — *the frontier item*
Sheath ion–neutral collision operator producing the tail **and** the fast-neutral species:
- differential cross sections: Born–Mayer repulsive potential per Khrabrov & Kaganovich
  (arXiv:2604.04214) — **do not** use Phelps-1994 alone (no ion elastic tail) or
  Phelps–Greene–Burke-2000 alone (~2× too much 1–5° scattering);
- survival law calibrated to the measured exponential main/total-vs-pressure trend
  (Kim 2025, JJAP 64, 096002);
- emit `SpeciesBoundaryState(role="fast_neutral", charge_number=0, ...)` (G3) — this also unblocks
  honest AR-200 charging, since neutralised flux deposits no charge.
**Gate S2:** predict tail fraction vs pressure without fitting it; independently, reproduce the
consistency check found in this pass — the digitized Krüger tail fraction (0.65) equals
1 − exp(−s/λ) at s = 14 mm, λ = 13.8 mm (0.64).
*This is the deterministic/differentiable-transport moonshot in miniature: PPPL is doing this with
Monte Carlo and still under-predicts the measured tail.*

### S3 — Azimuth transport fix + angular AMR (1 week, concurrent with S2)
- Fix the compressed-azimuth wall-flux under-delivery using S0's AR-resolved measurement of it.
- **Angular AMR keyed to the live aspect ratio:** polar bin ≤ arctan(1/AR)/5 (0.057° at AR 200),
  refined only inside the acceptance cone; the tail (θ > a few × acceptance) needs *no* refinement,
  because those ions all terminate near the mouth. This makes the 5× refinement nearly free.
- Retire the fixed `250 eV × 0.25°` production bin in favour of the AMR schedule.
**Gate S3:** S0's convergence harness returns < 1% change at AR 200 under 2× refinement in both
polar and azimuth; the azimuth deficit collapses to < 5%.

### S4 — Tier-1 module B (sheath-voltage waveform) — *unchanged from the design doc, revised inputs*
`reactor_tier1/sheath_voltage.py` as designed (CCP voltage division → V_dc + harmonic
`PeriodicSheathVoltage`), but seeded with the **complete** Krüger recipe from §B.3 (P_lf = 8.0 kW,
P_hf = 2.5 kW, −500 V on the top electrode, r = 15 cm, gap 4 cm) and with a **sheath-edge** density
rather than n_e,max in the Child-law closure. Benchmark the sheath dynamics against Miller & Riley
(1997) and Edelberg & Aydil (1999) rather than only against Krüger.
**Gate S4 = the revised Gate 1a:** score against the digitized Fig-4 **for energy only** (mean energy
15%, horn split 25%) and against **Kim 2025 for angle** (§S1 gate). Explicitly declare that Fig-4's
0.25° grid cannot grade angle.

### S5 — Tier-1 module A (global chemistry), with T_gas promoted (1–2 weeks)
As designed (`global_model.py`, CF₄/O₂/Ar from Gudmundsson 2019 first), **plus a gas-temperature
energy balance** (G1) reported with a sensitivity band, because T_gas *is* the IADF core width.
Gate 2-lite unchanged. Add the `radial_position_m` / `mean_incidence_angle_deg` provenance fields
(G6) so edge-die claims are structurally refused.

### S6 — AR-200 capability demonstration (the payoff)
Compose: Tier-1 recipe → (S5) fluxes/T_e/T_gas → (S4) V_sh(t) → (S2) collisional sheath →
(S1) bi-component IADF + NEAD → (S3) AMR angular transport at 0.057° → feature engine at AR 200.
**Gate S6:** bottom-vs-sidewall flux split at AR 200 within the band set by the measured
bi-Gaussian (§A.7 table: 0.61–0.69 bottom fraction), converged under S0's harness, with the tail
fraction *predicted* by S2 rather than fitted.

### Ordering rationale (one line each)
- **S0 before everything** — you cannot debug a 30–70× discrepancy without knowing which factor is
  representation and which is physics.
- **S1 before S2** — you need a container that can hold a tail before you can compute one.
- **S3 concurrent with S2** — independent code paths, and S3's AMR is what makes S2's tail affordable.
- **S4/S5 after S1** — Tier-1's whole value is producing an IADF; producing it into a container that
  cannot express a tail would bake the current error into the reactor layer.
- **Do not start Tier-1 module A (chemistry) first** even though it is the largest chunk of the
  design doc — it is the piece with the *least* leverage on AR 200.

---

## Appendix — citation ledger (all DOIs resolved via api.crossref.org this session unless flagged)

**Measured sub-degree IADFs (the anchor set)**
- K. Ichikawa, M. H. Chu, M. Moriyama, N. Nakahara, H. Suzuki, D. Iino, H. Fukumizu, K. Kurihara [+ H. Toyoda] — *Angular distribution measurement of high-energy argon neutral and ion in a 13.56 MHz capacitively-coupled plasma* — **Appl. Phys. Express 14, 126001 (2021)** — DOI 10.35848/1882-0786/ac33c4
- D. Kim, S. Kawamura, M. Naito, D. Iino, H. Fukumizu, K. Kurihara, H. Suzuki, H. Toyoda — *Measurement of energy-resolved ion angular distribution in a dual-frequency capacitively coupled argon plasma* — **Jpn. J. Appl. Phys. 64, 05SP15 (2025)** — DOI 10.35848/1347-4065/adce84
- D. Kim, S. Kawamura, K. Fujitani, M. Naito, D. Iino, H. Fukumizu, K. Kurihara, H. Suzuki [+ H. Toyoda] — *Influence of sheath collisions on the ion angular distributions in a dual-frequency capacitively coupled plasma* — **Jpn. J. Appl. Phys. 64, 096002 (2025)** — DOI 10.35848/1347-4065/ae0105
- S. Kawamura, D. Kim, K. Ichikawa, H. Suzuki, D. Iino, H. Fukumizu, K. Kurihara, H. Toyoda — *Angular distribution measurement of high-energy ions and neutrals impinging on an RF electrode in a dual-frequency capacitively-coupled Ar plasma* — **Plasma Sources Sci. Technol. 34, 055006 (2025)** — DOI 10.1088/1361-6595/add321
- J. R. Woodworth, M. E. Riley, D. C. Meister, B. P. Aragon, M. S. Le, H. H. Sawin — **J. Appl. Phys. 80, 1304–1311 (1996)** — DOI 10.1063/1.362977
- J. R. Woodworth, M. E. Riley, P. A. Miller, G. A. Hebner, T. W. Hamilton — **J. Appl. Phys. 81, 5950–5959 (1997)** — DOI 10.1063/1.364383
- S. Sharma, D. Gahan, P. Scullin, S. Daniels, M. B. Hopkins — *Ion angle distribution measurement with a planar retarding field analyzer* — **Rev. Sci. Instrum. 86, 113501 (2015)** — DOI 10.1063/1.4934808

**IADF/sheath theory**
- A. V. Khrabrov, I. D. Kaganovich (PPPL) — *Ion-neutral and neutral-neutral scattering in argon at KeV energies and implications for high-aspect-ratio etching* — **arXiv:2604.04214v2** (5 Apr 2026, rev. 8 Apr 2026) — full text fetched and quoted verbatim above
- L. L. Raja, M. Linne — **J. Appl. Phys. 92, 7032–7040 (2002)** — DOI 10.1063/1.1524020
- P. A. Miller, M. E. Riley — *Dynamics of collisionless rf plasma sheaths* — **J. Appl. Phys. 82, 3689–3709 (1997)** — DOI 10.1063/1.365732
- E. A. Edelberg, E. S. Aydil — **J. Appl. Phys. 86, 4799–4812 (1999)** — DOI 10.1063/1.371446
- E. Kawamura, V. Vahedi, M. A. Lieberman, C. K. Birdsall — *Ion energy distributions in rf sheaths; review, analysis and simulation* — **Plasma Sources Sci. Technol. 8, R45–R64 (1999)** — DOI 10.1088/0963-0252/8/3/202
- K. Denpoh — **Jpn. J. Appl. Phys. 53, 080304 (2014)** — DOI 10.7567/jjap.53.080304
- *(global model + fluid sheath + MCC → IEDF/IADF in biased ICP)* — **Plasma Process. Polym. (2016)** — DOI 10.1002/ppap.201600100 — **[VERIFY authors/volume; Wiley 403]**

**Cross sections / transport coefficients**
- A. V. Phelps — *The application of scattering cross sections to ion flux models in discharge sheaths* — **J. Appl. Phys. 76, 747–753 (1994)** — DOI 10.1063/1.357820
- A. V. Phelps, C. H. Greene, J. P. Burke Jr. — **J. Phys. B 33, 2965–2981 (2000)** — DOI 10.1088/0953-4075/33/16/303
- T. Stefánsson, H. R. Skullerud — *Measurements of the ratio between the transverse diffusion coefficient and the mobility for argon ions in argon* — **J. Phys. B 32, 1057–1066 (1999)** — DOI 10.1088/0953-4075/32/5/001

**Reactor / feature context**
- F. Krüger, D. Zhang, P. Luan, M. Park, A. Metz, M. J. Kushner — *Autonomous hybrid optimization of a SiO₂ plasma etching mechanism* — **J. Vac. Sci. Technol. A 42, 043008 (2024)** — DOI 10.1116/6.0003554
- F. Krüger, H. Lee, S. K. Nam, M. J. Kushner — *Voltage waveform tailoring for high aspect ratio plasma etching of SiO₂ using Ar/CF₄/O₂ mixtures: Consequences of low fundamental frequency biases* — **Phys. Plasmas 31, 033508 (2024)** — DOI **10.1063/5.0189397** (not 5.0189675)
- F. Krüger — PhD thesis (Univ. Michigan), local copy `tmp/pdfs/krueger_thesis.txt`; reactor details §4.2, §6.3; HPEM-as-ground-truth admissions L4283, L4852; sheath thicknesses L2798–2800; hot-neutral/NEAD L2353–2356
- S. Huang, C. Huard, S. Shim, S. K. Nam, I.-C. Song, S. Lu, M. J. Kushner — *Plasma etching of high aspect ratio features in SiO₂ using Ar/C₄F₈/O₂ mixtures* — **J. Vac. Sci. Technol. A 37 (2019)** — DOI 10.1116/1.5090606
- P. Ducluzaux, D. Ristoiu, G. Cunge, E. Despiau-Pujo — **J. Vac. Sci. Technol. A 42, 013002 (2024)** — DOI 10.1116/6.0003291
- I. Seong, J. Lee, S. Kim, Y. Lee, C. Cho, J. Lee, W. Jeong, Y. You *et al.* — *Characterization of an Etch Profile at a Wafer Edge in Capacitively Coupled Plasma* — **Nanomaterials 12, 3963 (2022)** — DOI 10.3390/nano12223963
- R. A. Gottscho, C. W. Jurgensen, D. J. Vitkavage — *Microscopic uniformity in plasma etching* — **J. Vac. Sci. Technol. B 10, 2133–2147 (1992)** — DOI 10.1116/1.586180
