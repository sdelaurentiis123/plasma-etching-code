# Deterministic reactor-tier advance — 2026-08-09

## Bottom line

The reactor stack now reaches beyond 0-D chemistry. It has deterministic,
conservation-audited implementations for axisymmetric neutral/reaction
transport, quasineutral positive/negative-ion drift diffusion, harmonic ICP
electromagnetics, line-resolved resonance radiation, finite-wafer delivery,
RF sheath transfer, and evolving-feature coupling. The spatial charged and
electromagnetic blocks expose exact implicit or analytic JVPs and reuse their
sparse factorizations; no particle or photon Monte Carlo is used.

That advance does **not** make the current depth boards green. It produces a
stronger result: on the 11-condition Mahorowala chlorine board, physically
distinct radial source fields—including a solved complex-conductivity ICP
field—change evolving-feature MAPE by at most `0.68` percentage points. The
errors remain mixed sign (`6` over, `5` under). Missing radial plasma shape is
therefore not the dominant Mahorowala absolute-depth failure.

Krüger remains `346.833 nm` simulated versus `825 nm` measured under the
published aggregate boundary. The new ICP/chlorine tiers cannot be transplanted
to his dual-frequency Ar/C4F6/O2 CCP, and no solver tier can reconstruct his
unpublished species-resolved positive-ion composition or stable C4F6 wafer
flux from feature depth without turning the result into a fit.

## Landed rungs and their actual claim

| rung | implementation | present claim boundary |
|---|---|---|
| axisymmetric finite-volume reaction diffusion | `src/petch/reactor_global/axisymmetric_reaction_diffusion.py` | exact local particle ledgers and source JVP; spatial lift only when chemistry inventory is supplied |
| quasineutral charged drift diffusion | `src/petch/reactor_global/axisymmetric_drift_diffusion.py` | Scharfetter--Gummel ion flux, Boltzmann electrons, Bohm positive-ion wall, absorbing negative-ion wall, exact implicit inventory JVP; mutual neutralization is still volume-average linearized |
| reusable chlorine wafer provider | `src/petch/reactor_global/chlorine_axisymmetric_transport.py` | species-resolved finite-wafer Cl+/Cl2+ flux conditional on a solved 0-D state and declared source moment |
| harmonic ICP field | `src/petch/reactor_global/axisymmetric_inductive_field.py` | complex Drude conductivity, toroidal absorbed-power shape, exact complex-power ledger and conductivity/coil-field JVP; generator setpoint does not yet determine coil-side field amplitude |
| axisymmetric/line radiation | `axisymmetric_radiation.py`, `vuv_radiation.py`, `zonal_radiation.py`, `chlorine_resonance_boundary.py` | deterministic finite-wafer line transfer, optically thick resonance escape and partial redistribution; absolute upper-state populations and several quenching rates remain evidence-limited |
| sheath and feature coupling | `wafer_sheath_transfer.py` plus the released feature engine | reactor species flux and RF waveform can drive the same evolving geometry; no depth selects a boundary parameter |

## Independent physics checks

### Wise GEC chlorine spatial board

The original NIST hardware drawing was visually audited at 300 dpi. The
landed geometry receipt carries the `165.1 mm` lower electrode, `40.5 mm`
gap, `15 mm` probe plane, five-turn coil, and source hashes. Wise Figure 3 was
PIL-digitized into 21 direct markers—seven each for electron density,
temperature, and plasma potential.

Those three measured fields directly check the bulk closure
`grad(phi)=grad(n_e T_e)/n_e`. Aligning only the arbitrary potential gauge at
the axis, the measured-density/temperature integral reconstructs the measured
potential with `5.07%` unweighted MAPE. No slope or temperature is fitted.
The digitized density FWHM is `7.53 cm`, inside Miller's independent typical
`7--9 cm` interval. Source error bars are absent, so this supports the closure
but is not relabeled as a formal uncertainty-weighted pass.

Receipt:
`results/curated/reactor_to_feature_chlorine/wise_1996_boltzmann_closure_audit/`.

### Electromagnetic numerics

The azimuthal-field implementation passes a manufactured Bessel/skin mode
with second-order refinement, closes its complex power ledger to roundoff,
and matches centered differences with its exact parameter JVP. Across the
Mahorowala conditions the reconstructed complex conductivity has positive
dissipative part, and the solved power layer has a roughly `10.4--12.9 mm`
mean top depth. The coil-side field is a declared geometry sensitivity, not a
power-setpoint calibration.

### Radiation numerics and source boundary

The dominant calculated Cl I shortwave source is at `118.88`, `109.74`, and
`110.75 nm`, not in the `104.82--106.67 nm` band where the landed large
photoetch yield was measured. The 118.88-nm partial-redistribution transport
changes by `0.160%` from its 12x12 production rule to the 16x16 check, after
exactly splitting the finite-wafer radial discontinuity.

## Exact Mahorowala reachability result

The same independent 400 W / 100 sccm center-current datum normalizes each
source family once. Every other condition and every feature depth is held out.

| source boundary | evolving-feature MAPE | sign count |
|---|---:|---:|
| global edge | `19.47%` | 6 over / 5 under |
| Wise-class top annulus | `19.28%` | 6 / 5 |
| compact top-center | `20.15%` | 6 / 5 |
| solved inductive EM annulus | `19.59%` | 6 / 5 |

Receipt:
`results/curated/reactor_to_feature_chlorine/axisymmetric_depth_reachability/`.

The mixed signs prove that neither one flux multiplier nor one additive
photon rate can close the board. Runs 6, 9, and 13 already overpredict before
an additive photochannel. Spatial transport remains necessary for radial and
equipment transfer, but it is not an authorized absolute-depth knob.

## What is and is not solved

- **Solved as implementation physics:** deterministic conservation, spatial
  drift/diffusion, quasineutral field closure, ICP absorbed-power shape,
  finite-wafer line radiation, sheath interface, and feature propagation.
- **Supported directly:** Wise electron-pressure/Boltzmann closure and density
  width; Tian spatial resonance mechanism reproduction; numerical ledgers and
  derivatives.
- **Still conditional:** generator power -> coil field/absorbed power;
  full local electron energy/EEDF; local nonlinear mutual neutralization;
  equipment wall state; species-resolved wafer IEAD and VUV amplitude.
- **Not matched:** Krüger absolute depth; Mahorowala held-out feature-depth
  board. Formal blind held-out feature-depth passes remain zero.
- **Other chemistry depth evidence remains intact:** SF5+ direct-beam
  depth-per-dose is `5.88%` MAPE over four points; Si/Cl2/Ar+ ALE has `12.88%`
  maximum relative error over three points. Those are surface-depth boards,
  not reactor-to-feature validations.

## Next solver rung—without using depth

The next useful reactor build is a self-consistent deterministic block solve,
not another source-shape family:

1. solve the complex ICP field from a measured coil-current/voltage phase or
   a validated matching-network/electromagnetic boundary;
2. solve spatial electron energy/EEPF and electron-impact rates on that power
   field, with exact implicit JVPs through the coupled residual;
3. replace volume-average mutual neutralization by local nonlinear
   Cl+/Cl2+/Cl- loss and couple neutral Cl/Cl2 transport and wall state;
4. close species-resolved presheath/sheath transfer and propagate covariance
   to the feature plane;
5. grade reactor state and wafer flux before opening feature depths.

This rung can improve machine transfer and possibly the condition dependence,
but the present Mahorowala result forbids claiming in advance that it will
remove the absolute-depth residual.

## Minimal experiments that now block prediction

### Lam/Mahorowala chlorine ICP

Run at least the center condition plus two sweep endpoints, synchronously:

1. calibrated coil voltage/current phase and a B-dot or equivalent field map,
   or a direct absorbed-power deposition diagnostic;
2. radial/axial `n_e`, EEDF/`T_e`, and plasma potential;
3. mass-resolved Cl+/Cl2+ flux and IEAD at the wafer plane;
4. Cl density/dissociation and documented wall-conditioned recombination;
5. absolute 109.7, 110.8, and 118.9 nm wafer flux;
6. chlorinated-poly-Si photon yield at those wavelengths with and without the
   target RF ion waveform, plus feature-floor photon delivery for the 310 nm
   opening.

Items 1--4 make the reactor boundary out of sample. Items 5--6 decide whether
the remaining mixed-sign depth residual is surface/photo synergy rather than
plasma delivery.

### Krüger Ar/C4F6/O2 CCP

At the exact published recipe, measure:

1. species-resolved positive-ion flux and IEAD/IAD at the wafer;
2. stable C4F6 and CFx/F neutral flux at the wafer;
3. blanket SiO2 and mask rates on the same conditioned chamber;
4. preferably one in-feature floor flux or reaction-state measurement.

Without at least items 1--3, the `825 nm` endpoint is not identifiable from
the publication. A 0-D or 2-D reactor solver can consume those measurements
or predict them after an equipment-specific grade; it cannot infer them from
the endpoint without circularity.

## Version-control boundary

The repository-wide verification after this tier advance is **1,856 passed,
1 skipped** in `12m47s`. The focused new reactor/surface/radiation set is
**112 passed**. The library contains **147 source cards** after adding the
primary Miller GEC-ICP geometry record.

All source cards, digitized markers, code, tests, and derived receipts are
intended for the active `codex/validation-first-multiphysics` branch. Raw
publisher PDFs, rights-limited OPEN-ADAS/LXCat payloads, and local user files
remain outside Git. The protected mixed-layer log and mouth-equilibrium probe
directory are not part of this campaign and must remain untouched.
