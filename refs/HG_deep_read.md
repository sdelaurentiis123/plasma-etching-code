# Deep read: HG JVST B 15,70 (1997) methods paper — full extraction (2026-07-07)

Read page-by-page from refs/HG_jvstb97.pdf (visual). This resolves most of the residual ledger.

## The algorithm (verbatim structure)
- Eq 3.1: 2D (x,y) equations of motion; no transverse field in the sheath; E_x significant near/in
  the microstructure. Electrons inertialess (respond to instantaneous field).
- Ion entry: Gaussian parallel velocity, mean = Bohm (kTe/M)^1/2, variance (kTi/M)^1/2, Ti=0.5 eV
  ("collisional broadening... presheath"); **"ions enter the sheath with a ratio of directed energy
  to transverse temperature of 8"** — exactly our Te/2 ÷ Ti/2 convention. Entry physics MATCHES ours.
- **IEDF (Fig 4a): ASYMMETRIC bathtub — "The high energy peak has lower intensity than the low
  energy peak, as expected from the self-consistent treatment of the sheath at the low rf bias
  frequency."** Low horn ≈ 2× high horn. Our instantaneous-crossing bathtub has equal horns →
  we UNDER-weight low-energy ions. (Our old ied_bias=0.25 heuristic mimicked their sheath better
  than the pure instantaneous derivation.) THIS is the bottom-flux (+30%) driver.
- IADF 4.3° HWHM (consistent w/ Woodworth measurements); energy-dependent (wing ions slowest) ✓ ours.
- EADF Fig 5b: dotted cos^0.6 fit over a NOISY MC histogram ✓ (injection: "isotropic flux
  distribution" at the sheath top = uniform-in-angle; our theorem stands).
- Charging step = 50 ions + 50 electrons ≈ half rf cycle → Laplace ∇²V=0 (FDM), BCs: V=0 at the
  sheath lower boundary, ∇V=0 at left/right centerlines (mirror) ✓ ours.
- **Four-step poly-Si charge redistribution (§III.D, verbatim):** (1) potential at each poly surface
  cell = Coulomb contributions from ALL charged cells (PR, SiO2, incl. through-dielectric ε=1.6
  effects) + mirror images; (2) arithmetic mean over both sidewalls = the equipotential; deviation
  per cell; (3) add surface charge ∝ deviation; (4) repeat to equipotential. Charge PILES at the
  poly/SiO2 foot ("critically important" for notching). NOTE: their published poly "equipotential"
  labels are COULOMB-MAP READOUTS of the redistributed charge.

## Steady state + the JVST B numbers (their AR=2 structure)
- Currents balance at ~1500 steps (0.965 — **3.5% of launched particles never hit**; ours 0%).
- POTENTIAL steady state needs ~7000 steps (3500 rf cycles).
- **Fig 7b floor profile: plateau ~20–22 V with foot peak 58.7 V** ~40 nm from the outermost line's
  inner foot. Fig 8 poly labels: outermost (edge) 7.8, neighbor 19.8.
- Kinoshita contrast: symmetric bottom potential is "not physical"; ions steered toward the LOWER
  line (the edge) — the notch side.

## The reconciliation ledger (their AR2 vs our phys-mode AR2)
| quantity | HG (JVST B / JAP) | petch | verdict |
|---|---|---|---|
| floor plateau | ~20–22 | 23.7 | ✓ |
| foot peak | 58.7 | ~61 | ✓ |
| edge line | 7.8 (JVST B) but ≈4 (JAP Fig 6!) | 4.4 | ✓ vs JAP — **HG's own two papers differ ~2× on this label** |
| neighbor | 19.8 / ≈19 | 11.0 | the one real magnitude gap |
| bottom flux | ~0.33 | 0.46 | IEDF asymmetry (above) |
| survivors | 3.5% | 0% | minor |

## What remains, precisely
1. **IEDF asymmetry** — implement the nonlinear-sheath weighting (low horn ≈ 2× high). Expected to
   close the bottom-flux offset and pull floorV/foot-E endpoints in.
2. **Neighbor magnitude** — likely tied to their charge-space conductor treatment (Coulomb-map
   labels + foot charge pile). Implement §III.D as a dynamics mode.
3. Their own inter-paper label spread (edge 4 vs 7.8) bounds the achievable "exactness" of any
   label comparison at ~2× — the flux/energy/structure channels are the real gates.
