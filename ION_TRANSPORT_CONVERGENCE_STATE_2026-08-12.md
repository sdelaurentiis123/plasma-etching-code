# Ion transport convergence state — 2026-08-12

## Outcome

The deterministic pure-Ar reactor-to-wafer stack no longer truncates the
ion-neutral collision-order expansion.  It now solves the complete discrete
linear Boltzmann series as an absorbing sparse system on fixed
potential/energy/angle ordinates.  At the frozen 1 keV, 10 mm, 500 K audit
condition, both 1 mTorr and 10 mTorr pass the dual-grid gate.

This closes a numerical failure, not the Krüger fluorocarbon depth gap.  The
effective-static operator described here also is not the final sheath model.
The separate current-driven moving-sheath path landed on 2026-08-13 and is
documented in `MOVING_RF_SHEATH_REACTOR_BREAKTHROUGH_2026-08-13.md`.
Krüger Figure 4 is already a downstream HPEM wafer-plane combined-positive-ion
IEAD.  Applying this Ar sheath after Figure 4 would double-count upstream
transport, and Ar+/Ar cross sections do not identify the paper's unknown CFx+
mixture.  The certified Krüger result therefore remains 346.83264081620524 nm
from the published aggregate boundary versus 825 nm measured.

## What changed

- Added a bounded deterministic discrete-ordinates operator.  For each RF
  phase it solves `(I - Q.T) x = s`, summing all ion collision orders without
  particle Monte Carlo or an exponentially growing path tree.
- Added the analytic neutral-density JVP of the implicit solve:
  `(I - Q.T) dx = ds + dQ.T @ x`.
- Retained equal-mass elastic/CX kinematics and Phelps/Born-Mayer provenance.
- Followed backscattered ions through electrostatic deceleration and turning;
  only ions energetic enough to reach the plasma edge escape.
- Made the implicit solver the default pure-Ar reactor-to-wafer collision
  transfer.  The finite collision-order implementation remains an explicit
  reference mode.
- Added a dual-grid certification covering ion probability, mean energy, RMS
  angle, expected collision count, and the resolved fast-neutral lower-bound.

## Frozen high-pressure receipt

At 10 mTorr, the old three-collision reference left 0.0530040 ion probability
unresolved and its order-2 to order-3 mean-energy shift was 3.5855%.

The new implicit operator gives:

| Quantity | Coarse grid (7/7/10) | Fine grid (9/9/13) | Relative change |
|---|---:|---:|---:|
| unresolved ion probability | 0 | 5.55e-16 | — |
| mean ion energy | 758.4427 eV | 763.2048 eV | 0.62396% |
| RMS polar angle | 0.658295° | 0.663897° | 0.84391% |
| expected ion collisions | 1.421419 | 1.413540 | 0.55433% |
| resolved fast-neutral arrivals/source ion | 0.360589 | 0.360630 | 0.01138% |

The maximum sparse linear-system residual is 8.33e-17 and the maximum
per-state probability-ledger residual is 1.33e-15.  The frozen machine-readable
receipt is `results/curated/collisional_sheath_depth_scope/audit.json` schema
`petch.collisional-sheath-depth-scope.v2`.

## Remaining physical, not numerical, limits

1. Subsequent fast-neutral Ar--Ar collisions remain an explicit unresolved
   ledger because they require neutral--neutral differential scattering and
   loss channels, not the Ar+--Ar law used for the ion solve.  The exported
   neutral boundary is a lower bound and cannot yet certify depth.  The later
   2026-08-13 current-driven path separately closes low-energy *ion* angles
   with the Phelps/LXCat isotropic/backscatter decomposition.
2. The legacy default in this receipt remains an effective phase-conditioned
   Child profile.  The new current-driven path resolves the moving electron
   front and RF phase as a kinetic coordinate, conditional on a sheath-current
   waveform; it does not retroactively change this frozen receipt.
3. The Ar global stack consumes absorbed bulk and delivered bias power, not
   generator forward power; chamber/circuit coupling remains equipment data.
4. Krüger still lacks species-resolved positive-ion flux/IEAD, stable C4F6
   wafer flux, molecular-ion collision data, and a validated reactor waveform.

No feature depth was used to select any coefficient or grid.  No Krüger input,
yield, or profile was altered by this work.

The numbers in this document are a frozen 2026-08-12 receipt for the former
hybrid collision law.  They are retained as history; current predictions use
the 2026-08-13 Phelps low-energy decomposition documented in the moving-sheath
receipt.
