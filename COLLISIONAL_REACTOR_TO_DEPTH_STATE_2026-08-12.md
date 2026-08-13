# Deterministic collisional reactor-to-depth state - 2026-08-12

## Outcome

The pure-Ar stack now runs deterministically from absorbed reactor inputs to a
common wafer boundary:

`absorbed bulk power + pressure + geometry + gas temperature`
-> Lee-Lieberman particle/power balance
-> `ne, Te, Ar+ Bohm flux`
-> delivered bias-power RF closure
-> finite-transit RF sheath
-> source-backed Ar+--Ar elastic/CX collision-order transport
-> resolved Ar+ and unscattered-fast-Ar lower-bound wafer distributions.

No particle Monte Carlo is used. Phase, thermal velocity, collision position,
impact parameter, and collision azimuth are fixed quadratures. The retained
discrete operator has an analytic neutral-density JVP, including at zero gas
density. Initial global azimuth is eliminated by symmetry and restored only at
the output boundary; collision azimuth is retained in full because the wafer
normal breaks the apparent +/- pairing.

This path is opt-in. The certified collisionless transfer is unchanged.
Unsupported CFx+ use fails rather than borrowing Ar cross sections.

## Physics gates

- Collisionless density limit reproduces the finite-transit sheath plus the
  correct transverse thermal energy.
- Every equal-mass event closes kinetic energy to machine precision.
- Ion arrival, escape, and collision-order truncation close a probability
  ledger to machine precision.
- Every born fast-neutral lineage closes into wafer arrival, subsequent
  neutral-collision truncation, or escape.
- The density JVP agrees with finite differences and is defined at density 0.
- The Khrabrov-Kaganovich 1 keV geometry statements are reproduced: inferred
  0.1 degree COM angle at 5.4 a0 and independent Rcx within 4%.
- The charge-label swap reproduces the source's theta/2 laboratory angle for
  the born fast neutral.
- Existing RF-sheath tests pass unchanged.

The executable audit is
`scripts/audit_collisional_sheath_and_depth_scope.py`; its frozen receipt is
`results/curated/collisional_sheath_depth_scope/audit.json`.

## Convergence

For the audit's declared 1 keV, 10 mm, 500 K pure-Ar sheath:

- At 1 mTorr, collision order 3 leaves 1.44e-5 unresolved ion probability and
  changes mean resolved energy by 0.0156% from order 2: the declared gate
  passes.
- At 10 mTorr, collision order 3 leaves 5.30% unresolved and changes mean
  resolved energy by 3.59% from order 2: the declared gate fails.

The path-tree expansion is therefore certified only where its convergence
receipt passes. At optical depth near one, the correct next numerical rung is
a sparse deterministic discrete-ordinates linear Boltzmann solve, not a deeper
exponential tree. That operator can reuse the same cross sections, boundary
contract, ledgers, and derivative tests.

## Krueger depth verdict

This build does not change the frozen Krueger endpoint: 346.833 nm simulated
versus 825 nm experimental.

Applying the new sheath operator to Krueger's feature input would be physically
wrong. Figure 4(a) is already the HPEM result at the wafer plane, after reactor
and sheath transport. Post-processing it through another sheath would
double-count collisions. It is also one aggregate distribution of all positive
ions; Ar cross sections cannot be assigned to its unknown molecular-ion mix.

The Figure 4 digitization was visually rechecked from a 400 dpi render. More
decisively, the newly obtained author PDF and its embedded Figure 4 PNG match
the digitizer's archived SHA-256 values exactly. The mean 3465.11 eV and signed
angular standard deviation 0.8332 degrees are source-pixel-bound. The downloaded
thesis also matches the earlier Figure 6.17 vision/PIL audit exactly; that plot
still reports aggregate `Ions` and no stable C4F6 parent flux.

Therefore the remaining Krueger depth-identifying inputs are still:

1. species-resolved positive-ion wafer fluxes and IEADs;
2. stable C4F6 wafer flux;
3. molecular-ion/neutral collision cross sections for any independent reactor
   reproduction; and
4. the reactor voltage/current waveform or a validated circuit closure.

Selecting any of these by agreement with 825 nm would be a depth fit, not a
prediction.

## What the full stack now is, and is not

The new pure-Ar route closes absorbed knobs to a resolved collisional wafer
boundary. It does not yet close generator forward power through a matching
network, a self-consistent moving sheath, subsequent fast-neutral collision
trees at high optical depth, fluorocarbon plasma chemistry, or a validated
C4F6 surface boundary. Its immutable result therefore keeps both
`supports_equipment_prediction` and `supports_feature_depth` false.

That boundary is deliberate: the architecture needed for a fast,
differentiable knobs-to-depth predictor now exists, while the evidence gaps
that prevent a Krueger depth claim remain visible rather than hidden in a
calibration.
