# ViennaPS fluorocarbon chemistry integration audit

Date: 2026-07-15

## Result

Petch now has a selectable, dimensional implementation of the La Magna--Garozzo three-coverage
fluorocarbon model used by ViennaPS 4.6.1. It runs inside the existing 3-D boundary-transport,
surface-state, material-ledger, level-set, reflection, and charging architecture; it is not a forked
feature solver. The prior reduced SiO2 chemistry remains available and remains the default in existing
campaign adapters, so this increment does not silently change archived operators.

This is an implementation and engine-integration milestone, **not an experimental validation claim**.
The transferred ViennaPS parameter set is explicitly marked nonpredictive outside its sourced domain.

## Sources and frozen parity target

- Scientific model: A. La Magna and G. Garozzo, *Journal of The Electrochemical Society* 150 (2003),
  DOI `10.1149/1.1602084`.
- Behavioral source: `ViennaTools/ViennaPS`, file
  `include/viennaps/models/psFluorocarbonEtching.hpp`, commit
  `2956ed587984c6dc38be24c6e2390e10c9b2f0a7`.
- ViennaPS source license: GPL-3.0. The exact source commit and URL are embedded in every mechanism
  provenance block.

## One-engine data path

```text
Jeong plasma-model radical densities
             |
             +-- CF + CF2 + CF3 ----------> FC_etchant flux
             |
             +-- C2F4 + C3F6 + C4F7 ------> FC_polymer flux
                                               |
Ar+ boundary --> common 3-D transport ----------+
                                               v
                              La Magna algebraic coverages
                              theta_e, theta_p, theta_pe
                                      |
                       +--------------+---------------+
                       |                              |
                 SiO2 recession                polymer growth/removal
                       |                              |
                       +---------- signed velocity --+
                                      |
                         common level-set evolution
                                      |
                 conservative film/material ledgers + remap
```

The same transported Ar+ event measure can therefore carry common-engine charged deflection and
certified grazing reflection before chemistry evaluates its energy/angle-dependent yields.

## Implemented model

The new `LaMagnaGarozzoFluorocarbonMechanism` implements:

- etchant coverage, polymer coverage, and etchant-on-polymer coverage;
- thermal chemical removal, ion-enhanced removal, and physical sputtering;
- the Vienna square-root threshold energy laws and incidence-angle factors;
- polymer saturation, deposition, ion removal, and a finite areal film inventory;
- removal of an existing film before substrate recession;
- per-step SiO2 and fluorocarbon-film material ledgers;
- neutral/surface fixed-point participation through the existing radiosity loop;
- exact ViennaPS 4.6.1 neutral-transport behavior as a selectable compatibility mode, including the
  source's use of material `beta_e` for both neutral particle types;
- a separately declared species-specific transport mode that uses `beta_p` for polymer neutrals.

The ViennaPS normalized `1e-6` zero guards are represented by the exact dimensional zero limit in v1.
This only changes branch behavior at unresolved vanishing flux and is disclosed rather than called
bitwise source parity.

## Moving-surface contract repaired

The integration exposed a pre-existing type error in the generic remap assumption: algebraic coverage
fractions were treated like conserved areal inventories. A saturated coverage on a contracting mesh
could therefore exceed the new mesh's artificial integral "capacity" and refuse an otherwise resolved
step.

The common remapper now supports an explicit per-field mode:

- `intensive`: interpolate bounded algebraic coverages;
- `conservative`: preserve area-integrated physical inventory exactly.

Legacy states default to `conservative`, so their operator is unchanged. The new fluorocarbon state
declares its three coverages intensive and its polymer-film/removed-material inventories conservative.
Material-ID routing preserves the declaration. Mechanism-owned polymer growth is combined with etch
recession in the common signed face velocity, and external product redeposition continues to use the
same path.

## Jeong 2023 adapter

`build_jeong_2023_boundary_state` now has two explicit modes:

- `aggregate` (default): the archived `FC_total` development closure;
- `heavy_light`: `FC_etchant = CF + CF2 + CF3` and
  `FC_polymer = C2F4 + C3F6 + C4F7`.

Each channel is the sum of the species-specific one-way thermal fluxes. Its representative mass is
flux weighted. Tests require the two split channels to sum back to the aggregate flux. The source data
remain digitized outputs of Jeong's volume-averaged plasma model, not measured feature-boundary fluxes,
so the boundary still refuses a predictive-evidence label.

`scripts/jeong_2023_transfer.py --chemistry-model lamagna_garozzo` selects the new chemistry, split
boundary, and common neutral/surface fixed point. `--chemistry-model reduced_si_o2` preserves the old
path and remains the default.

## Verification evidence

- Manufactured equations reproduce all three analytic coverages and all removal-rate components.
- Ion-only, deposition-only, film-depletion, and undeclared-species limits pass.
- Polymer deposition and removal ledgers close exactly.
- Material routing preserves both recession and mechanism-owned growth velocities.
- A common-engine 3-D deposition test moves the shared level set in the growth direction.
- A manufactured contracting-surface remap keeps saturated coverage bounded while conserving finite
  film inventory.
- Safe charging checkpoints round-trip the new state.
- The legacy focused surface/material/feature gates remain passing.
- A coarse end-to-end Jeong adapter smoke completed 80 moving-profile steps in about 20 seconds.

The coarse smoke predicted net deposition (`-82 nm` etch depth) for the transferred ViennaPS defaults,
where the experimental anchor is about `1223 nm` etch depth. This is retained as a useful negative
result: at that Jeong support point the grouped polymer thermal flux is about 59 times the grouped
etchant flux, so the uncalibrated transferred parameter set saturates polymer coverage. It proves the
new path executes; it does **not** validate the parameter transfer, the radical role closure, or Jeong.
No held-out datum was tuned.

## What this closes, and what remains

Closed by this increment:

1. The generic ViennaPS La Magna fluorocarbon reaction structure is available in petch.
2. Polymer deposition can move a petch feature without a second geometry engine.
3. Jeong's resolved radical species can reach distinct etchant/polymer channels without flux loss.
4. The chemistry can compose with petch's richer charging and reflected-ion transport.

Still required before claiming a strict predictive superset:

1. Direct numerical parity fixtures against ViennaPS for matched geometry, flux normalization, and
   discretization; analytic equation parity alone is not profile parity.
2. A preregistered one-anchor calibration/development step followed by untouched Jeong held-out scoring,
   with grid/ray/timestep uncertainty.
3. Evidence for the radical role/weight closure. Species density is not automatically reactive or
   depositing flux; effective reaction probabilities may need bounded reactor/surface evidence.
4. Charging-on causality tests only where charge-off residual morphology earns them.
5. Other ViennaPS process models (for example its distinct CF4/O2 and HBr/O2 parameterizations) if
   "strict superset" is meant across the entire ViennaPS catalog rather than this FC model.

Accordingly, the earned claim is: **petch can run the ViennaPS-class generic fluorocarbon mechanism in
its unified, more physics-capable engine.** The unearned claim is: **petch has validated Jeong or is
already more accurate than ViennaPS on experimental profiles.**
