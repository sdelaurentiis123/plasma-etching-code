# ARDE legacy-to-unified reconciliation

Date: 2026-07-15. Scope: local git history and current worktree only. No remote operation was used.

## Why this audit was necessary

The repository contains two feature-scale engine generations. The legacy compatibility engine
(`petch.threed.run_etch_3d`) accumulated substantial SF6/O2, neutral-transport, reflection, and
profile-evolution work before the dimensional common engine (`petch.feature_step_3d`) existed. A new
ARDE implementation must therefore begin by reconciling those paths, not by reimplementing either one.

The audit also found an evidence error that changes how the historical ARDE campaign may be used. The
sequence `1.00 / 0.43 / 0.29 / 0.20` at aspect ratios `0 / 10 / 20 / 40` is an evaluated
Blauw/Clausing model curve, not a direct digitization of de Boer Figure 9. Commits `41dee3c`,
`9c3e181`, `43083b1`, `7ff2dcb`, `8e937dc`, and `0e61e76` contain useful engineering, but their
"wafer pass" or "held-out AR40" language is not experimental validation. The directly digitized,
checksummed Figure-9 markers in `data/experimental/deboer_2002/` are now development data; a different
experiment is required for a future validation claim.

## History inventory

| Capability | Historical implementation / commit | Current common-engine status | Decision |
| --- | --- | --- | --- |
| Coupled F/O Belen kinetics | `src/petch/belen.py`, introduced in `387a5f0` | Reimplemented dimensionally as `BelenSiliconSF6O2Mechanism`, with parameter evidence/bounds, exact energetic-event integration, neutral/surface fixed-point coupling, and a target-material ledger | Use the common implementation. Keep the legacy function for replay only |
| Evolving 3-D level-set profile | `src/petch/threed.py::run_etch_3d`; full evolving campaigns in `43083b1` | `solve_feature_3d` already performs multi-step material-local motion, conservative state remap, topology checks, and validity reporting | Do not write a new profile driver. Use `solve_feature_3d` / `PhysicalProcess` |
| Deep neutral transport | MC, DDA, radiosity, and 1-D Knudsen paths; numerical fixes through `f8dced1` and `43083b1` | Conservative material/state-dependent diffuse radiosity is already authoritative. It closes source = reacted + escaped and couples to the same surface state | Do not port the calibrated `knudsen_wall_loss_scale` as governing physics. Retain Knudsen as a legacy reduced-model/replay path; reuse its bottom-area and nearest-valid-slice lessons in tests and diagnostics |
| Rounded-front Knudsen sink fix | `43083b1`: bottom-classified, area-weighted sink; avoids applying a full duct loss in every slice crossed by a rounded front | Common radiosity uses actual triangle areas and per-face reaction probabilities, so it does not contain that slice-multiplication bug | Add a regression statement/test where appropriate; no new sink formula is needed |
| Passivated sidewalls | `43083b1` bottom-only sink; opt-in depth-shaped `knudsen_front_loss` in `7ff2dcb` | Belen available-site coverage controls the neutral reaction probability on every face without an aspect-ratio branch | Keep the state-derived common law. Do not port the AR-shaped front-loss closure |
| Energetic ion reflection | Legacy faithful coned-cosine path, wired into deterministic transport in `9c3e181` | The common engine now applies one certified, charge/energy-accounted, adaptive full-lineage cascade in both charged co-evolution and ordinary profile evolution. Response-enabled straight trajectories reuse the certified zero-field tracer so exact impact position/direction and selective float64 replay are preserved | Migration complete; retain the legacy kernel for replay only |
| Redeposition/passivation | Rate-suppression closure in `9f1c1cf` | Conservative emitted/deposited/escaped product transport and same-material growth already exist. Cross-material films deliberately refuse | Use the common ledger. Do not port the legacy velocity-suppression shortcut |
| Charging/profile co-evolution | Multiple legacy 2-D/closure paths | `solve_charging_coevolution_3d` already runs physical-time charge, certified reflection/SEE, charge remap, and common profile motion | Keep this as the charged path. Do not route new work through legacy charging closures |
| Photon-assisted etching | No historical feature-etch implementation was found in branches, reflogs, or reachable commits | Explicitly absent from the Belen mechanism and provenance | Not a lost feature. Add only when a case supplies a measured photon flux/spectrum and material-specific photon-assisted yield/state law |

## Object-store and branch check

The relevant history is linear and reachable from the current branch. The only unreachable commit found
by `git fsck --full --unreachable --no-reflogs` is `47e05b0`, the pre-amend version of `df28b29`; the
amendment changes only the float64 replay refinement ceiling from four to five halvings. It contains no
lost ARDE, chemistry, photon, or profile implementation. The `codex/adaptive-cascade-engine` worktree
contains the already-cherry-picked cascade/horizon work plus its terminal-window commit lineage; it does
not contain a separate ARDE implementation to merge.

## Bounded implementation result

1. `apply_charged_surface_response_to_transport_3d` now applies a declared response to one existing
   face-resolved primary measure and returns the completed cascade plus chemistry-facing reimpacts.
2. `advance_feature_step_3d` and `solve_feature_3d` expose that same operator. In a field-free case only
   charged species use the certified zero-field tracer; neutral species retain the cheaper first-hit path.
3. Manufactured and multi-step gates cover lineage, response causality, charge/energy closure, large-event
   reduction order, periodic Godunov advection, and duplicate-endpoint seam closure. The full repository
   suite passes (460 passed, one CUDA-only skip on the local CPU build).
4. The source-correct de Boer development diagnostic with the literature-bounded reflection sensitivity
   improved the twelve-point development RMSE from 3.545 to 2.557 micrometres. Response charge errors were
   at or below 2.9e-14; energy closed exactly. This is improvement, not validation or a complete fit.
5. A three-level moving-profile gate completed through the same operator. Its base/refined/fine final
   centerline depths were 0.218634/0.217613/0.217393 micrometres from a 0.200000 micrometre start. The
   static-rate counterfactual ended at 0.218577 micrometres. The latest sampling change was 0.000220
   micrometres (4.6x contraction from the preceding change), while moving-minus-static was -0.001184
   micrometres. All response, cascade-tail, and state-remap conservation gates passed.

The moving gate also exposed and closed a real common-engine defect: periodic transport had been coupled
to nonperiodic level-set motion. A reflected ion could wrap across a subcell mismatch between duplicate
endpoint planes and encounter the periodic image from its solid side. Periodic profile evolution now uses
wrapped Godunov neighbors, wrapped-padding redistancing, and an explicitly reported seam projection. No
visibility or lineage tolerance was weakened.

## Claim boundary

This migration can establish that the common engine contains the already-earned mechanisms and that the
direct de Boer development miss is not caused by an omitted reflection wire. It cannot turn Figure 9 into
a predictive validation because the paper does not declare O2 flow, electrode temperature, mask profile,
ion flux, IEDF/IADF, or statistical measurement uncertainty for that figure. Photon physics remains a
separate, evidence-triggered channel rather than an unmeasured repair knob.
