# Overnight plan 2026-07-25 — per-event chemistry (the ml8 blocker)

Doctrine: audit-first, one verified commit per step, gates before box, ml6
(depth 811/825 zero-anchor, commit a8189aa lineage) is the protected
checkpoint — nothing may regress it.

## Why
ml8 proved the adapter's flux-weighted mean (E, cos) compression breaks when
faces receive heterogeneous populations (primary ions + oblique hot
neutrals): one bad average poisons a face through two nonlinear energy laws.
Krüger evaluates chemistry per event. This was the audit's declared omission;
it is now the blocking item.

## Design: atom-resolved ion channels in the mixed layer
1. `ModuleFluxes` gains optional sparse ion atoms: `ion_atom_face`,
   `ion_atom_flux`, `ion_atom_energy_eV`, `ion_atom_cosine` (the adapter's
   event arrays passed through). When present, `step()` computes EVERY
   ion-driven term as a sum over atoms with full state coupling
   (attenuation exp(-d_FC/lambda(E_i)) per atom, film sputter law per atom,
   Kress per atom, crosslink (E_i - E_iface_i) per atom, layer channels per
   atom); scalar path unchanged when absent.
2. Cost: atoms ~ O(100) per face; numpy sum over a sparse (atom -> face)
   scatter via np.add.at or segment sums. Acceptable at 9k faces.
3. Gates (tests/test_mixed_layer_atoms.py):
   - single-atom == scalar path to 1e-12 (same E, cos);
   - two-atom Jensen gate: rate != rate-at-mean, sign matches concavity;
   - ledger closure with atoms < 1e-9;
   - adapter integration: two EnergeticFlux populations -> atom arrays,
     no mean compression.
4. Adapter: build atom arrays from FaceResolvedEnergeticFlux events
   (event_face/flux/energy/cosine) and EnergeticFlux (broadcast per face x
   spectrum row). Remove the "spectrum compressed" omission from validity
   when atoms are active.
5. Full suite → commit.
6. Box: rerun ml6-config (NO reflection) with atoms → depth must stay in
   [780, 870] (Jensen correction may shift it; record honestly). Then
   ml9 = atoms + reflection literature_v1 → gate: opening [40,50], depth
   825±5%, mask erosion < 40 nm. Archive audits, commit, KILL BOX.
7. If ml9 passes → 8-condition scorecard overnight (o05/o15/o25/p4/p8),
   commit scorecard table. If not → localize with per-face energetic
   receipts; do NOT stack further mechanisms without analysis.

## Research (Opus agent, parallel)
How MCFPM/Kuboi/industry codes couple event-level energetic distributions to
surface-site chemistry; segment-sum patterns for GPU atom-resolved surface
kinetics; any published Jensen-bias quantification for IEAD compression.
Deliverable: RESEARCH_EVENT_RESOLVED_CHEMISTRY_2026-07-26.md (do not commit).

## Not tonight
Multi-bounce reflection, crosslinked sputter resistance [VERIFY], activation
states, redeposition — all stay parked behind the audit's P1/P2 queue.
