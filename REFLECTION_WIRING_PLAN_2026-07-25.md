# Grazing-ion reflection wiring plan (single change; audit P0.2)

Goal: primary ion events at grazing incidence split off reflected weight as a
single-bounce hot-neutral population; residuals (opening 18.8 vs 45, mask-top
erosion 133 vs ~0) are the target observables. Default OFF (preserves replays);
mixed-layer campaign turns it on explicitly.

## Mechanism (GrazingSpecularIonReflection3D.literature_bounded_sensitivity)
- P_reflect(cos) = 0.95 * (1 - cos^3); energy retention 0.90; specular
  direction d' = d - 2(d.n)n about the LOCAL face normal.
- Single bounce v1: secondary hits deposit fully (declared omission:
  multi-bounce tail, ~0.95^2*grazing^2 measure, recorded not modeled).

## Implementation (post-processing, tracer untouched)
1. New function in boundary_transport_3d.py:
   `split_grazing_reflection(population, verts, faces, normals, model, *,
   periodic_lateral, domain)`
   - Input: primary ion FaceResolvedEnergeticFlux (has event_position and
     event_incident_direction — REQUIRED; refuse if absent).
   - w_refl = P_reflect(event_cosine) per event; primary event flux scaled by
     (1 - w_refl) in place (new population object).
   - Secondary rays: origin = event_position + eps*normal[event_face],
     direction = specular; energy = 0.90 * event_energy.
   - Trace with a plain warp raycast (wp.mesh_query_ray, same mesh build as
     _ions_deterministic) with periodic lateral wrap; hits become a second
     FaceResolvedEnergeticFlux named f"{name}:hot_neutral" (charge 0);
     flux = w_refl * event_flux * area[src]/area[dst] measure bookkeeping
     (event_flux is per-hit-face density: convert via event_rate = flux*area).
   - Escapes (upward through source plane): allowed, weight recorded in
     diagnostics dict returned alongside.
2. Call site: in the pilot path where ion populations are finalized in
   gather_boundary_state_ballistic_3d face_gather mode (after primary gather,
   before chemistry roles). Gate on new kwarg
   `grazing_ion_reflection_model=None`.
3. Pilot: flag --grazing-ion-reflection {off,literature_v1}; role map gains
   the hot-neutral name -> "energetic_bombardment"; config records model
   provenance.
4. Adapter: nothing — it already folds any FaceResolvedEnergeticFlux into the
   flux-weighted mean per face.
5. Gates (tests/test_grazing_reflection_split.py):
   - normal incidence (cos=1): zero reflected weight, primary unchanged bitwise.
   - grazing event on a vertical wall of a box trench: reflected weight
     appears on the opposite wall/floor with 0.9 energy; total particle rate
     conserved (primary+secondary+escape == original) to 1e-12.
   - specular direction correctness on an inclined face (analytic).
   - periodic wrap: reflection crossing the lateral boundary lands.
6. Run: box, ml7-base-refl 60s all-published + reflection on. Gate:
   opening in [40,50], mask-top erosion < 40 nm, depth within +/-5% of 825.
