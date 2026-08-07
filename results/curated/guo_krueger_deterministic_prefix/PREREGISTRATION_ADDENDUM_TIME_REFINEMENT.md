# Addendum: production time step for the spatial gate

Frozen on 2026-08-06 after `time_gate_final.json` passed and before any 5 nm
deterministic-extruded trajectory was launched or inspected.

## Why an addendum is required

`PREREGISTRATION.md` fixed the first 10 nm pair at 0.0625/0.03125 s and
specified 0.03125 s for the subsequent 10/5 nm spatial comparison. The
first pair failed one frozen early-trajectory depth criterion, and the next
0.03125/0.015625 s pair failed the same criterion. The authorized finer
0.015625/0.0078125 s pair passed every original time, conservation, solver,
and symmetry gate. Therefore 0.015625 s, the coarse member of the passing
pair, is the production time step for spatial refinement.

Using 0.03125 s after it failed its twofold time comparison would mix known
temporal error into the spatial error. Using 0.0078125 s would add cost
without being the production member selected by the frozen twofold rule.

## Frozen spatial comparison

Run fresh `t=0` trajectories through 0.5 s at:

- reference: 10 nm spacing, 32 nominal steps, `dt = 0.015625 s`;
- refined: 5 nm spacing, 32 nominal steps, `dt = 0.015625 s`.

Both use:

- `surface_model=guo_tml`;
- published Krüger boundary with ion normalization exactly 1;
- `radiosity_backend=deterministic_extruded_2d`;
- `exchange_method=analytic_occlusion`;
- `exchange_geometry_tolerance=1e-9`;
- `exchange_relative_tolerance=1e-5`;
- `radiosity_tolerance=1e-12`;
- `radiosity_max_iterations=2000`;
- `surface_state_remap_backend=common_refinement`;
- `topology_change_policy=continue_gas_cavity`;
- `n_position=16`.

The original spatial thresholds remain unchanged:

- terminal depth absolute relative difference at most 5%;
- terminal mask-opening difference at most 5 nm;
- maximum matched depth absolute relative difference at most 7.5%;
- maximum matched mask-opening difference at most 7.5 nm;
- exact material ledgers;
- neutral-radiosity relative balance residual at most `1e-9`;
- extrusion-projection deviation at most `1e-9` mesh units;
- raw-face/resolved-grid maximum-speed ratio at most 2;
- zero rejected trials, accepted topology events, and asymmetric cells;
- native-resolution profiles must show one connected exterior gas region, a
  connected floor/mouth, and no disconnected one-cell artifact.

This is a numerical-convergence amendment only. It was selected without the
825 nm endpoint, any target profile, a flux scale, a yield change, or a
surface-capacity change. A pass authorizes only a no-fit published-boundary
sensitivity forecast; it does not identify the missing C4F6 reactor boundary.
