# Guo/Kwon–Krüger deterministic-extruded prefix gate

Frozen on 2026-08-06 after the 10 nm, 0.0625 s-rung completed and before the
0.03125 s-rung completed. No experimental endpoint, yield, flux normalization,
surface capacity, or target-dependent scale may select a numerical rung.

## Authority under test

Krüger's trench and published source are translationally invariant in the
periodic direction. The candidate neutral-exchange authority is therefore
`deterministic_extruded_2d`: Hottel crossed strings with exact projective
occlusion, reciprocity by construction, and fail-closed extrusion
certification. The general-3-D eight-ray scrambled-QMC operator is excluded
from moving-profile authority because a frozen-checkpoint control isolated an
inverse-receiver-area flux-density spike on a vanishing mask triangle.

Both rungs must start at `t=0` with:

- `surface_model=guo_tml`;
- the published Krüger boundary with ion normalization exactly `1`;
- `radiosity_backend=deterministic_extruded_2d`;
- `exchange_method=analytic_occlusion`;
- `surface_state_remap_backend=common_refinement`;
- `topology_change_policy=continue_gas_cavity`;
- 10 nm spacing and a 0.5 s physical horizon.

The nominal time steps are 0.0625 and 0.03125 s. The refined trajectory is the
comparison clock; coarse metrics are linearly interpolated only for scoring,
never for state continuation.

## Time gates

Over 0.0625--0.5 s:

- terminal depth absolute relative difference: at most 2%;
- terminal mask-opening absolute relative difference: at most 2%;
- maximum matched depth absolute relative difference: at most 5%;
- maximum matched mask-opening difference: at most 5 nm;
- exact material ledgers;
- neutral-radiosity relative balance residual at most `1e-9`;
- extrusion-projection deviation at most `1e-9` mesh units;
- raw-face to grid-resolved maximum-speed ratio at most 2 at every step;
- zero rejected trials, accepted topology events, and asymmetric cells.

Passing authorizes the 0.0625 s rung as the 10 nm temporal production
candidate. Failure holds the endpoint and requires a finer rung or an
interface/remap repair.

## Space gates

Only after the time gate passes, compare fresh 10 and 5 nm trajectories using
the same 0.03125 s nominal step through 0.5 s:

- terminal depth absolute relative difference: at most 5%;
- terminal mask-opening difference: at most 5 nm;
- maximum matched depth absolute relative difference: at most 7.5%;
- maximum matched mask-opening difference: at most 7.5 nm;
- the conservation, radiosity, extrusion, speed-ratio, retry, topology, and
  symmetry gates above remain unchanged;
- native-resolution profile images must show one connected exterior gas
  region, a connected floor/mouth, and no disconnected one-cell artifact.

A passing short prefix authorizes a no-fit 60 s **published-boundary
sensitivity forecast**, not a Tier-A absolute-depth prediction. Krüger's
unpublished ion composition and stable C4F6 wafer flux, the transfer beyond
Guo's measured energy/chemistry board, and the mask material closures remain
physical—not numerical—uncertainties.
