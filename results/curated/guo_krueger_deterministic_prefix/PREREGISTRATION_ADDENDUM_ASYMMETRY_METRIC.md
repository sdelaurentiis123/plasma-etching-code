# Addendum: resolved symmetry observable

Frozen on 2026-08-06 after the first 10/5 nm comparison was graded and before
either trajectory is rerun. This addendum changes a measurement implementation,
not the solver, boundary, chemistry, time step, geometry, or zero-tolerance
symmetry threshold.

## Preserved first result

`space_gate_first_attempt.json` is the immutable first grade. It failed only
the quantity called `maximum_asymmetric_cell_count`, whose value was
`0.004005852860336501`. The implementation divided a continuous difference
between two linearly interpolated zero-crossing offsets by `dx`. That result is
a fractional-cell interface-location residual; it is not a count of asymmetric
cells. The result is not waived or overwritten.

The first pair otherwise agreed to:

- `0.0798607%` in terminal depth;
- `0.128620 nm` in terminal mask opening;
- `0.427606%` over the full matched depth trajectory;
- `0.134785 nm` over the full matched mask-opening trajectory.

It had exact material ledgers, maximum neutral-radiosity balance residual
`2.19e-12`, maximum extrusion deviation `6.94e-18` mesh units, maximum
raw/resolved speed ratio `1.0000133`, and zero rejected trials or topology
events. Native-resolution inspection found one connected exterior gas region,
a connected floor and mouth, and no disconnected one-cell artifact.

## Corrected observable, frozen before rerun

The preregistered zero-tolerance `asymmetric cells` gate means:

1. compute the phase at every physical finite-volume cell center as the mean of
   its eight nodal signed-distance values;
2. classify the center as solid for `phi >= 0` and gas otherwise;
3. compare each x-mirror cell pair once; and
4. count pairs whose phase classifications disagree.

`asymmetry_cell_count` is this integer count and its frozen limit remains
exactly zero. Two stricter resolved checks are also frozen at zero:

- mirrored nodal solid/gas sign-mismatch pair count;
- mirrored resolved material-label mismatch pair count, with gas normalized to
  material zero.

The old continuous quantity remains in every new audit as
`maximum_subcell_interface_asymmetry_cells`, but is diagnostic only. No
continuous level-set calculation can honestly be called a cell count.

## Rerun and provenance rule

Both the 10 and 5 nm trajectories must be rerun fresh from `t=0` with the
already frozen `0.015625 s` nominal step and all physical settings unchanged.
The corrected metrics must be recorded at every accepted step. The original
first-attempt gate and its source audit hashes remain committed beside the new
receipt.

The first 5 nm directory was created at
`2026-08-07T08:30:31-0400`, immediately after commit `559c857`, and completed
at `2026-08-07T10:56:25-0400`. Across `559c857..07054aa`, the pilot blob stayed
`1f37fe5c943848724329db5cdecff28622bfd747` and the entire `src/petch` tree
stayed `4804806cd9fa3d55cc7b3913b1cbcc614afbb28b`; the intervening commits
changed research evidence only. The source audit hashes are:

- 10 nm:
  `60fb3033247ed6866b6d7192f57ef5a6dc3864fbedc34c47c814762a3bc1c539`;
- 5 nm:
  `b5fe395f35837492d974f1c4fbc4abe09ef81df373e3a2ba01a2b14695ae3666`.

No experimental endpoint, yield, flux normalization, surface capacity, or
target-dependent scale may select or alter the corrected grade.
