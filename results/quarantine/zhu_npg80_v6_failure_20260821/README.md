# Oxford NPG80 v6 failure quarantine

This directory preserves the remote-only evidence copied from Vast instance
`48177892` after the clean v6 square-pillar board exited on 2026-08-21.

Source paths:

- `/root/zhu_v6_board_980800a.log`;
- `/root/petch-4b656fd/results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/`.

The remote trajectory directory contained 175 historical caches across model
revisions v3-v6. Only the 23 files whose embedded `job_spec.model_revision` is
`two-material-moving-tio2-cr-owner-projection-v6` are preserved here under
`v6_trajectories/`. They are quarantine evidence, not yet a complete or
gradeable 56-cell board. `SHA256SUMS.v6` covers those 23 files plus the complete
parent failure log.

The parent failed on the missing v6 cell with width 120 nm, selectivity
18.016664028610727, high ion energy, and zero angular tail. The terminal chain
is:

```text
advance_feature_step_3d
  -> _remap_surface_state_with_indexed_transfer
  -> build_surface_transfer_3d
  -> SurfaceTransferWeights3D.__post_init__
  -> ValueError: invalid sparse surface-transfer weights
```

The aggregate v6 parent did not write a failed-cell checkpoint. A narrow
reproducer was launched separately after commit `48b78d4` added diagnostics
that name the exact rejected invariant without changing acceptance criteria.
Do not promote these caches into the curated board until that reproduction
determines whether the numerical change is diagnostics-only/bitwise-equivalent
or requires a model-revision bump and full recomputation.

The unfiltered local copy of all 175 remote files is intentionally not part of
this evidence package because it duplicates already committed historical
caches. It may remain as local untracked working data during diagnosis.
