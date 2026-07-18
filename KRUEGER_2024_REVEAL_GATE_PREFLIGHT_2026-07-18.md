# Krüger 2024 reveal-gate preflight

Date: 2026-07-18  
Status: freeze tooling ready; held-out observations remain sealed

## What was stale

The reveal tool accepted the superseded `base-axisymmetric-secant.v1` derivation but did not
understand the current `development-trust-region-proposal.v1` artifact that selected the fixed
R1.9 pair. It also trusted completed endpoint JSON files without requiring the launch manifests
that identify the clean source archive used to produce them.

The archived trust proposal was internally checksum-valid and selected the correct pair, but its
protocol hash referred to a transient pre-commit copy of the R1.9 document. It was not silently
grandfathered. The existing base-only generator was rerun against the committed protocol and the
same sealed inputs. The physical result is unchanged:

```text
effective_mask_crosslinked_growth_fraction = 0.9004722559883319
oxide_etch_yield_scale                      = 0.5586489665864749
```

The regenerated proposal is
`results/curated/krueger_2024_r19_trust_proposal_2026-07-18.json`, with embedded proposal SHA-256
`e7d10d9d2497c7abbc350073d4a01f70b78c56cf3df8a898e04f6d375a08309f`. It read only the base
calibration target and the already sealed development receipts.

## New freeze contract

Revision `28ffd45` requires all of the following before producing a reveal artifact:

1. the current R1.9 trust proposal is internally checksum-valid, bound to the committed protocol
   and base target table, and selects the exact preregistered pair;
2. the 10 nm and 5 nm endpoint audits both pass their existing completion, conservation, operator,
   and base-tolerance gates;
3. explicit 10 nm and 5 nm launch manifests are supplied;
4. both launch manifests bind the same clean git revision and source-archive SHA-256;
5. the 10 nm companion points to the exact 5 nm authority launch hash;
6. the manifests declare the 10/5 nm R1.9 refinement relationship and common-refinement remap;
7. the launched pilot, boundary transport, feature step, and common-refinement source hashes match
   the code being frozen.

The CLI now requires `--launch-10nm` and `--launch-5nm`. The paired launch receipts and shared
source epoch are copied into the reveal artifact. Mixed source epochs, tampered manifests, missing
operator checksums, or any source change after launch refuse before held-out execution.

## Real-campaign consequence

The partial `efaa070` trajectories remain valuable diagnostics, but their launch manifests now
correctly fail the freeze preflight because `feature_step_3d.py` predates the local-gradient normal
repair. The authoritative base must therefore start from zero under a new clean source epoch. This
was already required scientifically; it is now enforced mechanically.

No held-out observation was opened, no simulation was launched, and no GPU was rented. Focused
freeze tests pass 15/15; the repository suite passes 976 tests with one expected skip.
