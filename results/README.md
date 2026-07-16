# Result artifact policy

Simulation outputs are local runtime artifacts by default. Git tracks source code, experimental
inputs, validation contracts, audit reports, and deliberately curated evidence—not every seed,
checkpoint, heartbeat, or retry.

## Local runs

Campaigns may continue writing their established paths under `results/`. Those paths are ignored,
so an exploratory or long-running campaign does not make the source worktree dirty.

Before moving or deleting a large local campaign, run:

```bash
python scripts/index_local_result_artifacts.py \
  --output LOCAL_RESULT_ARTIFACT_INDEX_YYYY-MM-DD.json
```

The index records a checksum tree, file count, and byte count for each top-level local result
directory. It is an inventory, not a validation claim.

## Promoting evidence

Only an audit may promote a result into version control. Copy the smallest sufficient evidence set
under:

```text
results/curated/<campaign>/<revision>/
```

A curated set normally contains:

- the complete input/configuration manifest;
- Git revision and source checksums;
- seeds and sampling mode;
- a compact machine-readable summary;
- the final gate values and uncertainty/refinement evidence;
- the smallest restart or field artifact required to reproduce the conclusion, when necessary.

Raw per-replicate arrays remain in local or external artifact storage and are referenced by their
checksum-tree index. Existing result files already tracked before this policy remain tracked until a
separate, evidence-reviewed migration replaces them.
