# Worktree cleanup map

Date: 2026-07-16

This is a non-destructive inventory of the current local work. It does not declare unfinished
science complete, and it does not authorize deleting raw campaign evidence.

## Git state

- Active branch: `codex/unified-engine-root-fixes`
- HEAD and upstream: `dcbfa5f`, exactly aligned with
  `origin/codex/unified-engine-root-fixes` before the uncommitted work below.
- The branch is 189 commits ahead of `origin/main`; do not use `main` as a cleanup scratchpad.
- A second clean worktree exists at `/private/tmp/petch-adaptive-cascade` on
  `codex/adaptive-cascade-engine`.
- That branch has three differently hashed commits whose equivalent changes are already present on
  the active branch (`04a97e7`, `9bb882c`, and `eda4e1c`). It contains no unique unmerged work.
- There are no stashes.

## Inventory

At inventory time:

| Class | Count | Approximate size | Meaning |
| --- | ---: | ---: | --- |
| Modified tracked files | 53 | +8,648 / -637 lines | Substantial engine, tests, scripts, and docs |
| New engine modules | 15 | included in 0.9 MB below | Durable implementation |
| New tests | 18 | included in 0.9 MB below | Durable verification |
| New scripts | 23 | included in 0.9 MB below | Campaign and audit tooling |
| New experimental-data files | 12 | included in 0.9 MB below | Checksum/protocol evidence |
| New root audit documents | 7 | about 96 KB | Scientific decisions and closure reports |
| New result files | 1,669 | about 174 MB | Mixed curated evidence, raw replicates, and checkpoints |
| Generated `build/` files | 67 | about 1.7 MB | Disposable package build output |

The repository occupies about 310 MB. `results/` accounts for about 199 MB; source code is about
11 MB. The mess is therefore mostly artifact lifecycle, not an intrinsically oversized engine.

## Durable implementation clusters

### 1. Charging, electrostatics, and unattended runtime

Primary tracked changes:

- `src/petch/boundary_transport_3d.py`
- `src/petch/charged_surface_cascade_3d.py`
- `src/petch/charged_surface_response_3d.py`
- `src/petch/charging_coevolution_3d.py`
- `src/petch/charging_coupled_3d.py`
- `src/petch/charging_poisson_3d.py`

New durable modules:

- `src/petch/charging_checkpoint_3d.py`
- `src/petch/charging_stationarity_3d.py`
- `src/petch/conductor_terminal_3d.py`

This cluster contains float64 certification/replay, adaptive trajectory horizons, compatible-Q1
charge projection, periodic field topology, stochastic terminal windows, checkpoint/restart,
external conductor terminals, and recovery/error accounting. Its matching tests and bounded audit
scripts belong with it.

### 2. Unified profile, material, and chemistry engine

Primary tracked changes:

- `src/petch/feature_step_3d.py`
- `src/petch/surface_kinetics.py`
- `src/petch/physical_api.py`
- `src/petch/surface_charge_remap_3d.py`

New durable modules:

- `src/petch/material_mechanism_3d.py`
- `src/petch/surface_product_redeposition_3d.py`
- `src/petch/chlorine_poly_si.py`
- `src/petch/fluorocarbon_lamagna.py`
- `src/petch/silicon_sf6o2.py`

This cluster contains signed material motion, material-ID routing, conservative state/product
ledgers, same-material redeposition, the common charging/profile API, and process-specific surface
mechanisms. Its matching tests belong with it.

### 3. Experimental boundaries, observables, and campaign contracts

Primary tracked changes:

- `src/petch/experimental_boundary.py`
- `src/petch/experimental_data.py`
- `src/petch/sheath.py`

New durable modules:

- `src/petch/reactor_boundary.py`
- `src/petch/notching_validation_3d.py`
- `src/petch/nozawa_replay_3d.py`
- `src/petch/hwang_giapis_scatter_3d.py`
- `src/petch/physical_arrivals_3d.py`
- `src/petch/profile_observables_3d.py`
- `src/petch/twist_campaign_3d.py`

This cluster contains provenance-bound reactor/diagnostic adapters, waveform sheath transport,
Nozawa/Hwang replay, held-out scoring contracts, stochastic arrivals, and geometry-native notch,
bow, and twist observables. The new experimental CSV/JSON files, validation scripts, and matching
tests belong with it.

### 4. Product surface and documentation

- `src/petch/__init__.py` exports all three implementation clusters.
- `pyproject.toml` installs the Nozawa command and its bundled evidence.
- `README.md`, capability maps, verification contracts, HTML documentation, and root audit reports
  describe the same work.

These files should be committed only after the underlying clusters and tests are fixed in place so
the published surface cannot get ahead of the engine.

## Result-artifact classification

The 1,669 untracked result files are not one commit.

1. **Curated evidence:** final `audit.json`, `summary.json`, configuration, small plots, and explicit
   validation-ledger outputs referenced by reports.
2. **Restart evidence:** the smallest checkpoint needed to reproduce or continue an accepted run.
3. **Raw campaign evidence:** per-seed/per-level current audits and intermediate checkpoints. Keep
   them until a checksum-bound external archive exists, but do not put all of them in ordinary Git.
4. **Operational ephemera:** heartbeats, logs, plotting caches, superseded retries, and generated
   build copies. These are not scientific deliverables.

No raw result is deleted during this cleanup. Promotion from class 3 to class 1 must be explicit and
must preserve the source configuration, Git revision, seeds, checksums, and the conclusion that used
the artifact.

## Safe checkpoint sequence

1. Run the focused subsystem tests, then the complete test suite on the untouched scientific state.
2. Commit the charging/electrostatics/runtime cluster with its tests and only its required small
   evidence.
3. Commit the unified material/chemistry/profile cluster with its tests.
4. Commit experimental boundary/data/campaign contracts with their tests and checksum-bound source
   data.
5. Commit product exports, packaging, documentation, audit reports, and a curated validation ledger.
6. Build a checksum index for raw campaign directories before moving any of them to external/local
   archival storage.
7. Start reactor-scale implementation on a new branch from that checkpoint. Keep the reactor solver
   upstream of `PlasmaBoundaryState`; do not merge chamber and feature meshes into one monolith.

Because several tracked files span clusters, use explicit path staging and, only where necessary,
interactive hunk staging. Every commit must pass `git diff --cached --check` and its relevant tests.

## Hygiene rules going forward

- `build/`, `dist/`, per-run `.plot-cache/`, and `results/_scratch/` are ignored.
- New exploratory jobs write below `results/_scratch/<campaign>/<run-id>/`.
- A result leaves `_scratch` only when an audit names the exact files that support a conclusion.
- Campaign runners write manifests/checksums; reports refer to those manifests rather than copying
  raw arrays into Git.
- One implementation campaign per branch/worktree. A second worktree is used only when the change
  can be merged as a bounded commit, not as an alternate accumulating universe.
- Do not delete the current raw result tree until its archive/index is verified.
