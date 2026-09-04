# Codex handoff — 2026-09-04

Branch: `codex/validation-first-multiphysics` at `1b1ea5e`, pushed. Working tree
clean except the untracked files listed in section 6.

## 1. Where we are, in one paragraph

Two partner cases are at the "recipient pages built, owner has not sent them"
stage. Freddie (Oxford NPG80, TiO2 square pillars) has a frozen blind prediction
that was scored against his SEMs after the freeze: the failure mechanism was
predicted correctly, the severity was over-predicted. Arun (Resona, printed
clock-gate mask) has an exact-geometry entrance-transport result and a clean
explorer page, but that page is a depth heightfield, not a profile simulation,
and the owner has asked why it shows no sidewall or charging physics. The next
engineering step on both cases is the same kind of thing: run real feature-scale
profile physics through the common engine (`src/petch/feature_step_3d.py`)
where the current deliverables stop at surface laws.

## 2. Freddie / Oxford TiO2 — state

Recipe: CHF3/SF6/O2 55/5/1 sccm, 30 mTorr, 150 W, 20 C, 1200 s, 700 nm ALD TiO2
under 45 nm Cr, square pillars at 350 nm pitch, nine dose blocks.

Frozen boards (`FROZEN_OXFORD_PREDICTION_2026-08-31.md`, commit 378f3fa):

- `results/curated/zhu_npg80_moving_cr_profiles_v1/` — v7, 112 endpoints,
  audit.json check PASS.
- `results/curated/zhu_npg80_gds_square_profiles_v1/` — 368 endpoints,
  widths 105–250 nm, audit.json check PASS, envelopes rendered.
- `target_sem_used: false` in both audits. Nothing was re-run after the SEMs
  were opened.

Answer key (`data/experimental/zhu_2026_tio2_npg80/sem_0817/`, commit 82fad37):
80 Zeiss TIFs, manifest, per-dose measurements. Scripts:
`scripts/digitize_zhu_npg80_sem_0817.py`, `scripts/measure_zhu_npg80_sem_lattice.py`,
`scripts/measure_zhu_npg80_sem_height.py`, `scripts/measure_zhu_npg80_sem_occlusion_height.py`.
None of these four scripts has a test yet.

Score (`showcase/oxford_scorecard/oxford_tio2_blind_scorecard.html`):

| Quantity | Predicted | Measured | Verdict |
|---|---|---|---|
| Pitch | 350 nm design | 334–346 nm | consistent |
| Mid-height CD | 168–191 nm (w185/w195) | 181–236 nm eq-square by dose | consistent |
| Cr cap | exhausted at every width | loss zones in doses 3, 7, 8, 9; intact elsewhere | mechanism hit, severity over-predicted |
| Capillary collapse | none (8–11x margin) | none observed | hit |
| Pillar height | 372–406 nm | 0.25–0.55 um from tilt occlusion (dose 9: 533 nm, dose 7: 248 nm) | geometry estimate only |
| Trench depth | 673–676 nm | no cross-section | unscored |

Physics gaps this exposed (detail in `PHYSICS_AND_REPO_AUDIT_2026-09-01.md`):

1. Cr removal is a rate-normalized law, not surface chemistry. The board says
   the cap is gone by ~11 min (ion-enhanced, from the frozen endpoint's top
   loss); the blanket Cr rate says 14.5–23.8 min. The SEMs say it depends on
   dose block. A literature-first Cr/CHF3-SF6-O2 surface deck validated on
   published blanket Cr rates is the fix. Do not fit to Freddie's SEMs.
2. Reactor-scale uniformity was assumed; the dose-block dependence of Cr loss
   is the first sign the answer is not uniform across the 4-inch wafer.
3. Absolute TiO2 rate is interval-valued (34.1–43.5 nm/min from Janissen
   witnesses) until Freddie supplies a same-run blanket loss and self-bias.

Recipient pages (all self-contained HTML, safe to email):

- `showcase/oxford_client/oxford_tio2_prediction_vs_sem.html` — one-to-one
  predicted pillar vs SEM.
- `showcase/oxford_scorecard/oxford_tio2_blind_scorecard.html` — full scorecard
  with TL;DR.
- `showcase/oxford_interactive/oxford_tio2_interactive.html` — whole-board
  interactive: block picker, rotatable pillar mesh, ion rain with per-surface
  impact tallies, Cr cap thinning to t_c then top attack, scenario buttons
  (147/296 eV, tail 0/0.65, TiO2:Cr 14/18). Physics timeline is in
  `template.html`; `build.py` embeds all eight frozen endpoints per mapped
  width. `#static` and `#norain` URL flags plus `window.__ox` exist for
  headless QA.

Owner asks to Freddie (`ASK_FREDDIE_DATA_REQUEST_2026-08-20.md`): one
cross-section, dose-to-CD map, actual Cr thickness, self-bias reading,
rinse/dry protocol.

## 3. Arun / Resona clock-gate — state

Everything for this case lives in `partner-private/arun_resona_clockgate_2026/`
and must never reach a rented box (see section 7).

Done:

- Geometry pinned exactly: Nangate/FreePDK45 CLKGATE_X1 Metal-1, enlarged 40x,
  98.8 x 62.8 um plan, 30 um tall, 2.6 um minimum gap, 0.19 nm max mismatch.
  Topology repaired (nine zero-area facets). `README.md` has the full chain.
- Deterministic entrance transport through the exact ten-polygon footprint
  (`src/petch/extruded_mask_transport.py`, tests in
  `tests/test_extruded_mask_transport.py`): ions ~92% direct transmission at
  1.5 deg sigma, thermal radicals ~7%. The mask is a radical filter.
- Conditional depth map via Belen SF6/O2 surface law, anchored to Miao's
  3.5–3.9 um open-feature depth: 2.80 um mean, 2.53–2.95 um 5–95% spread at
  the central scenario. `RESULTS_EXACT_MASK_ENTRANCE_TRANSPORT_2026-08-30.md`.
- Explorer `explorer/arun_etch_explorer.html` rebuilt as a clean recipient
  version on 2026-09-02: badges, warning states, and fine print removed;
  substance unchanged. Build with `python explorer/build.py`.
- Draft email `DRAFT_EMAIL_ARUN_2026-09-02.md` (not sent).

The honest gap, raised by the owner on 2026-09-02: the explorer is a
heightfield. It has no sidewall evolution, no polymer mask erosion, no
charging, no micro-trenching. A smooth map is the correct output of what was
computed, and it undersells what the engine can do. The owner's instinct is to
hold the email until a real profile exists.

Next concrete step (owner agreed in principle, not yet started):

1. Take one representative cut: a 2.8 um track through the 30 um polymer mask,
   silicon below, cryo SF6/O2 at the Miao fidelity recipe.
2. Run it through the common engine with moving sidewalls, polymer mask
   erosion at ~15:1 selectivity, and the charging tracer toggled on and off.
   The insulating 30 um mask with a 2.6 um opening is a plausible notching
   geometry, so the on/off difference is the thing to show.
3. Report the profile difference as a cross-section Arun can compare to an
   SEM, then add it to the explorer.
4. Constraints: no target SEM exists, so nothing can be fitted. The charging
   tracer is the CPU numpy path (`charging-speed-route` memory); at 30 M voxels
   the full 3-D domain is not practical, which is why this is a 2-D cut.

Owner asks to Arun (in the draft email): polarity and target depth, etcher
model and whether it holds -110 C, self-bias/platen voltage and He backside
pressure, pre/post polymer height plus a blanket Si depth on the same run, and
a pre-etch image plus post-etch cross-section held until the prediction is
frozen.

## 4. Showcase inventory

- `showcase/index.html` — internal reactor-to-feature showcase with movies
  (`showcase/media/*.mp4`), built by `build_data.py`, `build_page.py`,
  `render_movies.py`. Published earlier as a claude.ai artifact; the artifact
  watch has since dropped and does not matter.
- `showcase/oxford_*` — the three Freddie pages above.

## 5. Test status

Full suite at 08be74d: 2274 passed, 7 skipped, 21m44s. Commits since then touch
only `showcase/` and `partner-private/` HTML, templates, build scripts, and
markdown. Re-run before any engine change.

## 6. Untracked files needing your decision

These are yours or pre-date this session; I did not touch them:

- `src/petch/reactor_global/c4f6_ion_neutral.py`,
  `tests/test_c4f6_ion_neutral.py`,
  `research_sources/library/morris-viggiano-paulson-1994-c4f6-ion-neutral.md`
  — C4F6 ion-neutral chemistry work. Commit or drop.
- `results/curated/mouth_equilibrium_probe_dx/` — probe outputs.
- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`.
- `results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/w200_s14.000_ion_low_tail_0p0_db9361cf72bb8203.json`
  — a duplicate-hash trajectory the audit did not need; the frozen audit does
  not reference it.
- `scratch_ignore_calc.py`.

## 7. Rules that stay in force

- Partner content stays under `partner-private/`. `scripts/make_box_archive.sh`
  excludes `partner-private`, `partner-private/**`, `tests/test_arun_*`, and
  the PARTNER_*/RESONA_* rules, and fails closed on a word-boundary path scan.
  Never publish a partner-private page as a claude.ai artifact.
- Never fit to a target SEM. Never claim atomic accuracy or reactor validation
  without measurements. Frozen boards stay frozen; scoring is post-hoc only.
- No `git reset`, `git checkout --`, `git clean`; preserve untracked files;
  no amend or force-push on this branch.
- Rented boxes: pick many fast vCPUs, check the cgroup memory cap (70 GB vs
  183 GB bit us), destroy your own boxes when done, never touch the account's
  exp*-spectate instances.
- Headless browser QA uses gstack `browse` only; do not use Chrome MCP tools.
  `browse` cold-starts flaky on 3 MB inline pages: run `browse stop`, then keep
  per-command waits at or under ~1 s and use the page's `#static` / `__ox`
  hooks instead of long sleeps.

## 8. Physics queue, ranked

1. Arun representative cut through the common engine with sidewalls, mask
   erosion, and charging on/off (section 3). Unblocks the Arun send.
2. Cr surface-chemistry deck for CHF3/SF6/O2, literature-first, validated on
   published blanket Cr rates, then re-run the Oxford boards through it as a
   new versioned board. Do not re-score against the SEMs until frozen.
3. Reactive F/O wall return for extruded polymer masks (Arun's dominant
   unknown), as a deterministic diffuse-return operator next to
   `extruded_mask_transport.py`.
4. Tests for the four SEM measurement scripts.
5. Minimal C4F6 reactor on Benck (your untracked work in section 6 is the
   start).
