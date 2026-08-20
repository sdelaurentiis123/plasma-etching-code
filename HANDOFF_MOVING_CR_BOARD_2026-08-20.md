# HANDOFF — moving-Cr blind board, 2026-08-20

> **SUPERSEDED:** Read `CODEX_TAKEOVER_MOVING_CR_2026-08-20.md` instead.
> The debug-file/process/cost claims below were disproved on live inspection,
> the numerical singularity was identified and fixed in v4, and v3 caches are
> no longer eligible for the corrected board.

For: Codex (or any resuming session). Read this before touching the board.
State is honest, including one operational failure at the end.

## Where the board stands

**55 of 56 trajectories are computed, certified-cacheable, and committed** to
`results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/` (this commit).
They are v3-revision caches; the campaign script reuses them byte-for-byte
(`job_spec` match), so any relaunch only computes what is missing.

**The one missing cell is `w320 / s14.000 / ion_low_tail_0p0`** — and it does
not merely "take long": the worker computing it was observed spinning with
100% CPU, **frozen RSS (byte-identical over 5 min), and 0% GPU over 5+ min of
sampling**. That is a non-converging CPU-only iteration inside
`advance_feature_step_3d`, not slow physics and not the GPU path. Every other
cell of the board completed, including all of w280 and the other three w320
scenarios, so the loop is specific to this cell's geometry trajectory.

Ruled out already:
- Not the marching-cubes degeneracy (fixed in `4b656fd`, regression-tested).
- Not the mask-exhaustion crash (declared terminal event since `2ad547d`).
- Not OOM (183 GB box, zero oom_kills, RSS frozen at 3.9 GB).
- Bounded loops inspected: the remap bracketing loop raises past 1e300; the
  adaptive-stepping loops have retry budgets; the neutral-radiosity
  `while True` at feature_step_3d.py:1860 is NOT exercised (ion-only
  scenario). The spin site is therefore unidentified — capture it, don't
  guess it (recipe below).

## This session's commits (all pushed)

- `2ad547d` — v3 mask-exhaustion guard: trajectories end as declared
  `cr_mask_below_vertical_resolution_in_footprint`; enriched step-failure
  context (elapsed, steps, mask metrics in the raise).
- `ac3b2e2` — make_box_archive.sh: word-boundary path scan ("resonance"
  files no longer refuse the build).
- `1a39d99` — failed-run logs + root-cause note preserved.
- `4b656fd` — extraction retry-with-nudge for grid-aligned float-noise
  fields; regression test carries the exact captured production field
  (`tests/data/w80_step80_grid_aligned_zero_phi.npz`).
- `f6003d0` — C4F6 resolved-current summation order pinned (hash-seed
  nondeterminism; exact-replay now stable across seeds).
- `5810075` — ASK_FREDDIE_DATA_REQUEST_2026-08-20.md (draft, unsent).
- `b955ff0` — pattern-collapse criterion research (9 library entries): intact
  pillars have ~8–11x stiffness margin against drying collapse; fallen
  pillars must be etch-weakened. Sharpens the blind narrative.
- This commit — the 55 trajectory caches + this handoff.

Full suite at `4b656fd`+: 2125 passed. Physics headline so far: the Cr cap
fails **corner-first at roughly one-third of the etch** on the tight cells
(w80 low/s14 tripped the guard at reference t≈391 s with centre still 27 nm).

## The box (LIVE, billing)

- Vast instance **48177892**, label `petch-zhu-moving-cr-bigram`,
  RTX 3090, 183 GB cgroup, 28 vCPUs, **$0.176/hr**, `ssh -p 17892
  root@ssh6.vast.ai`.
- Tree: `/root/petch-4b656fd` (clean archive at 4b656fd). Venv:
  `/root/petch-venv`. The venv has NO editable petch — always run with
  `PYTHONPATH=/root/petch-4b656fd/src`.
- The 55 caches also live on the box in the tree's trajectories/ dir.
- The old box 48118439 is destroyed. **Destroy 48177892 when done**
  (`yes | vastai destroy instance 48177892`).

## How to capture the spin (the part that was NOT completed)

An instrumented solo driver was written to `/root/run_w320_debug.py` on the
box (per-step prints + `faulthandler.register(SIGUSR1)`), but its launch
**silently failed and was misread as running** — two traps for the next
operator:

1. `pgrep -f run_w320_debug` **matches the ssh probe's own command line**;
   verify liveness with `ps -o pid,etime,pcpu -p <pid>` on a pid taken from
   a previous command, or pgrep with a pattern that the probe doesn't carry.
2. Launch nohup in its own parenthesized group and confirm the LOG EXISTS
   AND GROWS before believing anything:

```
ssh -p 17892 root@ssh6.vast.ai
cd /root && source petch-venv/bin/activate
PYTHONPATH=/root/petch-4b656fd/src nohup python run_w320_debug.py \
  > /root/w320_debug.log 2>&1 < /dev/null &
sleep 10 && wc -c /root/w320_debug.log && tail -2 /root/w320_debug.log
```

Then: watch the log's `step N elapsed ...` lines. When they stall for
minutes with the process at 100% CPU:

```
kill -USR1 <pid>       # faulthandler prints the exact stack into the log
tail -50 /root/w320_debug.log
```

That stack names the spinning loop. Fix it in the engine with the same
discipline as 4b656fd (never weaken a certifier; bounded iteration with a
raise, or a declared physical event), add the failing configuration as a
regression test, and re-run the cell.

Note: ptrace is blocked in the container (yama scope 1 + read-only), so
py-spy on a running process does NOT work — the faulthandler hook is the
only stack window. The run is deterministic: the spin will reproduce at the
same step every time.

## Landing sequence once the cell completes

1. `PYTHONPATH=/root/petch-4b656fd/src PETCH_PROFILE_WORKERS=4 python
   scripts/audit_zhu_npg80_moving_cr_profiles.py --write
   --transport-device cuda:0` — everything cached, it assembles and writes
   `audit.json` quickly. (8 workers fit this box; 12 do NOT fit a 70 GB box —
   per-step transients of many GB per worker, growing with width.)
2. Pull `audit.json` + the last cache; run `--check` locally; commit.
3. Render the atlas (see `render_zhu_npg80_conditional_profile_atlas.py`
   and the moving-mask audit's own `_render`).
4. Destroy the box.
5. Freeze the prediction package against Freddie's SEMs; the data request
   draft is `ASK_FREDDIE_DATA_REQUEST_2026-08-20.md` (owner sends).

## Open questions for the fixer

- Why does ONLY w320/s14/low/tail0 spin? w320 has the widest cap (interior
  margin spans ±150 nm) and s14 the fastest Cr erosion; low energy means the
  gentlest ion redistribution. Candidate: an oscillating sub-cell cleanup or
  a CFL/interface iteration that alternates between two states — but that is
  a GUESS; the stack capture decides.
- The two 12-second v2 smoke caches were deleted from the committed set
  (revision mismatch would have excluded them anyway).
