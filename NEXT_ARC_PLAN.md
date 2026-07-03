# Next arc: charging → measured notching (the paper arc)

Hand-off plan. State as of `273f612` (solo repo) / `e39d583` (shared docs). Everything below is
gated; no claim ships without its gate number. Working rules (non-negotiable, proven this project):
**literature-first when stuck; every claim gets a reproducer script + numeric gate vs published data;
honest fails documented, never hidden; no unpublished fudge knobs.** Repos: `plasma-etching-code` is
solo (push freely); `/Users/stanislavdelaurentiis/chip-etch` is shared with Craig — add ONLY docs/
files there, always `git fetch && git rebase origin/main` before push. PETCH_DEVICE=cpu suffices for
all of this (charging solver is small grids; no ViennaPS needed → no GPU box).

## Where charging stands (read RECONCILIATION.md "charging" sections first)

Two solver configs in `src/petch/charging2d.py::solve_trench_charging`, both reachable:
- **Closure config** (`poly_um=0, rf_bursts=False`): floor-flux gate **PASS 0.039**; wrong foot
  energetics. Sources the production `charging_floor_profile` table used by
  `flags.surface_charging="hg"` in the knudsen path.
- **Full-geometry config** (conductor + first-order RF bursts): poly-line potential 6.1→44.2 V
  **PASS**, foot-flux ~AR-independent **PASS**, rising foot peak appears (31.9→66.8 V), V_c(AR1)=11.1
  on-anchor — but floor-flux RMSE **0.060** (regressed past 0.05) and foot-ion energy decays deep
  (15.2→18.0 to AR2, then →10.3; HG wants 15→27.5 rising).

**One diagnosed root cause for everything open: the deep-AR floor over-charges (V_c 53.6 vs HG 33 V
at AR4) because the first-order burst model under-supplies electrons deep.** Gates:
`scripts/charging_gate.py` (8-pt HG floor-flux curve, digitized in-script; Matsui 300 eV asymptote)
and `scripts/notching_gate.py` (foot energy/flux, poly-line potential). Results land in
`charging_gate_result.npz` / `notching_gate_result.npz`; figures `viz/charging_hg.png`,
`viz/notching.png` regenerate from them.

## Workstream A — full RF-phase electron trajectories (the root-cause fix)

**What:** replace the first-order burst weighting with phase-resolved electron transport. Per
electron: sample RF phase φ; the instantaneous sheath drop V_s(φ)=V_dc+V_rf·sin(ωt) filters
(Boltzmann exp(−eV_s/Te)) AND sets arrival energy/angle at the mouth plane; trace through the
in-feature field as now. Ion energies are already phase-correct (bimodal). Key refs (full-text
notes in RECONCILIATION + memory): HG JVST B 15,70 §II (method — electrons penetrate in bursts at
sheath minima); JAP 82,566 (the gate curves); APL 71,1942 (V_dc Lieberman form). Keep the annealed
relaxation + tail-averaged segment potentials (they fixed the ratchet); keep in-plane sampling of
the published distributions (3-D→2-D projection bug already fixed — don't reintroduce).

**Gates (all must hold simultaneously in the full-geometry config):**
1. Floor-flux 8-pt RMSE ≤ 0.05 (restore the pass WITH conductor on)
2. V_c(AR4) → 33 V ± 40% (currently 53.6)
3. Foot-ion energy RISES 15→27.5 eV over AR 1→4 (±20%) — the notch energetics
4. Poly-line 6→39 V (±30%) and foot-flux constancy must NOT regress
5. Matsui 300 eV: floor stays open at AR4-6

If a gate can't be met, bisect: bursts-only vs conductor-only configs to attribute, document the
plateau precisely (which term, which AR range), and stop rather than tune.

## Workstream B — multi-material etch-stop → the measured notch-depth gate

**What:** minimal multi-material support in the evolving engine (`src/petch/threed.py`): a material
id on the grid (PR mask / poly-Si / buried-oxide etch-stop), per-material rate multiplier at the
surface (oxide ≈ 0, PR ≈ 0, poly = 1). Then the notch sim: HG/Nozawa line/space geometry (0.3 µm
poly on oxide, PR mask, AR set by PR height), etch to the oxide, then **100% overetch** with
`surface_charging="hg"` — the deflected-ion foot flux (already wired, uses the published E_defl
table until Workstream A supplies our own) digs the notch at the poly/oxide junction. Measure notch
depth (lateral penetration at the foot).

**Gates:**
1. Notch depth vs AR: HG's validated 0.08/0.09/0.10/0.12/0.165/0.185/0.215/0.23 µm at AR 1.0→4.0,
   100% overetch (±30% — first wiring tolerance, state it)
2. Nozawa width dependence (JJAP 34,2107): notch depth vs open-area width saturates at W ≈ 2–5 µm
   (qualitative: monotone rise then plateau)
3. Fujiwara (JJAP 34,2095): monotone notch growth with AR over 0.7–2.8
This is **measured wafer data** — passing #2/#3 even qualitatively is the "no open code does this"
claim. Figure: notch profile snapshots + depth-vs-AR/width against the measured points.

## Workstream C (optional, tracked as #42) — DDA flux-conservative march operator

Per-iteration flux accounting in `src/petch/dda.py` so each bounce redistributes exactly (1−s) of
arriving flux (kills the compounding per-hop leak: 0.98^50 ≈ 0.36 = the ~2× high-albedo deficit).
Gate: `scripts/dda_static_gate.py` RMSE ≤ 0.05 vs the measured ViennaPS static curve (radiosity+GMRES
currently passes at 0.043 and remains solver of record; DDA is quarantined in that regime).

## Order & scope

A → B (B consumes A's foot energies; B is the paper's headline). C independent, do last or skip.
After each workstream: commit with measured numbers, push both repos per the rules, update
RECONCILIATION.md + `docs/what-petch-does.html` limits box + the Experiments page, regenerate
figures from the npz results. If everything passes: the writeup is
"open, GPU, differentiable feature-scale etcher — held-out wafer ARDE prediction + literature-gated
charging + measured notching" — assemble as a new docs page, then it's the paper skeleton.
