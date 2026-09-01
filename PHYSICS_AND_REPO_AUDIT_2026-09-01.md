# Physics and repository audit — 2026-09-01

Branch `codex/validation-first-multiphysics`, HEAD `ca54a1c`. Read-only audit;
every number below is quoted from the file named beside it. "Not found" means
exactly that.

---

## Part 1 — What the science is showing us

### 1.1 The thesis, and how far it has been earned

The repository's thesis is that a deterministic, differentiable, no-fitted-knob
chain — machine knobs → species-resolved plasma → wafer flux/IEAD → deterministic
feature transport → surface chemistry → moving 3-D profile — can predict etched
features. The evidence now splits cleanly into three classes.

**Established (validated against independent measurement, no fit to the compared
data):**

- *Surface physics at the beam level.* Tinacba SF5+ on Si/SiO2: "mean absolute
  depth-per-dose error is 5.88%; the largest error is 15.04%"
  (`results/curated/tinacba_2021_sf5_depth/README.md:12`). Vella–Hao Si-Cl ALE:
  predictions "0.4842/0.9453/1.4006 nm/cycle" vs "0.5558/1.0453/1.5427",
  "Maximum nominal relative error is 12.9%" (`vella_hao_ale_depth/README.md:18-20`).
  Humbird–Graves held-out MD: all gates pass, cumulative-etch normalized RMSE
  0.0325 against a 0.15 threshold (`humbird_graves_reduced_response/audit.json`).
  Karahashi SiO2 chemical channel 1.570 vs 1.5 molecules/ion at 1 keV (session
  record; see `RESULTS_RATE_GAP_CLOSURE_2026-08-06.md`). Krüger mask width
  850.2 vs 850 nm and aperture 41.8 vs 41.1 nm, "PASS — exact" / "PASS — 1.7 %"
  (`VALIDATION_DOSSIER_KRUEGER_2026-08-05.md:31-33`).
- *Transport.* de Boer SF6/O2 ARDE through the coupled engine predicts the
  held-out AR40 floor (0.166–0.169 vs 0.20; `AUTONOMOUS_PROGRESS.md:114-115,134`).
  The radical-only channel was **definitively refuted** as the explanation of the
  deep floor ("collapses to ~0.035 at AR40", `AUTONOMOUS_PROGRESS.md:105-109`);
  the floor is an ion-assisted second channel.
- *Reactor tier, piecewise.* BOLSIG+ swarm grade 0.122–0.634% residuals
  (`FULL_STACK_EXECUTION_PLAN_2026-08-08.md:74-76`); Lam chlorine ICP electron
  density transferred 300 W → 500 W within "+1.76%"
  (`FULL_STACK_EXECUTION_PLAN_2026-08-08.md:222`); moving RF sheath with global
  residual 3.394e-15, 453.165 V max sheath, 234.458 eV mean wafer ion energy,
  1.151° RMS angle (`MOVING_RF_SHEATH_REACTOR_BREAKTHROUGH_2026-08-13.md:60,68,70`).
- *One sealed machine-to-wafer depth win.* Bosch cyclic SF6/C4F8: "0.239 um
  wafer-mean depth MAE and 0.541% MAPE over 13 unseen wafers, each with 89 radial
  points" (`CODEX_TAKEOVER_REPORT_2026-08-21_1700UTC.md:50-51`) — wafer depth, not
  a feature profile.
- *Blind feature prediction — mechanism level.* Oxford/Freddie: the frozen board
  (`FROZEN_OXFORD_PREDICTION_2026-08-31.md`) said the Cr cap is consumed and pillar
  tops are attacked with spatial survive/fail zones, and that loss is not
  capillary collapse. The 80 SEMs show loss zones in doses 3/7/8/9 whose sites
  hold stubs, intact flat-topped pillars elsewhere, and zero fallen pillars
  (`showcase/oxford_scorecard/template.html`, scorecard rows 1–2, both HIT).
  Pitch 334–346 nm vs 350 nm GDS; CDs 181–236 nm inside the 105–250 nm ladder
  (`data/experimental/zhu_2026_tio2_npg80/sem_0817/dose_summary.json`).

**Conditional (correct physics, interval-valued inputs):**

- Oxford absolute depth/CD: the surface law is "cross-machine rate-normalized
  (Janissen witnesses); self-bias/IED are interval inputs"
  (`FROZEN_OXFORD_PREDICTION_2026-08-31.md`, claim boundary). The magnitude call
  "exhaustion in every cell" was over-predicted: 5 of 9 dose blocks survived.
  Tilt-geometry height 0.25–0.55 µm ±150 nm (`sem_0817/occlusion_height.json`);
  trench depth unscored.
- Arun clock-gate: exact-mask entrance transport gives ion transmission 0.9165 at
  1.5° and direct thermal-neutral transmission 0.06874
  (`partner-private/.../RESULTS_EXACT_MASK_ENTRANCE_TRANSPORT_2026-08-30.md:20,23`);
  the depth map (2.80 µm mean, 2.53–2.95 µm) is "a sensitivity example, not the
  frozen prediction" (same file, §Conditional surface-law transfer) because F/O
  wall return is unmeasured. Quadrature and grid sensitivities are 0.0125 and
  0.0146 absolute (lines 56-59).
- Reactor stacks: every audit still carries `supports_equipment_prediction:
  false` / `supports_feature_depth: false`
  (`results/curated/current_driven_argon_reactor_stack/audit.json:99-100`).

**Refuted or retracted:**

- The old 790–811 nm Krüger "match" — "it arose from canceling implementation
  errors"; corrected published-input endpoint "about 346.833 nm"
  (`CODEX_TAKEOVER_REPORT_2026-08-21_1700UTC.md:213-216`).
- The 13.75 nm/s blanket-anchor plan — "that number is a feature-average rate, not
  a published blanket measurement" (`PROGRESS_STATE_2026-08-19.md:179`).
- The LER slope-0.5 rule from static shadowing — "not reproducible from static
  shadowing, and the shortfall is mechanistic rather than numerical"
  (`LER_DEMONSTRATOR_RESULTS_2026-08-06.md:76-77`).
- Capillary collapse as the Oxford loss mechanism — 8–11× stiffness margin
  (`research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`) and zero
  fallen pillars in the SEMs.
- The 47× ViennaPS speedup (weak-vCPU artifact; honest ~14×, memory record).

### 1.2 What it jointly means

Three lessons recur across every arc, and together they are the scientific
result of the summer.

1. **The surface laws are right where they were measured, and depth misses come
   from the boundary, not the surface.** Every beam-level board passes at the
   5–15% level with nothing fitted. Krüger's 825 nm cannot be reached from his
   published aggregate ion flux without exceeding "the supported physical-yield
   ceiling" (`CODEX_TAKEOVER_REPORT…:216`); the missing quantity is a
   species-resolved wafer boundary. The repository's own board records "the count
   of completed formal held-out feature/profile predictions is **zero**"
   (`results/curated/cross_chemistry_depth_evidence/README.md:30-31`); Bosch has
   since added one at the wafer-depth level, still not a profile.
2. **Knife-edge processes are decided by inputs the literature cannot supply.**
   The Oxford mask-exhaustion prediction was mechanistically right and
   magnitude-wrong because a 45 nm Cr cap against a ~49 nm/20-min Cr loss is
   decided by ±10% in a selectivity we could only take from other tools. The
   PROGRESS_STATE enumerates nine TiO2/Cr coefficients that "cross-process plasma
   recipes … do not uniquely identify" (`PROGRESS_STATE_2026-08-19.md:139-151`).
   The SEMs converted that abstract list into one dominant lever: Cr chemistry.
3. **Deterministic transport is now a settled instrument.** Direct-characteristic
   quadrature (Arun), Clausing/cascade hole delivery (`HOLE_STUDY_RESULTS_PHASE2`),
   analytic 2-D exchange (`deterministic_exchange_2d.py`), and the moving RF
   sheath all close conservation to roundoff and are graded against independent
   references. What remains open in transport is *reactive* return (F/O on
   polymer walls, Arun) and sub-degree IADF/charging at extreme AR (de Boer
   floor; hole study §Honesty appendix line 163).

The one-sentence synthesis: **the engine predicts the right mechanisms blind, and
its remaining absolute-depth errors are traceable, one by one, to a named missing
boundary measurement or a named placeholder surface law — never to the
transport or numerics.**

---

## Part 2 — Missing physics inventory

| # | Physics item | Where it bites | Current treatment (file) | Real treatment needs | Impact | Effort |
|---|---|---|---|---|---|---|
| 1 | **Cr surface chemistry** (CrOx/CrFx conversion, ion-assisted removal) | Oxford/Freddie mask exhaustion — the single lever behind the over-prediction | Rate-normalized dimensional law from selectivity; router carries `"nguyen_2021_topology_warning": "real Cr removal in F/O plasma resolves CrOx, CrFx, neutral conversion, and ion-assisted inhibitor removal"` (`scripts/audit_zhu_npg80_moving_cr_profiles.py:199-202`) | Mixed-layer deck (Cr–F, Cr–O sticking; CrFx/CrOx sputter yields vs energy; ion-assisted removal); sourced constants **not found** in the library — Nguyen 2021 supplies rates, not a mechanism. Validate on published blanket Cr rates in CHF3/SF6/O2 before any board | Decides every survive/fail block; ±10% flips 5/9 blocks | High: literature hunt + new deck + witness validation |
| 2 | **TiO2 surface coefficients** (9 items) | Oxford absolute rate, CD bias, top loss | Fail-closed deck with "deliberately no default surface probabilities or yield curves" (`src/petch/tio2_surface_deck.py` docstring); yield/selectivity are "caller-supplied evidence axes" (`tio2_square_pillar.py`) | Site density, ALD density, fluorination probability, passivation sticking/density, O cleanup/blocking, species/energy-resolved yields, passivation sputter yields, roughness form (`PROGRESS_STATE_2026-08-19.md:141-149`) | Sets the depth interval width; today spans 0.82–1.63 units/ion | High; needs Freddie's blanket loss + self-bias to close most of it |
| 3 | **Diffuse/reactive F and O wall return** | Arun printed-polymer mask; every deep feature's neutral budget | Only the direct no-return component: "diffuse/reactive wall return remains a separate declared operator" (`src/petch/extruded_mask_transport.py` docstring). Radiosity solver exists without species/material law (`neutral_radiosity_3d.py`); 2-D analytic exchange is authority for line trenches only (`deterministic_exchange_2d.py`) | Deterministic source-iteration on the extruded mask with a sourced polymer-wall F/O loss law; conservation + convergence on analytic geometries | Arun depth interval is dominated by this (the amber slider) | Medium-high: operator exists in pieces; the wall law needs one blanket measurement |
| 4 | **Moving polymer/silicon sidewalls for Arun** | Arun absolute 3-D profile | Heightfield transfer only; "does not yet evolve full silicon/polymer lateral sidewalls" (session takeover brief) | Route the extruded-mask boundary into `MaterialMechanismRouter3D` + level-set with polymer erosion law | Taper, bowing, mask erosion unrepresented | Medium: engine supports it; needs polymer law |
| 5 | **C4F6 reactor closure** for Krüger | Krüger 825 nm boundary | Electron collisions "failing closed on a C4F6 reactor state" (`c4f6_electron_collisions.py`); ion sources are a "branching *prior*, not an energy-resolved branching measurement" (`c4f6_ion_sources.py`); ion–neutral topology has "deliberately *no numerical defaults*" and the Morris 1994 table is "NOT RECOVERED" (`c4f6_ion_neutral.py` [untracked]; `research_sources/library/morris-viggiano-paulson-1994-c4f6-ion-neutral.md` [untracked]) | Smallest atom/charge-conserving global model retaining C4F6, C3F3, CF, CF2, CF3 + ions, secondary ionization, ion-neutral conversion, residence, Bohm/wall loss; graded on Benck feed/pressure board (`CODEX_TAKEOVER…:249-258`) | Only honest route to the 825 nm question | High |
| 6 | **Sub-degree IADF / in-feature charging at high AR** | de Boer AR40 floor (0.166 vs 0.20), 200:1 hole study | Reduced ion channel with σ≈1° (`AUTONOMOUS_PROGRESS.md:116`); charging modules exist but not coupled to these boards; hole study lists "In-feature charging (bounds every deep-AR claim…)" as open (`HOLE_STUDY_RESULTS_PHASE2_2026-08-06.md:163`) | Couple `charging_nodal_*`/`charging_coupled_3d` to the ARDE boards; sheath-MC-free IADF source (W3 in `CHARGING_PHYSICS_PLAN.md`) | Residual 0.03 at AR40; dominant at AR>50 | Medium-high |
| 7 | **LER dynamic mask-erosion rung** | Partner LER modality | Static shadowing only; slope ≈1.0 vs published 0.5 (`LER_DEMONSTRATOR_RESULTS_2026-08-06.md:44,76`) | Resist-erosion foot sweep through the 3-D engine, preregistered | Any LER-reduction claim blocked until built | Medium |
| 8 | **Per-block litho/Cr-thickness variation** | Oxford: why doses 1/2/4/5/6 survived and 7/8/9 did not | Not representable: one periodic cell, reactor radial variation ≈0.6% rules out the smooth explanation (`PROGRESS_STATE_2026-08-19.md:49-56`) | Per-block inputs (dose→CD, lift-off Cr thickness) — a data problem, not a solver problem; optionally a small multi-pillar patch run | Explains the survive/fail pattern | Low once Freddie supplies dose map |
| 9 | **Pillar footprint shape** (4-lobed clover vs square) | Oxford CD bias, corner-first cap failure | Square footprint (`tio2_square_pillar.py`, GDS squares) | Use the measured top footprint (lattice-averaged cell) or litho-proximity model as the mask prior | Corner budget changes exhaustion timing | Low-medium |
| 10 | **Charging on insulating stacks** (TiO2 on fused silica) | Oxford foot notching at over-etch; loss-zone bases | Not in the Oxford boards; charging validated only on Hwang–Giapis notch demonstrations (`hwang_giapis_notch_validation_2d.py`) | Add the D1–D5 charging ladder to the pillar board | Unknown until a cross-section shows the bases | Medium |
| 11 | **Krüger species-resolved ion boundary** | Krüger 825 nm | "Krueger reports an aggregate positive-ion population rather than a species-resolved ion flux" (`CODEX_TAKEOVER…:238-239`); three no-fit sensitivity branches declare `parameter_evidence_supports_prediction=false` (line 233) | Authors' HPEM/PCMCM flux + IEAD, or item 5 | Blocks the flagship benchmark | Data request / item 5 |
| 12 | **Fluorocarbon passivation on TiO2 sidewalls; roughness evolution** | Oxford sidewall angle, CD | Passivation inventory in the reduced kernel; "chemistry-dependent roughness evolution" listed missing (`PROGRESS_STATE…:149`) | Sourced CHF3 polymer deposition/sputter law on TiO2 | Sidewall angle claim (85–86.5°) unscored | Medium |

Code-wide honesty markers (grep counts, `src/petch`): "refuse" 39 files,
"sensitivity" 34, "conditional" 14, `supports_feature_depth` 24, "fail(s) closed"
7. Docs: "refuse" 59, "sensitivity" 60. There is no `unquoted` marker in code;
it lives in `research_sources/library/*` rows.

---

## Part 3 — Repository state

`git status --short` at HEAD `ca54a1c` (local == remote after the last push):

| Path | Provenance | Disposition |
|---|---|---|
| `src/petch/reactor_global/c4f6_ion_neutral.py` | Codex in-flight, C4F6 ion–neutral topology module (Morris 1994), no numerical defaults | **Codex's decision** — coherent module, but not referenced by any tracked file (`grep` found no importer); needs its test to pass before commit |
| `tests/test_c4f6_ion_neutral.py` | Codex in-flight companion test | Commit together with the module once green |
| `research_sources/library/morris-viggiano-paulson-1994-c4f6-ion-neutral.md` | Codex library entry; status line says table "NOT RECOVERED" | Safe to commit with the module (library convention: extract + entry same commit) |
| `scratch_ignore_calc.py` | Collapse-criterion arithmetic scratch (Tanaka/Chandra formulas) — matches the pattern-collapse research | Scratch; do not commit (or fold into `research_sources/` if kept) |
| `results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/w200_s14.000_ion_low_tail_0p0_db9361cf72bb8203.json` | v7 revision, `duration_s = 12.0` — a smoke probe, not a board cell | Scratch; harmless, exclude |
| `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log` | Pre-August user artifact; explicitly "never edit, stage, remove" (`CODEX_TAKEOVER…:66-80`) | Leave |
| `results/curated/mouth_equilibrium_probe_dx/` | Pre-August probe (`probe.json`, two PNGs) | Leave |

`git stash list`: empty. Nothing unpushed.

**Test state.** `python -m pytest --collect-only -q` → **2281 tests collected**.
Latest recorded full-suite run: "2229 passed, 7 skipped" (`CODEX_TAKEOVER_REPORT_2026-08-21_1700UTC.md`); the 2125-pass run recorded at
`4b656fd` predates it. **A full suite has not been recorded at HEAD** — the 52
newer collected tests (SEM scripts have none; the GDS/STL/extruded-mask/explorer
tests do) were run focused only. Focused files for the newest modules:
`tests/test_extruded_mask_transport.py`, `tests/test_arun_partner_etch_explorer.py`,
`tests/test_mask_footprints.py`, `tests/test_gds_import.py`,
`tests/test_zhu_npg80_gds_square_profiles.py`,
`tests/test_summarize_zhu_npg80_gds_square_profiles.py`,
`tests/test_closeout_zhu_npg80_profile_caches.py`, `tests/test_audit_stl_geometry.py`.
The three SEM digitization scripts added 2026-09-01 have **no tests**.

**Partner-safety.** `scripts/make_box_archive.sh` now excludes
`partner-private/**` and `tests/test_arun_*` (commits `887b18a`, `87e5b3b`);
the last archive reported "0 sensitive paths, 0 content mentions".

---

## Part 4 — Ranked next steps

**(a) Recipient-driven**

1. **Freddie: one cross-section through an intact block.** Scores trench depth
   and sidewall angle (656–684 nm, 85–86.5°) — the two largest frozen numbers
   still unscored; removes the ±150 nm tilt-geometry slack.
   *Send:* `ASK_FREDDIE_DATA_REQUEST_2026-08-20.md` §6 + scorecard §4.
2. **Freddie: dose→CD map, per-block Cr thickness, DC self-bias, rinse/dry.**
   Converts the survive/fail pattern from "consistent" to "explained" and collapses
   the yield interval. *Send:* same draft, items 1–5.
3. **Arun: one blanket Si etch + DC self-bias on his tool, pre/post polymer
   height, one blind cross-section.** Pins the F/O wall-return slider; nothing
   else turns the conditional map into a prediction.
   *Send:* `partner-private/arun_resona_clockgate_2026/explorer/arun_etch_explorer.html`
   with its "Minimal experimental closure" table.

**(b) Physics to build**

4. **Cr surface deck** (item 1). Literature-first: Nguyen 2021, Cr fluorination
   XPS, CrF3/Cr2O3 sputter yields; validate on published blanket Cr rates; then a
   labeled post-SEM Oxford v2. *Start:* `research_sources/library/` search for
   `nguyen-2021-cr-sf6-o2`, then a new `src/petch/cr_surface_deck.py` shaped like
   `tio2_surface_deck.py`.
5. **Deterministic reactive wall return on the extruded mask** (item 3).
   *Start:* `src/petch/extruded_mask_transport.py` + `neutral_radiosity_3d.py`;
   gate on an analytically checkable slot before Arun.
6. **Minimal C4F6 global reactor** graded on Benck (item 5). *Start:* land
   `c4f6_ion_neutral.py` with its test; follow `CODEX_TAKEOVER…:249-258`.
7. **Charging + sub-degree IADF into the ARDE boards** (item 6).
   *Start:* `CHARGING_PHYSICS_PLAN.md` W1–W3 against `scripts/deboer_feature3d.py`.
8. **Measured-footprint mask prior for Oxford** (item 9): feed the
   lattice-averaged pillar top (`data/experimental/…/sem_0817/`) as the mask
   footprint in a labeled sensitivity run. *Start:*
   `scripts/audit_zhu_npg80_gds_square_profiles.py` geometry path.

**(c) Hygiene**

9. **Run the full suite at HEAD and record the count** (last recorded 2229/7
   skipped predates 15+ commits): `python -m pytest -q 2>&1 | tail -3`, then note
   it in the next state doc. Add tests for the three SEM scripts (metadata parse,
   lattice pitch on a synthetic array, occlusion solve on a synthetic profile).
10. **Resolve the untracked C4F6 files with Codex**, commit or drop
    `scratch_ignore_calc.py`, delete the v7 12-second smoke cache, and update
    `FROZEN_OXFORD_PREDICTION_2026-08-31.md` with a one-line pointer to the
    scorecard (`showcase/oxford_scorecard/`) so the freeze → reveal → score chain
    is one hop from the freeze doc.
