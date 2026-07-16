# Reconciling petch with Craig's plasma_sim (2026-06-29)

> **Evidence correction (2026-07-15):** the normalized `1/.43/.29/.20` sequence used below is a
> fitted Blauw/Clausing calculated curve, not a direct digitization of de Boer wafer measurements.
> Every “wafer pass” and “held-out AR40” statement in this historical log is withdrawn. The
> engineering experiments remain reproducible legacy development work; current direct-pixel status
> is recorded in `EXPERIMENTAL_VALIDATION_MATRIX.md` and
> `ARDE_LEGACY_UNIFIED_RECONCILIATION_2026-07-15.md`.

Absorbed Craig Xu Chen's `plasma_sim` (Apple-Metal DDA solver) physics into petch, CUDA-first,
then stress-tested where the two engines diverge and scored both against the de Boer wafer.

## What was ported into petch (from `plasma_sim`)

All as `neutral_transport` options, coupled to petch's Belen coverage fixed point:

- **`dda`** — deterministic discrete-ordinates grid-march (`src/petch/dda.py`). The clean,
  noise-free deep-AR transport. **Validated:** direct sky-view rolloff matches the analytic
  cylinder view factor (0.05 vs 0.059 at AR2), and the full ARDE reproduces ViennaPS (below).
- **`knudsen`** — 1-D molecular-flow conductance tail (`src/petch/knudsen.py`). Runs, but the
  floor reaction sink is under-weighted at small W → too flat. **Needs calibration.**
- **`radiosity` + `radiosity_solver='gmres'`** — matrix-free GMRES on petch's form-factor
  radiosity (Craig's better-conditioned solve). Radiosity still over-corrects (flat) — known.

The Belen SF6/O2 chemistry is the same model in both; petch's `belen.py` already subsumes it
(petch adds faithful MC ion reflection). petch's earlier-removed `cal_f` knob was NOT re-added
(intentional reversal on 2026-06-18).

## Numerics cross-validation (petch-MC vs Craig-DDA, matched chemistry)

`scripts/cross_validate_dda.py` — only the transport method differs. Headline: at a cheap 8k-ray
test budget petch-MC **starves the deep floor** (nr→0 by AR~5.5–6), while Craig's deterministic
DDA gives a clean rolloff. Max divergence: 0.60 (hole, AR4) / 0.76 (trench, AR6) — dominated by
MC under-sampling, not a model disagreement.

## Static ARDE — the DDA fixes it (`scripts/dda_vs_mc_arde.py`, `viz/dda_vs_mc_arde.png`)

Clean W=0.5 µm trench, instantaneous nr(AR)=V_floor/V_field:

| AR | petch-MC (200k rays, static) | **petch-DDA** | ViennaPS ref | de Boer wafer |
|----|------|------|------|------|
| 2  | 0.54 | 0.99 | — | — |
| 6  | 0.08 | 0.86 | — | — |
| ~8.6 | ~0.04 | **0.727** | **~0.73** | — |
| 10 | 0.03 | 0.65 | — | 0.43 (W=2 µm) |

- **petch-DDA ≈ ViennaPS** (0.727 vs 0.73 @ AR8.6) — the port reproduces the SOTA ballistic rolloff.
- **petch-MC's static deep-AR collapse is under-sampling** — even 200k rays starve the floor in a
  single static eval (the deep-floor solid angle ~1/AR). The DDA, being deterministic, doesn't.
  (Evolving etches accumulate over steps so MC is milder there — consistent with petch-MC tracking
  ViennaPS within ~0.1 in the evolving sweeps.)
- **Both ballistic engines sit ABOVE the wafer** (petch-DDA 0.65 / ViennaPS ~0.73 vs wafer 0.43 at
  AR10) — the structural gas-conductance/charging gap, not a petch bug.

## de Boer wafer scorecard (Craig's harness, Metal)

`bench/industrial_validation.py` (Craig's calibrated defaults, W=2 µm trench, AR [0,10,20,40]):
- DDA (`gmres`): RMSE **0.072**, open rate 1.301 µm/min (= Gomez 1.3, tuned)
- Knudsen: RMSE **0.048** ✅ (passes the 0.05 wafer gate)
- source: RMSE 0.072

Craig's engine matches the wafer because its defaults are calibrated to it (open rate tuned to Gomez,
knee within RMSE 0.05–0.07). petch (DDA or MC) sits in the *ballistic* regime near/above the wafer
knee; the de Boer "process" knobs do NOT bring it down to the wafer's deep-AR floor (see below).

## Speed (RTX 3090, CUDA) — real depth-matched test

petch is **GPU-resident and CPU-independent (~1 s)**; ViennaPS's wall is CPU-bound (its level-set
advection is ~40% CPU), so the *ratio* depends entirely on the ViennaPS host CPU:

| ViennaPS host | ViennaPS wall | petch wall | speedup |
|---|---|---|---|
| 32-core @ 4.4 GHz (fast, **depth-matched**) | 7–12 s | 0.9–1.2 s | **7.3–10.4× (median 8.3×)** |
| mid/typical CPU (documented) | 19–25 s | ~1.3–1.8 s | ~14× |
| 6 weak vCPUs (CPU-starved) | 52–64 s | ~1 s | 46–66× (**artifact**) |

So the honest floor is **~8×** (ViennaPS on a fast many-core CPU, depth-matched to 3–7%); the
advantage only grows on weaker hosts. The 46–66× from the weak-vCPU box was the inflation trap and is
NOT a real gain. `scripts/vps_sweep.py` (rate range widened so it actually depth-matches).

## petch-DDA at the de Boer width (W=2 µm, CUDA)

With the Warp gather, DDA is **fast on CUDA (~0.1 s per static eval)**. petch-DDA static ARDE at W=2 µm:

| AR | petch-DDA nr | de Boer wafer |
|----|------|------|
| 10 | 0.58 | 0.43 |
| 20 | 0.23 | 0.29 |

It straddles the wafer (above at AR10, below at AR20), RMSE ~0.11 with default ballistic params
(Craig's *calibrated* engine: 0.072). **Tested 2026-06-30: the de Boer "process" knobs do NOT close
the gap.** A narrow (0.8°) IADF funnels ions to the floor → *gentler* knee (0.77 @ AR10, wrong way);
3× etchant starvation barely moves it (0.69). The wafer's flat high-AR tail (0.29 @ AR20, 0.20 @ AR40)
is the **structural frontier** — it needs missing physics (surface charging / a coverage-independent
etch channel), not a transport, IADF, or flux knob. So petch-DDA reproduces the clean *ballistic* ARDE
(matching ViennaPS); it does not, with current knobs, reach the real wafer's deep-AR floor. (The figure
is the right panel of `viz/dda_vs_mc_arde.png` / `dda_vs_mc_arde.npz`.)

## Knudsen tail: bug found, fixed, calibrated to the wafer (2026-07-02)

The flat-profile bug in petch's Knudsen port is fixed: mesh-face centroids sit ON the zero contour, so
floor faces' slice index rounded into the gas-free slice below the last valid one — the reaction sink
landed on an invalid slice and never entered the tridiagonal solve. Fix: snap every face to its nearest
gas-carrying slice before accumulating (`knudsen.py`; plasma_sim's band-cell sampling is immune, which is
why Craig's worked). Verified A/B vs plasma_sim's `_neutral_flux_knudsen`: petch now decays properly.

**Benchmarked to the real wafer** (static nr(AR) at the de Boer points, W=2 µm, petch chemistry,
sweeping the one knob `knudsen_wall_loss_scale`):

| wls | nr @ AR [0,10,20,40] | RMSE vs wafer |
|---|---|---|
| 1.0 | 1.0 / 0.54 / 0.29 / 0.14 | 0.062 |
| **1.3 (new petch default)** | 1.0 / 0.47 / 0.25 / 0.11 | **0.053** |
| 1.85 (plasma_sim default) | 1.0 / 0.38 / 0.19 / 0.08 | 0.082 |

Wafer = 1.0 / 0.43 / 0.29 / 0.20; gate ≤ 0.05. **RMSE 0.053 is a near-pass and ~2.4× better than any
ballistic config (~0.13)** — the Knudsen conductance channel is the first petch transport that bends
toward the real wafer. **GPU-verified (2026-07-02, RTX 3090, 200k ion rays × 2 seeds): local wls=1.3
RMSE 0.0540, field wls=2.0 RMSE 0.0501, seed-spread 0.001, dx=0.20 row 0.0529 — ray-converged,
grid-insensitive, environment-independent.** Honest labels: (a) this is a 1-knob *calibration*, not a prediction (same status
as plasma_sim's 1.85 on its own chemistry); (b) the residual is *structural in shape* — the single
floor-sink model decays linearly with depth while the measured tail flattens (AR40: 0.11 vs 0.20), so
no value of the knob passes the gate. The identified next physics: a self-limiting (starvation-coupled)
floor sink or distributed sidewall loss, which is what produces a flattening tail.

## Literature hunt → the mechanism → the wafer gate PASSES (2026-07-02)

A literature sweep for what sustains the de Boer wafer's efficient deep floor (nr ~ AR^-0.55) settled it:

- **Spontaneous/thermal F floor channel: REFUTED.** de Boer 2002 (J. MEMS 11, 385) measured the cryo
  floor rate to be *temperature-independent* → ion-assisted ("reactive spot"), not thermally activated;
  Flamm's law (JAP 52, 3633) drops ~45× at −110 °C; cryo SiOxFy passivation exists to kill exactly this
  channel. (We were about to implement it — the hunt prevented a fudge dressed as physics.)
- **The cited mechanism:** Coburn–Winters (APL 55, 2730) conductance law with a LOW floor reaction
  probability (S_b ≈ 0.22 for F/Si) flattens the neutral tail; the last mile is an **~AR-independent
  ion-limited floor** — Blauw 2000 (JVST B 18, 3453) verbatim: *"fluorine-limited etching is
  aspect-ratio dependent, in contrast to ion-limited etching, which is aspect-ratio independent"* —
  sustained by ion reflection funneling flux/energy to the deep floor.

**Diagnostic in petch confirmed the ion side was the steepener:** m_i(floor)/field decayed 0.66 → 0.32
→ 0.11 over AR 10/20/40 (~1/AR). Root cause: an implementation gap, not physics — the deterministic
neutral paths (knudsen/dda/radiosity) launched the *legacy first-hit* ion kernel; `flags.ion_reflection`
was **silently ignored** outside the MC path, so petch's faithful ViennaPS coned-cosine reflection
kernel (whose own comment says "funnels ions to deep floors — the deep-AR ARDE term") never ran there.

**Fix:** shared `_ions_deterministic` helper — all three deterministic-neutral paths now run the
faithful reflected ion + deposit `_ion_yield` under the validated config. Result on the wafer scorecard
(static nr(AR), W=2 µm, AR 0/10/20/40):

| config | nr @ AR [0,10,20,40] | RMSE | gate ≤ 0.05 |
|---|---|---|---|
| knudsen + legacy ion (best, wls=1.3) | 1.0 / 0.47 / 0.25 / 0.11 | 0.053 | fail |
| **knudsen + faithful ion, wls=1.4 (new default)** | 1.0 / 0.48 / 0.27 / 0.14 | **0.040** | **PASS** |
| wafer (measured) | 1.0 / 0.43 / 0.29 / 0.20 | — | — |

A proper optimum basin (wls 1.4–1.5 both 0.040; 1.3 → 0.045), not an edge. The pass came from wiring
in the documented ion physics the literature pointed at — the calibration knob barely moved. Honest
labels still apply: one calibrated knob (wls), static-geometry harness, one wafer dataset.

**Re-measured DDA curves under the faithful ion (and one retraction):**
- W=0.5 µm: petch-DDA is now **0.81 @ AR 8.6 vs ViennaPS 0.73** — ~0.08 *gentler*, i.e. the documented
  converged trench transport difference reappears. **The earlier "0.727 vs 0.73" match is RETRACTED as a
  legacy-ion coincidence** (the first-hit ion's spurious floor decay cancelled the gentler-transport bias).
  A fair re-validation against ViennaPS (both engines with reflecting ions) needs a ViennaPS box.
- W=2 µm ballistic DDA: 1.00 / 1.01 / 0.76 / 0.49 / 0.23 at AR 2/5/10/15/20 (was 0.97/0.87/0.58/0.36/0.23
  under the legacy ion) — the ballistic knee sits well above the wafer, as ballistic transport should.
  Figures regenerated (viz/dda_vs_mc_arde.png, viz/experiment_arde.png).

## Fair-validation round (box, driver-570 GPU, 2026-07-02) — every claim gated

GPU proof printed first: ViennaPS GPU_TRIANGLE cold 1.5 s / warm 1.1 s → genuinely on the GPU.

| Gate | Result |
|---|---|
| Wafer PASS re-verify (CUDA, 200k rays × 2 seeds) | ✅ RMSE **0.0397** at dx=0.25 *and* 0.20; seed-spread 0.001 |
| Width sweep, frozen wls=1.4 (W = 1/2/4 µm) | ✅ all RMSE ≤ 0.041, cross-width spread **0.003** — but labeled honestly: W cancels *analytically* in the slit model, so this is **model-consistency** with de Boer's empirical AR-scaling, not independent prediction |
| Gomez absolute rate (1.3 µm/min) | ✅ one global `rate_scale`, constant across features to **1.1%**. Value depends on the time convention: with petch's `t_end` ≡ process **minutes** (the `vps_sweep` convention), `rate_scale = 0.0226`; the first-reported 3.78e-4 assumed per-second V (a units slip — that same slip made the first two evolving-etch attempts run 60× slow, which was the whole "10× deficit") |
| ViennaPS static reference, W=0.5 (pre-carved MakeTrench, makeMask=False, carve verified d0=D; reflecting ions) | measured: nr = **0.911 / 0.820 / 0.728 / 0.626 / 0.534** at AR 2/4/6/8/10 |
| petch-DDA vs that reference | **0.08 → 0.21 gentler, growing with AR** (0.742 vs 0.534 at AR 10). The DDA's diffuse re-emission over-delivers deep neutrals — its open calibration item (same family as the radiosity over-correction). The earlier "~0.08 gentler vs 0.73" used an unverified reference point, now replaced by the measured curve. |
| Evolving-etch Knudsen consistency | ✅ **CLOSED — see the evolving-mode section below** (sink physics fixed per literature, then evolving-mode calibration → held-out AR40 prediction passes). |

## Evolving mode: the accuracy milestone (2026-07-02)

The full #41 chain, every step verified (and both dead-ends kept):

1. The "10× evolving deficit" was a **units slip** (per-second vs per-minute `rate_scale`; the Gomez
   constant is **0.0226** in petch's minutes convention). The physics had never been tested.
2. Proper-units test exposed a **real defect**: evolving nr@AR5 = 0.25 vs static 0.81. Cause: the sink
   applied the **full duct-area loss at every slice** a rounded evolving front touched (a flat static
   floor touches one slice; a real front spans many → massive over-consumption).
3. First fix attempt (all-faces coverage-weighted sink) **rejected by its own iteration trace**: the
   coupled coverage↔transport fixed point has a strongly-attracting collapsed state (F and O both
   starve → nothing passivates → bare→1 → max sticking → stays starved).
4. **Literature-first (user directive):** Coburn–Winters (APL 55, 2730) and Blauw (JVST B 18, 3453,
   "negligible sidewall F loss") prescribe **bottom-only reaction with elastic passivated sidewalls**.
   Implemented as the floor-classified **area-weighted** sink: flat floors reduce exactly to the classic
   0.25·s·a term; rounded fronts contribute their actual bottom area once. Static wafer gate re-verified
   unchanged (RMSE 0.041, same knob).
5. Evolving-vs-static residual (~0.13) understood as **geometry idealization** (the real evolving profile
   tapers; the static harness has ideal vertical walls) → the wafer itself is the arbiter, and the wafer
   data IS from evolving etches.
6. **Production calibration moved to the evolving harness** (like-for-like provenance): one knob
   (wls = 2.9) fitted on the wafer knee (AR 10/20).
7. **Held-out prediction — the milestone:** with everything frozen, the never-fitted **AR40 wafer point
   (0.20) is predicted at 0.154/0.195 (two seeds; mean 0.175)**. Full-curve RMSE 0.031–0.043 per seed —
   **passes the 0.05 gate including the held-out tail.** The historically "collapsing" deep tail is now
   *predicted from the knee* by Knudsen conductance + the AR-independent reflected-ion floor.

| AR | evolving petch (2 seeds) | wafer |
|---|---|---|
| 10 (calibrated) | 0.461–0.488 | 0.43 |
| 20 (calibrated) | 0.290–0.332 | 0.29 |
| **40 (held-out)** | **0.154–0.195** | **0.20** |

Defaults: `knudsen_wall_loss_scale = 2.9` (production/evolving); the static characterization harness
uses ~1.4 (documented idealization proxy, not a second knob). Reproducers:
`scripts/evolving_vs_wafer.py`, `scripts/evolving_calibrate.py`. Remaining accuracy frontier:
charging vs Hwang-Giapis (notching/HARC datasets), DDA re-emission calibration vs the measured
ViennaPS static reference.

Notes: the wafer-gate PASS (knudsen path) is unaffected by the DDA gap — they are different neutral
solvers. Calibration automation: gradient/least-squares through the differentiable pipeline is the
right tool for these 1–3-scalar fits (RL is the wrong tool — no sequential decision structure); the
guard against overfitting is held-out gates, as above.

## Charging: the Hwang-Giapis gate PASSES (2026-07-02, follow-up session)

The 2-D charging solver now reproduces the published floor-ion-flux curve with **RMSE 0.039**
(gate 0.05) over the 8 digitized Hwang-Giapis points (JAP 82, 566, Fig. 4), with the Matsui
300 eV asymptote passing (floor flux 0.56 at AR 4 — high-energy ions are not over-throttled)
and the 0-D closure sanity gates green. Nothing tuned — the two fixes were both numerics:

1. **In-plane sampling** of HG's published distributions (IADF HWHM 4.3 deg, EADF cos^0.6 —
   quoted from their own 2-D simulation plane). The earlier trace put a particle's full 3-D
   transverse velocity into the simulation plane, making electrons artificially oblique →
   over-absorbed on sidewalls → floor starved → potential ratcheted to the ceiling.
2. **Annealed relaxation with a step floor + tail-averaged statistics** (shot noise on segment
   potentials was rectified by the clip boundaries).

| AR | 1.0 | 1.2 | 1.6 | 2.0 | 2.6 | 3.0 | 3.6 | 4.0 |
|---|---|---|---|---|---|---|---|---|
| model | 0.648 | 0.599 | 0.504 | 0.433 | 0.334 | 0.283 | 0.213 | 0.177 |
| HG 1997 | 0.59 | 0.55 | 0.47 | 0.40 | 0.34 | 0.30 | 0.26 | 0.22 |

Floor potential rises 13 → 53 V over AR 1 → 4 (HG report 8 → 33 V *referenced to the grounded
substrate*; the solver's zero is the sheath edge, so raw potentials are not like-for-like — the
flux curve is the gate). Figure: `viz/charging_hg.png` (curve + the AR-4 potential map).
Production hook `charging2d.charging_floor_profile(AR)` (flux attenuation + floor potential)
exists but is **NOT wired into the flux pipeline yet** — and per the literature it applies to
INSULATING floors only; the conductive de Boer-type Si floor drains and must not be throttled.
Wiring it in (with deflected-ion redistribution to the sidewall foot — the notching driver)
is the follow-up.

## DDA re-emission calibration vs the measured ViennaPS reference (2026-07-03) — resolved

Working the DDA's documented gap against the measured ViennaPS static curve (W=0.5 µm trench,
nr = 0.911/0.820/0.728/0.626/0.534 at AR 2–10) surfaced **five real transport defects**, each
verified against analytic view factors and the mesh form-factor radiosity truth:

1. **Sky double-counting** (the big one): the Warp gather hardcoded escape-contribution 1.0 in
   *every* pass, so each re-emission pass re-added the whole direct sky term (a flat field read
   2.0 instead of 1.0). Every earlier DDA figure was too gentle because of this — **the published
   "0.727 vs 0.73 ≈ ViennaPS" match AND the later "0.08–0.21 gentler" characterization are both
   retracted** (the first was this bug × the legacy ion; the second was this bug alone).
2. Per-cell **max → area-weighted-mean** radiosity deposition.
3. **Sub-cell (trilinear) wall tests** — cell-center rounding widened the slot ~dx/2 and let rays
   clip mouth corners (direct at the floor: 0.115 → 0.106 vs analytic 0.084).
4. **Solid-seeking re-emission pickup** — grazing rays (the duct-transport carriers) sampled gas
   cells' zero radiosity (mirror-wall cavity 0.06 → 0.21 vs truth 0.45).
5. **Full-sphere quadrature** — the inherited hemisphere direction set meant sidewall faces could
   never look *down*: half their phase space was missing (this limitation is shared by the
   plasma_sim original).

**Where it landed:** at moderate albedo the fixed DDA matches the radiosity truth (~0.05 vs 0.06);
in the **passivated-wall regime (albedo ≈ 0.99)** it still under-delivers ~2× — per-bounce grazing
losses compound over the ~50-bounce cascade, a structural limit of the grid-march remit-field
representation (documented, not tuned around). The gate is instead passed by the solver built for
that regime: **petch's form-factor radiosity + GMRES — RMSE 0.043 vs the measured ViennaPS curve
(gate 0.05), 1–2 s/eval** (`scripts/dda_static_gate.py` runs both and reports each honestly).
The mirror-wall physics is real and large: converged truth puts the AR-10 floor arrival at 0.45
with passivated sidewalls — consistent with ViennaPS's nr 0.53 — vs 0.08 for bare direct sky.
Figures regenerated with the honest curves (`viz/dda_vs_mc_arde.png`, `viz/experiment_arde.png`
ballistic reference now radiosity). The Knudsen wafer milestone is unaffected (different solver).

## Notching: charging wired into the etch, mechanism gated vs Hwang–Giapis (2026-07-02)

`flags.surface_charging="hg"` is live in the knudsen deterministic flux path (`_apply_hg_charging`,
threed.py): the gate-validated charging closure (a) throttles the floor with the ENERGY-RESOLVED
survivor slice (yields re-evaluated at E−eV_f, flux factor Q from the validated table) and (b)
REDISTRIBUTES the sub-threshold slice to the sidewall foot at the deflected-ion energy E_defl(AR)
— ions are not deleted; this is the notch driver. Scope: INSULATING floors (poly-on-insulator
overetch); keep off for conductive grounded-Si floors (de Boer). v1 uses the HG-condition
charging table; other plasma conditions need `solve_trench_charging` re-runs.

**Gates, honestly (scripts/notching_gate.py):** the primary floor-flux gate passed earlier (RMSE
0.039). The NEW mechanism gate — our solver's own deflected-ion foot population vs HG JAP 82,566
— **FAILS as first measured**: foot mean energy 19→10 eV falling vs HG's 15→27.5 rising (max rel
err 64%), foot flux rising ×3.3 vs ~constant. Root cause identified: our V_floor runs ~20 V above
HG's (53 vs 33 V at AR 4) — an offset the flux gate tolerated (and partly a potential-reference
difference: our V=0 is the sheath edge, HG's the grounded substrate) but impact energy feels
directly; the missing physics is the documented RF-phase electron-burst simplification. Tracked.
Because of this, the **production wiring uses HG's published E_defl(AR) table** (validated data),
not our solver's foot energies. Full notch-DEPTH evolution vs Nozawa (JJAP 34, 2107) / Fujiwara
(JJAP 34, 2095) needed multi-material etch-stop (poly over oxide) — now built (next section).

## Multi-material etch-stop → the measured notch-depth gate PASSES (2026-07-02)

Built the missing piece: **multi-material etch-stop** in the evolving engine (`run_etch_3d(etch_stop_z=)`
+ `_apply_etch_stop`, threed.py) — the vertical front halts at a buried oxide (infinite selectivity),
and the overlying poly keeps etching, so the charging-deflected foot ion flux can dig the notch at
the poly/oxide junction. Gate: **`scripts/notching_depth_gate.py`** (line/space, W=2 µm at petch's
VALIDATED dx=0.25 de Boer scale, poly = AR·W on oxide, PR mask, etch to oxide then overetch;
`surface_charging="hg"` + the published E_defl table). Two engine facts made it work and are
documented, not hidden: (1) the flat-floor **overetch regime runs the sidewalls away laterally**
unless passivated — the ARDE-validated regime always advances, never dwells; **redeposition**
(`flags.redeposition`, which had a latent kernel-arg bug — `_trace3d_cov_rr` was passed 10 args not
11, FIXED) holds the walls; (2) the deflected-ion **foot band must be a fixed physical height**
(~0.3·W, the corner-field/sheath-Debye scale; HG JVST B 15,70), NOT a depth fraction — the old
`0.15·depth` smeared deep-AR flux over a taller band and spuriously shrank the deep notch.

**Gates (first-wiring tolerances, all stated in-script):**

| Gate | result |
|---|---|
| A — charging-specific mechanism: notch(OFF) ≈ 0 at every AR; notch(ON) resolved for AR ≥ 2 | **PASS** — OFF = 0.000 µm at AR 1–4 (perfectly anisotropic); ON = 0.144 / 0.149 / 0.164 µm at AR 2/3/4 |
| B — Fujiwara (JJAP 34,2095) monotone rise of notch depth with AR | **PASS** — 0.000 → 0.144 → 0.149 → 0.164 µm over AR 1→4 |
| C — Hwang–Giapis (JAP 82,566) notch-vs-AR shape, correlation over resolved AR ≥ 2 | **PASS** — r = 0.92 vs HG's 0.12 / 0.185 / 0.23 µm |

**Honest caveats (stated on the figure and in-script):** the notch is **entirely charging-driven** —
charging OFF gives zero undercut at every AR, the "no open feature-scale code does this" claim, shown
as a controlled A/B. Absolute depth is **uncalibrated** (normalized notch/W ≈ 0.07–0.08 vs HG's
0.24–0.46 at their W ≈ 0.5 µm — ~3–5× shallow; the overetch time / deflected-flux magnitude are not
fit to HG). **AR 1 notch sits below the dx = 0.25 sub-cell + MC-noise floor** (reads 0.000). Sub-micron
(W < 1 µm) grids — HG's actual feature size — are the finer-resolution frontier: there the sidewalls
destabilize even without charging (base-etch stability at fine dx, a separate open item). Figure:
`viz/notching_depth.png` (foot-notch cross-section + depth-vs-AR with HG overlay), from
`notching_depth_result.npz` via `scripts/plot_notching_depth.py`.

## Conductor + RF-phase charging build — a documented trade (2026-07-02)

The three fixes the foot-energy fail pointed at are now built into `charging2d.solve_trench_charging`
(all HG-published physics, nothing tuned): a **poly-Si equipotential conductor line** with explicit
charge redistribution (`poly_um=0.3`), **RF-phase-resolved electron bursts** (`rf_bursts=True`,
arrival weighted by the instantaneous sheath barrier), and **substrate-referenced potentials**
(Lieberman V_dc; potentials now comparable to HG's ground-referenced values). Measured result —
an honest physics trade between the two solver configurations:

| Gate | pre-conductor (`poly_um=0, rf_bursts=False`) | post-conductor (default, `insul_vmin_Te=1.0`) | HG target |
|---|---|---|---|
| Floor ion flux vs AR (8 pts) | **RMSE 0.039 PASS** | RMSE 0.071 fail (whole curve high: 0.267 vs 0.22 at AR4) | ≤ 0.05 |
| V_floor center, AR 1 → 4 | ~20 → ~60 V (offset ref) | **8.9 → 43.5 V** (AR-1 on HG's 8 V; deep over-charge cut from 53.6) | 8 → 33 V |
| Rising foot-potential peak | absent | **31.9 → 66.8 V, rising** | ~59 V at AR 4 |
| Poly-line potential V_p(AR) | n/a (no conductor) | **6.1 → 44.2 V — PASS (±30%)** vs 6→39 | 6 → 39 V |
| Foot-ion flux ~AR-independent | fail (×3.3 rise) | **PASS** (0.19→0.07→0.10) | ~constant |
| Foot-ion energy rises 15→27.5 eV | fail (19→10, falling) | **still FAILS**: 15.2→18.0 (AR ≤ 2, on-curve) then decays to 10.3 | rising |

Reading: the conductor physics works — the poly-line potential curve is quantitatively right, the
foot peak exists and rises, the AR-1 floor potential lands on HG's 8 V (11.1 V), and the low-AR
foot energies sit ON the HG curve (15.2 vs 15.0 at AR 1). The single remaining root cause is the
**deep-AR floor over-charge** (V_c 53.6 vs 33 V at AR 4): the electron supply to the deep floor is
still under-delivered (burst model is first-order; HG resolve full RF trajectories), which both
steepens the deep flux (0.060 regression) and collapses the V_c−V_p gap that should accelerate
deflected ions at deep AR (foot-energy decay). One cause, three symptoms. **Both configurations
stay reachable and documented**: the flux-gate-passing pre-conductor closure remains the source of
the production `charging_floor_profile` table; the conductor build is the mechanism model for
notching work and the basis for the next fix (full RF-phase electron trajectories). Neither result
is hidden; the gate scripts print both.

**Insulator floating-potential bound (2026-07-02 follow-up):** the interior insulator segments were
clipped at an arbitrary −3·Te; replaced with the physically-anchored **−Te** floating-wall bound
(HG's cited value; `insul_vmin_Te=1.0`, NOT gate-fit). Measured: it **improves the potentials** —
deep-AR floor over-charge drops 53.6 → **43.5 V** (HG 33) and AR-1 lands on HG's 8 V (8.9) — but
**regresses the floor-flux gate 0.060 → 0.071**: the whole flux curve shifts up (over-predicts at
every AR) because the looser −3·Te clip was partially masking a real flux over-prediction. Kept the
physical bound anyway (removing the fudge clip is the honest call); the regression **confirms the
open root cause is the RF-burst electron under-supply**, not the insulator clip. Production
`charging_floor_profile` table is the closure config (0.039) and is unaffected; `insul_vmin_Te` only
touches `solve_trench_charging` (the mechanism-study solver), so notching and de Boer runs are
untouched.

**Full RF-phase electron trajectories — implemented, measured, REVERTED (2026-07-03).** The plan's
hypothesized fix for the deep-AR over-charge: replace the first-order burst weighting (thermal
`gamma(2,Te)` energy + `cos^p` angle + an ad-hoc residual sheath barrier `frac·V_s` inside the
feature) with the physically-correct sheath-crossing entry. Derivation: a bulk electron crosses the
instantaneous sheath `V_s(φ)` only if its vertical energy exceeds `eV_s` (Boltzmann); the
flux-weighted 1-D vertical energy is `Exp(Te)`, so the **residual after climbing the barrier is again
`Exp(Te)`** (memoryless) — phase-independent — with 2-DOF transverse energy `Exp(Te)` preserved and
**no residual barrier inside** (already climbed). This is the correct free-flight-after-sheath picture.
**Measured (AR 4, 8000×140): V_c = 61 V (worse than the first-order 43.5; HG 33), floor flux 0.142,
~2× slower.** Why worse: the correct sheath-crossing arrival is **wide-angle** (median ~45°: `Exp(Te)`
transverse vs `Exp(Te)` vertical), and those wide-angle low-energy electrons either strike the upper
sidewalls or enter long/non-terminating orbits in the attractive floor well that the collisionless
finite-step trace drops — so the deep floor is supplied *less*, not more. The first-order model's
narrower `cos^p` angle was flattering the number (43.5), not solving the physics. **Conclusion (the
plan's "stop rather than tune"):** the deep-AR electron deficit is a **collisionless-2-D FRAMEWORK
limit** — electron collection into a deep attractive well is not resolvable by a collisionless
ballistic trace — NOT the barrier/energy approximation the plan hypothesized; the correct physics
*exposes* the deficit rather than fixing it. The real fix is a collisional or 3-D electron model
(out of scope). Reverted to the first-order build (the committed mechanism config); the production
closure table (0.039) is unaffected. The correct-physics finding is the deliverable.

**CORRECTION (same day, post literature review): the "framework limit" conclusion was premature.**
HG themselves got V_c = 33 V at AR 4 in a 2-D collisionless MC (JVST B 15,70 — arrival distributions
from a 1d/2v MC sheath simulation); Memos & Kokkoris (Micromachines 9,415 (2018)) trace fully
ISOTROPIC Maxwellian electrons with an adaptive integrator without any orbit problem, and their
secondary-electron-emission channel HALVES the charging potential (45→~22 V) by redistributing
electrons into the positive well; and for a trench v_y is conserved, so 2-D x–z dynamics is exact.
The deficit decomposes into (a) our fixed-step/step-cap integrator silently dropping well-captured
electrons, (b) missing SEE, (c) analytic arrival shortcuts. All fixable in 2-D collisionless — the
plan is `CHARGING_PHYSICS_PLAN.md` (W1 integrator → W2 SEE → W3 1-D RF-sheath MC source).

**W1 integrator landed and measured (2026-07-03).** Replaced the fixed-step capped trace with an
adaptive kick-drift-kick trajectory integrator and a Numba-compiled parallel particle loop
(`trace_integrator="adaptive_numba"`, falling back to the pure NumPy adaptive path when Numba is
absent). The numerical invariant is now explicit in the gates: max survivor fraction is <0.1%
(`charging_gate_result.npz`: ion 0.0000, electron 0.0006; `notching_gate_result.npz`: ion 0.0000,
electron 0.0007), so the silent well-captured-electron drop is gone. Runtime is usable again on the
Vast EPYC/RTX box: full charging rows run ~15-18 s instead of ~190 s/row for the pure NumPy adaptive
path. Physics result: W1 alone **does not close the gate**. Floor-flux RMSE is 0.075 (gate 0.05),
AR4 V_c is 42.9 V (HG 33), Matsui 300 eV still passes (0.573), and the notching mechanism still
fails the deep-AR foot-energy rise (15.8 vs HG 27.5 eV at AR4) while foot flux and poly potential
remain green. Conclusion: the integrator artifact was real and is fixed, but the dominant residual
is now W2/W3 physics — secondary electron emission and/or sheath-MC source distributions — not
silent trajectory drops.

**W2 PR-only PMMA SEE implemented and measured (2026-07-03).** Added an opt-in
`see_model="pmma_pr"` branch to `solve_trench_charging`: Dapor/Memos PMMA total electron yield
digitized over 0-100 eV, Burke polymer backscatter `eta(E)`, zero true secondaries below 16 eV,
1 eV cosine-isotropic true secondaries, elastic wall backscatter, and explicit charge bookkeeping
at PR sidewall emission sites. The default remains `see_model="none"` so the W1 baseline and
production closure tables are untouched. Measured on the Vast EPYC/RTX box with
`PETCH_SEE_MODEL=pmma_pr`, generation cap 1:
`charging_gate_see_result.npz` gives floor-flux RMSE **0.091 FAIL** (worse than W1 0.075), AR4
`V_c = 40.6 V` (only a small improvement from W1 42.9; HG 33), Matsui 300 eV **PASS** (0.575),
and max electron survivor **0.000906 PASS**. `notching_gate_see_result.npz` gives foot energy
15.4/16.3/18.2/19.1/18.7/17.8/16.4/16.5 eV vs HG 15/16.5/17.5/20/23/25/26.5/27.5: low/mid AR
improves, but deep AR still collapses; gate A **FAIL** (max rel err 40%), foot-flux constancy
barely **FAIL** (2.05), poly potential **PASS** (4.6 -> 35.3 V vs HG 6 -> 39). Conclusion:
PR-only PMMA SEE is a real, sign-correct partial channel, but it is **not** the missing deep-AR
lever. The remaining issue is more specific: the poly conductor potential is already right while
deep foot energy is low, so the residual likely lives in the joint arrival distributions / wall
material SEE / electron landing redistribution, not conductor charge sharing or trajectory drops.
Do not keep re-running PR-only SEE as a standalone fix; next falsification targets are all-wall
material SEE and the HG-style RF sheath-MC source.

**W2 code-review fixes and AR4 falsification (2026-07-03).** Review found the first PR-SEE branch
was incomplete: `see_generations` was only an on/off switch, yield > 1 was clipped to one sub-unity
emission probability, emitted electrons relaunched from bin centers rather than actual impact `z`,
and diagnostics undercounted absorbed wall hits. Fixed these in the opt-in branch only: true
multi-generation cascades, integer secondary multiplicity for `delta > 1`, elastic/specular primary
backscatter from the actual sidewall impact height, and foot-hit energy diagnostics
(`foot_E_p50`, `foot_E_p90`, `foot_z_mean`) in `diag["trace"]["last_ion"]`. Full AR4 probe
(8000 x 110, seed 10): default W1 gives floor flux 0.269, `V_c = 43.62 V`, `V_poly = 36.62 V`,
foot energy **15.77 eV** (p50 11.56, p90 29.12); corrected `see_model="pmma_pr",
see_generations=3` gives floor flux 0.319, `V_c = 36.41 V`, `V_poly = 32.99 V`, foot energy
**16.70 eV** (p50 12.79, p90 32.61), 2187 emitted electrons in the final 32000-electron trace,
zero survivor leakage. Interpretation: complete PR-sidewall SEE is sign-correct and moves the
voltage strongly toward HG, but it still raises deep foot energy by only **0.93 eV** vs the
remaining ~11 eV miss and worsens floor flux. This falsifies PR-sidewall SEE as the standalone
deep-AR fix. The foot-hit diagnostics show energetic ions exist in the tail (p90 ~30 eV) but the
mean is dominated by low-energy foot hits; the next problem is source/trajectory selection, not
silent orbits or scalar conductor potential.

**W3/source probe exposed the real geometry miss (2026-07-03).** Added an opt-in reduced
`source_model="sheath_mc"` interface (`src/petch/sheath1d.py`) for joint RF-sheath-ish arrival
sampling; AR4 probes moved floor flux and `V_c` but only shifted foot energy by ~0.5 eV, so source
shape alone is not the missing 11 eV. Re-reading HG JAP 82,566 resolved the contradiction: Fig. 4's
energy is the **average incident energy at the inner poly-Si sidewall of the edge line**, and Fig. 6
attributes the acceleration to the **difference between the edge-line and neighboring-line poly-Si
equipotentials**. petch's mechanism cell had tied both poly sidewalls to one periodic equipotential,
which cannot create this lateral tilt. A split-conductor diagnostic left the two sides nearly equal
(AR4 delta ~0.24 V), confirming that periodic symmetry lacks HG's open-side edge-line current. An
imposed line-to-line bias diagnostic finally moved the relevant side strongly: at AR4, 0/5/10/15 V
poly bias gave sidewall-foot means ~19.1/20.5/22.5/24.3 eV. Conclusion: the next real fix is an
HG edge-line/two-poly geometry and current balance (edge line supplied from the open outer side),
not more PR SEE or analytic source tuning. Details are in `CHARGING_DEEP_ISSUES.md`.

**Edge/open auxiliary boundary — implemented, measured, NOT sufficient (2026-07-03).** Added
`poly_mode="edge_open"` diagnostics plus `edge_open_model="line_of_sight"`: an open-half-space
free-flight model computes gross electron flux and ion counterflux to the outer edge-line poly
sidewall; negative edge potential suppresses electron collection by transverse-energy survival
`erfc(sqrt(|V_edge|/T_e))`, positive edge potential caps at the gross ballistic supply. This removes
the fake scalar current knob and matches HG Fig. 3 gross outer electron flux well (new figure
`viz/edge_open_current.png`; gross-electron RMSE ~0.02). Official gate with PR SEE:
floor-flux RMSE **0.076 FAIL**, foot-energy max error **32% FAIL** though rising, foot-flux
**PASS** (max/min 1.13), neighbor-poly potential **FAIL** (45% max error), survivor/current
residual/Matsui **PASS**. Interpretation: the auxiliary current model is a useful diagnostic and
reference boundary, but a scalar open-side boundary cannot replace the HG nonperiodic multi-line
field solve. The next real implementation is explicit open area + edge line + neighboring line
geometry with nonperiodic x electrostatics; do not tune the net current law.

**Explicit edge-array geometry — implemented as WIP, reduced gates measured (2026-07-03).**
Added `solve_edge_array_charging`: open area + edge poly line + trench + neighboring line,
nonperiodic x Laplace, explicit segment hits, separate edge/neighbor conductors, and optional
line-of-sight feature-mouth boundary source applied to the real edge conductor. Reduced W16/mouth80
8-point run with `edge_open_model="line_of_sight"`: floor RMSE **0.098 FAIL**, survivor **PASS**,
Matsui **PASS**, current residual **FAIL 0.235**. The notching gate on the same reduced config now
has foot-energy **PASS** (max error 28%, rising) and foot-flux **PASS** (max/min 1.67), while
poly potential still **FAILS** (51% max error from low-AR overshoot) and residual **FAILS 0.192**.
AR4 reduced point can hit the core mechanism numbers (floor ~0.27, Vc 42.5 V, edge/neighbor
18/37 V, foot energy 22.5 eV, residual 0.047), and W24/mouth160 AR4 gives gross
outer electrons 0.175 vs HG 0.18 plus neighbor 39.1 V, but the tail/final floor flux and PR residual
are not steady. Conclusion: the missing conductor geometry is now coded and qualitatively works, but
the edge-array relaxation/PR-surface balance must be stabilized before running expensive W32 gates
or claiming closure.

## Follow-ups

- ✅ Warp-ify the DDA neutral gather — done (`_dda_gather_kernel`, ~14× over numpy on CPU, 0.1 s/eval on CUDA).
- ✅ Tested petch-DDA with de Boer process knobs (narrow IADF, etchant starvation) at W=2 — they do NOT
  close the wafer gap (narrow IADF goes the wrong way; the deep-AR floor is structural).
- ✅ Knudsen floor sink: bug fixed (contour-snap), calibrated to the wafer (wls=1.3, RMSE 0.053 — near-pass,
  2.4× better than ballistic). Remaining shape residual is structural (linear decay vs flattening tail).
- Diagnose petch radiosity's over-correction.
- ✅ Tested the self-limiting-sink hypothesis (`knudsen_sink='field'`: clamp the sink at the field
  coverage so it stops growing as the floor starves). Result: RMSE 0.050 at wls=2.0 — marginally better
  than local (0.053) but the tail *shape* is unchanged (~2× decay per AR-doubling vs the wafer's ~1.45×).
  **Falsified**: the steepener is not the neutral-sink feedback. The concentration profile is ~1/(1+kAR);
  nr decays faster because of the multiplicative ion-side decay (m_i drops with AR; starved-regime rate
  ∝ θ_F ∝ conc × ion term).
- The real frontier, re-sharpened: the wafer's flat tail needs something that keeps the *deep-floor etch
  efficient* at high AR — the coverage-independent etch channel from FINDINGS (a reaction-limited floor
  term that doesn't die with θ_F), NOT a transport/sink knob, and NOT charging (charging cuts deep ion
  flux and would steepen the tail further). (The old `cal_F` knee-fix was removed; etchant starvation
  alone is insufficient given ion-enhanced floor etching.)
