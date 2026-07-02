# Reconciling petch with Craig's plasma_sim (2026-06-29)

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
| Gomez absolute rate (1.3 µm/min) | ✅ one global `rate_scale = 3.78e-4`, constant across features to **1.1%** |
| ViennaPS static reference, W=0.5 (pre-carved MakeTrench, makeMask=False, carve verified d0=D; reflecting ions) | measured: nr = **0.911 / 0.820 / 0.728 / 0.626 / 0.534** at AR 2/4/6/8/10 |
| petch-DDA vs that reference | **0.08 → 0.21 gentler, growing with AR** (0.742 vs 0.534 at AR 10). The DDA's diffuse re-emission over-delivers deep neutrals — its open calibration item (same family as the radiosity over-correction). The earlier "~0.08 gentler vs 0.73" used an unverified reference point, now replaced by the measured curve. |
| Evolving-etch Knudsen consistency | ❌ **OPEN** — as-configured the evolving run underperforms the static rate ~10× (AR 0.4 in 45 min at a nominal 1.3 µm/min field). Root cause undiagnosed (suspects: mask/sink interaction in `knudsen_face_flux` during evolution, or coverage collapse in the evolving fixed point). Two earlier attempts were invalid for script reasons (unphysical rate_scale; CFL-capped). Do not claim evolving-knudsen fidelity until this is closed. |

Notes: the wafer-gate PASS (knudsen path) is unaffected by the DDA gap — they are different neutral
solvers. Calibration automation: gradient/least-squares through the differentiable pipeline is the
right tool for these 1–3-scalar fits (RL is the wrong tool — no sequential decision structure); the
guard against overfitting is held-out gates, as above.

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
