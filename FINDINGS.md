# Phase 0 — Findings log

Running record of measured results as we retire the bias contributors. Numbers are on this
M1 Mac (CPU), 2D, vs the cached ViennaPS reference (width-8 depth 10.05 µm; normalized ARDE
`[0.932, 0.969, 0.985, 1.000]` for widths `[4, 6, 8, 12]`).

---

## Step 1 — Baseline parity
The `petch` package reproduces the original PoC exactly: width-8 center depth **9.18 µm** with
`rate_scale=0.29` (cached PoC value 9.40 — the 0.2 µm difference is MC run-to-run variance).

## Step 2 — Convergence harness
- **Grid:** depth plateaus by `dx=0.25` (0.25→0.125 changes depth 0.018 µm, below the
  seed-noise 0.063 µm). So the benchmark grid is essentially converged; contributor #5's
  *grid-resolution* part is small at `dx=0.25`. The −8.7% gap is therefore **not grid**.
- **Rays:** MC noise (across-seed std of depth) falls 0.25 → 0.12 µm over N = 5k → 40k,
  roughly `1/√N`. This is the noise floor QMC (Step 4) must beat.

## Step 3 — Contributor #1 (Belen coupled-coverage chemistry)
Implemented the exact ViennaPS Belen model (`src/petch/belen.py`): two coupled coverages with
`k_sigma=300`/`beta_sigma=0.04`, the standalone chemical-etch term `k_sigma·θ_F/4`, and
ViennaPS sticking `β_E=0.7 / β_O=1.0` (transport re-emission made consistent in the driver).

| model | width-8 depth | normalized ARDE | ARDE rmse |
|---|---|---|---|
| ViennaPS (ref) | 10.05 | 0.93, 0.97, 0.99, 1.00 | — |
| Langmuir @ rate_scale=0.29 | 9.18 (−8.7%) | 0.88, 0.93, 0.95, 1.00 | 0.038 |
| **Belen @ rate_scale=0.10** | **10.24 (+1.8%)** | 0.75, 0.86, 0.92, 1.00 | **0.110** |

**Result — two-sided:**
1. **Belen fixes the absolute depth** (+1.8% vs −8.7%) at its calibrated scale.
2. **Belen makes the ARDE lag *steeper* (worse)** — rmse 0.110 vs 0.038.

**Diagnosis (the valuable part):** Belen adopts ViennaPS sticking `β_E=0.7` (vs the PoC's 0.2).
Higher sticking consumes F near the trench top, so narrow/deep trenches get F-starved → steeper
lag. But **ViennaPS uses 0.7 too and still gets a flat ARDE**, which means our **neutral
transport over-depletes narrow trenches** (re-emission is insufficient — only 12 truncated
bounces, pseudo-random). The PoC's artificially-low `s_F=0.2` was *masking* this transport
deficiency by letting more F survive to the bottom.

**Conclusions:**
- Contributor #1 (chemistry) is **entangled with neutral transport** via the sticking
  coefficient; it cannot be cleanly retired in isolation.
- The real lever for ARDE is **better neutral re-emission** — the radiosity neutral solve
  (Step 4). That step now does double duty: variance-reduction speedup **and** the ARDE fix
  (restore flux to narrow-trench bottoms → flatten ARDE back toward ViennaPS).
- `rate_scale` (0.10 for Belen, 0.29 for Langmuir) is confirmed a **units / flux-normalization**
  issue (our `mc_flux` emits dimensionless open-field multipliers, not absolute ViennaPS fluxes
  in 1e15 cm⁻²s⁻¹), **separate from chemistry form**. Killing it needs flux-unit reconciliation,
  not a better coverage model.

**Next:** Step 4 — implement the neutral radiosity solve (all-bounces, no truncation) and
re-measure ARDE for Belen. Hypothesis: ARDE rmse drops back toward / below the Langmuir 0.038
once narrow-trench flux is restored.
