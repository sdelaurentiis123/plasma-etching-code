# plasma-etching-code

A **differentiable, GPU-accelerated, ViennaPS-calibrated** feature-scale plasma-etch simulator
(2D + 3D). Level-set surface evolution + Monte-Carlo / radiosity flux transport + SF₆/O₂ surface
chemistry, built so the flux kernel is a [NVIDIA Warp](https://github.com/NVIDIA/warp) kernel —
RT-core ray tracing **and** autodiff in one substrate.

Full research log with every measured number: [`FINDINGS.md`](FINDINGS.md). Design explainer:
[`docs/`](docs/) (open `docs/index.html`).

## What's built & verified

| | result |
|---|---|
| **2D fidelity** | matches real ViennaPS ground truth **parameter-free**: width-8 depth **+3%**, ARDE rmse **0.016**. |
| **3D accuracy** | **all of ViennaPS's mechanisms** (Belen coupled coverages, exact Russian-roulette weighted transport, coverage-dependent sticking, 1-neighbor flux smoothing) → 3D ARDE **rmse ~0.05–0.08** vs depth-resolved ViennaPS-3D, *bracketing* it. |
| **Speed (same-engine)** | **faster than ViennaPS-GPU**: matched d=6 hole **10.18 s vs ViennaPS-GPU 11.6 s** (was 22.6 s = 2.2× self-speedup). Shipped **narrow-band reinit** (exact, ~3–5×) + **GPU advection** (Warp kernel, on-device CFL substeps). Root cause of the old gap was dense-vs-sparse, not CPU-vs-GPU. |
| **Differentiable** | `wp.Tape` autodiff through the flux kernel; **inverse design** recovers a recipe from a target (E=139.9 vs 140 eV). ViennaPS is not differentiable. |
| **Physics beyond ViennaPS** | **full ion-energy-distribution** yield integration (mean / Gaussian=ViennaPS / **bimodal** sheath); **etch-product redeposition** (sidewall passivation → taper). Both omitted by ViennaPS. |
| **Real-wafer validation** | benchmark targets identified — Gomez/Belen 2004 SF₆/O₂, de Boer 2002 cryo RIE-lag + Blauw Knudsen model, Hoekstra-Kushner microtrench. See `docs/experimental-validation.html`. |

**The three together — faster than ViennaPS-GPU, as accurate, AND differentiable — is the edge**
(ViennaPS is neither GPU-resident in the level set nor differentiable).

**Honest open frontier** (also missing in ViennaPS, the path to beating *real wafers*): **surface
charging** (notching, AR-dependent ion deflection — biggest gap), then full bimodal-IEDF-tail
sensitivity and flux-dependent reaction probability. See `docs/` and `FINDINGS.md` for the roadmap.

## Layout

```
src/petch/        2D: params, geometry, transport (numba), chemistry (langmuir|belen),
                  belen, levelset, driver, metrics. 3D: threed.py (Warp flux + level set).
harness/          benchmark.py (2-axis scorecard), convergence.py (grid/ray + Richardson),
                  reference/ (cached + real ViennaPS 2D & 3D ground truth).
scripts/          run_phase0.py (2D scorecard), warp_spike / warp_3d_flux (kernel+autodiff),
                  gpu_benchmark.py, inverse_design.py, viennaps_*calibrate / groundtruth.py.
tests/            smoke + PoC parity + 3D loop.
```

## Run it

```bash
pip install -e .                      # numpy scipy scikit-fmm numba scikit-image matplotlib
pip install warp-lang                 # CPU on Apple Silicon; CUDA/RT-cores on NVIDIA
python -m pytest tests/ -q            # parity + 3D smoke
python scripts/run_phase0.py          # 2D scorecard vs ViennaPS
python scripts/warp_3d_flux.py        # 3D differentiable flux kernel (shadowing + autodiff)
python scripts/inverse_design.py      # target outcome -> recipe by autodiff
```

The **default config is the calibrated accurate model** (belen + ViennaPS angular yield +
`cal_F=12`). `Flags(chemistry="langmuir", yield_angular="cosine")` + `cal_F=1` recovers the
original 2D proof-of-concept (pinned by tests).

## Status & next

2D is parameter-free matched; 3D is fast + differentiable but the deep-HARC neutral-transport
fidelity gap is open. The companion docs repo [`plasma-etching`](https://github.com/sdelaurentiis123/plasma-etching)
holds the shareable explainer. Phase 2 (learned operators) is intentionally not started.
