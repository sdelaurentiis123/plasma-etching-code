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
| **2D fidelity** | matches real ViennaPS ground truth **parameter-free**: width-8 depth **+3%**, ARDE rmse **0.016** (from −8.7% / 0.110 that needed a per-case knob). Root cause was **flux normalization**, calibrated to `cal_F=12`. |
| **3D etcher** | full 3D loop (level set → marching cubes → Warp flux → chemistry → advection). Etches trenches + contact holes. |
| **GPU speed** | flux kernel **198 M rays/s on RT cores — 137× CPU**; production GPU reinit (Russo-Smereka) → fully GPU-resident loop; sub-second 3D etch. |
| **Differentiable** | `wp.Tape` autodiff through the flux kernel; **inverse design** demo recovers a recipe from a target (E=139.9 vs 140 eV). |
| **Speedups** | QMC (Sobol) source sampling (exact, fewer rays); ion/neutral split with deterministic radiosity neutral solve. |

**Honest open gap:** deep 3D HARC holes do **not yet** quantitatively match ViennaPS ARDE — our
MC neutral transport over-depletes deep narrow holes more than ViennaPS (which stays F-saturated).
A named model-form limit needing a better neutral-transport model. See FINDINGS §"3D fidelity gap".

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
