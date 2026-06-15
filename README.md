# petch — a fast, differentiable 3D plasma-etch simulator

**The only open-source GPU-accelerated, differentiable feature-scale plasma-etch simulator** —
and on the same GPU it runs **~14× faster than ViennaPS** (the open-source SOTA) at matched accuracy.

Level-set surface evolution + Monte-Carlo / radiosity flux transport + SF₆/O₂ surface chemistry,
with the flux and level-set kernels written in [NVIDIA Warp](https://github.com/NVIDIA/warp) — so the
whole pipeline is GPU-resident **and** autodifferentiable in one substrate. Runs on an NVIDIA GPU
(fast) or on CPU / Apple Silicon (portable, slower) with the *same code*.

## Headline numbers

**Speed — same RTX 3090, both engines warmed, matched etch depth, swept across aspect ratios:**

| hole Ø | aspect ratio | ViennaPS-GPU (OptiX RT-core) | **petch** | speedup |
|---|---|---|---|---|
| 4 µm | 2.1 | 19.4 s | **1.30 s** | **14.9×** |
| 6 µm | 1.6 | 22.7 s | **1.61 s** | **14.1×** |
| 8 µm | 1.2 | 25.8 s | **1.82 s** | **14.1×** |

(Conservative — petch etched slightly *deeper*, i.e. more work. Reproduce: `scripts/vps_sweep.py`.)

- **As accurate as ViennaPS**: replicates every ViennaPS mechanism (Belen coupled coverages, exact
  Russian-roulette weighted neutral transport, coverage-dependent sticking, 1-neighbour flux
  smoothing); 3D ARDE rmse **0.05–0.08** vs depth-resolved ViennaPS, 2D rmse **0.016** parameter-free.
- **Differentiable** end-to-end (`wp.Tape`): gradient-based inverse design of process recipes — which
  ViennaPS cannot do.
- **Portable & modular**: every GPU speedup auto-enables on CUDA and falls back to numpy/skimage on
  CPU. The loop is clean swappable stages (`mesh → flux → chemistry → advect → reinit`).

## Install

```bash
pip install -e .            # core: numpy scipy scikit-fmm scikit-image warp-lang
# optional, for the last bit of GPU speed (edge-adjacency sort):
pip install cupy-cuda12x
```

Warp runs on CPU (Apple Silicon included) and on NVIDIA CUDA. No GPU needed to try it.

## Quickstart

```python
import petch
from petch import threed as t3

# THE fast + accurate config — one flag. On a GPU this auto-enables the full GPU pipeline
# (GPU marching cubes, on-device flux, GPU warm-start, ...) and runs ~14x faster than ViennaPS.
flags = petch.Flags(coverage_sticking=True, warm_start_coverage=True, sampling="sobol")

geo = t3.run_etch_3d(trench_width=6.0, dx=0.25, n_steps=40,
                     Lx=14, Ly=14, Lz=24, mask_th=2, sub_top=18, hole=True, t_end=3.0,
                     par=dict(petch.PAR), flags=flags, n_ion=30000, n_neu=30000,
                     reinit_method="fsm")
print("center depth:", t3.center_depth_3d(geo), "µm")
```

On CPU the *same call* works (drops to the numpy/skimage path automatically). See
[`examples/`](examples/) for a runnable script.

## How it works (and how to extend it)

Each step is a swappable stage, so adding physics is local:

```
make_trench_3d → [ marching cubes → MC/radiosity flux → SF6O2 chemistry → upwind advect → reinit ] × N
```

- **Flux** (`mc_flux_3d_coupled`): Warp ray-traced ions + Russian-roulette neutrals on a BVH, with a
  flux↔coverage fixed point. `neutral_transport="radiosity"` swaps in a deterministic conductance solve.
- **Chemistry** (`chemistry.py`): Belen/Ertl coupled F/O coverages + Steinbrüchel √E yields; swap via
  `Flags(chemistry=...)`.
- **Level set** (`reinit_method`): GPU Jacobi Godunov-Eikonal reinit (`fsm`), or CPU skfmm narrow-band.
- **GPU speedups** are independent flags (`gpu_mesh`, `gpu_source`, `gpu_warmstart`, `device_flux`,
  `flux_smooth_gpu`) — all auto-on under CUDA, each individually overridable.

## Beyond ViennaPS (experimental)

petch includes physics ViennaPS omits, behind flags, **clearly marked experimental and not yet
calibrated to wafer data**:
- `redeposition` — etch-product redeposition → sidewall passivation/taper.
- `surface_charging` — differential electron/ion charging (Hwang–Giapis). *Currently over-predicts the
  AR rolloff; needs a better electron angular model. Off by default.*
- bimodal IEDF yield integration (`ied_mode`).

## Honest limitations

- **As-accurate-as-ViennaPS, not yet validated to real wafers.** Matching ViennaPS ≠ matching a fab.
  The real-wafer gap (surface charging, true Knudsen molecular-flow transport) is research-grade and
  needs experimental data to calibrate.
- ViennaPS is the only runnable open-source peer; other benchmarks are published experimental curves
  (Belen 2005, de Boer/Blauw, Gomez 2004, Hoekstra–Kushner).

## Repo

```
src/petch/       2D (params/geometry/transport/chemistry/levelset/driver) + 3D (threed.py)
scripts/         vps_sweep.py (the 14x benchmark), validate_*.py, inverse_design.py, profilers
tests/           parity + 3D smoke
docs/            design explainers (open docs/index.html)
FINDINGS.md      full research log with every measured number
```

License: MIT (see `LICENSE`). Benchmarks compare against ViennaPS (GPL-3.0) but ship no ViennaPS code.
