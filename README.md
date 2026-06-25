# petch — a fast, differentiable 3D plasma-etch simulator

**The only open-source GPU-accelerated, differentiable feature-scale plasma-etch simulator** —
and on the same GPU it runs **~14× faster than ViennaPS** (the open-source SOTA), tracking its ARDE
within ~0.1 at low/mid aspect ratio (see *Accuracy* below for the honest deep-AR difference).

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

- **Tracks ViennaPS**: replicates every ViennaPS mechanism (Belen coupled coverages, exact
  Russian-roulette weighted neutral transport, coverage-dependent sticking, faithful ion reflection,
  1-neighbour flux smoothing). 3D trench ARDE: petch and ViennaPS agree within **~0.1** at low/mid aspect
  ratio and diverge by **~0.1–0.2 at the deepest AR** (petch is slightly gentler — it delivers a bit more
  flux to the deep floor). This is a genuine *converged* difference between two independent ballistic codes
  (different source/transport discretization), not a calibration gap: both run the same documented model,
  and **both sit ~0.3 above the real de Boer wafer** (ballistic transport omits gas-conductance / charging).
  2D ARDE-shape rmse **0.016**. The ARDE *curve* is parameter-free; one global `rate_scale` sets only the
  absolute etch rate (the analog of ViennaPS's `unitConversion`), not the shape.
- **Differentiable** end-to-end: the flux and level-set kernels are written in Warp, so the pipeline
  carries gradients via `wp.Tape` — a substrate for gradient-based recipe optimization that ViennaPS
  cannot do. Demonstrated on a single-parameter inverse-design recovery in `scripts/inverse_design.py`
  (a one-call `petch.inverse_design(...)` API is not yet exposed — it's a demo script, not a product API).
- **Portable & modular**: every GPU speedup auto-enables on CUDA and falls back to numpy/skimage on
  CPU. The loop is clean swappable stages (`mesh → flux → chemistry → advect → reinit`).

## Install

```bash
pip install -e .            # core: numpy scipy scikit-fmm scikit-image warp-lang
# optional, for the last bit of GPU speed (edge-adjacency sort):
pip install cupy-cuda12x
```

Warp runs on CPU (Apple Silicon included) and on NVIDIA CUDA. No GPU needed to try it.

Verify the install:

```bash
pytest tests/            # 4 smoke/regression tests (2D parity, 3D engine, high-level API)
```

## Quickstart

The high-level API mirrors ViennaPS (`Domain` / process model / `Process`), with the full faithful
config (Belen coverages + ion reflection) built in. On a GPU it auto-enables the GPU pipeline and runs
~14× faster than ViennaPS; on CPU the same code drops to numpy/skimage.

```python
import petch

dom    = petch.Domain.hole(extent=14, dx=0.25, diameter=6, mask=2, depth=18)
model  = petch.SF6O2()                      # faithful SF6/O2; rate_scale calibrates absolute rate
result = petch.Process(dom, model, duration=3.0).run(steps=40)

print(f"depth {result.depth:.1f} µm   aspect ratio {result.aspect_ratio:.1f}   ({result.wall_time:.1f}s)")
result.save("etch.vtk")                     # ParaView / ViennaPS-readable surface mesh
```

For full control, the low-level `petch.run_etch_3d(...)`, `petch.Flags`, and `petch.PAR` stay public.
See [`examples/`](examples/) for runnable scripts.

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
- `surface_charging` — differential electron/ion charging (Hwang–Giapis). *Honest status: the reduced
  differential-shadowing model (electrons more HARC-shadowed than ions) does NOT reproduce the HG
  floor-current rolloff at any electron angular spread — it over-throttles. HG needs the self-consistent
  floor potential re-deflecting ions (a PIC-class field solve), not geometric shadowing. The
  infrastructure (electron trace, charge factor, `e_ang_sigma`) is in place for that future work. Off by
  default; the effect direction is right (throttles the floor) but it is not quantitatively calibrated.*
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
