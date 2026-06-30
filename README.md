# petch — a fast, differentiable 3D plasma-etch simulator

**The only open-source GPU-accelerated, differentiable feature-scale plasma-etch simulator** —
and on the same GPU it runs **~8–15× faster than ViennaPS** (the open-source SOTA), tracking its ARDE
within ~0.1 at low/mid aspect ratio (see *Accuracy* below for the honest deep-AR difference). petch is
GPU-resident (~1 s, **CPU-independent**) while ViennaPS's level-set advection is CPU-bound, so the
*ratio* depends on the host CPU — ~8× when ViennaPS gets a fast many-core CPU, more on typical/weak ones.

Level-set surface evolution + Monte-Carlo / radiosity flux transport + SF₆/O₂ surface chemistry,
with the flux and level-set kernels written in [NVIDIA Warp](https://github.com/NVIDIA/warp) — so the
whole pipeline is GPU-resident **and** autodifferentiable in one substrate. Runs on an NVIDIA GPU
(fast) or on CPU / Apple Silicon (portable, slower) with the *same code*.

📖 **[Read the docs →](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/index.html)** — the full technical explainer set (physics, numerics, transport, the ViennaPS head-to-head, and experimental validation), rendered with math and figures. Full index in [Documentation](#documentation) below.

## Headline numbers

**Speed — same RTX 3090, both engines warmed, depth-matched, swept across aspect ratios.**
Measured on a **32-core / RTX 3090** box (so ViennaPS gets a *fast* CPU — the conservative case):

| hole Ø | aspect ratio | ViennaPS-GPU (OptiX RT-core) | **petch** | speedup |
|---|---|---|---|---|
| 4 µm | 2.1 | 7.4 s | **0.90 s** | **8.3×** |
| 6 µm | 1.6 | 7.4 s | **1.01 s** | **7.3×** |
| 8 µm | 1.2 | 12.2 s | **1.17 s** | **10.4×** |

Depth-matched within 3–7%. **The ratio is CPU-dependent** because ViennaPS's level-set advection runs
on the CPU (~40% of its time): on this fast 32-core box ViennaPS is 7–12 s → ~8×; on a mid/typical CPU
it's ~19–25 s → ~14×; on a weak few-core box it balloons to 50–64 s → 40–60× (an *artifact*, not a real
gain). **petch's ~1 s is GPU-resident and steady regardless** — so ~8× is the honest floor, and the
advantage only grows on weaker hosts. Reproduce: `scripts/vps_sweep.py`.

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

## Documentation

A cross-linked technical explainer set (real math via MathJax, figures, derivations) lives in [`docs/`](docs/). GitHub renders those `.html` files as raw source when clicked, so the links below open them rendered through [githack](https://raw.githack.com) — math, figures, and cross-page nav all work in the browser:

**▶︎ [Open the docs site](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/index.html)** — start at the Overview; every page is cross-linked from there.

| Page | What |
|---|---|
| [Overview](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/index.html) | The thesis, the noise-vs-bias framing, the full stack top-to-bottom |
| [Physics](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/physics.html) | Ion–neutral synergy, competitive coverage, sputter/ion-enhanced yields |
| [Numerics](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/numerics.html) | Level-set Hamilton–Jacobi PDE, WENO + TVD-RK, reinitialization |
| [Flux &amp; Transport](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/flux-transport.html) | Ballistic ion/neutral transport, multi-bounce radiosity |
| [Acceleration](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/acceleration.html) | Variance reduction, GPU traversal, sparse data structures |
| [Differentiable &amp; ML](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/differentiable-ml.html) | Autodiff, inverse design, learned operators |
| [Multiscale](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/multiscale.html) | Reactor → feature-scale coupling |
| [**ViennaPS Validation**](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/viennaps-validation.html) | Head-to-head vs ViennaPS-GPU — ARDE agreement + the ~47× real-time speed race |
| [**Experiments**](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/experimental-validation.html) | ARDE vs the real de Boer/Blauw cryo wafer |
| [Performance](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/performance.html) | Speed breakdown and where the time goes |
| [SOTA &amp; Plan](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/sota-and-plan.html) | Where this sits vs the field, and the roadmap |
| [Physics Grounding](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/physics-grounding.html) | First-principles basis and the Kushner/Graves lineage |
| [References](https://raw.githack.com/sdelaurentiis123/plasma-etching-code/main/docs/references.html) | 52 BibTeX entries, human-readable |

> Prefer a permanent URL without the proxy? Enable GitHub Pages on this repo (**Settings → Pages → Source: `main` / `/docs`**) and the same set is served at `https://sdelaurentiis123.github.io/plasma-etching-code/`.

## Repo

```
src/petch/       2D (params/geometry/transport/chemistry/levelset/driver) + 3D (threed.py)
scripts/         vps_sweep.py (the 14x benchmark), validate_*.py, inverse_design.py, profilers
tests/           parity + 3D smoke
docs/            technical explainer set — see Documentation above for rendered links
FINDINGS.md      full research log with every measured number
```

License: MIT (see `LICENSE`). Benchmarks compare against ViennaPS (GPL-3.0) but ship no ViennaPS code.
