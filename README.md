# plasma-etching-code

Working repo for a **differentiable, GPU-accelerated 3D feature-scale plasma-etch simulator**. This is where the implementation lives; the shareable docs/PoC mirror is the [`plasma-etching`](https://github.com/sdelaurentiis123/plasma-etching) repo.

The thesis in one line: keep the *geometry* physical (level set + ballistic transport + radiosity, all differentiable) and make the *physics inputs* learned operators (plasma → fluxes, and an MD-grounded surface-response operator). The result is one fast, accurate, **invertible** model from atoms to reactor — defensible on three axes nobody ships together: equal-fidelity-at-far-fewer-rays, differentiability (inverse design), and learned operators.

## Contents

| Path | What |
|---|---|
| [`docs/`](docs/) | **Technical explainer set** — open [`docs/index.html`](docs/index.html). 9 cross-linked pages: physics, numerics, flux/transport, acceleration, differentiable & ML, multiscale, SOTA, references. Real math (MathJax), derivations, Kushner/Graves lineage. |
| [`docs/references.bib`](docs/references.bib) | 45 BibTeX entries (human-readable: `docs/references.html`). |
| [`Plasma Chemistry Interaction Proof/`](Plasma%20Chemistry%20Interaction%20Proof/) | The 2D proof-of-concept: `feature_etch.py` (level-set + Monte-Carlo flux SF₆/O₂ etcher), `run_benchmark.py` (head-to-head vs ViennaPS `SF6O2Etching`), `summary.json` + `etch_benchmark.png` (cached results). The starting point we refactor and extend. |

## The two problems (kept separate)

1. **Statistical noise** — Monte-Carlo sampling roughness (`~1/√N`). Cured by *counting*: QMC, control variates, radiosity. Exact, no fidelity cost.
2. **Model-form & numerical bias** — the `~7%` depth offset currently absorbed by one global `rate_scale = 0.29`. Cured by *modeling*: exact Belen coverage, yield-over-IED, structured angular yield, ion reflection, higher-order advection.

Six identified bias contributors stand behind `rate_scale`; retiring each one (with a measured number, no knob) is the point of Phase 0.

## Phase plan

- **Phase 0** — 2D fidelity/convergence harness (CPU / Warp-CPU). Retire the 6 bias contributors; stand up DDA + QMC + control-variate + ion/neutral solver split; benchmark against cached ViennaPS numbers.
- **Phase 1** — 3D + NVIDIA Warp. DDA-on-level-set vs BVH-triangle, per species. Head-to-head vs live ViennaPS (Chamfer / depth / sidewall / ARDE + rays/sec).
- **Phase 2** — learned operators for the stiff chemistry source + multiple-scattering integral, trained on the validated baseline. The plasma-chemistry foundation-model convergence.

## Dev environment

Local box is Apple-silicon (no CUDA/ViennaPS). Phase 0 runs on CPU / Warp-CPU; GPU validation (Phase 1+) is pending. To run the existing PoC:

```bash
pip install warp-lang scikit-fmm numba scikit-image numpy scipy matplotlib
```

`viennaps` is only needed to regenerate ground truth (requires a GPU box); the cached reference is in `Plasma Chemistry Interaction Proof/summary.json`.
