# TorchSim atomistic feasibility

This directory records a deliberately bounded evaluation of TorchSim 0.6.1 on the
multiphysics branch. It asks whether TorchSim can serve as an offline atomistic
surface-data provider for petch. It does **not** promote a new etch mechanism or alter
any certified depth result.

## Reproduction

TorchSim 0.6.1 requires Python 3.12, while petch supports Python 3.10. The scripts use
PEP 723 dependency blocks so the environments remain separate:

```bash
WARP_CACHE_PATH=/tmp/petch-warp-cache uv run --python 3.12 \
  scripts/torchsim_atomic_feasibility.py \
  --device cpu --steps 6000 --timestep-ps 0.00025 \
  --energies-eV 0.1 0.5 1.0 2.0 \
  --output results/curated/torchsim_feasibility/lj_argon_cpu_v1.json

WARP_CACHE_PATH=/tmp/petch-warp-cache MPLCONFIGDIR=/tmp/petch-mpl-cache \
  uv run --python 3.12 scripts/torchsim_mace_si_feasibility.py \
  --device cpu --steps 20 \
  --output results/curated/torchsim_feasibility/mace_mp_si_cpu_v1.json
```

The MACE script downloads the named public checkpoint. Its URL and SHA-256 are stored
in the output. Runtime fields are observations from this host, not portable benchmarks.

## Findings

- Batched deterministic NVE execution works on CPU.
- For the analytic Lennard-Jones model, the force equals the negative autograd energy
  gradient to the recorded numerical precision. That is only a static derivative.
- The stricter trajectory gate fails: after 50 NVE steps, the final separation has a
  nonzero finite-difference sensitivity to the Lennard-Jones length parameter
  (-2.758e-3 in the recorded units), while autograd reports no parameter gradient.
  TorchSim's pair-potential force path uses a detached force derivative in 0.6.1, so
  this configuration is not end-to-end differentiable through the trajectory.
- The toy impact ensemble serializes directly to the existing
  `SurfaceInteractionTable` payload and round-trips through petch.
- The toy potential produces an apparent energy-dependent ejection response, but the
  unshifted finite-cutoff Lennard-Jones dynamics lose 2.2--5.4% of incident energy in
  this sweep. The result is a software/contract probe and is forbidden as etch data.
- The real MACE-MP silicon batch places the lowest energy of the three sampled diamond
  lattices at 5.43 angstrom. The 20-step, 300 K NVE relative drift is
  4.1e-6--7.2e-6 for the three 64-atom systems.
- MACE 0.3.16's TorchSim adapter returns detached energy and force tensors. The tested
  MACE path is therefore deterministic batched inference, not an autograd-capable
  atomistic segment.
- On this Apple host, PyTorch exposes MPS, but TorchSim's default
  nvalchemiops/Warp neighbor list rejects the MPS device. CPU works; a CUDA host is the
  relevant acceleration target unless a different verified neighbor-list backend is
  selected.

## Physics boundary

MACE-MP is a useful integration check, not an etch potential. This evaluation contains
no evidence for energetic ion/radical impacts, reaction probabilities, product
branching, sputter yield, charging, SiO2, fluorocarbon films, masks, or reactor
conditions. Consequently neither generated artifact supports a petch prediction.

The repository already contains stronger etch-specific atomistic evidence than this
generic checkpoint: the Kounis-Melas Si-Cl-Ar DeepMD tables and the An et al. SiO2/Si3N4
NNP/ZBL beam comparison. TorchSim becomes scientifically useful only after it is paired
with a potential whose training support and direct-beam transfer are at least as well
audited. The An potential cannot presently be imported because its pinned repository
has no compatible license.

## Architectural verdict

Keep TorchSim outside the Python 3.10 production runtime as a versioned offline provider:

```text
DFT / licensed etch potential -> TorchSim impact ensembles
                               -> SurfaceInteractionTable artifact
                               -> petch feature calculation
```

The next meaningful gate is not a longer MACE-MP run. It is one exact-overlap replay
against a direct beam board using a licensed potential with explicit short-range
collision physics, followed by a held-out energy/species test. Independently, any claim
of differentiable impact dynamics must pass an autodiff-versus-finite-difference
trajectory gate; the current pair-potential path does not, and the current MACE adapter
detaches its outputs. Until those gates pass, TorchSim changes no reactor or depth
verdict.
