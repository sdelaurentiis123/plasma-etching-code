# Krüger 2024 bounded CUDA and refinement report

Date: 2026-07-17
Status: complete bounded preflight; no endpoint calibration authority and no held-out reveal

## Question answered

Before another long Krüger run, this audit asked three narrower questions:

1. Does putting both particle transport and level-set work on one CUDA device change the answer?
2. Where does a warmed 5 nm feature step actually spend wall time?
3. Which initial 10 nm versus 5 nm observables are already stable, and which remain grid-sensitive?

It reused the already-complete 0.5 s paired refinement run. No experimental held-out profile was
opened, no calibration parameter was changed, and no additional endpoint trajectory was launched.

Machine-readable evidence and the static plot are in
`results/krueger_2024_cuda_profile_summary/`.

## CUDA-device correction and parity

The first bounded profile set `transport_device=cuda:0` but left `PETCH_DEVICE` unset. Transport
therefore ran on the RTX 4090 while level-set redistance followed the module-import default on CPU.
The worker launcher now exports both declarations from one `--device` input and the profiler refuses
mixed devices.

The corrected run used:

- transport device: `cuda:0`;
- level-set device: `cuda:0`;
- one positive warm-up step and one profiled positive step;
- the same seed, 5 nm grid, operator, calibration pair, and total 0.05 s physical trajectory as the
  mixed-device reference;
- exact-zero material-ledger residual;
- maximum diffuse-radiosity balance residual `1.90e-12`;
- no topology event.

Endpoint parity is substantially tighter than any physical or spatial uncertainty in this campaign:

| Quantity | Absolute mixed-vs-unified difference |
| --- | ---: |
| Etch depth | `0 nm` |
| Mask opening | `8.59e-9 nm` |
| Maximum/top feature width | `1.02e-7 nm` |
| Mask-top height | `9.42e-9 nm` |
| Maximum per-y opening | `2.53e-8 nm` |

The maximum common-unit difference is `1.02e-7 nm`, below the report's `1e-6 nm` roundoff-scale
parity check. Unified CUDA device selection therefore passes this bounded numerical-equivalence gate.

## Warmed step cost

The corrected one-step wall time was `7.780 s`; the cumulative time inside
`advance_feature_step_3d` was `7.762 s`.

| Mutually exclusive top-level work | Time (s) | Share |
| --- | ---: | ---: |
| Chemistry and material routing | 2.292 | 29.5% |
| Ballistic boundary transport | 2.233 | 28.8% |
| Diffuse neutral exchange | 1.184 | 15.3% |
| Legacy surface-state remap | 1.064 | 13.7% |
| Remaining extraction/advection/topology/overhead | 0.931 | 12.0% |
| Level-set redistance | 0.059 | 0.8% |

This changes the optimization order. Moving redistance to CUDA was necessary to make device
provenance honest, but redistance is not the runtime bottleneck. The large remaining costs are mostly
Python/NumPy chemistry construction, boundary-transport setup and gathering, diffuse view-factor
work, and the legacy point-to-triangle remap certification.

A purely arithmetic projection of `2400` fixed 0.025 s steps gives `5.19 h` for 60 s physical time.
That is **not** an endpoint runtime prediction: it is based on one warmed initial-geometry step, while
surface size, access, topology, and adaptive timestep schedules change during a real process. It is
only enough to reject another blind long run before the dominant kernels are addressed.

## Initial 10/5 nm refinement

The existing 0.5 s paired run used the same physical cell, calibration pair, seeds, boundary
quadratures, and the 5 nm-owned timestep schedule.

| Rate observable | 10 nm | 5 nm | Absolute relative difference |
| --- | ---: | ---: | ---: |
| Etch depth (nm/s) | 13.6573 | 13.6672 | 0.073% |
| Definition-correct mask opening (nm/s) | -4.4285 | -4.3659 | 1.433% |
| Remaining mask thickness (nm/s) | 0.5023 | 0.4972 | 1.033% |
| Maximum feature width (nm/s) | -4.4493 | -15.4678 | 71.235% |
| Top feature width (nm/s) | -4.4493 | -15.4678 | 71.235% |

The scientifically narrow conclusion is:

- depth, mask opening, and mask-thickness rates are sufficiently close to construct a
  coarse-to-fine calibration-discrepancy proposal;
- maximum/top width is strongly grid-sensitive at this shallow state and cannot support a local
  profile-shape or matched-error superiority claim;
- this preflight does not certify a 60 s endpoint, a topology-changing late state, or AMR.

## Remap finding

The Krüger worker did not explicitly select a surface-state transfer backend, so it inherited
`surface_state_remap_backend="legacy_knn"`. That path accounted for `1.064 s` and its
point-to-triangle maximum-distance check accounted for almost all of that cost. The newer indexed,
partitioned-overlap, and common-refinement operators exist and pass their manufactured suites, but
this run did not use them. They must be selected explicitly in the next short paired operator test;
their correctness and cost may not be inferred from the legacy profile.

## Earned next work

In order:

1. expose the remap backend as a checksum-bound Krüger worker input and run a bounded legacy versus
   common-refinement/indexed comparison on the same short state;
2. remove repeated construction and subsetting in the chemistry router, preserving byte/roundoff
   parity and every material ledger;
3. retain device geometry/sample arrays across ballistic and diffuse calls where the surface has not
   changed beyond a declared reuse tolerance;
4. use the paired 10/5 evidence to implement the multi-fidelity trust-region *controller* while
   keeping every held-out outcome sealed;
5. launch a new long base trajectory only after the selected transfer operator and short-step
   performance/refinement gates close.

This is a performance and numerical-authority report, not experimental validation.
