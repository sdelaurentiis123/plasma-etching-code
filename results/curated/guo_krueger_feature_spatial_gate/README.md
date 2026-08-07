# Guo/Kwon–Krüger production-order spatial gate

Status: **early-time numerical gate passed; absolute-depth prediction gate
still failed/unidentified**.

This pair advances the no-depth-fit Guo/Kwon surface closure through the
unchanged common 3-D feature engine at 10 nm and 5 nm uniform spacing.  Both
use the same production-order boundary operator and conservative
common-refinement surface-state remap:

- 16 source-position samples;
- three face quadrature points;
- 8×16 neutral direction quadrature;
- 16 ion azimuths;
- 8 radiosity rays/face with 1e-12 balance tolerance;
- the compressed joint IEAD quadrature already endpoint-certified against the
  exact digitized table;
- published aggregate ion flux with normalization 1.0;
- Guo/Kwon surface parameters unchanged and no oxide-yield scale.

Implementation revision: `e3b9257`.

## Paired result at 0.5 s

| spacing | accepted steps | depth (nm) | mean rate (nm/s) | mask opening (nm) | wall time (s) |
|---:|---:|---:|---:|---:|---:|
| 10 nm | 2 | 5.30891814 | 10.6178363 | 85.1018661 | 67.60 |
| 5 nm | 3 | 5.06933169 | 10.1386634 | 84.6348241 | 786.82 |

Fine relative to coarse:

- depth: -4.513%
- mean depth rate: -4.513%
- mask opening: -0.549%

The 5 nm controller rejected the second 0.25 s proposal and used
0.157331037 s plus 0.092668963 s; the 10 nm run accepted two 0.25 s steps.
The comparison therefore includes converged adaptive time integration rather
than replaying a coarse schedule that violates the fine displacement gate.

Maximum recorded neutral-radiosity balance residual was 1.28e-12; material
ledgers closed exactly in both runs.

## Interpretation

The early 5/10 nm depth difference is below 5%, so the spatial gate passes at
this time and operator.  The cheap 10 nm scout used earlier
(`n_position=1`, one face point, 2×4 neutral directions, four ion azimuths and
four radiosity rays) gave 7.26233068 nm at the same time.  Production order
therefore lowers the 10 nm result by 26.90%; the scout must not be used for an
absolute-depth claim.

The fine production-order rate is 10.1387 nm/s.  Linear extrapolation is not a
prediction, but its 60 s scale is 608.3 nm, well below Krüger's 825 nm.
Nonlinear feature evolution and mouth equilibrium still require a long
production-order trajectory.  This gate therefore improves the numerical
credibility of the short-time operator while **not** establishing an absolute
depth match.

## Artifact receipts

- 10 nm config hash:
  `5f0013caa158f2656c63a9d3e409027a163b50c0b4758a66e2f780f5474b00b1`
- 10 nm `audit.json` SHA-256:
  `a90f88a006e5c4076bb4a05a91be0c0de7140f67e4c0bad1ff6567bf55b5e200`
- 5 nm config hash:
  `9b08f94847fb8bdc415e6bc7354f973c5696f1175c4cfe973cbccb2e36445562`
- 5 nm `audit.json` SHA-256:
  `21275ae08f2d3fcbe2a7c40e4496645c8c24cf9cd647fa156e82d55406acf195`

The full local artifacts are in
`/private/tmp/krueger_guo_refine_dx10` and
`/private/tmp/krueger_guo_refine_dx5`.

## Evidence boundary

The pair remains an out-of-board transfer audit, not Tier-A prediction:
Krüger's ion composition is aggregate; C2F3/C3F4 are mapped through Guo's
generic-neutral topology; most IEAD support exceeds Guo/Yin's 370 eV
regression board; and the off-normal physical angular polynomial includes the
declared source repair.  The transferred amorphous-carbon mask closure also
contains unmeasured parameters.  Spatial convergence cannot identify those
physical inputs.
