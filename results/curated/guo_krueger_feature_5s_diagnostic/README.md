# Guo/Kwon–Krüger five-second feature diagnostic

Status: **completed development diagnostic; not a predictive validation
claim**.

This receipt pins the first uninterrupted feature trajectory using the
atom-balanced Guo/Kwon translating-layer closure, the unchanged common 3-D
transport/level-set engine, Krüger's published aggregate ion flux, and the
digitized published IEAD.  No oxide yield scale, ion-flux normalization, or
feature-depth calibration was used.

The constitutive implementation is commit `9540970`.  The run was made from
the same implementation immediately before that commit; the only worktree
changes were the four files recorded by `9540970`.

## Result

- Physical time: 5.000000 s
- Accepted adaptive steps: 38
- Etch depth: 66.0488716 nm
- Mean depth rate: 13.2097743 nm/s
- Paper-implied 60 s mean rate: 13.7500000 nm/s
- Rate difference over this five-second prefix: -3.92891%
- Initial / final mask opening: 90.0000000 / 57.1298271 nm
- Maximum recorded radiosity relative balance error: 1.94673e-12
- Maximum material ledger residual: 0 units/m²
- Terminal topology event: none
- Pilot status: complete

Selected trajectory points:

| time (s) | depth (nm) | mask opening (nm) |
|---:|---:|---:|
| 0.000000 | 0.000000 | 90.000000 |
| 0.500000 | 7.262331 | 80.280060 |
| 1.492211 | 21.063334 | 76.381379 |
| 2.403572 | 31.844828 | 72.999605 |
| 3.233428 | 44.368878 | 69.427047 |
| 3.889184 | 52.174325 | 64.531978 |
| 4.456653 | 60.333112 | 61.285383 |
| 5.000000 | 66.048872 | 57.129827 |

## Reproduction

```text
python scripts/krueger_2024_trench_pilot.py \
  --output /private/tmp/krueger_guo_5s_bdf \
  --duration-s 5 --n-steps 10 --dx-um 0.01 --n-position 1 \
  --surface-model guo_tml --face-quadrature-points 1 \
  --radiosity-rays 4 --radiosity-max-iterations 500 \
  --neutral-transverse-order 2 --neutral-direction-polar-order 2 \
  --neutral-direction-azimuthal-order 4 --ion-azimuthal-order 4 \
  --maximum-accepted-steps 40 \
  --topology-change-policy continue_gas_cavity --max-wall-s 300
```

Configuration hash:
`fe25f1856df422dda9614ebe4fee273e0b66f2f314acaf9ca5af1af939f1fea7`.

Original local artifacts and SHA-256:

- `audit.json`: `92be6c757efe8b8be8bbbafc889162b468e7b5911e9882c18a4305df3b79689e`
- `checkpoint.npz`: `77b009dbd3f8767156e12abc2e5812287d6b77bf650eb90a68a81822a97ce805`
- `profile.png`: `a59a7e43de37cde040f2a2d13eaa4b3a1f22c6a7fcb4adc3097fb171c71ba3c6`
- `metric_trajectory.png`: `a40a7eb4fd056032e9bfbc3aedc15debac75cbfcc468a1e9b454a249e0213ce6`

The local full artifacts remain under
`/private/tmp/krueger_guo_5s_bdf`; this repository receipt deliberately
records the compact, reviewable result rather than committing a
worktree-dirty checkpoint.

## Claim boundary

This prefix is encouraging evidence for the absolute scale, not a 60 s depth
match.  It is out of Guo's declared board because Krüger supplies an aggregate
ion population, C2F3/C3F4 are transferred through the generic-neutral
topology, most of the IEAD exceeds 370 eV, and the source's repaired physical
angular polynomial is used off-normal.  The amorphous-carbon mask closure also
contains transferred and unmeasured quantities.

The transport quadrature here is intentionally low order
(`n_position=1`, one face point, four ion azimuths, 2×4 neutral directions,
four radiosity rays).  A depth claim requires:

1. continuation to 60 s or a declared physical/topological terminal event;
2. spatial, angular, IEAD, and radiosity refinement;
3. separation of mask narrowing from oxide-side polymer growth;
4. measured/species-resolved boundary data or an independently validated
   reactor closure;
5. held-out profile/depth transfer rather than agreement of this one prefix.
