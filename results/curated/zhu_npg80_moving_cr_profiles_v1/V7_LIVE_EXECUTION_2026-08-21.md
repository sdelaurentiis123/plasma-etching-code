# Oxford moving-Cr v7 live execution receipt

Status at 2026-08-21 16:52 UTC: **11 of 56 trajectories complete and retrieved**. This is a
recoverable execution checkpoint, not a complete board and not a grade against
Freddie's withheld SEM.

## Frozen implementation

- Git commit: `84a786e` (`Close sparse surface transfer without negative roundoff`)
- model revision: `two-material-moving-tio2-cr-positive-row-closure-v7`
- Vast instance: `48177892`
- remote tree: `/root/petch-4b656fd`
- parent process: `33309`
- worker count: 8
- transport device: `cuda:0`, RTX 3090
- log: `/root/zhu_v7_board_84a786e.log`

Deployed file hashes matched the pushed local sources byte for byte:

```text
c634011c047697615b0ec38cc351d409bd403ba0687ffa6e9ec399ade45f4276  src/petch/surface_transfer_3d.py
7d7cb39fc854421a7b8a45223bbffe5216d922cd71a29479395ee7ec34436973  scripts/audit_zhu_npg80_moving_cr_profiles.py
02277fbc38e63b4267374ab740d873394ac9a728cbb2d0eb3ad3150bffe24281  scripts/reproduce_zhu_npg80_moving_cr_cell.py
```

## Captured-failure regression

The exact v6 production checkpoint advances under v7 on both local arm64 CPU
and remote x86_64 CUDA. Both executions produce:

- surface-transfer fingerprint
  `64dafad5621689e785800294ec3408842939d24095c7c85d23c54e1abf498d96`;
- next-surface mesh fingerprint
  `950c8f9a1176e4184020ef9f6a72d9e4ff6c61d201ba58efa4351a6478264bac`;
- zero unresolved-node reassignment.

The maximum state-remap conservation residual is `1.371e-16` on local CPU and
`2.742e-16` on remote CUDA. The former `-2.2204460492503131e-16` sparse weight
does not recur, all weights satisfy the unchanged nonnegative validator, and
the validator tolerance was not relaxed.

## Certification run

The full local repository suite at this commit passed:

```text
2229 passed, 7 skipped in 1232.94s (0:20:32)
```

The focused common-remap/feature/Oxford suite passed 149 tests before the full
suite. The v7 model revision deliberately invalidates all v6 caches because
closing through the largest coefficient changes valid stored weights and
fingerprints at roundoff scale.

## Completed v7 caches

```text
2ab98c4c82856f6113372740b584e56041e3f2cb3e73bc45b82ae0646414de61  trajectories/w080_s14.000_ion_high_tail_0p0_42692f797926cb19.json
abbdf2093b32420cfa7bb9f982876c9eb9e26240d68d999a07f29ee78bbfce53  trajectories/w080_s14.000_ion_high_tail_0p65_f9ba2d90f86e4da1.json
bcb9e0b727dcc946ab282a2131158a545dd3fe932ca9d75695c0a759d61ff75a  trajectories/w080_s14.000_ion_low_tail_0p0_5060598c59a0533b.json
147ce4b0f0511ca60f6d2661e9aebd3ceedbcc7f9f1eda1f2a0065ad446f1aac  trajectories/w080_s14.000_ion_low_tail_0p65_c4467951bc65af3f.json
aa9410a9265e2d2b098f33d57ddd44f3ef9db38b060449ff5c5e4b86b5deae58  trajectories/w080_s18.017_ion_high_tail_0p0_7a5bb429a9a4abf2.json
771fc2bccff7e49964bf09804b5f081c310ff3e36c7bf9e1f0cbacab1be89428  trajectories/w080_s18.017_ion_high_tail_0p65_ba43b0e393d3c961.json
b63d9a3255b2b884ee26c0166d40ff23953c05aa4cdddfeab462006d242090be  trajectories/w080_s18.017_ion_low_tail_0p0_3191039eb71f0d8d.json
ba7c5a0825bc7b7413dce2e42fe1af3c89b352771cf188972ec2d1d54198db0b  trajectories/w080_s18.017_ion_low_tail_0p65_ae05c672dd49a0f9.json
17cb8f90bc96f099a23f0b61212b14c1901359342e8aff526c4729cfdfd64a63  trajectories/w120_s14.000_ion_high_tail_0p0_634c0d52394b18f0.json
ee67f8e967bed4cdfc8b22ec019721f1812688221ca083bd269bf4338833c2b4  trajectories/w120_s14.000_ion_low_tail_0p0_9713fd1a430912a4.json
28326afdce64da32ba366d4bcf2f7411e8646011e49e33530019777e8ce62c45  trajectories/w120_s18.017_ion_low_tail_0p0_f7ad6305552d158c.json
```

All 11 retrieved cells declare the frozen v7 model revision and preregistration
hash, contain both rate-normalized endpoints, have exact zero particle-balance
error, and have maximum surface-state remap residual below `1.47e-15`. The
width-80 endpoints finish at approximately 683.71-684.09 nm etched depth. The
three retrieved width-120 endpoints finish at approximately 684.05-684.28 nm.
These are blind conditional results, not observed SEM matches.

Do not assemble or publish `audit.json` until all 56 v7 cache keys exist and
the production audit passes locally in check mode.
