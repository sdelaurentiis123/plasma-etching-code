# Guo/Kwon–Krüger nonlinear 4 s time gate

Status: **0.125 s nominal step failed the preregistered trajectory gate;
0.0625 s and 0.03125 s comparison in progress**.

The criteria in `PREREGISTRATION.md` were committed and pushed before the
0.0625 s trajectory reached the one-cell mouth transition.

## Medium versus fine result

| quantity at 4.0 s | 0.125 s | 0.0625 s | relative difference |
|---|---:|---:|---:|
| depth (nm) | 47.6387010 | 47.4659419 | +0.3640% |
| mask opening (nm) | 64.0491885 | 64.5263752 | -0.7395% |

The terminal pair agrees well, but two precommitted maximum-trajectory gates
fail:

- maximum absolute depth difference over 0.5–4.0 s: **5.4233%**, above 5%;
- maximum mask-opening difference: **8.4949 nm**, above 5 nm.

The latter is a phase shift in a one-row, 10 nm throat kink. Neighboring
aperture rows differ by roughly 0–2 nm, the kink recovers in both trajectories,
and the terminal profiles agree. It is nevertheless a real accepted-geometry
time sensitivity and remains a failed gate rather than being relabeled after
inspection.

All other gates pass: terminal depth and mouth, exact material ledgers,
neutral-radiosity balance (`2.167e-12` maximum), and zero topology events.

`audit_medium_vs_fine.json` contains all 58 matched samples. The complete
0.0625 s artifact was copied before any horizon extension to
`/private/tmp/krueger_guo_transient_dt0625_dx10_4s_receipt`.

## Consequence

The running 0.125 s, 60 s trajectory is a scout, not the certified endpoint.
A fresh 0.03125 s trajectory is testing the 0.0625 s path against the same
unchanged gates. If that rung passes, the 0.0625 s checkpoint is the numerical
production path. If it fails, no endpoint depth will be promoted until the
time-sensitive interface/remap dynamics are fixed or a still-finer rung
passes.

This numerical result does not change the physical evidence boundary:
Krüger's aggregate ion identity and stable C4F6 wafer flux remain unpublished.
