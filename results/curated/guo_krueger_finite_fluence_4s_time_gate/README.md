# Guo/Kwon–Krüger nonlinear 4 s time gate

Status: **both twofold QMC time pairs fail the preregistered moving-profile
trajectory gate; QMC is not numerical authority**.

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

## Fine versus ultrafine result

The fresh 0.03125 s trajectory completed after the preceding receipt was
written.

| quantity at 4.0 s | 0.0625 s | 0.03125 s | relative difference |
|---:|---:|---:|---:|
| depth (nm) | 47.4659419 | 46.7541506 | +1.5224% |
| mask opening (nm) | 64.5263752 | 62.5920786 | +3.0903% |

The terminal pair and maximum depth trajectory now pass: the largest depth
difference over 0.5--4.0 s is 4.2837%. The maximum mouth difference remains
**7.4488 nm**, above the frozen 5 nm gate. Exact ledgers, radiosity balance,
and topology gates pass.

The complete 114-sample receipt is
`audit_fine_vs_ultrafine.json`. Its input audit hashes are:

- 0.0625 s: `bc2fccbc3466282c5bbe79e1b6d7e9259147765bd881891674ba24b5cec94c19`;
- 0.03125 s: `9f8d5c6240ddef1fd63c95912fec46aaf94bc4d7121ecf30d378e0bd2718f63a`.

## Consequence

Neither QMC time pair is a certified production path. More QMC time
refinement is not authorized because the independent frozen-checkpoint paired
control in `results/curated/krueger_qmc_receiver_singularity/` additionally
isolates a 149.6x inverse-receiver-area speed spike. Krueger's
translationally invariant geometry is instead being graded with the
deterministic-extruded exchange operator.

This numerical result does not change the physical evidence boundary:
Krüger's aggregate ion identity and stable C4F6 wafer flux remain unpublished.
