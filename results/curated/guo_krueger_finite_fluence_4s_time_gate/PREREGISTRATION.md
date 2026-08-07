# Guo/Kwon–Krüger nonlinear 4 s time gate

Frozen before the 0.0625 s trajectory reached the 2.9 s mouth transition.

The existing 60 s trajectory uses a 0.125 s nominal step with the adaptive
displacement controller. A cloned 0.0625 s checkpoint is being advanced to
4.0 s on the same 10 nm mesh and identical physical operator. The comparison
begins at 0.5 s, after both trajectories share the already certified prefix.

The 0.125 s trajectory is authorized numerically for the endpoint only if all
of these gates pass:

- terminal depth absolute relative difference at 4.0 s: at most 5%;
- terminal mask-opening absolute relative difference: at most 5%;
- maximum matched depth absolute relative difference over 0.5–4.0 s: at most
  5%;
- maximum matched mask-opening difference: at most 5 nm;
- no accepted topology event in either prefix;
- exact material ledgers;
- maximum neutral-radiosity relative balance residual at most `1e-9`.

Medium values are sampled at the fine accepted times by piecewise-linear
interpolation between accepted medium checkpoints. This interpolation is a
comparison operator only; no interpolated state is used to continue either
trajectory.

The script is
`scripts/audit_guo_krueger_matched_time_gate.py`. Passing this gate can certify
the time integration through the nonlinear prefix. It cannot identify
Krüger's unpublished ion mixture or stable C4F6 boundary and cannot promote
the 60 s result to Tier-A prediction.
