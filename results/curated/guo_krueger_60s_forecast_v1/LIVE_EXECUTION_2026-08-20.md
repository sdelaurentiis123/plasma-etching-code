# Live execution receipt — Guo/Kwon to Krueger no-fit forecast

Snapshot: 2026-08-20 15:42 EDT

This is an operational receipt, not a completed scientific result.  The
frozen authority and claim boundary are in `PREREGISTRATION.md`.

## Reproducible source

- local branch: `codex/validation-first-multiphysics`;
- source commit: `d852a1f`;
- remote source tree: `/root/petch-d852a1f`;
- source archive SHA-256:
  `153b416c236ab1493389fd050181d5afdb1279bdb719d522fa8c03eefff524e0`;
- authorized Vast instance: `48177892`;
- SSH endpoint: `ssh6.vast.ai:17892`;
- venv: `/root/petch-venv`.

The archive hash matched before extraction.  Six focused remote tests passed
before launch.

## Live processes

| Case | Supervisor PID | Child PID | Output | Log |
| --- | ---: | ---: | --- | --- |
| `nominal_unresolved` | `10284` | `10288` | `/root/krueger_guo_60s_nominal` | `/root/krueger_guo_60s_nominal.log` |
| `all_cf2` | `10285` | `10289` | `/root/krueger_guo_60s_cf2` | `/root/krueger_guo_60s_cf2.log` |
| `all_cf3` | `10286` | `10290` | `/root/krueger_guo_60s_cf3` | `/root/krueger_guo_60s_cf3.log` |

Every supervisor has a corresponding `/root/krueger_guo_60s_*.pid` file and
is configured for 1,800-second wall-budget checkpoints with bounded automatic
resume.  Use those exact PIDs for monitoring; `pgrep -f` is unsafe because it
can match the probe itself.

## First independently inspected checkpoint

At 0.09375 seconds of physical time, after six accepted steps:

| Case | Etch depth | Mask opening | Runtime target read |
| --- | ---: | ---: | --- |
| `nominal_unresolved` | 0.9830 nm | 88.9619 nm | false |
| `all_cf2` | 0.8547 nm | 88.9618 nm | false |
| `all_cf3` | 0.8278 nm | 88.9619 nm | false |

All cases reported one CFL substep, no rejected trial, no unresolved material
volume, zero material-ledger residual, and neutral-radiosity balance errors
between `1.0e-12` and `1.5e-12` at the inspected point.

The early transient ordering is not a final scientific result.  The CF2+ and
CF3+ cases are declared composition endpoints; no mixture fraction may be
selected from Krueger's measured depth.

## Runtime expectation and landing rule

Each trajectory contains 3,840 nominal 0.015625-second steps.  The observed
early shared-host rate is approximately 45 seconds per accepted step, which
would imply roughly 48 hours for a full trajectory if the cost remained
linear.  A physical `terminal_feature_clogged` event may end it earlier.

Accept a case only after it reaches 60 seconds or a declared terminal event,
then freeze and verify the hashes of `audit.json`, `checkpoint.npz`, and
`forecast_execution.json`.  Only a separate scorer may open the 825 nm target
after those runtime artifacts are sealed.
