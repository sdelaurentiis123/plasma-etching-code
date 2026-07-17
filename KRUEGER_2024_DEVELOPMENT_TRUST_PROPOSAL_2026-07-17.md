# Krueger 2024 bounded development proposal

Date: 2026-07-17
Status: development decision artifact only; held-out profiles remain sealed

## Why this exists

The completed 5 nm trajectory is valuable but not authoritative: the periodic remap repair entered
at `56.482184 s`, and gas-cavity continuation entered at `56.920695 s`.  It therefore supplies a
base-case residual and a calibration direction, not a parameter freeze or validation result.

The only response model reused here is the checksum-bound 10 nm base-only Jacobian.  The two inputs
are intentionally kept separate:

```text
old 10 nm base-only response matrix       completed mixed-operator 5 nm base residual
                  \                         /
                   \                       /
                    full local secant step
                              |
                    explicit trust safeguard
                              |
                   one development proposal
                              |
             actual-versus-predicted base comparison
```

No oxygen-ratio profile, power-sweep profile, or transfer observation enters this calculation.

## Current response and full local step

The evaluated parameter pair is

```text
crosslinked-film fraction = 0.8934059741411972
oxide-yield scale         = 0.5667632723491973
```

Its complete development endpoint is `39.466366780154374 nm` opening and
`900.8570308260097 nm` depth, compared with calibration targets `45 nm` and `825 nm`.
Using the archived response matrix (condition number `12.8484`) gives the full local step

```text
delta fraction    = +0.04215303421546874
delta yield scale = -0.04840489182999453
```

and the full, unclipped candidate

```text
fraction    = 0.9355590083566659
yield scale = 0.5183583805192028
```

The linear model necessarily predicts the two targets exactly at that full step.  That algebra is
not evidence that the nonlinear, fine-grid engine will do so.

## Trust-radius choice

R1.9 and WP5 require a bounded local step and one proposal at a time, but they do **not** declare a
numeric radius. R1.9 authorizes exactly this fixed 10 nm model check and closes the 10 nm sequence
afterward. A radius must therefore be identified as a numerical safeguard, never as sourced physics.
The implemented initial trust box uses the last proposal that was actually evaluated:

```text
coordinate scale in fraction    = 0.027611078244628235
coordinate scale in yield scale = 0.008114305762722413
normalized L-infinity radius    = 1
```

The full step has normalized coordinates `(1.52667, 5.96538)`.  It is therefore scaled uniformly by
`0.16763400259670264`, preserving its direction.  The one safeguarded development candidate is

```text
fraction    = 0.9004722559883319
yield scale = 0.5586489665864749
```

The frozen response predicts `40.39399186569918 nm` opening and `888.1408131235443 nm` depth.  This
small first move is intentionally information-seeking: it tests whether the old 10 nm response
direction transfers to the repaired 5 nm operator before another larger move is authorized.

This radius is conservative, not uniquely ordained.  If the actual-versus-predicted comparison is
good, the recorded trust policy may grow it.  If it is poor, the radius shrinks or the response is
rebuilt.  There is no held-out-dependent branch.

## Evidence binding and execution gate

Generator: `scripts/krueger_2024_development_trust_proposal.py`

Generated artifact:
`results/krueger_2024_base_calibration_r18/development_trust_region_proposal.json`

The artifact binds:

- the old response-model file and its embedded proposal checksum;
- the exact completed endpoint and continuation receipt;
- the base calibration-target table;
- the current R1.9 protocol and governing campaign;
- both the full and safeguarded parameter steps;
- `authority=false` and `held_out_profile_data_read=false`.

It authorizes no long run by itself.  WP3/WP4 numerical gates still precede the one development
evaluation.  Whatever the outcome, an eventual final pair must be confirmed from `t=0` by one clean
uniform-5-nm or certified-AMR-equivalent operator before held-out profiles are revealed.
