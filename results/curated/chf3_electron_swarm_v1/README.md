# CHF3 electron-swarm source replay

| transport deck | all-point median | all-point maximum | 40--100 Td maximum |
|---|---:|---:|---:|
| kushner_zhang_working_set | 25.01% | 82.19% | 13.42% |
| nist_evaluated_constant_join_ratio | 13.69% | 29.61% | 9.59% |
| nist_evaluated_linear_return_to_working_set_at_120eV | 13.69% | 29.61% | 9.59% |

The original working set reproduces the source paper's own low-field
overprediction. An independent BOLOS solve and petch agree on that diagnosis.
The NIST-evaluated elastic/momentum backbone materially improves transport in
the declared 40--100 Td engineering band while keeping the original inelastic
chemistry fixed.

This is deliberately a source-replay receipt. The working set was regressed
against swarm behavior, the NIST curve is not labeled bulk-versus-flux, and
the neutral-dissociation branches remain weakly constrained. It authorizes a
transport input for the next reactor rung, not a unique reactor state, wafer
flux, or feature depth.

Maximum representative-grid drift change: `0.1529%`.
