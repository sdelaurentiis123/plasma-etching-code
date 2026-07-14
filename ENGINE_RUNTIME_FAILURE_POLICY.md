# Unified engine runtime-failure policy

Revision: 2026-07-14-r1

Scope: production 3-D charged transport, surface response, physical-time charging, and charging
co-evolution. Configuration validation is distinguished from failure during a valid run.

## Product contract

An unattended run is not promised to be literally unstoppable. It is promised to recover from
every **predictable, bounded work-limit condition** without human intervention, to checkpoint
often enough that a process or host loss does not erase the campaign, and to stop rather than
invent physics when an invariant or declared validity domain is violated.

No recovery may change the physical operator, sampling epoch, accepted state, or convergence
gate. Every retry and every bounded approximation is recorded. Exact signed charge conservation
remains a per-step invariant; therefore Russian-roulette particle termination is not an allowed
default in the B4-certified charging path because it conserves charge only in expectation.

## Classes

| Class | Meaning | Engine response |
| --- | --- | --- |
| P — preflight refusal | Invalid, incomplete, or contradictory declared inputs | Refuse before state evolution; no checkpoint is needed |
| A — inline recovery | A correct result exists after more bounded numerical work | Extend/refine inside the same physical evaluation; keep the same state and sample epoch; count the recovery |
| B — priced closure | A declared deterministic approximation has a rigorous error bound | Continue; add the bound to the run error budget and final audit |
| C — hard integrity stop | Continuing could violate conservation, geometry, state integrity, or model validity | Atomically retain the last certified state and stop |
| D — operational recovery | Process, host, or storage interruption outside the physical operator | Resume the last atomic checkpoint with the exact configuration and epoch |

## Charging/transport classification

| Condition | Class | Current mechanism | Required evidence |
| --- | --- | --- | --- |
| Initial charged-cascade bounce budget exhausted while explicit lineage remains | A | Budget doubles inline up to the declared emergency ceiling | Initial/final budget, extension count, derived bound, exact charge ledger |
| Reflection model supplies uniform absolute-charge contraction `rho < 1` | A | Minimum sufficient response budget is derived from `rho` and tail tolerance; for `rho=0.95`, tolerance `1e-10`, the budget is 450 response evaluations | Model provenance, `rho`, tolerance, derived budget |
| Cascade remainder below declared relative absolute-charge tolerance | B | Deterministic absorption on the currently reached faces | Exact signed charge; remaining absolute-charge fraction; spatial-current relative L1 bound `<= 2*tail_fraction` |
| Cascade remains explicit at the emergency ceiling, or response has no subunit contraction and does not close | C | Refuse the evaluation | Full unresolved signed/absolute charge and lineage; replayable pre-step state |
| Float32 shared-edge or gas-side hit fails geometric certification | A | Replay only invalid lineages in float64; halve flight timestep up to the bounded ladder while preserving physical horizon | Replay eligible/count/fraction; edge-inset count; final gas-side certification |
| Float64 replay still cannot certify a path after the ladder | C | Refuse; do not soften visibility or discard the path | Exact ray/face/origin/velocity diagnostics and last certified state |
| Primary, adjoint, or re-impact trajectory exhausts its declared physical horizon | A target; C today | Current engine refuses incomplete transport. The next migration is deterministic horizon extension at fixed `dt`, bounded by a declared emergency horizon; diagnostic truncation remains non-production only | Truncated rate/weight, initial/final horizon, extension count, endpoint agreement with a larger strict horizon |
| Unknown trajectory termination, zero impact speed, inward surface emission, or inconsistent event lineage | C | Refuse | Offending event lineage and state checkpoint |
| Local surface-response charge/energy ledger fails | C | Refuse | Local and global signed inventories; offending response provenance |
| Step charge deposition, face-to-node projection, or remap ledger exceeds roundoff contract | C | Refuse | Both inventories, residual, scale, mesh/config hashes |
| Fresh-scramble SER requested | P | Refuse before evolution | Sampling and timestep policy in manifest |
| Frozen estimator method map differs from separately certified pilot | P/C | Refuse before evolution when known; hard-stop if discovered during evaluation | Pilot map checksum and sampling provenance |
| Safeguarded deterministic SER exhausts minimum timestep | C for that schedule | Retain state and stop; this does not prove equilibrium nonexistence | Complete attempted-step history and physical-time reference comparison |
| Heartbeat/checkpoint write fails | C operational | Stop with the in-memory certified state rather than continue without promised durability | I/O error and last successful checkpoint hash |
| Worker/process dies but host and atomic checkpoint survive | D | Supervisor resumes the progress checkpoint; no sample epoch is skipped or reused incorrectly | Config/source hashes, checkpoint hash, old/new process IDs, restart count |
| Host loss/preemption | D | External instance monitor restarts on another accuracy-certified device from the last synchronized checkpoint | Device parity gate, synchronized checkpoint hash, restart record |
| Non-finite/corrupt checkpoint or checksum mismatch | C | Refuse | Expected and observed hashes; no fallback to an older state unless explicitly selected and logged |

## Error-budget ledger

Every C3 summary carries separate, non-interchangeable entries:

- exact signed-charge conservation residual;
- cascade tail absolute-charge fraction and spatial-current L1 bound;
- float64 lineage replay count/fraction and edge-inset count;
- initial/final/derived/emergency bounce budgets and extension count;
- trajectory-horizon extension/truncation counts once that A-class migration lands;
- timestep, sample, grid, and CPU/GPU refinement evidence;
- retained per-node RMS/worst diagnostics and B1/B2 gate values.

Bounded numerical work is not a physics change. A run may claim the exact hard-visibility
operator only when every B-class bound is below its declared tolerance and the independent B5
audit uses the same operator. A C-class event can never be converted to success by relabeling it
as uncertainty.

## Operational lifecycle

1. Preflight validates model provenance, mesh/checkpoint/method-map hashes, sampling policy,
   numerical stack, and device parity.
2. Each certified evaluation invokes the progress hook. Heartbeats are atomically replaced and
   face-authoritative checkpoints are written at the declared accepted-step cadence.
3. Inline A-class work happens before the physical update is accepted.
4. A segment summary is written only after its checkpoint and current audit are durable; the
   summary is the segment commit marker.
5. The supervisor advances only from that committed checkpoint. It may automate legacy fixed-cap
   recovery, but the production engine handles cascade extension inline.
6. Completion still means the B gates plus independent high-sample exact-operator audit—not
   merely reaching a time, iteration, or work ceiling.

## Explicit non-goals

- no silent particle deletion;
- no conversion of truncation into escape;
- no soft visibility fallback;
- no stochastic charge conservation in place of exact B4 accounting;
- no automatic addition of SEE, leakage, or conduction to make a numerical failure disappear;
- no claim that a long-running engine should ignore corrupted state or failed invariants.
