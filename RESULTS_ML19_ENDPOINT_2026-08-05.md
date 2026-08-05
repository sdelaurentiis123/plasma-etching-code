# ml19 — 60 s endpoint with deposition-crosslinking (stopped t = 46.2 s)

Configuration: identical to ml18 (campaign 7) at `--duration-s 60`, dx = 10 nm,
mixed layer, `literature_v1` grazing reflection, deposition-driven crosslinking
(`3a931b1`), corrected axisymmetric lift (`6e97ef3`), O-channel normalisation
(`63cfefa`). Every mechanism in the stack is source-verbatim,
measurement-gated, or declared-open. Audit + checkpoint + step log:
`results/curated/mixed_layer_feature_v1/ml19-depxl-60s/`.

**The run stopped at t = 46.231 s** on a surface-topology refusal, not a
physics failure — see "Robustness stop" below. 77 % of the closure budget and
the entire mouth-equilibrium window are inside the recorded history.

## Headline: the mouth equilibrates for the first time

| t (s) | 34.1 | 37.2 | 40.4 | 43.6 | 46.2 |
|---|---|---|---|---|---|
| mask aperture (nm) | 50.58 | 51.30 | 50.92 | 50.84 | 50.37 |

Over t = 34–46.2 s (83 records): mean **50.92 nm**, spread 2.52 nm, drift
**−16.8 pm/s**. At that drift the aperture would reach t = 60 s at ≈ 50.1 nm.

Every prior run sealed monotonically (ml16a 11.1, ml16b 12.3, ml15 7.2). This
is the first configuration in the campaign whose mouth reaches a stable
equilibrium — the qualitative behaviour Krüger's feature has and ours never
did. The mechanism is the deposition-crosslinking correction: the lip film
hardens as it deposits (his §2.2.3), so growth self-limits instead of running
away.

## Gates

| gate | target | ml19 | verdict |
|---|---|---|---|
| mask constriction CD | 45 (his `w_m`) | **50.9** (equilibrium) | +13 %, near-miss |
| — vs his simulated neck | 38.8 | 50.9 | +31 % |
| constriction depth below mask top | 200 (SEM) / 271 (sim) | **170.2** | −15 % vs SEM |
| mask remaining | 850 | **850.2** | **PASS** (exact) |
| etch depth | 825 ± 5 % | 835.9 **at t = 46.2 s** | see below — **MISS** |
| closure/etch, 1–4 s | 0.0310 | 0.0566 (1.83×) | above band |
| closure/etch, 4–8 s | 0.0310 | 0.0490 (1.58×) | above band |
| closure/etch, 8–12 s | 0.0310 | 0.0257 (0.83×) | in band |
| closure/etch, 12–20 s | 0.0310 | 0.0153 (0.49×) | in band |
| closure/etch, 20–40 s | 0.0310 | 0.0160 (0.52×) | in band |
| closure/etch, 40–46 s | 0.0310 | 0.0037 (0.12×) | in band |

`neck_cd_nm` = 25.58 at 1680 nm below the mask top is the advancing etch-front
taper, not a mask neck (same caveat as ml18); the mask throat is the graded
object.

### The depth gate is a MISS, and the "PASS" above is an artifact

835.9 nm at t = 46.2 s sits inside 825 ± 5 %, but the experiment reaches 825 nm
at t = 60 s. Late etch rate is **16.74 nm/s against Krüger's 13.75 average
(1.22×)**; linear extrapolation gives **1066 nm at 60 s, +29 %**.

This inverts the previous failure and is the most useful new finding.
ml16a undershot depth (590 vs 825) *because* its sealed mouth throttled the
flux. With the mouth held open at ~51 nm the floor now over-etches. **Two
errors were compensating**: a mouth that closed too fast and a floor that
etches too fast. The first is fixed; the second is now exposed and is the next
target. It is a floor-side question (oxide removal per ion at normal
incidence), independent of everything audited at the lip.

## Trajectory vs the ml16a baseline

| t (s) | ml19 aperture | ml16a | ml19 depth | ml16a |
|---|---|---|---|---|
| 4 | 78.49 | 60.61 | 88.2 | 87.4 |
| 8 | 70.35 | 44.79 | 171.3 | 159.9 |
| 12 | 66.45 | 38.25 | 247.3 | 222.1 |
| 20 | 62.00 | 29.61 | 392.7 | 317.3 |
| 30 | 52.54 | 23.46 | 568.4 | 413.6 |
| 40 | 51.15 | 16.77 | 732.3 | 490.8 |

Closure budget spent by t = 46.2 s: **39.6 nm = 77 %** of Krüger's full-run
budget, against the baseline's 101 % by t = 12 s.

## Robustness stop (follow-up item, not a physics result)

```
petch.feature_step_3d.SurfaceTopologyChangeError: surface topology changed under
periodic_xy_component_cavity_breakthrough_v1 ...
periodic material-component sizes changed from ((1, (3444,)), (2, (940,)))
to ((1, (3440, 2)), (2, (940,)))
```

A **2-cell fragment detached from the mask** (component 1: 3444 → 3440 + 2).
The policy `continue_gas_cavity` covers gas-cavity topology changes, not
material-component splits, so the step is refused; the adaptive controller then
shrank dt 0.1625 → 7.6e-05 s over six steps with simulated time frozen, which
is a terminal stall rather than a recovery. Two-cell slivers are a
remeshing/level-set artifact, not physics — the mask is armoured and cannot
genuinely fragment. Needed: either a sliver-reabsorption rule (fragments below
a declared cell count merge back into their parent component) or an explicit
topology-event path for material splits. Until then, 60 s runs at dx = 10 nm
are not reliably completable.

## Standing

The mouth arc is closed at the mechanism level: the defect was the crosslink
inversion, the fix is verbatim from the source, and the mouth now equilibrates
at 50.9 nm against a 45 nm target with a residual bounded by the untabulated
bond-multiplicity integer (`[VERIFY]`; m = 3 → ≈ 1.4×, m = 8 → ≈ 1.0× on the
first-window ratio). Mask survival is exact. The open item is the floor etch
rate at 1.22×, newly visible now that the mouth no longer throttles it.
