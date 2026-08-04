# The mouth defect is an early transient, not a late-time runaway (2026-08-04)

Reproduce: `python scripts/early_transient_analysis.py`.
Inputs: archived pilot audits in `results/curated/mixed_layer_feature_v1/`
(`ml9a`, `ml13`, `ml16a`, `ml16b`) and the corrected-transport band audit
`results/curated/lip_deposition_audit/audit_neck45_dx0.01.json` (`cbbd2d6`).

## Why this pass was commissioned, and what it overturns

`RESULTS_DX5_RESOLUTION_VERDICT_2026-08-04.md` (`ab5b880`) reported that the
throat descends normally to 180 nm by t = 30 s and then *reverses* to 130 nm
while the aperture collapses — "a second, higher constriction overtaking the
descending neck" — and posed the remaining defect as a late-time
profile-evolution runaway in the t = 20–40 s window.

**Both halves of that framing are wrong.**

### 1. The "reversal" is a metric artifact

Stepping through the window step by step (not at 2 s samples), the reported
throat depth hops 160 → 170 → 110 → 180 → 120 → 180 → 130 nm across
consecutive steps while the aperture itself declines *smoothly and
monotonically* through every one of those hops. The neck is a long, nearly
flat constriction spanning roughly 110–180 nm below the mask top, and
`mask_opening_throat_z` is an argmin over near-degenerate minima. There is no
takeover event and no second constriction.

### 2. Closure decelerates monotonically — there is no runaway

Per-side closure rate in `ml16a`, by window:

| window (s) | per-side closure | etch rate | closure/etch |
|---|---|---|---|
| 1–4 | 3.391 nm/s | 21.57 nm/s | **0.1572** |
| 4–8 | 1.978 | 18.12 | **0.1091** |
| 8–12 | 0.818 | 15.56 | 0.0526 |
| 12–20 | 0.540 | 11.90 | 0.0454 |
| 20–30 | 0.308 | 9.63 | **0.0319** |
| 30–40 | 0.335 | 7.72 | 0.0433 |
| 40–50 | 0.176 | 5.50 | **0.0320** |
| 50–60 | 0.106 | 4.47 | 0.0237 |

Krüger's run-average, from his own endpoint (90 → 38.8 nm neck, 825 nm depth,
60 s): per-side 0.427 nm/s, etch 13.75 nm/s, **closure/etch = 0.0310**.

From t ≈ 20 s onward petch's closure/etch ratio *is* Krüger's ratio (0.032,
0.032). **The late-time behaviour is correct.** The excess is entirely in the
first ~10 s, where the ratio runs 5.1× / 3.5× / 1.7× his value.

## Where the closure budget is actually spent

Krüger's entire 60 s closure budget is 51.2 nm of aperture (90 → 38.8).

| by t = | petch aperture lost | % of Krüger's *full-run* budget |
|---|---|---|
| 2 s | 16.9 nm | 33 % |
| 4 s | 29.4 nm | 57 % |
| 6 s | 39.3 nm | 77 % |
| **8 s** | **45.2 nm** | **88 %** |
| 12 s | 51.8 nm | **101 %** |
| 60 s | 78.9 nm | 154 % |

petch spends Krüger's complete 60-second closure budget in **12 seconds**.
The much-quoted "38.4 nm at t = 12 s versus the experimental 39.0 nm neck" is
therefore *not* evidence of a correct early trajectory — it is the aperture
crossing his final value five times too early and continuing.

### Universal across every archived configuration

| run | t=2 | t=4 | t=8 | t=12 | t=30 | final | budget spent by t=8 |
|---|---|---|---|---|---|---|---|
| ml9a-base-atoms | 75.9 | 63.9 | 47.2 | 42.4 | 29.1 | 22.5 | 84 % |
| ml13-base-cascade | 75.9 | 63.9 | 47.2 | 42.3 | 28.0 | 24.8 | 84 % |
| ml16a-verbatim-lift | 73.1 | 60.6 | 44.8 | 38.2 | 23.5 | 11.1 | 88 % |
| ml16b-ml13c-lift | 72.8 | 60.7 | 44.6 | 38.8 | 25.0 | 12.3 | 89 % |

Reflection off/on, paper vs Table-6.5 constants, pre- and post-√2 lift: the
early transient is identical to within 5 %. It is **structural, not a
constant-set or beam question**. (The configurations differ only in the *late*
tail, which is where the √2 lift acts.)

## Both commissioned candidates are falsified by that timing

- **Source-flux competition as the trench deepens.** 88 % of the damage is
  done by t = 8 s, when the trench is ~160 nm deep out of an eventual 590 nm.
  A mechanism that grows with trench depth cannot cause a defect that is
  essentially complete before the trench is deep.
- **Mask-corner rounding over time.** `top_feature_width_nm` moves 85.7 → 76.0
  nm over the full 60 s (9.7 nm), and only 4.3 nm during the window in which
  the neck loses 45 nm. The corner is nearly static while the neck collapses.

## The term that erodes the margin: none — the top band never had margin

On the initial geometry the mask wall is vertical. Krüger's own setup is the
same, verbatim (thesis §6.4, line 4373): *"A SiO2 substrate is covered by a
850 nm thick AC film patterned to contain an ideal, straight walled, opening
with an initial width of 90 nm. The etch was performed for 60 seconds."*

On a vertical wall the ion removal channel vanishes structurally: areal flux
carries `cos θ` and the sputter yield carries it again, so removal collapses
~200× while the *isotropic* channels (depositor flux, O flux) are barely
attenuated. The audited top band bears this out exactly:

| quantity (top band, 0–50 nm, tilt 0.47°) | value |
|---|---|
| removal / deposition, measured | 0.216 |
| ...of which the O channel (geometry-free) | 0.180 |
| ...of which ions | 0.037 |

So the film grows at ≈ 0.78 × deposition from t = 0, and keeps growing until
shadowing throttles the delivered depositor flux. That is precisely the
observed shape: fast early closure, monotone deceleration, no equilibrium.
The arithmetic closes — the probe's top-band net at a 45 nm neck is
−1.77 nm/s per side (3.5 nm/s of aperture) against 8.5 nm/s observed at t < 2 s
with the aperture still at 90 nm and correspondingly less shadowed.

## The inversion, and why it is a mechanism-level statement

Deposition and O-removal are both thermal and isotropic (measured isotropy
ratio 1.0000 face by face), so delivery cancels in their ratio and the balance
is geometry-free:

```
removal/deposition  =  p_ox * J_O / (s_eff * J_dep)
```

With Krüger's converged `p_ox = 0.0423` and his own flux table
(`J_O = 7.70e20`, `J_dep = 3.094e21`):

| | |
|---|---|
| our effective lip sticking `s_eff` | 0.0586 (fresh/crosslinked blend at x_xl = 0.163) |
| `s_eff` required for balance | **0.0105** |
| published crosslinked row | 0.020 |
| removal/deposition if the lip were **fully** crosslinked | **0.526** |

**No state reachable within Krüger's published sticking rows balances the
near-vertical lip.** Even a completely crosslinked film — the lowest-sticking
state his mechanism offers — reaches only 0.53. This is consistent with his own
sensitivity statement that `p_ox = 0.005` produces *"a complete clog at the top
of the feature"*: in his model too, the top-of-feature balance is set by the O
channel and sits near the clog boundary. It also means the gap cannot be closed
by any crosslinking-state fix, which is the hypothesis this inversion was built
to test.

## Two candidate explanations tested in this pass

**Units-per-deposition-event — falsified.** The O-channel fix (`63cfefa`)
found that oxidation removed one *atom* where Krüger's row removes one
*polymer unit*, a 2.69× throttle. The deposition side was checked here for the
mirror-image error and is **correct**: `mixed_layer_mechanism.py:459-462`
accumulates `p * flux * c_atoms` and `p * flux * f_atoms`, so each sticking
event adds the species' full stoichiometry (CF2 → 1 C + 2 F), exactly one
polymer unit per collision, matching his `AC(s) + CF2 → AC(s) + CF2(s)` row.
No factor lives here.

**Absolute rate — the discrepancy is entirely the lip film's crosslink state.**
From the published rows and flux table alone (no transport, no geometry):

| unshadowed blanket rate | value |
|---|---|
| gross deposition, **fresh** polymer rows (0.1 / C2F3 0.03) | 6.74 nm/s |
| gross deposition, **crosslinked** row (0.02) | 1.98 nm/s |
| O removal (geometry-free) | 1.17 nm/s |

On a vertical wall, where ion removal vanishes, net film growth per side is
therefore **5.57 nm/s if the film is fresh** and **0.81 nm/s if it is
crosslinked** — a 6.9× lever. Scaled by the audited top-band delivery (0.372):

| lip film state | predicted per-side closure |
|---|---|
| fresh (our measured x_xl = 0.163 → mostly fresh) | ≈ 2.07 nm/s |
| fully crosslinked | ≈ 0.30 nm/s |
| **Krüger's run-average (his endpoint)** | **0.427 nm/s** |

petch's measured early per-side rate is 3.39 nm/s (t = 1–4 s, aperture still
near 90 nm and less shadowed than the 0.372 probe geometry) — consistent with
a mostly-fresh lip. **A fully crosslinked lip would reproduce Krüger's closure
rate to within a factor of 1.4.** The entire 5× defect is spanned by this one
state variable.

## What this leaves: one sharply-posed question with a preregistered target

> Our lip film sits at crosslinked fraction **x_xl = 0.163**. Matching Krüger's
> closure requires **x_xl ≳ 0.9**. Both models drive crosslinking with the same
> ion-dose row and both see the same collapsed ion flux on a vertical wall
> (cos 88.7° ≈ 0.023). Why is his lip film crosslinked and ours not?

The mechanism is self-consistent in our model — fast fresh deposition outruns
the ion dose available at grazing incidence, so the film never converts — which
is exactly why it is worth auditing rather than assuming: the crosslinking
kernel is the one channel on the lip that has *never* been dimensionally
audited, and it is the only remaining term with the required leverage.

Next steps, in order, all local and free:

1. **Crosslinking kernel — started here, no unit error, but a live
   inconsistency found.** `kernel_xl = flux * max(E_ion − E_interface, 0)` and
   `xl_rate = kernel_xl / 25 eV * (1 − x_xl)`
   (`mixed_layer.py:273, 525`): one *atom* converted per 25 eV of
   film-absorbed energy, and `n_xl_film` is counted in atoms against
   `film_total` in atoms — dimensionally self-consistent, no atom-vs-unit
   throttle of the kind the O channel had.

   But the 0-D estimate does not reconcile with the measured steady state. The
   lip film is thick (15.6 nm), so `E_interface → 0` and each ion deposits
   essentially its full energy in the film: at keV energies that is ~10²
   atoms crosslinked per ion, which even against the collapsed grazing ion
   flux should out-run the ~1.9e20 atoms/m²/s of deposition and drive
   `x_xl → O(1)`. The audit measures **x_xl = 0.163**. Either the delivered
   lip ion flux is far below the cos-projection estimate, or a suppression
   term (the `(1 − x_xl)` factor, `_overdraw_scale`, the `xl_share` draw-down,
   or the de-crosslink channel) is dominating. **Resolving that one number is
   now the top of the queue** — it is the only remaining term with the
   required 5× leverage, and it is measurable per-face in the probe.
   Also check whether the grazing slant path (a grazing ion traverses a longer
   path *in the film* and should deposit **more** energy per ion, not less) is
   represented, since that is precisely the lip geometry.
2. **Sweep x_xl at the lip in the probe** (frozen 45 nm geometry, corrected
   transport) and confirm the predicted 6.9× closure-rate lever end to end.
3. Only then a run — and it needs **12 s, not 60 s**, since 88 % of the defect
   is complete by t = 8 s. That makes the confirmation ~5× cheaper than the
   standing spec, graded on `closure/etch` by window against 0.0310.

## Method note

Grading a closure defect by its *endpoint* hid a 5× rate error behind a
plausible-looking intermediate value (38.4 nm at t = 12 s). Trajectory-level
grading against a run-average ratio exposed it immediately. Future mouth work
should grade `closure/etch` by window, not aperture at 60 s.
