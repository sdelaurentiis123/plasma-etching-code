# The two ion channels, measured against Gray's absolute yields

The joint solve (`RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`) graded RATIOS —
dynamic range, peak/normal, depth factor.  Gray 1993 also publishes the
ABSOLUTE yields, which grade each channel's magnitude separately.

Source: Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), 350 eV Ar+ on
SiO2, replotted as Kwon (ScD, MIT DMSE 2004) Fig. 3.4 p. 76 — floor **0.28**
(F/Ar+ -> 0, pure physical sputter), plateau **1.10** (F-saturated).

| quantity | petch | Gray 1993 | petch / measured |
|---|---|---|---|
| floor, F/Ar+ -> 0 (bare row alone) | 0.341 | **0.28** | **1.22x too strong** |
| plateau, F-saturated (complex row) | 0.390 | **1.10** | **0.35x — 2.8x too weak** |

The floor is exactly the bare row evaluated at 350 eV:
`0.0852 * (350-70)/(140-70) = 0.341`, so this is a direct read of that row's
magnitude, with no coverage or transport interpretation in between.

## Why this matters

**The two channels are wrong in opposite directions**, which is why every
single-knob attempt failed and why the model is stuck ion-limited:

- The **physical sputter** row over-delivers by **1.22x**.  That is the fourth
  independent route to the same term: the cascade audit put the coupled floor
  rate ~1.4x high; Krueger states the overestimate himself (thesis L4884-4888,
  *"the effect of ion energy ... might be overestimated in the mechanism"*); the
  joint solve's beam-selected assignment undershoots the depth gate by 1.24-1.35x;
  and now the bare row reads 1.22x against a published absolute yield.
- The **chemically-enhanced** row under-delivers by **2.8x**.  This is the
  regime defect stated as a magnitude: with the chemical channel that weak, it
  can never dominate the fluorine-free sputter channel, so coverage collapses
  12x under starvation while the rate moves 0.7% (`RESULTS_LIMITING_REGIME` §2).
  Neutral-limited behaviour is not reachable by re-weighting; the chemical
  channel simply is not strong enough to be the limiting one.

Together they explain the dynamic range miss exactly: measured 0.255, petch
0.873, because our floor is high and our plateau is low from both sides.

## Status: a receipt, not yet a change

These are two absolute measurements at ONE energy (350 eV) against a feature
that runs at keV, so they cannot be applied as magnitude corrections until the
energy scaling that carries them from 350 eV to 3406 eV is settled — and that
scaling is exactly the open ZBL-vs-published-linear question
(`RESULTS_LIMITING_REGIME` §3).  Landing a 350 eV correction under an unsettled
keV scaling would be fitting at the reference point and hoping.

What they do settle is the SHAPE of the departure the depth gate needs.  The
joint solve showed no sourced combination reaches it; this shows why, per
channel, with a published number on each.  The correction is a beam
measurement, not a tuned constant — which makes departing from Appendix B's
magnitudes a receipted move rather than a knob.
