# The ion-channel model space, solved jointly

Every prior pass changed one assignment and graded it against one observable.
That is why each fix over-corrected another channel — `RESULTS_ANGULAR_CONVENTION`
§5.2 put it exactly: *"any class-1 normalisation moves both at once, in opposite
senses, which is why the change that fixes one over-corrects the other. They
cannot be graded separately."*

This pass grades them together.  `scripts/ion_channel_model_solve.py` enumerates
the whole discrete, sourced space and scores every combination against every
measured constraint at once.  Nothing is fitted: each option carries a citation,
and combinations that would apply a measurement outside the material it was
taken on are not enumerated at all.

## 1. The space (citation-gated)

| axis | options | admissibility |
|---|---|---|
| polymer class-1 shape | Kress B=9.3 `f(0)=1`; Kress B=9.3 peak-normalised; Barklund yield reading; Barklund raw-rate reading | Kress is Krueger's own class-1 citation; Barklund & Blom measured an **FC film**, so it is admissible here and nowhere else |
| oxide/mask class-1 shape | Kress B=9.3 `f(0)=1`; Kress B=9.3 peak-normalised; B=1.7 `f(0)=1`; B=1.7 peak-normalised | Cho 2000 / Schaepkens 1998 measured **SiO2 in fluorocarbon**, so B=1.7 is admissible here and nowhere else |
| complex-channel energy form | ZBL `eps(E)/eps(140)`; published `(E-35)/105`, n=1 | ZBL is the K24-DEKNOB result that retired the fitted yield-scale knob; the linear form is Krueger Appendix B row `0.1471 35 1 140 2` |
| thermal-F sticking | 1.0; 0.0 | 1.0 is petch's standing Langmuir sticking; 0.0 follows Krueger having **no** thermal-F-on-bare-oxide row (L5905-5909; thermal F only fluorinates an already-complexed site at 0.1, L6548-6555) |

4 x 4 x 2 x 2 = **64 combinations**.

Two hooks were added to `mixed_layer.py` to make the space patchable —
`_complex_energy_factor` and `_THERMAL_F_STICKING`.  Both are pure refactors:
the defaults reproduce the previous expressions exactly, and the mixed-layer
gates pass unchanged.

## 2. The constraints

| id | quantity | band | source |
|---|---|---|---|
| C1a | Gray beam dynamic range `Y(F/Ar+ -> 0)/Y(sat)` | 0.20-0.30 | Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), 350 eV Ar+ on SiO2; replotted Kwon (ScD, MIT 2004) Fig. 3.4 p.76 — floor 0.28, plateau 1.10 |
| C1b | Gray half-rise `F/Ar+` | 19-35 | same source, half-rise ~27 |
| C2 | oxide total angular peak/normal | 1.28-1.36 | Cho, *JVST A* **18**, 2705 (2000) ~1.30; Schaepkens, *JVST A* **16**, 3281 (1998) ~1.33 |
| C3 | FC-film angular peak/normal | 1.30-3.50 | Barklund & Blom, *JVST A* **10**, 1212 (1992): 1.448 at 65 deg (yield reading) to 2.70 (raw-rate reading); the band spans both because the flux convention is `[VERIFY]` |
| C4 | coupled blanket rate at Krueger Table 6.1 fluxes | 4.5-14.0 nm/s | his 825 nm / 60 s = 13.75 nm/s at 1-3x funnelled delivery |
| C5 | coupled depth factor at the measured 3406 eV front | 0.735-0.812 | the [784, 866] nm gate against the +29% baseline |
| C6 | lip film growth | 0.30-0.60 nm/s | 0.7-1.4x Krueger's 0.427 nm/s per-side closure |
| C7 | ARDE sign, rate(AR 16)/rate(AR 0) | < 1 | must fall; real HARC loses ~80% by AR 40 (Huang L5430-5478) |
| C8 | film thickness at the front condition | <= 0.15 nm | validated ~0.085 nm; 4x thickening throttles interface energy |

C5 is a *ratio to the current tree*, so it measures the change each combination
makes rather than an absolute the 0-D reduction cannot claim.  Every coupled
column is evaluated in **coupled mode** at Krueger's Table 6.1 fluxes — the
lesson of `RESULTS_ANGULAR_CONVENTION` §4, where a beam-mode forecast said 0.78
and the coupled truth was 0.23.

## 3. C1b is unsatisfiable, and it names a missing constant

Before the sweep, the half-rise was traced analytically
(`scripts/gray_half_rise_scan.py`).  It is a function of the thermal-F sticking
magnitude **alone**:

| thermal-F sticking `s` | dynamic range | half-rise `F/Ar+` |
|---|---|---|
| 1.00 (petch) | 0.873 | **1.94** |
| 0.30 | 0.874 | 5.88 |
| 0.10 | 0.876 | 17.87 |
| **0.06** | 0.878 | **25.88** |
| 0.03 | 0.883 | 54.30 |
| Gray 1993 | **0.20-0.30** | **27 +/- 8** |

Two things fall out, and they are independent:

1. **The half-rise position is set by `s` and nothing else in this space.**
   Gray's measured knee implies `s ~ 0.06` on bare SiO2 — a factor **~17 below**
   petch's 1.0, and not reachable by either enumerated option (1.0 or 0.0, the
   latter removing the rise altogether).  So the survivor set is empty *by
   construction*, and the constraint that empties it names its own fix: the
   site-limited adsorption coefficient of Kwon/Sawin element **E1**, which is a
   constant petch does not carry.  It is **not** invented here — `s ~ 0.06` is
   what inverting Gray's published curve yields, and it is recorded `[VERIFY]`
   against the Gray/Kwon body until the number is read directly rather than
   inverted.
2. **`s` does not fix the dynamic range** (0.873 -> 0.878 across a 33x sweep).
   C1a and C1b are controlled by different terms: the knee position by the
   adsorption coefficient, the floor-to-plateau ratio by the bare/complex row
   balance.  Any pass that tried to fix both with one change was mis-specified.

C1b is therefore lifted out of the per-combination sweep and reported here, so
the remaining eight constraints can discriminate the space they actually span.

## 4. The matrix

**Tier 0 — free shape test (C3).**  The polymer peak/normal ratio costs nothing
to evaluate and is scale-invariant, so peak-normalisation cannot move it:

| polymer option | peak/normal | C3 (1.30-3.50) |
|---|---|---|
| Kress B=9.3, `f(0)=1` | 4.172 | FAIL |
| Kress B=9.3, peak-normalised | 4.172 (identical — ratio is scale-invariant) | FAIL |
| **Barklund yield reading** | **1.448** | **PASS** |
| Barklund raw-rate reading | 5.131 | FAIL |

**48 of 64 combinations are eliminated here, at zero cost, on a measurement.**
The measured FC-film curve selects the Barklund yield reading uniquely, which
also settles the `[VERIFY]` flux-convention question in the only way consistent
with the data: the raw-rate reading is *more* peaked than the Kress form it was
meant to bound.

**The axes are separable**, which the tier-2 rows then confirm, so the remaining
matrix is reported as axis scans rather than a 16-fold product (partial product
rows preserved in `results/curated/ion_channel_solve/sweep_partial.log`; they
agree row-for-row with the scans below):

| oxide assignment | C1a dyn. range (0.20-0.30) | C2 peak/normal (1.28-1.36) | C5 depth factor, ZBL / linear (0.735-0.812) |
|---|---|---|---|
| Kress B=9.3, `f(0)=1` | 0.873 FAIL | 1.257 FAIL | 1.000 / 1.087 FAIL |
| **Kress B=9.3, peak-normalised** | **0.210 PASS** | 1.176 FAIL | 0.577 / 0.628 FAIL |
| B=1.7, `f(0)=1` | 0.873 FAIL | 1.161 FAIL | 1.000 / 1.087 FAIL |
| B=1.7, peak-normalised | 0.667 FAIL | 1.158 FAIL | 0.868 / 0.945 FAIL |

| polymer assignment | C6 lip growth (0.30-0.60 nm/s) |
|---|---|
| Kress B=9.3, `f(0)=1` | 0.673 FAIL |
| Kress B=9.3, peak-normalised | 0.673 FAIL |
| Barklund yield | 0.673 FAIL |
| Barklund raw-rate | 0.670 FAIL |

Best row: **`barklund_yield | kress9.3_peaknorm | zbl | unity`**, 3 failures
(C2, C5, C6).  Every other combination is worse.  Thermal-F sticking `zero`
additionally breaks C1a (no F uptake, so the yield curve is flat), which is the
sweep's own check on that axis.

## 5. Verdict: the survivor set is empty, and four constraints say why

| constraint | satisfiable? | what the empty cell names |
|---|---|---|
| C1a Gray dynamic range | **YES** — uniquely by oxide peak-normalisation (0.210) | the reading Huang L2290-2296 states in words is the one the beam data selects |
| C3 FC-film shape | **YES** — uniquely by the Barklund yield reading | settles the flux-convention `[VERIFY]` |
| C1b Gray half-rise | **NO** — needs `s ~ 0.06`, space holds {1.0, 0.0} | the site-limited adsorption coefficient (Kwon/Sawin E1), a constant petch does not carry |
| C2 oxide peak/normal | **NO** — best 1.257 against a 1.28 floor (2% short) | definition-sensitive; see the caveat below |
| C5 depth factor | **NO** — required 0.735-0.812 falls in the **gap** between adjacent published options (0.628 \| 0.868) | the depth defect is not a model-selection question in this space |
| C6 lip growth | **NO** — 0.670-0.673 across the *entire* space | the lip is deposition-dominated (removal is ~200x collapsed there), so its lever is the crosslink channel, not the ion channel |

The central result is the pair C1a and C5.  **The assignment the measurements
select is the one the depth gate rejects.**  Oxide peak-normalisation is the
only reading that reproduces Gray's measured dynamic range — and it lands the
depth factor at 0.577-0.628 where the gate needs 0.735-0.812, i.e. it removes
**1.24-1.35x less** than required.

That factor is now the third independent route to the same number.  The cascade
audit put the coupled floor rate ~1.4x high; Krueger states it himself
(thesis L4884-4888) — *"the effect of ion energy (for example in sputter yield
or related processes) might be overestimated in the mechanism"*; and the joint
solve now shows the beam-selected assignment undershoots the depth gate by the
same 1.2-1.4x.  Three routes, one term: the ion-energy channel magnitude, which
petch inherits from Appendix B row-for-row.

**Caveat on C2, recorded rather than smoothed.**  This pass measures "oxide
total angular peak/normal" at the Gray beam condition (350 eV, saturated F,
all channels), giving 1.16-1.26.  An earlier pass reported 1.34 for the same
assignment under a different condition.  The quantity is definition-sensitive
(coverage, energy, which channels are summed), and Cho/Schaepkens measured at
their own biases, not at 350 eV.  C2 is therefore the weakest of the four
"unsatisfiable" verdicts and should be re-posed against the source's stated
measurement condition before it is used to reject anything.

## 6. Status

Two hooks landed as pure refactors (`_complex_energy_factor`,
`_THERMAL_F_STICKING`); defaults reproduce the previous expressions exactly and
the suite is green at **1174 passed, 1 skipped**.  **No model assignment was
changed on this pass** — the beam-selected assignment fails the depth gate, so
landing it would trade one miss for another, and the forecast-before-spend rule
says do not buy a run for that.

What the pass produces is the structure: two constraints select assignments
uniquely, four cannot be satisfied by any sourced combination, and each of the
four names a different missing term.  None of them is a parameter to move.
