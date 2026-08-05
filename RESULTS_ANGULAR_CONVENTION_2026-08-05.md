# The class-1 angular normalisation: measured against three datasets, NOT landed

Task: find the term that keeps the surface model ion-limited where the process
is neutral-limited, per `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`.

## 1. The question

Krueger Eq. (2.40) (thesis L2399) is

    p(E_i, t) = p0 * (E_i^n - Eth^n)/(Er^n - Eth^n) * f(t)

and neither he nor Huang (L2287, "p0 is the yield at the reference energy")
states the **scale** of `f`.  petch had been reading `f(0) = 1` for both
angular classes, so `p0` was the normal-incidence probability.  The research
pass identified this as the single remaining lever on neutral sensitivity at
HARC energies: with the energy term capped, the only thing separating the
chemical row from the physical row at a floor that sees near-normal ions is
their angular class.

## 2. The discriminating verbatim

Huang L2290-2296 describes both classes in one sentence:

> "For physical sputtering, f(theta) is an empirical function with a maximum at
> 60 deg, **reduced probability at normal incidence** and zero probability at
> grazing incidence. For chemically enhanced etching, f(theta) is **unity for
> normal incidence** and angles up to 45 deg, with a monotonic roll-off to zero
> probability at grazing incidence."

The sentence only carries information if the two classes share one scale:
class 2 is stated unity at normal, class 1 is stated *reduced* there.  Under a
common `f(0) = 1` the phrase is false for class 1 — it would equal class 2 at
normal incidence, which is precisely the condition at a HARC etch front.

That is a reading, not a proof, so it was measured rather than asserted.

## 3. The 2x2, against three independent published observables

`scripts/angular_convention_discriminator.py` sweeps convention x shape.  No
constant is fitted anywhere in the table; `B = 9.3` is Krueger's cited source
(Kress, JVST A 17, 2819) and `B = 1.7` was the substitution landed in `830e5c5`.

| convention | f(0) | Gray dynamic range | total peak/normal | depth factor |
|---|---|---|---|---|
| f(0)=1, B=1.7 (previous) | 1.000 | 0.873 | 1.11 | 1.000 |
| f(0)=1, B=9.3 | 1.000 | 0.873 | 2.04 | 1.000 |
| peak-normalised, B=1.7 | 0.764 | 0.667 | 1.09 | 0.932 |
| **peak-normalised, B=9.3** | **0.240** | **0.210** | **1.34** | **0.780*** |
| **measured / required** | — | **0.255** (band 0.20-0.30) | **1.30-1.33** | **0.78** |

`*` this depth column is the flawed forecast that §4 overturns — see below.

Sources for the three columns:

- **Gray dynamic range** `Y(F/Ar+ -> 0) / Y(saturated)`, 350 eV Ar+ on SiO2:
  Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), replotted as Kwon
  (ScD, MIT DMSE 2004) Fig. 3.4 p. 76 — floor 0.28, plateau 1.10.
- **Total oxide yield peak/normal**, all channels at HARC coverage: Cho,
  *JVST A* **18**, 2705 (2000) ~1.30; Schaepkens, *JVST A* **16**, 3281 (1998)
  ~1.33.
- **Depth factor**: the multiplier on removal per ion at the measured ml19
  etch-front energy (3406 eV) needed to bring the extrapolated 60 s depth from
  ~1066 nm into the 825 +/- 5% gate.  This column is a *forecast*, graded by
  its own run; the first two are direct measurements.

One cell appeared to hit all three.  It does not, and the next section is why.

## 4. The depth column was wrong, and it refutes the change

The depth forecast above was computed in **beam mode** — an F beam with no
carbon flux, so no polymer film.  The feature has one.  Re-forecasting at
Krueger's Table 6.1 wafer fluxes with the precursor and O channels live, at the
measured 3406 eV etch-front energy (`recession_velocity_m_s`, steady state):

| configuration | rate (nm/s) | film (nm) | relative | 60 s depth |
|---|---|---|---|---|
| previous (f(0)=1 both, oxide B=1.7) | 12.057 | 0.085 | 1.000 | ~1066 (+29%) |
| peak-normalised, all class-1 rows | 2.719 | 0.343 | 0.226 | **~240** |
| peak-normalised on oxide rows only | 6.522 | 0.085 | 0.541 | **~577** |
| **required** | — | — | **0.78** | **[784, 866]** |

Neither lands in the gate band, and the uniform application is catastrophic.
The mechanism of the miss is a feedback the beam-mode forecast could not see:
`f(0) = 0.24` on the **polymer** row cuts film sputter 4.17x, the film thickens
4x (0.085 -> 0.343 nm), and the thicker film throttles interface energy to the
oxide underneath — a compounding loss on top of the direct one.

**The change was therefore not landed and no run was bought.**  This is the
forecast-before-spend rule doing its job: the beam-mode number (0.78) would
have justified a 60 s run that the coupled number (0.23-0.54) says would have
produced a 30-75% under-etch — a worse miss than the +29% it was meant to fix.

### What survives, and what it now localises

The two **direct measurements** still select peak-normalisation for the
**oxide** rows: dynamic range 0.210 against Gray's 0.255, and total angular
dependence 1.34 against the measured 1.30-1.33.  Those are not overturned by a
depth forecast; they are measured, and the previous convention misses both
(0.873 and 1.11).  What the depth column adds is that the *same* reading
applied to the **polymer** row is refuted.

Three readings survive that tension, and the sources in hand cannot separate
them:

1. **The class label compresses two different curves.** Krueger's legend gives
   one "class 1" for oxide, mask and polymer alike, but the measured angular
   yields differ by material: SiO2 in fluorocarbon peaks 1.30-1.33 (Cho 2000,
   Schaepkens 1998) while an FC *film* peaks **1.448 at 65 deg** under Ar+
   (Barklund & Blom, JVST A 10, 1212 (1992), via RESEARCH_VERIFY_HUNT).  If
   each row carries its own measured curve, the polymer row is nearly flat and
   the oxide row is strongly peaked — which is exactly the split the data
   above wants.  This is the most likely resolution and it is a *measurement*
   change, not a convention change.
2. **The depth channel has a second defect of its own**, so no single
   normalisation can satisfy all three columns.  The depth wants a factor
   ~0.78 where the oxide-only convention delivers 0.54; the residual would sit
   in the complex-row energy form (the ZBL-vs-published-linear question left
   open in `RESULTS_LIMITING_REGIME` §3, worth 2.05/1.72 = 1.19x in the right
   direction) or in the ion-energy overestimate Krueger flags himself
   (L4884-4888).
3. **The convention is f(0) = 1 after all** and the two measured agreements are
   coincidence.  Least likely — two independent datasets, one of them the
   in-chemistry measurement that the withdrawn B=1.7 substitution was chosen
   to match — but it cannot be excluded from text alone.

Reading 1 is the next change to try, and it is cheap: digitize Barklund & Blom
for the polymer row, keep the oxide row peak-normalised, and re-forecast the
coupled rate before any run.

## 5. What the surviving evidence does resolve

Nothing in `src/` changed this pass — the tree is back at `830e5c5` behaviour.
What the pass produced is measurement, and three results outlive the refuted
change:

1. **`OXIDE_B = 1.7` is on weaker ground than when it landed.** The
   in-chemistry measurement of ~1.3 that motivated it bounds the **total**
   angular dependence of the oxide etch yield, which at HARC coverage is
   dominated by the chemically-enhanced (class 2) channel, not by the physical
   shape.  Reproducing it by lowering the *physical* shape gives 1.11, i.e. it
   does not actually reproduce the measurement it was chosen for; the
   peak-normalised Kress form gives 1.34, which does.  Whatever replaces it,
   `B = 1.7` should not be defended on the Cho/Schaepkens ratio alone.
2. **The timestep and the depth channel are coupled through the same term.**
   The stiffness driver behind the ml21 dt collapse is off-normal removal
   running above `p0`; the depth miss is removal at *normal* incidence.  Any
   class-1 normalisation moves both at once, in opposite senses, which is why
   the change that fixes one over-corrects the other.  They cannot be graded
   separately.
3. **The polymer row is the binding constraint**, not the oxide row.  The
   oxide-only variant misses the depth band by 26% (577 vs 784 low edge); the
   uniform variant misses it by 70%.  Every future class-1 change should be
   forecast on the polymer row *first*, in coupled mode, because that is where
   the 4x film-thickening feedback lives.

## 6. Scope item 2 — site turnover, checked and already present

Verbatim, first passivation of pristine oxide (Krueger L6006-6012):

    SiO2(s) + CF  -> SiO2CF(s)    0.278
    SiO2(s) + CF3 -> SiO2CF3(s)   0.2

against further passivation of an already-complexed site (Krueger L6556-6564;
Huang L10214-10222 carries 1e-4 for the same rows):

    SiO2CF(s)  + CF -> SiO2C2F3(s)   0.0002
    SiO2CF2(s) + CF -> SiO2C2F3(s)   0.0002
    SiO2CF3(s) + CF -> SiO2C2F4(s)   0.0002

a 1390x asymmetry.  petch already carries it, as Langmuir site blocking:
chemisorption enters through `site_open = (1 - theta_film) * (1 - theta_f_layer)`,
so an occupied layer accepts essentially nothing.  Measured in
`tests/test_neutral_limited_gates.py::test_repassivation_asymmetry_is_present`:
carbon uptake into a fluorine-saturated layer is more than 1390x below uptake
into a pristine one.  **No change needed** — E1 and E3 of the research doc's fix
list were already in the model, which is why the flat response could not be
fixed on the supply side.

## 7. An unrelated defect found while gating, not fixed here

Krueger's mechanism has **no thermal-F-on-bare-oxide row at all**.  Every
`SiO2(s) + F` entry is for the *ion* F+ (activation 0.9, sputter 0.0852,
L5905-5909).  Thermal F only fluorinates an **already-complexed** site, at
p = 0.1 (L6548-6555, "Fluorination of passivated surface"):

    SiO2CF(s)  + F -> SiO2CF2(s)   0.1
    SiO2CF2(s) + F -> SiO2CF3(s)   0.1

petch's `f_direct` admits thermal F into the bare layer at sticking **1.0**.
For the Krueger feature this is inert — his gas mix (Table 6.1) has no atomic
F — but it is wrong by 10x against the closest source row for any chemistry
that does deliver F, which includes the SF6/O2 silicon arm and the Gray beam
case gated here.  Recorded as the next chemistry-side item; it is the lever on
the N1 *half-rise* position (petch 1.6 against Gray's 27), which this pass does
not address.

## 8. E8 — declared next transport item, not implemented

Huang L5716-5727: above AR 10, **>95%** of the radicals reaching the etch front
are neutralised and thermalised CFx+/CxFy+ ions, not conducted thermal
radicals; a 15%->60% CFx+ sweep gives ~3x floor radical flux and **-30%** etch
rate (L5742-5754).  petch's cascade thermalises spent energetic species out of
the ledger entirely.  The spec:

1. A cascade particle dropping below `E_c` (10 eV) with a fluorocarbon identity
   re-enters the **thermal radical** ledger at its current position rather than
   being discarded.
2. Non-fluorocarbon partners (Ar) stay inert — verbatim, they "diffuse out of
   the feature with no surface reactions".
3. The returned population is then subject to the same conductance and sticking
   as any thermal radical, so it does not arrive free.
4. Gate: floor `radical / (ion + hot neutral)` at AR 10-20 against Huang's
   published **0.33-0.5** (GATE N3), and the sign check that raising the
   fluorocarbon ion fraction *lowers* the etch rate.

This is transport, not surface kinetics, and belongs to its own change.

## 9. What a graded run must show, when one is warranted

Preregistered, against `scripts/grade_ml21_final.py`:

| gate | target | prior (ml19) |
|---|---|---|
| etch depth @ 60 s | 825 +/- 5% = [784, 866] | ~1066 extrapolated (+29%) |
| mask aperture | 45 nm | 50.9 equilibrium (+13%) |
| constriction depth | 200 (SEM) / 271 (sim) nm | 170 |
| mask remaining | 850 nm | 850.2 (exact) |
| closure/etch, t >= 8 s | 0.0310 | in band |
| run completes 60 s | no dt collapse | blocked at t=46.2 (ml19), dt/2.6 (ml21) |

No run is warranted yet: the coupled forecast (§4) puts every sourced
convention outside the band.  The next candidate is reading 1 of §4 (the
per-material measured angular curves), and it must clear the *coupled*
forecast, not the beam-mode one, before a box is provisioned.

## 10. Reproducing this pass

    python scripts/gate_n1_beam.py                    # GATE N1/N2 baseline
    python scripts/angular_convention_discriminator.py # the 2x2, beam mode
    pytest tests/test_neutral_limited_gates.py         # the gates, pinned

Artifacts in `results/curated/neutral_limited_gates/`:
`gate_n1.json` (the beam curve against Gray 1993) and `angular_convention.json`
(the 2x2; convention-independent, since the script patches the shape internally).

The coupled forecast of §4 is not a committed script — it is four
`steady_state` calls at Krueger's Table 6.1 fluxes with the class-1 shapes
patched, and its numbers are tabulated above.  If a future pass repeats it,
promote it to a script: **every class-1 forecast must be run coupled**, because
beam mode systematically overstates the surviving rate (0.78 against 0.54 for
the identical change).
