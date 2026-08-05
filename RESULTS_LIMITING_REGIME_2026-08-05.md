# Why the floor is ion-limited, and what a unity yield cap would have cost

Task: find the term that puts our surface model in the ion-limited regime where
the process is neutral-limited, and land the class-1 angular discrimination.

## 1. The unity-cumulative hypothesis: tested, refuted, not landed

Huang's thesis states the MCFPM selection rule verbatim (L2313-2320):

> "Since the probability of a particle striking the surface upon arrival is, by
> definition, unity, the following procedure is followed to normalize selection
> of reaction probabilities. The cumulative yield of all allowed processes for
> the energy and angle of incidence is computed ... If the cumulative yield is
> greater than unity, the elastic yield is reduced so that the cumulative yield
> is unity. If after scaling the elastic yield to zero, the cumulative yield is
> still greater than unity, then the yields of all processes are scaled to
> provide a unity cumulative yield."

Our rows are unbounded in energy, and the oxide pair is far past unity at the
measured etch-front energy:

| E (eV) | bare SiO2 | complex | cumulative |
|---|---|---|---|
| 140 | 0.085 | 0.147 | 0.232 |
| 500 | 0.523 | 0.565 | 1.088 |
| 1500 | 1.741 | 1.721 | 3.462 |
| 3406 (measured front) | 4.060 | 3.889 | **7.949** |

Implemented the normalisation, then checked it against the source's own rate
before landing it — and it fails:

    Krueger 825 nm / 60 s = 13.75 nm/s over SiO2 at 2.2e28 m^-3
    => 3.03e20 formula units/m2/s against his 9.6e19 m^-2s^-1 ion flux
    => 3.15 units per incident ion (1.09 even at 2.9x funnelled delivery)

A unity cap allows at most 1.00 unit/ion, i.e. **4.36 nm/s at 1x delivery** —
below his 13.75 nm/s at any plausible funnelling. So the rule cannot be a cap
on removal per ion in a continuum layer. It is a *selection* rule in a
cell-based Monte Carlo: a pseudoparticle selects at most one process, but it
carries a numerical weight of many real ions and removes a whole mesh cell, so
units-removed-per-real-ion is set by that weight calibration, not by the
probability array. **Reverted; the change would have made the model 3x too
slow.** The receipt is kept because the arithmetic is what licenses our
continuum kernels to exceed unity at all.

## 2. The regime term, measured

The oxide has two removal channels with complementary coverage weights:

    sif4 = capacity(E) * theta_F        (complex; consumes 4 F)
    bare = kernel_bare(E) * (1 - theta_F)   (physical sputter; consumes NO F)

Their yields are nearly equal at feature energies:

| E (eV) | bare | complex | bare share |
|---|---|---|---|
| 140 (reference) | 0.085 | 0.147 | 36.7% |
| 500 | 0.523 | 0.565 | 48.1% |
| 1500 | 1.741 | 1.721 | 50.3% |
| 3406 | 4.060 | 3.889 | 51.1% |

so the total `capacity*theta + bare*(1-theta)` becomes **independent of
theta_F** — the two channels trade off exactly. Measured directly:

| neutral delivery | theta_F | rate |
|---|---|---|
| 0.02x | 0.056 | 7.59 nm/s |
| 0.10x | 0.201 | 7.58 nm/s |
| 1.00x | 0.664 | 7.54 nm/s |

**The fluorine coverage responds correctly to starvation — a 12x collapse — and
the rate does not move.** That is the whole ion-limited regime, and it is not a
bug in the F budget: supply, coverage and clamps all behave. It is a
consequence of two published rows whose energy scalings cross at the feature
energy. At Krueger's 140 eV reference the same expression is strongly
theta-dependent (0.085 + 0.062*theta, 73% swing); the model is neutral-limited
there and ion-limited at 1.5-3.4 keV.

Named consequence: ARDE from radical starvation cannot appear at these
energies regardless of transport, which closes the question opened by the
cascade audit (`fd61bb7`) — the cascade is right, the F budget is right, and
the regime is set by the energy scaling of the two removal rows.

## 3. Open, with a source line on each side

Our complex channel uses the DEKNOB-derived ZBL deposited-energy shape,
`0.1471 * eps(E)/eps(140)`, while the published row reads `0.1471 35 1 140 2`,
i.e. `(E-35)/(140-35)` with n=1. At 1500 eV that is 2.05 against our 1.72, and
it would restore an 18% theta-dependence (complex above bare rather than equal
to it). Both choices are sourced — the ZBL form is the K24-DEKNOB result that
retired a fitted knob and is validated against the power sweeps — so this is a
model-choice question needing its own graded run, not a transcription fix.
Recorded as the next depth-channel candidate.

## 4. Class-1 angular discrimination: landed

The ml21 timestep collapse is explained. Krueger's cited class-1 source is
Kress et al. (1999), a molecular-dynamics study of "Cu and Ar ion sputtering of
Cu(111) surfaces" — its shape gives peak/normal 4.17. Applying it to the oxide
rows (`4c66df1`) amplifies off-normal removal by **2.54x** at ~50 deg; the
measured ml21/ml19 timestep ratio is **2.63x**. Agreement to 4%.

The only angular sputter measurements on SiO2 in fluorocarbon bound
peak/normal at 1.30 (Cho 2000, JVST A 18, 2705) and 1.33 (Schaepkens 1998,
JVST A 16, 3281). The oxide/mask rows now use B = 1.7 (peak 1.31, inside the
measured band); the polymer row keeps Krueger's cited B = 9.3, so every
validated lip/mouth result is untouched. `f(0) = 1` for any B, so all
normal-incidence and blanket results are bitwise unchanged. Off-normal
amplification drops 2.54x -> 1.16x, which should restore the timestep.

Flagged follow-up with a receipt already in hand: the polymer row has a
measured in-chemistry counterpart too — Barklund & Blom, JVST A 10, 1212
(1992), Ar+ on a fluorocarbon film, peak 1.448 at 65 deg — 2.9x below the
Kress form now carried there. Changing it needs its own graded run against the
lip results.

## 5. No run spent

Forecast: the etch front sits at 0.92 deg median tilt, where both class-1
shapes are 1.000 to four digits, so the depth gate is **unchanged** by section
4 and still misses. Per forecast-before-spend, no confirmation run was bought;
the change's value is the restored timestep, which unblocks the 60 s endpoint
whenever a depth-channel change is ready to grade with it.
