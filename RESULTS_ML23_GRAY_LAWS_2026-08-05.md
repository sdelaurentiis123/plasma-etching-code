# ml23 — Gray's measured laws at feature scale (§5.2 path executed)

Confirms `2f1e218`, which replaced Krüger's two-row `n = 1` anomaly on the SiO₂
ion channels with Gray's own measured √E laws. The dossier (`866e69a`) named
this as the single path to closing the depth gate; this run executes it and
grades the result.

Run: `ml23-gray-12s`, dx = 10 nm, deterministic extruded 2-D transport,
`literature_v1` cascade, mixed layer, volatilization 1.0 — identical to
ml18/ml19/ml21 in every respect except the two SiO₂ rows. Cut at t = 2.80 s
(14 steps) once the trajectory was linear and the projection unambiguous;
see §6.

## 1. What the laws do at the feature front

The etch front runs at a *measured* 3406 eV, mean cosine 0.768
(`RESULTS_FLOOR_DELIVERY_2026-08-05`), not the 140 eV the previous factors were
normalised at. Per incident ion at that energy:

| channel | pre-Gray | Gray | ratio |
|---|---|---|---|
| complex / chemically-enhanced (needs F) | 3.889 | 2.987 | 0.77× |
| bare / physical sputter (F-free) | 4.060 | 0.752 | **0.19×** |

The physical channel was **5.4× too strong**, exactly as the row-vs-yield
measurement predicted (`74eb2fa`). That one number explains the campaign's
central pathology: a *fluorine-free* channel supplying half the removal at
4.06 units/ion etches at full speed however starved the surface is, so the
model could not be neutral-limited, could not produce ARDE, and over-etched.

## 2. The regime is fixed

0-D at the front energy, shadowing the radical supply at fixed ion flux:

| radical supply | etch response |
|---|---|
| open field | 1.00× |
| 20× shadowed | 0.057× |
| 100× shadowed | 0.010× |

Pre-Gray the identical sweep moved the rate **under 1 %** across a 50× range
(`RESULTS_LIMITING_REGIME_2026-08-05`). The floor now rides its radical supply,
which is the precondition for ARDE to exist at all.

## 3. Feature scale: matched simulated time

| t (s) | ml23 (Gray) | ml19 (pre-Gray) | ml21 (pre-Gray) |
|---|---|---|---|
| 0.50 | 3.46 nm / 87.87 nm | 11.00 nm / 87.47 nm | 11.03 nm / 87.92 nm |
| 1.00 | 6.95 nm / 86.75 nm | 22.03 nm / 86.00 nm | 22.99 nm / 86.84 nm |

(depth / aperture)

- **Floor rate 6.97 nm/s**, against ~22 nm/s pre-Gray — **3.2× slower**.
- **Aperture is untouched**: 86.75 vs 86.00 nm at t = 1.0. The Gray laws move
  the oxide rows only and the mouth is polymer-governed, so this separation is
  expected and confirms nothing leaked across channels.
- Projected 60 s depth **418 nm against the band [784, 866] — −49 %**,
  reproducing the 0-D forecast (`b7cdc32`) to within a percent. Feature
  coupling did **not** rescue it, so the 60 s endpoint was not bought.

The depth gate moves from **+29 % over to −49 % under**: a sign reversal, not a
closure. The truth lies between two sourced options that the campaign's
constants cannot reach from either side.

## 4. The ARDE sign

Two independent readings, both reported with their limits.

**Trajectory (rate vs accumulated depth), fitted slope:**

| run | rate span | slope | sign |
|---|---|---|---|
| ml23 (Gray) | 6.97 → 6.94 nm/s over 5.5–19.5 nm | **−0.0033** | flat / weakly decelerating |
| ml21 (pre-Gray, comparable early span) | 22.00 → 39.63 nm/s over 8.1–38.0 nm | **+0.7530** | **accelerating (anti-ARDE)** |
| ml19 (pre-Gray, full 46 s) | 21.97 → 9.00 nm/s over 16.4–835.9 nm | −0.0065 | decelerating |

Over comparable early spans the sign **flips**: the pre-Gray configuration
accelerated as it deepened (the anti-ARDE pathology), the Gray configuration
does not. Two caveats, stated rather than smoothed: ml23's span covers only
14 nm of oxide depth (mask-dominated AR 9.44 → 9.6), so the magnitude is not
resolvable; and ml19's full-run deceleration is **confounded** — its aperture
closed 88 → 50 nm over the same span, so that curve mixes true ARDE with mouth
throttling. Only the frozen-geometry route separates them.

**Frozen geometry (delivery, `scan_gray.json`, straight walls, fixed aperture):**

| oxide AR | direct ion | hot neutral | total energetic | CF₂ (radical) |
|---|---|---|---|---|
| 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| 2 | 0.975 | 1.214 | 0.996 | 0.946 |
| 4 | 0.950 | 1.430 | 0.992 | 0.915 |
| 8 | 0.899 | 1.857 | **0.984** | **0.897** |

Energetic delivery is flat to **−1.6 %** at AR 8 — the cascade still almost
exactly cancels ion shadowing, unchanged by this work, since the Gray laws are
chemistry and this is transport. What changed is which of these curves the rate
follows: radicals fall **−10 %** over the same span, and the model is now
coupled to *them* rather than to the flat energetic curve. That is the
mechanism by which ARDE becomes possible; its magnitude at these aspect ratios
is small because the radical decline is small.

## 5. Where the residual now lives

Gray defines β_e as *"the number of SiF₄ molecules removed from fluorine
saturated surface regions per incoming ion"* (Eq. 5-35). At the front energy
that is **2.99 SiF₄/ion**, against the **≈3.15 units/ion** Krüger's own blanket
arithmetic requires (`74eb2fa`) — **95 % agreement**.

So the chemical channel alone carries the observed etch **if the floor is
fluorine-saturated**. Ours is not: with the physical channel correctly weak,
the model's floor is F-starved and etches at half the observed rate. The depth
residual is now a **fluorine-delivery-and-coverage question at the etch front**,
not an ion-yield question — better posed than anything the campaign carried
before. Three candidates, all declared, none adopted here:

1. **Thermal-F sticking.** `_THERMAL_F_STICKING = 1.0` in petch; Gray's thesis
   prints **0.02**, and Krüger carries *no* thermal-F-on-bare-oxide row at all
   (only F onto an already-complexed site, at 0.1). Sourced, deliberately not
   landed — a second physics change would have confounded this grading. Note
   its sign: lowering it reduces F further.
2. **Thermalised-ion return to the radical ledger** (element E8 of
   `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05`, specified, unimplemented).
   Huang measures **> 95 % of the floor's radicals above AR 10** to be
   thermalised CF_x⁺ rather than conducted thermal radicals. petch delivers no
   such population, so a floor that should be F-fed is F-starved by
   construction. Largest named candidate, and a *transport* feature.
3. **The √E extrapolation.** The front is 1.7× above Gray's top measured point
   and ZBL `s_n ∝ √ε` holds only for `ε ≪ 0.01` where `ε = 0.048` here, so √E
   is the **upper** end of the defensible band — the true yield at 3406 eV is at
   or below what was used, which widens rather than closes the gap.

## 6. Gate table

| gate | target | ml23 | verdict |
|---|---|---|---|
| regime: rate responds to radicals | must respond | 100× shadow → 0.010× | **PASS** (was < 1 %) |
| channel magnitudes | Gray absolute | 0.752 / 2.987 per ion | **PASS** — sourced |
| beam dynamic range (unfitted) | 0.255 measured | 0.227 | **PASS** |
| aperture vs pre-Gray | unchanged | 86.75 vs 86.00 @ t = 1 | **PASS** — no leakage |
| ARDE sign | must not accelerate | −0.0033 (was +0.7530) | **PASS** — sign flipped |
| depth @ 60 s | 784–866 nm | ~418 nm projected | **MISS** — −49 % |
| closure/etch | 0.0310 | 0.1478 | **MISS** — arithmetic, see below |

The closure/etch miss is arithmetic rather than new physics: absolute mouth
closure is **unchanged** from the pre-Gray runs (1.17 vs 1.35 nm/s per side at
the same instant), so cutting the etch rate 3.2× inflates the ratio by the same
factor. It will move with the depth residual, not independently.

**Run cut at t = 2.80 s.** The trajectory was linear to four digits (rate 6.97 →
6.94 nm/s across 14 steps), the 60 s projection was already unambiguous at
−49 %, and the directive gated the 60 s endpoint on the projection being
band-plausible. Buying 50 more minutes of GPU to re-measure a settled number
was the wrong trade; the endpoint comparison against ml18's t = 12 s state is
therefore **not** in this record.

## 7. Verdict

The §5.2 path is **executed and partially successful**. It fixed what it was
predicted to fix — the limiting regime, removing a 5.4× magnitude error against
absolutely-calibrated published data, and flipping the ARDE sign — and it did
not close depth; it reversed the error's sign. The stop rule stands: depth is a
declared, bounded, now **better-localised** open item, and the named path
forward is fluorine delivery at the front (candidate 2), not further ion-yield
selection.
