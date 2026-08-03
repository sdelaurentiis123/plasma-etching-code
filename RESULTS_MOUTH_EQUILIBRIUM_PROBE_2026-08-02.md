# Frozen-geometry interrogation of the trench neck (2026-08-02)

`scripts/mouth_equilibrium_probe.py`, artefacts in
`results/curated/mouth_equilibrium_probe/`.

The 60 s evolution runs answer "where does the mouth end up" only after hours of
wall clock, and they conflate transport, chemistry and level-set kinematics.
This probe asks the balance question directly: **at a prescribed neck aperture,
does the neck film grow or erode?** One transport gather per geometry (~60 s
CPU) replaces a campaign, and the aperture where the net normal velocity changes
sign is the equilibrium the evolution run would converge to.

Method: geometry built to the digitised Fig. 7 profile (90 nm at the mask top, a
smooth constriction to the prescribed aperture 250 nm into the 850 nm mask,
re-opening below, over a 300 nm etched oxide trench). The transport gather runs
at `duration_s=0` — no geometry motion, no topology risk — and the surface
mechanism is then relaxed on that frozen flux field until face velocities stop
changing. Chemistry is the ml16a (fig-7/Table-6.5) set with the P1a-corrected
angular lift and `grazing_ion_reflection=literature_v1`.

## 1. Sweep: there is no equilibrium aperture

| prescribed neck | net normal velocity | film growth | neck ion flux | mean incidence |
|---|---|---|---|---|
| 45 nm | **−7.71e−4 nm/s** | 7.71e−4 | 6.50e18 m⁻²s⁻¹ | cos 0.120 (83.1°) |
| 39 nm | **−6.66e−4 nm/s** | 6.66e−4 | 7.30e18 | cos 0.137 (82.1°) |
| 33 nm | **−5.78e−4 nm/s** | 5.78e−4 | 8.09e18 | cos 0.154 (81.1°) |
| 27 nm | **−5.19e−4 nm/s** | 5.19e−4 | 8.88e18 | cos 0.146 (81.6°) |
| 21 nm | **−4.69e−4 nm/s** | 4.69e−4 | 9.63e18 | cos 0.159 (80.8°) |

`equilibrium_aperture_nm: None` — **the net velocity is negative at every
aperture and never changes sign.** The neck closes monotonically, and the
closing rate *decreases* as it narrows (−7.7e−4 → −4.7e−4 nm/s). This is an
asymptotic approach to closure, not a balance: petch's neck has no stable
aperture, which is the frozen-geometry statement of the ml13–ml16 mouth
residual.

Note the substrate recession term (`etch_velocity_m_s`) is **identically zero**
at every neck face — correct, not a bug: the a-C mask is sputter-armoured
(0.001 @ 200 eV, O-inert 1e-5) exactly as Krüger's converged mechanism
specifies, so the mask itself does not recede. The entire neck balance lives in
the fluorocarbon film.

## 2. Lip budget at the 39 nm (experimental) neck

Area-weighted over the 144 selected neck faces, per relaxation step:

| term | units m⁻² | share |
|---|---|---|
| polymer deposited | 1.3818e17 | 100% |
| polymer removed | 1.1820e17 | **85.5%** |
| bare-substrate removed | 0.0 | — |
| complex removed | 0.0 | — |
| **net** | **−0.00067 nm/s** | 14.5% deposition surplus |

**Removal is 85.5% of deposition.** The neck is not starved of removal — it is
within 15% of balance and loses. This is quantitatively the regime Krüger
describes ("a steady state polymer thickness occurs when these contributions
balance"), and it means the mouth residual is a *~15% budget error*, not a
missing mechanism. Any of the following would close it: 15% more grazing
removal, 15% less lip deposition, or the combination.

That sensitivity is the headline for what to do next. Two measured anchors from
`RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md` both point at this margin and
both are larger than 15%:

* **You et al., *Coatings* 13, 1452 (2023)**: measured normalised etch rates lie
  *above* the cosine curve out to 50–60° incidence. If petch's grazing removal
  sits at or below cosine, the deficit is of exactly the needed size.
* **Izawa et al., *JJAP* 46, 7870 (2007)**: sidewall CF*x* sticking coefficient
  **0.004**, against the ~0.094 class petch applies on the mask. If the
  definitions are commensurate this is a 20× over-deposition on the lip — far
  more than enough. Definition check required before acting (`[VERIFY]`).

### Angle test: inconclusive from this geometry

The 80 neck faces span only **78–80° incidence** — a vertical wall in a
near-vertical beam samples one narrow angular band. 70% of them sit above the
(self-normalised) cosine reference, but with a 2° lever arm that number cannot
discriminate the You-2023 law. Testing the cosine claim needs faces across
30–70°, i.e. a tapered or faceted lip geometry, not this one. Recorded as an
explicit negative result rather than a passing grade.

## 3. Depth-resolved closure — petch pinches at the wrong place

Net normal velocity for every mask face at the 45 nm geometry, binned by depth
below the mask top (`depth_profile_45nm.json`):

| depth below mask top | faces | mean net (nm/s) | min net | mean ion flux (m⁻²s⁻¹) |
|---|---|---|---|---|
| **0–50 nm** | 56 | **−0.03418** | **−0.05037** | 3.14e19 |
| **50–100 nm** | 40 | **−0.01559** | −0.04012 | 4.54e18 |
| 100–150 nm | 40 | −0.00190 | −0.00434 | 1.42e19 |
| 150–200 nm | 48 | −0.00018 | −0.00045 | 1.87e19 |
| **200–250 nm** (SEM/MCFPM neck) | 48 | −0.00115 | −0.00316 | 9.96e18 |
| 250–300 nm | 48 | −0.00095 | −0.00305 | 6.40e15 |
| 400–700 nm | 240 | −0.001 to −0.002 | — | 1e14–7e17 |

**Closure is 30–50× faster at the mask top than at the neck.** The band that
actually pinches in petch is 0–100 nm below the mask top; Krüger's simulated
minimum is at 271 nm and the SEM's at 200 nm. The 200–250 nm band — where the
experiment necks — is nearly the *most stable* part of our mask wall
(−0.0012 nm/s, and the 150–200 nm band is essentially at balance, −0.00018).

Two consequences, and they reframe the whole residual:

1. **`mask_opening` is not measuring Krüger's `w_m`.** Our metric reports the
   throat wherever it is; in petch that throat forms within 100 nm of the mask
   top, so ml13/ml16's 11–25 nm number is a *top-pinch* aperture being compared
   against his *mid-mask neck*. This is quantitatively the "difference in the
   vertical position of the minimum in necking" he flags as unreproduced
   ([T] 4757–4774) — we have now measured it on our side.
2. **The physics error is at the top-of-mask lip, not in the neck.** The top
   band sees the *highest* ion flux in the feature (3.14e19, ~25% of incident)
   and still loses to deposition. Removal there is not shadow-starved; the
   deposition/removal ratio at the exposed lip is simply wrong.

## 4. Resolution

The dx = 0.005 µm repeat of the 39 nm evaluation was launched
(`results/curated/mouth_equilibrium_probe_dx/`) and had not converged when this
document was written; it is the one open item. Note the depth profile already
weakens the pure-discretisation hypothesis: a grid floor would bias every band
equally, whereas the measured closure is 30–50× concentrated in the top two
bands. `[PENDING]`

## Implied fix, in order

1. **Report the neck metrics, not one aperture.** Adopt Top CD / Necking CD /
   `z_neck` (Kwon 2024 definitions, already recommended by the literature pass)
   in `measure_krueger_metrics`. Until `z_neck` is reported we cannot tell a
   mouth error from a *location* error, and this probe says we have both — the
   mid-mask band nearest his neck sits at −0.0012 nm/s, far closer to balance
   than the 11–25 nm endpoint suggests.
2. **Attack the top-of-mask lip balance, not the neck.** The 0–100 nm band
   carries 30–50× the closure rate and the highest ion flux in the feature.
   Two measured anchors bracket the ~15% budget error, both larger than needed:
   the You-2023 above-cosine grazing removal, and Izawa's 0.004 sidewall CF*x*
   sticking against our 0.094 mask class (`[VERIFY]` the definition first — a
   20× over-deposition at the lip would dominate everything else here).
3. **Finish the dx = 0.005 repeat** to close the discretisation hypothesis
   quantitatively, though the depth localisation already argues against it.
4. Do **not** add faceting, charging, or redeposition on this evidence: the
   budget is 15% from balance, all three are refuted as the mouth mechanism in
   both research passes, and redeposition has the wrong sign.

## Cost

Five apertures + one depth profile = six transport gathers, ~60 s each on CPU,
about 25 min wall clock end to end — against ~5 h for a single 60 s evolution
run on a rented GPU. The probe is the cheap instrument this question needed.
