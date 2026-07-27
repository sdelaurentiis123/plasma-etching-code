# Mouth-equilibrium research — why our lip sits at 22.5 nm and Krüger's at 45 nm

Date: 2026-07-27
Question: our zero-fitted-constant mixed layer reproduces Krüger 2024 base depth
(−4.1%), the O2=0.5 clog (exact), and O2=2.5 necking-absence, but the mask mouth
equilibrates at ~22.5 nm vs the experimental 45 nm, and mask-top erodes 133 nm vs
~0. This one residual drives all five graded misses. What compensating
mouth-region film **removal** (or suppressed deposition) does Krüger's MCFPM
carry that we lack, and what moves our lip from 22.5 to 45 while keeping
mask-top erosion ~0 and not breaking depth?

Sources actually read this session:
- **Krüger, Zhang, Luan, Park, Metz, Kushner, JVST A 42, 043008 (2024)** — read
  in full from the local text `tmp/pdfs/krueger-2024.txt` (Sec. IV mechanism,
  Table II reaction list, Table III tuning ranges, Table IV target metrics,
  **Table V final tuned values**, Sec. V metric definitions, Sec. VIII
  necking/clogging coupling, Sec. IX O2 + power transferability).
- Repo transcriptions of the surface constants: `src/petch/surface_kinetics.py`
  (`krueger_2024_reduced_projection`, `huang_kushner_2019_reduced_projection`),
  `src/petch/mixed_layer.py`, `src/petch/mixed_layer_mechanism.py`,
  `src/petch/boundary_transport_3d.py::split_grazing_ion_reflection`.
- **Krüger PhD thesis (2024), DOI 10.7302/23106** and **Huang JVST A 37, 031304
  (2019)** — retrieval **BLOCKED this session**: the only live copies are on
  `cpseg.eecs.umich.edu`, which returns a TLS error ("unable to verify the first
  certificate") to WebFetch; the Scribd mirror is JS/paywalled. Every number
  attributed below to the thesis Appendix B / Huang Table I is taken from the
  repo's own verbatim transcription and is marked **[VERIFY vs thesis]** where the
  article PDF does not independently confirm it.

---

## 1. What Krüger's article pins down exactly (no verification needed)

### 1a. The five tuned parameters — Table V (verbatim from krueger-2024.txt L830–840)

| Symbol | Description | Single feature | **4-feature array** |
|---|---|---|---|
| ps,SiO2 | Physical sputter probability of bare SiO2 | 0.0852 | **0.0909** |
| ps,SiO2CFXY | Sputter probability of SiO2·CxFy complex | 0.1471 | **0.1384** |
| pp,SiO2 | SiO2·CxFy complex-formation probability | 0.278 | **0.2729** |
| pe,poly | **O-based polymer etch probability** | 0.0423 | **0.0628** |
| pd,poly-AC | **Polymer deposition probability on mask** | 0.094 | **0.0842** |

Every constant in `mixed_layer.py` traces to the **4-feature-array** column:
`sputter (bare) 0.0909`, `complex 0.1384`, `complex-formation 0.2729`,
`oxidation_probability=0.0628` = **pe,poly**, and — the one worth noticing —
`sticking_probability=0.0842` **is pd,poly-AC**, the mask-region polymer
deposition probability, not the generic on-polymer 0.1. So on the deposition
side we are already using Krüger's *optimized mask* value. Good.

**Consequence for the mouth:** the deposition constant at the lip and the O-etch
constant at the lip both match Krüger's optimized numbers exactly. The 22.5-vs-45
gap therefore is **not** a wrong deposition or O-etch constant. It is either (a) a
missing *removal channel* at the lip, or (b) mis-*delivery* (transport/angular)
of the removal we already have. Sections 3–4 localize which.

### 1b. Polymer (P) sputter law — NOT tuned; inherited fixed

Polymer sputtering `P(s)+M(g)→P(g)+M(g)` is a Table-II reaction but is **not** in
Krüger's tuned set (Table III). Its yield is inherited from the parent Huang
mechanism and held fixed. Repo transcription (`surface_kinetics.py:854–857`,
Appendix B): **p0 = 0.9, ε_th = 20 eV, q = 0.5, ε0 = 500 eV, Kress B = 9.3, at
incident energy.** `mixed_layer.py` uses exactly this (L256, L281). This is the
audit-corrected law (the earlier 0.1384 mis-lift is fixed). **[VERIFY vs thesis]**
on the 0.9/20/0.5/500 quadruple — the article does not print it; only the repo
transcription of Appendix B does.

### 1c. Necking/clogging is governed by the O-etch balance, and Krüger says so explicitly

Sec. VIII, verbatim: *"A steady state polymer thickness occurs when [deposition
and removal] balance. Polymer removal occurs by sputtering and O-radical based
etching. **In oxygen rich gas mixtures, polymer removal can be dominated by
O-radical based etching**, mostly ground state atomic oxygen … necking and
clogging … can ultimately be controlled by the reaction probability of the O
based polymer etching."* Fig. 6: pe,poly=0.005 → full clog; 0.02 → necking but
open. Base is 0.0628, comfortably in the "necking but open" regime.

So in Krüger's own words the mouth equilibrium in the O-rich base case is set
**primarily by O-radical polymer etch, with sputter secondary**. This reframes
the campaign's "deposition vs ion sputter of the lip film" balance: the dominant
*remover* Krüger relies on is O, not ions. Our model has O-etch at the identical
0.0628 — so the question sharpens to **is our O-etch actually being delivered to
the lip at Krüger's rate, or is it under-delivered there?** (Section 4, receipt
candidate #4.)

### 1d. The mask stays intact because polymer covers it

Table IV target: initial AC mask 850 nm thick (Sec. IV L330), target **hm
(remaining mask thickness) = 850 nm** — i.e. essentially **zero net mask
erosion**. Sec. IX O2 sweep, verbatim: *"**The lack of polymer film on the mask
leads to increased mask erosion** with increasing O2 inflow … direct oxidation of
the AC mask by O2 was not included in the mechanism."* So in Krüger the mask top
is protected **by the polymer film that sits on it**; erosion only appears where
that film is stripped (high O2). Our base-case 133 nm mask loss therefore means
**our mask-top polymer film is not surviving** (or our AC sputter is not gated by
film coverage). This is a second, independent symptom of the same lip/top
polymer-balance error and a hard joint constraint (Section 5).

### 1e. wm definition and time behaviour (task items 4, 5)

- **Definition (Sec. V, Table IV):** wm = *"width of the mask opening **including
  deposition**"*, the **minimum** aperture at the necking location, extracted from
  the **final feature at 60 s** (the base-case SEM, Fig. 7). It is a minimum
  aperture, not an aperture at a fixed height. **Confirm our "opening" metric is
  the minimum aperture over height at 60 s** — if we currently report the aperture
  at the mask-plane height rather than the global minimum, part of the 22.5-vs-45
  gap is a metric mismatch. (Our scorecard labels it "opening (nm)"; the harness
  definition should be checked against "global minimum aperture".)
- **Time behaviour:** the article gives metric evolution over optimization
  *epochs* (Fig. 9b/14b/19b), **not** a wm(t) etch-time trace, so an
  overshoot-and-reopen curve is not shown. But Sec. VIII states the physics
  directly: necking reaches *"a steady state polymer thickness … when
  [deposition and removal] balance."* → Krüger's mouth **narrows monotonically to
  a steady 45 nm; it does not overshoot-and-reopen.** A per-etch-time wm(t) trace
  and the neck *height* are thesis-only (Fig. 7 of the thesis chapter). **[VERIFY
  vs thesis]** for the exact narrowing timescale.

---

## 2. What our controlled runs already proved (from the campaign doc)

- ml9a (atoms, Kress angular polymer sputter, **no reflection**): mouth 22.5, open,
  depth 790.8 (−4.1%). Checkpoint of record.
- ml9b (atoms + our `GrazingSpecularIonReflection3D`): mouth **0.0 (sealed)**.
- ml4 (weak sputter law, full grazing sputter, no Kress): mouth 30.
- Pattern the campaign extracted: **near-grazing wall-film sputter governs mouth
  survival**; suppress it (Kress-only *or* 95 % specular reflection) → seals;
  give it full grazing sputter → holds 22–30. Krüger holds **45** with **both**
  Kress **and** reflection active.

The key inference: **turning our reflection ON makes the mouth worse (seals it),
which is the opposite of Krüger, where reflection is part of what keeps 45 open.**
That inversion is the smoking gun and is explained in Section 3.

---

## 3. The mechanism we get wrong: reflection SIGN (subtractive vs additive)

### 3a. Krüger's reflection is a *redistribution that keeps sputtering*

Krüger (Sec. IV L344, Sec. V L458–460, Sec. V L478): ions that strike a surface
are *"neutralized during collisions"* into **hot neutrals** *"that retain energy
+ angle and continue to sputter/etch"*; *"the curvature of the polymer affects
the angle of the trajectory with which ions (hot neutrals) reflect from its
surface as they progress deeper into the surface"*; and the mask *"shadows and
reflects"* the broad-angle ions *"prior to entering the feature."* Crucially,
**the primary grazing impact still sputters the lip** (that is what the
Kress B=9.3 angular law is *for* — it enhances the sputter yield ~10× at grazing
vs normal), **and then the neutralized remnant reflects and sputters again**
downstream. Sputtering and reflection are **not** mutually exclusive at a grazing
impact; they are the same event's two outputs.

### 3b. Our implementation *subtracts* the reflected share from the lip

`boundary_transport_3d.py::split_grazing_ion_reflection` (L748–754):

```python
weight  = grazing_reflection_probability * (1 - cos^angular_exponent)   # →0.95 at grazing
primary = Population(..., flux * (1.0 - weight), ...)   # PRIMARY LIP FLUX SCALED DOWN
reflected_rate = weight * flux * areas[face]            # re-cast elsewhere
```

The docstring says it outright: *"the reflected share **leaves the primary
event** (its flux is scaled down in place)."* So at the grazing lip (cosθ→0)
**95 % of the flux is removed from the lip's own sputter event** and re-emitted.
The lip loses almost all of its sputter dose → the film there is no longer
removed → the mouth seals (ml9b). This is a modelling defect, not physics: a
grazing ion does **not** stop sputtering the lip because it is about to reflect.

### 3c. Why the correct (additive) sign fixes the mouth AND spares the mask top

Make the primary keep its **full** flux (its sputter yield already carries the
Kress grazing enhancement at the true per-event angle — ml9a machinery), and add
the reflected hot neutral as an **additional** population that sputters where it
lands. Then:

- **Lip:** gets its own grazing-enhanced primary sputter (the 22.5 baseline)
  **plus** reflected hot-neutral flux cross-bombarding from the opposite lip /
  mask shoulder → **more** removal → mouth opens 22.5 → toward 45.
- **Mask top (flat, horizontal):** sees **normal** incidence, Kress factor
  `(1+B(1−cos²θ))·cosθ = 1` at cosθ=1 (no grazing enhancement) and receives
  **negligible reflected flux** (reflection weight ∝ (1−cosⁿθ) → 0 at normal).
  So the additive reflected population is **geometrically concentrated at the
  grazing mouth and starved on the flat top** — it opens the mouth **without**
  adding mask-top removal. This is exactly the differential Krüger's "shadow and
  reflect off the mask" language describes.
- **Depth:** reflected hot neutrals *"progress deeper into the surface"* / funnel
  to the floor — additive reflection **feeds** the etch front rather than
  starving it, so it should not hurt (and may help) the −4.1 % depth result.

This single change flips reflection from "seals the mouth" (our ml9b) to "helps
hold it open" (Krüger). It is the **highest-leverage, best-evidenced** fix and it
is the only mechanism that simultaneously (i) widens the mouth, (ii) spares the
mask top, (iii) does not degrade depth.

**Uncertainty flag:** how much the reflected flux lands on the *opposite lip*
(widening) vs the *floor* (depth only) depends on lip **curvature**, which
Krüger calls out explicitly and which our voxel/level-set lip may under-resolve.
So additive reflection is directionally certain to widen but its **magnitude**
(does it reach 45, overshoot, or only reach ~30) is curvature-resolution
dependent. Receipt in 3d.

### 3d. Receipt that would confirm #3

Per-face energy/removal budget at the mouth lip voxel, ml9a vs ml9a+additive-
reflection: the lip face should gain a **`:hot_neutral` sputter term** roughly
comparable in magnitude to its primary grazing sputter term (removal at the lip
≈ doubles), while the **flat mask-top face gains ≈0** reflected term. Mouth width
should rise monotonically toward 45 and settle (Sec. 1e: steady state, no
reopen). If the lip removal doubles but the mouth only reaches ~30, the shortfall
is curvature-resolution (reflected flux going to floor not opposite lip), not a
missing chemical channel.

---

## 4. The remover Krüger leans on hardest: O-radical polymer etch delivery

Per Sec. 1c, Krüger says O-etch **dominates** polymer removal in the O-rich base
case and is the necking control. Our pe,poly=0.0628 and O-flux are Krüger's, and
in `mixed_layer.py` the film O-etch is
`ox_c = p_ox · O_flux · θ_film · x_c` (L338) — isotropic, no angular dependence,
applied wherever film exists. If our **transport under-delivers the O-atom flux
to the lip** (e.g. the already-narrowing mouth shadows O, a physical feedback
Krüger names: *"the narrowed opening can impede neutral gas transport into the
feature"*), the lip loses its dominant remover and over-narrows — a runaway that
seals faster than Krüger's. This is **mis-delivery of a channel we already have**,
not a missing channel.

**Receipt that would confirm #4:** at the mouth lip voxel, compare the local
O-atom flux used in `ox_c` against the incident (open-field) O flux, and compare
the lip's deposition rate against its (sputter + O-etch) removal rate. Krüger's
steady state requires removal ≈ deposition at 45 nm; if our lip shows removal ≪
deposition **and** local O-flux ≪ incident O-flux, the fix is in neutral
transport/re-emission to the mouth, not chemistry. Expected direction: restoring
lip O-flux widens the mouth. This should be checked **before** adding any new
chemistry, because if #4 is the cause, #3 alone will over-correct.

---

## 5. The mask-top constraint (hm≈850) and where crosslinking actually belongs

Our mask erodes 133 nm; Krüger's target is ~0 because polymer covers it (Sec.
1d). Two receipt-checkable causes, with fixes:

**5a. AC/mask sputter may not be gated by film coverage.** In `mixed_layer.py`
substrate removal is correctly gated by `open_area = max(1−θ_film,0)` (L380–393):
oxide is only removed where the film is open. **Verify the mask/AC arm does the
same** — if AC is sputtered even under polymer, the mask erodes through its own
protective film. Fix: gate AC sputter by exposed (unfilmed) AC area exactly as
oxide removal is gated. Direction: mask-top erosion → 0. No depth impact, no
mouth impact. **This is a pure-correctness fix and should be verified regardless.**

**5b. Crosslink SPUTTER-resistance — the one place it helps, and the sign trap.**
Krüger: deposited monomer *"can subsequently crosslink to create an actual
polymeric material that is **more resistant to sputtering**"* (Sec. IV L378), and
crosslinking changes *"sputtering probability"* (Fig. 5d). Our model already
lowers crosslinked *sticking* (`crosslinked_deposition` blend, `n_xl_film`) but
**PC and P sputter identically** — the resistance half is absent (audit class
#7). The **exact PC sputter factor is thesis-Appendix-B only and I could not
retrieve it → [VERIFY vs thesis].**

The sign trap: crosslinking is **ion-dose driven**, and *both* the flat mask top
and the lip receive heavy dose, so a naive PC-resistance term hardens **both** →
protects the mask top (good, helps hm) **but also hardens the lip** (bad, narrows
the mouth — the wrong direction for 22.5→45). It is therefore **only safe to add
PC sputter-resistance once #3 (additive reflection) provides grazing+reflected
removal at the lip strong enough to overwhelm the hardened lip film.** The
coherent joint state is: mask-top PC film survives the *weak, normal-incidence
(Kress≈1)* sputter → hm≈850; lip PC film, though equally hardened, is overwhelmed
by *grazing-enhanced (Kress~10×) primary + reflected* sputter → mouth opens to 45.
PC-resistance alone, without #3, will make the mouth worse.

---

## 6. Mechanisms to EXCLUDE (wrong direction for 22.5→45)

The audit's P2 queue names these as "complete-the-mechanism" additions, but for
our *specific* miss (mouth too **narrow** = need **more** removal / **less**
deposition at the lip) they push the wrong way:

- **Low-energy activation of polymer (Huang 5–30 eV, p0=0.3).** Task asked
  whether activation lowers sticking or enables removal — the repo transcription
  is unambiguous (`surface_kinetics.py:698`): activated sites take *"P* sticking
  probabilities, **ten times** the unactivated P values."* Activation **raises**
  sticking 10× → **more** lip deposition → **narrower** mouth. Adding it makes
  22.5 worse. **[VERIFY vs Huang]** on the 10× and the window, but the sign is
  clear from the transcription. Exclude as a mouth fix.
- **Energetic FC-ion direct polymer deposition (Huang p0=0.1, 5–70 eV, FC+
  family).** Adds deposition precisely in the ion-rich mouth → **narrower**.
  Exclude.
- **PC sputter-resistance applied uniformly** (see 5b): narrows the lip unless
  paired with #3. Not a standalone mouth fix.

These belong in a completeness pass, but none of them opens the mouth; they were
mis-queued as mouth candidates because "add the missing Huang states" was read as
"add removal," whereas these three are all **deposition/retention** channels.

---

## 7. VERDICT — ranked, with direction, magnitude, source, and confirming receipt

Diagnosis in one line: **our lip has the right chemistry constants but is missing
the additive reflected-hot-neutral grazing removal, and possibly under-delivered
O-flux; the mask top erodes because its polymer film isn't protected. Krüger's
45 nm is a steady state set by O-etch-dominated removal + grazing (direct +
reflected) sputter, geometrically differentiated from the normal-incidence mask
top.**

| Rank | Change | Direction / magnitude | Source (exact) | Confirming receipt |
|---|---|---|---|---|
| **1 (P0)** | **Flip reflection from subtractive to additive**: primary grazing event keeps full flux + Kress-enhanced sputter; reflected hot neutral is an *added* secondary population (retain 0.90 E, angle) that sputters where it lands. | Widens mouth 22.5 → toward 45; spares flat mask top (reflection weight→0 at normal); feeds floor → depth neutral/positive. Magnitude curvature-limited. | Krüger JVST A 42 Sec. IV L344 (hot-neutral neutralization), Sec. V L458/L478 (curvature-dependent reflection, mask shadow+reflect). Defect at `boundary_transport_3d.py:748–754`. Kernel constants p_grazing 0.95 / exp 3 / retention 0.90 already present. | Lip face gains `:hot_neutral` sputter term ≈ its primary term (lip removal ~2×); flat mask-top face gains ≈0; wm rises monotonically to steady value (no reopen). |
| **2 (P0-verify)** | **Audit O-atom flux delivery to the lip** before adding chemistry. | If lip O-flux ≪ incident, restoring it widens mouth; if it's fine, #1 carries it (else #1 over-corrects). | Krüger Sec. VIII (O-etch dominates removal; narrowed mouth impedes neutral transport — the feedback). `mixed_layer.py:338`. | At lip voxel: local O-flux vs incident O-flux, and deposition rate vs (sputter+O-etch) removal rate. Krüger steady state ⇒ removal≈deposition at 45. |
| **3 (P1)** | **Gate AC/mask sputter by exposed (unfilmed) area** (mirror the oxide `open_area` gate). | Mask-top erosion 133 → ~0; no mouth/depth impact. Pure correctness. | Krüger Sec. IX (mask erosion gated by polymer film coverage). `mixed_layer.py:380–393` shows the correct oxide gate to mirror. | Mask-top face shows persistent polymer film at 60 s and AC removal only on exposed area; hm → ~850. |
| **4 (P1, after #1)** | **Add PC sputter-resistance** (PC harder to sputter than P): raise threshold / lower p0 for crosslinked film. | Protects normal-dose mask-top film (helps hm); at the lip only safe because #1's grazing+reflected removal overwhelms it. Adding it *without* #1 narrows the mouth. | Krüger Sec. IV L378 ("more resistant to sputtering"), Fig. 5d. Numeric PC factor **[VERIFY vs thesis Appendix B]** — not retrievable this session. | With #1 on: mask-top film survives, lip film still removed; wm holds 45 and hm holds 850 jointly. Without #1: wm drops — the sign check. |
| **5 (P2)** | **Confirm wm metric = global minimum aperture at 60 s** (not aperture at mask-plane height). | Could recover part of the gap for free if we currently measure at fixed height. | Krüger Table IV + Fig. 7: wm is the minimum opening including deposition. | Re-extract ml9a opening as min-over-height; compare to 22.5. |
| **— (exclude)** | Low-energy polymer activation (10× sticking), energetic FC-ion deposition, uniform PC-resistance. | All **narrow** the mouth (deposition/retention channels). | `surface_kinetics.py:698` (activation = 10× sticking), Huang Table I. | n/a — do not add as mouth fixes. |

**Recommended order of operations:** (2) receipt-audit O-flux and lip balance
first — it's free and tells you whether the deficit is transport or removal.
Then (1) additive reflection — the structural fix, highest leverage, self-consistent
with the mask-top constraint. Then (3) mask-sputter gating — independent
correctness fix for hm. Only then (4) PC sputter-resistance, guarded by the sign
check. (5) is a cheap metric-definition confirmation to run alongside.

**What remains genuinely unverifiable without the thesis (all [VERIFY vs thesis]):**
the exact PC-vs-P sputter factor; the hot-neutral generation probability and any
threshold reflection angle off polymer vs mask; the wm(t) narrowing timescale and
neck height. The umich host (`cpseg.eecs.umich.edu`) is the only live source and
returned a TLS certificate error to every fetch this session — retry from a
network that trusts its chain, or pull the local thesis PDF into `tmp/pdfs/`.
