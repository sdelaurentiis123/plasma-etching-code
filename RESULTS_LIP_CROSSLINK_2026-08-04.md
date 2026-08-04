# The lip crosslink inversion: creation is a deposition event, not an ion dose

Reproduce the 0-D comparison: `scripts/lip_crosslink_check.py`.
Inputs: `RESULTS_EARLY_TRANSIENT_2026-08-04.md` (`43a2116`), the corrected
band audit `results/curated/lip_deposition_audit/audit_neck45_dx0.01.json`
(`cbbd2d6`), and Krüger's thesis (`tmp/pdfs/krueger_thesis.txt`).

## The question this pass was given

> Our lip film sits at crosslinked fraction **x_xl = 0.163**. Matching Krüger's
> closure requires **x_xl ≳ 0.9**. Both models drive crosslinking with the same
> ion-dose row and both see the same collapsed ion flux on a vertical wall.
> Why is his lip film crosslinked and ours not?

The premise in the second sentence is false, and that is the answer.

## 1. Our kernel is right, and dilution explains 0.163 exactly

`kernel_xl = flux x max(E_ion - E_interface, 0)` with
`xl_rate = kernel_xl / 25 eV x (1 - x)`: one atom converted per 25 eV absorbed
in the film. On the thick lip film (15.6 nm) essentially the whole ion energy
is absorbed, so ~60 atoms crosslink per 1500 eV ion — but the *areal* ion flux
on a 0.47°-tilted wall is what multiplies it.

Steady state of a growing film (fresh arrives, ions convert) is
`xl_rate = x·D + B`, so with `B ≈ 0`:

| term | value |
|---|---|
| deposited atoms `D` (audited: 6.77e19 units/m²/s x 2.69 atoms/unit) | 1.82e20 atoms/m²/s |
| required creation at x = 0.163 (`x/(1-x)·D`) | 3.55e19 atoms/m²/s |
| implied areal ion flux at 1500 eV | 5.9e17 m⁻²s⁻¹ |
| actual: 9.6e19 x sin(0.472°) x visibility 0.742 | **5.8e17 m⁻²s⁻¹** |

The measured 0.163 **is** the dilution steady state of ion-driven creation
against the collapsed grazing flux. No unit error, no suppression bug: the
kernel does exactly what it says. The defect is that this channel should not
be the dominant one.

## 2. Krüger creates crosslinks at deposition; ions only break them

Thesis §2.2.3, verbatim:

> "**Crosslinking occurs during the deposition of eligible materials** (Figure
> 2.2a) and b). Each material has a maximum number of crosslink partners
> associated with it, which is based on the number of available bonds (3 in the
> example depicted in Figure 2.2). **During deposition bonds to random eligible
> cell neighbors can be formed**, increasing the respective crosslink number...
> **Crosslinks can be broken by impinging particles such as ions, hot neutrals
> and photons.**"

Chapter 6, verbatim: *"The polymer, as deposited, consists of individual CxFy
radicals... **This radical based film can subsequently crosslink**... **Ion
bombardment can then break bonds (chain scission)**... Since these energetic
particles are typically delivered anisotropically to the surface, this
**spatially discriminate activation can result in shaping of the polymer
deposition**."*

Table 6.2 rows:

```
P(s) + P(s)         ->  PC(s) + PC(s)          Crosslinking
P(s) + P(s) + M(g)  ->  P(s)  + P(s)  + M(g)   Breaking of Crosslinking
```

The crosslinking row carries **no M(g)**. Appendix B confirms it from the other
side: every `(xs)` row is a *consumer* of the crosslinked state —

```
CF(xs) + Ar+  ->  EP     + Ar#    0.6  50  0.5  500  1     (crosslinked sputter)
CF(xs) + Ar+  ->  CF(s)  + Ar#    0.3   8  0.5  500  1     (breaking)
CF(xs) + CF   ->  CF(xs) + CF(s)  0.02                     (crosslinked sticking)
CF(xs) + O    ->  EP             0.0423                    (O etch)
```

— and **no row in the 1229-row table creates `(xs)`**. Creation lives in the
module, at deposition time.

**The geometry dependence is therefore inverted between the two models.**
Deposition is isotropic and does not collapse on a near-vertical wall; the ion
flux that *breaks* crosslinks collapses ~200x through the double cosine. So in
his model the lip is the *most* crosslinked surface in the feature — which is
his stated purpose ("spatially discriminate activation... shaping of the
polymer deposition"). In ours the lip was the *least* crosslinked, because we
made creation ride on the channel that vanishes there.

## 3. Implementation (zero new constants)

`mixed_layer.py`, added to `xl_rate`:

```
xl_rate += 2.0 * (dep_c + dep_f) * fresh_fraction * theta_film
```

- **2** is the printed row's own stoichiometry: `P(s)+P(s) -> PC(s)+PC(s)`
  converts *both* partners, so each deposition event converts the arriving
  material plus one eligible neighbour.
- **`fresh_fraction`** is his "eligible" gate — an already-crosslinked cell has
  no free bond to offer.
- **`theta_film`** is the neighbour-exists gate — nothing to bond to on a bare
  surface.
- The ion-dose channel (Bruce/Graves ion-processed skin) is retained; it is
  independently sourced and now sub-dominant where it should be.

## 4. Measured effect (0-D, audited top-band delivery 0.372, tilt 0.472°)

| condition | x_xl before | x_xl after | film growth before | after |
|---|---|---|---|---|
| **lip** (grazing ions) | 0.176 | **0.703** | 1.759 nm/s | **0.831 nm/s** |
| blanket (normal incidence) | 0.618 | 0.741 | 0 (steady) | 0 (steady) |

The 0-D lip reproduces the feature runs' audited 0.163 before the change, which
is what licenses the diagnosis. Against Krüger's run-average per-side closure
of **0.427 nm/s**, the lip goes from **4.1x too fast to 1.95x too fast** — the
single largest correction found in this investigation, and the blanket film
moves only 17 % in thickness (4.02e19 -> 3.33e19 atoms/m²), so no regime change.

## 5. What remains, and why it is not tuned here

The implemented form is the **minimal** transcription: one bond per deposition
event, two cells converted. Its exact ion-free steady state is `2(1-x) = x`,
i.e. **x = 2/3** (gated). His module allows *multiple* bonds per event — "a
maximum number of crosslink partners... (3 in the example depicted)" — which
would give `(1+m)(1-x) = x`, i.e. x = 0.8 at m = 3 and x -> 1 for large m. The
x ≈ 0.9 that reproduces his closure exactly sits in that family.

**The per-material maximum bond number is never tabulated in the thesis** — it
appears only as "3 in the example depicted in Figure 2.2". Adopting it would be
a judgement call on an undeclared constant, so it is recorded as `[VERIFY]`
rather than fitted. The residual 1.95x is bounded by this one undetermined
integer.

## 6. Gates

`tests/test_deposition_crosslinking.py` (5):

| gate | result |
|---|---|
| ion-free steady state equals the row stoichiometry, x = 2/3 | within 5 % |
| crosslinking survives grazing incidence (lip x_xl > 0.6) | 0.70 |
| crosslinked lip grows slower than fresh at identical flux | holds |
| `n_xl` never exceeds the film inventory (200 steps, aggressive dt) | holds |
| no creation without deposition or film | exactly 0 |

`tests/test_mixed_layer_mechanism.py::test_crosslink_brake_is_ion_dose_differential`
**repaired, not weakened**: it asserted `x_lip > 2·x_shadow` — that ion dose is
what crosslinks. That is the falsified direction. It now asserts the corrected
signature (shadowed surface crosslinks, `x_shadow > 0.6`; bombardment cannot
exceed it by more than 15 %) and keeps the growth differential unchanged.

## 7. Next

A confirmation run is now justified, and per the early-transient finding it
needs **12 s, not 60 s** (88 % of the defect is complete by t = 8 s), graded on
`closure/etch` by window against Krüger's 0.0310 — not on aperture at 60 s.
