# Cascade funnelling grade — the depth gate is not an over-delivery (2026-08-05)

Instrument: `scripts/cascade_funnelling_scan.py`. Straight-walled trenches
(0.09 µm opening, 0.85 µm mask — the Krüger cell) at a series of etch depths,
one frozen-geometry transport gather each with the cascade active, measuring
delivery to the etch-front band split into direct ions and cascade hot
neutrals. Chemistry evaluated separately at the measured delivery.

Motivating claim (overnight, `RESULTS_FLOOR_DELIVERY_2026-08-05.md`): oxide
recession is exactly proportional to ion flux and flat across a 50× range of
neutral delivery, so ml19's 21.9 / 16.7 nm/s implied **2.9× / 2.2× the source
ion flux** reaching the floor — an energetic over-delivery by the funnelled
cascade. This document tests that claim by measuring the delivery directly.

## Measured etch-front delivery (dx = 10 nm, CPU)

| oxide AR | direct ion (m⁻²s⁻¹) | hot neutral | ion/AR0 | hot/AR0 | hot/ion |
|---|---|---|---|---|---|
| 0 | 8.181e19 | 7.962e18 | 1.000 | 1.000 | 0.097 |
| 2 | 7.978e19 | 9.666e18 | 0.975 | 1.214 | 0.121 |
| 4 | 7.770e19 | 1.138e19 | 0.950 | 1.430 | 0.147 |
| 8 | 7.357e19 | 1.479e19 | 0.899 | 1.857 | 0.201 |

Boundary source ion flux (measured from the same boundary object): 1.20e20
m⁻²s⁻¹. Floor neutral delivery at AR 0 is 7.3 % of source (uniform across
CF/CF₂/C₂F₃/O — conductance through the 0.85 µm mask), 6.6 % at AR 8.

**The cascade does not over-deliver.** The etch front receives 0.68× the source
ion flux at AR 0 falling to 0.61× at AR 8, plus a hot-neutral component that
rises from 0.066× to 0.123×. Total energetic delivery is *flat to −1 %* across
AR 0→8.

### The over-delivery claim is refuted, and so was my first reading of it

Feeding those measured fluxes into a **1500 eV, normal-incidence** 0-D
evaluation gives 6.46 nm/s at AR 0 — far below ml19's 21.9 nm/s, which looks
like a factor ≈ 3.4 unaccounted for. That reading is wrong, and the reference
energy is why.

Relaxing the mixed layer on the *same gather* through the engine's own
per-atom path returns **19.15 nm/s at AR 0** — consistent with ml19. The
etch-front population measured on that band has **mean energy 3406 eV and mean
incidence cosine 0.768**, not 1500 eV at normal incidence: the published IEAD
carries substantial flux well above its nominal energy, and the band mixes the
flat floor with the corner fillet. 0-D checks against a single nominal energy
are therefore not a valid reference for this feature, and the overnight
inference of "2.9× the source ion flux" — which used exactly that reference —
does not survive.

| evaluation | rate at AR 0 |
|---|---|
| 0-D, 1500 eV, cos 1.0 (the overnight reference) | 6.46 nm/s |
| 0-D, measured E and cos (3406 eV, 0.768) | 53.5 nm/s |
| engine per-atom relaxation on the measured gather | **19.15 nm/s** |
| ml19 feature run | ~20–22 nm/s |

**Transport and chemistry are mutually consistent.** There is no missing factor
and no internal contradiction: the delivery we measure, run through the
chemistry we ship, reproduces the rate the feature runs produce. (The 0-D
single-energy columns bracket it, as Jensen requires — the per-atom path is the
one with the validated contract.)

## Grade against Huang's published funnelling numbers

Huang thesis (`tmp/pdfs/huang_thesis.txt`), verbatim:

- L5405–5407: *"This shadowing contributes to a decrease in ion flux to the etch
  front from 2.0 × 10¹⁵ to 0.3 × 10¹⁵ cm⁻²s⁻¹."* (oxide AR 0 → 40)
- L5408–5413: *"The flux of hot neutrals to the etch front increases from
  3.1 × 10¹⁵ to 8.0 × 10¹⁵ cm⁻²s⁻¹ as the etch depth increases from 0 to 480 nm
  (AR = 4). … decreases to 1.1 × 10¹⁵ cm⁻²s⁻¹ … due to diffusive scattering
  from the sidewalls and thermalization."*
- L5399–5402: *"Each strike of the etch front by a neutral particle (either hot
  or thermal) increments the flux-count."* (re-arrival counting)

| observable | Huang | petch | note |
|---|---|---|---|
| hot/ion at etch front, AR 0 | 1.55 | **0.097** | 16× low |
| hot-neutral rise AR 0→4 | 2.58× | **1.43×** | 1.8× low |
| ion flux AR 0→8 | strongly falling | **0.90** | nearly flat |

**The comparison does not transfer, and the reason is the source, not the
cascade.** Huang's ion flux decays because his beam is wide: he states
(L5403–5405) that *"the maximum incident angle of an ion which can directly hit
the etch front … decreases from about 4° to 1°"*, and his etch-front ion flux at
AR 0 is 2.0e15 of a 4.2e15 open-field flux (L5414–5415) — i.e. **≈ 52 % of his
ions arrive beyond 4.4°**, the geometric acceptance of his AR-13 photoresist.
The Krüger IEAD we consume is narrow: signed angles span ±2.86°, polar rms
0.84° after the P1a lift. A beam that narrow is barely shadowed at these aspect
ratios and produces few wall strikes, hence few hot neutrals. Both discrepancies
are one fact about the source distribution.

Two different reactors: Huang's funnelling curve is **not** a valid grading
reference for the Krüger cell. It remains valid for grading the cascade *given
a matched beam*, which is the open test (§ next steps).

## What the depth overshoot actually implicates

Krüger's own assessment of his mechanism, verbatim
(`tmp/pdfs/krueger_thesis.txt` L4884–4888):

> "These trends indicate that above a certain threshold energy the etch
> progression and the mask removal process are not ion starved, but rather
> **limited by neutral gas transport**. To some degree this trend is reproduced
> by the simulations where etch depth does increase with increasing low
> frequency power however the rate of increase is substantially sublinear.
> These outcomes indicate that **the effect of ion energy (for example in
> sputter yield or related processes) might be overestimated in the
> mechanism.**"

Our floor is **ion-limited** (rate ∝ ion flux, flat across 50× neutral
delivery). The experiment he describes is **neutral-transport limited**.

With the internal inconsistency removed, the depth gate decomposes into two
measured statements, neither of which is a transport defect:

1. **Absolute rate.** We etch ≈ 19 nm/s at low AR where his experiment averages
   825 nm / 60 s ≈ 13.8 nm/s including all attenuation — roughly 1.4× high, in
   the channel his own text says is likely over-estimated.
2. **Too little attenuation.** Our total energetic delivery is flat to −1 % over
   AR 0→8 because the beam we consume is ±2.86° and is barely shadowed at these
   aspect ratios. A process that is ion-limited (ours) with an unshadowed beam
   produces almost no ARDE; the experiment he describes is neutral-limited,
   where conductance falls with AR.

Both are **chemistry- and source-side**, reached by eliminating transport:
delivery is faithful, the cascade is faithful to its own published rule, and the
engine's chemistry reproduces its own feature runs from the measured delivery.

## Standing [VERIFY]

- Whether the Krüger Fig-4 IEAD digitisation is truncated at ±3° (a real
  collisional tail beyond the plotted axis would change the shadowing and is not
  visible in the digitised table). The figure axis range is not stated in the
  text.
- Matched-beam cascade grade: run this scan with a Huang-width beam and check
  whether hot/ion → 1.55 at mask AR 13 and the AR-4 peak appears at 2.58×. This
  tests the cascade mechanism itself rather than the source.

## Next

1. Matched-beam cascade grade (above) — free, decides whether the cascade
   mechanism is right independent of the source width.
2. The ion-limited vs neutral-limited crossover in our mixed layer: locate the
   neutral scale at which the floor becomes neutral-limited and compare with the
   flux ratios Krüger's Table I implies. If our crossover sits far from his, the
   ion-energy channel magnitude is the depth gate's owner, as his own text
   suggests.
3. Retire 0-D single-energy reference checks from the diagnostic toolkit, or
   force them to consume the measured (E, cos) population. The 1500 eV reference
   produced a 3.4x phantom in this pass and a "2.9x over-delivery" conclusion in
   the previous one.

## Addendum — completed sweep to AR 16, and the rate rises with depth

The full frozen-geometry sweep (`results/curated/cascade_funnelling/scan.json`):

| oxide AR | ion/AR0 | hot/AR0 | total/AR0 | hot/ion |
|---|---|---|---|---|
| 0 | 1.000 | 1.000 | 1.000 | 0.097 |
| 4 | 0.950 | 1.430 | 0.992 | 0.147 |
| 8 | 0.899 | 1.857 | 0.984 | 0.201 |
| 12 | 0.850 | 2.284 | 0.977 | 0.262 |
| 16 | 0.798 | 2.713 | 0.968 | 0.331 |

Over AR 0 -> 16 the direct ion flux falls 20 % and the cascade hot-neutral flux
rises 171 %, so **total energetic delivery falls 3 %**. The rising cascade term
very nearly cancels the shadowing term across the entire measured range.

Engine per-atom relaxation on those gathers gives floor recession
**19.15 nm/s at AR 0 and 24.30 nm/s at AR 4** — the etch rate *increases* 27 %
with depth. These are straight-walled trenches, so this is not taper focusing
(Huang's "small initial increase … due to tapering", L5578); it is the
hot-neutral term outrunning shadowing on a vertical wall.

This is the quantitative form of the ARDE statement in the main text: the
measured delivery produces no aspect-ratio slowdown at all over the range the
validation case spans, and the depth channel therefore has nothing to attenuate
it. Anti-ARDE over AR 0 -> 4 is a sharper falsifier than the +29 % endpoint
miss, and it is free to re-measure after any change to the beam or the
ion-energy channel.
